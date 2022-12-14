**In [another post](https://rolansen.github.io/radiometric_correction/local_histogram_correction.html), 
I discuss a few approaches for radiometrically correcting aerial imagery and evaluate their accuracy using real imagery. 
I figured it'd be nice to test them out on simulated landscapes to see how characteristics of the imagery, 
like the presence of objects like trees and houses, affect these methods' performances.  
Here I talk about how we'd begin to approach doing that sort of thing
using Gaussian random fields and random sets.
Simulations are ran with R and the libraries NLMR, spatstat, raster, sf, and fasterize.**

-----

The first thing I'd like to look at are the effects of spatial (auto)covariance. 
One way to do this is to repeatedly simulate a Gaussian random field (GRF) and adjust the mean and covariance function's parameters (sill, range, nugget). 

I'll talk about GRF's briefly--the next couple of paragraphs should be familiar to those who have worked with geostatistics, and much more detail can be found in Cressie (1993). Say we have some process *Y*(**s**), where the *d*-dimensional location vector **s** is in some study area **D** ⊂ ℝ*ᵈ* (usually *d* = 2, for the geography-minded). It's often advantageous to assume *Y*(**s**) has a Gaussian distribution, because this implies it can be completely described by its mean and covariance. If this is the case and *Y*(**s**) is defined for all **s** &isin; **D**, then *Y*(**s**) is a GRF. Furthermore, if for two arbitrary locations **sᵢ** and **sⱼ** we have cov(*Y*(**sᵢ**), *Y*(**sⱼ**)) = cov(*h*), where *h* = \|\|**sᵢ** - **sⱼ**\|\|, then *Y*(**s**) is stationary (properties don't depend on location) and isotropic (direction between two locations doesn't matter).

<div style="text-align: center">
  <figure>
      <img
       src="/assets/variogram.png"
       width="662"
       height="414"
     />
     <figcaption>Exponential variogram model. Sill = 1, nugget = 0.2, range = 10.</figcaption>
  </figure>
</div>\

For this type of GRF the covariance structure is often described with the (semi)variogram *&gamma;*(*h*) = cov(0) - cov(*h*), where cov(0) is the sill of the variogram, representing its upper bound as *h* approaches infinity. One common model for *&gamma;*(*h*) is the exponential model, which for *h* > 0 is defined as: 

&gamma;(*h*) = *a* + *&sigma;*²(1 - exp(-*h*/*r*)), 

and is equal to 0 otherwise. *a* describes variation at *h* = 0 and is called the nugget. *r* is the (effective) range, a value of *h* at which &gamma;(*h*) exceeds 95% of the sill cov(0) = *a* + *&sigma;*². *&sigma;*² is referred to as the partial sill.

<div style="text-align: center">
  <figure>
      <img
       src="/assets/gfs_with_different_params.png"
       width="621"
       height="540"
     />
     <figcaption>Simulated Gaussian fields all with exponential covariance functions and the same mean (10) and nugget (0.2), but varying partial sill and range. Each image is 100 by 100 pixels.</figcaption>
  </figure>
</div>\

GRF's are often invoked when working with spatial data--for example, kriging methods assume the surface of interest is a GRF. 

Let's try simulating a stationary, isotropic GRF. We can do this very easily with NLMR. Below I stick with image dimensions of 100 X 100 and the following parameters:
* mean = 0.22
* partial sill = 0.007
* nugget = 0.002 
* range = 10 units
{% highlight R %}
set.seed(1)

lindim <- 100 #width/height of image
gf_range <- 10
gf_partial_sill <- 0.007 
gf_nug <- 0.002
gf_mean <- 0.22

gfield <- nlm_gaussianfield(ncol=lindim, nrow=lindim, autocorr_range=gf_range, mag_var=gf_partial_sill, nug=gf_nug, mean=gf_mean, rescale=FALSE)
gfield <- reclassify(gfield, matrix(c(-Inf,0,0, 1,Inf,1), ncol=3, byrow=TRUE)) #ensures values are between 0 and 1
{% endhighlight %}
I use a low mean and variance and force values to be greater than 0 and less than 1 because I'd like this field to be similar to visible range surface reflectance measurements. 

We're using an exponential variogram/covariance model, which is the only kind allowed by NLMR's *nlm_gaussianfield()* function. Under the hood *nlm_gaussianfield()* is using the RandomFields library, and if we'd wanted to use a different covariance model we could have just as easily used RandomFields to simulate the field, converted the resulting object to a matrix, truncated values so they're between 0 and 1, and then converted the matrix to a raster. For now, though, I'll just use *nlm_gaussianfield()* to keep things simple.  

-----

Most aerial images include clearly distinguisable "objects" such as trees and buildings. Intuitively, we know these objects are not randomly distributed. [In the post I mentioned at the beginning of this one](https://rolansen.github.io/radiometric_correction/local_histogram_correction.html), I found that the accuracy of each radiometric correction method varied with the fine-scale heterogeneity of each image--where these kinds of objects tended to be present, every method seemed to do worse. 

How can we simulate objects like this? The most straightforward way to do so would be to generate a random set, more specifically a Boolean model. Informally, this means we'd like to simulate a point process, making up the "germs" of the set, and then for each germ create a corresponding geometry, or "grain." This set will be used to mask out our GRF and substitute other values in its place. I don't expect the shapes of the grains to be very relevant, so to keep things as simple as possible they'll just be disks with a random radius centered on the germs. For more detailed and formal descriptions of random sets, see Chiu et. al. (2013) and once again Cressie (1993). 

The first step is to decide how to vary the intensity of the germ-generating process throughout our study area. One way we could do this is to simulate an inhomogenous Poisson process: 

*Pr*(*X* = *k*) = *&lambda;*(**s**)*ᵏ* exp(*-&lambda;*(**s**)) / *k*!, 

where *k* is some number of events (points) and *&lambda;*(**s**) is our spatially varying intensity (mean rate of points per unit area). 
*&lambda;*(**s**) could be some function of space that's independent of our GRF, e.g. *&lambda;*(**s**) is proportional to the y coordinate. It might make sense to have the intensity depend on the GRF, though. Doing this would be sort of like having the intensity depend on land cover type (we'll divide our landscape into discrete patches as we go forward, which may make this analogy work better). 

To accomplish this, we can make the intensity a Gaussian function of the GRF. Adjusting the mean will control the GRF value corresponding to the maximum intensity, and adjusting the variance will control how clustered the process is.

<div style="text-align: center">
  <figure>
      <img
       src="/assets/intensity_dists.png"
       width="662"
       height="414"
     />
     <figcaption>Two Gaussian probability density functions. Left: mean=0.3, standard deviation=0.15. Right: mean=0.6, standard deviation=0.03.</figcaption>
  </figure>
</div>\

Here's how we can make a raster of the intensity (note that we use a factor to keep the expected number of points in each pixel from being too high):
{% highlight R %}
max_intensity_val <- 0.3
stdev_intensity <- 0.15
intensity_factor <- 0.005
intensity_image <- gfield
values(intensity) <- intensity_factor * dnorm(as.vector(gf_model), mean=max_intensity_val, sd=stdev_intensity)
{% endhighlight %}

<div style="text-align: center">
  <figure>
      <img
       src="/assets/intensity_raster.png"
       width="621"
       height="540"
     />
     <figcaption>Intensity raster.</figcaption>
  </figure>
</div>\

Next we'll use spatstat to simulate the germs as a Poisson process:
{% highlight R %}
intensity_im <- as.im(as.matrix(intensity_raster))
germs <- rpoispp(lambda = intensity_im)
{% endhighlight %}

<div style="text-align: center">
  <figure>
      <img
       src="/assets/germs.png"
       width="621"
       height="540"
     />
     <figcaption>Germs (black) overlain on the GDF. Note how the germs tend to occur on pixels where the value is around 0.3.</figcaption>
  </figure>
</div>\

Since the intensity surface is a random field derived from our first GRF, we can think of our Poisson process as a Cox process, which is just a Poisson process where the intensity is random.

Now we can create the grains. I'll keep the grain radius distribution uniform, so any value between 2 units and 4 units is equally likely for all germs. We'll sample from this distribution so we have one value for each germ, then create the grains using *sf::st_buffer()*:
{% highlight R %}
min_radius <- 2
max_radius <- 4
radii <- runif(n=germs$n, min=min_radius, max=max_radius)
sf_germs <- st_as_sf(germs)[1:germs$n+1,] #st_as_sf() will add the window as the first row
sf_germs$radius <- radii
sf_grains <- st_buffer(sf_germs, dist=sf_germs$radius)
{% endhighlight %}

<div style="text-align: center">
  <figure>
      <img
       src="/assets/grains.png"
       width="621"
       height="540"
     />
     <figcaption>Grains.</figcaption>
  </figure>
</div>\

Finally, we can use the random set to replace values in the original GRF. The replacement values will be drawn from a Gaussian noise image, i.e. a GRF with no spatial covariance: 
{% highlight R %}
noise_dist_mean <- 0.05
noise_dist_sd <- 0.02
noise_image <- gfield
values(noise_image) <- rnorm(ncell(noise_image), mean=noise_dist_mean, sd=noise_dist_sd)

grain_raster <- fasterize(sf_grains, noise_image)
gfield_masked <- mask(gf_model, grain_raster, updatevalue=0, inverse=TRUE)
noise_image_masked <- mask(noise_image, grain_raster, updatevalue=0)
landscape <- gfield_masked + noise_image_masked
{% endhighlight %}

<div style="text-align: center">
  <figure>
      <img
       src="/assets/landscape.png"
       width="621"
       height="540"
     />
     <figcaption>GRF with values replaced by Gaussian noise where pixel centroids are in the random set.</figcaption>
  </figure>
</div>\

Note that we use fasterize to convert the grains to a raster sampled to the resolution of the other images, which speeds up raster::mask(), and that the noise GRF has little variation compared to our first one.

-----

There's still some work to do before we can test out our radiometric correction methods. 
We're in a good position now, though, to take a step back and see how efficiently we can run these simulations. 
I plan on doing at least a few hundred iterations for each combination of parameters that I decide to work with
and working with much larger rasters than the ones I did for this post,
so computation time will be a relevant issue.

I ran 100 simulations for images with 100-by-100, 300-by-300, and 1000-by-1000 pixels
Parameters such as the Gaussian field mean, etc., are kept constant here, but in practice I'll vary them.

The table below shows, for each image size, the elapsed time various steps took in total. 
For example, the first cell of the first row shows how many seconds GRF simulation took, in total, 
for all 100-by-100-pixel runs.

|   | Gaussian field simulation | Germ intensity image creation | Germ simulation | Grain simulation | Image masking |
| -------------  | -------------  | -------------  | -------------  |-------------  |-------------  |
| 100-by-100 | 4.16 s | 0.34 s | 0.45 s | 2.75 s | 1.90 s |
| 300-by-300 | 29.86 s | 1.44 s | 1.96 s | 13.68 s | 4.11 s |
| 1000-by-1000 | 323.00 s | 8.55 s | 14.11 s | 86.03 s | 16.06 s |

\
While all steps seem to take longer as image size increases, GRF simulation and grain simulation stand out. It's known that GRF's can be more efficiently approximated with Gaussian Markov random fields (GMRF's), and it's possible that the grain raster could be generated more efficiently by running *raster::cellFromXY()* and then repeatedly running *raster::adjacent()*. I plan on trying both of these approaches going forward. 

-----

My next post on this topic will go over how we can divide our landscape into something like a patch mosaic with each patch containing values drawn from one of several GRF's. Then, I'll address the efficiency issues discussed above and also introduce nonstationarity to our GRF's by using the popular INLA approach for approximating GRF's with GMRF's. After all this we'll be in a good position to start evaluating our radiometric correction methods.

-----

References

Chiu, S. N., Stoyan, D., Kendall, W. S., & Mecke, J. (2013). *Stochastic geometry and its applications*. John Wiley & Sons.

Cressie, N. (1993). *Statistics for spatial data*. John Wiley & Sons.

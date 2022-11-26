Aerial photographs can be helpful sources of environmental data when fine detail is needed or long-term change is of interest. 
Many aerial images date back to the early 20th century, and their grain size isn't rivaled by any satellite imagery that can be worked with for free 
(or without a security clearance). 
Aerial imagery has a reputation for being outmoded (Morgan et. al., 2010), but nevertheless is still being used in industry and research today. 
An interesting example is given by Suraci et. al. (2020); the authors used telemetry data and classified aerial imagery to show how patch size, 
patch configuration, and building density drives puma movement in developed landscapes. 

Unfortunately, digital numbers (DN's) of digital aerial images are generally physically meaningless, 
so they can't be used to accurately estimate things like spectral indices. 
Also, there are usually systematic differences in DN's between images making up a time series due to contrasting illumination conditions 
and image processing procedures, limiting aerial imagery's capacity for supporting long-term, fine-scale landscape change monitoring. 

We can solve these problems by radiometrically correcting the aerial images we're working with, 
so that pixels represent (surface or top-of-atmosphere) reflectance estimates. 
Histogram matching is one simple and intuitive method for doing this sort of thing. 
Given some image with values we wish to change the distribution of, say an aerial image, 
and another image with a distribution that we want to map the other image's values to, like a Landsat image, 
we can map the values of the first image with the equation 

*y* = *Cₛ⁻¹*(*Cₐ*(*x*)),

where *x* is an input value from the first image, *Cₐ* is the cumulative distribution function (CDF) of the first image, 
*Cₛ⁻¹* is the inverse cumulative distribution function (iCDF) of the second image, and *y* is a target value (Richards, 2016).

Other simple empirical methods for radiometric correction, like dark pixel subtraction and pseudo-invariant pixel-based methods, 
assume a linear relationship between reflectance values and the DN's to be corrected, 
which may not be realistic ([Campbell & Wynne?]), 
and also require the presence of bodies of water and pseudo-invariant surfaces in the scene of interest, respectively. 
Histogram matching doesn't assume a linear relationship between DN's and reflectances, and only requires an image pair. 
Empirical methods in general also have the advantage of not requiring a unique training set for each unique aerial imagery collection/satellite combination, 
e.g. NAIP and Landsat versus Vexcel and Sentinel-2. 
These kinds of training sets may be necessary to get the best performance out of machine learning-based approaches. 
Several other methods for radiometric correction, particularly physical methods, require the use of field data, 
but this is not always practical depending on the user's available resources, especially if the aerial images were collected a few decades ago.

Histogram matching does have some obvious drawbacks. It works well if the relationship between aerial DN's and satellite reflectances is constant throughout the area of interest. Despite common measures such as BRDF correction, there are several reasons why this may not be true. One is that aerial images are often delivered as mosaics, and seamline data which allows us to separate the original images making up the mosaic may not always be available. Often these images will have been collected at different times of day, not to mention different dates (more on that below). Color balancing of the original images while mosaicing exacerbates things.  

Can we come up with a way to do histogram matching that accounts for these inconsistences? I think so. A starting point would be to just define a grid over the area of interest, then within each grid cell match our aerial image to a multispectral satellite image captured around the same time. I'll call this approach "localized traditional histogram matching" from here onwards.

We'll work with imagery taken over Fort Worth, Texas. The aerial imagery I want to correct is a subset of a NAIP mosaic which can be downloaded [here](https://nrcs.app.box.com/v/naip/file/769545426773). The images making up this mosaic were recorded in late 2020, and it's effectively 6-bit due to SID compression. We'll use a Sentinel-2 surface reflectance image that has a similar vintage (L2A_T14SPB_A027708_20201011T171836) as our reference. I'll be working with the red, green and blue bands of each image, i.e. the first three of the aerial and bands 4, 3, and 2 from the Sentinel-2 image. Each corresponding band pair will be treated separately from the others. 

*fig1: 300 m grid overlaying NAIP subset*

It seems very possible that we'd see some clear breaks in the image along the borders of many grid cells, if we applied localized traditional histogram matching. Trying it out on my study area with a length of 300 meters for the width and height of each grid cell shows that this is the case. 

*fig2: discontinuities due to differences in local histograms*

Increasing the granularity of the grid lessens these effects, as I discuss below. 

How could we address these discontinuities? Adaptive histogram equalization, a popular technique for image contrast enhancement, addresses this sort of thing by finding the 4 nearest (2 nearest for pixels near edges, nearest for pixels near corners) grid cells to each pixel, applying histogram equalization functions defined for each grid cell to the pixel, then interpolating the outputs based on the pixels' proximities to each grid cell center. It's pretty straightforward to do something analogous with histogram matching--let's call this approach "adaptive histogram matching." [Here's an R script](https://github.com/rolansen/rolansen.github.io/blob/main/code/ahm_no_subgrid.R) defining a function that will do this for us. It requires the libraries [terra](https://cran.r-project.org/web/packages/terra/index.html) and [sf](https://cran.r-project.org/web/packages/sf/index.html). After running the script, the user can correct their imagery with a command like this: 
{% highlight R %}
ahm(aoi_poly, aerial_path, sat_paths, c_region_size, grid_lindim_length, out_name)
{% endhighlight %}
The arguments are as follows:
* *aoi_poly* is an sf-tibble or sf-data.frame containing one polygon which contains each pixel from the aerial image we'd like to correct
* *aerial_path* is the path to an aerial image file. This and each file listed in *sat_paths* can take any GDAL-readable format
* *sat_paths* is a character vector listing the paths for each satellite band file. Each path should be listed in the same order as the corresponding bands from the aerial image file
* *c_region_size* is the length of the width and height for the "computational region" of each grid cell, in the units of the aerial image's coordinate reference system. Pixels falling in the computational region, which is centered on the corresponding grid cell's centroid, will be used to find the distribution of pixel values for that cell.
* *grid_lindim_length* is the length of each grid cell's width and height, in the units of the aerial image's coordinate reference system.
* *out_name* is the path of the output file, without the extension. The output file will always be a GeoTiff. 

In the examples below, *c_region_size* and *grid_lindim_length* are set to be identical.

In addition to writing a file for the corrected image, the function will output two other files: an Rds file containing a DataFrame that lists the satellite iCDF and aerial CDF corresponding to each grid cell, and a GeoPackage containing a layer for the grid used to find the aforementioned iCDF's/CDF's as well as a layer representing a separate grid that the function uses to identify, for each aerial pixel, the neighboring cells from the first grid.

*fig3: the two grids* 

The script is very much a work in progress, and likely has several bugs. I'll improve flexibility and efficiency going forward. Things I plan on doing in the future include making the aoi_poly and output files optional, putting in an option for returning a SpatRaster object, allowing both the aerial imagery and the satellite imagery to have their respective bands placed in one or multiple files, adding behavior allowing for a number of output bands other than 3, and making a faster interpolation step. I'll update this post as I do these things. I'd also like to try implementing this method with Google Earth Engine. 

Let's try applying adaptive histogram matching to the scene I described above. I used a grid cell and computational region width/height of 300 meters. We can see that the discontinuities have been smoothed over.

*fig 4: fig2 but from adaptive histogram matching*

Tuominen & Pekkinaren (2004) introduced another simple, empirical approach for radiometrically correcting aerial imagery. They used the formula
[corrected aerial = (mean sat in group) / (mean aerial in group) * (uncorrected aerial)], 
where [uncorrected aerial] is some pixel in the aerial image, [mean sat in group] is the mean value of pixels in the satellite image which are part of some "group" defined for the pixel of interest, [mean aerial in group] is the mean value of aerial pixels in the group, and [corrected aerial] is the corrected value of the pixel of interest.
They defined a "group" to be any one of the following: a) a satellite pixel and all coincident aerial pixels, b) all pixels within some radius of the pixel of interest, c) segments of the aerial image. They found the first method didn't work well. The process of image segmentation would also introduce its own set of problems, since the best parameters for most segmentation algorithms will likely vary from image to image. However, we can easily try out something like the radius-based method--I did this here, although rather than working with a circle, as they suggest, I just used the focal() function from the terra package to run a moving window across the imagery. 

...

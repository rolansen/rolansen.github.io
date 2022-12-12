**In another post, 
I discuss a few approaches for radiometrically correcting aerial imagery and evaluate their accuracy using real imagery. 
I figured it'd be nice to work with simulated landscapes to see how characteristics of the imagery, 
like the presence of objects like trees and houses, affect these methods' performances.  
Here I talk about how we'd begin to approach doing that sort of thing
using Gaussian random fields and random sets.
Simulations are ran with R and the libraries NLMR, spatstat, raster, and sf.**

-----

The first thing I'd like to look at are the effects of spatial (auto)covariance. 
One way to do this is to repeatedly simulate a Gaussian random field (GRF) and adjust the mean and covariance function's parameters (sill, range, nugget). 

I'll talk about GRF's briefly--the next couple of paragraphs should be familiar to those who have worked with geostatistics, and much more detail can be found in Cressie (1993). Say we have some process *Y*(**s**), where the *d*-dimensional location vector **s** is in some study area **D** ⊂ ℝ*ᵈ* (usually *d* = 2, for the geography-minded). It's often advantageous to assume *Y*(**s**) has a Gaussian distribution, because this implies it can be completely described by its mean and covariance. If this is the case and *Y*(**s**) is defined for all **s** &isin; **D**, then *Y*(**s**) is a GRF. Furthermore, if for two arbitrary locations **sᵢ** and **sⱼ** we have cov(*Y*(**sᵢ**), *Y*(**sⱼ**)) = cov(*h*), where *h* = \|\|**sᵢ** - **sⱼ**\|\|, then *Y*(**s**) is stationary (properties don't depend on location) and isotropic (direction doesn't matter).

*plot of exponential variogram model.*

For this type of GRF the covariance structure is often described with the (semi)variogram &gamma;(*h*) = &sigma;² - cov(*h*), where &sigma;² is the variance. One common model for &gamma;(*h*) is the exponential model: 

&gamma;(*h*) = &sigma;²(1 - exp(-3*h*/*r*)) + *a*. 

*a* describes variation at *h* = 0 and is called the nugget. *r* is the range, a value of *h* at which &gamma;(*h*) approaches an upper bound referred to as the sill. The sill is equal to *a* + &sigma;².

GRF's are often invoked when working with spatial data--for example, kriging methods assume the surface of interest is a GRF. 

Let's try simulating a stationary, isotropic GRF. We can do this very easily with NLMR:

*code, plots*

I use this low mean and variance and force values to be greater than 0 and less than 1 because I'd like this field to be similar to typical visible range surface reflectance measurements. Below I stick with image dimensions of 100 X 100 and the following parameters:
* mean = 0.1
* variance = 0.25
* nugget = 0.02 
* range = 10 units
We're using an exponential variogram/covariance model, since that's the only kind allowed by NLMR's nlm_gaussianfield() function. Under the hood nlm_gaussianfield() is using the RandomFields library, and if we'd wanted to use a different covariance model we could have just as easily used RandomFields to simulate the field, convert the resulting object to a matrix, truncated values so they're between 0 and 1, and then converted the matrix to a raster. For now, though, I'll just use nlm_gaussianfield() to keep things simple.  

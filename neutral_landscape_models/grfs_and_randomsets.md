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

*plot of exponential variogram model. List the parameter values.*

For this type of GRF the covariance structure is often described with the (semi)variogram &gamma;(*h*) = &sigma;² - cov(*h*), where &sigma;² is the variance. One common model for &gamma;(*h*) is the exponential model &gamma;(*h*) = &sigma;²(1 - exp(-3*h*/*r*)) + *a*. *a* describes variation at *h* = 0 and is called the nugget. *r* is the range, a value of *h* at which &gamma;(*h*) approaches an upper bound referred to as the sill. The sill is equal to *a* + &sigma;².

GRF's are often invoked when working with spatial data--for example, kriging methods assume the surface of interest is a GRF. 

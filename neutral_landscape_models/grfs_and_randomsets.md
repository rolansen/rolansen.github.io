**In (this post)[https://github.com/rolansen/rolansen.github.io/blob/main/radiometric_correction/local_histogram_correction.md], 
I discuss a few approaches for radiometrically correcting aerial imagery and evaluate their accuracy using real imagery. 
I figured it'd be nice to work with simulated landscapes to see how characteristics of the imagery, 
like the presence of objects like trees and houses, affect these methods' performances.  
Here I talk about how we'd begin to approach doing that sort of thing
using Gaussian random fields and random sets.
Simulations are ran with R and the libraries NLMR, spatstat, raster, and sf.**

-----

The first thing I'd like to look at are the effects of spatial (auto)covariance. 
We can do this by repeatedly simulating a Gaussian random field (GRF) and adjusting the mean covariance parameters (sill, range, nugget). 

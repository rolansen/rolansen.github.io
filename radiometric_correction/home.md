Here I'll write a few posts describing a method I made for the radiometric correction of aerial imagery.
The goal is to do something about as conceptually simple and easily applicable as histogram matching, but which accounts for inconsistencies in the relationship between aerial DN's and measured reflectances due to BRDF and other effects.
We achieve this basically by deriving a set of histogram matching functions throughout the scene of interest, then applying a few to each pixel and interpolate the results.
A script which demonstrates this procedure can be found [here](/code/ahm_no_subgrid.R). The script is very much a work in progress, but the functions it contains work. So far, the method has compared very favorably to traditional histogram matching. 
I'll give much more detail in the very, very near future.

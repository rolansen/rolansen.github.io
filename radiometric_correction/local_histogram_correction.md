Aerial photographs can be helpful sources of environmental data when fine detail or [long-term records] is needed. 
Many images date back to the early 20th century, and their grain size isn't rivaled by any satellite imagery that can be worked with for free 
(or without a security clearance). 
Aerial imagery has a reputation for being outmoded (Morgan et. al., 2010), but nevertheless is still being used in industry and research today. 
An interesting example is given by Suraci et. al. (2020); the authors used telemetry data and classified aerial imagery to show how patch size, 
patch configuration, and building density drives puma movement. 

Unfortunately, digital numbers (DN's) of digital aerial images are generally physically meaningless, 
so they can't be used to accurately estimate things like spectral indices. 
There are also usually systematic differences in DN's between recorded over the same area due to contrasting illumination conditions 
and image processing procedures, limiting  aerial imagery's capacity for supporting long-term, fine-scale landscape change monitoring. 

We can solve these problems by radiometrically correcting the aerial images we're working with, 
so that pixels represent (surface or top-of-atmosphere) estimates. 
Histogram matching is one simple and intuitive method for doing this sort of thing. 
Given some image with values we wish to change the distribution of, say an aerial image, 
and another image with a distribution that we want to map the other image's values to, like a Landsat image, 
we can map the values of the first image with the equation 
[y = iCDF_sat(CDF_aerial(x))],
where [x] is an input value from the first image, [CDF_aerial] is the cumulative distribution function (CDF) of the first image, 
[iCDF_sat] is the inverse cumulative distribution function (iCDF) of the second image, and [y] is a target value (Richards, 2016).

Other simple empirical methods for radiometric correction, like dark pixel subtraction and pseudo-invariant pixel-based methods, 
assume a linear relationship between reflectance values and the DN's to be corrected, 
which may not be realistic ([Campbell & Wynne?]), 
and also require the presence of bodies of water and pseudo-invariant surfaces, respectively. 
Histogram matching doesn't assume a linear relationship between DN's and reflectances and only requires an image pair. 
Empirical methods in general also have the advantage of not requiring a unique training set for each unique aerial imagery collection/satellite combination, 
e.g. NAIP and Landsat versus Vexcel and Sentinel-2. 
These kinds of training sets may be necessary to get the best performance out of machine learning-based approaches. 
Several other methods for radiometric correction require the use of field data, 
but this is not always practical depending on the user's available resources, especially if the aerial images were collected a few decades ago.

"Title"
"brief summary. Caveat about “lidar data” == discrete return aerial lidar point clouds"

In the past few years free, high-quality lidar datasets have become much more common than they were a decade ago. A good portion of the planet now has at least one representative lidar survey. Surveys will become more widespread and regular, and quality will keep improving.
This increase in the availability of lidar data has caused it to become more "widely used" in the GIS world. LAStools and PDAL are "in wide use by GIS professionals", and ArcGIS and QGIS have both been adding more options for visualizing and processing lidar point clouds. 

"screenshot of recent QGIS release showing lidar support. Give credit to source, probably Lutra"

However, tools for integrating lidar data with other GIS data types (vector & raster) can still feel somewhat limited. Say we wanted to identify the returns in or around features from a polygon layer, then compute a new field for this layer based on the returns corresponding to each polygon? What if we wanted to ""raster example. Limit returns to just those hitting certain classes in a land cover raster?""? 

It seems relatively difficult to do these things with the tools I listed above. 
R could be helpful here via the package [lidR](https://cran.r-project.org/web/packages/lidR/index.html). Another option is Python, which is perhaps more ubiquitous at this point and has the huge advantage of letting us work with point clouds as numpy arrays. Here, I’ll go over a hypothetical use case which will show how easy, flexibly, and efficiently we can work with lidar point clouds using tools from the Python ecosystem. 

The libraries that we’ll need are laspy, the Python bindings for the Rust crate laz-rs (see the [laspy installation instructions](https://laspy.readthedocs.io/en/latest/installation.html), NumPy, SciPy, and GeoPandas. This may be easier to follow along with if the reader has worked with NumPy and GeoPandas before, but I don’t think they have to be too familiar with those packages. 

$$
E = mc^2
$$

<<[old]
Tools for implicating lidar data with other GIS data types (vector & raster) are still relatively limited. ArcGIS & QGIS have stuff, but … . LasTools and PDAL are great, but. (come up with a few hypothetical scenarios). What if we want to do ___, or ____?
Luckily, we have Python is at this point universal among those who may not be professional software developers (analysts, technicians, data scientists…). Here I’ll show how to use it to do lidar stuff. Laspy (+ lazrs python port), numpy, scipy, and geopandas make working with vector/lidar easy, flexible, and fast. >>

-----

Let’s say we’ve got a set of polygons representing tree canopies and we want to compute each tree’s maximum height and its leaf area index (LAI). Aerial lidar can help us easily compute both of these. 

For this post I drew a set of polygons over tree canopies around Soldier’s Circle [link], a parkway in Buffalo, New York. The majority of these trees are American elms. I used Bing imagery as a reference.

[Google StreetView screenshot, hollow-fill polygons over Bing imagery].

You can download these polygons in GeoJSON format here [link]. 

The lidar we’ll work with is from a 2019 survey. Here we’ll just use a single tile. USGS provides this data in laz format, and you can download it here [link].

Both the lidar and the polygons use the coordinate reference system specified by EPSG:6350. 

If we were to do something like this in practice, we’d probably be segmenting trees rather than working with polygons. I’ll go over how we can do this in a future post. We should also keep in mind that some of the polygons I drew will include some returns from rootftops. I’ll discuss how to handle that in the future, as well.

{explain hypothetical scenario. Mention how tree polygons are unrealistic, but we’ll go over segmentation in a future post. Show pictures. Link to data}

-----

First, let’s change directories to the folder in which we’re keeping the data, read our polygon layer as a GeoDataFrame, and instantiate the fields we’d like to calculate with the point cloud:
{os.chdir(), read polygons, create columns.}

Next, we’ll read our lidar dataset and convert it to an ndarray. laspy makes it very easy to read lidar data, and while laspy’s LasData objects already represent point clouds as numpy arrays (see the point.array attribute), I find that it’s often easier to just convert the dimensions of interest directly, since so many libraries are able to work with ndarrays. Since we want the scaled [link {to laspy page explaining what scale means}] versions of the x, y, and z dimensions, which aren’t included in las.point.array, we’ll make a new array which includes these dimensions, the class code dimension (which we’ll use to find heights above ground elevation), the scan angle (which we’ll use to calculate LAI), and an empty column that we’ll fill in with heights. 
[code]
In practice, we might want to make an array from more than one lidar file, in which case something like this could be helpful: np.append(np_las_1, np_las_2, axis=0)
{read lidar and convert to numpy. As an aside, mention that in practice we might want to make an array from more than one lidar file, in which case np.append can be helpful: np.append(np_las_1, np_las_2, axis=0). Explain why we “convert” to numpy array (laspy dimensions are numpy arrays, but it’ll be easier to just work with an ndarray directly).} 

Next, we’ll make a spatial index. k-D trees are popular spatial indexes for point clouds because they’re good for nearest neighbor queries (what are the n nearest neighbors?) and range queries (what are the neighbors within a distance r from some point?) on point datasets. We’ll use scipy’s implementation [link]. Since we’re only concerned with finding which returns have (x,y) coordinates that fall in polygons, we’ll just build the tree for those dimensions. For more on k-D trees, Xiao (2016) gives a great introduction for GIS folks. 
[code]
[{picture of k-d tree. Take from Xiao, p 81?}]
{kd-trees. Very briefly explain the idea, mention they’re popular for point clouds b/c of what they’re good at, mention Xiao (2016) as a good reference. Credit laspy documentation here. For this section, mention how we can use multithreading to make things faster for both kd-tree construction and numpy stuff, although we don’t do that here}

Now we’re ready to calculate our variables of interest. For each polygon:
1) Do a range query to find the points with (x,y) coordinates near the polygon. I like to use a range that’s slightly larger than the minimum enclosing circle of the minimum bounding rectangle of the polygon, since it’s guaranteed to include each point in the polygon and also makes it more likely that we’ll include many ground returns. Here, we’ll use a buffer of 2 meters for the minimum enclosing circle.
[code]
[picture of buffered MEC, MEC, MBR, and polygon]
2) Next, we’ll interpolate ground returns and use the derived surface to compute height above ground for each point satisfying the range query. USGS lidar datasets usually have good ground classifications, so we’ll just use theirs. We’ll use TIN interpolation, which is fast if we’re only considering a few hundred points at a time like we are here, and will accurate enough for our purposes (separating returns near or on the ground from those which probably correspond to tree canopies). Again, we’ll take advantage of scipy for this.
[code]
3) Now we’re ready to find the points which fall in our polygon. One approach for doing this that seems appealing is the shapely.contains() method [link shapely.contains — Shapely 2.0b2 documentation], since our polygons are already represented as shapely objects by GeoPandas. However, I’ve found that it takes too long to convert an ndarray to a list of shapely points. What works better is to convert the polygon to a list or ndarray of coordinates, then use some other function for doing point-in-polygon tests. Here I’ll use a function mostly lifted from Xiao (2016), but which doesn’t consider returns intersecting polygons’ line segments or vertices. This works well enough for our purposes, but if you’d like more speed you can check out this implementation of the winding number algorithm [link] which leverages numba—as shown by this stackexchange post [link], this can make a huge difference. For more on point-in-polygon tests, Xiao (2016) is again a good reference.
[code]
[picture, if we have time]
4) From here, we can easily find the maximum height above ground.
[code]
5) LAI can also be easily calculated now. Often the distribution of light throughout canopies is described in a method similar to the Beer-Lambert law for light attenuation through a homogenous medium (Jones, 2013). More specifically, this looks something like:
$$
E = mc^2
$$
Where I is the irradiance (i.e., radiant flux) at the ground surface [W m^-2], I_0 is the irradiance at the top of the canopy, L is the leaf area index, and k is an “extinction coefficient” representing the ratio of the area of shadows cast by leaves to the actual area of the leaves (Jones, 2013). Solving for LAI gives us
L = -1/k * ln(I / I_0).
For a spherical leaf angle distribution, meaning all leaves have a uniform probability for all angles from the horizontal plane (zenith angles %theta%) %in% [0°, 90°], we have k = 0.5. Following Richardson et. al. (2009), we’ll substitute I with the total number of ground returns R_g and I_0 with the total number of returns in the polygon R_t, and also model the effects of lidar scanning angle using Lambert’s cosine law I = I_0*cos(%theta%). This gives us the model:
L = -cos(%mean_theta_lidar%) / 0.5 * ln(R_g/R_t)
Where %mean_theta_lidar% is the mean lidar scanning angle of all returns in the polygon. Here’s the code for implementing this model:
[code]
Note we consider a return to be a ground return if its height above ground is less than or equal to 5 cm. Also note that the scan angle dimension doesn’t actually provide us with angles, and we multiply it by a factor to get the actual angle. (see this [link https://github.com/ASPRSorg/LAS/issues/41#issuecomment-344300998] comment). {TODO: do this calculation beforehand, and copy/paste this part above}.
 

{Explain procedure for each polygon:
	1) do range query around polygon. Explain why I like MEC + buffer of MER.
	2) do TIN interpolation of ground returns. Mention recently collected USGS point clouds generally have good ground classifications. Explain why (fast if we just look at a few points at a time, and should be accurate enough for our purposes.
	3) do point-in-polygon test. Mention we’re avoiding shapely methods because converting numpy points to shapely objects takes too long. Very briefly explain ray tracing method, again mention Xiao (2016) as a good reference and mention that below we use a modified version of his code, but don’t consider points on line segments or verticies of polygons. Mention stackexchange post which goes over quickest ways to do point-in-polygon tests with python + numpy + numba (best is winding number + numba?).
	4) max height this is just the maximum height above ground with (x,y) coordinates falling in the polygon 
	5) calculate LAI. Go in deep here.
Then show the code}

{Explain results (run time, etc), make some comments, show some pictures.}

If you’d like to reproduce exactly what we did here, here’s a script [link] which puts together everything we did above. Again, the polygons can be found here [link] and the laz file here [link].
{Here’s the whole code, if you’d like. Here’s the geojson and here’s the laz files, if you’d like to completely reproduce}

{Note “horizontal error” of polygons and lidar?}

{One thing I didn’t go over here is multithreading.}

-----

Next post: rasterio, vegetation and building classification…

-----

References

Jones (2013). Plants and microclimate, Chapter 2.

Richardson, Moskal, Kim (2009). Modeling approaches to estimate effective leaf area index from aerial discrete return LIDAR.

Xiao (2016).

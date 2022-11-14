**Title**
**brief summary. Caveat about “lidar data” == discrete return aerial lidar point clouds**

-----

Free, high-quality lidar datasets have become much more common in the last few years than they were a decade ago. A good portion of the planet now has at least one representative lidar survey. In the future surveys will become more widespread and regular, and quality will keep improving.
This increase in the availability of lidar data has caused it to become more commonly utilized in the GIS world. **The trend** has been **bolstered** by the awesome tools that are available nowadays. LAStools and PDAL are often used by GIS professionals, WhiteboxTools has great support for lidar data, and ArcGIS and QGIS have both been adding more options for visualizing and processing lidar point clouds. 

**TODO: screenshot of recent QGIS release showing lidar support. Give credit to source, probably Lutra**

However, at the moment tools for integrating lidar data with other GIS data types (vector & raster) can still feel somewhat limited. What if we wanted to do something like identify the returns in or around features from a polygon layer, then compute a new field for this layer based on the returns corresponding to each polygon? What if we wanted to limit our data to returns which fall in certain classes in a land cover raster? 

It seems relatively difficult to do these things with the software I listed above. R could be helpful here via the [lidR package](https://cran.r-project.org/web/packages/lidR/index.html). Another option is Python, which is perhaps more ubiquitous at this point and has the huge advantage of letting us work with point clouds as numpy arrays. Here, I’ll go over an example which will show how easily, flexibly, and efficiently we can work with lidar point clouds using tools from the Python ecosystem. 

The libraries that we’ll need are laspy, the Python bindings for the Rust crate laz-rs (see the [laspy installation instructions](https://laspy.readthedocs.io/en/latest/installation.html), NumPy, SciPy, and GeoPandas. This may be easier to follow along with if the reader has worked with NumPy and GeoPandas before, but I don’t think they have to be too familiar with those packages. 

$\sqrt{3x-1}+(1+x)^2$

$$
  \tilde{y} \sim N(\tilde{\mu}, Q^{-1}),
$$

```math
\sqrt{3x-1}+(1+x)^2
```

test1

$$
\sqrt{3}
$$

test2

$$
E = mc^2
$$

-----

Let’s say we’ve got a set of polygons representing tree canopies and we want to compute each tree’s maximum height and its leaf area index (LAI). For this post I drew a set of polygons over tree canopies around [Soldier’s Circle](https://www.tclf.org/landscapes/soldiers-circle), a parkway in Buffalo, New York. The majority of these trees are American elms. I used Bing imagery as a reference; in the image below you can see the polygons in green overlaying Bing imagery. You can download the polygons in GeoJSON format [here](https://google.com). **TODO: make sure geojson's good, upload it to an "assets" folder or something, and actually link to it**

<div style="text-align: center">
  <img
    src="/assets/tree_polygons_over_bing_imagery.png"
    width="1000"
    height="700"
  />
</div>

Aerial lidar can help us easily compute both of the variables we're interested in. The lidar we’ll work with is from a 2019 survey. Here we’ll just use a single tile. USGS provides this data in laz format, and you can download it [here](https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/LPC/projects/NY_3County_2019_A19/NY_3County_2019/LAZ/USGS_LPC_NY_3County_2019_A19_e1382n2339_2019.laz).

Both the lidar and the polygons use the coordinate reference system specified by EPSG:6350. 

If we were to do something like this in practice, we’d probably be segmenting trees rather than working with polygons. I’ll go over how we can do this in a future post. Also, several of the polygons I drew will include some returns from rootftops. I’ll discuss how to handle that in the future, as well, but for now LAI estimates will be affected.

-----

First, let’s change directories to the folder in which we’re keeping the data, read our polygon layer as a GeoDataFrame, and instantiate the fields we’d like to calculate with the help of the point cloud:
{% highlight Python %}
#change directories and read polygons
import geopandas as gpd
import os

os.chdir(r'D:\HeDemo\data')
tree_polys_filename = 'soldiers_circle_tree_polys_epsg6350.geojson'

tree_polys = gpd.read_file(tree_polys_filename)
tree_polys['max_height'] = np.nan
tree_polys['lai'] = np.nan
{% endhighlight %}
Next, we’ll read our lidar dataset and convert it to an ndarray. laspy makes it very easy to read lidar data, and while laspy’s LasData objects already represent point clouds as mdarrays (see the point.array attribute), I find that it’s often easier to just convert the dimensions of interest directly, since so many libraries are able to work with ndarrays. We want the [scaled versions of the x, y, and z dimensions](https://laspy.readthedocs.io/en/latest/intro.html#point-records), which aren’t included in LasData.point.array, so we’ll make a new ndarray which includes these dimensions, the class code dimension (which we’ll use to find heights above ground elevation), the scan angle (which we’ll use to calculate LAI), and an empty column that we’ll fill in with heights. 
{% highlight Python %}
import laspy
import numpy as np

lidar_filename = r'D:\HeDemo\data\lidar\USGS_LPC_NY_3County_2019_A19_e1382n2339_2019.laz'
scan_angle_factor = 0.006 #to convert from value given by las object, see https://github.com/ASPRSorg/LAS/issues/41#issuecomment-344300998.

las = laspy.read(fpath)
np_las = np.transpose(np.vstack([las.x, las.y, las.z, las.classification, las.scan_angle*scan_angle_factor, np.zeros(len(las))]))
{% endhighlight %}
Note that the scan angle dimension doesn’t actually provide us with angles, and we need to multiply it by a factor to get the actual angle. See [this comment](https://github.com/ASPRSorg/LAS/issues/41#issuecomment-344300998) for an in-depth explanation.
In practice, we might want to make an array from more than one lidar file, in which case something like this could be helpful: 
{% highlight Python %}
np.append(np_las_1, np_las_2, axis=0)
{% endhighlight%}

Next, we’ll make a spatial index. k-D trees are popular spatial indexes for point clouds because they’re good for nearest neighbor queries (what are the n nearest neighbors?) and range queries (what are the neighbors within a distance r from some point?) on point datasets. Here we’ll use [scipy’s implementation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html). Since we’re only concerned with finding which returns have (x,y) coordinates that fall in polygons, we’ll just build the tree for those dimensions. For more on k-D trees, Xiao (2016) gives a great introduction for GIS folks. 
{% highlight Python %}
from scipy.spatial import KDTree
kdtree_las_xy = KDTree(np_las[:, :2])
{% endhighlight %}
**TODO: picture of k-d tree. Take from Xiao, p 81?**

Now we’re ready to calculate our variables of interest. For each polygon:

* Do a range query to find the points with (x,y) coordinates near the polygon. I like to use a range that’s slightly larger than the minimum enclosing circle of the minimum bounding rectangle of the polygon, since it’s guaranteed to include each point in the polygon and also makes it more likely that we’ll include many ground returns. Here, we’ll use a buffer of 2 meters for the minimum enclosing circle.
{% highlight Python %}
#"poly" is a row returned by GeoDataFrame.iterrows()
search_radius_buffer = 2
poly_mbr = np.array(poly['geometry'].bounds) #minx, miny, maxx, maxy
poly_mbr_centroid = np.array([(poly_mbr[2] + poly_mbr[0]) / 2, (poly_mbr[3] + poly_mbr[1]) / 2])
poly_search_radius = search_radius_buffer + ((poly_mbr[2] - poly_mbr_centroid[0])**2 + (poly_mbr[3] - poly_mbr_centroid[1])**2) ** 0.5
near_poly_np_las = np_las[kdtree_las_xy.query_ball_point(poly_mbr_centroid, poly_search_radius)]
{% endhighlight %}
**TODO: picture of buffered MEC, MEC, MBR, and polygon**
* Next, we’ll interpolate ground returns and use the derived surface to compute height above ground for each point satisfying the range query. USGS lidar datasets usually have good ground classifications, so we’ll just work with theirs and not try to come up with our own. We’ll use TIN interpolation, which is fast if we’re only considering a few hundred points at a time like we are here, and will accurate enough for our purposes (separating returns near or on the ground from those which probably correspond to tree canopies). Again, we’ll take advantage of scipy for this.
{% highlight Python %}
def get_height_above_ground_for_points(points):
    '''Estimate height above ground for a (ground-classified) point cloud, using TIN interpolation.
    see https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html'''
    ground_las_class_code = 2
    ground_points = points[points[:,3] == ground_las_class_code]
    ground_interp = LinearNDInterpolator(points = ground_points[:,:2], values = ground_points[:,2])
    ground_elevations = ground_interp(points[:,:2])
    return points[:,2] - ground_elevations

near_poly_np_las[:, 5] = get_height_above_ground_for_points(near_poly_np_las)
{% endhighlight %}
* Now we’re ready to find the points which fall in our polygon. One approach for doing this that seems appealing is the [shapely.contains() method](https://shapely.readthedocs.io/en/latest/reference/shapely.contains.html#shapely.contains), since our polygons are already represented as shapely objects by GeoPandas. However, I’ve found that in practice it takes too long to convert an ndarray to a list of shapely points. What works better is to convert the polygon to a list or ndarray of coordinates, then run some other function for doing point-in-polygon tests. Here I'll use a function mostly lifted from Xiao (2016), but which doesn’t consider points intersecting polygons’ line segments or vertices. This works well enough for our purposes, but if you’d like more speed you can check out [this implementation of the winding number algorithm](https://github.com/sasamil/PointInPolygon_Py/blob/master/pointInside.py) which leverages Numba—as shown by [this stackexchange post](https://stackoverflow.com/a/48760556), Numba can make a huge difference here. For more on point-in-polygon tests, Xiao (2016) is again a good reference.
{% highlight Python %}
def point_in_polygon(point, polygon):
    '''Tests whether (x,y) coordinates of point is in polygon using ray tracing algorithm
    Argument point is numpy array with shape (m,) where m >= 2, polygon is a numpy array with shape (n, 2) where n-1 is the number of vertices in the polygon and polygon[0] == polygon[n]
    Adopted from Xiao (2016), see https://github.com/gisalgs/geom/blob/master/point_in_polygon.py'''
    point_is_in_polygon = False
    for i in range(len(polygon)-1):
        poly_vertex1, poly_vertex2 = polygon[i], polygon[i+1]
        yside1 = poly_vertex1[1] >= point[1]
        yside2 = poly_vertex2[1] >= point[1]
        xside1 = poly_vertex1[0] >= point[0]
        xside2 = poly_vertex2[0] >= point[0]
        if yside1 != yside2: #if point's y coord is between y coords of line segment
            if xside1 and xside2: #if point is to the left of both polygon vertices:
                point_is_in_polygon = not point_is_in_polygon
            else:
                m = poly_vertex2[0] - (poly_vertex2[1]-point[1])*(poly_vertex1[0]-poly_vertex2[0])/(poly_vertex1[1]-poly_vertex2[1]) #intersection point
                if m > point[0]:
                    point_is_in_polygon = not point_is_in_polygon
    return point_is_in_polygon
    
poly_vertex_array = np.array(poly['geometry'].exterior.coords)
in_poly_np_las = near_poly_np_las[[point_in_polygon(las_point, poly_vertex_array) for las_point in near_poly_np_las]]
{% endhighlight %}
**TODO?: picture, if we have time**
* From here, we can easily find the maximum height above ground.
{% highlight Python %}
#tree_polys.at[index, 'max_height'] = np.max(in_poly_np_las[~np.isnan(near_poly_np_las[:, 5]), 5]) #see text below for why we use ~np.isnan()
tree_polys.at[index, 'max_height'] = np.max(in_poly_np_las[near_poly_np_las[:, 5], 5]) #see text below for why we use ~np.isnan()
{% endhighlight %}
* LAI can also be easily calculated now. Often the distribution of light throughout canopies is described in a method similar to the Beer-Lambert law for light attenuation through a homogenous medium (Jones, 2013). More specifically, this looks something like:
$$
E = mc^2
$$
Where I is the irradiance (i.e., radiant flux) at the ground surface [W m^-2], I_0 is the irradiance at the top of the canopy, L is the leaf area index, and k is an “extinction coefficient” representing the ratio of the area of shadows cast by leaves to the actual area of the leaves (Jones, 2013). Solving for LAI gives us
L = -1/k * ln(I / I_0).
For a spherical leaf angle distribution, meaning all leaves have a uniform probability for all angles from the horizontal plane (zenith angles %theta%) %in% [0°, 90°], we have k = 0.5. Following Richardson et. al. (2009), we’ll substitute I with the total number of ground returns R_g and I_0 with the total number of returns in the polygon R_t, and also model the effects of lidar scanning angle using Lambert’s cosine law I = I_0*cos(%theta%). This gives us the model:
L = -cos(%mean_theta_lidar%) / 0.5 * ln(R_g/R_t)
Where %mean_theta_lidar% is the mean lidar scanning angle of all returns in the polygon. Here’s the code for implementing this model:
{% highlight Python %}
lambert_beer_extinction_coefficient_when_scan_angle_is_0 = 0.5 #assumes spherical leaf angle distribution.
ground_elev_threshold = 0.05 #in m

total_number_of_returns_in_poly = len(in_poly_np_las)
number_of_ground_returns_in_poly = np.sum(in_poly_np_las[:,5] <= ground_elev_threshold)
mean_lidar_scanning_angle = np.mean(in_poly_np_las[:,4])
tree_polys.at[index, 'lai'] = -np.cos(mean_lidar_scanning_angle) / lambert_beer_extinction_coefficient_when_scan_angle_is_0 * np.log(number_of_ground_returns_in_poly / total_number_of_returns_in_poly)
{% endhighlight %}
Note we consider a return to be a ground return if its height above ground is less than or equal to 5 cm.

**TODO: Explain results (run time, etc), make some comments, show some pictures.**

If you’d like to reproduce exactly what we did here, here’s a [script](google.com) which puts together everything we did above. **TODO: actually link to script** I tried to keep things in the order we discuss them in this post, in practice it would definitely make more sense to import our libraries at the start of the script, etc. Again, the polygons can be found here [link] and the laz file here [link]. 
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

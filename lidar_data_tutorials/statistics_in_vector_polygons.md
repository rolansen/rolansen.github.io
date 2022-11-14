**This is the first of a series of tutorials I'm writing on using Python to do interesting things that mix lidar data and other types of GIS datasets. By "lidar data" I mean discrete return aerial lidar point clouds.**

**In this post, I describe a use case demonstrating how to identify lidar returns which have some spatial relationship to vector (polygon) data. Once we do this, it's easy to compute any sort of statistic or metric, then assign the values to some new field in the vector dataset.**

-----

Free, high-quality lidar datasets have become much more common in the last few years than they were a decade ago. A good portion of the planet now has at least one representative lidar survey. In the future, surveys will become more widespread and regular, and quality will keep improving.

This increase in the availability of lidar data has caused it to become more commonly utilized in the GIS world, a trend which has been bolstered by the great tools that are available to the community nowadays. LAStools and PDAL are specialized software suites that are now often used by GIS professionals, WhiteboxTools has great support for lidar data, and ArcGIS and QGIS have both been adding more options for visualizing and processing lidar point clouds. 

<div style="text-align: center">
  <figure>
      <img
       src="https://www.lutraconsulting.co.uk/img/posts/pointcloud_copc_melbourne.gif"
       width="600"
       height="397"
     />
     <figcaption>Cloud Optimized Point Cloud rendered by QGIS. Image by Lutra Consulting.</figcaption>
  </figure>
</div>\

However, at the moment tools for integrating lidar data with other GIS data types (vector & raster) can still feel somewhat limited. What if we wanted to do something like identify the returns in or around features from a polygon layer, then compute a new field for this layer based on these returns? What if we wanted to limit our data to returns which fall in certain classes of a land cover image? 

It seems relatively difficult to do these things with the software I listed above. The R package [lidR](https://cran.r-project.org/web/packages/lidR/index.html) is one helpful option. Another is Python, which is perhaps more ubiquitous at this point, and has the huge advantage of letting us work with point clouds in the form of NumPy arrays. Here, I’ll go over an example which will show how easily, flexibly, and efficiently we can work with lidar point clouds using tools from the Python ecosystem. 

To reproduce the code on this page, you'll need a Python environment which includes laspy, the Python bindings for the Rust crate laz-rs (see the [laspy installation instructions](https://laspy.readthedocs.io/en/latest/installation.html)), NumPy, SciPy, and GeoPandas. Things may be easier to follow along with if the reader has worked with NumPy and GeoPandas before, but I don’t think they have to be too familiar with those packages. 

-----

Let’s consider a set of polygons representing tree canopies, and say that we want to compute each tree’s maximum height and its leaf area index (LAI). For this post I drew a set of polygons over tree canopies around [Soldier’s Circle](https://www.tclf.org/landscapes/soldiers-circle), a parkway in Buffalo, New York. The majority of these trees are American elms. I used Bing imagery as a reference. You can download the polygons in GeoJSON format [here](https://github.com/rolansen/rolansen.github.io/blob/main/data/canopy_polygons.geojson).

<div style="text-align: center">
  <figure>
      <img
       src="/assets/tree_polygons_over_bing_imagery.png"
       width="500"
       height="350"
     />
     <figcaption>Polygons in green overlaying Bing imagery</figcaption>
  </figure>
</div>\

Aerial lidar can help us easily compute both of the variables we're interested in. We’ll work with data from a single tile, collected during a 2019 survey. USGS provides this dataset in laz format, and you can download it [here.](https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/LPC/projects/NY_3County_2019_A19/NY_3County_2019/LAZ/USGS_LPC_NY_3County_2019_A19_e1382n2339_2019.laz)

Both the lidar and the polygons use the coordinate reference system specified by EPSG:6350. The code below assumes both datasets are in the same directory.

If we were to do something like this in practice, we’d probably be segmenting trees rather than working with polygons. I’ll go over how we can do this in a future post. Also, several of the polygons I drew will include some returns from rootftops. I’ll discuss how to handle that in the future, as well, but for now a few LAI estimates will be affected.

-----

First, let’s change directories to the folder in which we’re keeping the data, read our polygon layer as a GeoDataFrame, and instantiate fields for the we’re interested in:
{% highlight Python %}
import geopandas as gpd
import numpy as np
import os

os.chdir(r'your\directory')
tree_polys_filename = 'canopy_polygons.geojson'

tree_polys = gpd.read_file(tree_polys_filename)
tree_polys['max_height'] = np.nan
tree_polys['lai'] = np.nan
{% endhighlight %}
Next, we’ll read our lidar dataset and convert it to an ndarray. laspy makes it very easy to read lidar data, and while laspy’s LasData objects already represent point clouds as mdarrays (see the point.array attribute), I find that it’s often easier to just convert the dimensions of interest directly, since so many libraries are able to work with ndarrays. We need the [scaled versions of the x, y, and z dimensions](https://laspy.readthedocs.io/en/latest/intro.html#point-records), which aren’t included in LasData.point.array, so we’ll make a new ndarray which includes these dimensions, the class code dimension (which we’ll use to find heights above ground elevation), the scan angle (which we’ll use to calculate LAI), and an empty column that we’ll fill in with heights. 
{% highlight Python %}
import laspy

lidar_filename = 'USGS_LPC_NY_3County_2019_A19_e1382n2339_2019.laz'
scan_angle_factor = 0.006 #to convert from value given by las object, see https://github.com/ASPRSorg/LAS/issues/41#issuecomment-344300998.

las = laspy.read(lidar_filename)
np_las = np.transpose(np.vstack([las.x, las.y, las.z, las.classification, las.scan_angle*scan_angle_factor, np.zeros(len(las))]))
{% endhighlight %}
Note that we need to multiply the scan angle dimension by a factor to get actual scan angles. See [this comment](https://github.com/ASPRSorg/LAS/issues/41#issuecomment-344300998) for an in-depth explanation.
In practice, we might want to make an array from more than one lidar file, in which case something like this could be helpful: 
{% highlight Python %}
np.append(np_las_1, np_las_2, axis=0)
{% endhighlight%}

Next, we’ll make a spatial index. *k*-D trees are popular spatial indexes for point clouds because they can efficiently answer nearest neighbor queries (what are the *n* nearest neighbors?) and range queries (what are the neighbors within a distance *r* from some point?) on point datasets. Let's use [SciPy’s implementation.](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html) Since we’re only concerned with finding which returns have (x,y) coordinates that fall in polygons, we’ll build the tree for just those dimensions. Xiao (2015) gives a nice introduction to *k*-D trees that's geared towards GIS folks, if you'd like to learn more. 
{% highlight Python %}
from scipy.spatial import KDTree
kdtree_las_xy = KDTree(np_las[:, :2])
{% endhighlight %}

<div style="text-align: center">
  <figure>
      <img
       src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/bf/Kdtree_2d.svg/1024px-Kdtree_2d.svg.png"
       width="400"
       height="400"
     />
     <figcaption>Depiction of k-d tree for a two-dimensional dataset. From Wikipedia.</figcaption>
  </figure>
</div>\

Now we’re ready to calculate the variables of interest. For each polygon:

* Do a range query to find the points with (x,y) coordinates near the polygon. I like to use the centroid of the polygon's minimum bounding rectangle as the reference point, and a range that’s slightly larger than the minimum enclosing circle of the rectangle. This way the query's guaranteed to include each point in the polygon and also makes it less likely that we won't exclude many ground returns. Below we use a buffer of 2 meters for the minimum enclosing circle.
{% highlight Python %}
#"poly" is a row returned by GeoDataFrame.iterrows()
search_radius_buffer = 2
poly_mbr = np.array(poly['geometry'].bounds) #minx, miny, maxx, maxy
poly_mbr_centroid = np.array([(poly_mbr[2] + poly_mbr[0]) / 2, (poly_mbr[3] + poly_mbr[1]) / 2])
poly_search_radius = search_radius_buffer + ((poly_mbr[2] - poly_mbr_centroid[0])**2 + (poly_mbr[3] - poly_mbr_centroid[1])**2) ** 0.5
near_poly_np_las = np_las[kdtree_las_xy.query_ball_point(poly_mbr_centroid, poly_search_radius)]
{% endhighlight %}
<div style="text-align: center">
  <figure>
      <img
       src="/assets/bounding_geometries.png"
       width="500"
       height="350"
     />
     <figcaption>Polygons (green), their envelopes (blue), the envelopes' minimum enclosing circles (red), and buffered circles (pink).</figcaption>
  </figure>
</div>\

* Next, we’ll apply TIN interpolation to ground returns and use the derived surface to compute height above ground for each point satisfying the range query. USGS lidar datasets usually have good ground classifications, so we’ll just work with theirs. TIN interpolation is fast if we’re only considering a few hundred points at a time, like we are here, and should be accurate enough for our purposes. Again, we’ll take advantage of SciPy, this time using their [LinearNDInterpolator.](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html)
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
* Now we’re ready to find the points which fall in our polygon. One approach that seems appealing is the [shapely.contains() method](https://shapely.readthedocs.io/en/latest/reference/shapely.contains.html#shapely.contains), since our polygons are already represented as Shapely objects by GeoPandas. However, I’ve found that in practice it takes too long to convert an ndarray to a list of Shapely points. What works better is to convert the polygon to a list or ndarray of coordinates, then run some other function for doing point-in-polygon tests. Here I'll use a function mostly lifted from Xiao (2015), but which for simplicity's sake doesn’t consider points intersecting polygons’ line segments or vertices. This works well enough, but if you’d like more generality and speed you can check out [this implementation of the winding number algorithm which leverages Numba](https://github.com/sasamil/PointInPolygon_Py/blob/master/pointInside.py)—as shown by [this stackexchange post](https://stackoverflow.com/a/48760556), Numba can make a huge difference here. For more on point-in-polygon tests, Xiao (2015) is again a good reference.
{% highlight Python %}
def point_in_polygon(point, polygon):
    '''Tests whether (x,y) coordinates of point is in polygon using ray tracing algorithm
    Argument point is numpy array with shape (m,) where m >= 2, polygon is a numpy array with shape (n, 2) where n-1 is the number of vertices in the polygon and polygon[0] == polygon[n]
    Adopted from Xiao (2015), see https://github.com/gisalgs/geom/blob/master/point_in_polygon.py'''
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
* From here, we can easily find the maximum height above ground.
{% highlight Python %}
#"index" is a row index from GeoDataFrame.iterrows()
tree_polys.at[index, 'max_height'] = np.max(in_poly_np_las[in_poly_np_las[:, 5], 5])
{% endhighlight %}
* LAI can also be easily calculated at this point. Often the distribution of light throughout canopies is described in a method similar to the Beer-Lambert law for light attenuation through a homogenous medium (Jones, 2013). More specifically, this looks something like:

*I* = *I*₀ exp(-*kL*),

Where *I* is the irradiance (i.e., radiant flux) at the ground surface [W m^-2], *I*₀ is the irradiance at the top of the canopy, *L* is LAI, and *k* is an “extinction coefficient” representing the ratio of the area of shadows cast by leaves to the actual area of the leaves (Jones, 2013). Solving for LAI gives 

*L* = -ln(*I* / *I*₀) / *k*.

For a spherical leaf angle distribution, meaning all leaves have a uniform probability for any zenith angle *&theta;* &isin; [0°, 90°], we have *k* = 0.5. Following Richardson et. al. (2009), we’ll substitute *I* with the total number of returns reaching the ground surface *R*ₛ and *I*₀ with the total number of returns in the polygon *R*ₜ, and also account for the effects of lidar scanning angle using Lambert’s cosine law *I* = *I*₀cos(*&theta;*). This gives us the model:

L = -cos(&Theta;) ln(*R*ₛ/*R*ₜ) / 0.5,

Where &Theta; is the mean lidar scanning angle of all returns in the polygon. Here’s some code which implements this model:
{% highlight Python %}
lambert_beer_extinction_coefficient_when_scan_angle_is_0 = 0.5 #assumes spherical leaf angle distribution.
ground_elev_threshold = 0.05 #in m

total_number_of_returns_in_poly = len(in_poly_np_las)
number_of_ground_returns_in_poly = np.sum(in_poly_np_las[:,5] <= ground_elev_threshold)
mean_lidar_scanning_angle = np.deg2rad(np.mean(in_poly_np_las[:,4]))
tree_polys.at[index, 'lai'] = -np.cos(mean_lidar_scanning_angle) / lambert_beer_extinction_coefficient_when_scan_angle_is_0 * np.log(number_of_ground_returns_in_poly / total_number_of_returns_in_poly)
{% endhighlight %}
Note we consider a return to be a ground return if its height above ground is less than or equal to 5 cm.

When we're finished, we can write the polygons with their new fields to a file:
{% highlight Python %}
output_path = 'canopy_polygons_with_height_and_lai.geojson'
tree_polys.to_file(output_path)
{% endhighlight %}
-----

On my mid-end home desktop this all takes about a minute. We would have seen some performance gains from multithreading, since NumPy operations aren't constrained by Python's global interpreter lock. The same can be said of the SciPy methods used here. Multithreading would have been especially advantageous if we had to read multiple laz files. As mentioned before, we could have made a more efficient point-in-polygon function, too.

If you’d like to reproduce exactly what we did here, here’s a [script](https://github.com/rolansen/rolansen.github.io/blob/main/code/find_height_and_lai_in_canopy_polys.py) which puts together everything we did above. I tried to keep things in the order we discuss them in this post, but in practice it would definitely make more sense to import our libraries at the start of the script, etc. I also put in some logic to deal with cases where polygons don't intersect any lidar returns. Again, the polygons can be found [here](https://github.com/rolansen/rolansen.github.io/blob/main/data/canopy_polygons.geojson) and the laz file [here](https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/LPC/projects/NY_3County_2019_A19/NY_3County_2019/LAZ/USGS_LPC_NY_3County_2019_A19_e1382n2339_2019.laz). 

-----

In future posts, I'll show how rasterio can help you can do some cool stuff with lidar and raster data, and how we can segment individual trees and buildings.

-----

References

Jones, H. G. (2013). *Plants and microclimate: a quantitative approach to environmental plant physiology*. Cambridge university press.

Richardson, J. J., Moskal, L. M., & Kim, S. H. (2009). Modeling approaches to estimate effective leaf area index from aerial discrete-return LIDAR. *Agricultural and Forest Meteorology, 149*(6-7), 1152-1160.

Xiao, N. (2015). *GIS algorithms*. Sage.

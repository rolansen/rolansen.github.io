**This is the second post in a series on using the geospatial Python ecosystem to analyze spatial relationships between lidar point clouds and other types of geospatial datasets. This will be the first of maybe two or three on working with raster data. This post in particular will discuss how to query [Entwine Point Tile (EPT) datasets](https://entwine.io/en/latest/entwine-point-tile.html) and efficiently create raster data from the results.**

-----

Today desktop GIS software generally makes it easy to create digital elevation/surface models (DEM’s/DSM’s) from lidar data. While that’s actually more or less what we’re going to do here, doing this with Python gives us some obvious advantages, perhaps the most important being flexibility.

Also, the existence of a [small library written by the developers of the EPT format](https://github.com/hobu/ept-python) for easily querying EPT services and returning [laspy](https://laspy.readthedocs.io/en/latest/index.html) objects, as well as the current lack of EPT support in most desktop GIS software (with the notable exception of [QGIS](https://docs.qgis.org/testing/en/docs/user_manual/working_with_point_clouds/point_clouds.html)) makes Python an ideal way to take advantage of the massive amount of [EPT data USGS is hosting](https://usgs.entwine.io/).

In addition to the ept-python library mentioned above, to run the code shown in this post you'll need rasterio, geopandas (optionally pyogrio), shapely, nest_asyncio, and some of the libraries typically included with scientific Python distributions: numpy, scipy, and matplotlib. We'll also be using asyncio syntax that will only work for Python >= 3.7.

For convenience, here’s the full list of imports:

{% highlight Python %}
import rasterio
from rasterio.windows import Window
from rasterio.io import MemoryFile
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.sample import sample_gen
from rasterio.mask import mask
from pyogrio import read_dataframe #delete this line if pyogrio is not installed, and instead use geopandas
from shapely.geometry import Polygon
import ept
from scipy.ndimage import label, convolve
from scipy.interpolate import griddata
from scipy.stats import f_oneway
import numpy as np
import matplotlib.pyplot as plt
import nest_asyncio
from aiohttp import ClientSession
import asyncio
from urllib.request import urlopen
import os
{% endhighlight %}

-----

Our primary goal will be to create a digital height model (DHM) representing height above ground elevation. Before downloading the lidar we’ll use to do this, we need to get the following::
* the boundary of the lidar survey of interest
* a "template" raster with a known CRS, spatial resolution, and extent that we'll sample the DHM to
* a grid that we'll use to make requests to the EPT resource

We'll start by extracting the survey boundary from USGS's [WESM layer](https://www.usgs.gov/ngp-standards-and-specifications/wesm-data-dictionary), a source of spatial metadata for the lidar datasets it distributes. Many of these are [available in EPT format](https://usgs.entwine.io). We'll be working with the [work unit](https://www.usgs.gov/ngp-standards-and-specifications/wesm-data-dictionary-general-attributes#workunit) "NC_Phase4_Rowan_2017." I picked this one because it's relatively small and includes a diverse mix of cultivated, impervious, and forested land cover (see the second plot below).

Below we use pyogrio to read the boundary to a single-rowed geopandas GeoDataFrame. This is fast and doesn't require us to write much code, but if you'd prefer you could do this [directly with geopandas if you have a recent version of Fiona](https://geopandas.org/en/stable/docs/user_guide/io.html#sql-where-filter), read the whole WESM layer and then subset it, or since we'll only be working with the work unit's geometry from here on out read the row of interest with Fiona and convert the geometry to a Shapely polygon with shapely.wkb.loads().

{% highlight Python %}
srid = 3857 #EPT service uses EPSG:3857, and we'll project workunit boundary to and download NLCD with this CRS
workunit_boundary_not_proj = read_dataframe('path\\to\\WESM.gpkg', layer='WESM', where="workunit = 'NC_Phase4_Rowan_2017'") #this reads to GeoDataFrame, use geopandas instead if you don't have pyogrio
workunit_boundary = workunit_boundary_not_proj.to_crs(srid)
{% endhighlight %}

<div style="text-align: center">
  <figure>
      <img
       src="/assets/workunit_boundary.png"
       width="482"
       height="413"
     />
     <figcaption>Boundary for work unit "NC_Phase4_Rowan_2017." CRS here and for plots below is EPSG:3857</figcaption>
  </figure>
</div>\

Now we'll use this boundary to define the template raster. We'll be using rasterio to read and write our raster data, allowing us to leverage numpy and scipy for most of our computations.

We could just manually define our DHM's extent, CRS, and spatial resolution. However, often in practice we'll want to compare elevation models to other raster datasets, in which case it'll make sense to sample the former to the spatial resolution, etc. of the latter. For this post I'd like to compare the DHM with [National Land Cover Database (NLCD) data](https://www.mrlc.gov/), so I can look a bit into any relationships that might exist between height and land cover (LC) class. We'll work with NLCD data for 2016, since that's the closest available year to the vintage of the lidar work unit (2017). If you'd like, [you could download the ERDAS .img file from the Multi-Resolution Land Characteristics Consortium](https://www.mrlc.gov/data/nlcd-2016-land-cover-conus), which will take up about 30 GB after being unzipped. I’ll make a request from [their WMS site](https://www.mrlc.gov/geoserver/mrlc_display/wms?service=WMS&request=GetCapabilities), though, so we can get into a quick tangent about reading WMS data to rasterio.

Either way, we’ll want to limit the data we read to the work unit’s bounding box. Below I use the geometry of the work unit to set the parameters we’ll use to make the WMS request–-namely, the width, height, and bounding box. The CRS will be EPSG:3857, which we projected the work unit geometry to. The raster's width and height are set so its spatial resolution is 30 meters and it covers the work unit’s bounding box. Note the WMS service will resample the image, so it won't exactly match the results of cropping and reprojecting the .img file ourselves. This is because we aren’t considering the CRS’s origin. After setting the parameters, we request a GeoTIFF and just read it as a bytes object.

{% highlight Python %}
#set parameters
#we increase bounding box and calculate image width/height to get the smallest raster with 30 m pixels which contains the workunit's bbox
pixel_lindim_length = 30 #in meters
workunit_bbox = np.asarray(workunit_boundary.bounds)[0] #xmin, ymin, xmax, ymax
workunit_width = workunit_bbox[2] - workunit_bbox[0]
workunit_height = workunit_bbox[3] - workunit_bbox[1]
img_bbox = [workunit_bbox[0], workunit_bbox[1], \
                     workunit_bbox[2] + pixel_lindim_length - workunit_width % pixel_lindim_length, \
                     workunit_bbox[3] + pixel_lindim_length - workunit_height % pixel_lindim_length] #this variable and others here appear later
img_width = int(workunit_width // pixel_lindim_length) + 1 #workunit_width is np.float64, need to cast
img_height = int(workunit_height // pixel_lindim_length) + 1
img_bbox_str = ','.join([str(_) for _ in img_bbox])

#request NLCD 2016 from WMS
nlcd_url = f'https://www.mrlc.gov/geoserver/mrlc_display/wms/v2?service=WMS&version=1.3.0&request=GetMap&layers=NLCD_2016_Land_Cover_L48&format=image/geotiff&width={img_width}&height={img_height}&bbox={img_bbox_str}&styles=&srs=EPSG:{srid}'
nlcd_bytes = urlopen(nlcd_url).read()
{% endhighlight %}

Next, we’ll take this bytes object, open it as a [rasterio MemoryFile](https://rasterio.readthedocs.io/en/stable/topics/memory-files.html), and save its cell values (a 2-D numpy array) and [Affine transformation matrix](https://pypi.org/project/affine/) as variables for later use. [This stackexchange post explains how to do this smoothly](https://gis.stackexchange.com/a/363923), and we mostly copy the first workflow it discusses. However, we care about the values assigned to cells--ideally they'd match NLCD's [Level II LC class codes](https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description). Since we’re using a GetMap request rather than GetFeatureInfo requests, the ndarray won't use these class codes for cell values. To get them from the values returned by the GetMap request, we'll use a hard-coded dictionary I wrote by comparing [the WMS layer’s legend](https://www.mrlc.gov/geoserver/mrlc_display/ows?service=WMS&request=GetLegendGraphic&format=image%2Fpng&width=20&height=20&layer=NLCD_2016_Land_Cover_L48) with [MLRC's NLCD legend](https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description). Afterwards, we convert the majority of Level II class codes to Level I class codes, mostly just to keep patches from being too small. However, to avoid massive, inconvenient developed patches, class codes 21 are 22 are consolidated into one class while 23 and 24 are consolidated into another. Note also that we mask the image using the work unit’s geometry.

{% highlight Python %}
nlcd_dn_dict = {0:np.nan, 1:11, 2:12, 3:21, 4:22, 5:23, 6:24, 7:31, 8:32, 9:41, 10:42, 11:43, 12:51, 13:52, 14:71, 15:72, 16:73, 17:74, 18:81, 19:82, 20:90, 21:95} #0 comes from rasterio.mask.mask()
with rasterio.Env(), MemoryFile(nlcd_bytes) as memfile:
    with memfile.open() as nlcd_not_masked:
        nlcd, nlcd_transform = mask(nlcd_not_masked, workunit_boundary['geometry'])
        nlcd = nlcd[0] #reshape from (1, row, col) to (row, col)
        nlcd_profile = nlcd_not_masked.profile
        nlcd_ndarray_level_2 = np.vectorize(nlcd_dn_dict.get)(nlcd)
        nlcd_ndarray_level_1 = nlcd_ndarray_level_2 // 10 #2/22/23: use level 1 for non-urban, modify urban a bit
        nlcd_ndarray = np.where(nlcd_ndarray_level_1 != 2, nlcd_ndarray_level_1, 21 + 2*(nlcd_ndarray_level_2 > 22)) #assign classes 21 & 22 to value 21 and 23 & 24 to 22
{% endhighlight %}

<div style="text-align: center">
  <figure>
      <img
       src="/assets/nlcd_consolidated.png"
       width="557"
       height="300"
     />
     <figcaption>NLCD 2016 within workunit boundary, consolidated as described above.</figcaption>
  </figure>
</div>\

When we start downloading lidar, we could make one request per pixel, which would keep our code simple since we wouldn't have to consider lidar returns' (x, y) coordinates after making requests. This would be very inefficient, though-–at our desired EPT resolution of 2 (see below) it takes my computer about 5 seconds to get a laspy object for a single 30-by-30 meter pixel. Making one request for the entire point cloud at this resolution also seems impractical unless one has a lot of RAM to work with.

For these reasons, we’ll make requests using a grid defined over the extent of the raster. The code below implements this grid as a list of ept.Bounds objects ([see ept-python's GitHub page](https://github.com/hobu/ept-python)). Grid cell boundaries are aligned with pixel boundaries, and cells which don’t intersect the work unit geometry are discarded.

{% highlight Python %}
grid_lindim_size = 66 #in pixels 
grid_cols = img_width // grid_lindim_size + 1
grid_cols_remainder = img_width % grid_lindim_size
grid_rows = img_height // grid_lindim_size + 1
grid_rows_remainder = img_height % grid_lindim_size

grid_row_length = img_bbox[2] - img_bbox[0]

cell_x_length = grid_lindim_size * pixel_lindim_length
cell_y_length = grid_lindim_size * pixel_lindim_length
cell_ll = np.array([img_bbox[0], img_bbox[1]])
cell_lr = cell_ll + [cell_x_length, 0]
cell_ur = cell_ll + [cell_x_length, cell_y_length]
cell_ul = cell_ll + [0, cell_y_length]
cell = np.array([cell_ll, cell_lr, cell_ur, cell_ul])

#define grid as list of ept.Bounds objects which intersect the workunit's boundary
grid = []
for row in range(grid_rows):
    for col in range(grid_cols):
        #append cell to list if it intersects workunit boundary
        if workunit_boundary['geometry'].intersects(Polygon(cell)).any():
            cell_ept_bounds = ept.Bounds(cell[0,0], cell[0,1], -np.inf,
                                         cell[2,0], cell[2,1], np.inf) #xmin, ymin, zmin, xmax, ymax, zmax
            grid.append(cell_ept_bounds)

        #shift cell along x-axis
        if col != grid_cols - 1:
            cell += [cell_x_length, 0]
        else:
            cell += [grid_cols_remainder * pixel_lindim_length, 0] #avoids making EPT requests outside workunit boundaries

    #shift cell back to min x
    cell -= [grid_row_length, 0]

    #shift cell along y-axis
    if row != grid_rows - 1:
        cell += [0, cell_y_length]
    else:
        cell += [0, grid_rows_remainder * pixel_lindim_length] #avoids making EPT requests outside workunit boundaries
{% endhighlight %}

The chosen grid resolution of 66 pixels (1980 meters) is pretty arbitrary. One could increase it if they'd like to speed up DHM creation, or decrease it if they'd rather use fewer machine resources.

<div style="text-align: center">
  <figure>
      <img
       src="/assets/classified_grid.png"
       width="557"
       height="300"
     />
     <figcaption>Grid defined over subsetted NLCD dimensions. Cell width & height are 990 meters. Cells are symbolized with the majority value of consolidated NLCD values within them.</figcaption>
  </figure>
</div>\

-----

Now we're almost ready to work with the lidar data. All that's left to do is set a couple of parameters and a few other small, miscellaneous things.

The first will involve asyncio. As shown below, lidar data will be requested and processed asynchronously. 
Since ept_python will be internally running its own event loop, making our own would throw an error. 
The [nest_asyncio library](https://pypi.org/project/nest-asyncio/) conveniently allows us to make our own event loop, which we'll use to 
run coroutines that in turn call coroutines from ept_python--all we need to do is import the library and run *nest_asyncio.apply()*.
We also make a [semaphore](https://docs.python.org/3/library/asyncio-sync.html#asyncio.Semaphore), since if we don't we'll run out of stack space while making requests. 
I use a value of 15, but similar to the granularity of the grid this is somewhat arbitrary and can likely be adjusted for your hardware.

{% highlight Python %}
nest_asyncio.apply() #so ept event loop doesn't interfere with our own
sem_value = 15 #arbitrary, could be adjusted to account for machine's capabilities
sem = asyncio.Semaphore(sem_value) 
{% endhighlight %}

We also need to tell ept_python where to make requests and at what resolution. Below *ept_site* points to the URL of the lidar survey's corresponding EPT dataset. 
This is done because, unfortunately, WESM doesn't have an attribute for EPT URL's, although it's probably possible to construct them from other attributes.
As discussed [here](https://pdal.io/en/2.4.3/stages/readers.ept.html), the resolution value for an EPT request roughly specifies a cutoff for octree leaf sizes, with at most one lidar return given per leaf. 
Once again, the chosen value here is a bit subjective.

{% highlight Python %}
ept_site = 'https://s3-us-west-2.amazonaws.com/usgs-lidar-public/USGS_LPC_NC_Phase4_Rowan_2017_LAS_2019' #unfortunately is not given in WESM
ept_resolution = 2 #compromise between speed, accuracy, and reliably having enough ground returns in NLCD pixels to calculate ground elevation
{% endhighlight %}

One of the key steps for assigning values to our raster data will be converting coordinates for the projection (EPSG:3857) to image coordinates. 
When we made the NLCD raster we also got an affine transformation matrix, for converting image coordinates to projection coordinates, as an Affine object. 
We could use this and *rasterio.transform.rowcol()* to get the image coordinates corresponding to our lidar returns' (x, y) coordinates, but in practice that's slow.
A faster option is to use the ~ operator to invert the affine transformation matrix, then represent the matrix with an ndarray.
That way we can just multiply an ndarray containing points' (x, y) coordinates by this transformation matrix, as will be shown below. 

{% highlight Python %}
inv_transform = ~nlcd_transform
affine_mat = np.transpose([[inv_transform.a, inv_transform.b, inv_transform.c],
                           [inv_transform.d, inv_transform.e, inv_transform.f],
                           [0, 0, 1]]) #transposing so results of transform have shape (number of returns, 3), or (number of returns, 2) after removing column of 1's
{% endhighlight %}

Finally, let's initialize the ndarray representing our elevation data. To make the DHM, we'll create a DEM from lidar ground returns and a DSM from first returns, then subtract the former from the latter.
As mentioned, the elevation rasters will be sampled to the resolution, etc. of the NLCD raster. Each pixel value will represent lidar returns' average z coordinate values.
For no reason in particular, I'll keep the DSM and DEM in the same ndarray.

{% highlight Python %}
dsm_and_dem = np.full(shape = (img_height, img_width, 2), fill_value = np.nan, dtype=np.float32) #DSM at [:,:,0] and DEM at [:,:,1]
{% endhighlight %}

The time it takes for each tile to have its DEM/DSM values filled in is largely limited
by how long it takes to download requested lidar data.
To compensate, we'll use a set of coroutines to do this work.
The main coroutine, shown below, awaits two others.
It may seem odd to do this instead of just awaiting one coroutine, but I've found that in practice this is faster.

{% highlight Python %}
async def set_elev_values_main(grid):
    epts = await asyncio.gather(*(make_ept_query(req_bounds) for req_bounds in grid)) #can cause stack overflow if semaphore value too great
    await asyncio.gather(*(download_las_and_assign_elev_values(ept_query) for ept_query in epts))
{% endhighlight %}

The first coroutine ran by *set_elev_values_main()*, *make_ept_query()*, will make an EPT request for the given EPT dataset URL, resolution, and tile bounds.  

{% highlight Python %}
async def make_ept_query(req_bounds):
    async with sem:
        return ept.EPT(ept_site, bounds=req_bounds, queryResolution=ept_resolution)
{% endhighlight %}

Next, the requested point cloud data is downloaded and used to set elevation raster values.

{% highlight Python %}
async def download_las_and_assign_elev_values(grid):
    las = await download_las(eqt_query)
    await find_heights_for_tile(las)
{% endhighlight %}

*download_las()* simply runs a function from ept_python, *ept.EPT.as_laspy()*, to create a laspy object using the request.

{% highlight Python %}
async def download_las(ept_query):
    async with sem:
        return ept_query.as_laspy()
{% endhighlight %}

*find_heights_for_tile()* finds image coordinates for each lidar return by using the transformation matrix made prior to running *set_elev_values_main()*. It then makes Boolean indexes for the first returns and ground returns in the tile, and runs a final coroutine for each index.

{% highlight Python %}
async def find_heights_for_tile(las):
    #return if las has no records (in case of imprecision with workunit boundary, etc)
    if las is None or len(las) == 0: 
        return
    
    #find image coordinates
    lidar_rowcols = np.floor(np.column_stack((las.x, las.y, np.ones(len(las.z)))) @ affine_mat).astype(np.int32)[:,(1,0)] #for convenience coordinates are set in order (image row, image column)
    
    #find indexes for first returns and ground returns
    first_returns_idx = las.return_number == 1
    ground_returns_idx = las.classification == 2
    
    #set raster cell values
    await assign_lidar_z_means_within_pixels(lidar_rowcols, las, first_returns_idx, 0)
    await assign_lidar_z_means_within_pixels(lidar_rowcols, las, ground_returns_idx, 1)
{% endhighlight %}

*find_heights_for_tile()* is more intricate than the others. z-coordinates for the lidar data are sorted by image coordinates, the unique image coordinates are identified, a [local reduce operation](https://numpy.org/doc/stable/reference/generated/numpy.ufunc.reduceat.html) is used to find the mean z values corresponding to each image coordinate. In order to efficiently find mean z coordinates for each of the identified image coordinates, and these means are set as values for the appropriate cells of one of the elevation rasters.

{% highlight Python %}
async def assign_lidar_z_means_within_pixels(lidar_rowcols, las, las_idx, image_idx):
    #convert image coordinates to 1-dimensional representation, so we can sort lidar
    lidar_rowcols_at_idx = lidar_rowcols[las_idx]
    rowcol_ids = lidar_rowcols_at_idx[:,1] + lidar_rowcols_at_idx[:,0] * img_width
    
    #sort returns by 1-dimensional image pixel coordinates. need to do this to reduce since np.unique() sorts values
    rowcol_ids_argsort = rowcol_ids.argsort()
    rowcol_ids = rowcol_ids[rowcol_ids_argsort]
    las_z_of_interest = las.z[las_idx][rowcol_ids_argsort]
    
    #identify unique 1D and 2D image pixel coordinates, as well as variables we'll use to find mean
    unique_rowcol_ids, unique_rowcol_ids_idx, counts_in_pixels = np.unique(rowcol_ids, return_index=True, return_counts=True)
    unique_rowcols = np.unique(lidar_rowcols_at_idx, axis=0)
    
    #use reducer to efficiently find mean z and assign to appropriate pixels in dsm/dem
    elev_means = np.add.reduceat(las_z_of_interest, unique_rowcol_ids_idx) / counts_in_pixels
    dsm_and_dem[unique_rowcols[:,0], unique_rowcols[:,1], image_idx] = elev_means
{% endhighlight %}

With the grid tile size and semaphore value I use (raster resolution doesn't matter as much), it takes me a little less than 3 hours to fill in the DEM & DSM. If tiles are not processed asynchronously, it takes maybe about a day. The large majority of this time is spent running *download_las()*, i.e. *ept.EPT.as_laspy()*. Everything before and after setting DEM/DSM values should take just a few seconds to run.

DSM cells should be approximately fully filled in within the workunit's boundaries (excepting pixels located over water bodies, see below), but DEM cells may not be filled in where no ground returns were downloaded. One way to partially resolve this issue to increase the requested EPT resolution, but this won't work for every case--for example, some pixels are completely contained by large rooftops. A more complete solution is to interpolate the DEM so missing values are filled in. Below, we use *scipy.interpolate.griddata()* to apply [cubic interpolation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CloughTocher2DInterpolator.html#scipy.interpolate.CloughTocher2DInterpolator) for the empty cells.  

{% highlight Python %}
#identify non-empty cells from both the DEM & DSM
dsm_not_nodata_idx = ~np.isnan(dsm_and_dem[:,:,0])
dem_not_nodata_idx = ~np.isnan(dsm_and_dem[:,:,1])

#identify cells with non-missing values for both the DEM and DSM
dsm_and_dem_filled_idx = dsm_not_nodata_idx & dem_not_nodata_idx
dsm_and_dem_filled_rowcols = np.column_stack(dsm_and_dem_filled_idx.nonzero())

#identify cells with non-missing values for the DSM but missing values for the DEM
dsm_filled_dem_not_idx = dsm_not_nodata_idx & ~dem_not_nodata_idx
dsm_filled_dem_not_rowcols = np.column_stack(dsm_filled_dem_not_idx.nonzero())

#interpolate DEM
dsm_and_dem[dsm_filled_dem_not_rowcols[:,0], dsm_filled_dem_not_rowcols[:,1], 1] = griddata(dsm_and_dem_filled_rowcols, 
                                                                                            dsm_and_dem[dsm_and_dem_filled_rowcols[:,0],
                                                                                                        dsm_and_dem_filled_rowcols[:,1]],
                                                                                            dsm_filled_dem_not_rowcols, method='cubic')
{% endhighlight %} 

Now we can create the DHM.  

{% highlight Python%}
#make DHM. np.maximum() is used to keep values from being < 0
heights = np.maximum(dsm_and_dem[:,:,0] - dsm_and_dem[:,:,1], 0) 
{% endhighlight %}

We'll also mask out pixels placed over bodies of water, according to NLCD 2016, since these can be misleading where values are greater than 0 (another option would have been to just set all water values to 0)

{% highlight Python%}
heights[nlcd_ndarray==1] = np.nan
{% endhighlight %}

<div style="text-align: center">
  <figure>
      <img
       src="/assets/DHM.png"
       width="482"
       height="413"
     />
     <figcaption>DHM. Values represent meters above the ground surface. </figcaption>
  </figure>
</div>\

In the near future I'll create a script which given a raster (ndarray and rasterio profile), a specified EPT dataset, and a set of desired pixel-level statistics to calculate from the EPT dataset returns a set of corresponding rasters, similar to what we've done here so far.

What's written below won't go over that sort of thing. Instead, it'll cover patch identification from the NLCD dataset, a technique for finding DHM pixels' surface areas and volumes, calculating the topographic roughness index, aggregating these measures of topography to the patch level, and how these measures seem to be related to LC class. The amount of code listed will be more limited than it is above, since this post is already on the lengthy side.

-----

Generally, patch identification is done using [connected component labeling algorithms](https://en.wikipedia.org/wiki/Connected-component_labeling). For example, [landscapemetrics](https://r-spatialecology.github.io/landscapemetrics/reference/get_patches.html) does this with a connected-component labeling algorithm that's also used in [watershed segmentation](https://people.cmm.minesparis.psl.eu/users/beucher/wtshed.html). 

[*scipy.ndimage.label()* implements one of these algorithms](https://scipy-user.scipy.narkive.com/V8Naooer/scipy-ndimage-measurements-label), and we'll take advantage of it here. The function will consider all non-zero values to be equivalent, so we'll have to apply it once per LC class. After initializing a patch ID raster, for each class we make a binary raster where True values correspond to the current class, get a patch ID raster for this class from *scipy.ndimage.label()*, then update the more general patch ID raster by adding the class-specific one to it. Note also that with each loop patch ID's are incremented with the current maximum value.

{% highlight Python%}
#label patches
max_patch_id = 0
patch_ids = np.zeros(shape=(img_height, img_width))
for lc_type in np.unique(nlcd_ndarray):
    if not np.isnan(lc_type):
        has_lc_type = nlcd_ndarray == lc_type
        patch_ids_for_lc_type, num_patches = label(has_lc_type)
        patch_ids += patch_ids_for_lc_type + has_lc_type * max_patch_id
        max_patch_id == num_patches
patch_ids[patch_ids == 0]
{% endhighlight %}

<div style="text-align: center">
  <figure>
      <img
       src="/assets/patch_ids.png"
       width="482"
       height="413"
     />
     <figcaption>Identified patches, each colored randomly. </figcaption>
  </figure>
</div>\

There are about 50,000 patches in the landscape, given the LC class definitions we set and using 4-neighborhood for each pixel. You might be able to tell by comparing the plot immediately above with the NLCD plot that the largest patches tend to be developed. If we liked, we could investigate further by comparing the distributions of classes' patch-level SHAPE index values. 

-----

It seems predictable that there would be clear relationships between LC type and height statistics--the most obvious such relationship being that forested and developed patches probably have greater heights above ground elevation, on average, than others. 

How can we investigate? One clear place to start is to just look at the distributions of corresponding DHM values for each LC class. Averages, of course, will suggest which types of LC tend to be shortest and tallest and so on, and if we'd like to measure how heights vary for each class we could just take variance. Another option for assessing height variability is the terrain ruggedness index (TRI), which is often used to describe pixel-level heterogeneity of DEM values (Riley et. al., 1999). 

The TRI value for a pixel can be expressed as: 

*TRI* = √[*&Sigma;*ᵢ₌₁ⁿ(*x*-*xᵢ*)²],

where *x* is the value of the pixel and *xᵢ* is the value of the pixel's ith neighbor. Of course the TRI wasn't designed for DHM data, but we'll compute it anyway just because I'd like to show how, and it should suit our purposes regardless. It's a pixel-level metric, so it can be aggregated to the patch or class level.  

Speaking of patches, *let's say we'd like to work at the patch level rather than the class level*...

-----
**References**

Jenness, J. S. (2004). Calculating landscape surface area from digital elevation models. Wildlife Society Bulletin, 32(3), 829-839.

Kedron, P., Zhao, Y., & Frazier, A. E. (2019). Three dimensional (3D) spatial metrics for objects. Landscape Ecology, 34, 2123-2132.

Riley, S. J., DeGloria, S. D., & Elliot, R. (1999). A terrain ruggedness index that quantifies topographic heterogeneity. intermountain Journal of sciences, 5(1-4), 23-27.

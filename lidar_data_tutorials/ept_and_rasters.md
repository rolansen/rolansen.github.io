**This is the second post in a series on using the geospatial Python ecosystem to analyze spatial relationships between lidar point clouds and other types of GIS datasets. This will be the first of maybe two or three on working with raster data. This post in particular will discuss how to query [Entwine Point Tile (EPT) datasets](https://entwine.io/en/latest/entwine-point-tile.html) and efficiently create raster data from the results.**

-----

Today most GIS software makes creating digital elevation/surface models (DEM's/DSM's) from lidar data relatively easy. While that's actually more or less what we're going to do here, doing this with Python gives us some obvious advantages, perhaps the most important being flexibility. 

Additionally, EPT support in most GIS software, with the notable exception of [QGIS](https://docs.qgis.org/testing/en/docs/user_manual/working_with_point_clouds/point_clouds.html), is for now nonexistent. Fortunately, the developers of the format have written a small Python library, [ept-python](https://github.com/hobu/ept-python), which allows us to conveniently create [laspy](https://laspy.readthedocs.io/en/latest/index.html) objects from queried EPT data. With the massive amount of [EPT data USGS is hosting](https://usgs.entwine.io/), there's a lot for us to take advantage of here.

In addition to ept-python, the user will need to install rasterio, geopandas (optionally pyogrio), shapely, [nest_asyncio,](https://pypi.org/project/nest-asyncio/) and some of the libraries typically included with scientific Python distributions: numpy, scipy, and matplotlib. We'll also be working with asyncio, so the user will need a recent version of Python. 

For convenience, ere's the full list of imports:

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

Our primary goal will be to create a digital height model (DHM) representing height above ground elevation. To prepare for downloading the lidar we'll use to create our DHM, we need to get the following:
* the boundary of the lidar survey of interest
* a "template" raster with a known CRS, spatial resolution, and extent that we'll sample the DHM to
* a grid that we'll use to make requests to the EPT resource

We'll need to start by extracting the survey boundary. USGS's [WESM layer](https://www.usgs.gov/ngp-standards-and-specifications/wesm-data-dictionary) is a source of spatial metadata for the lidar datasets it distributes, many of which are [provided in EPT format](https://usgs.entwine.io). We'll be working with the [work unit](https://www.usgs.gov/ngp-standards-and-specifications/wesm-data-dictionary-general-attributes#workunit) "NC_Phase4_Rowan_2017." I picked this one because it's relatively small and includes a diverse mix of cultivated, impervious, and forested land cover (see the second plot below).

The code below [pyogrio](https://pyogrio.readthedocs.io/en/latest/) to read the boundary to a single-row geopandas GeoDataFrame. This is fast and doesn't require us to write much code, but if you'd prefer you could do this [directly with geopandas if you have a recent version of Fiona](https://geopandas.org/en/stable/docs/user_guide/io.html#sql-where-filter), read the whole WESM layer and then subset it, or since we'll only be working with the work unit's geometry from here out read the row of interest with Fiona and convert the geometry to a Shapely polygon using shapely.wkb.loads().

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
     <figcaption>Boundary for work unit "NC_Phase4_Rowan_2017." CRS is EPSG:3857</figcaption>
  </figure>
</div>\

Now we'll use this boundary to define our raster of interest. We'll be using rasterio to read and write our raster data, allowing us to take advantage of the always-helpful numpy and scipy packages for	our image processing tasks. 

We could just manually define our raster's extent, CRS, and spatial resolution. However, it often it makes sense to sample a DHM to another raster of interest so they can be compared. Let's work with [National Land Cover Database (NLCD) data](https://www.mrlc.gov/), so later we can look a bit into any relationships between height and land cover (LC). We'll go with NLCD data from 2016, since it's close to the vintage of the lidar work unit. If you like you could [download the ERDAS .img file from the Multi-Resolution Land Characteristics Consortium](https://www.mrlc.gov/data/nlcd-2016-land-cover-conus), which will take up about 30 GB after being unzipped. I'll make a request from [their WMS site](https://www.mrlc.gov/geoserver/mrlc_display/wms?service=WMS&request=GetCapabilities), though, so we can get into a quick tangent about reading WMS data to rasterio.

Either way, we'll want to limit the data we read to the work unit's bounding box. Below I set parameters we'll use to make our WMS request--namely, the width, height, and bounding box--using the geometry of the work unit we just read. The CRS will be that which we projected this geometry to, EPSG:3857. The width and height are set so the spatial resolution is 30 meters and the work unit's bounding box is covered. Note the image will still be resampled because it's being reprojected and we aren't considering the CRS's origin. After setting the parameters, we request a GeoTIFF and just read it as a bytes object.

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
**This is the second post in a series on using the geospatial Python ecosystem to analyze spatial relationships between lidar point clouds and other types of GIS datasets. This will be the first of maybe two or three on working with raster data. This post in particular will discuss how to query [Entwine Point Tile (EPT) datasets](https://entwine.io/en/latest/entwine-point-tile.html) and efficiently create raster data from the results.**

-----

Today most GIS software makes creating digital elevation/surface models (DEM's/DSM's) from lidar data relatively easy. While that's actually more or less what we're going to do here, doing this with Python gives us some obvious advantages, perhaps the most important being flexibility. 

Additionally, EPT support in most GIS software, with the notable exception of [QGIS](https://docs.qgis.org/testing/en/docs/user_manual/working_with_point_clouds/point_clouds.html), is for now nonexistent. Fortunately, the developers of the format have written a small Python library, [ept-python](https://github.com/hobu/ept-python), which allows us to conveniently create [laspy](https://laspy.readthedocs.io/en/latest/index.html) objects from queried EPT data. With the massive amount of [EPT data USGS is hosting](https://usgs.entwine.io/), there's a lot for us to take advantage of here.

In addition to ept-python, the user will need to install rasterio, geopandas (optionally pyogrio), shapely, [nest_asyncio,](https://pypi.org/project/nest-asyncio/) and some of the libraries typically included with scientific Python distributions: numpy, scipy, and matplotlib. We'll also be working with asyncio, so the user will need a recent version of Python. 

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

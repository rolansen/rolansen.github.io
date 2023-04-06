**This is the second post in a series on using the geospatial Python ecosystem to analyze spatial relationships between lidar point clouds and other types of GIS datasets. This will be the first of maybe two or three on working with raster data. This post in particular will discuss how to query [Entwine Point Tile (EPT) datasets](https://entwine.io/en/latest/entwine-point-tile.html) and efficiently create raster data from the results.**

-----

Today most GIS software makes creating digital elevation/surface models (DEM's/DSM's) from lidar data relatively easy. While that's actually more or less what we're going to do here, doing this with Python gives us some obvious advantages, perhaps the most important being flexibility. 

Additionally, EPT support in most GIS software, with the notable exception of [QGIS](https://docs.qgis.org/testing/en/docs/user_manual/working_with_point_clouds/point_clouds.html), is for now nonexistent. Fortunately, the developers of the format have written a small Python library, [ept-python](https://github.com/hobu/ept-python), which allows us to conveniently create [laspy](https://laspy.readthedocs.io/en/latest/index.html) objects from queried EPT data. With the massive amount of [EPT data USGS is hosting](https://usgs.entwine.io/), there's a lot for us to take advantage of here.

In addition to ept-python, the user will need to install rasterio, geopandas (optionally pyogrio), shapely, [nest_asyncio,](https://pypi.org/project/nest-asyncio/) and some of the libraries typically included with scientific Python distributions: numpy, scipy, and matplotlib. We'll also be working with asyncio, so the user will need a recent version of Python. 

-----

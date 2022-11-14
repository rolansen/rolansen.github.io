#change directories and read polygons
import geopandas as gpd
import numpy as np
import os

os.chdir(r'D:\HeDemo\data\final')
tree_polys_filename = 'canopy_polygons.geojson'

tree_polys = gpd.read_file(tree_polys_filename)
tree_polys['max_height'] = np.nan
tree_polys['lai'] = np.nan

#---------------------------------------------------------------------------------------------------------------------------
#read lidar data
import laspy

lidar_filename = 'USGS_LPC_NY_3County_2019_A19_e1382n2339_2019.laz'
scan_angle_factor = 0.006 #to convert from value given by las object, see https://github.com/ASPRSorg/LAS/issues/41#issuecomment-344300998.

las = laspy.read(lidar_filename)
np_las = np.transpose(np.vstack([las.x, las.y, las.z, las.classification, las.scan_angle*scan_angle_factor, np.zeros(len(las))]))

#------------------------------------------------------------------------------------------------------------------------------------------------
#make KDTree
from scipy.spatial import KDTree
kdtree_las_xy = KDTree(np_las[:, :2])

#---------------------------------------------------------------------------------------------------------------------------
#find maximum heights and LAI's within polys
#ie for each tree polygon, do range query, calculate HAG's, find points in polygon, calculate LAI
from scipy.interpolate import LinearNDInterpolator

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

def get_height_above_ground_for_points(points):
    '''Estimate height above ground for a (ground-classified) point cloud, using TIN interpolation.
    see https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html'''
    ground_las_class_code = 2
    ground_points = points[points[:,3] == ground_las_class_code]
    ground_interp = LinearNDInterpolator(points = ground_points[:,:2], values = ground_points[:,2])
    ground_elevations = ground_interp(points[:,:2])
    return points[:,2] - ground_elevations

search_radius_buffer = 2
lambert_beer_extinction_coefficient_when_scan_angle_is_0 = 0.5 #assumes spherical leaf angle distribution.
ground_elev_threshold = 0.05 #in m
for index, poly in tree_polys.iterrows():
    print(index)
    #do range query
    poly_mbr = np.array(poly['geometry'].bounds) #minx, miny, maxx, maxy
    poly_mbr_centroid = np.array([(poly_mbr[2] + poly_mbr[0]) / 2, (poly_mbr[3] + poly_mbr[1]) / 2])
    poly_search_radius = search_radius_buffer + ((poly_mbr[2] - poly_mbr_centroid[0])**2 + (poly_mbr[3] - poly_mbr_centroid[1])**2) ** 0.5
    near_poly_np_las = np_las[kdtree_las_xy.query_ball_point(poly_mbr_centroid, poly_search_radius)]

    #find heights above ground
    if len(near_poly_np_las) > 0:
        near_poly_np_las[:, 5] = get_height_above_ground_for_points(near_poly_np_las)
    else:
        continue

    #find points in polygon
    poly_vertex_array = np.array(poly['geometry'].exterior.coords)
    in_poly_np_las = near_poly_np_las[[point_in_polygon(las_point, poly_vertex_array) for las_point in near_poly_np_las]]

    #assign max height in polygon
    if len(in_poly_np_las) > 0:
        tree_polys.at[index, 'max_height'] = np.max(in_poly_np_las[~np.isnan(in_poly_np_las[:, 5]), 5])
    else:
        continue

    #calculate LAI
    total_number_of_returns_in_poly = len(in_poly_np_las)
    number_of_ground_returns_in_poly = np.sum(in_poly_np_las[:,5] <= ground_elev_threshold)
    mean_lidar_scanning_angle = np.deg2rad(np.mean(in_poly_np_las[:,4]))
    tree_polys.at[index, 'lai'] = -np.cos(mean_lidar_scanning_angle) / lambert_beer_extinction_coefficient_when_scan_angle_is_0 * np.log(number_of_ground_returns_in_poly / total_number_of_returns_in_poly)

#---------------------------------------------------------------------------------------------------------------------------
#write results
output_path = 'canopy_polygons_with_height_and_lai.geojson'
tree_polys.to_file(output_path)

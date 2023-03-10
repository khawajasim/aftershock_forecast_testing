#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 15:20:42 2022

@author: khawaja
"""
# import json
# import matplotlib.pyplot as plt
import numpy as np
import numpy
from shapely.geometry import Point, Polygon, box #GeometryCollection, LineString, mapping, 
# from shapely.ops import cascaded_union, polygonize, unary_union
import copy

def get_circumference_points(radius, number_of_points):
    """
    Get the radius of circle, and provive the number of points along circumference. 
    """
    
    radians_between_each_point = 2*np.pi/number_of_points
    list_of_points = []
    for p in range(0, number_of_points):
        list_of_points.append([radius*np.cos(p*radians_between_each_point),radius*np.sin(p*radians_between_each_point)])
    return list_of_points


# def write_geojson(grid_polygons,filename):
#     """
#     This functions takes in Shape.Geometery type list of shapes and writes
#     them in the format of geojason file. The geojason file will be then used to plot in QGIS

#     Parameters
#     ----------
#     shape_list : List of Shapely.geometery
#         Bounding boxes or points or lines that are used for display in QGIS
#     filename : string
#         Intended name of GEOJSON file

#     Returns 
#     -------
#     GEOJSON file

#     """
#     shape_list=[] #Converting to the proper format. 
#     for poly in grid_polygons:
#         shape_list.append(mapping(poly))
        
        
#     geo_shape = {"type":"FeatureCollection",
#                   "features":[]}
#     for bb in shape_list:
# #        print(bb.centroid)
#         feature = {"type":"Feature",
#               "properties":{"some_parameter": "XYZ",
#                             },
#                   "geometry":bb};
#         geo_shape['features'].append(feature)

# #   poly['features'][0]['geometry']['coordinates']

#     with open(filename+'.geojson', 'a') as f:
#         json.dump(geo_shape, f)
    
#     return geo_shape

# def write_geojson_feature(grid_polygons,feature, filename):
#     """
#     Write coordsniates along with feature
#     This functions takes in Shape.Geometery type list of shapes (polygon or point) and writes
#     them in the format of geojason file. The geojason file will be then used to plot in QGIS

#     Parameters
#     ----------
#     shape_list : List of Shapely.geometery
#         Bounding boxes or points or lines that are used for display in QGIS
#     filename : string
#         Intended name of GEOJSON file

#     Returns 
#     -------
#     GEOJSON file

#     """
#     shape_list=[] #Converting to the proper format. 
#     for poly in grid_polygons:
#         shape_list.append(mapping(poly))
        
        
#     geo_shape = {"type":"FeatureCollection",
#                   "features":[]}
#     for bb, ff in zip(shape_list, feature):
#         feature = {"type":"Feature",
#               "properties":{"feature1": ff,
#                             },
#                   "geometry":bb};
#         geo_shape['features'].append(feature)

# #   poly['features'][0]['geometry']['coordinates']

#     with open(filename+'.geojson', 'w') as f:
#         json.dump(geo_shape, f)
    
#     return geo_shape



def create_full_arc(radius,start_angle, end_angle, origin=(0,0), numsegments = 100):
    """
    Create a polygon circular arc for radius, between start_angle (degrees) to end_angle (degrees).
    Inputs:
        radius of circle
        start_angle, end_angle: in degrees
        origin: tuple the center of cirlce. Default (0,0)
    """
    theta = np.radians(np.linspace(start_angle, end_angle, numsegments))
    x = origin[0] + radius * np.cos(theta) #I think There is a potential round off error here. 
    y = origin[1] + radius * np.sin(theta)
#    xx = np.append(np.append(origin[0],x),origin[0])
#    yy =  np.append(np.append(origin[1],y),origin[1])
#    arc = Polygon(np.column_stack([xx, yy]))
    arc = Polygon(np.column_stack([x, y]))
    return arc



def create_circular_grid(max_radius, r_segments, a_segments, origin = (0,0) ):
    """
    Get circular grid, that divides "max_radius" into "r_segments" concentric circles and different 
    angles.
    Inputs: 
        max_radius: Maximum radius of circular grid
        r_segments: Number of concentric circles, or the number of segments in circle.
        a_segments: The number circular segments. Or segments in the angle.
        origin = the center point of the circular grid
    Returns:
        A list of polygon for each circular cell. 
    """

    radius = np.linspace(0,max_radius, r_segments+1)[1:] #,4,6,8,10]
    angles = np.linspace(0,360, a_segments+1)
    #angles = [0,90]
    polygon_list = [] # list of al the arcs.
    
    for i in range(len(angles)-1):
#        print('Arc between Angle :', angles[i],' and ', angles[i+1])
    
        #Frst take one angle segment and find arcs for all the radius.
        #the grid cell at every radius is the difference between this segment and all smaller radius segments
        for r in np.sort(radius):
#            print('The radius of Arc : ',r)
            arc = create_full_arc(r,angles[i], angles[i+1], origin=origin)
            if r == min(radius):
                coords = arc.exterior.coords[:]
                coords.pop(-1)
                polygon_list.append(Polygon([origin]+coords))
                arc_tmp = copy.copy(arc)
            else: #elif r == max(radius)
    #            diff = arc.difference(arc_tmp)
    #            polygon_list.append(diff)
                coords1 = arc_tmp.exterior.coords[:]
                coords2 = arc.exterior.coords[:]
                coords1.pop(-1)
                coords2.pop(-1)
                coords2.reverse()
                polygon_list.append(Polygon(coords1+coords2))

                arc_tmp = copy.copy(arc)
    return polygon_list        



def bounds_to_polygons(bounding_coords):
    """
    taking corner points and convert them to polygons boxes
    Inputs:
        List or Arrany of points: nx4.
        [bottom_left_x, bottom_left_y, top_right_x, top_right_y]
    """
    bounding_polygons = []
    for bb in bounding_coords:
        bounding_polygons.append(Polygon(box(bb[0], bb[1], bb[2], bb[3])))
    return bounding_polygons


def create_square_grid_bounds(bottom_left_corner, top_right_corner, ncells_1D):
    """
    Create a square grid with starting from bottom_left_corner to top_right_corner with dh_1D cells in one side
    dh_1D = sqrt(Total_cells_in_grid)

    """
    xx = np.linspace(bottom_left_corner[0], top_right_corner[0],ncells_1D+1)
    yy = np.linspace(bottom_left_corner[1],top_right_corner[1], ncells_1D+1)
    
    bottom_left_coords = []
    for i in range(len(xx)-1):
        for j in range(len(yy)-1):
            bottom_left_coords.append([xx[i],yy[j]])
    
    bottom_left_coords = np.array(bottom_left_coords)
    dh = np.unique(np.diff(xx))[0]
    top_right_coords = bottom_left_coords + dh
    
    bounding_coords = np.column_stack((bottom_left_coords, top_right_coords))
    return bounding_coords
    

def forecast_aggregation_mean(test_bounds, model_loc, stress_value):
    """
    Forecast mapping from "model_locations" to the "test_bounds"
    Take the mean of all the points inside the "test_bounds", and then multiply by the area. 
    Input:
        test_bounds: A list of polygons of test grid
        model_loc: A list of polygons of forecast grid
        stress_values: forecast
    """
#    model_bounds = [[gg.bounds[0], gg.bounds[1]] for gg in model_grid]
#    model_bounds = numpy.array(model_bounds)
    test_stress = [] #numpy.zeros(len(test_bounds))
    
    for bounds in test_bounds:
        inside_locs = numpy.logical_and(numpy.logical_and(bounds[0] <= model_loc[:,0],model_loc[:,0] < bounds[2]),
                            numpy.logical_and(bounds[1] <= model_loc[:,1], model_loc[:,1] < bounds[3]))
        
#        len(model_loc[inside_locs])
        cell_area = abs((bounds[0]-bounds[2]) * (bounds[1] - bounds[3]))
        test_stress.append(sum(stress_value[inside_locs])/len(stress_value[inside_locs]) * cell_area)
     
    return test_stress


def forecast_aggregation(test_grid, model_grid, stress_value, by_addition=False):
    """
    Forecast mapping from one grid "model_grid" to the "test_grid"
    Input:
        test_grid: A list of polygons of test grid
        model_grid: A list of polygons of forecast grid
        stress: forecast
    Process: If by_addition=True, then add the fraction of stress shared by the common areas.
            If by_addition=False, then average stress of all the points touching a cell. 
            (Since, the stress is computed for the points, so I find the latter better for aggregation.)
    """
    test_stress = numpy.zeros(len(test_grid))
#    points_taken = 0
    if by_addition:
        for i, p1 in enumerate(test_grid):
            # print('Cell: ',i)
            for j, p2 in enumerate(model_grid):
                if p1.intersects(p2):
                   test_stress[i] = test_stress[i]+ ((p1.intersection(p2).area /p2.area) * stress_value[j] )
    else:
        model_points = [[bb.bounds[0], bb.bounds[1]] for bb in model_grid]
        for i, p1 in enumerate(test_grid):
#            print('Cell: ',i)
            l = 0
            for j, p2 in enumerate(model_points):
                if p1.intersects(Point(p2)):
                   l=l+1 
                   test_stress[i] = test_stress[i]+ stress_value[j] #((p1.intersection(p2).area /p2.area) * stress[j] )
            test_stress[i] = test_stress[i]/l * p1.area
            # print('Cell ',i, ' took ',l, ' Points')
#            points_taken = points_taken+l
        
          
#    for n,pn in enumerate(test_grid):  #Normalizing the aggregated forecasts
#        test_stress[n] = test_stress[n] / pn.area
               
    return test_stress

def forecast_aggregation_bounds(grid_bounds, stress_loc,stress_data):
    """
    Aggregating the stress into the grid_bounds
    Inputs:
        grid_bounds[lon1,lat1, lon2,lat2] nx4
        stress_loc [lon,lat] nx2
        stress_data nx1: MAS or OOP. 
        
    Outputs:
        Aggregated stress: nx1
    
    """
    aggregated_forecast = numpy.zeros(len(grid_bounds))
    for i, bb in enumerate(grid_bounds):
        inside_points = numpy.logical_and(numpy.logical_and(bb[0]<=stress_loc[:,0],bb[1]<=stress_loc[:,1]),
                                numpy.logical_and(stress_loc[:,0]<bb[2], stress_loc[:,1]<bb[3]))
        aggregated_forecast[i] = sum(stress_data[inside_points])
    
    return aggregated_forecast    

def locate_eq_to_grid(grid, eq_loc):
    """
    Map xy location of earthquakes onto a grid. Convert catalog into gridded catalog
    Inputs:
        grid: List of shapely.polygon of grid cells
        eq_loc: List of xy location of earthqyuakes
    Outputs:
        6
    """
    gridded_cat = numpy.zeros(len(grid))
    for xy in eq_loc:
        for i, poly in enumerate(grid):
            if poly.contains(Point(xy[0], xy[1])):
                gridded_cat[i] = gridded_cat[i]+1
                break
    return gridded_cat    

def locate_eq_to_grid3D(stress_coords, cat): 
    """
    Take catalog and Stress bounds in 3D (Lon,lat,depth). Locate each earthquake of catalog to the respective 3D grid.
    Inputs: 
        stress_coords 
            nx6: [lon1,lat1,lon2,lat2,d1,d2]
        cat
            nx4: [lat,laon,mag, depth]
    """
    lon1 = stress_coords[:,0]
    lon2 = stress_coords[:,2]
    lat1 = stress_coords[:,1]
    lat2 = stress_coords[:,3]
    depth1 = stress_coords[:,4]
    depth2 = stress_coords[:,5]
    
    grid_eq = []
    # cat = AScataBuffer[ AScataBuffer[:,2] >= Mcut ]
    for i in range(len(stress_coords)):
        flag_coords = numpy.logical_and(
            numpy.logical_and(cat[:,1]>=lon1[i], cat[:,0]>=lat1[i]), 
            numpy.logical_and(cat[:,1] < lon2[i], cat[:,0] < lat2[i]))
        flag_coords_depth = numpy.logical_and( flag_coords, numpy.logical_and(depth1[i]<= cat[:,3], cat[:,3]<depth2[i]))
        grid_eq.append(len(cat[flag_coords_depth]))
    
    grid_eq = numpy.array(grid_eq) 
    return grid_eq
# Robust Normal Estimation in point clouds

## Paper

* Title: 	Iterative weighted PCA for robust and edge-aware normal vector estimation.
* Authors: 	Julia Sanchez, Florence Denis, David Coeurjolly, Laurent Trassoudaine, Paul Checchin.
* Year: 	2019.

## Dependencies

This filter have a Flann dependency, it means that it can be applied to any data structure but need to be compiled with PCL. 
The original algorithm was built with OpenMp. It was disabled during the port.

### Parameters

#### Number of neighbours

Number of points taken into account for local surface fitting. 
Low value leads to a faster computation time but can decrease quality in case of important sampling anisotropy or noise level. 

Typical: 100 
For scattered point clouds : 50 
For very noisy point clouds : 300 

-------------------------------------------------------------------- 

#### Noise

The sensor noise, in meters. 

For bad sensor: 0.03 
For good sensor : 0.005 
For artificial point cloud: 0 

-------------------------------------------------------------------- 

#### Curvature

Maximum curvature radius to define a sharp feature, in meters. 

For curved object: 0.1 
For piecewise planar objects: inf 

### Suggested improvement 

Actually the original algorithm is just "plugged" into a generic MEPP2 filter, there is a duplication of data between the generic cloud write and read by MEPP and the one used by the computation ( this data structure is named "raw cloud").

A suggested improvement can be to propagate the generic cloud through all the code in order to convert the full algorithm into a generic one, and avoid data duplication just as flann dependency.

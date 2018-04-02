# SurfaceFromCloud
Surface generarating from pointcloud and texture with camera images

# Dependencies:
+ Boost 1.66.0
+ Qt 5.10.1
+ CGAL (this has too much dependency)
+ Eigal (not necessary)

# Progress
This project is able to make surface from point clouds using poisson surface reconstructor by CGAL and own implemented Power Crust algorithm
The result is a triangle mesh and you can check the result of Power Crust steps as well like display Voronoi Diagramm / Cell,
Power Diagram / Cell or power_crust surface or the power_shape. 
There are some simplifier, outlier removal and smoother filters as well what can be used on the cloud.
Right now we cannot use mine kdtree implementation.

The cloud files are in .ply format
ply, vertex number and ending terms are needed only
data format: x -y z r g b
                ?? yes, I know
                
Incoming:
+ Make better performance for Power Crust (multithread / GPU)
+ Make better memory usage for Power Crust (as well we can display show many things...)
+ Add texture coordinates to the result (using pictures provided by the kinect)
+ Add one more surface reconstructor provided by CGAL
+ Refactoring some things...
+ Let the user to use my kdtree
+ adding config file
+ ...

# SurfaceFromCloud
Surface generarating from pointcloud and texture with camera images

# Dependencies:
+ Boost 1.66.0
+ Qt 5.10.1
+ CGAL (this has too much dependency)
+ Eigal (not necessary)

# Progress
This project just have started.
Right now you can use the CGAL functions to get the surface, you can add / remove filters and change their properties whenever you want. I added one more outlier "algorithm". I also added one kd tree implementation but that currantly cant be used.

The cloud files are in .ply format
ply, vertex number and ending terms are needed only
data format: x -y z r g b
                ?? yes, I know
                
Incoming:
+ Voronoi Diagram calculation and display
+ http://web.cs.ucdavis.edu/~amenta/pubs/sm.pdf proceed
+ Let the user to use my kdtree
+ Linking CGAL with TBB to get better performance (just test needed)
+ adding config file
+ ...

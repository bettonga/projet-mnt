# DEM Project
A program that creates a raster of a Digital Elevation Model given a series of measurement.

## Requirements
* CMake
* Proj 8.1.1 (download packet [here](https://download.osgeo.org/proj/proj-8.1.1.tar.gz), installation guide for CMake [here](https://proj.org/install.html))
* CGAL (```sudo apt-get install libcgal-dev```)

## How to use it
Once cloned, open a terminal in the repository and execute ```./build.sh```, then ```cd build/```.

The program can now be executed
```
./create_raster <your_file> <img_width>
```

The image will be created in the ```/render``` folder.

## Note
The file must be of this format :

```
#latitude   longitude     depth
48.29762887 -004.41737340 14.334
48.29762971 -004.41735997 14.379
48.29763698 -004.41738809 14.452
48.29763783 -004.41737466 14.370
48.29763867 -004.41736124 14.376
48.29763951 -004.41734781 14.463
48.29764425 -004.41741620 14.302
48.29764510 -004.41740278 14.341
48.29764594 -004.41738935 14.459
48.29764678 -004.41737593 14.402
```

With WGS84 coordinates (standard system used by GPS).

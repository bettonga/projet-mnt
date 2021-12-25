#include <iostream>
#include <fstream>
#include <proj.h>
#include <cstdio>
#include <math.h>
#include <chrono>
#include <ctime>
#include <string>
#include <array>
#include <thread>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;

typedef Delaunay::Face_handle     Face_handle;
typedef Delaunay::Vertex_handle   Vertex_handle;
typedef K::Point_3                Point_3;
typedef K::Vector_3               Vector_3;

using namespace std;


#define MAX_AREA  10.0      // define the maximum area of a triangle
#define MAX_PGM   255       // grey level coded on 1 byte
#define MIN_PGM   10        // >0 so that the background does not blend with the MNT
#define BIN       true      // true for binary image file output
#define NB_THREAD 4         // number of threads creating the image (4 seems to be the most effective, even for larger file | set at 1 for non-parallel processing)

// boundaries of the input terrain
double min_x, max_x, min_y, max_y, min_z, max_z;


// project GPS datas into cartesian coordinates, and insert them into a Delaunay 2D triangulation
int proj93(Delaunay& dt, char* file_name);

// return the depth z according to a (x,y) point and the triangulation
double get_z(const Delaunay& dt, double x, double y, Face_handle& old_fh);

// updates min_x, max_x, min_y, max_y, min_z, max_z according to the boundaries of the input terrain
void update_maxmin(double& min_x, double& max_x, double& min_y, double& max_y, double& min_z, double& max_z, double new_x, double new_y, double new_z);

// function and thread associated to create the PGM image
void create_pgm(const int n, const Delaunay& dt, const int& img_width, const int& img_height, const double& density, const bool& binary);
void thread_pgm(const int n, const int k, const Delaunay& dt, const int& img_width, const int& img_height, const double& density, const bool& binary);



int main(int argc, char *argv[]) {
  // chrono start
  auto start = std::chrono::system_clock::now();

  // projection GPS WGS84 -> cartesian xyz + Delaunay 2D triangulation
  Delaunay dt;
  proj93(dt, argv[1]);

  // calculation of the final image dimension and pixel density (pixel/m)
  const int img_width = strtol(argv[2], NULL, 10);
  const double density = img_width / (max_x - min_x);
  const int img_height = ceil(density * (max_y - min_y));

  // PGM image creation, result will be stored in /img_output
  create_pgm(NB_THREAD, dt, img_width, img_height, density, BIN);

  // chrono end
  auto end = std::chrono::system_clock::now();
  chrono::duration<double> elapsed_seconds = end-start;
  cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";

  return EXIT_SUCCESS;
}






int proj93(Delaunay& dt, char* file_name) {
  //// project GPS datas into cartesian coordinates, and insert them into a Delaunay 2D triangulation

  // projection initialization, will project GPS WGS84 -> cartesian xyz
  PJ *P;
  PJ_COORD c, c_out;
  P = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
                             "+proj=longlat +datum=WGS84",                     // (x_0, y_0) is the offset point
                             "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
                             NULL);
  PJ* P_for_GIS = proj_normalize_for_visualization(PJ_DEFAULT_CTX, P);
  proj_destroy(P);
  P = P_for_GIS;


  // projection and triangulation
  fstream data;
  data.open(file_name, fstream::in);
  double lat, lon, z;
  bool first = true;
  while (data >> lat >> lon >> z){          // read through all lines
    // projection of the points
    c.lpz.lam = lat;
    c.lpz.phi = lon;
    c.lpz.z = z;
    c_out = proj_trans(P, PJ_FWD, c);
    if (first == true) {                    // min_x, max_x, min_y, max_y, min_z, max_z initialization on the first point
      min_x = c_out.xyz.x; max_x = c_out.xyz.x; min_y = c_out.xyz.y; max_y = c_out.xyz.y; min_z = c_out.xyz.z; max_z = c_out.xyz.z;
      first = false;
    }
    update_maxmin(min_x, max_x, min_y, max_y, min_z, max_z, c_out.xyz.x, c_out.xyz.y, c_out.xyz.z);
    // triangulation insertion
    Point_3 p(c_out.xyz.x, c_out.xyz.y, c_out.xyz.z);
    dt.insert(p);
  }

  // close file and destroy projection
  proj_destroy(P);
  data.close();
  return 0;
}


double get_z(const Delaunay& dt, double x, double y, Face_handle& old_fh){
  //// return the depth z according to a (x,y) point and the triangulation

  // location of a valid triangle containing the (x,y) point
  Point_3 p(x,y,0.0);
  Face_handle fh = dt.locate(p, old_fh);                                    // old_fh is the starting point of locate method
  if (fh==nullptr || dt.is_infinite(fh)) {throw invalid_argument("");}      // exit if infinite or null triangle (outside of the terrain)
  old_fh = fh;                                                              // old_fh is updated with the new triangle, next locate will start from this triangle
  Point_3 a = fh->vertex(0)->point();
  Point_3 b = fh->vertex(1)->point();
  Point_3 c = fh->vertex(2)->point();
  Vector_3 u(a,b);
  Vector_3 v(a,c);
  if (0.5*(cross_product(u,v).squared_length()) > MAX_AREA) {throw invalid_argument("");}     // exit if triangle area is too large (not really on the terrain)

  // interpolation : the (x,y,0) point is projected according to the vertical vector (0,0,1) on the face defined by the triangle located before
  double denom = (c[0]-a[0])*(b[1]-a[1]) - (c[1]-a[1])*(b[0]-a[0]);
  double mu = (  (b[1]-a[1])*(x-a[0]) - (b[0]-a[0])*(y-a[1]) ) / denom;
  double la = ( -(c[1]-a[1])*(x-a[0]) + (c[0]-a[0])*(y-a[1]) ) / denom;
  double z = a[2] + la*(b[2]-a[2]) + mu*(c[2]-a[2]);
  return z;
}


void update_maxmin(double& min_x, double& max_x, double& min_y, double& max_y, double& min_z, double& max_z, double new_x, double new_y, double new_z){
  //// updates min_x, max_x, min_y, max_y, min_z, max_z according to the boundaries of the input terrain
  if (new_x<min_x) {min_x = floor(new_x);}
  else if (new_x>max_x) {max_x = ceil(new_x);}
  if (new_y<min_y) {min_y = floor(new_y);}
  else if (new_y>max_y) {max_y = ceil(new_y);}
  if (new_z<min_z) {min_z = new_z;}
  else if (new_z>max_z) {max_z = new_z;}
}


void create_pgm(const int n, const Delaunay& dt, const int& img_width, const int& img_height, const double& density, const bool& binary) {
  //// function to create the PGM image

  // generate file_name automatically according to date and hour
  time_t rawtime;
  struct tm * ptm;
  time ( &rawtime );
  ptm = gmtime ( &rawtime );
  string file_name = "../img_output/" + to_string(ptm->tm_mday) + "-" + to_string(ptm->tm_mon) + "_" + to_string(ptm->tm_hour) + "-" + to_string(ptm->tm_min) + "_render.pgm";

  // writes PGM file header (P5 for binary, P2 for ASCII)
  fstream img;
  if (binary) {img.open(file_name, fstream::out | fstream::binary); img<<"P5 ";}
  else {img.open(file_name, fstream::out); img << "P2 ";}
  img << img_width << " "
      << img_height << " "
      << (int) MAX_PGM << endl;

  // start image creation, each thread process (1/n)th of the image, then all tmp file are concatenated
  array<thread,NB_THREAD> threads;

  for (int k=0; k<n; k++) {   // starts the threads
    threads[k] = thread( [=] { thread_pgm(n, k, dt, img_width, img_height, density, binary); } ); }

  for(int k=0; k<n; k++) {    // stops them
    threads[k].join(); }

  for(int k=0; k<n; k++) {   // concatenation of all n tmp files
    fstream img_mrg;
    file_name = "/tmp/MNT_tmp" + to_string(k);
    img_mrg.open(file_name, fstream::in);

    img << img_mrg.rdbuf();
    img_mrg.close();
  }

  // image is created, stream can be closed
  img.close();
}


void thread_pgm(const int n, const int k, const Delaunay& dt, const int& img_width, const int& img_height, const double& density, const bool& binary) {
  //// thread called to create the PGM image
  //cout << "Thread " + to_string(k) + " starting" << endl;

  // opening of a tmp file for writing
  fstream img_tmp;
  string file_name;
  file_name = "/tmp/MNT_tmp" + to_string(k);    // each file_name is different and placed in /tmp folder (deleted at each startup)
  if (binary) {img_tmp.open(file_name, fstream::out | fstream::binary);}
  else {img_tmp.open(file_name, fstream::out);}

  // calculation of the section the k-th thread has to process
  int len = ceil(img_height / n);
  int i_min = k*len;
  int i_max;
  if (k==n-1) {i_max = img_height;}   // last thread takes the rest
  else {i_max = (k+1)*len;}

  // calculation of transformation parameters from z -> grey level
  const double coeff = (MAX_PGM-MIN_PGM)/(max_z-min_z);
  const double offset = MIN_PGM - coeff*min_z;

  // browse pixel and write color correspondance on the image file
  double img_z;
  int pgm;
  const int pgm_null = 1;             // background color, 0 cause write method to skip value
  Face_handle old_fh = NULL;

  for (int i=i_min; i<i_max; i++) {
    for (int j=0; j<img_width; j++) {
      // color calculation if pixel point on a valid triangle (non-infinite, non-null, area < MAX_AREA)
      try {
        img_z = get_z(dt, min_x+j/density, max_y-i/density, old_fh);
        pgm = round(coeff*img_z+offset);
        if (binary) {img_tmp << reinterpret_cast<char*>(&pgm);}
        else {img_tmp << pgm << " ";}
      }
      catch (invalid_argument& e) {
        if (binary) {img_tmp << reinterpret_cast<char*>(&pgm_null);}
        else {img_tmp << pgm_null << " ";}
      }
      }
      if (!binary) {img_tmp << endl;}
    }

  // image is created, stream can be closed
  img_tmp.close();

  //cout << "Thread " + to_string(k) + " finished" << endl;
}

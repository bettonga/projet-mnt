/**
 * \file main.cpp
 * \brief DEM Project
 * \author Gabriel Betton
 * \date 9/01/2022
 *
 * A program that creates a raster of a Digital Elevation Model given a series of measurement.
 *
 */

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
#include <mutex>
#include <cstddef>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;

typedef Delaunay::Face_handle     Face_handle;
typedef Delaunay::Vertex_handle   Vertex_handle;
typedef Delaunay::Line_face_circulator   Line_face_circulator;
typedef K::Point_3                Point_3;
typedef K::Vector_3               Vector_3;

using namespace std;


/** \brief define the maximum area of a triangle */
#define MAX_AREA 10.0
/** \brief grey level coded on 1 byte */
#define MAX_PGM 255
/** \brief >0 so that the background does not blend with the minimal depth | can be 0 for PPM with colormap */
#define MIN_PGM 0
/** \brief sun altitude, 0 for sun on the horizon, 90 for vertical */
#define ALTITUDE 45
/** \brief sun azimut, 0 for North, 90 for East, 180 for South, 270 for West */
#define AZIMUT 315
/** \brief true for PPM, false for PGM */
#define PPM true
/** \brief true for binary file, false for ASCII */
#define BIN true
/** \brief true for hill-shadind */
#define SHD true
/** \brief number of threads creating the image */
#define NB_THREAD 4



// boundaries of the input terrain
double min_x, max_x, min_y, max_y, min_z, max_z;
mutex foo, bar;

// Haxby colormap
double tab_pixm [256][3]= {{0.1451,0.22353,0.68627},{0.14556,0.23429,0.69796},{0.14602,0.24506,0.70965},{0.14648,0.25582,0.72134},{0.14694,0.26659,0.73303},{0.1474,0.27735,0.74471},{0.14787,0.28812,0.7564},{0.14833,0.29889,0.76809},{0.14879,0.30965,0.77978},{0.14925,0.32042,0.79146},{0.14971,0.33118,0.80315},{0.15017,0.34195,0.81484},{0.15063,0.35271,0.82653},{0.1511,0.36348,0.83822},{0.15156,0.37424,0.8499},{0.15202,0.38501,0.86159},{0.15248,0.39577,0.87328},{0.15294,0.40654,0.88497},{0.1534,0.4173,0.89666},{0.15386,0.42807,0.90834},{0.15433,0.43883,0.92003},{0.15479,0.4496,0.93172},{0.15525,0.46036,0.94341},{0.15571,0.47113,0.95509},{0.15617,0.48189,0.96678},{0.15663,0.49266,0.97847},{0.15763,0.50288,0.98462},{0.15917,0.51257,0.98524},{0.16071,0.52226,0.98585},{0.16225,0.53195,0.98647},{0.16378,0.54164,0.98708},{0.16532,0.55133,0.9877},
{0.16686,0.56101,0.98831},{0.1684,0.5707,0.98893},{0.16993,0.58039,0.98954},{0.17147,0.59008,0.99016},{0.17301,0.59977,0.99077},{0.17455,0.60946,0.99139},{0.17609,0.61915,0.992},{0.17762,0.62884,0.99262},{0.17916,0.63852,0.99323},{0.1807,0.64821,0.99385},{0.18224,0.6579,0.99446},{0.18378,0.66759,0.99508},{0.18531,0.67728,0.99569},{0.18685,0.68697,0.99631},{0.18839,0.69666,0.99692},{0.18993,0.70634,0.99754},{0.19146,0.71603,0.99815},{0.193,0.72572,0.99877},{0.19454,0.73541,0.99938},{0.19608,0.7451,1},{0.20469,0.75202,1.},{0.2133,0.75894,1.},{0.22191,0.76586,1.},{0.23053,0.77278,1.},{0.23914,0.7797,1.},{0.24775,0.78662,1.},{0.25636,0.79354,1.},{0.26498,0.80046,1.},{0.27359,0.80738,1.},{0.2822,0.8143,1.},{0.29081,0.82122,1.},{0.29942,0.82814,1.},
{0.30804,0.83506,1.},{0.31665,0.84198,1.},{0.32526,0.8489,1.},{0.33387,0.85582,1.},{0.34248,0.86275,1.},{0.3511,0.86967,1.},{0.35971,0.87659,1.},{0.36832,0.88351,1.},{0.37693,0.89043,1.},{0.38554,0.89735,1.},{0.39416,0.90427,1.},{0.40277,0.91119,1.},{0.41138,0.91811,1.},{0.41815,0.92165,0.99377},{0.42307,0.9218,0.98131},{0.42799,0.92195,0.96886},{0.43291,0.92211,0.9564},{0.43783,0.92226,0.94394},{0.44275,0.92241,0.93149},{0.44767,0.92257,0.91903},{0.4526,0.92272,0.90657},{0.45752,0.92288,0.89412},{0.46244,0.92303,0.88166},{0.46736,0.92318,0.8692},{0.47228,0.92334,0.85675},{0.4772,0.92349,0.84429},{0.48212,0.92364,0.83183},{0.48704,0.9238,0.81938},{0.49196,0.92395,0.80692},{0.49689,0.92411,0.79446},{0.50181,0.92426,0.78201},{0.50673,0.92441,0.76955},{0.51165,0.92457,0.75709},
{0.51657,0.92472,0.74464},{0.52149,0.92488,0.73218},{0.52641,0.92503,0.71972},{0.53133,0.92518,0.70727},{0.53626,0.92534,0.69481},{0.54118,0.92549,0.68235},{0.55148,0.92841,0.68051},{0.56178,0.93133,0.67866},{0.57209,0.93426,0.67682},{0.58239,0.93718,0.67497},{0.5927,0.9401,0.67313},{0.603,0.94302,0.67128},{0.6133,0.94594,0.66943},{0.62361,0.94887,0.66759},{0.63391,0.95179,0.66574},{0.64421,0.95471,0.6639},{0.65452,0.95763,0.66205},{0.66482,0.96055,0.66021},{0.67512,0.96348,0.65836},{0.68543,0.9664,0.65652},{0.69573,0.96932,0.65467},{0.70604,0.97224,0.65283},{0.71634,0.97516,0.65098},{0.72664,0.97809,0.64913},{0.73695,0.98101,0.64729},{0.74725,0.98393,0.6454},{0.75755,0.98685,0.6436},{0.76786,0.98977,0.64175},{0.77816,0.9927,0.63991},{0.78847,0.99562,0.63806},{0.79877,0.99854,0.63622},{0.80661,0.99854,0.63214},{0.812,0.99562,0.62584},{0.81738,0.9927,0.61953},{0.82276,0.98977,0.61323},{0.82814,0.98685,0.60692},{0.83353,0.98393,0.60062},{0.83891,0.98101,0.59431},{0.84429,0.97809,0.588},{0.84967,0.97516,0.5817},{0.85506,0.97224,0.57539},{0.86044,0.96932,0.56909},{0.86582,0.9664,0.56278},{0.8712,0.96348,0.55648},{0.87659,0.96055,0.55017},{0.88197,0.95763,0.54387},{0.88735,0.95471,0.53756},{0.89273,0.95179,0.53126},{0.89812,0.94887,0.52495},{0.9035,0.94594,0.51865},{0.90888,0.94302,0.51234},{0.91426,0.9401,0.50604},{0.91965,0.93718,0.49973},{0.92503,0.93426,0.49343},{0.93041,0.93133,0.48712},{0.93579,0.92841,0.48082},{0.94118,0.92549,0.47451},{0.94348,0.91826,0.46928},{0.94579,0.91103,0.46405},{0.9481,0.90381,0.45882},{0.9504,0.89658,0.45359},{0.95271,0.88935,0.44837},{0.95502,0.88212,0.44314},{0.95732,0.87489,0.43791},{0.95963,0.86767,0.43268},{0.96194,0.86044,0.42745},{0.96424,0.85321,0.42222},{0.96655,0.84598,0.41699},{0.96886,0.83875,0.41176},{0.97116,0.83153,0.40654},{0.97347,0.8243,0.40131},{0.97578,0.81707,0.39608},{0.97809,0.80984,0.39085},{0.98039,0.80261,0.38562},{0.9827,0.79539,0.38039},{0.98501,0.78816,0.37516},{0.98731,0.78093,0.36993},{0.98962,0.7737,0.36471},{0.99193,0.76647,0.35948},{0.99423,0.75925,0.35425},{0.99654,0.75202,0.34902},{0.99885,0.74479,0.34379},{1.,0.73902,0.33972},{1.,0.73472,0.33679},{1.,0.73041,0.33387},{1.,0.72611,0.33095},{1.,0.7218,0.32803},{1.,0.71749,0.32511},{1.,0.71319,0.32218},{1.,0.70888,0.31926},{1.,0.70458,0.31634},{1.,0.70027,0.31342},{1.,0.69596,0.3105},{1.,0.69166,0.30757},{1.,0.68735,0.30465},{1.,0.68304,0.30173},{1.,0.67874,0.29881},{1.,0.67443,0.29589},{1.,0.67013,0.29296},{1.,0.66582,0.29004},{1.,0.66151,0.28712},{1.,0.65721,0.2842},{1.,0.6529,0.28128},{1.,0.6486,0.27835},{1.,0.64429,0.27543},{1.,0.63998,0.27251},{1.,0.63568,0.26959},{1.,0.63137,0.26667},{1.,0.63522,0.27666},{1.,0.63906,0.28666},{1.,0.64291,0.29666},{1.,0.64675,0.30665},{1.,0.6506,0.31665},{1.,0.65444,0.32664},{1.,0.65829,0.33664},{1.,0.66213,0.34664},{1.,0.66597,0.35663},{1.,0.66982,0.36663},{1.,0.67366,0.37662},{1.,0.67751,0.38662},{1.,0.68135,0.39662},{1.,0.6852,0.40661},{1.,0.68904,0.41661},{1.,0.69289,0.42661},{1.,0.69673,0.4366},{1.,0.70058,0.4466},{1.,0.70442,0.45659},{1.,0.70827,0.46659},{1.,0.71211,0.47659},{1.,0.71596,0.48658},{1.,0.7198,0.49658},{1.,0.72364,0.50657},{1.,0.72749,0.51657},{1.,0.73472,0.53095},{1.,0.74533,0.54971},{1.,0.75594,0.56847},{1.,0.76655,0.58724},{1.,0.77716,0.606},{1.,0.78777,0.62476},{1.,0.79839,0.64352},{1.,0.809,0.66228},{1.,0.81961,0.68105},{1.,0.83022,0.69981},{1.,0.84083,0.71857},{1.,0.85144,0.73733},{1.,0.86205,0.75609},{1.,0.87266,0.77486},{1.,0.88328,0.79362},{1.,0.89389,0.81238},{1.,0.9045,0.83114},{1.,0.91511,0.8499},{1.,0.92572,0.86867},{1.,0.93633,0.88743},{1.,0.94694,0.90619},{1.,0.95755,0.92495},{1.,0.96817,0.94371},{1.,0.97878,0.96248},{1.,0.98939,0.98124},{1.,1.,1.}};



// project GPS data into cartesian coordinates, and insert them into a Delaunay 2D triangulation
void proj93(Delaunay& dt, char* file_name);

// return the depth z according to a (x,y) point and the triangulation
double get_z(const Delaunay& dt, double x, double y, Face_handle& old_fh, const Vector_3& sun, double& shadow);
// if PGM: val1 is set to match the grey value | if PPM: val1, val2, val3 are respectively set to match R, G and B value
void z_to_color(const double z, int& val1, int& val2, int& val3, const double shadow, int (&rgbMax)[2]);
// writes to file in binary or ASCII
void write_val(fstream& data, const int val, const bool bin);
// reads file in binary or ASCII
int read_val(fstream& data, int& val);
// writes to file according to PPM or PGM
void write_img(fstream& data, const int val1, const int val2, const int val3, const bool bin);

// updates min_x, max_x, min_y, max_y, min_z, max_z according to the boundaries of the input terrain
void update_maxmin(double& min_x, double& max_x, double& min_y, double& max_y, double& min_z, double& max_z, double new_x, double new_y, double new_z);

// function and thread associated to create the PGM image
void create_img(const Delaunay& dt, const int& img_width, const int& img_height, const double& density, const Vector_3& sun);
void thread_img(const int k, const Delaunay& dt, const int& img_width, const int& img_height, const double& density, const Vector_3& sun, int (&rgbMax)[2]);

// brightens the image after the hill shading process
void thread_lift(const int k, const int (&rgbMax)[2]);





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

  // calculation for hill-shading
  double altitude_rad = (ALTITUDE*M_PI)/180;
  double azimut_rad   = (AZIMUT*M_PI)/180;
  Vector_3 sun( cos(azimut_rad)*cos(altitude_rad) , sin(azimut_rad)*cos(altitude_rad) , -sin(altitude_rad) );

  // PGM image creation, result will be stored in /renders
  create_img(dt, img_width, img_height, density, sun);

  // chrono end
  auto end = std::chrono::system_clock::now();
  chrono::duration<double> elapsed_seconds = end-start;

  cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";

  return EXIT_SUCCESS;
}





/**
 * \fn int proj93(Delaunay& dt, char* file_name)
 * \brief projects GPS data into cartesian coordinates, and insert them into a Delaunay 2D triangulation
 * \param dt edited Delaunay_triangulation_2 given the file file_name
 */
void proj93(Delaunay& dt, char* file_name) {
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
    c.lpz.z = abs(z);                       // useful because of z axis can be ascending or descending
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
}


/**
 * \fn double get_z(const Delaunay& dt, double x, double y, Face_handle& old_fh, const Vector_3& sun, double& shadow)
 * \brief computes the depth z of a (x,y) point and its shadow coefficient
 * \param dt Delaunay_triangulation_2, x and y coordinates of the point, old_fh keeps the face of the last point located, the sun vector is compared to the normale of the face to compute shadow
 * \return the depth z of the (x,y) point
 */
double get_z(const Delaunay& dt, double x, double y, Face_handle& old_fh, const Vector_3& sun, double& shadow){
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

  if (SHD) {
    Vector_3 n = normal(a, b, c);
    if (n.z() < 0) {n = n*(-1);}
    shadow = - scalar_product(n,sun) / (n.squared_length()*sun.squared_length());
    if (shadow < 0) {shadow = 0;}
    if (shadow > 1) {shadow = 1;}
  }

  return z;
}

/**
 * \fn void z_to_color(const double z, int& val1, int& val2, int& val3, const double shadow, int (&rgbMax)[2])
 * \brief computes the color corresponding to the depth z
 * \param given the depth z and the shadow coefficient, if PGM: val1 is set to match the grey value, if PPM: val1, val2, val3 are respectively set to match R, G and B value, additionaly rgbMax keeps the maximum value of R, G and B for post lifting
 */
void z_to_color(const double z, int& val1, int& val2, int& val3, const double shadow, int (&rgbMax)[2]) {
  int val = round( (MAX_PGM-MIN_PGM)/(max_z-min_z)*(z-min_z) + MIN_PGM );
  if (!PPM) {
    val1 = val*shadow;
    if (val1==0) {val1 = 1;}
    else if (val1>rgbMax[0]) {
      lock (foo,bar);
      rgbMax[0]=val1; rgbMax[1]=0;
      bar.unlock();
      foo.unlock();}
  }
  else {
    val1 = round(tab_pixm[val][0] * MAX_PGM * shadow);
    if (val1==0) {val1 = 1;}
    else if (val1>rgbMax[0]) {
      lock (foo,bar);
      rgbMax[0]=val1; rgbMax[1]=0;
      bar.unlock();
      foo.unlock();}

    val2 = round(tab_pixm[val][1] * MAX_PGM * shadow);
    if (val2==0) {val2 = 1;}
    else if (val1>rgbMax[0]) {
      lock (foo,bar);
      rgbMax[0]=val2; rgbMax[1]=1;
      bar.unlock();
      foo.unlock();}

    val3 = round(tab_pixm[val][2] * MAX_PGM * shadow);
    if (val3==0) {val3 = 1;}
    else if (val1>rgbMax[0]) {
      lock (foo,bar);
      rgbMax[0]=val3; rgbMax[1]=2;
      bar.unlock();
      foo.unlock();}
  }
}

/**
 * \fn void write_val(fstream& data, int val, const bool bin)
 * \brief writes a value in a stream, in binary or ASCII
 * \param writes val in data, in binary if bin=true, in ASCII if bin=false
 */
void write_val(fstream& data, int val, const bool bin) {
  if (bin) {data << reinterpret_cast<char*>(&val);}
  else {data << val << " ";}
}

/**
 * \fn void write_img(fstream& data, const int val1, const int val2, const int val3, const bool bin)
 * \brief writes color values in a stream, for PGM or PPM
 * \param writes values in data, val1 if PPM=true, val1,val2,val3 if PPM=false, in binary if bin=true, in ASCII if bin=false
 */
void write_img(fstream& data, const int val1, const int val2, const int val3, const bool bin) {
  if (!PPM) {
    write_val(data, val1, bin);
  }
  else {
    write_val(data, val1, bin);
    write_val(data, val2, bin);
    write_val(data, val3, bin);
  }
}

/**
 * \fn void update_maxmin(double& min_x, double& max_x, double& min_y, double& max_y, double& min_z, double& max_z, double new_x, double new_y, double new_z)
 * \brief updates the min and max values of x, y and z of the terrain
 * \param given new_x, new_y and new_z, updates min_x, max_x, min_y, max_y, min_z, max_z
 */
void update_maxmin(double& min_x, double& max_x, double& min_y, double& max_y, double& min_z, double& max_z, double new_x, double new_y, double new_z){
  if (new_x<min_x) {min_x = floor(new_x);}
  else if (new_x>max_x) {max_x = ceil(new_x);}
  if (new_y<min_y) {min_y = floor(new_y);}
  else if (new_y>max_y) {max_y = ceil(new_y);}
  if (new_z<min_z) {min_z = new_z;}
  else if (new_z>max_z) {max_z = new_z;}
}

/**
 * \fn void create_img(const Delaunay& dt, const int& img_width, const int& img_height, const double& density, const Vector_3& sun)
 * \brief creates the image
 * \param creates the image of the dt triangulation given the sun vector, image size is defined by img_width, img_height and density
 */
void create_img(const Delaunay& dt, const int& img_width, const int& img_height, const double& density, const Vector_3& sun) {
  // generate file_name automatically according to date and hour
  time_t rawtime;
  struct tm * ptm;
  time ( &rawtime );
  ptm = gmtime ( &rawtime );
  string file_name = "../renders/" + to_string(ptm->tm_mday) + "-" + to_string(ptm->tm_mon) + "_" + to_string(ptm->tm_hour) + "-" + to_string(ptm->tm_min) + "_render";
  if (PPM) {file_name = file_name + ".ppm";}
  else {file_name = file_name + ".pgm";}

  // writes PGM file header
  fstream img;
  if        (!BIN && !PPM)    {img.open(file_name, fstream::out);                     img<<"P2 ";}
  else if   (!BIN &&  PPM)    {img.open(file_name, fstream::out);                     img<<"P3 ";}
  else if   ( BIN && !PPM)    {img.open(file_name, fstream::out | fstream::binary);   img<<"P5 ";}
  else if   ( BIN &&  PPM)    {img.open(file_name, fstream::out | fstream::binary);   img<<"P6 ";}
  img << img_width << " "
      << img_height << " "
      << (int) MAX_PGM << endl;

  // start image creation, each thread process (1/n)th of the image, then all tmp file are concatenated
  array<thread,NB_THREAD> threads_img;
  array<thread,NB_THREAD> threads_lift;
  int rgbMax [2] = {0,0};    // {val, i} : val = max value of R G and B after hill shading | i=0 for R, i=1 for G, i=2 for B

  for (int k=0; k<NB_THREAD; k++) {   // starts the threads
    threads_img[k] = thread( [&, k] { thread_img(k, dt, img_width, img_height, density, sun, rgbMax); } ); }

  for(int k=0; k<NB_THREAD; k++) {    // stops them
    threads_img[k].join(); }

  if (SHD) {
    for (int k=0; k<NB_THREAD; k++) {   // starts the threads
      threads_lift[k] = thread( [&, k] { thread_lift(k, rgbMax); } ); }

    for(int k=0; k<NB_THREAD; k++) {    // stops them
      threads_lift[k].join(); }
  }

  for(int k=0; k<NB_THREAD; k++) {   // concatenation of all n tmp files
    fstream img_mrg;
    if (SHD) {file_name = "/tmp/MNT_tmplift" + to_string(k);}
    else {file_name = "/tmp/MNT_tmp" + to_string(k);}
    img_mrg.open(file_name, fstream::in);

    img << img_mrg.rdbuf();
    img_mrg.close();
  }

  // image is created, stream can be closed
  img.close();
}

/**
 * \fn void thread_img(const int k, const Delaunay& dt, const int& img_width, const int& img_height, const double& density, const Vector_3& sun, int (&rgbMax)[2])
 * \brief thread of create_img()
 * \param create_img() parameters are passed, computes the k th part of the dt triangulation, keeps rgbMax updated
 */
void thread_img(const int k, const Delaunay& dt, const int& img_width, const int& img_height, const double& density, const Vector_3& sun, int (&rgbMax)[2]) {
  //// thread called to create the PGM image

  // opening of a tmp file for writing
  fstream img_tmp;
  string file_name;
  file_name = "/tmp/MNT_tmp" + to_string(k);    // each file_name is different and placed in /tmp folder (deleted at each startup)
  if (BIN) {img_tmp.open(file_name, fstream::out | fstream::binary);}
  else {img_tmp.open(file_name, fstream::out);}

  // calculation of the section the k-th thread has to process
  int len = ceil(img_height / NB_THREAD);
  int i_min = k*len;
  int i_max;
  if (k==NB_THREAD-1) {i_max = img_height;}   // last thread takes the rest
  else {i_max = (k+1)*len;}

  // browse pixel and write color correspondance on the image file
  double img_z;
  Face_handle old_fh = NULL;
  int val1, val2, val3;
  double shadow = 1.0;

  for (int i=i_min; i<i_max; i++) {
    for (int j=0; j<img_width; j++) {
      // color calculation if pixel point on a valid triangle (non-infinite, non-null, area < MAX_AREA)
      try {
        img_z = get_z(dt, min_x+j/density, max_y-i/density, old_fh, sun, shadow);    // will throw exception if non-valid triangle
        z_to_color(img_z, val1, val2, val3, shadow, rgbMax);
      }
      catch (invalid_argument& e) {
        val1 = 1;
        val2 = 1;
        val3 = 1;
      }
      write_img(img_tmp, val1, val2, val3, BIN);
      }
    }

  // image is created, stream can be closed
  img_tmp.close();
}

/**
 * \fn void thread_lift(const int k, const int (&rgbMax)[2])
 * \brief lifts the color of the image after the shadowing
 * \param computes the k th part of the dt triangulation given rgbMax
 */
void thread_lift(const int k, const int (&rgbMax)[2]) {
  // opening of a tmp file for reading
  fstream img_tmp;
  string file_name;
  file_name = "/tmp/MNT_tmp" + to_string(k);    // each file_name is different and placed in /tmp folder (deleted at each startup)
  if (BIN) {img_tmp.open(file_name, fstream::in | fstream::binary);}
  else {img_tmp.open(file_name, fstream::in);}

  // opening of a tmp file for writing
  fstream img_tmplift;
  string file_namelift;
  file_namelift = "/tmp/MNT_tmplift" + to_string(k);    // each file_name is different and placed in /tmp folder (deleted at each startup)
  if (BIN) {img_tmplift.open(file_namelift, fstream::out | fstream::binary);}
  else {img_tmplift.open(file_namelift, fstream::out);}

  int val;
  double coeff = (MAX_PGM-MIN_PGM)/rgbMax[0];

  while(read_val(img_tmp, val)) {
    if (val==1) {
      write_val(img_tmplift, 1, BIN);
    }
    else {
      val = round(val*coeff);
      write_val(img_tmplift, val, BIN);
    }
  }
  img_tmp.close();
  img_tmplift.close();
}

/**
 * \fn int read_val(fstream& data, int& val)
 * \brief reads int from stream
 * \param updates val given the stream data, reads in binary or ASCII
 * \return the depth z of the (x,y) point
 */
int read_val(fstream& data, int& val) {
  if (BIN) {
    char x ;
    data.read(&x, 1);
    val = static_cast<uint8_t>(x);
    if (data.eof()) {return 0;}
    return 1;
  }
  else {
    if(data>>val) {return 1;}
    return 0;
  }
}

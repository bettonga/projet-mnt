#include <iostream>
#include <fstream>
#include <proj.h>
#include <cstdio>
#include <math.h>
#include <chrono>
#include <ctime>
#include <string>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;

typedef Delaunay::Face_handle Face_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef K::Point_3   Point_3;
typedef K::Vector_3   Vector_3;

using namespace std;

#define MAX_AREA 10.0
#define MAX_PGM 255
#define MIN_PGM 10


int proj93(Delaunay& dt, char* file_name);
void update_maxmin(double& min_x, double& max_x, double& min_y, double& max_y, double& min_z, double& max_z, double new_x, double new_y, double new_z);
double get_z(const Delaunay& dt, double x, double y, Face_handle& old_fh);
void create_pgm(const Delaunay& dt, const int& img_width, const int& img_height, const double& density);
void create_pgm_bin(const Delaunay& dt, const int& img_width, const int& img_height, const double& density);

double min_x, max_x, min_y, max_y, min_z, max_z;




int main(int argc, char *argv[]) {
  // début du chrono
  auto start = std::chrono::system_clock::now();

  // projection des coordonnées GPS, WGS84 -> xyz cartésien et triangulation Delaunay avec la librairie CGAL
  Delaunay dt;
  proj93(dt, argv[1]);

  // calculs de la taille de l'image, de la densité de pixel (pixel/m) et des coeff et offset pour pgm
  const int img_width = strtol(argv[2], NULL, 10);
  const double density = img_width / (max_x - min_x);
  const int img_height = ceil(density * (max_y - min_y));

  // création de l'image
  create_pgm_bin(dt, img_width, img_height, density);
  // create_pgm(dt, img_width, img_height, density);

  // fin du chrono
  auto end = std::chrono::system_clock::now();
  chrono::duration<double> elapsed_seconds = end-start;
  cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

  return EXIT_SUCCESS;
}




int proj93(Delaunay& dt, char* file_name) {
  fstream data;
  double lat;
  double lon;
  double z;

  data.open(file_name, fstream::in);

  // initialisation de la projection
  PJ *P;
  PJ_COORD c, c_out;
  // x_0 et y_0 première mesure, par défaut +x_0=-6.74423e+06 +y_0=4.38037e+06
  P = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
                             "+proj=longlat +datum=WGS84",
                             "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
                             NULL);
  if (P==0){return 1;}
  else {
    PJ* P_for_GIS = proj_normalize_for_visualization(PJ_DEFAULT_CTX, P);
    if( 0 == P_for_GIS ) {proj_destroy(P); return 1;}
    proj_destroy(P);
    P = P_for_GIS;

    //lecture du fichier de données
    bool first = true;
    while (data >> lat >> lon >> z){
      // projection de chaque coordonnées
      c.lpz.lam = lat;
      c.lpz.phi = lon;
      c.lpz.z = z;
      c_out = proj_trans(P, PJ_FWD, c);
      if (first == true){
        min_x = c_out.xyz.x;
        max_x = c_out.xyz.x;
        min_y = c_out.xyz.y;
        max_y = c_out.xyz.y;
        min_z = c_out.xyz.z;
        max_z = c_out.xyz.z;
        first = false;
      }
      update_maxmin(min_x, max_x, min_y, max_y, min_z, max_z, c_out.xyz.x, c_out.xyz.y, c_out.xyz.z);
      // insertion dans la triangulation
      Point_3 p(c_out.xyz.x, c_out.xyz.y, c_out.xyz.z);
      dt.insert(p);
    }
    proj_destroy(P);
  }
  data.close();
  return 0;
}

double get_z(const Delaunay& dt, double x, double y, Face_handle& old_fh){
  // renvoie la profondeur z correspondant à un point (x,y)
  Point_3 p(x,y,0.0);
  Face_handle fh = dt.locate(p, old_fh);      // old_fh est le point de départ de locate (améliore fortement la rapidité du programme)
  if (fh==nullptr || dt.is_infinite(fh)) {throw invalid_argument("");}
  old_fh = fh;
  Point_3 a = fh->vertex(0)->point();
  Point_3 b = fh->vertex(1)->point();
  Point_3 c = fh->vertex(2)->point();
  Vector_3 u(a,b);
  Vector_3 v(a,c);
  if (0.5*(cross_product(u,v).squared_length()) > MAX_AREA) {throw invalid_argument("");}
  // on projette le point (x,y) sur le plan définie par le triangle
  double denom = (c[0]-a[0])*(b[1]-a[1]) - (c[1]-a[1])*(b[0]-a[0]);
  double mu = (  (b[1]-a[1])*(x-a[0]) - (b[0]-a[0])*(y-a[1]) ) / denom;
  double la = ( -(c[1]-a[1])*(x-a[0]) + (c[0]-a[0])*(y-a[1]) ) / denom;
  double z = a[2] + la*(b[2]-a[2]) + mu*(c[2]-a[2]);
  return z;
}

void update_maxmin(double& min_x, double& max_x, double& min_y, double& max_y, double& min_z, double& max_z, double new_x, double new_y, double new_z){
  // appelée pendant la projection, stocke les valeurs maximales et minimales de x, y et z
  if (new_x<min_x) {min_x = floor(new_x);}
  else if (new_x>max_x) {max_x = ceil(new_x);}
  if (new_y<min_y) {min_y = floor(new_y);}
  else if (new_y>max_y) {max_y = ceil(new_y);}
  if (new_z<min_z) {min_z = new_z;}
  else if (new_z>max_z) {max_z = new_z;}
}

void create_pgm(const Delaunay& dt, const int& img_width, const int& img_height, const double& density) {
  const double coeff = (MAX_PGM-MIN_PGM)/(max_z-min_z);
  const double offset = MIN_PGM - coeff*min_z;
  double img_z;
  int pgm;
  Face_handle old_fh = NULL;
  fstream img;

  // génération du file_name selon la date et l'heure
  time_t rawtime;
  struct tm * ptm;
  time ( &rawtime );
  ptm = gmtime ( &rawtime );
  string file_name = "../img_output/" + to_string(ptm->tm_mday) + "-" + to_string(ptm->tm_mon) + "_" + to_string(ptm->tm_hour) + "-" + to_string(ptm->tm_min) + "_render.pgm";

  // initialisaition de l'image pgm (P2 pour pgm ASCII)
  img.open(file_name, fstream::out);
  img << "P2" << " " << img_width << " " << img_height << " " << MAX_PGM << endl;

  // parcours des pixels et correspondance pixel <-> coordonnées cartésiennes
  for (int i=0; i<img_height; i++) {
    for (int j=0; j<img_width; j++) {
      // calcul du niveau de gris si le pixel pointe sur une face acceptable (finie, non nulle, et d'aire < MAX_AREA)
      try{
        img_z = get_z(dt, min_x+j/density, max_y-i/density, old_fh);
        pgm = round(coeff*img_z+offset);
      }
      catch (invalid_argument& e) {pgm = 0;}
      img << pgm << " ";
    }
    img << endl;
  }
}

void create_pgm_bin(const Delaunay& dt, const int& img_width, const int& img_height, const double& density) {
  const double coeff = (MAX_PGM-MIN_PGM)/(max_z-min_z);
  const double offset = MIN_PGM - coeff*min_z;
  double img_z;
  int pgm;
  int pgm_null = 1;             // 0 cause skipping
  Face_handle old_fh = NULL;
  fstream img;

  // génération du file_name selon la date et l'heure
  time_t rawtime;
  struct tm * ptm;
  time ( &rawtime );
  ptm = gmtime ( &rawtime );
  string file_name = "../img_output/" + to_string(ptm->tm_mday) + "-" + to_string(ptm->tm_mon) + "_" + to_string(ptm->tm_hour) + "-" + to_string(ptm->tm_min) + "_render.pgm";

  // initialisaition de l'image pgm (P2 pour pgm ASCII)
  img.open(file_name, fstream::out | fstream::binary);
  img << "P5" << " "
      << img_width << " "
      << img_height << " "
      << (int) MAX_PGM << endl;


  // parcours des pixels et correspondance pixel <-> coordonnées cartésiennes
  for (int i=0; i<img_height; i++) {
    for (int j=0; j<img_width; j++) {
      // calcul du niveau de gris si le pixel pointe sur une face acceptable (finie, non nulle, et d'aire < MAX_AREA)
      try{
        img_z = get_z(dt, min_x+j/density, max_y-i/density, old_fh);
        pgm = round(coeff*img_z+offset);
        img << reinterpret_cast<char*>(&pgm);
      }
      catch (invalid_argument& e) {img << reinterpret_cast<char*>(&pgm_null);}
    }
    //img << endl;
  }
}

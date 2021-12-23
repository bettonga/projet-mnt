#include <iostream>
#include <fstream>
#include <proj.h>
#include <cstdio>
#include <math.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;

typedef Delaunay::Face_handle Face_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef K::Point_3   Point_3;

using namespace std;



int proj93(Delaunay& dt);
void update_maxmin(double& min_x, double& max_x, double& min_y, double& max_y, double& min_z, double& max_z, double new_x, double new_y, double new_z);
double get_z(Delaunay& dt, double x, double y, Face_handle& old_fh);

double min_x, max_x, min_y, max_y, min_z, max_z;

int main() {
  int img_width = 800;
  int img_height;
  int pgm_max = 256;
  double density, coeff, offset;
  double img_z;
  int pgm;

  // projection des coordonnées GPS, WGS84 -> xyz cartésien et triangulation Delaunay avec la librairie CGAL
  Delaunay dt;
  proj93(dt);

  // calculs de la taille de l'image, de la densité de pixel (pixel/m) et des coeff et offset pour pgm
  density = img_width / (max_x - min_x);
  img_height = ceil(density * (max_y - min_y));
  coeff = pgm_max/(max_z-min_z);
  offset = -(min_z*pgm_max)/(max_z-min_z);

  // création de l'image
  Face_handle old_fh = NULL;
  fstream img;
  img.open("../img_output/yoyo.pgm", fstream::out);
  img << "P2" << endl;
  img << img_width << " " << img_height << endl;
  img << pgm_max << endl;
  for (int i=0; i<img_height; i++) {
    for (int j=0; j<img_width; j++) {
      try{
        img_z = get_z(dt, min_x+j/density, max_y-i/density, old_fh);
        pgm = round(coeff*img_z+offset);
      }
      catch (invalid_argument& e){pgm = 0;}
      img << pgm << " ";
    }
    img << endl;
  }

  cout << min_x << " " << max_x << " " << min_y << " " << max_y << " " << min_z << " " << max_z << endl;

  return EXIT_SUCCESS;

}



int proj93(Delaunay& dt) {
  fstream data;
  double lat;
  double lon;
  double z;

  data.open("../datas/rade_1m_IM.txt" , fstream::in);

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

double get_z(Delaunay& dt, double x, double y, Face_handle& old_fh){
  Point_3 p(x,y,0.0);
  Face_handle fh = dt.locate(p, old_fh);
  if (fh==nullptr || dt.is_infinite(fh)) {throw invalid_argument("Infinite or Null or  face");}
  old_fh = fh;
  Point_3 a = fh->vertex(0)->point();
  Point_3 b = fh->vertex(1)->point();
  Point_3 c = fh->vertex(2)->point();
  double denom = (c[0]-a[0])*(b[1]-a[1]) - (c[1]-a[1])*(b[0]-a[0]);
  double mu = (  (b[1]-a[1])*(x-a[0]) - (b[0]-a[0])*(y-a[1]) ) / denom;
  double la = ( -(c[1]-a[1])*(x-a[0]) + (c[0]-a[0])*(y-a[1]) ) / denom;
  double z = a[2] + la*(b[2]-a[2]) + mu*(c[2]-a[2]);
  return z;
}

void update_maxmin(double& min_x, double& max_x, double& min_y, double& max_y, double& min_z, double& max_z, double new_x, double new_y, double new_z){
  if (new_x<min_x) {min_x = floor(new_x);}
  else if (new_x>max_x) {max_x = ceil(new_x);}
  if (new_y<min_y) {min_y = floor(new_y);}
  else if (new_y>max_y) {max_y = ceil(new_y);}
  else if (new_z<min_z) {min_z = new_z;}
  if (new_z>max_z) {max_z = new_z;}
}

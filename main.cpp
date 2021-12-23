#include <iostream>
#include <fstream>
#include <proj.h>
#include <cstdio>

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



int proj93(fstream& output);

double get_z(Delaunay& dt, double x, double y, Face_handle& old_fh);

int main() {
  // projection des coordonnées GPS, WGS84 -> xyz cartésien
  fstream data_proj;

  data_proj.precision(11);    // precision of raw datas
  data_proj.open("../datas/MNT_lilrade_PROJ.txt" , fstream::out);
  proj93(data_proj);
  data_proj.close();

  ifstream in("../datas/MNT_lilrade_PROJ.txt");
  istream_iterator<Point_3> begin(in);
  istream_iterator<Point_3> end;

  Delaunay dt(begin, end);
  cout << dt.number_of_vertices() << endl;

  Face_handle old_fh = NULL;
  double z = get_z(dt, 15.0, 2.0, old_fh);
  cout << z << endl;

  return EXIT_SUCCESS;

}



int proj93(fstream& output) {
  fstream data;
  double lat;
  double lon;
  double z;

  data.open("../datas/MNT_lilrade.txt" , fstream::in);

  // initialisation de la projection
  PJ *P;
  PJ_COORD c, c_out;
  // x_0 et y_0 première mesure, par défaut: +x_0=700000 +y_0=6600000
  P = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
                             "+proj=longlat +datum=WGS84",
                             "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=-6744040 +y_0=4381740 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
                             NULL);
  if (P==0){return 1;}
  else {
    PJ* P_for_GIS = proj_normalize_for_visualization(PJ_DEFAULT_CTX, P);
    if( 0 == P_for_GIS ) {proj_destroy(P); return 1;}
    proj_destroy(P);
    P = P_for_GIS;

    //lecture du fichier de données
    while (data >> lat >> lon >> z){
      // projection de chaque coordonnées
      c.lpz.lam = lat;
      c.lpz.phi = lon;
      c.lpz.z = z;
      c_out = proj_trans(P, PJ_FWD, c);
      // écriture du fichier output
      output << c_out.xyz.x << " " << c_out.xyz.y << " " << c_out.xyz.z << endl;
    }
    proj_destroy(P);
  }
  data.close();
  return 0;
}

double get_z(Delaunay& dt, double x, double y, Face_handle& old_fh){
  Point_3 p(x,y,0.0);
  Face_handle fh = dt.locate(p, old_fh);
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

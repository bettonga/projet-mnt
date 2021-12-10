#include <iostream>
#include <fstream>
#include <proj.h>
#include

using namespace std;

int proj93(fstream& output);

int main() {
  // projection des coordonnées GPS, WGS84 -> xyz cartésien
  fstream data_proj;
  data_proj.precision(11);    // precision of raw datas
  data_proj.open("../datas/MNT_lilrade_PROJ.txt" , fstream::out);
  proj93(data_proj);




  data_proj.close();
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

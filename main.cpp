#include <iostream>
#include <fstream>
#include <stdio.h>
#include <proj.h>
#include <string.h>

using namespace std;

int proj93(fstream& output);

int main() {
  fstream data_proj;
  data_proj.open("../projet-mnt/datas/MNT_lilrade_PROJ.txt" , fstream::out);
  proj93(data_proj);
  data_proj.close();
}



int proj93(fstream& output) {
  fstream data;
  char line [100];
  double lat;
  double lon;
  double z;

  data.open("../projet-mnt/datas/MNT_lilrade.txt" , fstream::in);

  // initialisation de la projection
  PJ *P;
  PJ_COORD c, c_out;
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
  }
  for (int i=0; i<5; i++){
    data >> lat >> lon >> z;
    std::cout << lat << " " << lon << " " << z << std::endl;
  }


  // projection de chaque coordonnées
  // while ( fgets (line , sizeof(line) , data) != NULL ) {       // lit la ligne du fichier de départ
  //   char * lineToken = strtok(line, " ");
  //   c.lpz.lam = atof(lineToken.c_str);
  //   lineToken = strtok(NULL, " ");
  //   c.lpz.phi = atof(lineToken.c_str);
  //   lineToken = strtok(NULL, " ");
  //   c.lpz.z = atof(lineToken.c_str);
  //   c_out = proj_trans(P, PJ_FWD, c);
  //   printf("%.2f\t%.2f\n", c_out.xy.x, c_out.xy.y);
  //   fprintf(output, "% % %", c_out.xyz.x, c_out.xyz.y, c_out.xyz.z, 1);                                  // écrit les coordonnées projetées dans le fichier
  // }

  proj_destroy(P);
  data.close();
  return 0;
}

#include <iostream>
#include <stdio.h>
#include <proj.h>
#include <string.h>

int proj93(FILE * output);

PJ init_proj();

int main() {
  FILE * data_proj;
  data_proj = fopen("../projet-mnt/datas/MNT_lilrade_PROJ.txt" , "w");
  proj93(data_proj);
  fclose(data_proj);
}

PJ* init_proj() {
  PJ *P;
  P = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
                             "+proj=longlat +datum=WGS84",
                             "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
                             NULL);
  if (P==0){return 1;}

  {
      PJ* P_for_GIS = proj_normalize_for_visualization(PJ_DEFAULT_CTX, P);
      if( 0 == P_for_GIS )  {
          proj_destroy(P);
          return 1;
      }
      proj_destroy(P);
      P = P_for_GIS;
  }
  return *P
}

int proj93(FILE * output) {
  FILE * data;
  char line [100];

  data = fopen("../projet-mnt/datas/MNT_lilrade.txt" , "r");
  if (data == NULL) perror ("Error opening file");
  else {
    PJ_COORD c, c_out;
    while ( fgets (line , sizeof(line) , data) != NULL ) {       // lit la ligne du fichier de départ
      char * lineToken = strtok(line, " ");
      c.lpz.lam = lineToken;
      lineToken = strtok(NULL, " ");
      c.lpz.phi = lineToken;
      lineToken = strtok(NULL, " ");
      c.lpz.z = lineToken;
      c_out = proj_trans(P, PJ_FWD, c);
      printf("%.2f\t%.2f\n", c_out.xy.x, c_out.xy.y);
      fprintf(output, "% % %", c_out,xyz.x, c_out.xyz.y, c_out.xyz.z, 1);                                  // écrit les coordonnées projetées dans le fichier
      proj_destroy(P);
    }
    fclose (data);
  }
  return 0;
}

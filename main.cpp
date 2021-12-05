#include <iostream>
#include <stdio.h>
#include <proj.h>
#include <string.h>

int proj93(FILE * output);

int main() {
  FILE * data_proj;
  data_proj = fopen("../projet-mnt/datas/MNT_lilrade_PROJ.txt" , "w");
  proj93(data_proj);
  fclose(data_proj);
}


int proj93(FILE * output) {
  FILE * data;
  char line [100];

  data = fopen("../projet-mnt/datas/MNT_lilrade.txt" , "r");
  if (data == NULL) perror ("Error opening file");
  else {
    while ( fgets (line , sizeof(line) , data) != NULL ) {       // lit la ligne du fichier de départ
      fprintf(output, line, 1);                                  // écrit les coordonnées projetées dans le fichier
    }
    fclose (data);
  }
  return 0;
}

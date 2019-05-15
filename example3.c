#include <stdio.h>
#include "IERS.h"
#include <stdlib.h>

double **eopdata = NULL;
size_t eopsize = 0;
int main()
{
    FILE *fp = fopen("eop19620101.txt","r");
    if(!fp)
    {
        fprintf(stderr, "ERROR: No se ha podido abrir eop19620101.txt\n");
        exit(1);
    }

    fseek(fp, 0, SEEK_END);
    long filesize = ftell(fp);
    if(filesize == -1)
    {
        fprintf(stderr, "ERROR: Archivo inválido.\n");
        fclose(fp);
        exit(2);
    }

    rewind(fp);
    eopsize = (size_t) filesize;
    eopdata = malloc(eopsize * sizeof(double*));

    char line[128];
    int c1 = 0, c2 = 0, c3 = 0, c4 = 0, c13 = 0;
    float c5 = 0.F, c6 =0.F, c7= 0.F, c8 =0.F, c9=0.F, c10 = 0.F, c11 = 0.F, c12 = 0.F;
    for(size_t i=0; i<eopsize; i++)
    {
        eopdata[i] = malloc(13 * sizeof(double));
        fgets(line, 255, fp);
        sscanf(line, "%d %d %d %d %f %f %f %f %f %f %f %f %d ", &c1,  &c2,  &c3,  &c4,  &c5, &c6,
                 &c7,  &c8,  &c9,  &c10,  &c11,  &c12,  &c13);

        eopdata[i][0] = c1;
        eopdata[i][1] = c2;
        eopdata[i][2] = c3;
        eopdata[i][3] = c4;
        eopdata[i][4] = (double) c5;
        eopdata[i][5] = (double) c6;
        eopdata[i][6] = (double) c7;
        eopdata[i][7] = (double) c8;
        eopdata[i][8] = (double) c9;
        eopdata[i][9] = (double) c10;
        eopdata[i][10] = (double) c11;
        eopdata[i][11] = (double) c12;
        eopdata[i][12] = c13;
    }

    fclose(fp);

    fp = fopen("sat3.txt","r");
    if(!fp)
    {
        fprintf(stderr, "ERROR: No se ha podido abrir sat1.txt\n");
        exit(1);
    }
    int Y = 0;
    int M = 0;
    int D = 0;
    int h = 0;
    int m = 0;
    float s = 0.F;
    float rtasc = 0.F;
    float decl = 0.F;

    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    rewind(fp);

    double **obs = calloc(3, sizeof(double*));
    for(long i=0; i<fsize; i++)
    {
        fgets(line, 128, fp);
        sscanf(line, "%d/%d/%d %d:%d:%f %f %f", &Y, &M, &D, &h, &m, &s, &rtasc, &decl);
    }

}

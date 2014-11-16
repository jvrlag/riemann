#ifndef GRAPHICS_H
#define GRAPHICS_H

#include"easyx.h"
#include"filam.h"
//#include"measures_conos.h"

class Graphics
{
public:
     int xsize, ysize;
     palette amarillo, blanco;
     bool show_points;
     double fx, fy, yav;

     Graphics();
     ~Graphics() {};
     void Start(int xs, int ys);
     void Show(const Filament &F);
     void Show_Segment(point *P1, point *P2);
     void Mark_Point(point *P);
     void Pixel(double x, double y);
     void Line(double x0, double y0, double x1, double y1);
};

#endif

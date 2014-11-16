#include"graphics.h"

Graphics::Graphics() { }

void Graphics::Start(int xs, int ys)
{
     xsize=xs; ysize=ys;
     EXStart(100,100,xsize,ysize);
     EXEnableBuffer();
     fx=10.0, fy=10.0;
     yav=0.0;
     amarillo=EXAllocNamedColor("yellow");
     blanco=EXAllocNamedColor("white");
     show_points=false;
}

void Graphics::Show(const Filament &F)
{
     point *P=F.root;
     for (long i=1;i<=F.N;i++)
     {
	  Show_Segment(P,P->next);
	  if (show_points) Mark_Point(P);
	  P=P->next;
     }
}

// it shows just a segment
void Graphics::Show_Segment(point *P1, point *P2)
{
     double x1=P1->x, y1=P1->y;
     double x2=P2->x, y2=P2->y;
     int ix=xsize/2+(int)floor(fx*x1);
     int iy=ysize/2-(int)floor(fy*y1);
     int ix2=xsize/2+(int)floor(fx*x2);
     int iy2=ysize/2-(int)floor(fy*y2);
     EXLine(ix,iy,ix2,iy2);
}

void Graphics::Mark_Point(point *P)
{
     int ix=xsize/2+(int)floor(fx*P->x);
     int iy=ysize/2-(int)floor(fy*P->y);
     EXFillCircle(ix,iy,2);
}

void Graphics::Pixel(double x, double y)
{
     int ix=xsize/2+(int)floor(fx*x);
     int iy=ysize/2-(int)floor(fy*y);
//     printf("ix: %d, iy: %d\n",ix,iy);
     EXPixel(ix,iy);
}

void Graphics::Line(double x0, double y0, double x1, double y1)
{
     int ix0=xsize/2+(int)floor(fx*x0);
     int iy0=ysize/2-(int)floor(fy*y0);
     int ix1=xsize/2+(int)floor(fx*x1);
     int iy1=ysize/2-(int)floor(fy*y1);
     EXLine(ix0,iy0,ix1,iy1);
}

#ifndef FILAM_H
#define FILAM_H

#include"geom.h"

extern double Lx;
extern double Ly;

typedef struct 
{
     bool crumb; // Whether the cut is a miguita (crumb) or a bubble
     double length; // How much length was removed
}Cut_Info;


class Filament
{
 public:
     long N;       // number of points
     point *root;  // entry
     double maxl;    // distances should be in the interval [minl, maxl]
     double minl; 
     long Ncuts;    // Number of times we've had to cut
     double length_cut; // total length we've cut so far
//     long overrel; // Number of "over-relaxations"
     Metric &G;

     Filament(Metric &g);
     ~Filament();
     
     void Destroy();
     
     point* Add_Point(point*); // add a P betw a given P and its successor
     point* Add_Point(point*,double,double); 
     void Remove_Point(point*);
     void Cut(point *P1, point *P2);

     void Create_Backbone(long); 

     void Check_Resolution();
     void Check_Consistency();
     void Check_Selfintersect();
     bool Remove_Selfintersect(Cut_Info &);
     void Update(); // update all points in the filament

     void Write() const;
     void Write(char *filename) const;

     double curv_distance(point *P, point *Q) const;

     void Read(char *filename);
};


void filament_random_circle(Filament &F, double R0, long n);
void filament_circle(Filament &F, double R0, long n);
void filament_otra_cerrada_suave(Filament &F, double R0, double A, double k, long n);
void Copy(Filament &F2, const Filament &F);

#endif

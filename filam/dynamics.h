// 080823
#ifndef DYNAMICS_HEADER
#define DYNAMICS_HEADER
#include"filam.h"
#include"geom.h"

class Dynamics
{
public:
     double Dt;
     Metric &G;
     Dynamics(Metric &g) : G(g) {};
     ~Dynamics() {};
};

void dynamics(point *P, Dynamics &D);
void evolstep(Filament &F, Dynamics &D);

#endif

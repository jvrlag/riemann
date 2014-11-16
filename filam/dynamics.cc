// 080823-140428
#include"dynamics.h"

// execute the dynamics on point P
void dynamics(point *P, Dynamics &D)
{
     Vector vn=normal_vector(P, D.G);
     P->x2+=D.Dt*vn(1);
     P->y2+=D.Dt*vn(2);
}

void evolstep(Filament &F, Dynamics &D)
{
     point *P=F.root;
     do
     {
	  dynamics(P,D);
	  P=P->next;
     }while(P!=F.root);
     F.Update();
}




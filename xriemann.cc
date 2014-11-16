// Program to compute a geodesic in any metric and show it
// 140509
#include"geom.h"
#include"graphics.h"
#include"filam.h"
#include"dynamics.h"


Matrix cosines(double x, double y, const Vector &D)
{
     Matrix g(2);
     g(1,1)=0.5+5.0*sqr(cos(x));
     g(2,2)=0.5+8.0*sqr(cos(y));
     return g;
}

Matrix eggcrate(double x, double y, const Vector &D)
{
     Matrix g(2);
     g(1,1)=1.0+sqr(cos(x));
     g(2,2)=1.0+sqr(cos(y));
     g(1,2)=cos(x)*cos(y);
     g(2,1)=cos(x)*cos(y);
     return g;
}

Matrix paraboloid(double x, double y, const Vector &D)
{
     Matrix g(2);
     g(1,1)=1.0+4.0*x*x;
     g(2,2)=1.0+4.0*y*y;
     g(1,2)=4.0*x*y;
     g(2,1)=4.0*x*y;
     return g;
}

Matrix saddle(double x, double y, const Vector &D)
{
     Matrix g(2);
     g(1,1)=1.0+4.0*x*x;
     g(2,2)=1.0+4.0*y*y;
     g(1,2)=-4.0*x*y;
     g(2,1)=-4.0*x*y;
     return g;
}

Matrix euclidean(double x, double y, const Vector &D)
{
     return 2.0*Unit(2);
}


int main()
{
     Metric metrica(cosines);
     double dt=1e-2;
     long ntimes=1600;

     Graphics G;
     G.Start(800,800);
     G.fx=G.fy=40.0;

     Dynamics D(metrica);
     D.Dt=dt;

     
     for (double alpha=0.0;alpha<=1.0;alpha+=0.01)
     {
	  Vector X(2), V(2);
	  X(1)=0.0; X(2)=0.0; 
	  V(1)=cos(2.0*M_PI*alpha);
	  V(2)=sin(2.0*M_PI*alpha);
	  // Normalize according to the metric!!
	  Matrix g=metrica(X(1),X(2));
	  double norma=sqrt(g.Elem(V,V));
	  V*=1.0/norma;
	  
	  for (long it=1;it<=ntimes;it++)
	  {
	       point P, P2;
	       P.x=X(1);
	       P.y=X(2);
	       P2.x=X(1)+dt*V(1);
	       P2.y=X(2)+dt*V(2);
	       Vector Vp=parallel_transport(V,&P,&P2,metrica);
	       X=X+dt*V;
	       V=Vp;
	       double factor=(double)it/(double)ntimes;
	       EXColor(1.-factor,factor,0);
	       G.Pixel(X(1),X(2));
	  }
	  EXFlush();
     }


     Filament F(metrica);
     filament_circle(F,0.01,100);
     F.maxl=0.1;
     F.minl=0.01;

     for (long it=1;it<=ntimes;it++)
     {
	  evolstep(F,D);
	  if (it%20==0)
	  {
	       F.Check_Selfintersect();
	       F.Check_Resolution();
	  }
	    
	  if (it%200==0) 
	  {
	       double factor=(double)it/(double)ntimes;
	       EXColor(1.-factor,factor,0.0);

	       F.Check_Selfintersect();
	       F.Check_Resolution();
	       G.Show(F);
	       EXFlush();
	  }
     }


     EXReadKey();

     
}

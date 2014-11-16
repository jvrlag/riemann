// Geometry file for filament structures
// It computes distances between points, angles, curvatures
// Fits a set of points to a straight line (diagonally) and
// finds the deviation of the points to the line...
// NO PBC, 110506
// METRIC, 140428
#ifndef GEOMETRY
#define GEOMETRY

#include"geom.h"

void write(point *P)
{
     Vector V=get_vector(P);
     V.Write();
}

Vector get_vector(point* P)
{
     Vector V(2);
     V(1)=P->x; V(2)=P->y;
     return V;
}

// Update point: x<-x2, y<-y2 
void update(point *P)
{
     P->x=P->x2; P->y=P->y2;
}

// Returns P->Q
Vector get_vector(point *P, point *Q, bool future)
{
     Vector U(2);
     if (future)
     {
	  U(2)=(Q->y2)-(P->y2);
	  U(1)=(Q->x2)-(P->x2);
     }
     else
     {
	  U(2)=(Q->y)-(P->y);
	  U(1)=(Q->x)-(P->x);
     }
     return U;
}

/////////////////////////
// METRIC FUNCTIONS
/////////////////////////

double distance(point *c1, point *c2, const Metric &G)
{
     Vector v=get_vector(c1,c2);
     Vector x=mid_point(c1,c2);
     Matrix g=G(x(1),x(2));
     return sqrt(g.Elem(v,v));
}

double distance(point *c1, point *c2)
{
     Vector v=get_vector(c1,c2);
     return v.Norm();
}

double cos_angle(point *P, const Metric &G)
{
     point *Pb=P->prev;
     point *Pf=P->next;
     Matrix g=G(P->x,P->y);
     Vector u=get_vector(Pb,P);
     Vector v=get_vector(P,Pf);
     double modu=g.Elem(u,u);
     double modv=g.Elem(v,v);
     double dotprod=g.Elem(u,v);
     return dotprod/sqrt(modu*modv);
}

double cos_angle(point *P)
{
     point *Pb=P->prev;
     point *Pf=P->next;
     Vector u=get_vector(Pb,P);
     Vector v=get_vector(P,Pf);
     double modu=u.Norm();
     double modv=v.Norm();
     double dotprod=Dot(u,v);
     return dotprod/(modu*modv);
}

// Formula for the curvature: 2 sin(x_1-x_2-x_3)/|x_1-x_3|
// If middle point is left of the tangent line, then curvature is +
double curv(point *P)
{
      point *Pb=P->prev;
      point *Pf=P->next;
      Vector V1=get_vector(Pb,P);
      Vector V2=get_vector(Pf,P);
      Vector V3=get_vector(Pf,Pb);
      double crossp=det2(V1(1),V1(2),V2(1),V2(2));
      double K=2.0*crossp/(V1.Norm()*V2.Norm()*V3.Norm());
      return K;
}

double turning_angle(point *P)
{
     double sina=sin_angle(P);
     double angle=asin(sina);
     double cosa=cos_angle(P);
     if (sina>0.0 && cosa<0.0) angle=M_PI-angle;
     if (sina<0.0 && cosa<0.0) angle=-(M_PI+angle);
     return angle;
}

double sin_angle(point *P)
{
      point *Pb=P->prev;
      point *Pf=P->next;
      Vector V1=get_vector(Pb,P);
      Vector V2=get_vector(P,Pf);
      double crossp=det2(V1(1),V1(2),V2(1),V2(2));
      double sina=crossp/(V1.Norm()*V2.Norm());
      return sina;
}

// Estimate the coordinates of the normal vector at the point
Vector normal_vector(point *P)
{
     Vector v=get_vector(P->prev,P->next);
     Vector vn(2);
     vn(1)=v(2);
     vn(2)=-v(1);
     vn.Normalize();
     return vn;
}

// Vector normal_vector(point *P, const Metric &G)
// {
//      Vector v=get_vector(P->prev,P->next);
//      Vector vn(2);
//      vn(1)=v(2);
//      vn(2)=-v(1);
// //     vn.Normalize();
//      Matrix g=G(P->x,P->y);
//      // double det=g.Det();
//      // if (fabs(det)<1e-10) merror("Metric is singular\n");
//      g.Solve(vn);
//      double norma=sqrt(g.Elem(vn,vn));
//      vn*=1.0/norma;
//      return vn;
// }

Vector normal_vector(point *P, const Metric &G)
{
     Vector t=get_vector(P->prev,P->next);
     Matrix g=G(P->x,P->y);
     Vector tg(2);
     tg(1)=g(1,1)*t(1)+g(1,2)*t(2);
     tg(2)=g(2,1)*t(1)+g(2,2)*t(2);
     
     Vector ng(2);
     ng(1)=tg(2);
     ng(2)=-tg(1);
          
     double norma2=g(1,1)*sqr(ng(1))+g(2,2)*sqr(ng(2))+2.0*g(1,2)*ng(1)*ng(2);
     
     return (1.0/sqrt(norma2)*ng);
}

Vector parallel_transport(const Vector &V, point *P1, point *P2, 
			  const Metric &G)
{
     Matrix Gp=G(P1->x,P1->y); // métrica en el punto 1
     
     Matrix IGp=Invert(Gp);  //inversa de la métrica en el punto 1
     
     //derivadas de los elementos de la métrica
     double h=1e-9;
     //respecto x
     Matrix Gpx2=G(P1->x+h,P1->y);
     Matrix dGpx=Gpx2-Gp;
     dGpx*=(1.0/h);
     //respecto y 
     Matrix Gpy2=G(P1->x,P1->y+h);
     Matrix dGpy=Gpy2-Gp;
     dGpy*=(1.0/h);
     
     //los simbolos
     Matrix Ch1(2);  // Gamma^1_{k,l}
     
     Ch1(1,1)=(IGp(1,1)*dGpx(1,1)+IGp(1,2)*(2.0*dGpx(2,1)-dGpy(1,1)))/2.0;
     Ch1(1,2)=(IGp(1,1)*dGpy(1,1)+IGp(1,2)*dGpx(2,2))/2.0;
     Ch1(2,1)=Ch1(1,2);
     Ch1(2,2)=(IGp(1,1)*(2.0*dGpy(1,2)-dGpx(2,2))+IGp(1,2)*dGpy(2,2))/2.0;
     
     Matrix Ch2(2);  // Gamma^2_{i,j}
     Ch2(1,1)=(IGp(2,1)*dGpx(1,1)+IGp(2,2)*(2.0*dGpx(2,1)-dGpy(1,1)))/2.0;
     Ch2(1,2)=(IGp(2,1)*dGpy(1,1)+IGp(2,2)*dGpx(2,2))/2.0;
     Ch2(2,1)=Ch2(1,2);
     Ch2(2,2)=(IGp(2,1)*(2.0*dGpy(1,2)-dGpx(2,2))+IGp(2,2)*dGpy(2,2))/2.0;

     // el transporte: Vt^i = V^i - Gamma^i_{kl} V^k V^l
     Vector Vd=get_vector(P1,P2);
     Vector Vt(2);
     
     Vt(1)=V(1)-Ch1(1,1)*V(1)*Vd(1)-Ch1(1,2)*V(1)*Vd(2)
	  -Ch1(2,1)*V(2)*Vd(1)-Ch1(2,2)*V(2)*Vd(2);
     Vt(2)=V(2)-Ch2(1,1)*V(1)*Vd(1)-Ch2(1,2)*V(1)*Vd(2)
	  -Ch2(2,1)*V(2)*Vd(1)-Ch2(2,2)*V(2)*Vd(2);
     return Vt;
}

double geodesic_curvature(point *P, const Metric &G)
{
      point *Pb=P->prev;
      point *Pf=P->next;
      Vector V1=get_vector(Pb,P);
      Vector V2=get_vector(P,Pf);
      Vector V1t=parallel_transport(V1,Pb,P,G);
      Matrix g=G(P->x,P->y);
      double length=sqrt(g.Elem(V1t,V1t));
      double cosa=g.Elem(V2,V1t)/sqrt(g.Elem(V2,V2)*g.Elem(V1t,V1t));
      double alpha;
      if (fabs(cosa)>=1.0) 
      {
	   if (cosa>0.0) alpha=0.0; 
	   else alpha=M_PI;
      }
      else alpha=acos(cosa);
      double crossprod=det2(V1t(1),V1t(2),V2(1),V2(2));
      return SIGN(1,crossprod)*alpha/length;
}

double geodesic_curvature_solo_angulo(point *P, const Metric &G)
{
     point *Pb=P->prev;
     point *Pf=P->next;
     Vector V1=get_vector(Pb,P);
     Vector V2=get_vector(P,Pf);
     Vector V1t=parallel_transport(V1,Pb,P,G);
     Matrix g=G(P->x,P->y);
     double cosa=g.Elem(V2,V1t)/sqrt(g.Elem(V2,V2)*g.Elem(V1t,V1t));
     if (cosa>1.0) cosa=1.0;
     double alpha=acos(cosa);
     double crossprod=det2(V1t(1),V1t(2),V2(1),V2(2));
     return SIGN(1,crossprod)*alpha;
}

///////////////////////////////////
// AFFINE FUNCTIONS
///////////////////////////////////

// return the mid point... affine, it's just an approximation!
Vector mid_point(point *p1, point *p2)
{
     Vector V(2);
     V(1)=.5*(p1->x+p2->x);
     V(2)=.5*(p1->y+p2->y);
     return V;
}

// returns +1 if point c2 is at the left side of the c1->c3 line
bool left_q(point *c1, point *c2, point *c3)
{
     Vector v=get_vector(c2,c3);
     Vector u=get_vector(c2,c1);
     return (det2(v(1),v(2),u(1),u(2))<0.0);     
}

// checks whether the segments intersect
// perhaps there is a more elegant way, but...
bool intersect_q(point *P1, point *P2, point *Q1, point *Q2)
{
     Vector p1=get_vector(P1);
     Vector p2=get_vector(P2);
     Vector q1=get_vector(Q1);
     Vector q2=get_vector(Q2);
     return intersect_q(p1,p2,q1,q2);
}

// euclidean checker
bool intersect_q(Vector &P1, Vector &P2, Vector &Q1, Vector &Q2)
{
     if (!overlap(P1(1),P2(1),Q1(1),Q2(1))) return false;
     if (!overlap(P1(2),P2(2),Q1(2),Q2(2))) return false;

     double detA=det2(P2(1)-P1(1), Q2(1)-Q1(1), P2(2)-P1(2), Q2(2)-Q1(2));
     double L=det2(Q1(1)-P1(1), Q2(1)-Q1(1), Q1(2)-P1(2), Q2(2)-Q1(2))/detA;
     double M=-det2(P2(1)-P1(1), Q1(1)-P1(1), P2(2)-P1(2), Q1(2)-P1(2))/detA;
     
     if (L>=0.0 && L<=1.0 && M>=0.0 && M<=1.0) return true;
     else return false;
}

/////////////////////////////
// AUXILIARY ROUTINES
/////////////////////////////

double det2(double a11, double a12, double a21, double a22)
{
     return a11*a22-a12*a21;
}

// return true if x is in [a,b], even if b<a
bool is_between(double x, double a, double b)
{
     if (a<b) return ( (a<x) && (x<b) );
     else return ( (b<x) && (x<a) );
}

// return true if there is overlap between [x1,x2] and [y1,y2]...
// Caution! they need not be in order!
bool overlap(double x1, double x2, double y1, double y2)
{
     return (is_between(x1,y1,y2) || is_between(x2,y1,y2) ||
	     is_between(y1,x1,x2) || is_between(y2,x2,x2) );
}

// void make_straight_line(point *P, bool nuevo)
// {
//      Vector u=get_vector(P->prev,P->next);
//      P->X=P->prev->x+0.5*u(1);
//      P->Y=P->prev->y+0.5*u(2);
// }

// Finds out the distance from point P to the tangent line P->prev, P->next
// with sign!! 
// double distance_to_baseline(point *P, bool nuevo)
// {
//      Vector u=get_vector(P->prev,P->next,nuevo);
//      Vector v=get_vector(P->prev,P,nuevo);
//      double crossprod=cross_prod(u(1),u(2),v(1),v(2));
//      return crossprod/u.Norm();
// }

#endif

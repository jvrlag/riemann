// Geometry file header
// 040525-071130-140428
#ifndef GEOMETRY_HEADER
#define GEOMETRY_HEADER
#include"common.h"
#include"matrix.h"

////////////////////////////////////////////////////////////////
// Point struct and referred functions
////////////////////////////////////////////////////////////////
typedef struct point
{
     double x, y;      // current position
     double x2, y2;   // future position, being computed right now
     point *next;    // next point in filament
     point *prev;    // previous point in filament
} point;
// write
void write(point *);
// Get the coordinates of the given point
Vector get_vector(point*);
// Get the vector joining the two points (if future, then use future components)
Vector get_vector(point*, point*, bool future=false);
// (x2,y2) -> (x,y)
void update(point *P);


// AFFINE FUNCTIONS: DO NOT NEED METRIC
// Find whether B is at the left or right of A->C segment
bool left_q(point *A, point *B, point *C);
// wrapper for euclidean checker
bool intersect_q(point *P1, point *P2, point *Q1, point *Q2);
// euclidean checker
bool intersect_q(Vector &P1, Vector &P2, Vector &Q1, Vector &Q2);
// mid point
Vector mid_point(point *, point*);


//////////////////////////////////////////////////////////////////
// Metric structure
//////////////////////////////////////////////////////////////////

typedef Matrix (*Metric_F)(double x, double y, const Vector &param);

class Metric
{
public:
     Metric_F G;
     Vector Data; // required parameters
     Metric(Metric_F Gq) { G=Gq; Data.Start(); }
     Metric(Metric_F Gq, const Vector &D) { G=Gq; Data=D; };
     ~Metric() {};
     Matrix operator()(double x, double y) const { return G(x,y,Data); }
};

// METRIC FUNCTIONS
// Estimate normal vector to Filament at point P
Vector normal_vector(point *P); // EUCLIDEAN
Vector normal_vector(point *P, const Metric &G);
// Estimate distance between two (near) points
double distance(point *c1, point *c2, const Metric &G);
double distance(point *c1, point *c2); // EUCLIDEAN
// Find the cosine of the angle at a point
double cos_angle(point *P, const Metric &G);
double cos_angle(point *P); // EUCLIDEAN
// Find the sine of the angle at a point
double sin_angle(point *P);
// Find the actual angle at a point
double turning_angle(point *P, const Metric &G);
// Extrinsic curvature - from the osculating circle
double curv(point *P); 
// Geodesic curvature 
double geodesic_curvature(point *P, const Metric &G);
double geodesic_curvature_solo_angulo(point *P, const Metric &G);
// parallel transport
Vector parallel_transport(const Vector &, point *p1, point *p2, 
			  const Metric &G);

// AUXILIARY FUNCTIONS
// Find the determinant of the 2x2 matrix
double det2(double a11, double a12, double a21, double a22);
// Find whether x belongs to [a,b] // no need for a<b
bool is_between(double x, double a, double b);
// Find whether [x1,x2] and [y1,y2] overlap // no need for x1<x2, y1<y2
bool overlap(double x1, double x2, double y1, double y2);

#endif

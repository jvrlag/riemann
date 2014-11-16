#include"filam.h"

Filament::Filament(Metric &g) :  G(g)
{
     N=1;
     root=(point*)malloc(sizeof(point));
     root->next=root;
     root->prev=root;
     Ncuts=0;
     length_cut=0.0;

     maxl=1e-3;
     minl=1e-4;
}

Filament::~Filament()
{
     Destroy();
     free(root);
}

void Filament::Destroy()
{
     point *P=root, *Q;
     do
     {
          Q=P->next;
          Remove_Point(P);
          P=Q;
     }while(P!=P->next);
//     free(P);
     root=P;
}


// CHECKS. Consistency, Self-intersection and resolution

void Filament::Check_Consistency()
{
     point* P=root;
     int n=0;
     if (!P || !(P->next)) merror("No points!");
     if (P!=(P->prev)->next) merror("Open filament!");
     do
     {
	  if (isnan(P->x+P->x2+P->y+P->y2)) merror("Nan error!");
          if ((P->next)->prev!=P)
               merror("Linking error!\n");
          P=P->next;
          n++;
     }while(P!=root);
}

void Filament::Check_Resolution()
{
     if (N<10) merror("Filament has collapsed!!!!\n");
     bool done;
     do
     {
	  done=true;
	  point *P=root;
	  do
	  {
	       point *Q=P->next;
	       double dist=distance(P,Q,G);
	       if (dist<minl)
	       {
		    Remove_Point(P);
		    done=false;
		    P=Q;
	       }
	       if (dist>maxl)
	       {
		    point *P2=Add_Point(P);
		    Vector v=get_vector(P,Q);
		    P2->x2=P->x+0.5*v(1);
		    P2->y2=P->y+0.5*v(2);
		    update(P2);
		    P=P2;
		    done=false;
	       }
	       P=P->next;
	  }while(P!=root);
     }while(!done);
}


// To cut from point P1 to point P2 out of a filament (both of them survive!)
// Ah, the sense is important! P1->....->P2
void Filament::Cut(point *P1, point *P2)
{
     point *P=P1->next;
     do
     {
	  point *Q=P->next;
	  free(P);
	  N--;
	  P=Q;
     }while(P!=P2);

     P1->next=P2;
     P2->prev=P1;
}


// return true if, when going from P to Q, we find root
bool is_root_between(const Filament &F, point *P, point *Q)
{
     // we're assuming that both P and Q belong to F... otherwise, disaster!
     point *A=P;
     do
     {
	  if (A==F.root) return true;
	  A=A->next;
     }while(A!=Q);
     return false;
}

// Distance, ALONG the curve, between points P and Q
double Filament::curv_distance(point *P, point *Q) const
{
     point *A=P;
     double d=0.0;
     while(A!=Q) 
     {
	  d+=distance(A,A->next,G);
	  A=A->next;
     }
     return d;
}

void Filament::Check_Selfintersect()
{
     Cut_Info CI;
     do{} while(Remove_Selfintersect(CI));
}

bool Filament::Remove_Selfintersect(Cut_Info &K)
{
     point *P=root, *Q;
     bool done=false; // true if we've cut something: only one cut per check!
     for (long i=1;i<N;i++)
     {
	  Q=P->next->next;
	  for (long j=1;j<N/2;j++)
	  {
	       if (intersect_q(P,P->next,Q,Q->next))
	       {
		    K.length=curv_distance(P,Q);
		    K.crumb=left_q(P,P->next,Q); // it's a crumb if P->next is on the left of P->Q
		    Cut(P,Q);
		    root=Q;
		    done=true;
		    break;
	       }
	       Q=Q->next;
	  }
	  P=P->next;
	  if (done) break;
     }
     return done;
}

// void Filament::Check_Selfintersect()
// {
//      point *P=root;
//      do
//      {
// 	  point *Q=P->next->next;
// 	  long d=2; // distance from P to Q
// 	  do
// 	  {
// 	       if (intersect_q(P,P->next,Q,Q->next))
// 	       {
// 		    Q=Q->next;
// 		    Ncuts++;
// 		    // Find out if root is going to disappear!
// 		    if (is_root_between((*this),P,Q)) 
// 			 root=Q; 
// 		    if (d<N/2)
// 			 Cut(P,Q);
// 		    else Cut(Q,P);
// 		    break;
// 	       }	       
// 	       Q=Q->next;
// 	       d++;
// 	  }while(Q->next->next!=P);
// 	  P=P->next;
//      }while(P->next!=root);
// }



// ADDITION AND REMOVAL OF POINTS

// Create an "empty" filament of n points, all of them linked
void Filament::Create_Backbone(long n)
{
     point *P=root;
     for (long i=2;i<=n;i++)
	  P=Add_Point(P);
} 

// Add a point, between "P" and "P->next".
point* Filament::Add_Point(point* P)
{
     point* Q=(point*)malloc(sizeof(point));
     Q->next=P->next;
     Q->prev=P;
     P->next->prev=Q;
     P->next=Q;
     N++;
     return Q;
}

point* Filament::Add_Point(point* P, double x, double y)
{
     point* P2=Add_Point(P);
     P2->x=P2->x2=x;
     P2->y=P2->y2=y;
     return P2;
}

// remove a point and re-link the chain
void Filament::Remove_Point(point *P)
{
     if (!N) merror("No more points in filament!\n");
     if (P==root) root=P->next;

     point* P0=P->prev;
     point* P1=P->next;
     free(P);
     P0->next=P1;
     P1->prev=P0;
     N--;
}


void Filament::Write() const
{
     Write("stdout");
     // printf("# Output information about filament\n");
     // printf("# Number of points: %ld\n",N);
     // point *P=root;
     // for (long i=1;i<=N;i++)
     // {
     // 	  printf("%15.12g %15.12g\n",P->x,P->y);
     // 	  P=P->next;
     // }
}

void Filament::Write(char *filename) const
{
     FILE *fout=fopen(filename,"wt");
     fprintf(fout,"# %ld\n",N);
     point *P=root;
     for (long i=1;i<=N;i++)
     {
	  fprintf(fout,"%10g %10g\n",P->x,P->y);
	  P=P->next;
     }
     fclose(fout);
}

void Filament::Read(char *filename)
{
     FILE *fin=fopen(filename,"rt");
  
     long n;
     fscanf(fin,"# %ld\n",&n);
     // Now, read all the points
     Create_Backbone(n);

     point *P=root;
     for (long i=1;i<=N;i++)
     {
	  float x, y;
	  fscanf(fin,"%g %g\n",&x,&y);
	  P->x2=x; P->y2=y;
	  P=P->next;
     }
     Update();
     fclose(fin);
}

void Copy(Filament &F2, const Filament &F)
{
//      printf("Entramos Copy\n");
     F2.G=F.G; // copy the metric
     F2.Destroy();
     F2.Create_Backbone(F.N);
     point *P2=F2.root, *P=F.root;
     for (long i=1;i<=F2.N;i++)
     {
	  P2->x2=P->x2; 
	  P2->y2=P->y2;
	  P=P->next;
	  P2=P2->next;
     }
     F2.Update();
//      printf("Salimos copy\n");
}

// Create a circle, with points randomly scattered around it
void filament_random_circle(Filament &F, double R0, long n)
{
     Vector Angle(n);
     for (long i=1;i<=n;i++)
	  Angle(i)=rand_double(0.0,2.0*M_PI);
     Angle.Order();
     
     point* P=F.root;
     P->x2=R0*cos(Angle(1));
     P->y2=R0*sin(Angle(1));
     for (long i=2;i<=n;i++)
     {
// 		 double R=R0+0.5*cos(Angle(i));
//		 double R=R0*(1.0+0.2*cos(Angle(i)*5.0));
	  double x=R0*cos(Angle(i)), y=R0*sin(Angle(i)); 
          P=F.Add_Point(P,x,y);
     }
     F.Update();
}

void filament_circle(Filament &F, double R0, long n)
{
     Vector Angle(n);
     for (long i=1;i<=n;i++)
	  Angle(i)=(i-1)*2.0*M_PI/(double)n;
     Angle.Order();
     
     point* P=F.root;
     P->x2=R0*cos(Angle(1));
     P->y2=R0*sin(Angle(1));
     for (long i=2;i<=n;i++)
     {
	  // 		 double R=R0+0.5*cos(Angle(i));
	  //		 double R=R0*(1.0+0.2*cos(Angle(i)*5.0));
	  double x=R0*cos(Angle(i)), y=R0*sin(Angle(i)); 
	  P=F.Add_Point(P,x,y);
     }
          F.Update();
}

void filament_otra_cerrada_suave(Filament &F, double R0, double A, double k, long n)
{
     Vector Angle(n);
     for (long i=1;i<=n;i++)
	  Angle(i)=(i-1)*2.0*M_PI/(double)n;
     Angle.Order();
     
     point* P=F.root;
     P->x2=(R0+A*cos(k*Angle(1)))*cos(Angle(1));
     P->y2=(R0+A*cos(k*Angle(1)))*sin(Angle(1));
     for (long i=2;i<=n;i++)
     {
	  // 		 double R=R0+0.5*cos(Angle(i));
	  //		 double R=R0*(1.0+0.2*cos(Angle(i)*5.0));
	  double x=(R0+A*cos(k*Angle(i)))*cos(Angle(i)), y=(R0+A*cos(k*Angle(i)))*sin(Angle(i)); 
	  P=F.Add_Point(P,x,y);
     }
               F.Update();
}

void Filament::Update()
{
     point *P=root;
     do
     {
	  update(P);
	  P=P->next;
     }while(P!=root);
}


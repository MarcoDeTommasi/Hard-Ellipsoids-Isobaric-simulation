#include<stdlib.h>
#include<stdio.h>
#include<math.h>
/* perram wertheim overlap of hard ellipsoids */
#define Sqr(x) ((x)*(x))
double scalProd(double *A, double *B)
{
  int kk;
  double R=0.0;
  for (kk=0; kk < 3; kk++)
    R += A[kk]*B[kk];
  return R;
}

void tRDiagRpw(double M[3][3], double D[3], double Ri[3][3])
{
  int k1, k2, k3;
  double Di[3][3];
  double Rtmp[3][3];
  /* calcolo del tensore d'inerzia */ 
  Di[0][0] = D[0];
  Di[1][1] = D[1];
  Di[2][2] = D[2];
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	if (k1 != k2)
	  Di[k1][k2] = 0.0;
      } 
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	Rtmp[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    if (Di[k1][k3] == 0.0)
	      continue;
	    Rtmp[k1][k2] += Di[k1][k3]*Ri[k3][k2];
	  }
      }
  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      {
	M[k1][k2] = 0.0;
	for (k3=0; k3 < 3; k3++)
	  {
	    M[k1][k2] += Ri[k3][k1]*Rtmp[k3][k2];
	  }
      }
}
void xlambda(double lambda, double rA[3], double A[3][3], double rB[3], double B[3][3], double x[3])
{
  double lamA[3][3], onemlamB[3][3], ABL[3][3], invABL[3][3];
  double x1[3], x2[3], x3[3], detinvABL;
  int k1, k2;
  /* calculate xlambda, see L. Paramonov and S. N. Yaliraki J. Chem. Phys. 123, 194111 (2005) */
  for (k1=0; k1 < 3; k1++)
    {
      for (k2=0; k2 < 3; k2++)
	{
	  lamA[k1][k2] = lambda*A[k1][k2];
	  onemlamB[k1][k2] = (1.0-lambda)*B[k1][k2];
 	  ABL[k1][k2] = lamA[k1][k2] + onemlamB[k1][k2];
	}
    }
  for (k1=0; k1 < 3; k1++)
    {
      x1[k1]=0;
      x2[k1]=0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  x1[k1] += lamA[k1][k2]*rA[k2];
	  x2[k1] += onemlamB[k1][k2]*rB[k2];
	}
      x3[k1] = x1[k1] + x2[k1];
    }
  detinvABL=-ABL[0][2]*ABL[1][1]*ABL[2][0] + ABL[0][1]*ABL[1][2]*ABL[2][0] + 
    ABL[0][2]*ABL[1][0]*ABL[2][1] - ABL[0][0]*ABL[1][2]*ABL[2][1] - 
    ABL[0][1]*ABL[1][0]*ABL[2][2] + ABL[0][0]*ABL[1][1]*ABL[2][2]; 

  invABL[0][0] = -ABL[1][2]*ABL[2][1] + ABL[1][1]*ABL[2][2];
  invABL[0][1] =  ABL[0][2]*ABL[2][1] - ABL[0][1]*ABL[2][2];
  invABL[0][2] = -ABL[0][2]*ABL[1][1] + ABL[0][1]*ABL[1][2];
  invABL[1][0] =  ABL[1][2]*ABL[2][0] - ABL[1][0]*ABL[2][2]; /* a12 a20 - a10 a22 */
  invABL[1][1] = -ABL[0][2]*ABL[2][0] + ABL[0][0]*ABL[2][2]; /* -a02 a20 + a00 a22 */
  invABL[1][2] =  ABL[0][2]*ABL[1][0] - ABL[0][0]*ABL[1][2]; /* a02 a10 - a00 a12 */
  invABL[2][0] = -ABL[1][1]*ABL[2][0] + ABL[1][0]*ABL[2][1]; /* -a11 a20 + a10 a21 */
  invABL[2][1] =  ABL[0][1]*ABL[2][0] - ABL[0][0]*ABL[2][1]; /* a01 a20 - a00 a21 */
  invABL[2][2] = -ABL[0][1]*ABL[1][0] + ABL[0][0]*ABL[1][1]; /* -a01 a10 + a00 a11 */

  for (k1 = 0; k1 < 3; k1++)
    for (k2 = 0; k2 < 3; k2++)
      invABL[k1][k2] /= detinvABL;

  for (k1 = 0; k1 < 3; k1++)
    {
      x[k1] = 0.0;
      for (k2 = 0; k2 < 3; k2++)
	{
	  x[k1] += invABL[k1][k2]*x3[k2];
	}
    }
}

double Slam(double lambda, double rA[3], double A[3][3], double rB[3], double B[3][3])
{
  int k1, k2;
  double xlam[3], fA[3], fB[3], SA, SB;

  xlambda(lambda, rA, A, rB, B, xlam);

  for (k1=0; k1 < 3; k1++)
    {
      fA[k1] = 0;
      fB[k1] = 0;
      for (k2=0; k2 < 3; k2++)
	{
	  fA[k1] += A[k1][k2]*(xlam[k2]-rA[k2]);
	  fB[k1] += B[k1][k2]*(xlam[k2]-rB[k2]);
	}
    }

  SA = SB = 0.0;
  for (k1=0; k1 < 3; k1++)
    {
      SA += lambda*(xlam[k1]-rA[k1])*fA[k1];
      SB += (1.0-lambda)*(xlam[k1]-rB[k1])*fB[k1];
    }
  /* ho messo un - così la funzione ha un minimo invece
     che un massimo e questo minimo viene trovato dalla funzione brentPW */
  //printf("C SA+SB=%.15G\n", SA+SB);
  return -(SA+SB);
}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); 
int brentPWTooManyIter=0;
double brentPW(double ax, double bx, double cx, double tol, double *xmin, double rA[3], double A[3][3], double rB[3], double B[3][3])
/*Given a function f, and given a bracketing triplet of abscissas ax, bx, cx 
 * (such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)),
 * this routine isolates the minimum to a fractional precision of about tol using Brent's
 * method. The abscissa of the minimum is returned as xmin, and the minimum function value 
 * is returned as brent, the returned function value. */
{ 
  int iter, ITMAXBR=100;
  const double CGOLD=0.3819660;
  const double ZEPSBR=1E-20;
  double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0, fuold;
  brentPWTooManyIter=0;
  /* This will be the distance moved on the step before last.*/
  a=(ax < cx ? ax : cx); /*a and b must be in ascending order, 
			   but input abscissas need not be.*/
  b=(ax > cx ? ax : cx);
  x=w=v=bx; /*Initializations...*/
  fw=fv=fx=Slam(x, rA, A, rB, B); 

  //printf("C fw=%.15G\n", fw);
  if (fw < -1.0)
    {
      /* non-overlap! */
      *xmin=x;
      return -100.0;
    }
  fuold = fv;
  for (iter=1;iter<=ITMAXBR;iter++)
    { 
      /*Main program loop.*/
      //printf("C iter=%d fx=%.15G\n", iter, fx);
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPSBR); 
      if (fabs(x-xm) <= (tol2-0.5*(b-a)))
	{ /*Test for done here.*/
	  *xmin=x;
	  return fx;
	} 
      if (fabs(e) > tol1) 
	{ /*Construct a trial parabolic fit.*/
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if (q > 0.0)
	    p = -p; 
	  q=fabs(q);
	  etemp=e; 
	  e=d; 
	  if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=CGOLD*(e=(x >= xm ? a-x : b-x)); 
	    /*The above conditions determine the acceptability of the parabolic fit.
	     * Here we take the golden section step into the larger of the two segments.*/
	  else
	    {
	      d=p/q; /* Take the parabolic step.*/
	      u=x+d; 
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-x); 
	    }
	}
      else
	{
	  d=CGOLD*(e=(x >= xm ? a-x : b-x));
	} 
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      fu=Slam(u, rA, A, rB, B); /*This is the one function evaluation per iteration.*/
      if (fu < -1.0)
	{
	  /* non overlap! */
	  *xmin=x;
	  return -100.0;
	}
      fuold = fu;
      if (fu <= fx)
	{ /*Now decide what to do with our function evaluation.*/
	  if (u >= x) 
	    a=x;
	  else
	    b=x;
	  SHFT(v,w,x,u); /* Housekeeping follows:*/
	  SHFT(fv,fw,fx,fu); 
	} 
      else
	{ 
	  if (u < x) 
	    a=u; 
	  else 
	    b=u; 
	  if (fu <= fw || w == x)
	    {
	      v=w; w=u; fv=fw; fw=fu;
	    }
	  else if (fu <= fv || v == x || v == w)
	    { 
	      v=u; fv=fu;
	    }
	} /* Done with housekeeping. Back for another iteration.*/
    }
  printf("Too many iterations in brent!\n");
  brentPWTooManyIter=1;
  //nrerror("Too many iterations in brent"); 
  *xmin=x; /*Never get here.*/
  return fx;
}
/* 
 * saxi[3]: semiaxes of first hard ellipsoid
 * saxj[3]: semiaxes of second hard ellipsoid
 * rA[3]: center of mass of first ellipsoid
 * rB[3]: center of mass of second ellipsoid
 * RA[3][3]: orientation matrix of first ellipsoid, where each row is a unit vector of hard body reference system,
 * i.e. RA[0][0..2] is the unit vector of x-axis of hard body reference system
 *      RA[1][0..2] is the unit vector of y-axis of hard body reference system
 *      RB[2][0..2] is the unit vector of z-axis of hard body reference system
 * Since the hard ellipsoids are uniaxial, you can generate randomly a unit vector and then build two vectors perpendicular
 * to the first one by the Gram-Schmidt orthogonalization procedure.      
 * RB[3][3]: orientation matrix of second
 * */
double check_overlap_pw_c(double saxi[3], double rA[3], double RA[3][3], double saxj[3], double rB[3], double RB[3][3])
{
  const double tolPW=1.0E-12;
  double res, A[3][3], B[3][3], xmin; 
  double  DA[3], DB[3]; 
  int k1;
  for (k1=0; k1 < 3; k1++)
    {
      DA[k1]= 1.0/Sqr(saxi[k1]);
      DB[k1]= 1.0/Sqr(saxj[k1]);
    }
  
  tRDiagRpw(A, DA, RA);
  tRDiagRpw(B, DB, RB);

  res =  -brentPW(0, 0.5, 1.0, tolPW, &xmin, rA, A, rB, B);
  if (brentPWTooManyIter)
    {
      printf("res=%f xmin=%f\n", res, xmin);
      exit(-1);
    }
  return res - 1.0;
}

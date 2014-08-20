
#include "euclideanvector.h"

#include "linalg.h"

#include <iostream>
#include <stdio.h>


#include <gmpxx.h>

using namespace std;


bool computesubfactors(EuclideanVector *orthobasis, int ambientdimension)
{

  
  if(ambientdimension == 2) {
    orthobasis[1].data[0] = orthobasis[0].data[1];
    orthobasis[1].data[1] = -1 * orthobasis[0].data[0];
    return true;
  }

  
   if(ambientdimension == 3) {
    orthobasis[2].data[0] = orthobasis[0].data[1]*orthobasis[1].data[2] -  orthobasis[0].data[2]*orthobasis[1].data[1] ;    

    orthobasis[2].data[1] = -1*orthobasis[0].data[0]*orthobasis[1].data[2] +  orthobasis[0].data[2]*orthobasis[1].data[0] ;    

    orthobasis[2].data[2] = orthobasis[0].data[0]*orthobasis[1].data[1] -  orthobasis[0].data[1]*orthobasis[1].data[0] ;    

    if((orthobasis[2].data[0]*orthobasis[2].data[0]) == 0)
    if((orthobasis[2].data[1]*orthobasis[2].data[1]) == 0)
    if((orthobasis[2].data[2]*orthobasis[2].data[2]) == 0)
      return false;

    return true;
  } 

  if(ambientdimension == 4) {

   mpq_class a,b,c,d,e,f,g,h,i,j,k,l;
   
   a = orthobasis[0].data[0];
   b = orthobasis[0].data[1];
   c = orthobasis[0].data[2];
   d = orthobasis[0].data[3];

   e = orthobasis[1].data[0];
   f = orthobasis[1].data[1];
   g = orthobasis[1].data[2];
   h = orthobasis[1].data[3];

   i = orthobasis[2].data[0];
   j = orthobasis[2].data[1];
   k = orthobasis[2].data[2];
   l = orthobasis[2].data[3];

   orthobasis[3].data[0] = b*g*l + f*k*d + c*h*j - b*k*h - f*c*l - j*g*d;
   orthobasis[3].data[1] = -1*(a*g*l + e*k*d + c*h*i - a*k*h - e*c*l - i*g*d);

   orthobasis[3].data[2] = (a*f*l + e*j*d + b*h*i - a*j*h - e*b*l - i*f*d);
   orthobasis[3].data[3] = -1*(a*f*k + e*j*c + b*g*i - a*j*g - e*b*k - i*f*c);

   if((orthobasis[3].data[0]*orthobasis[3].data[0]) == 0)
   if((orthobasis[3].data[1]*orthobasis[3].data[1]) == 0)
   if((orthobasis[3].data[2]*orthobasis[3].data[2]) == 0)
   if((orthobasis[3].data[3]*orthobasis[3].data[3]) == 0)
     return false;

   return true;

  }

  if(ambientdimension == 5) {

   mpq_class a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t;
   
   a = orthobasis[0].data[0];
   b = orthobasis[0].data[1];
   c = orthobasis[0].data[2];
   d = orthobasis[0].data[3];
   e = orthobasis[0].data[4];

   f = orthobasis[1].data[0];
   g = orthobasis[1].data[1];
   h = orthobasis[1].data[2];
   i = orthobasis[1].data[3];
   j = orthobasis[1].data[4];

   k = orthobasis[2].data[0];
   l = orthobasis[2].data[1];
   m = orthobasis[2].data[2];
   n = orthobasis[2].data[3];
   o = orthobasis[2].data[4];

   p = orthobasis[3].data[0];
   q = orthobasis[3].data[1];
   r = orthobasis[3].data[2];
   s = orthobasis[3].data[3];
   t = orthobasis[3].data[4];
   
   orthobasis[4].data[4] = fourbyfourdet(a,b,c,d, f,g,h,i, k,l,m,n, p,q,r,s);
   orthobasis[4].data[3] = -1*fourbyfourdet(a,b,c,e, f,g,h,j, k,l,m,o, p,q,r,t);
   orthobasis[4].data[2] = fourbyfourdet(a,b,d,e, f,g,i,j, k,l,n,o, p,q,s,t);
   orthobasis[4].data[1] = -1*fourbyfourdet(a,c,d,e, f,h,i,j, k,m,n,o, p,r,s,t);
   orthobasis[4].data[0] = fourbyfourdet(b,c,d,e, g,h,i,j, l,m,n,o, q,r,s,t);

   if((orthobasis[4].data[0]*orthobasis[4].data[0]) == 0)
   if((orthobasis[4].data[1]*orthobasis[4].data[1]) == 0)
   if((orthobasis[4].data[2]*orthobasis[4].data[2]) == 0)
   if((orthobasis[4].data[3]*orthobasis[4].data[3]) == 0)
   if((orthobasis[4].data[4]*orthobasis[4].data[4]) == 0)
     return false;

   return true;
  }

  //if(gmp_rationals || ambientdimension > 5) 
  {
  
    mpq_class matrix[ambientdimension - 1][ambientdimension - 1];
    for(int i = 0; i < ambientdimension-1; i++)
    for(int j = 0; j < ambientdimension-1; j++)
      matrix[i][j] = 0;

  
    int sign = 1;
    for(int k = 0; k < ambientdimension; k++) {
      // make kth subfactor matrix
      for(int i = 0; i < ambientdimension - 1; i++) {
      // leave out kth column
        for(int j = 0; j < ambientdimension - 1; j++) {
          if(j < k)
            matrix[i][j] = orthobasis[i].data[j];
          else
            matrix[i][j] = orthobasis[i].data[j+1];
        }
      }

      #ifdef DEBUG2
      cout << "\nbefore triangularize...\n";
      printmatrix((mpq_class *)(matrix),ambientdimension-1,ambientdimension-1);
      #endif

      uppertriangular((mpq_class *)(matrix), ambientdimension-1, ambientdimension-1);

      
      #ifdef DEBUG2
      cout << "\n...and after\n";
      printmatrix((mpq_class *)(matrix),ambientdimension-1,ambientdimension-1);
      #endif

      orthobasis[ambientdimension-1].data[k] = sign * diagonalproduct((mpq_class *)(matrix), ambientdimension-1); 

      sign *= -1;
    } 

  
    #ifdef DEBUG2
    cout << "\n Computed normal vector was:\n";
    orthobasis[ambientdimension-1].Print();
    #endif  

    for(int i = 0; i < ambientdimension; i++)
      if(orthobasis[ambientdimension-1].data[i] != 0)
        return true;

    return false;

  }

  return true;
}

inline mpq_class threebythreedet(mpq_class a, mpq_class b, mpq_class c, mpq_class d, mpq_class e, mpq_class f, mpq_class g, mpq_class h, mpq_class i) {
  return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
}


mpq_class fourbyfourdet(mpq_class a, mpq_class b, mpq_class c, mpq_class d, mpq_class e, mpq_class f, mpq_class g, mpq_class h, mpq_class i, mpq_class j, mpq_class k, mpq_class l, mpq_class m, mpq_class n, mpq_class o, mpq_class p) 
{

  mpq_class ans = 0;

  ans += a * threebythreedet(f,g,h, j,k,l, n,o,p);
  ans -= b * threebythreedet(e,g,h, i,k,l, m,o,p);
  ans += c * threebythreedet(e,f,h, i,j,l, m,n,p);
  ans -= d * threebythreedet(e,f,g, i,j,k, m,n,o); 

  return ans;

}

void EuclideanVector::Constructor(int d) {
  dimension = d;
  data = new mpq_class[d];

  for(int i = 0; i < d; i++)
    data[i] = 0;
}

void EuclideanVector::Print(){
  #ifdef DEBUG2
  cout << "Euclidean Vector " << this << "with data located at " << data << ":  ";
  #endif

  for(int i = 0; i < dimension; i++) {
    printNumber(data[i]);
    cout << " ";
  } 

  cout << "\n";
}

mpq_class dotproduct(EuclideanVector *v, EuclideanVector *w)
{
  if(v->dimension != w->dimension) {
    cerr << "Dot product of two vectors of unequal dimension " << v->dimension <<"," << w->dimension <<  "!";
    return (mpq_class)(0);
  }

  mpq_class ans;
  ans = 0;

  for(int i = 0; i < v->dimension; i++)
    ans += v->data[i] * w->data[i];

  return ans;

}

void EuclideanVector::timesequals(mpq_class c) {
  for(int i = 0; i < dimension; i++)
    data[i] *= c;

  return ;
}

void EuclideanVector::minusequals(EuclideanVector w) {
  for(int i = 0; i < dimension; i++)
    data[i] -= w.data[i];

  return;
}


void EuclideanVector::plusequals(EuclideanVector w) {
  for(int i = 0; i < dimension; i++)
    data[i] += w.data[i];

  return ;
}

void EuclideanVector::Negation() {
  for(int i =0; i < dimension; i++)
    data[i] = -1 * data[i];

  return;
} 

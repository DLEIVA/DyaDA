/* C routines for the DyaDA package */

/******************************************************************************************************************************/
/******************************************************************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>


/* Declaration of functions */

  void linearh(double *X, int *nrow, int *rep, double *pvhright, double *pvhleft);
  double getlandau(double *mat_V, int *maxrow, int *maxcol);
  void statsigZ(double *X, double *Y, int *nrow, int *perm, double *pvZright, double *pvZleft);
  void rectangularstatsigZ(double *X, double *Y, int *nrow, int *ncol, int *perm, double *pvZright, double *pvZleft);
  void swap(int *setn, const int i, const int j);
  void permuteZ(int *setn,double *matperm,double *mat_X, double *mat_Y, double *results, double z, const int start, 
                int maxrow, int maxcol);
  void rotateLeft(int *setn, const int start, const int maxrow);
  void exactsigZ(double *X, double *Y, int *nrow, double *pvZright, double *pvZleft);
  void statsigR(double *X, double *Y, int *nrow, int *perm, double *pvRright, double *pvRleft);
  double getR(double *mat_X, double *mat_Y, int *maxrow, int *maxcol);
  void rectangularstatsigR(double *X, double *Y, int *nrow, int *ncol, int *perm, double *pvRright, double *pvRleft);
  void permuteR(int *setn,double *matperm,double *mat_X, double *mat_Y, double *results, double R,
                const int start, int maxrow, int maxcol);
  void exactsigR(double *X, double *Y, int *nrow, double *pvRright, double *pvRleft);
  double getKr(double *mat_X, double *mat_Y, int *maxrow, int *maxcol);
  void statsigKr(double *X, double *Y, int *nrow, int *perm, double *pvKrright, double *pvKrleft);
  void rectangularstatsigKr(double *X, double *Y, int *nrow, int *ncol, int *perm, double *pvKrright, double *pvKrleft);
  void permuteKr(int *setn,double *matperm,double *mat_X, double *mat_Y, double *results, double Kr,
                 const int start, int maxrow, int maxcol);
  void exactsigKr(double *X, double *Y, int *nrow, double *pvKrright, double *pvKrleft);
  double getZr(double *mat_X, double *mat_Y, double *wi,double *ri, double *Xterm, double *Yterm, int *maxrow, int *maxcol);
  void statsigZr(double *X, double *Y, int *nrow, int *perm, double *pvZrright, double *pvZrleft);
  void rectangularstatsigZr(double *X, double *Y, int *nrow, int *ncol, int *perm, double *pvZrright, double *pvZrleft);
  void permuteZr(int *setn,double *matperm,double *mat_X, double *mat_Y, double *wi, double *ri, double *Xterm,
                 double *Yterm, double *results, double Zr, const int start, int maxrow, int maxcol);
  void exactsigZr(double *X, double *Y, int *nrow, double *pvZrright, double *pvZrleft);
  double getRr(double *mat_X, double *mat_Y, double *wi,double *ri, double *Xterm, double *Yterm, int *maxrow, int *maxcol);
  void statsigRr(double *X, double *Y, int *nrow, int *perm, double *pvRrright, double *pvRrleft);
  void rectangularstatsigRr(double *X, double *Y, int *nrow, int *ncol, int *perm, double *pvRrright, double *pvRrleft);
  void permuteRr(int *setn,double *matperm,double *mat_X, double *mat_Y, double *wi, double *ri, double *Xterm,
                 double *Yterm, double *results, double Rr, const int start, int maxrow, int maxcol);
  void exactsigRr(double *X, double *Y, int *nrow, double *pvRrright, double *pvRrleft);
  double getrrw(double *mat_X, double *mat_Y, double *wi,double *ri, double *X1term, double *X2term, int *maxrow, int *maxcol);
  double getrrwxyz(double *mat_X, double *mat_Y, double *mat_Z, double *wi,double *ri, double *X1term, double *X2term, 
                   int *maxrow, int *maxcol);
  void partialstatsigZr(double *X, double *Y, double *Z, int *nrow, int *perm, double *pvpartialZrright,
                        double *pvpartialZrleft);
  void rectangularstatsigpartZr(double *X, double *Y, double *Z, int *nrow, int *ncol, int *perm, double *pvpartialZrright,
                                double *pvpartialZrleft);
  void partialexactsigZr(double *X, double *Y, double *Z, int *nrow, double *pvpartialZrright, double *pvpartialZrleft);
  void permuterrwxyz(int *setn,double *matperm,double *mat_X, double *mat_Y, double *mat_Z, double *wi, double *ri, 
                     double *X1term,double *X2term, double *results, double rrwxyz, const int start, int maxrow, int maxcol);
  double getrhorw(double *mat_X, double *mat_Y, double *wi,double *ri, double *X1term, double *X2term, int *maxrow,
                  int *maxcol);
  double getrhorwxyz(double *mat_X, double *mat_Y, double *mat_Z, double *wi,double *ri, double *X1term, double *X2term,
                     int *maxrow, int *maxcol);
  void partialstatsigRr(double *X, double *Y, double *Z, int *nrow, int *perm, double *pvpartialRrright,
                        double *pvpartialRrleft);
  void permuterhorwxyz(int *setn,double *matperm,double *mat_X, double *mat_Y, double *mat_Z, double *wi, double *ri,
                       double *X1term, double *X2term, double *results, double rhorwxyz, const int start, int maxrow, 
                       int maxcol);
  void partialexactsigRr(double *X, double *Y, double *Z, int *nrow, double *pvpartialRrright, double *pvpartialRrleft);
  double gettaurw(double *mat_X, double *mat_Y, double *wi,double *taui, double *X1term, double *X2term, int *maxrow,
                  int *maxcol);
  double gettaurwxyz(double *mat_X, double *mat_Y, double *mat_Z, double *wi,double *taui, double *X1term, double *X2term,
                     int *maxrow, int *maxcol);
  void partialstatsigKr(double *X, double *Y, double *Z, int *nrow, int *perm, double *pvpartialKrright,
                        double *pvpartialKrleft);
  void rectangularstatsigpartKr(double *X, double *Y, double *Z, int *nrow, int *ncol, int *perm, double *pvpartialKrright,
                                double *pvpartialKrleft);
  void permutetaurwxyz(int *setn,double *matperm,double *mat_X, double *mat_Y, double *mat_Z, double *wi, double *taui,
                       double *X1term, double *X2term, double *results, double taurwxyz, const int start, int maxrow,
                       int maxcol);
  void partialexactsigKr(double *X, double *Y, double *Z, int *nrow, double *pvpartialKrright, double *pvpartialKrleft);
  void getphirmat(double *mat_X, double *dyadc, double *phirmat, int *maxrow, int *maxcol);
  double getPHIr(double *mat_X, double *dyadc, double *phirmat, int *maxrow, int *maxcol);
  void getphii (double *mat_X,double *dyadc,double *phirmat, double *phii, int *maxrow, int *maxcol);
  void getphij (double *mat_X,double *dyadc,double *phirmat, double *phij, int *maxrow, int *maxcol);
  void getPHIij(double *mat_X, double *dyadc, double *phirmat,double *PHIij, int *maxrow, int *maxcol);
  void recip(double *X, double *pi, int *nrow, int *rep, double *pvPHIrright, double *pvPHIrleft, double *pvphirmatright,
             double *pvphirmatleft, double *pvphiiright, double *pvphiileft, double *pvphijright, double *pvphijleft,
             double *pvPHIijright, double *pvPHIijleft);
  void recip0(double *X, double *pi, int *nrow, int *rep, double *pvphi, double *pvpsi, double *pvdc, double *pvdelta, double *pvepsilon,
              double *pvkappa, double *pvnuj, double *pvlambdaj, double *pvrationu, double *pvratiolambda, double *pvomega);  
  double getpsi(double *mat_X, double *tras_X, double *mat_S, double *tras_S, double *mat_XTX,
                double *mat_STS, int *maxrow, int *maxcol);
  double getphi(double *mat_X, double *tras_X, double *mat_K, double *tras_K, double *mat_XTX,
                double *mat_KTK, int *maxrow, int *maxcol);
  double getdelta(double *mat_X, double *tras_X, double *mat_S, double *tras_S, double *mat_XTX, double *mat_STS,
                  double *mat_K, double *tras_K, double *mat_KTK, int *maxrow, int *maxcol);
  double getdc(double *mat_X, double *tras_X, int *maxrow, int *maxcol);
  double getlambda(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *lambda,
                   int *maxrow, int *maxcol);
  double getlambdaj(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *lambda,
                    double *lambdaj, int *maxrow, int *maxcol);
  double getratiolambda(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *lambda,
                        double *lambdaj, double *ratiolambda, int *maxrow, int *maxcol);
  double getnu(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *nu,
               int *maxrow, int *maxcol);
  double getnuj(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *nu,
                double *nuj, int *maxrow, int *maxcol);
  double getrationu(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *nu,
                    double *nuj, double *rationu, int *maxrow, int *maxcol);
  double getomega(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *lambda,
                  double *nu, double *omega, int *maxrow, int *maxcol);
  double getkappa(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *nu,
                  double *nuj, double *rationu, double *difrationu, int *maxrow, int *maxcol);  
  double getepsilon(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *nu,
                    double *nuj, double *nui, int *maxrow, int *maxcol);  
  void DminScomp(double *mat_X, double *minS, int*maxrow);
  int Icomp(double *mat_X, int *maxrow);
  int SIcomp(double *mat_X, int *maxrow);
  void swaprowcol (double *mat_X, int *k, int *m, int *maxrow);
  void swapnames (char **namesOrd, int *k, int *m, char *temp);
  void ISImethod (double *X, int *nrow, char **vecNames, int *tries, double *matord,
                  char **namord);  

/******************************************************************************************************************************/
/******************************************************************************************************************************/

/* Set of functions to generate random sociomatrices under the null hypothesis and compute Landau's h statistics */

/* Function to compute Landau's linear hierarchy index */

double getlandau(double *mat_V, int *maxrow, int *maxcol)
{
  int i, j;
  double *Vvector,term,sumabil,h;
  
  /* Allocate in memory size of vector Vvector */
  
  Vvector = calloc(*maxrow, sizeof(double));   
  
  for (i= 0; i< *maxrow; i++)
    Vvector[i] = 0.;
  for (i= 0; i< *maxrow; i++)
    for (j= 0; j< *maxcol; j++)
      Vvector[i] += mat_V[i**maxcol+j];
  
  sumabil = 0.;
  for (j= 0; j< *maxcol; j++){
    term = Vvector[j] - ((*maxrow) - 1.)/2.;
    sumabil += pow(term,2);
  }
  h = 12./(pow(*maxrow,3)-(*maxrow))*sumabil;
  
  /* Deallocate the matrices and vectors used in the routine */
  
  free(Vvector);   
  
  return(h);
}

void linearh(double *X, int *nrow, int *rep, double *pvhright, double *pvhleft)
{
  int i, j, m, maxrow, maxcol, iter;
  double *mat_X, *mat_V, numb, *matgen, h, hsim, epsilon;
  
  GetRNGstate();
  maxrow = maxcol = *nrow;
  
  /* Allocate in memory size of original matrix */
  
  mat_X = malloc(maxrow*maxcol*sizeof(double));
  
  if (mat_X == NULL)
  {
    error("Null dimension");
  }
  
  m = 0.;
  for (i= 0; i< maxrow; i++)
    for (j= 0; j< maxcol; j++)
    {   
      mat_X[i*maxcol+j] = X[m];
      m++;
    }
    
    /* Allocate in memory size of matrix of abilities */
    
    mat_V = malloc(maxrow*maxcol*sizeof(double));
  
  for (i= 0; i< (maxrow-1); i++)
    for (j= (i+1); j< maxcol; j++)
    {
      if (X[i*maxcol+j] > X[j*maxcol+i])
      {
        mat_V[i*maxcol+j] = 1.;
        mat_V[j*maxcol+i] = 0.;
      }
      if (X[i*maxcol+j] < X[j*maxcol+i])
      {
        mat_V[j*maxcol+i] = 1.;
        mat_V[i*maxcol+j] = 0.;
      }
      if ((X[i*maxcol+j] == X[j*maxcol+i]) & (X[i*maxcol+j] != 0.))
      {
        mat_V[i*maxcol+j] = 1./2.;
        mat_V[j*maxcol+i] = 1./2.;
      }
      if ((X[i*maxcol+j] == X[j*maxcol+i]) & (X[i*maxcol+j] == 0.))
      {
        numb = runif(0,1);
        if (numb >= 0.5)
          mat_V[i*maxcol+j] = 1.;
        else X[i*maxcol+j] = 0.;
        mat_V[j*maxcol+i] = 1-X[i*maxcol+j];
      }
    }
    
    /* Compute linear hierarchy measure for the original sociomatrix */
    
    h = getlandau(mat_V,&maxrow,&maxcol);
  
  matgen = malloc(maxrow*maxcol*sizeof(double));
  
  *pvhright = 0.;
  *pvhleft = 0.;
  
  for (iter= 0; iter< *rep; iter++) 
  {
    for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
        {  if (i < j){
          numb = runif(0,1);   
          if (numb >= 0.5){matgen[i*maxcol+j]=1;}
          else matgen[i*maxcol+j]=0;
          matgen[j*maxcol+i]=1-matgen[i*maxcol+j];}
        if (i == j){matgen[i*maxcol+j]=0.;}
        }
      
      hsim = getlandau(matgen,&maxrow,&maxcol);
    
    epsilon = fabs(hsim - h);
    if ( (hsim > h) & (epsilon > 0.000001) ) 
    {*pvhright += 1.;}
    if ( (hsim < h) & (epsilon > 0.000001) ) 
    {*pvhleft += 1.;}
    
    if ((hsim == h) | (epsilon <= 0.000001))
    {*pvhright += 1.;
      *pvhleft += 1.;}
  }
  
  /* Compute p values for the linear hierarchy index under the specified null hypothesis */
  
  *pvhright = (*pvhright+1.)/(*rep+1.);
  *pvhleft =(*pvhleft+1.)/(*rep+1.);
  
  PutRNGstate();
  free(matgen);
  free(mat_V);
  free(mat_X);
}


/******************************************************************************************************************************/
/******************************************************************************************************************************/

/* Set of functions to randomly permute sociomatrices and to estimate statistical significance for several matrix
   correlation statistics */

/* Function to obtain the Mantel's Z statistic */

   double getZ(double *mat_X, double *mat_Y, int *maxrow, int *maxcol)
   {
   int i,j;
   double z;	   
   z = 0.;
   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         {z += mat_X[i**maxcol+j]*mat_Y[i**maxcol+j];}
   return(z);
   }

/* Function to estimate statistical significance for Mantel's Z statistic for two square matrices */

  void statsigZ(double *X, double *Y, int *nrow, int *perm, double *pvZright, double *pvZleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, iter, *index, pointi, pointj;
   double *mat_X, *mat_Y, z, *matperm, *setn, zsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
  {
    error("Null dimension");
  }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Compute statistic value for the original sociomatrix */

   z = getZ(mat_X,mat_Y,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector setn and index */

   setn = calloc(maxrow, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices
   under the specified null hypothesis */

   *pvZright = 0.;
   *pvZleft = 0.;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow);
   for (i = 0; i< maxrow; i++){
      setn[i] = index[i];
   }
   for (i= 0; i< (maxrow-1); i++){
      pointi = setn[i];     
      for (j= i+1; j< maxcol; j++){
         pointj = setn[j];
         matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
         matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
      for (i = 0; i < maxrow; i++)
         for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}
   
   zsim = getZ(matperm,mat_Y,&maxrow,&maxcol);
   if (zsim >= z)
        {*pvZright += 1.;}
      if (zsim <= z)
        {*pvZleft += 1.;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvZright = (*pvZright+1.)/(*perm+1.);
   *pvZleft =(*pvZleft+1.)/(*perm+1.);
  }

/* Function to estimate statistical significance for Mantel's Z statistic for two rectangular matrices */

  void rectangularstatsigZ(double *X, double *Y, int *nrow, int *ncol, int *perm, double *pvZright, double *pvZleft)
  {
/* Declaration of local variables */

   int i, j, k, m, maxrow, maxcol, iter, *index, point;
   double *mat_X, *mat_Y, z, *matperm, *setn, zsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = *nrow;
   maxcol = *ncol;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
  {
  error("Null dimension");
  }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Compute statistic value for the original sociomatrix */

   z = getZ(mat_X,mat_Y,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector setn and index */

   setn = calloc(maxrow*maxcol, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow*maxcol, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices
   under the specified null hypothesis */

   *pvZright = 0.;
   *pvZleft = 0.;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow*maxcol; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow*maxcol);
   for (i = 0; i< maxrow*maxcol; i++){
      setn[i] = index[i];
   }

   k = 0.;
   for (i= 0; i< maxrow; i++){     
      for (j= 0; j< maxcol; j++){
         point = setn[k];
         matperm[i*maxcol+j] = mat_X[point];
         k += 1;}}
  
   zsim = getZ(matperm,mat_Y,&maxrow,&maxcol);
   if (zsim >= z)
        {*pvZright += 1.;}
      if (zsim <= z)
        {*pvZleft += 1.;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvZright = (*pvZright+1.)/(*perm+1.);
   *pvZleft =(*pvZleft+1.)/(*perm+1.);
}

/* Some functions required in the exact permutation test */

   void swap(int *setn, const int i, const int j)
   {
   int temp;
   temp = setn[i];
   setn[i] = setn[j];
   setn[j] = temp;
   } 

   void rotateLeft(int *setn, const int start, const int maxrow)
   {
   int temp = setn[start];
   for (int i = start; i < maxrow-1; i++) {
      setn[i] = setn[i+1];
   }
   setn[maxrow-1] = temp;
  }

/* A function required to estimate statistical significance for Mantel's Z statistic */

   void permuteZ(int *setn,double *matperm,double *mat_X, double *mat_Y, double *results, double z, 
	       const int start, int maxrow, int maxcol)
   {
    int i, j, pointi, pointj;
    double zsim;
    for (i= 0; i< (maxrow-1); i++){
       pointi = setn[i];     
       for (j= i+1; j< maxcol; j++){
          pointj = setn[j];
          matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
          matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
       for (i = 0; i < maxrow; i++)
          for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}
       zsim = getZ(matperm,mat_Y,&maxrow,&maxcol);
   if (zsim >= z)
   {results[0] += 1.;}
   if (zsim <= z)
   {results[1] += 1.;}
   if (start < maxrow) {
      int i, j;
      for (i = maxrow-2; i >= start; i--) {
         for (j = i + 1; j < maxrow; j++) {
	    swap(setn, i, j);
	    permuteZ(setn,matperm,mat_X,mat_Y,results,z,i+1,maxrow,maxcol);
         }
      rotateLeft(setn, i, maxrow);
    }
    }
   }

/* Function to compute exact statistical significance for Mantel's Z */
  
   void exactsigZ(double *X, double *Y, int *nrow, double *pvZright, double *pvZleft)
   {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, *setn, start, fact;
   double *mat_X, *mat_Y, z, *matperm, *results;

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Compute statistic value for the original sociomatrix */

   z = getZ(mat_X,mat_Y,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector setn */

   setn = calloc(maxrow, sizeof(int));

   if (setn == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector results and factorial number of n */

   fact = 1;
   for (i = maxrow;i >= 1; i--){
      fact *= i;}

   results = calloc(2,sizeof(double));

   if (results == NULL)
   {
     error("Null dimension");
   }

/* Carry out the exact permutation procedure */

   for (i = 0; i < maxrow; i++){
      setn[i] = i;}

   for (i = 0; i < 2; i++){
      results[i] = 0.;}

   start =0;
   permuteZ(setn,matperm,mat_X,mat_Y,results,z,start,maxrow,maxcol);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvZright = (results[0])/(fact);
   *pvZleft =(results[1])/(fact);

/* Deallocate all the matrices used in the routine */

   free(results);
   free(setn);
   free(matperm);
   free(mat_Y);
   free(mat_X);
}

/* Function to obtain Dietz's R statistic */

   double getR(double *mat_X, double *mat_Y, int *maxrow, int *maxcol)
   {
   int i,j;
   double R;	   
   R = 0.;
   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         {R += mat_X[i**maxcol+j]*mat_Y[i**maxcol+j];}
   return(R);
   }

/* Function to estimate statistical significance for Dietz's R statistic for two square matrices */

  void statsigR(double *X, double *Y, int *nrow, int *perm, double *pvRright, double *pvRleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, iter, *index, pointi, pointj;
   double *mat_X, *mat_Y, R, *matperm, *setn, Rsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Compute statistic value for the original sociomatrix */

   R = getR(mat_X,mat_Y,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices
   under the specified null hypothesis */

   *pvRright = 0.;
   *pvRleft = 0.;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow);
   for (i = 0; i< maxrow; i++){
      setn[i] = index[i];
   }
   for (i= 0; i< (maxrow-1); i++){
      pointi = setn[i];     
      for (j= i+1; j< maxcol; j++){
         pointj = setn[j];
         matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
         matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
      for (i = 0; i < maxrow; i++)
         for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}
   
   Rsim = getR(matperm,mat_Y,&maxrow,&maxcol);
   if (Rsim >= R)
        {*pvRright += 1.;}
      if (Rsim <= R)
        {*pvRleft += 1.;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvRright = (*pvRright+1.)/(*perm+1.);
   *pvRleft =(*pvRleft+1.)/(*perm+1.);
}

/* Function to estimate statistical significance for Dietz's R statistic for two rectangular matrices */

  void rectangularstatsigR(double *X, double *Y, int *nrow, int *ncol, int *perm, double *pvRright, double *pvRleft)
  {
/* Declaration of local variables */

   int i, j, k, m, maxrow, maxcol, iter, *index, point;
   double *mat_X, *mat_Y, R, *matperm, *setn, Rsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = *nrow;
   maxcol = *ncol;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
  {
  error("Null dimension");
  }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Compute statistic value for the original sociomatrix */

   R = getR(mat_X,mat_Y,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector setn and index */

   setn = calloc(maxrow*maxcol, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow*maxcol, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices 
   under the specified null hypothesis */

   *pvRright = 0.;
   *pvRleft = 0.;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow*maxcol; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow*maxcol);
   for (i = 0; i< maxrow*maxcol; i++){
      setn[i] = index[i];
   }
   k = 0.;
   for (i= 0; i< maxrow; i++){     
      for (j= 0; j< maxcol; j++){
         point = setn[k];
         matperm[i*maxcol+j] = mat_X[point];
         k += 1;}}
   Rsim = getR(matperm,mat_Y,&maxrow,&maxcol);
   if (Rsim >= R)
        {*pvRright += 1.;}
      if (Rsim <= R)
        {*pvRleft += 1.;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvRright = (*pvRright+1.)/(*perm+1.);
   *pvRleft =(*pvRleft+1.)/(*perm+1.);
  }

/* A function required in the exact permutation test */

  void permuteR(int *setn,double *matperm,double *mat_X, double *mat_Y, double *results, double R,
	       const int start, int maxrow, int maxcol)
  {
    int i, j, pointi, pointj;
    double Rsim;
    for (i= 0; i< (maxrow-1); i++){
       pointi = setn[i];     
       for (j= i+1; j< maxcol; j++){
          pointj = setn[j];
          matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
          matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
       for (i = 0; i < maxrow; i++)
          for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}
   
   Rsim = getR(matperm,mat_Y,&maxrow,&maxcol);
   if (Rsim >= R)
   {results[0] += 1.;}
   if (Rsim <= R)
   {results[1] += 1.;}
   if (start < maxrow) {
      int i, j;
      for (i = maxrow-2; i >= start; i--) {
         for (j = i + 1; j < maxrow; j++) {
	    swap(setn, i, j);
	    permuteR(setn,matperm,mat_X,mat_Y,results,R,i+1,maxrow,maxcol);
         }
      rotateLeft(setn, i, maxrow);
    }
    }
   }

/* Function to compute exact statistical significance for Dietz's R */

  void exactsigR(double *X, double *Y, int *nrow, double *pvRright, double *pvRleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, *setn, start, fact;
   double *mat_X, *mat_Y, R, *matperm, *results;

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Compute statistic value for the original sociomatrix */

   R = getR(mat_X,mat_Y,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector setn */

   setn = calloc(maxrow, sizeof(int));

   if (setn == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector results and compute factorial number of n */

   fact = 1;
   for (i = maxrow;i >= 1; i--){
      fact *= i;}

   results = calloc(2,sizeof(double));

   if (results == NULL)
   {
     error("Null dimension");
   }

/* Carry out the exact permutation procedure */

   for (i = 0; i < maxrow; i++){
      setn[i] = i;}

   for (i = 0; i < 2; i++){
      results[i] = 0.;}

   start =0;
   permuteR(setn,matperm,mat_X,mat_Y,results,R,start,maxrow,maxcol);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvRright = (results[0])/(fact);
   *pvRleft =(results[1])/(fact);

/* Deallocate all the matrices used in the routine */

   free(results);
   free(setn);
   free(matperm);
   free(mat_Y);
   free(mat_X);
}

/* Function to obtain Kendall's Kr statistic */

   double getKr(double *mat_X, double *mat_Y, int *maxrow, int *maxcol)
   {
   int i,j,k,sig;
   double Kr;	   
   Kr =0.;
   for (i = 0; i < *maxrow; i++)
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              sig = sign((mat_X[i**maxcol+j]-mat_X[i**maxcol+k])*(mat_Y[i**maxcol+j]-mat_Y[i**maxcol+k]));
              if ((i != j) & (i != k)) Kr = Kr + sig;
                else Kr = Kr;}}
   return(Kr);
   }

/* Function to estimate statistical significance for Kendall's Kr statistic for two square matrices */

  void statsigKr(double *X, double *Y, int *nrow, int *perm, double *pvKrright, double *pvKrleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, iter, *index, pointi, pointj;
   double *mat_X, *mat_Y, Kr, *matperm, *setn, Krsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Compute statistic value for the original sociomatrix */

   Kr = getKr(mat_X,mat_Y,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   } 
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices under the specified null hypothesis */

   *pvKrright = 0.;
   *pvKrleft = 0.;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow);
   for (i = 0; i< maxrow; i++){
      setn[i] = index[i];
   }
   for (i= 0; i< (maxrow-1); i++){
      pointi = setn[i];     
      for (j= i+1; j< maxcol; j++){
         pointj = setn[j];
         matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
         matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
      for (i = 0; i < maxrow; i++)
         for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}

   Krsim = getKr(matperm,mat_Y,&maxrow,&maxcol);
   if (Krsim >= Kr)
        {*pvKrright += 1.;}
      if (Krsim <= Kr)
        {*pvKrleft += 1.;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvKrright = (*pvKrright+1.)/(*perm+1.);
   *pvKrleft =(*pvKrleft+1.)/(*perm+1.);
  }

/* Function to estimate statistical significance for Kendall's Kr statistic for two rectangular matrices */

  void rectangularstatsigKr(double *X, double *Y, int *nrow, int *ncol, int *perm, double *pvKrright, double *pvKrleft)
  {
/* Declaration of local variables */

   int i, j, k, m, maxrow, maxcol, iter, *index, point;
   double *mat_X, *mat_Y, Kr, *matperm, *setn, Krsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Compute statistic value for the original sociomatrix */

   Kr = getKr(mat_X,mat_Y,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   } 
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow*maxcol, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow*maxcol, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices
   under the specified null hypothesis */

   *pvKrright = 0.;
   *pvKrleft = 0.;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow*maxcol; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow*maxcol);
   for (i = 0; i< maxrow*maxcol; i++){
      setn[i] = index[i];
   }

   k = 0.;
   for (i= 0; i< maxrow; i++){     
      for (j= 0; j< maxcol; j++){
         point = setn[k];
         matperm[i*maxcol+j] = mat_X[point];
         k += 1;}}
  
   Krsim = getKr(matperm,mat_Y,&maxrow,&maxcol);
   if (Krsim >= Kr)
        {*pvKrright += 1.;}
      if (Krsim <= Kr)
        {*pvKrleft += 1.;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvKrright = (*pvKrright+1.)/(*perm+1.);
   *pvKrleft =(*pvKrleft+1.)/(*perm+1.);
  }

/* A required function to compute exact statistical significance for Kendall's Kr */

   void permuteKr(int *setn,double *matperm,double *mat_X, double *mat_Y, double *results, double Kr,
		 const int start, int maxrow, int maxcol)
   {
    int i, j, pointi, pointj;
    double Krsim;
    for (i= 0; i< (maxrow-1); i++){
       pointi = setn[i];     
       for (j= i+1; j< maxcol; j++){
          pointj = setn[j];
          matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
          matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
       for (i = 0; i < maxrow; i++)
          for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}
   
   Krsim = getKr(matperm,mat_Y,&maxrow,&maxcol);
   if (Krsim >= Kr)
   {results[0] += 1.;}
   if (Krsim <= Kr)
   {results[1] += 1.;}
   if (start < maxrow) {
      int i, j;
      for (i = maxrow-2; i >= start; i--) {
         for (j = i + 1; j < maxrow; j++) {
	    swap(setn, i, j);
	    permuteKr(setn,matperm,mat_X,mat_Y,results,Kr,i+1,maxrow,maxcol);
         }
      rotateLeft(setn, i, maxrow);
    }
    }
   }

/* Function to compute exact statistical significance for Kendall's Kr */

  void exactsigKr(double *X, double *Y, int *nrow, double *pvKrright, double *pvKrleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, *setn, start, fact;
   double *mat_X, *mat_Y, Kr, *matperm, *results;

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Compute statistic value for the original sociomatrix */

   Kr = getKr(mat_X,mat_Y,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector setn */

   setn = calloc(maxrow, sizeof(int));

   if (setn == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector results and compute factorial number of n */

   fact = 1;
   for (i = maxrow;i >= 1; i--){
      fact *= i;}

   results = calloc(2,sizeof(double));

   if (results == NULL)
   {
     error("Null dimension");
   }

/* Carry out the exact permutation procedure */

   for (i = 0; i < maxrow; i++){
      setn[i] = i;}

   for (i = 0; i < 2; i++){
      results[i] = 0.;}

   start =0;
   permuteKr(setn,matperm,mat_X,mat_Y,results,Kr,start,maxrow,maxcol);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvKrright = (results[0])/(fact);
   *pvKrleft =(results[1])/(fact);

/* Deallocate all the matrices used in the routine */

   free(results);
   free(setn);
   free(matperm);
   free(mat_Y);
   free(mat_X);
  }

/* Function to obtain rowwise Mantel's Zr statistic */

   double getZr(double *mat_X, double *mat_Y, double *wi,double *ri, double *Xterm, double *Yterm, int *maxrow, int *maxcol)
   {
   int i,j,k;
   double Zr, cont;	   	   
   for (i = 0; i < *maxrow; i++){
      Xterm[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              if ((i != j) && (i != k)) {Xterm[i] += pow((mat_X[i**maxcol+j]-mat_X[i**maxcol+k]),2);}
            }
         }
      Yterm[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              if ((i != j) && (i != k)) {Yterm[i] += pow((mat_Y[i**maxcol+j]-mat_Y[i**maxcol+k]),2);}
            }   
         }
      wi[i] = sqrt(Xterm[i]*Yterm[i]);
   }
   Zr = 0.;
   for (i = 0; i < *maxrow; i++){
      ri[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
	      cont = wi[i];
	      if (cont != 0.){
                ri[i] += ((mat_X[i**maxcol+j]-mat_X[i**maxcol+k])*(mat_Y[i**maxcol+j]-mat_Y[i**maxcol+k])/wi[i]);}
	        else {ri[i] += 0.;}
              }
	 }
   Zr += wi[i]*ri[i];
   }
   return(Zr);
   }

/* Function to estimate statistical significance for rowwise Mantel's Zr statistic for square matrices */

  void statsigZr(double *X, double *Y, int *nrow, int *perm, double *pvZrright, double *pvZrleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, iter, *index, pointi, pointj;
   double *mat_X, *mat_Y, *wi, *ri, *Xterm, *Yterm, Zr, *matperm, *setn, Zrsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory size of vector wi, ri, Xterm and Yterm */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   ri = calloc(maxrow, sizeof(double));

   if (ri == NULL)
   {
     error("Null dimension");
   }

   Xterm = calloc(maxrow, sizeof(double));

   if (Xterm == NULL)
   {
     error("Null dimension");
   }
   
   Yterm = calloc(maxrow, sizeof(double));

   if (Yterm == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   Zr = getZr(mat_X,mat_Y,wi,ri,Xterm,Yterm,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   } 
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices under the specified null hypothesis */

   *pvZrright = 0.;
   *pvZrleft = 0.;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow);
   for (i = 0; i< maxrow; i++){
      setn[i] = index[i];
   }
   for (i= 0; i< (maxrow-1); i++){
      pointi = setn[i];     
      for (j= i+1; j< maxcol; j++){
         pointj = setn[j];
         matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
         matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
      for (i = 0; i < maxrow; i++)
         for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}

   Zrsim = getZr(matperm,mat_Y,wi,ri,Xterm,Yterm,&maxrow,&maxcol);
   if (Zrsim >= Zr)
        {*pvZrright += 1.;}
      if (Zrsim <= Zr)
        {*pvZrleft += 1.;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(Yterm);
   free(Xterm);
   free(ri);
   free(wi);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvZrright = (*pvZrright+1.)/(*perm+1.);
   *pvZrleft =(*pvZrleft+1.)/(*perm+1.);
}

/* Function to estimate statistical significance for rowwise Mantel's Zr statistic for rectangular matrices */

  void rectangularstatsigZr(double *X, double *Y, int *nrow, int *ncol, int *perm, double *pvZrright, double *pvZrleft)
  {
/* Declaration of local variables */

   int i, j, k, m, maxrow, maxcol, iter, *index, point;
   double *mat_X, *mat_Y, *wi, *ri, *Xterm, *Yterm, Zr, *matperm, *setn, Zrsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = *nrow;
   maxcol = *ncol;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory size of vector wi, ri, Xterm and Yterm */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   ri = calloc(maxrow, sizeof(double));

   if (ri == NULL)
   {
     error("Null dimension");
   }

   Xterm = calloc(maxrow, sizeof(double));

   if (Xterm == NULL)
   {
     error("Null dimension");
   }
   
   Yterm = calloc(maxrow, sizeof(double));

   if (Yterm == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   Zr = getZr(mat_X,mat_Y,wi,ri,Xterm,Yterm,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector setn and index */

   setn = calloc(maxrow*maxcol, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow*maxcol, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices under the specified null hypothesis */

   *pvZrright = 0.;
   *pvZrleft = 0.;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow*maxcol; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow*maxcol);
   for (i = 0; i< maxrow*maxcol; i++){
      setn[i] = index[i];
   }

   k = 0.;
   for (i= 0; i< maxrow; i++){     
      for (j= 0; j< maxcol; j++){
         point = setn[k];
         matperm[i*maxcol+j] = mat_X[point];
         k += 1;}}
  
   Zrsim = getZr(matperm,mat_Y,wi,ri,Xterm,Yterm,&maxrow,&maxcol);

   if (Zrsim >= Zr)
        {*pvZrright += 1.;}
      if (Zrsim <= Zr)
        {*pvZrleft += 1.;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(Yterm);
   free(Xterm);
   free(ri);
   free(wi);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvZrright = (*pvZrright+1.)/(*perm+1.);
   *pvZrleft =(*pvZrleft+1.)/(*perm+1.);
}

/* A required function to compute exact p value for Mantel's Zr statistic */

   void permuteZr(int *setn,double *matperm,double *mat_X, double *mat_Y, double *wi, double *ri, double *Xterm,
	       double *Yterm, double *results, double Zr, const int start, int maxrow, int maxcol)
   {
    int i, j, pointi, pointj;
    double Zrsim;
    for (i= 0; i< (maxrow-1); i++){
       pointi = setn[i];     
       for (j= i+1; j< maxcol; j++){
          pointj = setn[j];
          matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
          matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
       for (i = 0; i < maxrow; i++)
          for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}
   
   Zrsim = getZr(matperm,mat_Y,wi,ri,Xterm,Yterm,&maxrow,&maxcol);
   if (Zrsim >= Zr)
   {results[0] += 1.;}
   if (Zrsim <= Zr)
   {results[1] += 1.;}
   if (start < maxrow) {
      int i, j;
      for (i = maxrow-2; i >= start; i--) {
         for (j = i + 1; j < maxrow; j++) {
	    swap(setn, i, j);
	    permuteZr(setn,matperm,mat_X,mat_Y,wi,ri,Xterm,Yterm,results,Zr,i+1,maxrow,maxcol);
         }
      rotateLeft(setn, i, maxrow);
    }
    }
   }

/* Function to estimate exact statistical significance for rowwise Mantel's Zr statistic */

  void exactsigZr(double *X, double *Y, int *nrow, double *pvZrright, double *pvZrleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, *setn, start, fact;
   double *mat_X, *mat_Y, *wi, *ri, *Xterm, *Yterm, Zr, *matperm, *results;

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory size of vector wi, ri, Xterm and Yterm */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   ri = calloc(maxrow, sizeof(double));

   if (ri == NULL)
   {
     error("Null dimension");
   }

   Xterm = calloc(maxrow, sizeof(double));

   if (Xterm == NULL)
   {
     error("Null dimension");
   }
   
   Yterm = calloc(maxrow, sizeof(double));

   if (Yterm == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   Zr = getZr(mat_X,mat_Y,wi,ri,Xterm,Yterm,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector setn */

   setn = calloc(maxrow, sizeof(int));

   if (setn == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector results and compute factorial number of n */

   fact = 1;
   for (i = maxrow;i >= 1; i--){
      fact *= i;}

   results = calloc(2,sizeof(double));

   if (results == NULL)
   {
     error("Null dimension");
   }

/* Carry out the exact permutation procedure */

   for (i = 0; i < maxrow; i++){
      setn[i] = i;}

   for (i = 0; i < 2; i++){
      results[i] = 0.;}

   start =0;
   permuteZr(setn,matperm,mat_X,mat_Y,wi,ri,Xterm,Yterm,results,Zr,start,maxrow,maxcol);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvZrright = (results[0])/(fact);
   *pvZrleft =(results[1])/(fact);

/* Deallocate all the matrices used in the routine */

   free(results);
   free(setn);
   free(matperm);
   free(Yterm);
   free(Xterm);
   free(ri);
   free(wi);   
   free(mat_Y);
   free(mat_X);
}

/* Function to obtain rowwise Dietz's Rr statistic */

   double getRr(double *mat_X, double *mat_Y, double *wi,double *ri, double *Xterm, double *Yterm, int *maxrow, int *maxcol)
   {
   int i,j,k;
   double Rr,cont;	   	   
   for (i = 0; i < *maxrow; i++){
      Xterm[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              if ((i != j) && (i != k)) {Xterm[i] += pow((mat_X[i**maxcol+j]-mat_X[i**maxcol+k]),2);}
            }
         }
      Yterm[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              if ((i != j) && (i != k)) {Yterm[i] += pow((mat_Y[i**maxcol+j]-mat_Y[i**maxcol+k]),2);}
            }   
         }
      wi[i] = sqrt(Xterm[i]*Yterm[i]);
   }
   Rr = 0.;
   for (i = 0; i < *maxrow; i++){
      ri[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
	      cont = wi[i];
	      if (cont != 0.){
                ri[i] += ((mat_X[i**maxcol+j]-mat_X[i**maxcol+k])*(mat_Y[i**maxcol+j]-mat_Y[i**maxcol+k])/wi[i]);}
	        else {ri[i] += 0.;}
              }
      }
   Rr += wi[i]*ri[i];
   }
   return(Rr);
   }

/* Function to estimate statistical significance for rowwise Dietz's Rr statistic for square matrices */

  void statsigRr(double *X, double *Y, int *nrow, int *perm, double *pvRrright, double *pvRrleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, iter, *index, pointi, pointj;
   double *mat_X, *mat_Y, *wi, *ri, *Xterm, *Yterm, Rr, *matperm, *setn, Rrsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory size of vector wi, ri, Xterm and Yterm */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   ri = calloc(maxrow, sizeof(double));

   if (ri == NULL)
   {
     error("Null dimension");
   }

   Xterm = calloc(maxrow, sizeof(double));

   if (Xterm == NULL)
   {
     error("Null dimension");
   }
   
   Yterm = calloc(maxrow, sizeof(double));

   if (Yterm == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   Rr = getRr(mat_X,mat_Y,wi,ri,Xterm,Yterm,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   } 
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices 
   under the specified null hypothesis */

   *pvRrright = 0.;
   *pvRrleft = 0.;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow);
   for (i = 0; i< maxrow; i++){
      setn[i] = index[i];
   }
   for (i= 0; i< (maxrow-1); i++){
      pointi = setn[i];     
      for (j= i+1; j< maxcol; j++){
         pointj = setn[j];
         matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
         matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
      for (i = 0; i < maxrow; i++)
         for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}

   Rrsim = getRr(matperm,mat_Y,wi,ri,Xterm,Yterm,&maxrow,&maxcol);
   if (Rrsim >= Rr)
        {*pvRrright += 1.;}
      if (Rrsim <= Rr)
        {*pvRrleft += 1.;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(Yterm);
   free(Xterm);
   free(ri);
   free(wi);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvRrright = (*pvRrright+1.)/(*perm+1.);
   *pvRrleft =(*pvRrleft+1.)/(*perm+1.);
  }

/* Function to estimate statistical significance for rowwise Dietz's Rr statistic for rectangular matrices */

  void rectangularstatsigRr(double *X, double *Y, int *nrow, int *ncol, int *perm, double *pvRrright, double *pvRrleft)
  {
/* Declaration of local variables */

   int i, j, k, m, maxrow, maxcol, iter, *index, point;
   double *mat_X, *mat_Y, *wi, *ri, *Xterm, *Yterm, Rr, *matperm, *setn, Rrsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = *nrow;
   maxcol = *ncol;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory size of vector wi, ri, Xterm and Yterm */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   ri = calloc(maxrow, sizeof(double));

   if (ri == NULL)
   {
     error("Null dimension");
   }

   Xterm = calloc(maxrow, sizeof(double));

   if (Xterm == NULL)
   {
     error("Null dimension");
   }
   
   Yterm = calloc(maxrow, sizeof(double));

   if (Yterm == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   Rr = getRr(mat_X,mat_Y,wi,ri,Xterm,Yterm,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   } 
   
/* Allocate in memory size of vector setn and index */

   setn = calloc(maxrow*maxcol, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow*maxcol, sizeof(int));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices under the specified null hypothesis */

   *pvRrright = 0.;
   *pvRrleft = 0.;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow*maxcol; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow*maxcol);
   for (i = 0; i< maxrow*maxcol; i++){
      setn[i] = index[i];
   }

   k = 0.;
   for (i= 0; i< maxrow; i++){     
      for (j= 0; j< maxcol; j++){
         point = setn[k];
         matperm[i*maxcol+j] = mat_X[point];
         k += 1;}}

   Rrsim = getRr(matperm,mat_Y,wi,ri,Xterm,Yterm,&maxrow,&maxcol);
   if (Rrsim >= Rr)
        {*pvRrright += 1.;}
      if (Rrsim <= Rr)
        {*pvRrleft += 1.;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(Yterm);
   free(Xterm);
   free(ri);
   free(wi);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvRrright = (*pvRrright+1.)/(*perm+1.);
   *pvRrleft =(*pvRrleft+1.)/(*perm+1.);
}

/* A required function to estimate exact statistical significance for rowwise Mantel's Zr statistic */

  void permuteRr(int *setn,double *matperm,double *mat_X, double *mat_Y, double *wi, double *ri, double *Xterm,
   	        double *Yterm, double *results, double Rr, const int start, int maxrow, int maxcol)
   {
    int i, j, pointi, pointj;
    double Rrsim;
    for (i= 0; i< (maxrow-1); i++){
       pointi = setn[i];     
       for (j= i+1; j< maxcol; j++){
          pointj = setn[j];
          matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
          matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
       for (i = 0; i < maxrow; i++)
          for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}
   
   Rrsim = getRr(matperm,mat_Y,wi,ri,Xterm,Yterm,&maxrow,&maxcol);
   if (Rrsim >= Rr)
   {results[0] += 1.;}
   if (Rrsim <= Rr)
   {results[1] += 1.;}
   if (start < maxrow) {
      int i, j;
      for (i = maxrow-2; i >= start; i--) {
         for (j = i + 1; j < maxrow; j++) {
	    swap(setn, i, j);
	    permuteRr(setn,matperm,mat_X,mat_Y,wi,ri,Xterm,Yterm,results,Rr,i+1,maxrow,maxcol);
         }
      rotateLeft(setn, i, maxrow);
    }
    }
   }


/* Function to estimate exact statistical significance for rowwise Dietz's Zr statistic */

  void exactsigRr(double *X, double *Y, int *nrow, double *pvRrright, double *pvRrleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, *setn, start, fact;
   double *mat_X, *mat_Y, *wi, *ri, *Xterm, *Yterm, Rr, *results, *matperm;

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory size of vector wi, ri, Xterm and Yterm */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   ri = calloc(maxrow, sizeof(double));

   if (ri == NULL)
   {
     error("Null dimension");
   }

   Xterm = calloc(maxrow, sizeof(double));

   if (Xterm == NULL)
   {
     error("Null dimension");
   }
   
   Yterm = calloc(maxrow, sizeof(double));

   if (Yterm == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   Rr = getRr(mat_X,mat_Y,wi,ri,Xterm,Yterm,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector setn */

   setn = calloc(maxrow, sizeof(int));

   if (setn == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector results and compute factorial number of n */

   fact = 1;
   for (i = maxrow;i >= 1; i--){
      fact *= i;}

   results = calloc(2,sizeof(double));

   if (results == NULL)
   {
     error("Null dimension");
   }

/* Carry out the exact permutation procedure */

   for (i = 0; i < maxrow; i++){
      setn[i] = i;}

   for (i = 0; i < 2; i++){
      results[i] = 0.;}

   start =0;
   permuteRr(setn,matperm,mat_X,mat_Y,wi,ri,Xterm,Yterm,results,Rr,start,maxrow,maxcol);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvRrright = (results[0])/(fact);
   *pvRrleft =(results[1])/(fact);

/* Deallocate all the matrices used in the routine */

   free(results);
   free(setn);
   free(matperm);
   free(Yterm);
   free(Xterm);
   free(ri);
   free(wi);   
   free(mat_Y);
   free(mat_X);
}


/* Function to obtain rowwise Pearson statistic */

   double getrrw(double *mat_X, double *mat_Y, double *wi,double *ri, double *X1term, double *X2term, int *maxrow, int *maxcol)
   {
   int i,j,k;
   double Zr, sumX1term, sumX2term, cont, term, rrw;
   for (i = 0; i < *maxrow; i++){
      X1term[i] = 0.;
      X2term[i] = 0.;
      wi[i] = 0.;
      ri[i] = 0.;}
         
   for (i = 0; i < *maxrow; i++){
      X1term[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              X1term[i] += pow((mat_X[i**maxcol+j]-mat_X[i**maxcol+k]),2);}
         }
      X2term[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              X2term[i] += pow((mat_Y[i**maxcol+j]-mat_Y[i**maxcol+k]),2);}  
         }
      wi[i] = sqrt(X1term[i]*X2term[i]);
   }
   Zr = 0.;
   for (i = 0; i < *maxrow; i++){
      ri[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
	      cont = wi[i];
	      if (cont != 0.){
                ri[i] += ((mat_X[i**maxcol+j]-mat_X[i**maxcol+k])*(mat_Y[i**maxcol+j]-mat_Y[i**maxcol+k])/wi[i]);}
	        else {ri[i] += 0.;}
              }
   }	    
   Zr += wi[i]*ri[i];
   }
   sumX1term = 0.;
   sumX2term = 0.;
   for (i = 0; i < *maxrow; i++){
      sumX1term += X1term[i];
      sumX2term += X2term[i];
   }
   term = sqrt(sumX1term*sumX2term);
   if (term == 0.) {rrw = 0.;}
     else {rrw = Zr/term;}     
   return(rrw);
   }

/* Function to obtain partial rowwise Pearson statistic */

   double getrrwxyz(double *mat_X, double *mat_Y, double *mat_Z, double *wi,double *ri, double *X1term, double *X2term, int *maxrow, int *maxcol)
   {
   double rrwxy, rrwxz, rrwyz, rrwxyz, term;	   	   
   rrwxy = getrrw(mat_X,mat_Y,wi,ri,X1term,X2term,maxrow,maxcol);
   rrwxz = getrrw(mat_X,mat_Z,wi,ri,X1term,X2term,maxrow,maxcol);
   rrwyz = getrrw(mat_Y,mat_Z,wi,ri,X1term,X2term,maxrow,maxcol);
   term = sqrt((1-pow(rrwxz,2))*(1-pow(rrwyz,2)));
   if (term == 0.){rrwxyz = -5.;}
     else {
         rrwxyz = (rrwxy - rrwxz*rrwyz)/term;
     }           
   return(rrwxyz);
   }

/* Function to estimate statistical significance for rowwise Mantel's Z statistic for square matrices */

  void partialstatsigZr(double *X, double *Y, double *Z, int *nrow, int *perm, double *pvpartialZrright, double *pvpartialZrleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, iter, *index, pointi, pointj, valid;
   double *mat_X, *mat_Y, *mat_Z, *wi, *ri, *X1term, *X2term, rrwxyz, *matperm, *setn, rrwxyzsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory the size of original matrix Z, exit if error */

   mat_Z = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Z == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Z - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_Z[i*maxcol+j] = Z[m];
         m++;
   }

/* Allocate in memory size of vector wi, ri, X1term and X2term */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   ri = calloc(maxrow, sizeof(double));

   if (ri == NULL)
   {
     error("Null dimension");
   }

   X1term = calloc(maxrow, sizeof(double));

   if (X1term == NULL)
   {
     error("Null dimension");
   }
   
   X2term = calloc(maxrow, sizeof(double));

   if (X2term == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   rrwxyz = getrrwxyz(mat_X,mat_Y,mat_Z,wi,ri,X1term,X2term,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices under the specified null hypothesis */

   *pvpartialZrright = 0.;
   *pvpartialZrleft = 0.;
   valid = *perm;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow);
   for (i = 0; i< maxrow; i++){
      setn[i] = index[i];
   }
   for (i= 0; i< (maxrow-1); i++){
      pointi = setn[i];     
      for (j= i+1; j< maxcol; j++){
         pointj = setn[j];
         matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
         matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
      for (i = 0; i < maxrow; i++)
         for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}

   rrwxyzsim = getrrwxyz(matperm,mat_Y,mat_Z,wi,ri,X1term,X2term,&maxrow,&maxcol);
   if (rrwxyzsim != -5.){
   if (rrwxyzsim >= rrwxyz)
        {*pvpartialZrright += 1.;}
      if (rrwxyzsim <= rrwxyz)
        {*pvpartialZrleft += 1.;}
   }
     else { valid -= 1;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(X2term);
   free(X1term);
   free(ri);
   free(wi);
   free(mat_Z);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvpartialZrright = (*pvpartialZrright+1.)/(valid+1.);
   *pvpartialZrleft =(*pvpartialZrleft+1.)/(valid+1.);
}

/* Function to estimate statistical significance for rowwise Mantel's Z statistic for rectangular matrices */

  void rectangularstatsigpartZr(double *X, double *Y, double *Z, int *nrow, int *ncol, int *perm, double *pvpartialZrright, double *pvpartialZrleft)
  {
/* Declaration of local variables */

   int i, j, k, m, maxrow, maxcol, iter, *index, point, valid;
   double *mat_X, *mat_Y, *mat_Z, *wi, *ri, *X1term, *X2term, rrwxyz, *matperm, *setn, rrwxyzsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = *nrow;
   maxcol = *ncol;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory the size of original matrix Z, exit if error */

   mat_Z = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Z == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Z - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_Z[i*maxcol+j] = Z[m];
         m++;
   }

/* Allocate in memory size of vector wi, ri, X1term and X2term */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   ri = calloc(maxrow, sizeof(double));

   if (ri == NULL)
   {
     error("Null dimension");
   }

   X1term = calloc(maxrow, sizeof(double));

   if (X1term == NULL)
   {
     error("Null dimension");
   }
   
   X2term = calloc(maxrow, sizeof(double));

   if (X2term == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   rrwxyz = getrrwxyz(mat_X,mat_Y,mat_Z,wi,ri,X1term,X2term,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow*maxcol, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow*maxcol, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices under the specified null hypothesis */

   *pvpartialZrright = 0.;
   *pvpartialZrleft = 0.;
   valid = *perm;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow*maxcol; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow*maxcol);
   for (i = 0; i< maxrow*maxcol; i++){
      setn[i] = index[i];
   }
   k = 0;
   for (i= 0; i< maxrow; i++){     
      for (j= 0; j< maxcol; j++){
         point = setn[k];
         matperm[i*maxcol+j] = mat_X[point];
         k += 1;}}

   rrwxyzsim = getrrwxyz(matperm,mat_Y,mat_Z,wi,ri,X1term,X2term,&maxrow,&maxcol);
   if (rrwxyzsim != -5.){
   if (rrwxyzsim >= rrwxyz)
        {*pvpartialZrright += 1.;}
      if (rrwxyzsim <= rrwxyz)
        {*pvpartialZrleft += 1.;}
   }
     else { valid -= 1;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(X2term);
   free(X1term);
   free(ri);
   free(wi);
   free(mat_Z);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvpartialZrright = (*pvpartialZrright+1.)/(valid+1.);
   *pvpartialZrleft =(*pvpartialZrleft+1.)/(valid+1.);
}

/* A required function to estimate exact statistical significance for the partial rowwise Mantel's Zr statistic */

   void permuterrwxyz(int *setn,double *matperm,double *mat_X, double *mat_Y, double *mat_Z, double *wi, double *ri, double *X1term,
	       double *X2term, double *results, double rrwxyz, const int start, int maxrow, int maxcol)
   {
    int i, j, pointi, pointj;
    double rrwxyzsim;
    for (i= 0; i< (maxrow-1); i++){
       pointi = setn[i];     
       for (j= i+1; j< maxcol; j++){
          pointj = setn[j];
          matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
          matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
       for (i = 0; i < maxrow; i++)
          for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}
   
   rrwxyzsim = getrrwxyz(matperm,mat_Y,mat_Z,wi,ri,X1term,X2term,&maxrow,&maxcol);
   if (rrwxyzsim != -5){
   if (rrwxyzsim >= rrwxyz)
   {results[0] += 1.;}
   if (rrwxyzsim <= rrwxyz)
   {results[1] += 1.;}
   }
   else {results[2] += 1.;}
   if (start < maxrow) {
      int i, j;
      for (i = maxrow-2; i >= start; i--) {
         for (j = i + 1; j < maxrow; j++) {
	    swap(setn, i, j);
	    permuterrwxyz(setn,matperm,mat_X,mat_Y,mat_Z,wi,ri,X1term,X2term,results,rrwxyz,i+1,maxrow,maxcol);
         }
      rotateLeft(setn, i, maxrow);
    }
    }
   }

/* Function to estimate exact statistical significance for the partial rowwise Mantel's Zr statistic */

  void partialexactsigZr(double *X, double *Y, double *Z, int *nrow, double *pvpartialZrright, double *pvpartialZrleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, *setn, start, fact;
   double *mat_X, *mat_Y, *mat_Z, *wi, *ri, *X1term, *X2term, rrwxyz, *matperm, *results;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory the size of original matrix Z, exit if error */

   mat_Z = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Z == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Z - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_Z[i*maxcol+j] = Z[m];
         m++;
   }

/* Allocate in memory size of vector wi, ri, X1term and X2term */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   ri = calloc(maxrow, sizeof(double));

   if (ri == NULL)
   {
     error("Null dimension");
   }

   X1term = calloc(maxrow, sizeof(double));

   if (X1term == NULL)
   {
     error("Null dimension");
   }
   
   X2term = calloc(maxrow, sizeof(double));

   if (X2term == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   rrwxyz = getrrwxyz(mat_X,mat_Y,mat_Z,wi,ri,X1term,X2term,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow, sizeof(int));

   if (setn == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector results and compute factorial number of n */

   fact = 1;
   for (i = maxrow;i >= 1; i--){
      fact *= i;}

   results = calloc(3,sizeof(double));

   if (results == NULL)
   {
     error("Null dimension");
   }

/* Carry out the exact permutation procedure */

   for (i = 0; i < maxrow; i++){
      setn[i] = i;}

   for (i = 0; i < 3; i++){
      results[i] = 0.;}

   start =0;
   permuterrwxyz(setn,matperm,mat_X,mat_Y,mat_Z,wi,ri,X1term,X2term,results,rrwxyz,start,maxrow,maxcol);

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(results);
   free(setn);
   free(matperm);
   free(X2term);
   free(X1term);
   free(ri);
   free(wi);
   free(mat_Z);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   fact -= results[2];
   *pvpartialZrright = (results[0])/(fact);
   *pvpartialZrleft =(results[1])/(fact);
}

/* Function to obtain rowwise Spearman statistic */

   double getrhorw(double *mat_X, double *mat_Y, double *wi,double *ri, double *X1term, double *X2term, int *maxrow, int *maxcol)
   {
   int i,j,k;
   double Rr, sumX1term, sumX2term, cont, term, rhorw;
   for (i = 0; i < *maxrow; i++){
      X1term[i] = 0.;
      X2term[i] = 0.;
      wi[i] = 0.;
      ri[i] = 0.;}
         
   for (i = 0; i < *maxrow; i++){
      X1term[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              if ((i != j) && (i != k)) {X1term[i] += pow((mat_X[i**maxcol+j]-mat_X[i**maxcol+k]),2);}
            }
         }
      X2term[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              if ((i != j) && (i != k)) {X2term[i] += pow((mat_Y[i**maxcol+j]-mat_Y[i**maxcol+k]),2);}
            }   
         }
      wi[i] = sqrt(X1term[i]*X2term[i]);
   }
   Rr = 0.;
   for (i = 0; i < *maxrow; i++){
      ri[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              if ((i != j) && (i != k)){
		cont = wi[i]; 
	        if (cont != 0.)
	        {ri[i] += ((mat_X[i**maxcol+j]-mat_X[i**maxcol+k])*(mat_Y[i**maxcol+j]-mat_Y[i**maxcol+k])/wi[i]);}
		  else {ri[i] += 0.;}
	      }
            }
      }
   Rr += wi[i]*ri[i];
   }
   sumX1term = 0.;
   sumX2term = 0.;
   for (i = 0; i < *maxrow; i++){
      sumX1term += X1term[i];
      sumX2term += X2term[i];
   }
   term = sqrt(sumX1term*sumX2term);
   if (term != 0.){
     rhorw = Rr/term;}
     else {rhorw = 0.;}
   return(rhorw);
   }

/* Function to obtain partial rowwise Spearman statistic */

   double getrhorwxyz(double *mat_X, double *mat_Y, double *mat_Z, double *wi,double *ri, double *X1term, double *X2term, int *maxrow, int *maxcol)
   {
   double rhorwxy, rhorwxz, rhorwyz, rhorwxyz, term;	   	   
   rhorwxy = getrhorw(mat_X,mat_Y,wi,ri,X1term,X2term,maxrow,maxcol);
   rhorwxz = getrhorw(mat_X,mat_Z,wi,ri,X1term,X2term,maxrow,maxcol);
   rhorwyz = getrhorw(mat_Y,mat_Z,wi,ri,X1term,X2term,maxrow,maxcol);
   term = sqrt((1-pow(rhorwxz,2))*(1-pow(rhorwyz,2)));
   if (term == 0.){rhorwxyz = -5.;}
     else {
         rhorwxyz = (rhorwxy - rhorwxz*rhorwyz)/term;
     } 
   return(rhorwxyz);
   }

/* Function to estimate statistical significance for rowwise Dietz's Rr statistic for square matrices */

void partialstatsigRr(double *X, double *Y, double *Z, int *nrow, int *perm, double *pvpartialRrright, double *pvpartialRrleft)
{
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, iter, *index, pointi, pointj, valid;
   double *mat_X, *mat_Y, *mat_Z, *wi, *ri, *X1term, *X2term, rhorwxyz, *matperm, *setn, rhorwxyzsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory the size of original matrix Z, exit if error */

   mat_Z = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Z == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Z - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_Z[i*maxcol+j] = Z[m];
         m++;
   }

/* Allocate in memory size of vector wi, ri, X1term and X2term */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   ri = calloc(maxrow, sizeof(double));

   if (ri == NULL)
   {
     error("Null dimension");
   }

   X1term = calloc(maxrow, sizeof(double));

   if (X1term == NULL)
   {
     error("Null dimension");
   }
   
   X2term = calloc(maxrow, sizeof(double));

   if (X2term == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   rhorwxyz = getrhorwxyz(mat_X,mat_Y,mat_Z,wi,ri,X1term,X2term,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices under the specified null hypothesis */

   *pvpartialRrright = 0.;
   *pvpartialRrleft = 0.;
   valid = *perm;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow);
   for (i = 0; i< maxrow; i++){
      setn[i] = index[i];
   }
   for (i= 0; i< (maxrow-1); i++){
      pointi = setn[i];     
      for (j= i+1; j< maxcol; j++){
         pointj = setn[j];
         matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
         matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
      for (i = 0; i < maxrow; i++)
         for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}

   rhorwxyzsim = getrhorwxyz(matperm,mat_Y,mat_Z,wi,ri,X1term,X2term,&maxrow,&maxcol);
   if (rhorwxyzsim != -5.){
   if (rhorwxyzsim >= rhorwxyz)
        {*pvpartialRrright += 1.;}
      if (rhorwxyzsim <= rhorwxyz)
        {*pvpartialRrleft += 1.;}
   }
     else { valid -= 1;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(X2term);
   free(X1term);
   free(ri);
   free(wi);
   free(mat_Z);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvpartialRrright = (*pvpartialRrright+1.)/(valid+1.);
   *pvpartialRrleft =(*pvpartialRrleft+1.)/(valid+1.);
}

/* Function to estimate statistical significance for rowwise Dietz's Rr statistic for rectangular matrices */

void rectangularstatsigpartRr(double *X, double *Y, double *Z, int *nrow, int *ncol, int *perm, double *pvpartialRrright, double *pvpartialRrleft)
{
/* Declaration of local variables */

   int i, j, k, m, maxrow, maxcol, iter, *index, point, valid;
   double *mat_X, *mat_Y, *mat_Z, *wi, *ri, *X1term, *X2term, rhorwxyz, *matperm, *setn, rhorwxyzsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = *nrow;
   maxcol = *ncol;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory the size of original matrix Z, exit if error */

   mat_Z = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Z == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Z - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_Z[i*maxcol+j] = Z[m];
         m++;
   }

/* Allocate in memory size of vector wi, ri, X1term and X2term */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   ri = calloc(maxrow, sizeof(double));

   if (ri == NULL)
   {
     error("Null dimension");
   }

   X1term = calloc(maxrow, sizeof(double));

   if (X1term == NULL)
   {
     error("Null dimension");
   }
   
   X2term = calloc(maxrow, sizeof(double));

   if (X2term == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   rhorwxyz = getrhorwxyz(mat_X,mat_Y,mat_Z,wi,ri,X1term,X2term,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow*maxcol, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow*maxcol, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices under the specified null hypothesis */

   *pvpartialRrright = 0.;
   *pvpartialRrleft = 0.;
   valid = *perm;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow*maxcol; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow*maxcol);
   for (i = 0; i< maxrow*maxcol; i++){
      setn[i] = index[i];
   }
   k = 0;
   for (i= 0; i< maxrow; i++){     
      for (j= 0; j< maxcol; j++){
         point = setn[k];
         matperm[i*maxcol+j] = mat_X[point];
         k += 1;}}

   rhorwxyzsim = getrhorwxyz(matperm,mat_Y,mat_Z,wi,ri,X1term,X2term,&maxrow,&maxcol);
   if (rhorwxyzsim != -5.){
   if (rhorwxyzsim >= rhorwxyz)
        {*pvpartialRrright += 1.;}
      if (rhorwxyzsim <= rhorwxyz)
        {*pvpartialRrleft += 1.;}
   }
     else { valid -= 1;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(X2term);
   free(X1term);
   free(ri);
   free(wi);
   free(mat_Z);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvpartialRrright = (*pvpartialRrright+1.)/(valid+1.);
   *pvpartialRrleft =(*pvpartialRrleft+1.)/(valid+1.);
}

/* A required function to estimate exact statistical significance for partial rowwise Dietz's Rr statistic */

   void permuterhorwxyz(int *setn,double *matperm,double *mat_X, double *mat_Y, double *mat_Z, double *wi, double *ri, double *X1term,
	       double *X2term, double *results, double rhorwxyz, const int start, int maxrow, int maxcol)
   {
    int i, j, pointi, pointj;
    double rhorwxyzsim;
    for (i= 0; i< (maxrow-1); i++){
       pointi = setn[i];     
       for (j= i+1; j< maxcol; j++){
          pointj = setn[j];
          matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
          matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
       for (i = 0; i < maxrow; i++)
          for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}
   
   rhorwxyzsim = getrhorwxyz(matperm,mat_Y,mat_Z,wi,ri,X1term,X2term,&maxrow,&maxcol);
   if (rhorwxyzsim != -5){
   if (rhorwxyzsim >= rhorwxyz)
   {results[0] += 1.;}
   if (rhorwxyzsim <= rhorwxyz)
   {results[1] += 1.;}
   }
   else {results[2] += 1.;}
   if (start < maxrow) {
      int i, j;
      for (i = maxrow-2; i >= start; i--) {
         for (j = i + 1; j < maxrow; j++) {
	    swap(setn, i, j);
	    permuterhorwxyz(setn,matperm,mat_X,mat_Y,mat_Z,wi,ri,X1term,X2term,results,rhorwxyz,i+1,maxrow,maxcol);
         }
      rotateLeft(setn, i, maxrow);
    }
    }
   }

/* A function to estimate exact statistical significance for partial rowwise Dietz's Rr statistic */

  void partialexactsigRr(double *X, double *Y, double *Z, int *nrow, double *pvpartialRrright, double *pvpartialRrleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, *setn, start, fact;
   double *mat_X, *mat_Y, *mat_Z, *wi, *ri, *X1term, *X2term, rhorwxyz, *matperm, *results;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory the size of original matrix Z, exit if error */

   mat_Z = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Z == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Z - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_Z[i*maxcol+j] = Z[m];
         m++;
   }

/* Allocate in memory size of vector wi, ri, X1term and X2term */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   ri = calloc(maxrow, sizeof(double));

   if (ri == NULL)
   {
     error("Null dimension");
   }

   X1term = calloc(maxrow, sizeof(double));

   if (X1term == NULL)
   {
     error("Null dimension");
   }
   
   X2term = calloc(maxrow, sizeof(double));

   if (X2term == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   rhorwxyz = getrhorwxyz(mat_X,mat_Y,mat_Z,wi,ri,X1term,X2term,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow, sizeof(int));

   if (setn == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector results and compute factorial number of n */

   fact = 1;
   for (i = maxrow;i >= 1; i--){
      fact *= i;}

   results = calloc(3,sizeof(double));

   if (results == NULL)
   {
     error("Null dimension");
   }

/* Carry out the exact permutation procedure */

   for (i = 0; i < maxrow; i++){
      setn[i] = i;}

   for (i = 0; i < 3; i++){
      results[i] = 0.;}

   start =0;
   permuterhorwxyz(setn,matperm,mat_X,mat_Y,mat_Z,wi,ri,X1term,X2term,results,rhorwxyz,start,maxrow,maxcol);

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(results);
   free(setn);
   free(matperm);
   free(X2term);
   free(X1term);
   free(ri);
   free(wi);
   free(mat_Z);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   fact -= results[2];
   *pvpartialRrright = (results[0])/(fact);
   *pvpartialRrleft =(results[1])/(fact);
}

/* Function to obtain rowwise Kendall statistic */

   double gettaurw(double *mat_X, double *mat_Y, double *wi,double *taui, double *X1term, double *X2term, int *maxrow, int *maxcol)
   {
   int i,j,k;
   double Kr, sumX1term, sumX2term, cont, term, taurw;
   for (i = 0; i < *maxrow; i++){
      X1term[i] = 0.;
      X2term[i] = 0.;
      wi[i] = 0.;
      taui[i] = 0.;}
         
   for (i = 0; i < *maxrow; i++){
      X1term[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              if ((i != j) && (i != k)) {X1term[i] += sign(pow(mat_X[i**maxcol+j]-mat_X[i**maxcol+k],2));}
            }
         }
      X2term[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              if ((i != j) && (i != k)) {X2term[i] += sign(pow(mat_Y[i**maxcol+j]-mat_Y[i**maxcol+k],2));}
            }   
         }
      wi[i] = sqrt(X1term[i]*X2term[i]);
   }
   Kr = 0.;
   for (i = 0; i < *maxrow; i++){
      taui[i] = 0.;
      for (j = 0; j < *maxcol; j++)
         for (k = 0; k < *maxcol; k++){
            if (j > k){
              if ((i != j) && (i != k)){
		cont = wi[i]; 
	        if (cont != 0.)
	        {taui[i] += (sign(mat_X[i**maxcol+j]-mat_X[i**maxcol+k])*sign(mat_Y[i**maxcol+j]-mat_Y[i**maxcol+k])/wi[i]);}
		  else {taui[i] += 0.;}
	      }
            }
      }
   Kr += wi[i]*taui[i];
   }
   sumX1term = 0.;
   sumX2term = 0.;
   for (i = 0; i < *maxrow; i++){
      sumX1term += X1term[i];
      sumX2term += X2term[i];
   }
   term = sqrt(sumX1term*sumX2term);
   if (term != 0.){
     taurw = Kr/term;}
     else {taurw = 0.;}
   return(taurw);
   }

/* Function to obtain partial rowwise Kendall statistic */

   double gettaurwxyz(double *mat_X, double *mat_Y, double *mat_Z, double *wi,double *taui, double *X1term, double *X2term, int *maxrow, int *maxcol)
   {
   double taurwxy, taurwxz, taurwyz, taurwxyz, term;	   	   
   taurwxy = gettaurw(mat_X,mat_Y,wi,taui,X1term,X2term,maxrow,maxcol);
   taurwxz = gettaurw(mat_X,mat_Z,wi,taui,X1term,X2term,maxrow,maxcol);
   taurwyz = gettaurw(mat_Y,mat_Z,wi,taui,X1term,X2term,maxrow,maxcol);
   term = sqrt((1-pow(taurwxz,2))*(1-pow(taurwyz,2)));
   if (term == 0.){taurwxyz = -5.;}
     else {
         taurwxyz = (taurwxy - taurwxz*taurwyz)/term;
     }
   return(taurwxyz);
   }

/* Function to compute statistical significance for partial rowwise Kendall's Kr statistic for square matrices */

  void partialstatsigKr(double *X, double *Y, double *Z, int *nrow, int *perm, double *pvpartialKrright, double *pvpartialKrleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, iter, *index, pointi, pointj, valid;
   double *mat_X, *mat_Y, *mat_Z, *wi, *taui, *X1term, *X2term, taurwxyz, *matperm, *setn, taurwxyzsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory the size of original matrix Z, exit if error */

   mat_Z = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Z == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Z - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_Z[i*maxcol+j] = Z[m];
         m++;
   }

/* Allocate in memory size of vector wi, ri, X1term and X2term */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   taui = calloc(maxrow, sizeof(double));

   if (taui == NULL)
   {
     error("Null dimension");
   }

   X1term = calloc(maxrow, sizeof(double));

   if (X1term == NULL)
   {
     error("Null dimension");
   }
   
   X2term = calloc(maxrow, sizeof(double));

   if (X2term == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   taurwxyz = gettaurwxyz(mat_X,mat_Y,mat_Z,wi,taui,X1term,X2term,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices under the specified null hypothesis */

   *pvpartialKrright = 0.;
   *pvpartialKrleft = 0.;
   valid = *perm;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow);
   for (i = 0; i< maxrow; i++){
      setn[i] = index[i];
   }
   for (i= 0; i< (maxrow-1); i++){
      pointi = setn[i];     
      for (j= i+1; j< maxcol; j++){
         pointj = setn[j];
         matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
         matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
      for (i = 0; i < maxrow; i++)
         for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}

   taurwxyzsim = gettaurwxyz(matperm,mat_Y,mat_Z,wi,taui,X1term,X2term,&maxrow,&maxcol);
   if (taurwxyzsim != -5.){
   if (taurwxyzsim >= taurwxyz)
        {*pvpartialKrright += 1.;}
      if (taurwxyzsim <= taurwxyz)
        {*pvpartialKrleft += 1.;}
   }
     else { valid -= 1;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(X2term);
   free(X1term);
   free(taui);
   free(wi);
   free(mat_Z);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvpartialKrright = (*pvpartialKrright+1.)/(valid+1.);
   *pvpartialKrleft =(*pvpartialKrleft+1.)/(valid+1.);
}

/* Function to compute statistical significance for partial rowwise Kendall's Kr statistic for rectangular matrices */

  void rectangularstatsigpartKr(double *X, double *Y, double *Z, int *nrow, int *ncol, int *perm, double *pvpartialKrright, double *pvpartialKrleft)
  {
/* Declaration of local variables */

   int i, j, k, m, maxrow, maxcol, iter, *index, point, valid;
   double *mat_X, *mat_Y, *mat_Z, *wi, *taui, *X1term, *X2term, taurwxyz, *matperm, *setn, taurwxyzsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = *nrow;
   maxcol = *ncol;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory the size of original matrix Z, exit if error */

   mat_Z = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Z == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Z - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_Z[i*maxcol+j] = Z[m];
         m++;
   }

/* Allocate in memory size of vector wi, ri, X1term and X2term */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   taui = calloc(maxrow, sizeof(double));

   if (taui == NULL)
   {
     error("Null dimension");
   }

   X1term = calloc(maxrow, sizeof(double));

   if (X1term == NULL)
   {
     error("Null dimension");
   }
   
   X2term = calloc(maxrow, sizeof(double));

   if (X2term == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   taurwxyz = gettaurwxyz(mat_X,mat_Y,mat_Z,wi,taui,X1term,X2term,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow*maxcol, sizeof(double));

   if (setn == NULL)
   {
     error("Null dimension");
   }

   index = calloc(maxrow*maxcol, sizeof(double));

   if (index == NULL)
   {
     error("Null dimension");
   }

/* Carry out the permutation procedure: set several p values. Then randomly permuted sociomatrices under the specified null hypothesis */

   *pvpartialKrright = 0.;
   *pvpartialKrleft = 0.;
   valid = *perm;
   for (iter = 0; iter < *perm; iter++) {
      for (i = 0; i< maxrow*maxcol; i++){
         setn[i] = runif(0,1);
	 index[i] = i;
      }
   rsort_with_index(setn,index,maxrow*maxcol);
   for (i = 0; i< maxrow*maxcol; i++){
      setn[i] = index[i];
   }
   k = 0;
   for (i= 0; i< maxrow; i++){     
      for (j= 0; j< maxcol; j++){
         point = setn[k];
         matperm[i*maxcol+j] = mat_X[point];
         k += 1;}}

   taurwxyzsim = gettaurwxyz(matperm,mat_Y,mat_Z,wi,taui,X1term,X2term,&maxrow,&maxcol);
   if (taurwxyzsim != -5.){
   if (taurwxyzsim >= taurwxyz)
        {*pvpartialKrright += 1.;}
      if (taurwxyzsim <= taurwxyz)
        {*pvpartialKrleft += 1.;}
   }
     else { valid -= 1;}
   }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(index);
   free(setn);
   free(matperm);
   free(X2term);
   free(X1term);
   free(taui);
   free(wi);
   free(mat_Z);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   *pvpartialKrright = (*pvpartialKrright+1.)/(valid+1.);
   *pvpartialKrleft =(*pvpartialKrleft+1.)/(valid+1.);
}

/* A required function to compute exact statistical significance for partial rowwise Kendall's Kr statistic */

   void permutetaurwxyz(int *setn,double *matperm,double *mat_X, double *mat_Y, double *mat_Z, double *wi, double *taui, double *X1term,
	       double *X2term, double *results, double taurwxyz, const int start, int maxrow, int maxcol)
   {
    int i, j, pointi, pointj;
    double taurwxyzsim;
    for (i= 0; i< (maxrow-1); i++){
       pointi = setn[i];     
       for (j= i+1; j< maxcol; j++){
          pointj = setn[j];
          matperm[i*maxcol+j] = mat_X[pointi*maxcol+pointj];
          matperm[j*maxcol+i] = mat_X[pointj*maxcol+pointi];}}
       for (i = 0; i < maxrow; i++)
          for (j = 0; j < maxrow; j++)
             {if (i==j){matperm[i*maxrow+j]=0.;}}
   
   taurwxyzsim = gettaurwxyz(matperm,mat_Y,mat_Z,wi,taui,X1term,X2term,&maxrow,&maxcol);
   if (taurwxyzsim != -5){
   if (taurwxyzsim >= taurwxyz)
   {results[0] += 1.;}
   if (taurwxyzsim <= taurwxyz)
   {results[1] += 1.;}
   }
   else {results[2] += 1.;}
   if (start < maxrow) {
      int i, j;
      for (i = maxrow-2; i >= start; i--) {
         for (j = i + 1; j < maxrow; j++) {
	    swap(setn, i, j);
	    permutetaurwxyz(setn,matperm,mat_X,mat_Y,mat_Z,wi,taui,X1term,X2term,results,taurwxyz,i+1,maxrow,maxcol);
         }
      rotateLeft(setn, i, maxrow);
    }
    }
   }

/* Function to compute exact statistical significance for partial rowwise Kendall's Kr statistic */

void partialexactsigKr(double *X, double *Y, double *Z, int *nrow, double *pvpartialKrright, double *pvpartialKrleft)
{
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, *setn, start, fact;
   double *mat_X, *mat_Y, *mat_Z, *wi, *taui, *X1term, *X2term, taurwxyz, *matperm, *results;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

   mat_X = malloc(maxrow*maxcol*sizeof(double));
   if (mat_X == NULL)
   {
     error("Null dimension");
   }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of original matrix Y, exit if error */

   mat_Y = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Y == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Y - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_Y[i*maxcol+j] = Y[m];
      m++;
   }

/* Allocate in memory the size of original matrix Z, exit if error */

   mat_Z = malloc(maxrow*maxcol*sizeof(double));
   if (mat_Z == NULL)
   {
     error("Null dimension");
   }

/* Transform vector Z - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_Z[i*maxcol+j] = Z[m];
         m++;
   }

/* Allocate in memory size of vector wi, ri, X1term and X2term */

   wi = calloc(maxrow, sizeof(double));

   if (wi == NULL)
   {
     error("Null dimension");
   }

   taui = calloc(maxrow, sizeof(double));

   if (taui == NULL)
   {
     error("Null dimension");
   }

   X1term = calloc(maxrow, sizeof(double));

   if (X1term == NULL)
   {
     error("Null dimension");
   }
   
   X2term = calloc(maxrow, sizeof(double));

   if (X2term == NULL)
   {
     error("Null dimension");
   }

/* Compute statistic value for the original sociomatrix */

   taurwxyz = gettaurwxyz(mat_X,mat_Y,mat_Z,wi,taui,X1term,X2term,&maxrow,&maxcol);

/* Allocate in memory the size of permuted sociomatrix, exit if error */   

   matperm = malloc(maxrow*maxcol*sizeof(double));
   if (matperm == NULL)
   {
     error("Null dimension");
   }
   
/* Allocate in memory size of vector setn, index and res */

   setn = calloc(maxrow, sizeof(int));

   if (setn == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory size of vector results and compute factorial number of n */

   fact = 1;
   for (i = maxrow;i >= 1; i--){
      fact *= i;}

   results = calloc(3,sizeof(double));

   if (results == NULL)
   {
     error("Null dimension");
   }

/* Carry out the exact permutation procedure */

   for (i = 0; i < maxrow; i++){
      setn[i] = i;}

   for (i = 0; i < 3; i++){
      results[i] = 0.;}

   start =0;
   permutetaurwxyz(setn,matperm,mat_X,mat_Y,mat_Z,wi,taui,X1term,X2term,results,taurwxyz,start,maxrow,maxcol);

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */
      
   free(results);
   free(setn);
   free(matperm);
   free(X2term);
   free(X1term);
   free(taui);
   free(wi);
   free(mat_Z);
   free(mat_Y);
   free(mat_X);

/* Compute p values for the statistic under the specified null hypothesis */

   fact -= results[2];
   *pvpartialKrright = (results[0])/(fact);
   *pvpartialKrleft =(results[1])/(fact);
}

/******************************************************************************************************************************/
/******************************************************************************************************************************/

/* A set of functions to generate random sociomatrices under the null hypothesis and to estimate statistical significance 
   for several measures of social reciprocity*/

/* Function to obtain matrix of dyadic directional asymmetry - phirmat - */
   	 
   void getphirmat(double *mat_X, double *dyadc, double *phirmat, int *maxrow, int *maxcol)
   {
   int i, j;
   double z, max, oddterm, min;
 
   z = 0.;
   for (i= 0; i< *maxrow; i++)
      for (j= i+1; j< *maxcol; j++)
      {
	 if (dyadc[i**maxcol+j] != 0.)
           z += 1 ;
      }
   max = z;
   oddterm =0.;
   for (i= 0; i< *maxrow; i++)
      for (j= i+1; j< *maxcol; j++)
      { 
         if (fmod(dyadc[i**maxcol+j],2.) != 0.)
         {oddterm = oddterm + (1/(pow(dyadc[i**maxcol+j],2)));}
	 {oddterm = oddterm;}
      }
   oddterm = oddterm/2.;
   min = (z/2.)+oddterm;
   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
      {
	 if (i == j){phirmat[i**maxcol+j]=0.;}
	 else
	 {
	   if  ( (fmod(dyadc[i**maxcol+j],2.) != 0.) & (dyadc[i**maxcol+j] != 0.) )
	   {phirmat[i**maxcol+j]=(((4*pow(mat_X[i**maxcol+j],2)/pow(dyadc[i**maxcol+j],2))-1)-
	    (1/pow(dyadc[i**maxcol+j],2)))/(4*max-4*min);}
	   if  ( (fmod(dyadc[i**maxcol+j],2.) == 0.) & (dyadc[i**maxcol+j] != 0.) )
	     {phirmat[i**maxcol+j]=((4*pow(mat_X[i**maxcol+j],2)/pow(dyadc[i**maxcol+j],2))-1)/(4*max-4*min);}
	   if (dyadc[i**maxcol+j] == 0.) {phirmat[i**maxcol+j]=0.;}
         }
      }
   }

/* Function to obtain the overall asymmetry index - PHIr - */

   double getPHIr(double *mat_X, double *dyadc, double *phirmat, int *maxrow, int *maxcol)
   {
   int i, j;
   double PHIr;
   getphirmat(mat_X,dyadc,phirmat,maxrow,maxcol);		   
   PHIr = 0.;
   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
         if (i != j)
           {PHIr += phirmat[i**maxcol+j];}
   return(PHIr);
   }

/* Function to obtain vector of individual contribution to asymetry as actor - phii - */

   void getphii (double *mat_X,double *dyadc,double *phirmat, double *phii, int *maxrow, int *maxcol)
   {
   int i, j;

   getphirmat(mat_X,dyadc,phirmat,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      phii[i] = 0.;
   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)      
         phii[i] += phirmat[i**maxcol+j];
   }
   
/* Function to obtain vector of individual contribution to asymetry as partner - phij - */

   void getphij (double *mat_X,double *dyadc,double *phirmat, double *phij, int *maxrow, int *maxcol)
   {
   int i, j;

   getphirmat(mat_X,dyadc,phirmat,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      phij[i] = 0.;

   for (j= 0; j< *maxcol; j++)
      for (i= 0; i< *maxrow; i++)
         phij[i] += phirmat[j**maxcol+i];
   
   }

/* Function to obtain matrix of dyadic contribution to asymmetry - PHIij - */
   	 
   void getPHIij(double *mat_X, double *dyadc, double *phirmat,double *PHIij, int *maxrow, int *maxcol)
   {
   int i, j;

   getphirmat(mat_X,dyadc,phirmat,maxrow,maxcol);

   for (i= 0; i< *maxrow; i++)
      for (j= 0; j< *maxcol; j++)
	 PHIij[i**maxcol+j]=phirmat[i**maxcol+j]+phirmat[j**maxcol+i];

   }

/* Function to estimate statistical significance for social reciprocity measurements */

  void recip(double *X, double *pi, int *nrow, int *rep, double *pvPHIrright, double *pvPHIrleft, double *pvphirmatright,
  	     double *pvphirmatleft, double *pvphiiright, double *pvphiileft, double *pvphijright, double *pvphijleft,
	     double *pvPHIijright, double *pvPHIijleft)
  {
/* Declaration of local variables */

   int i, j, m, maxrow, maxcol, iter;
   double *mat_X, *mat_pi, *dyadc, *phirmat, *phii, *phij, *PHIij, PHIr1, *phirmat1, *phii1, *phij1, *PHIij1,
	  *matgen, PHIrsim, *phirmatsim, *phiisim, *phijsim, *PHIijsim;

/* Setting the random seed for the random number generation routine */

   GetRNGstate();

/* Setting size of matrices */

   maxrow = maxcol = *nrow;

/* Allocate in memory the size of original matrix X, exit if error */

  mat_X = malloc(maxrow*maxcol*sizeof(double));
  if (mat_X == NULL)
  {
  error("Null dimension");
  }

/* Transform vector X - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
   {   
         mat_X[i*maxcol+j] = X[m];
         m++;
   }

/* Allocate in memory the size of matrix pi, exit if error */

   mat_pi = malloc(maxrow*maxcol*sizeof(double));
   if (mat_pi == NULL)
   {
     error("Null dimension");
   }

/* Transform vector pi - from R - into a matrix */

   m = 0.;
   for (i= 0; i< maxrow; i++)
   for (j= 0; j< maxcol; j++)
   {   
      mat_pi[i*maxcol+j] = pi[m];
      m++;
   }

/* Allocate in memory the size of matrix with dyadic interaction frequencies dyadc, exit if error. Then define and obtain its elements */

   dyadc = malloc(maxrow*maxcol*sizeof(double));
   if (dyadc == NULL)
   {
   error("Null dimension");
   }

   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
         dyadc[i*maxcol+j] = mat_X[i*maxcol+j] + mat_X[j*maxcol+i];
   
/* Allocate in memory the size of matrix of dyadic directional asymmetry, exit if error */

   phirmat = malloc(maxrow*maxcol*sizeof(double));
   if (phirmat == NULL)
   {
     error("Null dimension");
   }

/* Allocate size of vector of individual contribution to asymmetry as actor, exit if error */

   phii = calloc(maxrow, sizeof(double));

   if (phii == NULL)
   {
     error("Null dimension");
   }

/* Allocate size of vector of individual contribution to asymmetry as partner, exit if error */

   phij = calloc(maxrow, sizeof(double));

   if (phij == NULL)
   {
     error("Null dimension");
   }

/* Allocate in memory the size of matrix of dyadic contribution to asymmetry, exit if error */

   PHIij = malloc(maxrow*maxcol*sizeof(double));
   if (PHIij == NULL)
   {
     error("Null dimension");
   }

/* Compute social reciprocity measures at different levels for the original sociomatrix */

   PHIr1 = getPHIr(mat_X,dyadc,phirmat,&maxrow,&maxcol);

   phirmat1 = malloc(maxrow*maxcol*sizeof(double));
   if (phirmat1 == NULL)
   {
     error("Null dimension");
   }
   getphirmat(mat_X,dyadc,phirmat1,&maxrow,&maxcol);

   phii1 = calloc(maxrow, sizeof(double));

   if (phii1 == NULL)
   {
     error("Null dimension");
   }
   getphii(mat_X,dyadc,phirmat1,phii1,&maxrow,&maxcol);

   phij1 = calloc(maxrow, sizeof(double));

   if (phij1 == NULL)
   {
     error("Null dimension");
   }
   getphij(mat_X,dyadc,phirmat1,phij1,&maxrow,&maxcol);

   PHIij1 = malloc(maxrow*maxcol*sizeof(double));
   if (PHIij1 == NULL)
   {
     error("Null dimension");
   }
   getPHIij(mat_X,dyadc,phirmat1,PHIij1,&maxrow,&maxcol);

/* Allocate in memory the size of random sociomatrix, exit if error */   

   matgen = malloc(maxrow*maxcol*sizeof(double));
   if (matgen == NULL)
   {
     error("Null dimension");
   }

/* Carry out the MC sampling: set several p values and allocate simulated matrices. 
   Then generate random sociomatrices under the specified null hypothesis */

   *pvPHIrright = 0.;
   *pvPHIrleft = 0.;
   
   phirmatsim = malloc(maxrow*maxcol*sizeof(double));
   if (phirmatsim == NULL)
   {
     error("Null dimension");
   }

   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {
         pvphirmatright[i*maxcol+j] = 0.;
         pvphirmatleft[i*maxcol+j] = 0.;
      }

   phiisim = calloc(maxrow, sizeof(double));

   if (phiisim == NULL)
   {
     error("Null dimension");
   }
   for (i= 0; i< maxrow; i++)
   {
      pvphiiright[i] = 0.;
      pvphiileft[i] = 0.;
   }

   phijsim = calloc(maxrow, sizeof(double));

   if (phijsim == NULL)
   {
     error("Null dimension");
   }
   for (i= 0; i< maxrow; i++)
   {
      pvphijright[i] = 0.;
      pvphijleft[i] = 0.;
   }

   PHIijsim = malloc(maxrow*maxcol*sizeof(double));
   if (PHIijsim == NULL)
   {
     error("Null dimension");
   }

   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {
         pvPHIijright[i*maxcol+j] = 0.;
         pvPHIijleft[i*maxcol+j] = 0.;
      }

   for (iter= 0; iter< *rep; iter++) 
   {
      for (i= 0; i< maxrow; i++)
         for (j= 0; j< maxcol; j++)
         {    
            if (i < j){matgen[i*maxcol+j]=rbinom(dyadc[i*maxcol+j],mat_pi[i*maxcol+j]);}
            if (i > j){matgen[i*maxcol+j]=dyadc[i*maxcol+j]-matgen[j*maxcol+i];}
            if (i == j){matgen[i*maxcol+j]=0.;}
         }
 
   PHIrsim = getPHIr(matgen,dyadc,phirmat,&maxrow,&maxcol);
   if (PHIrsim >= PHIr1)
     {*pvPHIrright += 1.;}
   if (PHIrsim <= PHIr1)
     {*pvPHIrleft += 1.;}

   getphirmat(matgen,dyadc,phirmatsim,&maxrow,&maxcol);

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {
	 if (i != j)
         { 
           if (phirmatsim[i*maxcol+j] >= phirmat1[i*maxcol+j])
             {pvphirmatright[m] += 1.;}
           if (phirmatsim[i*maxcol+j] <= phirmat1[i*maxcol+j])
             {pvphirmatleft[m] += 1.;}
	 }
	 m += 1;
      }

   getphii(matgen,dyadc,phirmatsim,phiisim,&maxrow,&maxcol);

   for (i= 0; i< maxrow; i++)
   {
      if (phiisim[i] >= phii1[i]) 
        {pvphiiright[i] += 1.;}
      if (phiisim[i] <= phii1[i])
        {pvphiileft[i] += 1.;}
   }

   getphij(matgen,dyadc,phirmatsim,phijsim,&maxrow,&maxcol);

   for (i= 0; i< maxrow; i++)
   {
      if (phijsim[i] >= phij1[i]) 
        {pvphijright[i] += 1.;}
      if (phijsim[i] <= phij1[i])
        {pvphijleft[i] += 1.;}
   }

   getPHIij(matgen,dyadc,phirmatsim,PHIijsim,&maxrow,&maxcol);

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {
	 if (i != j)
         { 
           if (PHIijsim[i*maxcol+j] >= PHIij1[i*maxcol+j])
             {pvPHIijright[m] += 1.;}
           if (PHIijsim[i*maxcol+j] <= PHIij1[i*maxcol+j])
             {pvPHIijleft[m] += 1.;}
	 }
	 m += 1;
      }

   }


/* Compute p values for the social reciprocity measures under the specified null hypothesis */

   *pvPHIrright = (*pvPHIrright+1.)/(*rep+1.);
   *pvPHIrleft =(*pvPHIrleft+1.)/(*rep+1.);

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {
	 if (i != j)
	 {
	   pvphirmatright[m] = (pvphirmatright[m]+1)/(*rep+1);
           pvphirmatleft[m] = (pvphirmatleft[m]+1)/(*rep+1);
         }
         else
	 {
            pvphirmatright[m] = 0.;
            pvphirmatleft[m] = 0.;
	 }
	 m += 1;
      }

   
   for (i= 0; i< maxrow; i++)
   { 
      pvphiiright[i] = (pvphiiright[i]+1)/(*rep+1);
      pvphiileft[i] = (pvphiileft[i]+1)/(*rep+1);
   }

   for (i= 0; i< maxrow; i++)
   {  
      pvphijright[i] = (pvphijright[i]+1)/(*rep+1);
      pvphijleft[i] = (pvphijleft[i]+1)/(*rep+1);
   }

   m = 0.;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {
	 if (i != j)
	 {
	   pvPHIijright[m] = (pvPHIijright[m]+1)/(*rep+1);
           pvPHIijleft[m] = (pvPHIijleft[m]+1)/(*rep+1);
         }
         else
	 {
            pvPHIijright[m] = 0.;
            pvPHIijleft[m] = 0.;
	 }
	 m += 1;
      }

/* Eliminate the random seed used in the routine */
	 	 
   PutRNGstate();

/* Deallocate all the matrices used in the routine */   

   free(PHIijsim);
   free(PHIij1);
   free(phijsim);
   free(phij1);
   free(phiisim);   
   free(phii1);
   free(phirmatsim);
   free(phirmat1);
   free(matgen);
   free(PHIij);
   free(phij);
   free(phii);
   free(phirmat);
   free(dyadc);
   free(mat_pi);
   free(mat_X);
   }


/******************************************************************************************************************************/
/******************************************************************************************************************************/
   
   /* Function to compute overall symmetry index - psi - */
   
   double getpsi(double *mat_X, double *tras_X, double *mat_S, double *tras_S, double *mat_XTX,
                 double *mat_STS, int *maxrow, int *maxcol)
   {
     int i, j, k;
     double trXTX, trSTS, psi;
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         tras_X[i**maxcol+j] = mat_X[j**maxcol+i];
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         mat_S[i**maxcol+j] = (mat_X[i**maxcol+j] + tras_X[i**maxcol+j])/2.;
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         tras_S[i**maxcol+j] = mat_S[j**maxcol+i];
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         mat_XTX[i**maxcol+j] = 0.;
     
     for (k= 0; k< *maxrow; k++)
       for (j= 0; j< *maxcol; j++)
         for (i= 0; i< *maxrow; i++)
           mat_XTX[k**maxcol+j] += (tras_X[k**maxcol+i] * mat_X[i**maxcol+j]);
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         mat_STS[i**maxcol+j] = 0.;
     
     for (k= 0; k< *maxrow; k++)
       for (j= 0; j< *maxcol; j++)
         for (i= 0; i< *maxrow; i++)
           mat_STS[k**maxcol+j] += (tras_S[k**maxcol+i] * mat_S[i**maxcol+j]);
     
     trXTX = 0.;
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         if (i == j)
           trXTX += mat_XTX[i**maxcol+j];
         trSTS = 0.;
         for (i= 0; i< *maxrow; i++)
           for (j= 0; j< *maxcol; j++)
             if (i == j)
               trSTS += mat_STS[i**maxcol+j];
             
             psi = trSTS / trXTX;
             
             return(psi);
   }
   
   /* Function to compute overall skew-symmetry index - phi - */
   
   double getphi(double *mat_X, double *tras_X, double *mat_K, double *tras_K, double *mat_XTX,
                 double *mat_KTK, int *maxrow, int *maxcol)
   {
     int i, j, k;
     double trXTX, trKTK, phi;
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         tras_X[i**maxcol+j] = mat_X[j**maxcol+i];
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         mat_K[i**maxcol+j] = (mat_X[i**maxcol+j] - tras_X[i**maxcol+j])/2.;
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         tras_K[i**maxcol+j] = mat_K[j**maxcol+i];
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         mat_XTX[i**maxcol+j] = 0.;
     
     for (k= 0; k< *maxrow; k++)
       for (j= 0; j< *maxcol; j++)
         for (i= 0; i< *maxrow; i++)
           mat_XTX[k**maxcol+j] += (tras_X[k**maxcol+i] * mat_X[i**maxcol+j]);
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         mat_KTK[i**maxcol+j] = 0.;
     
     for (k= 0; k< *maxrow; k++)
       for (j= 0; j< *maxcol; j++)
         for (i= 0; i< *maxrow; i++)
           mat_KTK[k**maxcol+j] += (tras_K[k**maxcol+i] * mat_K[i**maxcol+j]);
     
     trXTX = 0.;
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
       {   
         if (i == j)
           trXTX += mat_XTX[i**maxcol+j];
       }
       trKTK = 0.;
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
       {
         if (i == j)
           trKTK += mat_KTK[i**maxcol+j];
       }
       
       phi = trKTK / trXTX;
     
     return(phi);
   }
   
   /* Function to compute overall skew-symmetry/symmetry ratio - delta - */
   
   double getdelta(double *mat_X, double *tras_X, double *mat_S, double *tras_S, double *mat_XTX, double *mat_STS,
                   double *mat_K, double *tras_K, double *mat_KTK, int *maxrow, int *maxcol)
   {
     
     double psi, phi, delta;
     
     psi = getpsi(mat_X,tras_X,mat_S,tras_S,mat_XTX,mat_STS,maxrow,maxcol);
     
     phi = getphi(mat_X,tras_X,mat_K,tras_K,mat_XTX,mat_KTK,maxrow,maxcol);
     
     delta = phi /psi;
     
     return(delta);
   }
   
   /* Function to compute the directional consistency index -dc- */
   
   double getdc(double *mat_X, double *tras_X, int *maxrow, int *maxcol)
   {   
     
     int i, j;
     double totalN, sumdif, dc;
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         tras_X[i**maxcol+j] = mat_X[j**maxcol+i];
     
     totalN = 0;
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         totalN += mat_X[i**maxcol+j];
     
     sumdif = 0.;
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
       {
         if (j > i)	   
           sumdif += fabs(mat_X[i**maxcol+j]-tras_X[i**maxcol+j]);
       }
       
       dc = sumdif/totalN;
     return(dc);
   }
   
   /* Function to compute unweighted dyadic contributions to symmetry -lambda- */
   
   double getlambda(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *lambda,
                    int *maxrow, int *maxcol)
   {
     int i, j;  
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         tras_X[i**maxcol+j] = mat_X[j**maxcol+i];
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         mat_S[i**maxcol+j] = (mat_X[i**maxcol+j] + tras_X[i**maxcol+j])/2.;
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         mat_K[i**maxcol+j] = (mat_X[i**maxcol+j] - tras_X[i**maxcol+j])/2.;
     
     for (i= 0; i< *maxrow; i++)
       sumvect[i] = 0.;
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         sumvect[i] += ((mat_S[i**maxcol+j]*mat_S[i**maxcol+j]) + (mat_K[i**maxcol+j]*mat_K[i**maxcol+j]));
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
       {
         if (sumvect[i] != 0.)	      
           lambda[i**maxcol+j] = (mat_S[i**maxcol+j] * mat_S[i**maxcol+j])/sumvect[i];
         else
           lambda[i**maxcol+j] = 0.;
       }
       return(*lambda);
   }
   
   /* Function to obtain individuals' contribution to symmetry -lambdaj- */
   
   double getlambdaj(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *lambda,
                     double *lambdaj, int *maxrow, int *maxcol)
   {
     
     int i, j;
     
     *lambda = getlambda(mat_X,tras_X,mat_S,mat_K,sumvect,lambda,maxrow,maxcol);
     
     for (i= 0; i< *maxrow; i++)
       lambdaj[i] = 0.;
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         lambdaj[i] += lambda[i**maxcol+j];
     return(*lambdaj);
   }
   
   /* Function to obtain dyadic contributions to symmetry -ratiolambda- */
   
   double getratiolambda(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *lambda,
                         double *lambdaj, double *ratiolambda, int *maxrow, int *maxcol)
   {
     int i, j;
     
     *lambda = getlambda(mat_X,tras_X,mat_S,mat_K,sumvect,lambda,maxrow,maxcol);
     
     *lambdaj = getlambdaj(mat_X,tras_X,mat_S,mat_K,sumvect,lambda,lambdaj,maxrow,maxcol);
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
       {
         if (lambdaj[i] != 0.)	      
           ratiolambda[j**maxcol+i] = lambda[i**maxcol+j]/lambdaj[i];
         else
           ratiolambda[i**maxcol+j] = 0.;
       }
       return(*ratiolambda);
   }
   
   /* Function to obtain unweighted contributions to skew-symmetry -nu- */
   
   double getnu(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *nu,
                int *maxrow, int *maxcol)
   {
     int i, j;  
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         tras_X[i**maxcol+j] = mat_X[j**maxcol+i];
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         mat_S[i**maxcol+j] = (mat_X[i**maxcol+j] + tras_X[i**maxcol+j])/2.;
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         mat_K[i**maxcol+j] = (mat_X[i**maxcol+j] - tras_X[i**maxcol+j])/2.;
     
     for (i= 0; i< *maxrow; i++)
       sumvect[i] = 0.;
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         sumvect[i] += ((mat_S[i**maxcol+j]*mat_S[i**maxcol+j]) + (mat_K[i**maxcol+j]*mat_K[i**maxcol+j]));
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
       {
         if (sumvect[i] != 0.)	      
           nu[i**maxcol+j] = (mat_K[i**maxcol+j] * mat_K[i**maxcol+j])/sumvect[i];
         else
           nu[i**maxcol+j] = 0.;
       }
       return(*nu);
   }
   
   /* Function to obtain individuals' contribution to skew-symmetry -nuj- */
   
   double getnuj(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *nu,
                 double *nuj, int *maxrow, int *maxcol)
   {
     
     int i, j;
     
     *nu = getnu(mat_X,tras_X,mat_S,mat_K,sumvect,nu,maxrow,maxcol);
     
     for (i= 0; i< *maxrow; i++)
       nuj[i] = 0.;
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         nuj[i] += nu[i**maxcol+j];
     return(*nuj);
   }
   
   /* Function to obtain dyadic contributions to symmetry -rationu- */
   
   double getrationu(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *nu,
                     double *nuj, double *rationu, int *maxrow, int *maxcol)
   {
     int i, j;
     
     *nu = getnu(mat_X,tras_X,mat_S,mat_K,sumvect,nu,maxrow,maxcol);
     
     *nuj = getnuj(mat_X,tras_X,mat_S,mat_K,sumvect,nu,nuj,maxrow,maxcol);
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
       {
         if (nuj[i] != 0.)	      
           rationu[j**maxcol+i] = nu[i**maxcol+j]/nuj[i];
         else
           rationu[i**maxcol+j] = 0.;
       }
       return(*rationu);
   }
   
   /* Function to obtain matrix of dyadic balanced reciprocity -omega- */
   
   double getomega(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *lambda,
                   double *nu, double *omega, int *maxrow, int *maxcol)
   {
     int i, j;
     
     *lambda = getlambda(mat_X,tras_X,mat_S,mat_K,sumvect,lambda,maxrow,maxcol);
     
     *nu = getnu(mat_X,tras_X,mat_S,mat_K,sumvect,nu,maxrow,maxcol);
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
       {
         if (lambda[i**maxcol+j] != 0.)	      
           omega[i**maxcol+j] = nu[i**maxcol+j]/lambda[i**maxcol+j];
         else
           omega[i**maxcol+j] = 0.;
       }
       return(*omega);
   }   
   
   /* Function to obtain dyadic reciprocity index -kappa- */
   
   double getkappa(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *nu,
                   double *nuj, double *rationu, double *difrationu, int *maxrow, int *maxcol)
   {
     int i, j;
     double kappa;
     
     *rationu = getrationu(mat_X,tras_X,mat_S,mat_K,sumvect,nu,nuj,rationu,maxrow,maxcol);
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         difrationu[i**maxcol+j] = fabs(rationu[i**maxcol+j]-rationu[j**maxcol+i]);
     
     kappa = 0.;
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
       {	      
         if (j > i)
           kappa += difrationu[i**maxcol+j];
       }
       
       kappa = 1 -(kappa/(*maxrow));
     
     return(kappa);
   }
   
   /* Compute generalized reciprocity based on skewsymmetrical part */
   
   double getepsilon(double *mat_X, double *tras_X, double *mat_S, double *mat_K, double *sumvect, double *nu,
                     double *nuj, double *nui, int *maxrow, int *maxcol)
   {
     
     int i, j;
     double sumnuj, epsilon;
     
     *nu = getnu(mat_X,tras_X,mat_S,mat_K,sumvect,nu,maxrow,maxcol);
     
     *nuj = getnuj(mat_X,tras_X,mat_S,mat_K,sumvect,nu,nuj,maxrow,maxcol);
     
     for (i= 0; i< *maxrow; i++)
       nui[i] = 0.;
     
     for (i= 0; i< *maxrow; i++)
       for (j= 0; j< *maxcol; j++)
         nui[i] += nu[j**maxcol+i]; 
     
     epsilon = 0.;
     for (i= 0; i< *maxrow; i++)
       epsilon += fabs(nuj[i]-nui[i]);
     
     sumnuj = 0.;
     for (i= 0.; i< *maxrow; i++)
       sumnuj += nuj[i];
     
     if (sumnuj != 0.)
       epsilon = 1 - (epsilon/sumnuj);
     else
       epsilon = 1.;
     
     return(epsilon);
   }
   
   
/* Program to generate random sociomatrices under the null hypothesis and to estimate statistical significance for several measures of social reciprocity */


void recip0(double *X, double *pi, int *nrow, int *rep, double *pvphi, double *pvpsi, double *pvdc, double *pvdelta, double *pvepsilon,
           double *pvkappa, double *pvnuj, double *pvlambdaj, double *pvrationu, double *pvratiolambda, double *pvomega)
{
  
  /* Declaration of local variables */
  
  int i, j, m, maxrow, maxcol, iter;
  double *mat_X, *mat_pi, *dyadc, *tras_X, *mat_S, *mat_K, *tras_S, *tras_K, *mat_STS, *mat_KTK, *mat_XTX, phi, psi, delta, dc,
  epsilon,kappa, *nu, *lambda, *nuj, *nui, *lambdaj, *rationu, *ratiolambda, *difrationu, *sumvect, *omega, *matgen,
  phisim, psisim, deltasim, dcsim, epsilonsim, kappasim, *nusim, *lambdasim, *nujsim, *lambdajsim, *rationusim,
  *ratiolambdasim, *omegasim;
  
  /* Setting the random seed for the random number generation routine */
  
  GetRNGstate();
  
  /* Setting size of matrices */
  
  maxrow = maxcol = *nrow;
  
  /* Allocate in memory the size of original matrix X, exit if error */
  
  mat_X = malloc(maxrow*maxcol*sizeof(double));
  if (mat_X == NULL)
  {
    error("Null dimension");
  }
  
  /* Transform vector X - from R - into a matrix */
  
  m = 0.;
  for (i= 0; i< maxrow; i++)
    for (j= 0; j< maxcol; j++)
    {   
      mat_X[i*maxcol+j] = X[m];
      m++;
    }
    
    /* Allocate in memory the size of matrix pi, exit if error */
    
    mat_pi = malloc(maxrow*maxcol*sizeof(double));
  if (mat_pi == NULL)
  {
    error("Null dimension");
  }
  
  /* Transform vector pi - from R - into a matrix */
  
  m = 0.;
  for (i= 0; i< maxrow; i++)
    for (j= 0; j< maxcol; j++)
    {   
      mat_pi[i*maxcol+j] = pi[m];
      m++;
    }
    
    /* Allocate in memory the size of matrix with dyadic interaction frequencies dyadc, exit if error. Then define and obtain its elements */
    
    dyadc = malloc(maxrow*maxcol*sizeof(double));
  if (dyadc == NULL)
  {
    error("Null dimension");
  }
  
  for (i= 0; i< maxrow; i++)
    for (j= 0; j< maxcol; j++)
      dyadc[i*maxcol+j] = mat_X[i*maxcol+j] + mat_X[j*maxcol+i];
  
  /* Allocate in memory the transposed of original matrix */
  
  tras_X = malloc(maxrow*maxcol*sizeof(double));
  
  if (tras_X == NULL)
  {
    error("Null dimension");
  } 
  
  /* Allocate in memory the symmetrical part of matrix X */
  
  mat_S = malloc(maxrow*maxcol*sizeof(double));
  
  if (mat_S == NULL)
  {
    error("Null dimension");
  } 
  
  /* Allocate in memory the skewsymmetrical part of matrix X */
  
  mat_K = malloc(maxrow*maxcol*sizeof(double));
  
  if (mat_K == NULL)
  {
    error("Null dimension");
  }  
  
  /* Allocate in memory the transposed of matrix S */
  
  tras_S = malloc(maxrow*maxcol*sizeof(double));
  
  if (tras_S == NULL)
  {
    error("Null dimension");
  } 
  
  /* Allocate in memory the transposed of matrix K */
  
  tras_K = malloc(maxrow*maxcol*sizeof(double));
  
  if (tras_K == NULL)
  {
    error("Null dimension");
  } 
  
  /* Allocate in memory the matrix X'X */
  
  mat_XTX = malloc(maxrow*maxcol*sizeof(double));
  
  if (mat_XTX == NULL)
  {
    error("Null dimension");
  } 
  
  /* Allocate in memory the matrix S'S */
  
  mat_STS = malloc(maxrow*maxcol*sizeof(double));
  
  if (mat_STS == NULL)
  {
    error("Null dimension");
  } 
  
  /* Allocate in memory the matrix K'K */
  
  mat_KTK = malloc(maxrow*maxcol*sizeof(double));
  
  if (mat_KTK == NULL)
  {
    error("Null dimension");
  } 
  
  /* Allocate matrix of dyadic decomposition of symmetrical part -lambda- */
  
  sumvect = calloc(maxrow, sizeof(double));
  
  if (sumvect == NULL)
  {
    error("Null dimension");
  }
  
  lambda = malloc(maxrow*maxcol*sizeof(double));
  
  if (lambda == NULL)
  {
    error("Null dimension");
  }
  
  /* Allocate array of individuals contribution to symmetry -lambdaj- */
  
  lambdaj = calloc(maxrow, sizeof(double));
  
  if (lambdaj == NULL)
  {
    error("Null dimension");
  }
  
  /* Allocate matrix of dyadic contribution to symmetry -ratiolambda- */
  
  ratiolambda = malloc(maxrow*maxcol*sizeof(double));
  
  if (ratiolambda == NULL)
  {
    error("Null dimension");
  }
  
  /* Allocate matrix of unweighted dyadic decomposition of skew-symmetrical part -nu- */
  
  nu = malloc(maxrow*maxcol*sizeof(double));
  
  if (nu == NULL)
  {
    error("Null dimension");
  }
  
  /* Allocate array of individuals contribution to skew-symmetry -nuj- */   
  
  nuj = calloc(maxrow, sizeof(double));
  
  if (nuj == NULL)
  {
    error("Null dimension");
  }
  
  nui = calloc(maxrow, sizeof(double));
  
  if (nui == NULL)
  {
    error("Null dimension");
  }
  
  /* Allocate matrix of dyadic contribution to skew-symmetry -rationu- */
  
  rationu = malloc(maxrow*maxcol*sizeof(double));
  
  if (rationu == NULL)
  {
    error("Null dimension");
  } 
  
  /* Allocate matrix omega of ratios nu/lambda */
  
  omega = malloc(maxrow*maxcol*sizeof(double));
  
  if (omega == NULL)
  {
    error("Null dimension");
  }
  
  /* Allocate matrix necessary for computating kappa index */
  
  difrationu = malloc(maxrow*maxcol*sizeof(double));
  
  if (difrationu == NULL)
  {
    error("Null dimension");
  }
  

  /* Compute social reciprocity measures at different levels for the original sociomatrix */
  
  psi = getpsi(mat_X,tras_X,mat_S,tras_S,mat_XTX,mat_STS,&maxrow,&maxcol);
  
  phi = getphi(mat_X,tras_X,mat_K,tras_K,mat_XTX,mat_KTK,&maxrow,&maxcol);
  
  delta = getdelta(mat_X,tras_X,mat_S,tras_S,mat_XTX,mat_STS,mat_K,tras_K,mat_KTK,&maxrow,&maxcol);
  
  dc = getdc(mat_X,tras_X,&maxrow,&maxcol);
  
  *lambda = getlambda(mat_X,tras_X,mat_S,mat_K,sumvect,lambda,&maxrow,&maxcol);
  
  *lambdaj = getlambdaj(mat_X,tras_X,mat_S,mat_K,sumvect,lambda,lambdaj,&maxrow,&maxcol);
  
  *ratiolambda = getratiolambda(mat_X,tras_X,mat_S,mat_K,sumvect,lambda,lambdaj,ratiolambda,&maxrow,&maxcol);
  
  *nu = getnu(mat_X,tras_X,mat_S,mat_K,sumvect,nu,&maxrow,&maxcol);
  
  *nuj = getnuj(mat_X,tras_X,mat_S,mat_K,sumvect,nu,nuj,&maxrow,&maxcol);
  
  *rationu = getrationu(mat_X,tras_X,mat_S,mat_K,sumvect,nu,nuj,rationu,&maxrow,&maxcol);
  
  *omega = getomega(mat_X,tras_X,mat_S,mat_K,sumvect,lambda,nu,omega,&maxrow,&maxcol);
  
  kappa = getkappa(mat_X,tras_X,mat_S,mat_K,sumvect,nu,nuj,rationu,difrationu,&maxrow,&maxcol);
  
  epsilon = getepsilon(mat_X,tras_X,mat_S,mat_K,sumvect,nu,nuj,nui,&maxrow,&maxcol);
  
  /* Allocate in memory the size of random sociomatrix, exit if error */   
  
  matgen = malloc(maxrow*maxcol*sizeof(double));
  if (matgen == NULL)
  {
    error("Null dimension");
  }
  
  /* Carry out the MC sampling: set several p values and allocate simulated matrices. Then generate random sociomatrices under the specified null hypothesis */
  
  *pvphi = 0.;
  *pvpsi = 0.;
  *pvdelta = 0.;
  *pvdc = 0.;
  *pvkappa = 0.;
  *pvepsilon = 0.;
  
  lambdasim = malloc(maxrow*maxcol*sizeof(double));
  if (lambdasim == NULL)
  {
    error("Null dimension");
  }
  
  lambdajsim = calloc(maxrow, sizeof(double));
  
  if (lambdajsim == NULL)
  {
    error("Null dimension");
  }
  
  for (i= 0; i< maxrow; i++)
  {
    pvlambdaj[i] = 0.;
  }
  
  ratiolambdasim = malloc(maxrow*maxcol*sizeof(double));
  if (ratiolambdasim == NULL)
  {
    error("Null dimension");
  }
  
  for (i= 0; i< maxrow; i++)
    for (j= 0; j< maxcol; j++)
    {
      pvratiolambda[i*maxcol+j] = 0.;
    }
    
    
    nusim = malloc(maxrow*maxcol*sizeof(double));
  if (nusim == NULL)
  {
    error("Null dimension");
  }
  
  nujsim = calloc(maxrow, sizeof(double));
  
  if (nujsim == NULL)
  {
    error("Null dimension");
  }
  
  for (i= 0; i< maxrow; i++)
  {
    pvnuj[i] = 0.;
  }
  
  rationusim = malloc(maxrow*maxcol*sizeof(double));
  if (rationusim == NULL)
  {
    error("Null dimension");
  }
  
  for (i= 0; i< maxrow; i++)
    for (j= 0; j< maxcol; j++)
    {
      pvrationu[i*maxcol+j] = 0.;
    }
    
    omegasim = malloc(maxrow*maxcol*sizeof(double));
  if (omegasim == NULL)
  {
    error("Null dimension");
  }
  
  for (i= 0; i< maxrow; i++)
    for (j= 0; j< maxcol; j++)
    {
      pvomega[i*maxcol+j] = 0.;
    }
    
    for (iter= 0; iter< *rep; iter++) 
    {
      for (i= 0; i< maxrow; i++)
        for (j= 0; j< maxcol; j++)
        {    
          if (i < j){matgen[i*maxcol+j]=rbinom(dyadc[i*maxcol+j],mat_pi[i*maxcol+j]);}
          if (i > j){matgen[i*maxcol+j]=dyadc[i*maxcol+j]-matgen[j*maxcol+i];}
          if (i == j){matgen[i*maxcol+j]=0.;}
        }
        
        psisim = getpsi(matgen,tras_X,mat_S,tras_S,mat_XTX,mat_STS,&maxrow,&maxcol);
      
      if (psisim <= psi)
      {*pvpsi += 1.;}
      
      phisim = getphi(matgen,tras_X,mat_K,tras_K,mat_XTX,mat_KTK,&maxrow,&maxcol);
      
      if (phisim >= phi)
      {*pvphi += 1.;}
      
      deltasim = getdelta(matgen,tras_X,mat_S,tras_S,mat_XTX,mat_STS,mat_K,tras_K,mat_KTK,&maxrow,&maxcol);
      
      if (deltasim >= delta)
      {*pvdelta += 1.;}
      
      dcsim = getdc(matgen,tras_X,&maxrow,&maxcol);
      
      if (dcsim >= dc)
      {*pvdc += 1.;}
      
      *lambdasim = getlambda(matgen,tras_X,mat_S,mat_K,sumvect,lambdasim,&maxrow,&maxcol);
      
      *lambdajsim = getlambdaj(matgen,tras_X,mat_S,mat_K,sumvect,lambdasim,lambdajsim,&maxrow,&maxcol);
      
      for (i= 0; i< maxrow; i++)
      {
        if (lambdajsim[i] <= lambdaj[i])
        {pvlambdaj[i] += 1.;}
      }
      
      *nusim = getnu(matgen,tras_X,mat_S,mat_K,sumvect,nusim,&maxrow,&maxcol);
      
      *nujsim = getnuj(matgen,tras_X,mat_S,mat_K,sumvect,nusim,nujsim,&maxrow,&maxcol);
      
      for (i= 0; i< maxrow; i++)
      {
        if (nujsim[i] >= nuj[i])
        {pvnuj[i] += 1.;}
      }
      
      *ratiolambdasim = getratiolambda(matgen,tras_X,mat_S,mat_K,sumvect,lambdasim,lambdajsim,ratiolambdasim,&maxrow,&maxcol);
      
      m = 0.;
      for (i= 0; i< maxrow; i++)
        for (j= 0; j< maxcol; j++)
        {
          if (i != j)
          { 
            if (ratiolambdasim[i*maxcol+j] <= ratiolambda[i*maxcol+j])
            {pvratiolambda[m] += 1.;}
          }
          m += 1;
        }
        
        *rationusim = getrationu(matgen,tras_X,mat_S,mat_K,sumvect,nusim,nujsim,rationusim,&maxrow,&maxcol);
      
      m = 0.;
      for (i= 0; i< maxrow; i++)
        for (j= 0; j< maxcol; j++)
        {
          if (i != j)
          { 
            if (rationusim[i*maxcol+j] >= rationu[i*maxcol+j])
            {pvrationu[m] += 1.;}
          }
          m += 1;
        }
        
        *omegasim = getomega(matgen,tras_X,mat_S,mat_K,sumvect,lambdasim,nusim,omegasim,&maxrow,&maxcol);
      
      m = 0.;
      for (i= 0; i< maxrow; i++)
        for (j= 0; j< maxcol; j++)
        {
          if (i != j)
          { 
            if (omegasim[i*maxcol+j] >= omega[i*maxcol+j])
            {pvomega[m] += 1.;}
          }
          m += 1;
        }
        
        kappasim = getkappa(matgen,tras_X,mat_S,mat_K,sumvect,nusim,nujsim,rationusim,difrationu,&maxrow,&maxcol);
      
      if (kappasim <= kappa)
      {*pvkappa += 1.;}
      
      epsilonsim = getepsilon(matgen,tras_X,mat_S,mat_K,sumvect,nusim,nujsim,nui,&maxrow,&maxcol);
      
      if (epsilonsim <= epsilon)
      {*pvepsilon += 1.;}
      
    }
    
    /* Compute p values for the social reciprocity measures under the specified null hypothesis */
    
    *pvpsi = (*pvpsi+1.)/(*rep+1.);
    *pvphi = (*pvphi+1.)/(*rep+1.);
    *pvdelta = (*pvdelta+1.)/(*rep+1.);
    *pvdc = (*pvdc+1.)/(*rep+1.);
    *pvkappa = (*pvkappa+1.)/(*rep+1.);
    *pvepsilon = (*pvepsilon+1.)/(*rep+1.);
    
    for (i= 0; i< maxrow; i++)
    { 
      pvlambdaj[i] = (pvlambdaj[i]+1)/(*rep+1);
    }
    
    for (i= 0; i< maxrow; i++)
    { 
      pvnuj[i] = (pvnuj[i]+1)/(*rep+1);
    }
    
    m = 0.;
    for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {
        if (i != j)
        {
          pvratiolambda[m] = (pvratiolambda[m]+1)/(*rep+1);
        }
        else
        {
          pvratiolambda[m] = 0.;
        }
        m += 1;
      }
      
      m = 0.;
    for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {
        if (i != j)
        {
          pvrationu[m] = (pvrationu[m]+1)/(*rep+1);
        }
        else
        {
          pvrationu[m] = 0.;
        }
        m += 1;
      }
      
      m = 0.;
    for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {
        if (i != j)
        {
          pvomega[m] = (pvomega[m]+1)/(*rep+1);
        }
        else
        {
          pvomega[m] = 0.;
        }
        m += 1;
      }
      
      
      /* Eliminate the random seed used in the routine */
      
      PutRNGstate();
    
    /* Deallocate all the matrices used in the routine */   
    
    free(matgen);
    free(omegasim);   
    free(rationusim);
    free(ratiolambdasim);
    free(nujsim);
    free(lambdajsim);
    free(nusim);
    free(lambdasim);
    free(omega);   
    free(difrationu);
    free(rationu);
    free(ratiolambda);
    free(nui);
    free(nuj);
    free(lambdaj);
    free(nu);
    free(lambda);
    free(sumvect);
    free(mat_STS);
    free(mat_KTK);
    free(mat_XTX);
    free(tras_K);
    free(tras_S);
    free(mat_K);
    free(mat_S);
    free(dyadc);
    free(mat_pi);
    free(mat_X);
}


/******************************************************************************************************************************/
/******************************************************************************************************************************/


/* Functions for I&SI routine in C */

void DminScomp(double *mat_X, double *minS, int*maxrow)
{
  int i, j, N;
  N = *maxrow;
  for (i= 0; i< N; i++)
  {
    minS[i] = 0;
    for (j= 0; j< N; j++)
    {
      minS[i] += sign(mat_X[i*N+j]-mat_X[j*N+i]);
    }
    minS[i] = -minS[i];  
  }
}
 
int Icomp(double *mat_X, int *maxrow)
{
  int numbinc, i, j, N;
  
  numbinc = 0;
  N = *maxrow;
  for (j= 1; j< N; j++)
   for (i= 0; i< j; i++)
   {  
     if (mat_X[j*N+i] > mat_X[i*N+j])
     {
       numbinc += 1;
     }
   }
  return(numbinc); 
}

int SIcomp(double *mat_X, int *maxrow)
{
  int strinc, i, j, N;
  
  strinc = 0;
  N = *maxrow;
  for (j= 1; j< N; j++)
   for (i= 0; i< j; i++)
   {  
     if (mat_X[j*N+i] > mat_X[i*N+j])
     {
       strinc += j-i;
     }
   }
  return(strinc); 
}

void swaprowcol (double *mat_X, int *k, int *m, int *maxrow)
{
  int i, j, row, col, N;
  double tempX;
  
  N = *maxrow;
  row = *k;
  col = *m;
  for (i= 0; i< N; i++)
  {
    tempX = mat_X[i*N+row];
    mat_X[i*N+row] = mat_X[i*N+col];
    mat_X[i*N+col] = tempX;
  }
  for (j= 0; j< N; j++)
  {
    tempX = mat_X[row*N+j];
    mat_X[row*N+j] = mat_X[col*N+j];
    mat_X[col*N+j] = tempX;
  } 
}

void swapnames (char **namesOrd, int *k, int *m, char *temp)
{
  int row, col;
  
  row = *k;
  col = *m;
  temp = namesOrd[row];
  namesOrd[row] = namesOrd[col];
  namesOrd[col] = temp;
}

void ISImethod (double *X, int *nrow, char **vecNames, int *tries, double *matord,
                char **namord)
{
     
  int i, j, k, m, Imin, SImin, maxrow, maxcol, *order, netincs, newI, newSI,
      t, randi, DiMinSi, DjMinSj, stopiter1, stopiter2, stopiter3;
  double *DminS, *Xc, *bestXthusfar;
  char **names, **namesbest, *temp;
  
   maxrow = maxcol = *nrow;
    
/* Allocate in memory size of temporary matrix Xc */

   Xc = malloc(maxrow*maxcol*sizeof(double));
   if (Xc == NULL)
   {
     error("Null dimension");
   }
   
/* Allocate in memory size of temporary matrix names, namesbest and temp */

   names = (char **)calloc(maxrow, sizeof (char *)); 
   if (names == NULL)
   {
     error("Null dimension");
   }
   
   namesbest = (char **)calloc(maxrow, sizeof (char *)); 
   if (namesbest == NULL)
   {
     error("Null dimension");
   }    
   
   temp = calloc(maxrow, sizeof(char));

   if (temp == NULL)
   {
     error("Null dimension");
   }  

/* Allocate in memory size of temporary matrix DminS */

   DminS = malloc(maxrow*maxcol*sizeof(double));
   if (DminS == NULL)
   {
     error("Null dimension");
   }   
   
/* Allocate in memory size of temporary matrix BestXthusfar */

   bestXthusfar = malloc(maxrow*maxcol*sizeof(double));
   if (bestXthusfar == NULL)
   {
     error("Null dimension");
   }   
   
/* Allocate in memory size of vector order */

   order = calloc(maxrow, sizeof(int));

   if (order == NULL)
   {
     error("Null dimension");
   }

/* Assigning values to Xc and order */
 
   m = 0;
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      {  
         Xc[i*maxcol+j] = X[m];
         m++;
      }      
      
   for (i= 0; i< maxrow; i++)
   {
     order[i] = i;
   }  
    
/* Assigning values to temp and Xc according initial order*/

   /*DminScomp(Xc,DminS,&maxrow);   
   rsort_with_index(DminS,order,maxrow);*/
   /*for (i= 0; i< maxrow; i++)
   Rprintf("%d ",order[i]," ");*/
   
   for (i= 0; i< maxrow; i++)
   {     
     names[i] = vecNames[order[i]];
     namesbest[i] = vecNames[order[i]];
   }
   
/*   for (i= 0; i< maxrow; i++)
   Rprintf("%s ",names[i]," ");
   i=1;
   j=2;
   swapnames(names, &i,&j, temp);
   
      for (i= 0; i< maxrow; i++)
   Rprintf("%s ",names[i]," ");*/
   
  /* for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
      { 
        Xc[i*maxrow+j]=X[order[i]*maxrow+order[j]];
        Xc[j*maxrow+i]=X[order[j]*maxrow+order[i]];
      }*/
      
/*for (i= 0; i< maxrow; i++)
{
   for(j= 0; j< maxcol; j++)
   {
      Rprintf("%10.6lf ",Xc[i*maxcol+j]);
   }
   Rprintf("\n");
}
Rprintf ("\n");*/      

      
   Imin = Icomp(Xc, &maxrow);
   SImin = SIcomp(Xc, &maxrow);
      
   for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
         bestXthusfar[i*maxrow+j] = Xc[i*maxrow+j];
   
   t = 0;
   stopiter1 = 0;
   GetRNGstate();
   
   while (stopiter1 == 0)
   {
     stopiter2 = 0;
     while (stopiter2 == 0)
     {
        stopiter2 = 1;
        for (j= 1; j< maxrow; j++)
           for (i= 0; i< j; i++)
           {  
              if (Xc[j*maxrow+i] > Xc[i*maxrow+j])
              {
                netincs = 0;
                for (k= i; k< j; k++)
                {
                   netincs += sign(Xc[j*maxrow+k]-Xc[k*maxrow+j]);
                }
                if (netincs > 0)
                {
                  swaprowcol(Xc,&i,&j,&maxrow);
                  swapnames(names, &i,&j, temp); 
                  stopiter2 = 0;
                }  
              }       
           }
     }   

     newI = Icomp(Xc,&maxrow);
     newSI = SIcomp(Xc,&maxrow);
     if ( (newI < Imin) || ((newI == Imin) && (newSI < SImin)) )
    {
      Imin = newI;
      SImin = newSI;
      for (i= 0; i< maxrow; i++)
        for (j= 0; j< maxcol; j++)
          bestXthusfar[i*maxrow+j] = Xc[i*maxrow+j];  
      for (i= 0; i< maxrow; i++)
         namesbest[i] = names[i];         
      stopiter1 = 0;
    } else
    {
      t = t+1;
      if ( (SImin > 0) && (t < *tries) )
      {
        for (j=1; j< maxrow; j++)
        {
          k = -1;
          for (i= 0; i< j; i++)
          {
            if (Xc[j*maxrow+i] > Xc[i*maxrow+j]) k = i;
          }
          if (k > -1)
          {
            randi = floor(runif(0,j));
            swaprowcol(Xc,&j,&randi,&maxrow);
            swapnames(names,&j,&randi,temp);
            stopiter1 = 0;
          }
        }
      } else
      {
        stopiter1 = 1;
      }
    }
  }
    for (i= 0; i< maxrow; i++)
      for (j= 0; j< maxcol; j++)
        Xc[i*maxrow+j] = bestXthusfar[i*maxrow+j];   
    for (i= 0; i< maxrow; i++)
       names[i] = namesbest[i];
    stopiter3 = 0;
    while (stopiter3 == 0)
    {
      stopiter3 = 1;
      for (j=1; j< maxrow; j++)
      {
        i = j-1;  
        if (Xc[j*maxrow+i] == Xc[i*maxrow+j])
        { 
          DiMinSi = DjMinSj = 0;
          for (k= 0; k< maxrow; k++)
          {
            DiMinSi += sign(Xc[i*maxrow+k]-Xc[k*maxrow+i]);
            DjMinSj += sign(Xc[j*maxrow+k]-Xc[k*maxrow+j]);
          }

          if (DjMinSj > DiMinSi) 
          {
            swaprowcol(Xc,&i,&j,&maxrow);
            swapnames(names, &i,&j,temp);
            newSI = SIcomp(Xc,&maxrow);
            if (newSI > SImin)
            {
              swaprowcol(Xc,&i,&j,&maxrow);
              swapnames(names, &i,&j,temp);
            } else
            {   
              SImin = newSI;
              stopiter3 = 0;           
            }  
          }       
        }
      }
    }

  m=0;  
  for (i= 0; i< maxrow; i++)
    for(j= 0; j< maxcol; j++)
    {
      matord[m] = Xc[i*maxcol+j];
      m++;
    }

  for (i= 0; i< maxrow; i++)
      namord[i] = names[i];
      
   free(bestXthusfar);
   free(DminS);
   free(temp);
   free(namesbest);  
   free(names);
   free(order);
   free(Xc);    
}


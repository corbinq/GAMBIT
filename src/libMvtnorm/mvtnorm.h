#ifndef _MVT_H_
#define _MVT_H_

#ifdef __cplusplus
extern"C" {


extern void mvtdst_(int* n,
                    int* nu,
                    double* lower,
                    double* upper,
                    int* infin,
                    double* correl,
                    double* delta,
                    int* maxpts,
                    double* abseps,
                    double* releps,
                    double* error,
                    double* value,
                    int* inform);

#endif

  
double pmvnorm(int* n,
               int* nu,
               double* lower,
               double* upper,
               int* infin,      //   if INFIN(I) < 0, Ith limits are (-infinity, infinity);
                                //   if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
                                //   if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
                                //   if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
               double* correl,  //   the correlation coefficient in row I column J of the 
                                //     correlation matrixshould be stored in 
                                //    CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
               double* delta,   // non-central parameter
               int* maxpts,     // param
               double* abseps,  // param
               double* releps,  // param
               double* error,   // estimated abs. error. with 99% confidence interval
               double* value,      // results store here.
               int* inform);    // inform message goes here

double pmvnorm_P(int n,
                 double* bound,
                 double* correlationMatrix, // (2,1), (3,1), (3,2) .....
                 double* error);
double pmvnorm_Q(int n,
                 double* bound,
                 double* correlationMatrix, // (2,1), (3,1), (3,2) .....
                 double* error);

#ifdef __cplusplus
}
#endif
#endif /* _MVT_H_ */

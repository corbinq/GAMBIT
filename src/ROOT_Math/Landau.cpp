#include <math.h>

using namespace std;

// The landau_ccdf routine below was modified from the landau_cdf routine 
// within the CERN ROOT system (https://root.cern.ch)

double landau_ccdf(double x, double mu, double sigma)
{

	// The parameterization below mimics the behavior of
	// PAR = FMStable::setParam(alpha=1,location=mu,logscale=log(scale),pm=0);
	// std_landau_cdf = function(x) FMStable::pEstable(x, PAR);
	// as implemented in the harmonicmeanp R package.
	
	double xi = 0.636619772367581; //   =  pi/2
	double x0 = -0.288607832450766; //  =  digamma(1)/2
	
	// Other comments below are from the original ROOT system landau_cdf code.  
	
   // implementation of landau distribution (from DISLAN)
   //The algorithm was taken from the Cernlib function dislan(G110)
   //Reference: K.S.Kolbig and B.Schorr, "A program package for the Landau
   //distribution", Computer Phys.Comm., 31(1984), 97-111

   static double p1[5] = {0.2514091491e+0,-0.6250580444e-1, 0.1458381230e-1,-0.2108817737e-2, 0.7411247290e-3};
   static double q1[5] = {1.0            ,-0.5571175625e-2, 0.6225310236e-1,-0.3137378427e-2, 0.1931496439e-2};

   static double p2[4] = {0.2868328584e+0, 0.3564363231e+0, 0.1523518695e+0, 0.2251304883e-1};
   static double q2[4] = {1.0            , 0.6191136137e+0, 0.1720721448e+0, 0.2278594771e-1};

   static double p3[4] = {0.2868329066e+0, 0.3003828436e+0, 0.9950951941e-1, 0.8733827185e-2};
   static double q3[4] = {1.0            , 0.4237190502e+0, 0.1095631512e+0, 0.8693851567e-2};

   static double p4[4] = {0.1000351630e+1, 0.4503592498e+1, 0.1085883880e+2, 0.7536052269e+1};
   static double q4[4] = {1.0            , 0.5539969678e+1, 0.1933581111e+2, 0.2721321508e+2};

   static double p5[4] = {0.1000006517e+1, 0.4909414111e+2, 0.8505544753e+2, 0.1532153455e+3};
   static double q5[4] = {1.0            , 0.5009928881e+2, 0.1399819104e+3, 0.4200002909e+3};

   static double p6[4] = {0.1000000983e+1, 0.1329868456e+3, 0.9162149244e+3,-0.9605054274e+3};
   static double q6[4] = {1.0            , 0.1339887843e+3, 0.1055990413e+4, 0.5532224619e+3};

   static double a1[4] = {0              ,-0.4583333333e+0, 0.6675347222e+0,-0.1641741416e+1};
   static double a2[4] = {0              , 1.0            ,-0.4227843351e+0,-0.2043403138e+1};

   double v = (x - mu)/sigma;
   v = (v - x0)/xi;
   double u;
   double lan;

   if (v < -5.5)
   {
	  u   = exp(v+1);
	  lan = 0.3989422803*exp(-1./u)*sqrt(u)*(1+(a1[1]+(a1[2]+a1[3]*u)*u)*u);
   }
   else if (v < -1 )
   {
	  u   = exp(-v-1);
	  lan = (exp(-u)/sqrt(u))*(p1[0]+(p1[1]+(p1[2]+(p1[3]+p1[4]*v)*v)*v)*v)/
										(q1[0]+(q1[1]+(q1[2]+(q1[3]+q1[4]*v)*v)*v)*v);
   }
   else if (v < 1)
	  lan = (p2[0]+(p2[1]+(p2[2]+p2[3]*v)*v)*v)/(q2[0]+(q2[1]+(q2[2]+q2[3]*v)*v)*v);
   
   else if (v < 4)
	  lan = (p3[0]+(p3[1]+(p3[2]+p3[3]*v)*v)*v)/(q3[0]+(q3[1]+(q3[2]+q3[3]*v)*v)*v);
   
   else if (v < 12)
   {
	  u   = 1./v;
	  lan = (p4[0]+(p4[1]+(p4[2]+p4[3]*u)*u)*u)/(q4[0]+(q4[1]+(q4[2]+q4[3]*u)*u)*u);
   }
   else if (v < 50)
   {
	  u   = 1./v;
	  lan = (p5[0]+(p5[1]+(p5[2]+p5[3]*u)*u)*u)/(q5[0]+(q5[1]+(q5[2]+q5[3]*u)*u)*u);
   }
   else if (v < 300)
   {
	  u   = 1./v;
	  lan = (p6[0]+(p6[1]+(p6[2]+p6[3]*u)*u)*u)/(q6[0]+(q6[1]+(q6[2]+q6[3]*u)*u)*u);
   }
   else
   {
	  // modified to provide more precise tail probabilities
	  
	  //u   = 1./(v-v*log(v)/(v+1));
	  double log_u = (-1.0) * log(v - v*log(v)/(v+1));
	  lan = (a2[2]+a2[3]*exp(log_u));
	  lan = (1.0+lan*exp(log_u));
	  lan = 1.0 - lan*exp(log_u);
	  if( lan == 1.0 ){
		  return exp(log_u);
	  }
   }
   return 1 - lan;
}

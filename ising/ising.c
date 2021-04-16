/* 
   C program: 2d Ising model monte carlo simulations.
*/

#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdbool.h>

/* Set these things by hand */


#define INITIAL_TEMP 0 /*  0 => (T=0,  start)
                           1 => (T=oo, start)
						   2 => (T=0,  strip geometry, start)  */
#define MCS 550
#define LENGTH 37
#define Tc     2.2691853142130216092  /* Critical temperature */
#define N 10
#define STEP 0.005

float TEMP = 2.4;  	/* Temperature in units of interaction and k_B */
/* transition probabilities for the glauber rule, better to make them global, more efficient, don't have to compute tanh every time you update... */
float glauber_p1;// = 0.5*(1-tanh(4/TEMP));
float glauber_p2;// = 0.5*(1-tanh(2/TEMP));
float glauber_p3;// = 0.5;
float glauber_p4;// = 0.5*(1+tanh(2/TEMP));
float glauber_p5;// = 0.5*(1+tanh(4/TEMP));
                       
/* Subroutines for boundary conditions, initial conditions, 
   Monte Carlo moves, random number generator, and movie frames. */

void	initialize_boundary_conditions( int [] , int []);
void	initialize_spin_configuration( int [][LENGTH]);
void	mcmove( int[][LENGTH] , int [] , int [] );
void	mcmove_glauber( int[][LENGTH], int [] , int [] );
void	frame_xterm ( int [][LENGTH] );
float	magnetization_per_spin( int[LENGTH][LENGTH]);
double	ran3( long int *);

int main() 
{
	int itime, i;
	int spin[LENGTH][LENGTH];	// the 2d grid of spins
	int nbr1[LENGTH];		// used for boundary conditions
	int nbr2[LENGTH];		// used for boundary conditions

	float mperspin[MCS];	// magnetization per spin at all times of a single run
	float mperspin_data[N][MCS];	// large matrix of mag per spin
	float avg;				// dummy variable, used when computing average mperspin / time
	float mperspin_avg[MCS];	// mperspin_data averaged along 0 axis

	long int iseed;      // random number stuff
	float junk;
	double runtime_seconds; 
	int time_seed;		// random number stuff
	
	bool verbose=false; // if verbose is set to false, program only outputs useful data, to be piped to a data file for analysis
	// if verbose is set to true the script makes lots of print statements, for debugging purposes

	// initiate seed with time in micro seconds
	struct timeval start, stop;
	gettimeofday(&start,NULL); // for the seed but also time the runtime
	time_seed = (int)(-start.tv_usec) - (int)(1000000 * start.tv_sec);
	if (time_seed%2==0) time_seed += 1; // needs to be negative odd integer, not sure why

	//iseed = -12888333;   // any negative odd integer (deterministic seed)
	iseed = time_seed;   // any negative odd integer (time based seed)
	junk = ran3(&iseed);  

	/* Initialize cold (all magnetized T = 0), or, hot (random, T = oo), 
	   or, cold strip (magnetized T = 0 strip in middle, to look at 
	   interface fluctuations):  INITIAL_TEMP = 0 => (T=0)
								 INITIAL_TEMP = 1 => (T=oo)
								 INITIAL_TEMP = 2 => (T=0, strip) */
	initialize_spin_configuration(spin);		// initial temp is zero, all spins up 
	initialize_boundary_conditions(nbr1, nbr2);	// Boundary conditions (periodic) 

	/* Start */
	itime = 0;
	if (verbose){
		printf("START. itime is %d, Temp is %fTc \n", itime, TEMP/Tc);
	}

	for (TEMP=2.53 ; TEMP>2.35; TEMP-=STEP){ // we need epsilon=0.03 or blow up in tau
		/* glauber rule probability updates (global variables) */
		glauber_p1 = 0.5*(1-tanh(4/TEMP));
		glauber_p2 = 0.5*(1-tanh(2/TEMP));
		glauber_p3 = 0.5;
		glauber_p4 = 0.5*(1+tanh(2/TEMP));
		glauber_p5 = 0.5*(1+tanh(4/TEMP));
 
		// do N monte simulations 
		for (i=0; i<N; i++){
			initialize_spin_configuration(spin);
			mperspin_data[i][0] = magnetization_per_spin(spin);
			for (itime=1;itime<MCS; itime++){
				mcmove_glauber(spin, nbr1, nbr2);
				mperspin_data[i][itime] = magnetization_per_spin(spin);
				if (verbose && itime % 100==0 && i % 1 == 0) {
					printf("RUNNING...itime is %d, Temp is %fTc, @ %d'th simulation\n", itime, TEMP/Tc, i);
				}
			}
		}
		// average over N simulations
		for (itime=0; itime<MCS;itime++){
			avg = 0;
			for (i=0;i<N;i++){
				avg += mperspin_data[i][itime];
			}
			avg = avg/N;
			mperspin_avg[itime] = avg;
		}
		
		// display the magnetization per spin
		if (verbose){
			printf("Magnetization per spin, samples over time\n");
			for (itime=0; itime<MCS; itime+=15){
				printf("%f\t",mperspin_avg[itime]);
			}		
			gettimeofday(&stop,NULL);
			runtime_seconds = (double)(stop.tv_sec - start.tv_sec) + (double)(stop.tv_usec - start.tv_usec)/1000000;
			printf("\nRuntime in seconds %f\n",runtime_seconds);
		}
		else // just print the data, in csv format, first element of each row is the temp
		{
			printf("%f,",TEMP);
			for (itime=0;itime<MCS-1;itime++){
				printf("%f,",mperspin_avg[itime]);
			}
			printf("%f",mperspin_avg[MCS-1]);
		}
		printf("\n");
	}
}

// determine the magnetization per spin
float magnetization_per_spin ( int spin[LENGTH][LENGTH]){ 
	int i,j;
	int sum=0;
	for (i=0 ; i<LENGTH ; i++) {
		for (j=0 ; j<LENGTH ; j++) {
			sum += spin[i][j];
		}
	}
	return (float) sum / (LENGTH*LENGTH);
}

void initialize_boundary_conditions (int nbr1[], int nbr2[]){
	int i;
   // Periodic boundary conditions 
	for (i = 0 ; i < LENGTH ; i++) { 
       nbr1[i] = i - 1;
       nbr2[i] = i + 1;
       if (i == 0 )         nbr1[i] = LENGTH - 1;
       if (i == LENGTH - 1) nbr2[i] = 0;
     }

}

void initialize_spin_configuration( int spin[][LENGTH] ) {
	int i, ix, iy;
	long int idum;
	int zero_one, plus_minus;
	int bottom_of_strip , top_of_strip;
 
	if(INITIAL_TEMP==0) {
   // start magnetized all spins = 1 
     for (ix = 0 ; ix < LENGTH; ix++) {  
		 for (iy = 0 ; iy < LENGTH; iy++) {
			spin[ix][iy] = 1;
		 }
     }
 }
 else if(INITIAL_TEMP==1) {
   /* Start all spins random plus or minus.  a little awkward
      coding with plus_minus because ran3 is real */

     for (ix = 0 ; ix < LENGTH; ix++){    
       for (iy = 0 ; iy < LENGTH; iy++){       
			 zero_one     = 2*ran3(&idum);     // random 0 or 1  
			 plus_minus   = 2*zero_one - 1;    // now random +/-1
			 spin[ix][iy] = plus_minus;
	   }
	 }
 }
 else {  
   /*  Start magnetized strip of spins = 1 in middle only.
       Note this is what happens if the INITIAL_TEMP is 
	   not zero or one. */    

	bottom_of_strip = 0.333 * LENGTH;
	top_of_strip    = 0.666 * LENGTH;

     for (ix = 0 ; ix < LENGTH; ix++) {   
       for (iy = 0 ; iy < bottom_of_strip ; iy++) { 
	       spin[ix][iy] = -1;
       }
       for (iy = bottom_of_strip ; iy < top_of_strip ; iy++) {  
	       spin[ix][iy] =  1;
	   }
       for (iy = top_of_strip ; iy < LENGTH ; iy++) {  
	       spin[ix][iy] = -1;
	   }
     }
	}
}

void mcmove_glauber( int spin[][LENGTH], int nbr1[] , int nbr2[]) {
/*
   ONE MONTE CARLO STEP by Glauber:
   FLIP WITH prob1   prob2    0.5    prob4   prob5   (Below spins called)
               +       -       -       -       -            ss1
             + + +   + + +   + + -   + + -   - + -      ss2 ss0 ss4
               +       +       +       -       -            ss3  
*/

  int i, ix, iy;
  int ixpick, iypick;
  int ss0, ss1, ss2, ss3, ss4, de;
  double ran_no;
  long int idum;
  for (i = 1 ; i <= LENGTH*LENGTH ; i++) {
	  ixpick = LENGTH * ran3(&idum) ; 
	  iypick = LENGTH * ran3(&idum) ;  

      ss0 = spin [ixpick]       [iypick]       ;     
      ss1 = spin [nbr1[ixpick]] [iypick]       ;
      ss2 = spin [ixpick]       [nbr1[iypick]] ;
      ss3 = spin [nbr2[ixpick]] [iypick]       ;
      ss4 = spin [ixpick]       [nbr2[iypick]] ;

      de =  2*ss0*(ss1+ss2+ss3+ss4);

		ran_no = ran3(&idum);
		if ((de==8 && ran_no < glauber_p1) || 
			(de==4 && ran_no < glauber_p2) ||
			(de==0 && ran_no < glauber_p3) || 
			(de==-4 && ran_no < glauber_p4) ||
			(de==-8 && ran_no < glauber_p5))
		{
			spin[ixpick][iypick] *= -1;
		}

}
}


void mcmove( int spin[][LENGTH], int nbr1[] , int nbr2[]) {
/*
   ONE MONTE CARLO STEP by Metropolis: Flip probability 1 if Enew < Eold, 
   else prob is exp -(Enew-Eold)/T.  Simplified here since only there 
   are five cases in d=2 for external field = 0.
   FLIP WITH prob1   prob2    1.0     1.0     1.0   (Below spins called)
               +       -       -       -       -            ss1
             + + +   + + +   + + -   + + -   - + -      ss2 ss0 ss4
               +       +       +       -       -            ss3  
*/

  int i, ix, iy;
  int ixpick, iypick;
  int ss0, ss1, ss2, ss3, ss4, de;
  int flag;
  long int idum;
  double prob1 , prob2;
  prob1 = exp(-8.0/TEMP);
  prob2 = exp(-4.0/TEMP);
  for (i = 1 ; i <= LENGTH*LENGTH ; i++) {
	  ixpick = LENGTH * ran3(&idum) ; // this is Martin Grant's pseudo random stuff, always gives same output 
	  iypick = LENGTH * ran3(&idum) ;	  
	  //ixpick = LENGTH * rand();
	  //iypick = LENGTH * rand();

      ss0 = spin [ixpick]       [iypick]       ;     
      ss1 = spin [nbr1[ixpick]] [iypick]       ;
      ss2 = spin [ixpick]       [nbr1[iypick]] ;
      ss3 = spin [nbr2[ixpick]] [iypick]       ;
      ss4 = spin [ixpick]       [nbr2[iypick]] ;

	  //printf("ixpick=%d\tiypick=%d\n",ixpick,iypick);// trace

      de =  2*ss0*(ss1+ss2+ss3+ss4);

      flag = 1;                     // flip spin if flag = 1

             if ( ( (de == 8) && (ran3(&idum) > prob1) )
			||  
         	  ( (de == 4) && (ran3(&idum) > prob2) )  )   
				flag = 0;
	 
		spin[ixpick][iypick] = (1 - 2*flag )*spin[ixpick][iypick];
}
}

#include <stdlib.h>         
#define MBIG 1000000000
#define MSEED 161803398              // Portable random number generator
#define MZ 0                       
#define FAC (1.0/MBIG)
 
int time_usec(){ // returns the time in micro seconds mod 1000 seconds (i think)
	
}

double ran3(long *idum) {
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  
  if (*idum < 0 || iff == 0){
      iff=1;
      mj = labs(MSEED-labs(*idum));
      mj %= MBIG;
      ma[55]=mj;
      mk=1;
      for (i=1;i<=54;i++){
          ii=(21*i) % 55;
          ma[ii]=mk;
          mk=mj-mk;
          if (mk < MZ) mk += MBIG;
          mj =ma[ii];
      }
      for (k=1;k<=4;k++)
          for (i=1;i<=55;i++) {
              ma[i] -= ma[1+(i+30) % 55];
              if (ma[i] < MZ) ma[i] += MBIG;
          }
      inext=0;
      inext=31;
      *idum=1;
  }
  if (++inext == 56) inext=1; 
  if (++inextp == 56) inextp=1; 
  mj =ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext] =mj;
  return mj*FAC;
}

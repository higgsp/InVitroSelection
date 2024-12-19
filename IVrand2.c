#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define ROUNDS 20 //rounds of selection
#define N 50

int pois(double mu);
double binom(double qq,double SS);
double gaus(double mu,double sigma);

int iseed=4402;
double eps=0.0; //1e-4;
double xmin=1e-10;

int main(){
	int i,j,m,n,ntypes;
	double Q,xsum,delta;
	double x[ROUNDS+1][N+1][3];
	double xp[N+1][3],A[N+1][3],k[N+1][3];
	double tau[ROUNDS+1];
	double freq,rat,rat2,kmean,kA;
	int RR;

	srand48(iseed);
	A[0][0]=1.0;
	k[0][0]=0.001;
	for(i=1;i<=N;i++) {
		A[i][0]=0.8+drand48()*0.2;
		A[i][1]=A[i][0] - 0.2*drand48();
		A[i][2]=A[i][1] - 0.2*drand48();
		k[i][0]=0.5+0.5*drand48();
		k[i][1]=k[i][0]/(2+4*drand48());
		k[i][2]=k[i][1]/(2+4*drand48());
	}
	
	for(j=1;j<=ROUNDS;j++) {
		//tau[j]=1.0;
		tau[j]=5.0;
		if(j>5) tau[j]=tau[j-1]/5;
	}
	Q=0.2;
	
	for(i=1;i<=N;i++) {
		x[0][i][2]=1e-14;
		x[0][i][1]=0.0;
		x[0][i][0]=0.0;
	}
	x[0][0][0]=1-N*1e-14;
	x[0][0][1]=0.0;
	x[0][0][2]=0.0;
	
	printf("rounds ");
	for(i=0;i<=N;i++) printf("x%d ",i);
	printf("\n");
	j=0;
	printf("%d ",j);
	for(i=0;i<=N;i++) printf("%e ",x[j][i][0]+x[j][i][1]+x[j][i][2]);
	printf("\n");
	for(j=1;j<=ROUNDS;j++) {
		xsum=0.0;
		for(i=0;i<=N;i++) {
			xp[i][0]=x[j-1][i][0]*Q*A[i][0]*(1-exp(-k[i][0]*tau[j]));
			xp[i][1]=x[j-1][i][1]*Q*A[i][1]*(1-exp(-k[i][1]*tau[j]));
			xp[i][2]=x[j-1][i][2]*Q*A[i][2]*(1-exp(-k[i][2]*tau[j]));
			xsum+=(xp[i][0]+xp[i][1]+xp[i][2]);
		}

		for(i=0;i<=N;i++) {
			x[j][i][0]=xp[i][0]/xsum;
			x[j][i][1]=xp[i][1]/xsum;
			x[j][i][2]=xp[i][2]/xsum;
		}
		
		for(i=1;i<=N;i++) {
			if(x[j][i][2]>xmin) {
				delta=x[j][i][2]*eps;
				x[j][i][1]+=delta;
				x[j][i][2]-=delta;
			}
			if(x[j][i][1]>xmin) {
				delta=x[j][i][1]*eps;
				x[j][i][0]+=delta;
				x[j][i][1]-=delta;
			}
		}
		printf("%d ",j);
		for(i=0;i<=N;i++) printf("%e ",x[j][i][0]+x[j][i][1]+x[j][i][2]);
		printf("\n");
	}
	
	RR=12;
	printf("i k kmean kA freq rat rat2 R=%d\n",RR);
	for(i=1;i<=N;i++) {
		freq=x[RR][i][0]+x[RR][i][1]+x[RR][i][2];
		rat=freq/(x[RR-1][i][0]+x[RR-1][i][1]+x[RR-1][i][2]);
		rat2=(x[RR+1][i][0]+x[RR+1][i][1]+x[RR+1][i][2])/freq;
		kmean=(k[i][0]*x[RR][i][0]+k[i][1]*x[RR][i][1]+k[i][2]*x[RR][i][2])/freq;
		kA=(k[i][0]*A[i][0]*x[RR][i][0]+k[i][1]*A[i][1]*x[RR][i][1]+k[i][2]*A[i][2]*x[RR][i][2])/freq;
		printf("%d %f %f %f %e %f %f\n",i,k[i][0],kmean,kA,freq,rat,rat2);
	}
	RR=10;
	printf("i k kmean kA freq rat rat2 R=%d\n",RR);
	for(i=1;i<=N;i++) {
		freq=x[RR][i][0]+x[RR][i][1]+x[RR][i][2];
		rat=freq/(x[RR-1][i][0]+x[RR-1][i][1]+x[RR-1][i][2]);
		rat2=(x[RR+1][i][0]+x[RR+1][i][1]+x[RR+1][i][2])/freq;
		kmean=(k[i][0]*x[RR][i][0]+k[i][1]*x[RR][i][1]+k[i][2]*x[RR][i][2])/freq;
		kA=(k[i][0]*A[i][0]*x[RR][i][0]+k[i][1]*A[i][1]*x[RR][i][1]+k[i][2]*A[i][2]*x[RR][i][2])/freq;
		printf("%d %f %f %f %e %f %f\n",i,k[i][0],kmean,kA,freq,rat,rat2);
	}
	RR=8;
	printf("i k kmean kA freq rat rat2 R=%d\n",RR);
	for(i=1;i<=N;i++) {
		freq=x[RR][i][0]+x[RR][i][1]+x[RR][i][2];
		rat=freq/(x[RR-1][i][0]+x[RR-1][i][1]+x[RR-1][i][2]);
		rat2=(x[RR+1][i][0]+x[RR+1][i][1]+x[RR+1][i][2])/freq;
		kmean=(k[i][0]*x[RR][i][0]+k[i][1]*x[RR][i][1]+k[i][2]*x[RR][i][2])/freq;
		kA=(k[i][0]*A[i][0]*x[RR][i][0]+k[i][1]*A[i][1]*x[RR][i][1]+k[i][2]*A[i][2]*x[RR][i][2])/freq;
		printf("%d %f %f %f %e %f %f\n",i,k[i][0],kmean,kA,freq,rat,rat2);
	}
	exit(0);
}
//-----------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define LMAX 80
#define N 39

int motifhss3(int seq[LMAX+1],int len);
int motifhss1(int seq[LMAX+1],int len);

double eps = 0.0;
long npop=1e8;
int motseq[LMAX+1];
int mismatch[4][4];

int main(){
	int i,j,k,len,m1,m2,m,s,b,iseed;
	long iseq,nmatch,nsave;

	double fact[LMAX+1];
	double pobs[22],hobs[LMAX+1],hth[LMAX+1];
	int seq[LMAX+1],W;
	int saveseq[2000][LMAX+1];
	double q[22],ps[14],pb[11],G[22],PW[22];
	double sum;
	
	for(i=0;i<4;i++) {
	for(j=0;j<4;j++) {
		mismatch[i][j]=1;
	}}
	mismatch[0][3]=0;
	mismatch[1][2]=0;
	mismatch[2][1]=0;
	mismatch[3][0]=0;
	mismatch[2][3]=0;
	mismatch[3][2]=0;
	
	iseed = 997;
	srand48(iseed);
	fact[0]=1.0;
	for(m=1;m<=LMAX;m++) fact[m]=fact[m-1]*(double)m;		
	
	s=13;
	b=8;
	for(m1=0;m1<=s;m1++) ps[m1]=pow(0.25,s-m1)*pow(0.75,m1)*fact[s]/(fact[m1]*fact[s-m1]);
	for(m2=0;m2<=b;m2++) pb[m2]=pow(0.375,b-m2)*pow(0.625,m2)*fact[b]/(fact[m2]*fact[b-m2]);
	
	for(m=0;m<=s+b;m++) q[m]=0.0;
	for(m1=0;m1<=s;m1++) {
	for(m2=0;m2<=b;m2++) {
		q[m1+m2]+=(ps[m1]*pb[m2]);
	}}
	W=(LMAX-N+3)*(LMAX-N+2)*(LMAX-N+1)/6;
//	W=LMAX-N+1;
	
	G[s+b]=q[s+b];
	for(m=s+b-1;m>=0;m--) G[m]=G[m+1]+q[m];
	for(m=0;m<s+b;m++) PW[m]=pow(G[m],W)-pow(G[m+1],W);
	PW[s+b]=pow(G[s+b],W);
	
/*	sum=0.0;
	for(m=0;m<=s+b;m++) {
		sum+=PW[m];
		printf("%d %e %e\n",m,q[m],PW[m]);
	}
	printf("\nsum %e\n\n",sum);*/
	
	nsave=0;
	for(m1=0;m1<=s+b;m1++) pobs[m1]=0.0;
	for(iseq=1;iseq<=npop;iseq++) {
		for(i=1;i<=LMAX;i++) seq[i]=4*drand48();
		m1=motifhss1(seq,LMAX);
		pobs[m1]+=1.0;
		if(m1<=1&&nsave<1900) {
			nsave++;
			//printf("nsave= %d\n",nsave);
			for(i=1;i<=LMAX;i++) saveseq[nsave][i]=seq[i];
		}
	}
	printf("m q PW Pobs\n");
	for(m=0;m<=s+b;m++) {
		pobs[m]/=npop;
		printf("%d %e %e %e\n",m,q[m],PW[m],pobs[m]);
	}
	
	for(m=0;m<=LMAX;m++) {
		hobs[m]=0.0;
		hth[m]=pow(0.75,m)*pow(0.25,LMAX-m)*fact[LMAX]/(fact[m]*fact[LMAX-m]);
	}
	for(m1=1;m1<nsave;m1++) {
	for(m2=m1+1;m2<=nsave;m2++) {
		m=0;
		for(i=1;i<=LMAX;i++) if(saveseq[m1][i]!=saveseq[m2][i]) m++;
		hobs[m]+=1.0;
	}}
	printf("\nm hth hobs\n");
	for(m=0;m<=LMAX;m++) {
		hobs[m]/=(0.5*nsave*(nsave-1));
		printf("%d %e %e\n",m,hth[m],hobs[m]);
	}
	
	exit(0);
	
}

//-------------------------------------------------------------------------
int motifhss3(int seq[LMAX+1],int len) {
	int j1,j2,j3,m,mmin,nn;
	int wseq[LMAX+1];
	
	mmin=len;
	for(j1=0;j1<=len-N;j1++) {
	for(j2=0;j2<=len-N-j1;j2++) {
	for(j3=0;j3<=len-N-j1-j2;j3++) {
		for(nn=1;nn<=11;nn++) wseq[nn]=seq[nn+j1];
		for(nn=12;nn<=28;nn++) wseq[nn]=seq[nn+j1+j2];
		for(nn=29;nn<=N;nn++) wseq[nn]=seq[nn+j1+j2+j3];
		m=0;
		m+=mismatch[wseq[1]][wseq[39]];
		m+=mismatch[wseq[2]][wseq[38]];
		m+=mismatch[wseq[3]][wseq[37]];
		//m+=mismatch[wseq[4]][wseq[36]]; //extra pair
		m+=mismatch[wseq[6]][wseq[17]];
		m+=mismatch[wseq[7]][wseq[16]];
		m+=mismatch[wseq[8]][wseq[15]];
		m+=mismatch[wseq[9]][wseq[14]];
		//m+=mismatch[wseq[25]][wseq[32]]; //extra pair
		m+=mismatch[wseq[26]][wseq[31]];
		if(wseq[4]!=3) m++; //extra site
		if(wseq[18]!=1) m++;
		if(wseq[19]!=3) m++;
		if(wseq[20]!=2) m++;
		if(wseq[21]!=0) m++;
		if(wseq[23]!=2) m++;
		if(wseq[24]!=0) m++;
		if(wseq[25]!=2) m++; //extra site
		if(wseq[32]!=1) m++; //extra site
		if(wseq[33]!=2) m++;
		if(wseq[34]!=0) m++;
		if(wseq[35]!=0) m++;
		if(wseq[36]!=0) m++; //extra site
		if(m<mmin) mmin=m;
	}}}
	return mmin;
}
//-------------------------------------------------------------------------
int motifhss1(int seq[LMAX+1],int len) {
	int j1,j2,j3,m,mmin,nn;
	int wseq[LMAX+1];
	
	mmin=len;
	for(j1=0;j1<=len-N;j1++) {
		for(nn=1;nn<=N;nn++) wseq[nn]=seq[nn+j1];
		m=0;
		m+=mismatch[wseq[1]][wseq[39]];
		m+=mismatch[wseq[2]][wseq[38]];
		m+=mismatch[wseq[3]][wseq[37]];
		//m+=mismatch[wseq[4]][wseq[36]]; //extra pair
		m+=mismatch[wseq[6]][wseq[17]];
		m+=mismatch[wseq[7]][wseq[16]];
		m+=mismatch[wseq[8]][wseq[15]];
		m+=mismatch[wseq[9]][wseq[14]];
		//m+=mismatch[wseq[25]][wseq[32]]; //extra pair
		m+=mismatch[wseq[26]][wseq[31]];
		if(wseq[4]!=3) m++; //extra site
		if(wseq[18]!=1) m++;
		if(wseq[19]!=3) m++;
		if(wseq[20]!=2) m++;
		if(wseq[21]!=0) m++;
		if(wseq[23]!=2) m++;
		if(wseq[24]!=0) m++;
		if(wseq[25]!=2) m++; //extra site
		if(wseq[32]!=1) m++; //extra site
		if(wseq[33]!=2) m++;
		if(wseq[34]!=0) m++;
		if(wseq[35]!=0) m++;
		if(wseq[36]!=0) m++; //extra site
		if(m<mmin) mmin=m;
	}
	return mmin;
}
//-------------------------------------------------------------------------
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#define rand01 (0.9999999*double(rand())/RAND_MAX) 
				
#define getrandom(max1) (((rand())%(max1))) // random integer between 0 and max-1
	  


//There are two implementations of polychronous group search algorithm here.
// The old one is event-trigered (very fast). 
// The new one is activity-trigered (slow, but more precise).
// By default, the new definition is used. Uncomment the line below to use the old one.

//#define old_definition_of_polychronous 


//M=100;                 % the number of synapses per neuron
const	int		M = 100; 

//D=10;                  % maximal axonal conduction delay
const	int	D=20;		

//% excitatory neurons   inhibitory neurons   total number 
//Ne=800;                Ni=200;              N=Ne+Ni;

const	int		Ne =800;			
const	int		Ni =200;			 
const	int		N = Ne+Ni;	
		double	C_max=10;		


const	int	W=3;	// initial width of polychronous groups
int		min_group_path = 7;		// minimal length of a group
int		min_group_time = 40;	// minimal duration of a group (ms)



double	a[N], d[N];		//
int		post[N][M];		//
double	s[N][M], sd[N][M];	  //
short	delays[N][D][M];	//
short	delays_length[N][D];  //
int		N_pre[N], I_pre[N][3*M], D_pre[N][3*M];	   //
double	*s_pre[N][3*M], *sd_pre[N][3*M];
double	LTP[N][1001+D], LTD[N];		  //
double	v[N], u[N];	  //
int		N_firings;
const int N_firings_max=100000;
int		firings[N_firings_max][2];

void initialize()
{	int i,j,k,jj,dd, exists, r;

//	a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
	for (i=0;i<Ne;i++) a[i]=0.02;
	for (i=Ne;i<N;i++) a[i]=0.1;

//	d=[   8*ones(Ne,1);      2*ones(Ni,1)];
	for (i=0;i<Ne;i++) d[i]=8;
	for (i=Ne;i<N;i++) d[i]=2;


//	post=ceil([Ne+Ni*rand(Ne,M*Ni/N),Ne*rand(Ne,M*Ne/N);Ne*rand(Ni,M)]); 
	for (i=0;i<Ne;i++)
	{

		for (j=0;j<M;j++) 
		{
			do
			{
				exists = 0;
				r = int(floor(N*rand01));
				if (r == i) exists=1;  
				for (k=0;k<j;k++) if (post[i][k] == r) exists = 1;
			}
			while (exists == 1);
			post[i][j]=r;
		}

	}
	for (i=Ne;i<N;i++)
	{
		for (j=0;j<M;j++) 
		{
			do
			{
				exists = 0;
				r = int(floor(Ne*rand01));
				for (k=0;k<j;k++) if (post[i][k] == r) exists = 1;
			}
			while (exists == 1);
			post[i][j]=r; 
		}
	}

//	s=[6*ones(Ne,M);-5*ones(Ni,M)];     % initial synaptic weights
	for (i=0;i<Ne;i++) for (j=0;j<M;j++) s[i][j]=6; 

	for (i=Ne;i<N;i++) for (j=0;j<M;j++) s[i][j]=-5;

	
//	sd=zeros(N,M);                      % derivatives of synaptic weights
  	for (i=0;i<N;i++)for (j=0;j<M;j++) sd[i][j]=0;
	

//	for i=1:N
  	for (i=0;i<N;i++)
	{
		short ind=0;

//		if i<=Ne
		if (i<Ne)
		{

//			for j=1:D
//            delays{i,j}=M/D*(j-1)+(1:M/D);
//			end;
			for (j=0;j<D;j++) 
			{	delays_length[i][j]=M/D;
				for (k=0;k<delays_length[i][j];k++)
					delays[i][j][k]=ind++;
			}
		}
//		else
//        delays{i,1}=1:M;
//		end;
		else
		{
			for (j=0;j<D;j++) delays_length[i][j]=0;
			delays_length[i][0]=M;
			for (k=0;k<delays_length[i][0];k++)
					delays[i][0][k]=ind++;
		}
	}
	
//		pre{i} = find(post==i & s>0);     % Indeces of pre excitatory neurons
//		aux{i} = N*(D-1-ceil(ceil(pre{i}/N)/(M/D))) + 1+mod(pre{i}-1,N);
  	for (i=0;i<N;i++)
	{
		N_pre[i]=0;
		for (j=0;j<Ne;j++)
		for (k=0;k<M;k++)
		if (post[j][k] == i) 
		{
			I_pre[i][N_pre[i]]=j;
			for (dd=0;dd<D;dd++)
				for (jj=0;jj<delays_length[j][dd];jj++)
					if (post[j][delays[j][dd][jj]]==i) D_pre[i][N_pre[i]]=dd;
			s_pre[i][N_pre[i]]=&s[j][k];
			sd_pre[i][N_pre[i]++]=&sd[j][k];
		}

//	end;
	}

//	LTP = zeros(N,1001+D);
	for (i=0;i<N;i++)
		for (j=0;j<1+D;j++)
			LTP[i][j]=0;

//	LTD = zeros(N,1);
	for (i=0;i<N;i++)	LTD[i]=0;

//	v = -65+10*rand(N,1);               % initial values for v
	for (i=0;i<N;i++)	v[i]=-65+10*rand01;

//	u = 0.2.*v;                         % initial values for u
	for (i=0;i<N;i++)	u[i]=0.2*v[i];

//	firings=[-D 0];                     % spike timings
	N_firings=1;
	firings[0][0]=-D;
	firings[0][1]=0;
}



// ----------------------------------------------------------------------------------

void save_all(char fname[30])
{
	int		i,j,k;
	FILE	*fce;
	fce = fopen(fname,"w");
	std::cout << "\nsaving data \n";
	
	for (i=0; i < N; ++i)
	{
		fprintf(fce, "%5.3f %5.3f \n", v[i], u[i]);
		for (j=0; j < M; ++j)
			fprintf(fce, "%d %5.3f %6.5f \n", post[i][j], s[i][j], sd[i][j]);
		
		
		for (k=0; k < D; ++k)
		{
			fprintf(fce, "%d  ", delays_length[i][k]);
			for (j=0; j < delays_length[i][k]; ++j)
				fprintf(fce, "%d ", delays[i][k][j]);
			fprintf(fce, "\n");
		}
		
		fprintf(fce, "%d  ", N_pre[i]);
		for (j=0; j < N_pre[i]; ++j)
			fprintf(fce, "%d %d ", I_pre[i][j], D_pre[i][j]);
		
		fprintf(fce, "\n %5.4f ", LTD[i]);
		for (j=0; j < D+1; ++j)
			fprintf(fce, "%5.4f ", LTP[i][j]);
		fprintf(fce, "\n");
	}

	fprintf(fce, " %d", N_firings);
	for (i=0; i < N_firings; ++i)
		fprintf(fce, "%d %d ", firings[i][0],firings[i][1]);

	fclose(fce);
}

// ----------------------------------------------------------------------------------

void load_all(char fname[30])
{
//	std::cout << "\nloading data \n";
	int		i,j, k, Np;
	float  x;
	int		dd;
	
	FILE	*stream;
	stream = fopen( fname, "r" );
    if( stream == NULL )
		std::cout << " \n Error: The file " << fname << " cannot be opened \n";
    else
    {
	  /* Set pointer to beginning of file: */
      fseek( stream, 0L, SEEK_SET );
	  for (i=0; i < N; ++i)
	  {
		fscanf( stream, "%f", &x);
		v[i]=x;
		fscanf( stream, "%f", &x);
		u[i]=x;

		for (j=0; j < M; ++j)
		{
			fscanf( stream, "%d", &dd);
			post[i][j]=dd;
			fscanf( stream, "%f", &x);
			s[i][j]=x;
			fscanf( stream, "%f", &x);
			sd[i][j]=x;
		}
		for (k=0; k < D; ++k)
		{
			fscanf( stream, "%d", &dd);
			delays_length[i][k]=dd;
			for (j=0; j < delays_length[i][k]; ++j)
			{
				fscanf( stream, "%d", &dd);
				delays[i][k][j]=dd;
			}
		}

		fscanf( stream, "%d", &dd);
		N_pre[i] = dd;
	    for (j=0; j < N_pre[i]; ++j)
		{
		  fscanf( stream, "%d", &dd);
		  I_pre[i][j]=dd;
		  fscanf( stream, "%d", &dd);
		  D_pre[i][j]=dd;
		}

		fscanf( stream, "%f", &x);
		LTD[i]=x;
		for (j=0; j < D+1; ++j)
		{
			fscanf( stream, "%f", &x);
			LTP[i][j]=x;
		}
	  }
	  
	  fscanf( stream, "%d", &dd);
	  N_firings=dd;
	  for (i=0; i < N_firings; ++i)
	  {
			fscanf( stream, "%d", &dd);
			firings[i][0]=dd;
			fscanf( stream, "%d", &dd);
			firings[i][1]=dd;
	  }
	
	  fclose( stream );

	  for (i=0; i < N; ++i)
	  {
		for (Np=0;Np<N_pre[i];Np++)
		{
			j = I_pre[i][Np];
			for (k=0;k<M;k++)
			if (post[j][k] == i) 
			{
				s_pre[i][Np]=&s[j][k];
				sd_pre[i][Np++]=&sd[j][k];
			}
		}
	  }
	  
	}
}






//--------------------------------------------------------------
int			N_polychronous;


double		C_rel = 0.95*C_max;
const int	polylenmax = N;
int			N_postspikes[polylenmax], I_postspikes[polylenmax][N], J_postspikes[polylenmax][N], D_postspikes[polylenmax][N], L_postspikes[polylenmax][N];
double		C_postspikes[polylenmax][N];
int			N_links, links[2*W*polylenmax][4];
int			group[polylenmax], t_fired[polylenmax], layer[polylenmax];
int			gr3[W], tf3[W];
int			I_my_pre[3*M], D_my_pre[3*M], N_my_pre;
int			N_fired;


FILE		*fpoly;





#ifdef old_definition_of_polychronous


//This is the old definition of polychronous groups
//--------------------------------------------------------------
const	int	V=2;			// number of spikes needed to fire a cell
		int	latency = 4;	// latency from the moment spikes arrive to the moment the neuron fires
		int	tf_final=2;
		int	inh_span=5;
		int	slacke = 0;
		int	slacki = 0;
int			pre_list[3*M], del_list[3*M];
int			inh_list[3*Ni*M],t_inh_list[3*Ni*M], N_inh_list; 

//--------------------------------------------------------------
void	polychronous(int nnum)
{
	int	i,j, t, tf, p, k;
	int npre[W];
	int d, exc, NL, NL_max;
	int	t_fire, t_last, timing;
	int	N_fired, Dmax, already; 
	int	first, last;
	int	used[W], discard;

	

	N_my_pre = 0;
	for (i=0;i<N_pre[nnum];i++)
	if (*s_pre[nnum][i] > C_rel) 
	{
		I_my_pre[N_my_pre]=I_pre[nnum][i];
		D_my_pre[N_my_pre]=D_pre[nnum][i];
		N_my_pre++;
	}
	if (N_my_pre<W) return;

	for (i=0;i<W;i++)	npre[i]=i;

	while (0==0) 
	{
		Dmax=0;
		for (i=0;i<W;i++) if (Dmax < D_my_pre[npre[i]]) Dmax=D_my_pre[npre[i]];
		
		for (i=0;i<W;i++)
		{
			group[i]=I_my_pre[npre[i]];
			t_fired[i]= Dmax-D_my_pre[npre[i]];
			
			for (d=0; d<D; d++)		 
			for (j=0; j<delays_length[group[i]][d]; j++)
			{
				p = post[group[i]][delays[group[i]][d][j]];
				if (((s[group[i]][delays[group[i]][d][j]] > C_rel) & (d>=D_my_pre[npre[i]])) | (p >=Ne))
				{
					timing = t_fired[i]+d+1;
					J_postspikes[timing][N_postspikes[timing]]=group[i];				// presynaptic
					D_postspikes[timing][N_postspikes[timing]]=d;						// delay
					L_postspikes[timing][N_postspikes[timing]]=1;						// layer
					I_postspikes[timing][N_postspikes[timing]++]=p;// index of post target	
				}
			}
		}

		
		N_inh_list=0;
		NL_max = 1;
		N_links = 0;
		N_fired=W;
		t_last = D+D+latency+1;
		for (t=0;t<t_last;t++)
		while (N_postspikes[t] >0)
		{
			N_postspikes[t]--;

			already = 0;

			if (I_postspikes[t][N_postspikes[t]] <Ne)
			{
				for (j=0;j<N_fired;j++)
					if ((I_postspikes[t][N_postspikes[t]] == group[j]) & (t_fired[j] > t-20)) already = 2;
				for (j=0;j<N_inh_list;j++)
					if ((inh_list[j]==I_postspikes[t][N_postspikes[t]]) & (t_inh_list[j]<=t)&(t_inh_list[j]+inh_span>=t)) already++;
			}
			else
			{
				for (j=0;j<N_fired;j++)
					if ((I_postspikes[t][N_postspikes[t]] == group[j]) & (t_fired[j] > t-20)) already = 2;
			}
 
			if ((already<=0) & (N_fired < polylenmax))
			{
				NL = 0;
				
				exc=0;
				first = -1;
				t_fire = t+4;
				for (tf=0;tf<=tf_final;tf++)
				for (j=0;j<=N_postspikes[t+tf]-0.5*tf;j++)
				if (I_postspikes[t+tf][j]==I_postspikes[t][N_postspikes[t]]) 
				{
					if (first < 0) first = tf;
					last = tf;
					pre_list[exc]=J_postspikes[t+tf][j];
					del_list[exc]=D_postspikes[t+tf][j];
					exc++;
					t_fire = t+tf;
					if (NL <L_postspikes[t+tf][j]) NL=L_postspikes[t+tf][j];
				}
				
				if (((first+1>=last) & (exc>=V-slacke)) | (((I_postspikes[t][N_postspikes[t]]<Ne)&(exc>=V)) | ((I_postspikes[t][N_postspikes[t]]>=Ne)&(first+2>=last)&(exc>=V-slacki))))
				{
					i = N_fired++;
	 				group[i]= I_postspikes[t][N_postspikes[t]];
					t_fired[i]= t_fire+latency;
					if (exc==2) t_fired[i]=t_fired[i]+latency; // longer latencies for weaker inputs

					for (j=0;j<exc;j++)
					{
						links[N_links][0]=pre_list[j];
						links[N_links][1]=group[i];
						links[N_links][2]=del_list[j];
						links[N_links++][3]=NL;

					}

					if (group[i] <Ne)
					{
						for (d=0; d<D; d++)
						for (j=0; j<delays_length[group[i]][d]; j++)
						if (s[group[i]][delays[group[i]][d][j]] > C_rel) 
						{
							timing = t_fired[i]+d+1;
							J_postspikes[timing][N_postspikes[timing]]=group[i];				// presynaptic
							D_postspikes[timing][N_postspikes[timing]]=d;						// delay
							L_postspikes[timing][N_postspikes[timing]]=NL+1;					// layer
							I_postspikes[timing][N_postspikes[timing]++]=post[group[i]][delays[group[i]][d][j]];// index of post target	
						}
					}
					else
					{
						for (d=0; d<D; d++)
						for (j=0; j<delays_length[group[i]][d]; j++)
						{
							inh_list[N_inh_list]= post[group[i]][delays[group[i]][d][j]];
							t_inh_list[N_inh_list++]= t_fired[i]+d;
						}

						j=0;
						while (t_inh_list[j]+inh_span<t) j++;
						if (j>10)	// needs cleaning
						{
							for (k=j;k<N_inh_list;k++)
							{
								inh_list[k-j]= inh_list[k];
								t_inh_list[k-j]= t_inh_list[k];
							}
							N_inh_list=N_inh_list-j;
						}
					}
					
					if (NL > NL_max) NL_max = NL;
				
					if (t_last < timing+1) 
					{
						t_last = timing+1;
						if (t_last > polylenmax-D-latency-tf_final-1) t_last = polylenmax-D-latency-tf_final-1;
					}
				}
			}
		}
		if (t_last == polylenmax-D) for (d=t_last;d<polylenmax;d++) N_postspikes[d]=0;


		if ((NL_max>=min_group_path)) 
		{
			discard = 0;
			for (i=0;i<W;i++)
			{
				used[i]=0;
				for (j=0;j<N_links;j++) if ((links[j][0] == group[i]) & (links[j][1] < Ne)) used[i]++;
				if (used[i] == 1) discard = 1;
			}

			if (discard == 0)
			{
				N_polychronous++;
				std::cout << "\ni= " << nnum << ", N_polychronous= " << N_polychronous << ", N_fired = " << N_fired << ",  NL = " << NL_max << ", T=" << t_fired[N_fired-1];
				fprintf(fpoly, " %d  %d,       ", N_fired, NL_max);
				for (j=0; j<N_fired; j++)
					fprintf(fpoly, " %d %d, ", group[j], t_fired[j]);
				fprintf(fpoly, "        ");
				for (j=0; j<N_links; j++)
					fprintf(fpoly, " %d %d %d %d,  ", links[j][0], links[j][1], links[j][2]+1, links[j][3]);
				fprintf(fpoly, "\n");
			}
		}

		i=1;
		while (++npre[W-i] > N_my_pre-i) if (++i > W) return; 
		while (i>1) {npre[W-i+1]=npre[W-i]+1; i--;}
	}	
}

#else

// The new (default) algorithm to find polychronous groups

const	int	latency = D; // maximum latency 


//--------------------------------------------------------------
void	polychronous(int nnum)
{
	int	i,j, t, p, k;
	int npre[W];
	int dd;
	int	t_last, timing;
	int	Dmax, L_max; 
	int	used[W], discard;

	double v[N],u[N],I[N];
	

	N_my_pre = 0;
	for (i=0;i<N_pre[nnum];i++)
	if (*s_pre[nnum][i] > C_rel) 
	{
		I_my_pre[N_my_pre]=I_pre[nnum][i];
		D_my_pre[N_my_pre]=D_pre[nnum][i];
		N_my_pre++;
	}
	if (N_my_pre<W) return;

	for (i=0;i<W;i++)	npre[i]=i;

	while (0==0) 
	{
		Dmax=0;
		for (i=0;i<W;i++) if (Dmax < D_my_pre[npre[i]]) Dmax=D_my_pre[npre[i]];
		
		for (i=0;i<W;i++)
		{
			group[i]=I_my_pre[npre[i]];
			t_fired[i]= Dmax-D_my_pre[npre[i]];
			layer[i]=1;
			
			for (dd=0; dd<D; dd++)		 
			for (j=0; j<delays_length[group[i]][dd]; j++)
			{
				p = post[group[i]][delays[group[i]][dd][j]];
				if ((s[group[i]][delays[group[i]][dd][j]] > C_rel) & (dd>=D_my_pre[npre[i]]))
				{
					timing = t_fired[i]+dd+1;
					J_postspikes[timing][N_postspikes[timing]]=group[i];				// presynaptic
					D_postspikes[timing][N_postspikes[timing]]=dd;						// delay
					C_postspikes[timing][N_postspikes[timing]]=s[group[i]][delays[group[i]][dd][j]];	// syn weight
					I_postspikes[timing][N_postspikes[timing]++]=p;						// index of post target	
				}
			}
		}

		for (i=0;i<N;i++) {v[i]=-70; u[i]=0.2*v[i]; I[i]=0;};

		N_links = 0;
		N_fired=W;
		t_last = D+D+latency+1;
		t=-1;
		while ((++t<t_last) & (N_fired < polylenmax))
		{
			for (p=0;p<N_postspikes[t];p++) 
			  I[I_postspikes[t][p]]+=C_postspikes[t][p]; 
 
		  	for (i=0;i<N;i++)
			{
				v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);
				v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);
				u[i]+=a[i]*(0.2*v[i]-u[i]);
				I[i]=0;
			}

			for (i=0;i<N;i++) 
			if (v[i]>=30)
			{
				v[i] = -65;
				u[i]+=d[i];

				if (N_fired < polylenmax)
				{
					t_fired[N_fired]= t;
					group[N_fired++]=i;
					for (dd=0; dd<D; dd++)
					for (j=0; j<delays_length[i][dd]; j++)
					if ((s[i][delays[i][dd][j]] > C_rel) | (i>=Ne)) 
					{
						timing = t+dd+1;
						J_postspikes[timing][N_postspikes[timing]]=i;				// presynaptic
						D_postspikes[timing][N_postspikes[timing]]=dd;				// delay
//						L_postspikes[timing][N_postspikes[timing]]=NL+1;			// layer
						C_postspikes[timing][N_postspikes[timing]]=s[i][delays[i][dd][j]];	   // syn weight
						I_postspikes[timing][N_postspikes[timing]++]=post[i][delays[i][dd][j]];// index of post target	
					}
					if (t_last < timing+1) 
					{
						t_last = timing+1;
						if (t_last > polylenmax-D-1) t_last = polylenmax-D-1;
					}
				}
			}
		}
		
		if (N_fired>2*W)
		{
			N_links=0;
			L_max=0;
			for (i=W;i<N_fired;i++)
			{
				layer[i]=0;
				for (p=t_fired[i]; (p>t_fired[i]-latency) & (p>=0); p--)
				for (j=0;j<N_postspikes[p];j++)
				if ((I_postspikes[p][j]==group[i]) & (J_postspikes[p][j]<Ne)) 
				{
				   for (k=0;k<i;k++)
				   if ((group[k]==J_postspikes[p][j]) & (layer[k]+1>layer[i])) layer[i]=layer[k]+1;
				   {
					   links[N_links][0]=J_postspikes[p][j];
					   links[N_links][1]=I_postspikes[p][j];
					   links[N_links][2]=D_postspikes[p][j];
					   links[N_links++][3]=layer[i];
					   if (L_max < layer[i]) L_max = layer[i]; 
				   }
				}
			}
										 
			discard = 0;
			for (i=0;i<W;i++)
			{
				used[i]=0;
				for (j=0;j<N_links;j++) if ((links[j][0] == group[i]) & (links[j][1] < Ne)) used[i]++;
				if (used[i] == 1) discard = 1;
			}

//			if ((discard == 0) & (t_fired[N_fired-1] > min_group_time) )  // (L_max >= min_group_path))
			if ((discard == 0) & (L_max >= min_group_path))
			{

				for (i=0;i<W;i++) {gr3[i]=group[i]; tf3[i]=t_fired[i];};

				N_polychronous++;
				std::cout << "\ni= " << nnum << ", N_polychronous= " << N_polychronous << ", N_fired = " << N_fired << ", L_max = " << L_max << ", T=" << t_fired[N_fired-1];
				fprintf(fpoly, " %d  %d,       ", N_fired, L_max);
				for (i=0; i<N_fired; i++)
					fprintf(fpoly, " %d %d, ", group[i], t_fired[i]);
				fprintf(fpoly, "        ");
				for (j=0;j<N_links;j++)
				   fprintf(fpoly, " %d %d %d %d,  ", links[j][0], links[j][1], links[j][2], links[j][3]);
				fprintf(fpoly, "\n");
			}
		}

  		for (dd=Dmax;dd<t_last;dd++) N_postspikes[dd]=0;
		if (t_last == polylenmax-D) for (dd=t_last;dd<polylenmax;dd++) N_postspikes[dd]=0;

		i=1;
		while (++npre[W-i] > N_my_pre-i) if (++i > W) return; 
		while (i>1) {npre[W-i+1]=npre[W-i]+1; i--;}
	}	
}




#endif




 
//--------------------------------------------------------------
void	all_polychronous()
{
	int	i;
	N_polychronous=0;
	fpoly = fopen("..//polyall.dat","w");
   	for (i=0;i<polylenmax;i++) N_postspikes[i]=0;

	for (i=0;i<Ne;i++) polychronous(i);

	std::cout << "\nN_polychronous=" << N_polychronous << "\n";
	fclose(fpoly);
}




void shuffle()
{
	int i, j, ri, rj;
	double x,y;
	std::cout << "***** scrambling ****";
	for (i=0;i<Ne;i++)
	for (j=0;j<M;j++)
	if (post[i][j] < Ne)
	{
		ri = int(floor(rand01*Ne));
		do 
		{
			rj = int(floor(rand01*M));
		}
		while (post[ri][rj] >= Ne);	 
		x=s[ri][rj];
		y=sd[ri][rj];
		s[i][j]=s[ri][rj];
		sd[i][j]=sd[ri][rj];
		s[ri][rj]=x;
		sd[ri][rj]=y;
	}
}


			   

// --------------------------------------------------------------------------
int main()
{
	int		i,j,k;
	int		sec, t;
	double	I[N];
	FILE	*fs, *fx, *fd;
	
	
	srand(0);
	initialize();



//	for sec=1:60*60*5
	for (sec=0; sec<5*60*60; sec++)
	{
	


//		for t=1:1000                  % simulation of 1 sec
		for (t=0;t<1000;t++)
		{

			for (i=0;i<N;i++) I[i] = 0;
			I[int(floor(N*rand01))]=20;


			for (i=0;i<N;i++) 
//			fired = find(v>=30);          % indices of fired neurons
			if (v[i]>=30)
			{
        

//			    v(fired)=-65;
				v[i] = -65;

//	            u(fired)=u(fired)+d(fired);
				u[i]+=d[i];

//	            LTP(fired,t+D)=0.1;
				LTP[i][t+D]= 0.1;

//	            LTD(fired)=0.12;
				LTD[i]=0.12;

//				for k=1:length(fired)
//					sd(pre{fired(k)}) = sd(pre{fired(k)})+LTP(N*t+aux{fired(k)});
//				end;
				for (j=0;j<N_pre[i];j++) *sd_pre[i][j]+=LTP[I_pre[i][j]][t+D-D_pre[i][j]-1];

//	            firings=[firings; t+zeros(length(fired),1), fired];
				firings[N_firings  ][0]=t;
				firings[N_firings++][1]=i;
				if (N_firings == N_firings_max)
				{
					std::cout << "*** Two many spikes, t=" << t << "*** (ignoring)";
					N_firings=1;
				}
			}

//	        k=size(firings,1);
			k=N_firings-1;

//	        while t-firings(k,1)<D
			while (t-firings[k][0] <D)
			{

//	            del=delays{firings(k,2),t-firings(k,1)+1};
				for (j=0; j< delays_length[firings[k][1]][t-firings[k][0]]; j++)
				{
//					ind = post(firings(k,2),del);
					i=post[firings[k][1]][delays[firings[k][1]][t-firings[k][0]][j]]; 

//					I(ind)=I(ind)+s(firings(k,2), del)';
					I[i]+=s[firings[k][1]][delays[firings[k][1]][t-firings[k][0]][j]];

//					if firings(k,2) <=Ne
					if (firings[k][1] <Ne)
                
//						sd(firings(k,2),del)=sd(firings(k,2),del)-LTD(ind)';
						sd[firings[k][1]][delays[firings[k][1]][t-firings[k][0]][j]]-=LTD[i];
//					end;
				}

//				k=k-1;
				k--;

//		    end;
			}

			for (i=0;i<N;i++)
			{
//		        v = v + 0.5*((0.04*v+5).*v+140-u+I);    % for numerical stability
//			    v = v + 0.5*((0.04*v+5).*v+140-u+I);    % time step is 0.5 ms
				v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);
				v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);

//				u = u + a.*(0.2*v-u);
				u[i]+=a[i]*(0.2*v[i]-u[i]);

//				LTP(:,t+D+1)=0.95*LTP(:,t+D); % tau = 20 ms
				LTP[i][t+D+1]=0.95*LTP[i][t+D];

//				LTD=0.95*LTD;                 % tau = 20 ms
				LTD[i]*=0.95;
			}

//		end;
		}
	
//		frate(end+1)=sum(firings(:,2)<=Ne)/Ne;
		double	frate=0;
		for (i=1;i<N_firings;i++)
			if ((firings[i][0] >=0) && (firings[i][1] <Ne)) frate++;
		frate = frate/Ne;

//		str(end+1) = sum(sum(s(find(post<=Ne)) > 6.3))/Ne;
		double	str=0;
		for (i=0;i<Ne;i++)
			for (j=0;j<M;j++)
			if ((s[i][j] > 0.9*C_max) && (post[i][j] <Ne)) str++;
		str=100*str/Ne/M;

//		sec, [frate(end), str(end), sum(firings(:,2)>Ne)/Ni]
		double	ifrate=0;
		for (i=1;i<N_firings;i++)
			if ((firings[i][0] >=0) && (firings[i][1] >=Ne)) ifrate++;
		ifrate = ifrate/Ni;
		std::cout << "sec=" << sec << ", exc. frate=" << frate << ",    exc->exc str=" << str << ",    inh. frate=" << ifrate << ".\n";
		fx = fopen("..//dat.dat","a");
		fprintf(fx, "%d  %2.2f  %2.2f  %2.2f\n", sec, frate, str, ifrate);
		fclose(fx);

		
//		plot(firings(:,1),firings(:,2),'.');
//		axis([0 1000 0 N]); drawnow;
   		fs = fopen("..//spikest.dat","w");
		for (i=1;i<N_firings;i++)
			if (firings[i][0] >=0)
				fprintf(fs, "%d  %d\n", sec*1000+firings[i][0], firings[i][1]);
		fclose(fs);
		remove("..//spikes.dat"); 
		rename( "..//spikest.dat", "..//spikes.dat" );


//		LTP(:,1:D+1)=LTP(:,1001:1001+D);
		for (i=0;i<N;i++)
			for (j=0;j<D+1;j++)
			LTP[i][j]=LTP[i][1000+j];

//	    ind = find(1001-firings(:,1) < D);
		k=N_firings-1;
		while (1000-firings[k][0]<D) k--;
		
//		firings=[-D 0; firings(ind,1)-1000, firings(ind,2)];
		for (i=1;i<N_firings-k;i++)
		{
			firings[i][0]=firings[k+i][0]-1000;
			firings[i][1]=firings[k+i][1];
		}
		N_firings = N_firings-k;

//      sd=0.9*sd;                         % tau = 250 ms
//      s(1:Ne,:)=max(0,min(7, 0.01+s(1:Ne,:)+sd(1:Ne,:)));
		for (i=0;i<Ne;i++)
		for (j=0;j<M;j++)
		{
			sd[i][j]*=0.9;
			s[i][j]+=0.01+sd[i][j];
			if (s[i][j]>C_max) s[i][j]=C_max;
			if (s[i][j]<0) s[i][j]=0;
		}
    
//    if mod(sec,10)==0, 
//        save all; 
//    end;
//      if ( (sec%10==0) & (sec > 0))
//	  {
//		  save_all("..//all.dat");
//
//		  fs = fopen("..//s.dat", "w");
//		  for (i=0; i<Ne; i++)
//			  for (j=0;j<M; j++)
//				fprintf(fs, "%d %3.3f\n", post[i][j], s[i][j]);
//		  fclose(fs);
//	  }
			
//	end;

	}
    save_all("..//all.dat");

//fpoly = fopen("..//polyall.dat","w");
	all_polychronous(); k=N_polychronous;
//	fclose(fpoly);
//	shuffle();
//	all_polychronous(); 
//	std::cout << "ratio (true/shuffled) = " << double(k)/(N_polychronous+1) << "                                                                                     ";


}

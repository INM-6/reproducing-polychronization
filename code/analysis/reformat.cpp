#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "jsoncpp/json/json.h"
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

FILE	*fout;

Json::FastWriter writer;
Json::Value json_data(Json::arrayValue);

void load_custom_data(char fname[30])
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


void save_json_data(char fname[30]){

    fout = fopen(fname,"w");
    // iterate over all presynaptic neurons
    for (int n = 0; n < N; ++n)
    {
        // iterate over all connections
        for( int m = 0; m < M; ++m)
        {
            Json::Value json_connection;
            json_connection["pre"] = n + 1; //TODO can we remove the +1 ?
            json_connection["post"] = post[n][m] + 1;
            json_connection["weight"] = s[n][m];
            if (n<Ne)
            {json_connection["delay"] = m / (M/D) + 1;} // is that correct??
            else
            {json_connection["delay"] = 1;}
            json_data.append(json_connection);
        }

    }

    fprintf(fout, writer.write(json_data).c_str());
    fclose(fout);
}

int main(int argc, char *argv[])
{

  load_custom_data(argv[1]);

  save_json_data(argv[2]);


}


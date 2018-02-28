#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "jsoncpp/json/json.h"

const int Ne = 800;
const int Ni = 200;
const int N  = Ne+Ni;
const int M  = 100;
const int D  =  20;

const double C_max = 10.0;

int post[N][M];
double s[N][M];

short delays[N][D][M];

int I_pre[N][3*M];
double * s_pre[N][3*M];
int D_pre[N][3*M];

int N_pre[N];
int N_post[N];
short delays_length[N][D];

//--------------------------------------------------------------

double a[N], d[N];

const int W = 3; // initial width of polychronous groups
const int min_group_path = 7; // minimal length of a group

int N_polychronous;

double C_rel = 0.95*C_max;
const int polylenmax = N;

int N_postspikes[polylenmax];
int I_postspikes[polylenmax][N];
int J_postspikes[polylenmax][N];
int D_postspikes[polylenmax][N];
int L_postspikes[polylenmax][N];
double C_postspikes[polylenmax][N];

int N_links;
int links[2*W*polylenmax][4];

int group[polylenmax];
int t_fired[polylenmax];
int layer[polylenmax];

int gr3[W], tf3[W];

int I_my_pre[3*M];
int D_my_pre[3*M];
int N_my_pre;

int N_fired;
FILE *fpoly;

Json::FastWriter writer;
Json::Value json_data(Json::arrayValue);

//--------------------------------------------------------------
// The new (default) algorithm to find polychronous groups

const int latency = D; // maximum latency

void polychronous(int nnum)
{
  int i, j, t, p, k;
  int npre[W];
  int dd;
  int t_last, timing;
  int Dmax, L_max;
  int used[W], discard;

  double v[N], u[N], I[N];


  N_my_pre = 0;
  for ( i=0; i<N_pre[nnum]; i++ )
    if ( *s_pre[nnum][i] > C_rel )
    {
      I_my_pre[N_my_pre]=I_pre[nnum][i];
      D_my_pre[N_my_pre]=D_pre[nnum][i];
      N_my_pre++;
    }

  if ( N_my_pre<W )
    return;

  for ( i=0; i<W; i++ )
    npre[i] = i;

  while ( true )
  {
    Dmax=0;
    for ( i=0; i<W; i++ )
      if ( Dmax < D_my_pre[npre[i]] )
	Dmax=D_my_pre[npre[i]];

    for ( i=0; i<W; i++ )
    {
      group[i] = I_my_pre[npre[i]];
      t_fired[i] = Dmax-D_my_pre[npre[i]];
      layer[i] = 1;

      for ( dd=0; dd<D; dd++ )
	for ( j=0; j<delays_length[group[i]][dd]; j++ )
	{
	  p = post[group[i]][delays[group[i]][dd][j]];
	  if ( (s[group[i]][delays[group[i]][dd][j]] > C_rel) & (dd>=D_my_pre[npre[i]]) )
	  {
	    timing = t_fired[i]+dd+1;
	    J_postspikes[timing][N_postspikes[timing]]=group[i];			      // presynaptic
	    D_postspikes[timing][N_postspikes[timing]]=dd;				      // delay
	    C_postspikes[timing][N_postspikes[timing]]=s[group[i]][delays[group[i]][dd][j]];  // syn weight
	    I_postspikes[timing][N_postspikes[timing]++]=p;				      // index of post target
	  }
	}
    }

    for ( i=0; i<N; i++ )
    {
      v[i] = -70;
      u[i] = 0.2*v[i];
      I[i] = 0;
    };

    N_links = 0;
    N_fired = W;
    t_last = D+D+latency+1;
    t = -1;
    while ( (++t<t_last) & (N_fired < polylenmax) )
    {
      for ( p=0; p<N_postspikes[t]; p++ )
      {
	I[I_postspikes[t][p]]+=C_postspikes[t][p];
      }

      for (i=0;i<N;i++)
      {
	v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);
	v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);
	u[i]+=a[i]*(0.2*v[i]-u[i]);
	I[i]=0;
      }

      for ( i=0; i<N; i++ )
	if ( v[i]>=30 )
	{
	  v[i]  = -65;
	  u[i] += d[i];

	  if ( N_fired < polylenmax )
	  {
	    t_fired[N_fired] = t;
	    group[N_fired++] = i;
	    for ( dd=0; dd<D; dd++ )
	      for ( j=0; j<delays_length[i][dd]; j++ )
		if ( (s[i][delays[i][dd][j]] > C_rel) | (i>=Ne) )
		{
		  timing = t+dd+1;
		  J_postspikes[timing][N_postspikes[timing]]=i;				   // presynaptic
		  D_postspikes[timing][N_postspikes[timing]]=dd;			   // delay
		  C_postspikes[timing][N_postspikes[timing]]=s[i][delays[i][dd][j]];	   // syn weight
		  I_postspikes[timing][N_postspikes[timing]++]=post[i][delays[i][dd][j]];  // index of post target
		}
	    if ( t_last < timing+1 )
	    {
	      t_last = timing+1;
	      if ( t_last > polylenmax-D-1 )
		t_last = polylenmax-D-1;
	    }
	  }
	}
    }

    if ( N_fired>2*W )
    {
      N_links=0;
      L_max=0;
      for ( i=W; i<N_fired; i++ )
      {
	layer[i]=0;
	for ( p=t_fired[i]; (p>t_fired[i]-latency) & (p>=0); p-- )
	  for ( j=0; j<N_postspikes[p]; j++ )
	    if ( (I_postspikes[p][j]==group[i]) & (J_postspikes[p][j]<Ne) )
	    {
	      for ( k=0; k<i; k++ )
		if ( (group[k]==J_postspikes[p][j]) & (layer[k]+1>layer[i]) )
		  layer[i]=layer[k]+1;
	      {
		links[N_links][0]=J_postspikes[p][j];
		links[N_links][1]=I_postspikes[p][j];
		links[N_links][2]=D_postspikes[p][j];
		links[N_links++][3]=layer[i];
		if ( L_max < layer[i] )
		  L_max = layer[i];
	      }
	    }
      }

      discard = 0;
      for ( i=0; i<W; i++ )
      {
	used[i]=0;
	for ( j=0; j<N_links; j++ )
	  if ( (links[j][0] == group[i]) & (links[j][1] < Ne) )
	    used[i]++;
	if ( used[i] == 1 )
	  discard = 1;
      }

      if ( (discard == 0) & (L_max >= min_group_path) )
      {
	    for ( i=0; i<W; i++ )
	    {
	      gr3[i]=group[i];
	      tf3[i]=t_fired[i];
	    };

	    N_polychronous++;
        std::cout << "\ni= " << nnum << ", N_polychronous= " << N_polychronous << ", N_fired = " << N_fired << ", L_max = " << L_max << ", T=" << t_fired[N_fired-1];


        // save group in JSON format
        Json::Value json_group;
        json_group["N_fired"] = N_fired; 
        json_group["L_max"] = L_max; 

        Json::Value json_fired(Json::arrayValue);
	    for ( i=0; i<N_fired; i++ )
        {
            Json::Value json_fire;
            json_fire["neuron_id"] = group[i];
            json_fire["t_fired"] = t_fired[i];
            json_fired.append(json_fire);
        }
        json_group["fired"] = json_fired;

        Json::Value json_links(Json::arrayValue);
	    for ( j=0; j<N_links; j++)
        {
            Json::Value json_link;
            json_link["pre"] = links[j][0];
            json_link["post"] = links[j][1];
            json_link["delay"] = links[j][2];
            json_link["layer"] = links[j][3];
            json_links.append(json_link);
        }
        json_group["links"] = json_links;
        json_data.append(json_group);

        //fprintf(fpoly, "N_fired = %d, L_max = %d, ", N_fired, L_max);

	    //for ( i=0; i<N_fired; i++ )
	    //  fprintf(fpoly, "group[%d] = %d, t_fired[%d] = %d, ", i, group[i], i, t_fired[i]);
        //
	    //for ( j=0; j<N_links; j++)
	    //  fprintf(fpoly, "links[%d] = [%d, %d, %d, %d], ", j, links[j][0], links[j][1], links[j][2], links[j][3]);
	    //fprintf(fpoly, "\n");
      }
    }

    for ( dd=Dmax; dd<t_last; dd++ )
      N_postspikes[dd]=0;
    if ( t_last == polylenmax-D )
      for ( dd=t_last; dd<polylenmax; dd++ )
	N_postspikes[dd]=0;

    i=1;
    while ( ++npre[W-i] > N_my_pre-i )
      if ( ++i > W )
	return;
    while ( i>1 )
    {
      npre[W-i+1]=npre[W-i]+1;
      i--;
    }
  }
}

//--------------------------------------------------------------

void all_polychronous(char *argv[])
{
  int i;

  N_polychronous = 0;
  for ( i=0; i<polylenmax; i++ )
    N_postspikes[i] = 0;

  fpoly = fopen(argv[2],"w");
  for ( i=0; i<Ne; i++ )
    polychronous(i);

  //std::cout << writer.write(json_data).c_str() << std::endl;
  fprintf(fpoly, writer.write(json_data).c_str());

  std::cout << "\nN_polychronous=" << N_polychronous << "\n";
  fclose(fpoly);
}

//--------------------------------------------------------------

int main(int argc, char *argv[])
{
  // initialize neuron variables and counters
  for ( int i=0; i<N; i++ )
  {
    N_pre[i] = 0;
    N_post[i] = 0;
    if ( i<Ne )
    {
      a[i]=0.02;
      d[i]=8.0;
    }
    else
    {
      a[i]=0.1;
      d[i]=2.0;
    }
    for ( short j=0; j<D; j++ )
      delays_length[i][j] = 0;
  }

  int source;
  int target;
  double weight;
  short delay;

  // file should contain four tab-separated columns
  // 1st col: gid of presyn. neuron
  // 2nd col: gid of postsyn. neuron
  // 3rd col: weight of synapse
  // 4th col: delay of synapse
  //
  // gids of exc. neurons 1, .., 800
  // gids of inh. neurons 801, .., 1000
  // delays 1, .., 20
  std::ifstream fp_in;
  fp_in.open(argv[1], std::ios::in);

  Json::Value json_in_data(Json::arrayValue);

  fp_in >> json_in_data;

  for (int c = 0; c < json_in_data.size(); ++c)
  {
    
    source = json_in_data[c]["pre"].asInt();
    target = json_in_data[c]["post"].asInt();
    delay = json_in_data[c]["delay"].asInt();
    weight = json_in_data[c]["weight"].asDouble();

    // map to index
    // gids 0, .., 1000
    // delays 0, .., 19
    
    source--;
    target--;
    delay--;

    //std::cout << source << " " << target << " " << delay << " " << weight << std::endl;

    post[source][N_post[source]] = target;
    s[source][N_post[source]] = weight;

    delays[source][delay][delays_length[source][delay]] = N_post[source];

    // consider only exc. sources here
    if ( source < Ne )
    {
      I_pre[target][N_pre[target]] = source;
      s_pre[target][N_pre[target]] = &s[source][N_post[source]];
      D_pre[target][N_pre[target]] = delay;
      N_pre[target]++;
    }

    // increment counters
    N_post[source]++;
    delays_length[source][delay]++;
  }
  fp_in.close();

  all_polychronous(argv);
}

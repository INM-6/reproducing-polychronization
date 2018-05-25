#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "jsoncpp/json/json.h"

const int Ne = 800;
const int Ni = 200;
const int N = Ne + Ni;
const int M = 100;
const int D = 20;

const double C_max = 10.0;

int post[N][M]; // postsynaptic neuron id of the Mth connection from Nth neuron
double s[N][M]; // outgoing synpatic weights ""

short delays[N][D][M]; // mapping neuron ids and connection id to delays with
                       // delay = m / M / D + 1

int I_pre[N][3 *
             M]; // presynaptic neuron id of the Mth connection from Nth neuron
double *s_pre[N][3 * M]; // incoming synaptic weight
int D_pre[N][3 * M];     // incoming synaptic delay

int N_pre[N];              // number of incoming connections
int N_post[N];             // number of outgoing connections
short delays_length[N][D]; // 5

//--------------------------------------------------------------

double a[N], d[N]; // neuron parameter

const int W = 3;              // initial width of polychronous groups
const int min_group_path = 7; // minimal length of a group

int N_polychronous; // number of groups found

double C_rel =
    0.95 * C_max;         // relative synaptic weight considered to be 'strong'
const int polylenmax = N; // maximal timespan of group (ms) is the number of
                          // neurons in the network (1000)
                          // also maximal number of spikes in the group

int N_postspikes[polylenmax];       // number of spikes over time
int I_postspikes[polylenmax][N];    // postesynaptic neuron id of these spikes
int J_postspikes[polylenmax][N];    // presynaptic neuron id of these spikes
int D_postspikes[polylenmax][N];    // delay of these spikes
int L_postspikes[polylenmax][N];    // layer of this connection belonging to these
                                    // spikes
double C_postspikes[polylenmax][N]; // weight of these spikes`

int N_links;                      // number of links ??
int links[2 * W * polylenmax][4]; // links consisting of pre, post, delay, layer

int group[polylenmax];   // neuron id of spiking neuron
int t_fired[polylenmax]; // timepoint of spike
int layer[polylenmax];   // layer of the connection belonging to this spike

int gr3[W], tf3[W]; // unused??

int I_my_pre[3 * M]; // presynaptic neuron ids of the mother neuron with weight
                     // > 9.5
int D_my_pre[3 * M]; // delay between presynaptic neurons and the mother neuron
                     // with weight > 9.5
int N_my_pre; // number of presynaptic neurons of the mother neuron with weight
              // > 9.5

int N_fired; // number of spikes in the group
FILE *fpoly;

Json::FastWriter writer;
Json::Value json_data(Json::arrayValue);

//--------------------------------------------------------------
// The new (default) algorithm to find polychronous groups

const int latency = D; // maximum latency

void polychronous(int nnum) {

  int i, j, t, p, k;
  int npre[W];
  int dd;
  int t_last, timing;
  int Dmax, L_max;
  int used[W], discard;

  double v[N], u[N], I[N];

  N_my_pre = 0;
  for (i = 0; i < N_pre[nnum]; i++)
    if (*s_pre[nnum][i] > C_rel) // find presynaptic neurons of the mother
                                 // neuron with strong weight (> 9.5)
    {
      I_my_pre[N_my_pre] = I_pre[nnum][i];  // presynaptic neuron id of the mother neuron
      D_my_pre[N_my_pre] = D_pre[nnum][i];  // delay from the presynaptic neuron to mother neuron
      N_my_pre++;
    }

  if (N_my_pre < W) // if number of strong connections < 3, there are not enough
                    // anchor neurons therefore discard this mother neuron
    return;

  for (i = 0; i < W; i++) // permutation of anchor neurons ?
    npre[i] = i;

  while (true) {
    Dmax = 0;
    for (i = 0; i < W; i++)
      if (Dmax < D_my_pre[npre[i]])
        Dmax = D_my_pre[npre[i]]; // Dmax is the maximal delay between anchor neurons
                                  // and mother neuron

    for (i = 0; i < W; i++) {                               // loop over all anchor neurons
      group[i] = I_my_pre[npre[i]];                         // add them to the group
      t_fired[i] = Dmax - D_my_pre[npre[i]];                // set theit firing time relative to dmax. The neuron with the longest delay has t_fired = 0
      layer[i] = 1;                                         // the layer of anchor neurons is 1

      for (dd = 0; dd < D; dd++)                                            // loop over all possible delays
        for (j = 0; j < delays_length[group[i]][dd]; j++) {                 // iterate over the connections from current anchor neuron to its postsynaptic neurons

          p = post[group[i]][delays[group[i]][dd][j]];                      // get postsynaptic neuron for the current connection

          if ((s[group[i]][delays[group[i]][dd][j]] > C_rel) &              // check if the synaptic weight is larger than 9.5 
              (dd >= D_my_pre[npre[i]])) {                                  // and the delay of this connection is equal or lager than the delay from the current anchor neuron to the mother neuron
            timing = t_fired[i] + dd + 1;                                   // calculate the arrival time of the PSP

                                                                            // save the timepoints and weights when to stimulate the postsynaptic neurons 
                                                                            // note that the anchor neurons are not stimulated but only their postsynaptic targets
                                                                            // also note that connections having a delay shorter than the delay to the mother neurons are excluded.
                                                                            //
            J_postspikes[timing][N_postspikes[timing]] = group[i];          // set presynaptic neuron id
            D_postspikes[timing][N_postspikes[timing]] = dd;                // delay   
            C_postspikes[timing][N_postspikes[timing]] =                    // weight
                s[group[i]][delays[group[i]][dd][j]]; 
            I_postspikes[timing][N_postspikes[timing]++] = p;               // postsynaptic neuron id
                
          }
        }
    }

    for (i = 0; i < N; i++) {               // reset neurons
      v[i] = -70;
      u[i] = 0.2 * v[i];
      I[i] = 0;
    };

    N_links = 0;
    N_fired = W;                                                    // the anchor neurons have spiked already
    t_last = D + D + latency + 1;                                   // simulate at least until 61 ms
    t = -1;                                                         // current timestep 
    while ((++t < t_last) & (N_fired < polylenmax)) {               // start simulation until t_last or until number of spikes in the group exceeds 1000


      for (p = 0; p < N_postspikes[t]; p++) {                       // get all arriving spikes in this timestep
        I[I_postspikes[t][p]] += C_postspikes[t][p];                // and stimulate the postsynaptic neurons 
      }

      for (i = 0; i < N; i++) {                                         // simulate neuron dynamics
        v[i] += 0.5 * ((0.04 * v[i] + 5) * v[i] + 140 - u[i] + I[i]);
        v[i] += 0.5 * ((0.04 * v[i] + 5) * v[i] + 140 - u[i] + I[i]);
        u[i] += a[i] * (0.2 * v[i] - u[i]);
        I[i] = 0;                                                       // reset input to 0
      }

      for (i = 0; i < N; i++)                                           // loop over all neurons 
        if (v[i] >= 30) {                                               // check if a neuron has emmited a spike
          v[i] = -65;                                                   // reset neuron
          u[i] += d[i];

          if (N_fired < polylenmax) {                                   // check if group is already too large
            t_fired[N_fired] = t;                                       // add timepoint of spike to the group
            group[N_fired++] = i;                                       // and the id of the spiking neuron
            for (dd = 0; dd < D; dd++)                                  
              for (j = 0; j < delays_length[i][dd]; j++)
                if ((s[i][delays[i][dd][j]] > C_rel) | (i >= Ne)) {     // generate spikes and save for later delivery but only for strong and inhibitory connections 
                                                                        // note that N_postspikes can get larger than 1000 (polynmax)
                                                                        // this leads to an overflow of the arrays "J_postspikes", "D_postspikes", "C_postspikes" and "I_postspikes"
                                                                        // resulting of corruption of data

                  // WE ADDED A DEBUGGING LINE HERE
                  if ( N_postspikes[timing] >= polylenmax ) 
                  { 
                    std::cout << "N postspikes larger than polynmax (" << N_postspikes[timing] << ") writing data to next timestep. Data corrupted." << std::endl;
                  }
                  // UNTIL HERE
                  
                  
                  timing = t + dd + 1;                                  // calculate timepoint of arrival of this spike
                  J_postspikes[timing][N_postspikes[timing]] = i;       // presynaptic
                  D_postspikes[timing][N_postspikes[timing]] = dd;      // delay
                  C_postspikes[timing][N_postspikes[timing]] =
                      s[i][delays[i][dd][j]];                           // syn weight
                  I_postspikes[timing][N_postspikes[timing]++] =
                      post[i][delays[i][dd][j]];                        // index of post target
                }


            if (t_last < timing + 1) {                                  // set maximal simulation time to the last arriving spike + 1 ms.
              t_last = timing + 1;                                      // note that sometimes neurons spike up to 4 ms later than the arriving spike
                                                                        // and therefore sometimes the last spike might get lost
                    
              if (t_last > polylenmax - D - 1)                          // use polynmax for its second meaning of maximal timespan of the group
                t_last = polylenmax - D - 1;
            }
          }
        }
    }
                                                                                            // simulation is over. From here on the spiking activity is analized for group formation 

    if (N_fired > 2 * W) {                                                                  // a group must be larger than 6 neurons, why?
      N_links = 0;  
      L_max = 0;
      for (i = W; i < N_fired; i++) {                                                       // iterate over all spikes after the three spikes from the anchor neurons
        layer[i] = 0;                                                               

        for (p = t_fired[i]; (p > t_fired[i] - latency) & (p >= 0); p--)                    // iterate 20 ms back in time from current spike time ( or to spike-time of the first anchor neuron )
          for (j = 0; j < N_postspikes[p]; j++)                                             // iterate over all active connections at this timestep (in presynaptic neuron ids)

            if ((I_postspikes[p][j] == group[i]) & (J_postspikes[p][j] < Ne)) {             // discard inhibitory connections and consider only connection to the currently analyzed neuron

              for (k = 0; k < i; k++)                                                       // iterate over all neurons in the group which spiked earlier than the currently analyzed neuron
                if ((group[k] == J_postspikes[p][j]) &                                      // set layer of the current connection to the max layer of presynaptic neuron + 1
                    (layer[k] + 1 > layer[i]))
                  layer[i] = layer[k] + 1;                                                  // note that this can lead to a inconsistent definition of layer dependent on the neuron id and relative timing
              {
                links[N_links][0] = J_postspikes[p][j];                                     // generate layer with pre, post, delay and layer
                links[N_links][1] = I_postspikes[p][j];
                links[N_links][2] = D_postspikes[p][j];
                links[N_links++][3] = layer[i];
                if (L_max < layer[i])                                                       // update LMAX
                  L_max = layer[i];
              }
            }
      }

      discard = 0;                                                                          // viability test 
      for (i = 0; i < W; i++) {                                                             // count the excitatory connections from each anchor neuron in the group where the anchor neuron is presynaptic
        used[i] = 0;                                                                        // the group is discarded if at least one anchor neuron has only one of those connections. (this is always to the mother neuron)
        for (j = 0; j < N_links; j++)
          if ((links[j][0] == group[i]) & (links[j][1] < Ne))
            used[i]++;
        if (used[i] == 1)
          discard = 1;
      }

      if ((discard == 0) & (L_max >= min_group_path)) {                                     // drop groups with less than 7 layers, why?
        for (i = 0; i < W; i++) {                                                           // unused code. It just sets gr3 and tf3 without using these later
          gr3[i] = group[i];
          tf3[i] = t_fired[i];
        };

        N_polychronous++;                                                                   // group found. increase group counter by 1
        std::cout << "\ni= " << nnum << ", N_polychronous= " << N_polychronous              // print group to display
                  << ", N_fired = " << N_fired << ", L_max = " << L_max
                  << ", T=" << t_fired[N_fired - 1];

        // save group in JSON format                                                        // we added this part. It saves the group in JSON format
        Json::Value json_group;
        json_group["N_fired"] = N_fired;
        json_group["L_max"] = L_max;

        Json::Value json_fired(Json::arrayValue);
        for (int iiii = 0; iiii < N_fired; iiii++) {
          Json::Value json_fire;
          json_fire["neuron_id"] = group[iiii];
          json_fire["t_fired"] = t_fired[iiii];
          json_fired.append(json_fire);
        }
        json_group["fired"] = json_fired;

        Json::Value json_links(Json::arrayValue);
        for (int jjjj = 0; jjjj < N_links; jjjj++) {
          Json::Value json_link;
          json_link["pre"] = links[jjjj][0];
          json_link["post"] = links[jjjj][1];
          json_link["delay"] = links[jjjj][2];
          json_link["layer"] = links[jjjj][3];
          json_links.append(json_link);
        }
        json_group["links"] = json_links;
        json_data.append(json_group);                                                       // here our part ends

        // fprintf(fpoly, "N_fired = %d, L_max = %d, ", N_fired, L_max);                    // we commented out the old data format

        // for ( i=0; i<N_fired; i++ )
        //  fprintf(fpoly, "group[%d] = %d, t_fired[%d] = %d, ", i, group[i], i,
        //  t_fired[i]);
        //
        // for ( j=0; j<N_links; j++)
        //  fprintf(fpoly, "links[%d] = [%d, %d, %d, %d], ", j, links[j][0],
        //  links[j][1], links[j][2], links[j][3]);
        // fprintf(fpoly, "\n");                                                            // until here
      }
    }

    for (dd = Dmax; dd < t_last; dd++)                                                      // reset N_postspikes. Make it ready of the next simulation of the next group
      N_postspikes[dd] = 0;

    if (t_last == polylenmax - D)                                                           // if simulation took 980 ms, reset N_postspikes between 980 and 1000 ms
      for (dd = t_last; dd < polylenmax; dd++)
        N_postspikes[dd] = 0;

    i = 1;
    while (++npre[W - i] > N_my_pre - i)                                                    // iterate over all possible triplets of anchor neurons for current mother neuron
      if (++i > W)
        return;

    while (i > 1) {
      npre[W - i + 1] = npre[W - i] + 1;
      i--;
    }
  }
}

//--------------------------------------------------------------

void all_polychronous(char *argv[]) {                                                       // iterate over all excitatory neurons. run polychronous() for each (mother) neuron
  int i;

  N_polychronous = 0;
  for (i = 0; i < polylenmax; i++)
    N_postspikes[i] = 0;

  for (i = 0; i < Ne; i++)
    polychronous(i);
  fpoly = fopen(argv[2], "w");

  // std::cout << writer.write(json_data).c_str() << std::endl;
  fprintf(fpoly, writer.write(json_data).c_str());

  std::cout << "\nN_polychronous=" << N_polychronous << "\n";
  fclose(fpoly);
}

//--------------------------------------------------------------

int main(int argc, char *argv[]) {                                                              // read in connectivity data set neuron parameter and run all_polychronous()
  // initialize neuron variables and counters
  for (int i = 0; i < N; i++) {
    N_pre[i] = 0;
    N_post[i] = 0;
    if (i < Ne) {
      a[i] = 0.02;
      d[i] = 8.0;
    } else {
      a[i] = 0.1;
      d[i] = 2.0;
    }
    for (short j = 0; j < D; j++)
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

  for (int c = 0; c < json_in_data.size(); ++c) {

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

    // std::cout << source << " " << target << " " << delay << " " << weight <<
    // std::endl;

    post[source][N_post[source]] = target;
    s[source][N_post[source]] = weight;

    delays[source][delay][delays_length[source][delay]] = N_post[source];

    // consider only exc. sources here
    if (source < Ne) {
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


  //This eusures a file is created even when the job fails ( which sometimes happens)
  //mainly needed to keep snakemake happy because it works based on output/input files
   Json::Value json_error;
   json_error["Failed"] = 1;

    Json::Value json_fail(Json::arrayValue);
    json_fail.append(json_error);                                                       // here our part ends


  fpoly = fopen(argv[2], "w");
  fprintf(fpoly, writer.write(json_fail).c_str());
  fclose(fpoly);

  all_polychronous(argv);
}

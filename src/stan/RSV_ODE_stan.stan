functions {
  vector RSV_ODE_stan(real time,
                      vector states,
                      int nAges,
                      int Ss_index,
                      int Es_index,
                      int Is_index,
                      int Sp_index,
                      int Ep_index,
                      int Ip_index,
                      int R_index,
                      real b0,
                      real b1,
                      real phi,
                      real delta,
                      real gammas,
                      real gammap,
                      real nu,
                      vector omega_vect,
                      vector sigma_vect,
                      vector alpha_vect,
                      matrix matrix_mean){

    vector[nAges] Sp;
    vector[nAges] Ep;
    vector[nAges] Ip;
    vector[nAges] Ss;
    vector[nAges] Es;
    vector[nAges] Is;
    vector[nAges] R;
    vector[nAges] N;
    vector[nAges] temp;

    matrix[nAges, nAges] s_ij;
    vector[nAges] lambda;
    vector[nAges] infect_p;
    vector[nAges] infect_s;
    vector[(nAges*7)] deriv;

    // filling in the positions
    for (i in 1:nAges) {
      Sp[i] = states[(Sp_index * nAges) + i];
      Ep[i] = states[(Ep_index * nAges) + i];
      Ip[i] = states[(Ip_index * nAges) + i];
      Ss[i] = states[(Ss_index * nAges) + i];
      Es[i] = states[(Es_index * nAges) + i];
      Is[i] = states[(Is_index * nAges) + i];
      R[i] = states[(R_index * nAges) + i];
      N[i] = Ss[i] + Es[i] + Is[i] + Sp[i] + Ep[i] + Ip[i] + R[i];
      temp[i] = omega_vect[i] * (Is[i] + Ip[i]) / N[i];
    }

    for(i in 1:nAges){
      for(j in 1:nAges){
        s_ij[i, j] = matrix_mean[i, j] * temp[j];
      }
    }

    for(i in 1:nAges){
      lambda[i] = b0 * (1 + b1 * cos(2 * pi() / 365.25 * (time - phi))) * sum(s_ij[i,]);
      infect_p[i] = lambda[i] * sigma_vect[i] * Sp[i];
      infect_s[i] = lambda[i] * sigma_vect[i] * alpha_vect[i] * Ss[i];
      deriv[(Sp_index * nAges) + i] = -infect_p[i];
      deriv[(Ep_index * nAges) + i] = infect_p[i] - delta * Ep[i];
      deriv[(Ip_index * nAges) + i] = delta * Ep[i] - gammap * Ip[i];
      deriv[(Ss_index * nAges) + i] = -infect_s[i] + nu * R[i];
      deriv[(Es_index * nAges) + i] = infect_s[i] - delta * Es[i];
      deriv[(Is_index * nAges) + i] = delta * Es[i] - gammas * Is[i];
      deriv[(R_index * nAges) + i] = gammap * Ip[i] + gammas * Is[i] - nu * R[i];
    }

    return deriv;
  }

  array[] vector cohort_ageing_stan(int n_times,
                                    int n_steps,
                                    real t0,
                                    array[,] real times_array,
                                    array[] int n_times_array,
                                    array[] int cumn_times_array,
                                    int Ss_index,
                                    int Es_index,
                                    int Is_index,
                                    int Sp_index,
                                    int Ep_index,
                                    int Ip_index,
                                    int R_index,
                                    vector init_conds,
                                    int nAges,
                                    real total_population,
                                    real b0,
                                    real b1,
                                    real phi,
                                    real delta,
                                    real gammas,
                                    real gammap,
                                    real nu,
                                    vector omega_vect,
                                    vector sigma_vect,
                                    vector alpha_vect,
                                    matrix matrix_mean,
                                    vector transition_rate,
                                    vector rel_sizes
                                    ){

                                      array[(n_times + 1)] vector[nAges * 7] out;
                                      vector[nAges * 7] next_state = init_conds;
                                      real t0_in = t0;
                                      // initial time
                                      out[1] = init_conds;

                                      // subsequent times
                                      for(i in 1:n_steps){

                                        array[n_times_array[i]] vector[nAges * 7] states = ode_rk45(RSV_ODE_stan,
                                                                                                    next_state,
                                                                                                    t0_in,
                                                                                                    times_array[i],
                                                                                                    nAges, Ss_index, Es_index, Is_index, Sp_index, Ep_index, Ip_index, R_index,
                                                                                                    b0, b1, phi, delta, gammas, gammap, nu, omega_vect, sigma_vect, alpha_vect, matrix_mean);

                                       for(j in 1:n_times_array[i]){
                                          out[(j + cumn_times_array[i] + 1)] = states[j];
                                        }

                                        next_state = states[n_times_array[i]];

                                        next_state[(Sp_index * nAges) + 1] = next_state[(Sp_index * nAges) + 1] + rel_sizes[1] * transition_rate[1] * total_population - transition_rate[1] * next_state[(Sp_index * nAges) + 1];
                                        next_state[(Ep_index * nAges) + 1] = next_state[(Ep_index * nAges) + 1] - transition_rate[1] * next_state[(Ep_index * nAges) + 1];
                                        next_state[(Ip_index * nAges) + 1] = next_state[(Ip_index * nAges) + 1] - transition_rate[1] * next_state[(Ip_index * nAges) + 1];
                                        next_state[(Ss_index * nAges) + 1] = next_state[(Ss_index * nAges) + 1] - transition_rate[1] * next_state[(Ss_index * nAges) + 1];
                                        next_state[(Es_index * nAges) + 1] = next_state[(Es_index * nAges) + 1] - transition_rate[1] * next_state[(Es_index * nAges) + 1];
                                        next_state[(Is_index * nAges) + 1] = next_state[(Is_index * nAges) + 1] - transition_rate[1] * next_state[(Is_index * nAges) + 1];
                                        next_state[(R_index * nAges) + 1] = next_state[(R_index * nAges) + 1] - transition_rate[1] * next_state[(R_index * nAges) + 1];

                                        for(k in 2:nAges){
                                          next_state[(Sp_index * nAges) + k] = next_state[(Sp_index * nAges) + k] + next_state[(Sp_index * nAges) + k - 1] * transition_rate[k - 1] - next_state[(Sp_index * nAges) + k] * transition_rate[k];
                                          next_state[(Ep_index * nAges) + k] = next_state[(Ep_index * nAges) + k] + next_state[(Ep_index * nAges) + k - 1] * transition_rate[k - 1] - next_state[(Ep_index * nAges) + k] * transition_rate[k];
                                          next_state[(Ip_index * nAges) + k] = next_state[(Ip_index * nAges) + k] + next_state[(Ip_index * nAges) + k - 1] * transition_rate[k - 1] - next_state[(Ip_index * nAges) + k] * transition_rate[k];
                                          next_state[(Ss_index * nAges) + k] = next_state[(Ss_index * nAges) + k] + next_state[(Ss_index * nAges) + k - 1] * transition_rate[k - 1] - next_state[(Ss_index * nAges) + k] * transition_rate[k];
                                          next_state[(Es_index * nAges) + k] = next_state[(Es_index * nAges) + k] + next_state[(Es_index * nAges) + k - 1] * transition_rate[k - 1] - next_state[(Es_index * nAges) + k] * transition_rate[k];
                                          next_state[(Is_index * nAges) + k] = next_state[(Is_index * nAges) + k] + next_state[(Is_index * nAges) + k - 1] * transition_rate[k - 1] - next_state[(Is_index * nAges) + k] * transition_rate[k];
                                          next_state[(R_index * nAges) + k] = next_state[(R_index * nAges) + k] + next_state[(R_index * nAges) + k - 1] * transition_rate[k - 1] - next_state[(R_index * nAges) + k] * transition_rate[k];
                                          }

                                        t0_in =  times_array[i, n_times_array[i]];

                                      }

                           return out;
      }

      array[] vector calc_incidence_stan(array[] vector out,
                                         vector times_vec,
                                         int n_times_vec,
                                         int nAges,
                                         int Sp_index,
                                         int Ep_index,
                                         int Ip_index,
                                         int Ss_index,
                                         int Es_index,
                                         int Is_index,
                                         int R_index,
                                         real b0,
                                         real b1,
                                         real phi,
                                         vector omega_vect,
                                         vector sigma_vect,
                                         vector alpha_vect,
                                         matrix matrix_mean
      ){

    vector[nAges] Sp;
    vector[nAges] Ep;
    vector[nAges] Ip;
    vector[nAges] Ss;
    vector[nAges] Es;
    vector[nAges] Is;
    vector[nAges] R;
    vector[nAges] N;
    vector[nAges] temp;
    vector[nAges * 7] states;

    matrix[nAges, nAges] s_ij;
    vector[nAges] lambda;

    array[(n_times_vec)] vector[nAges] inc;


    // filling in the positions

    for(i in 1:n_times_vec){

      states = out[i];

      for (j in 1:nAges) {
        Sp[j] = states[(Sp_index * nAges) + j];
        Ep[j] = states[(Ep_index * nAges) + j];
        Ip[j] = states[(Ip_index * nAges) + j];
        Ss[j] = states[(Ss_index * nAges) + j];
        Es[j] = states[(Es_index * nAges) + j];
        Is[j] = states[(Is_index * nAges) + j];
        R[j] = states[(R_index * nAges) + j];
        N[j] = Ss[j] + Es[j] + Is[j] + Sp[j] + Ep[j] + Ip[j] + R[j];
        temp[j] = omega_vect[j] * (Is[j] + Ip[j]) / N[j];
        }

        for(l in 1:nAges){
          for(j in 1:nAges){
            s_ij[l, j] = matrix_mean[l, j] * temp[j];
            }
            }

            for(j in 1:nAges){
              lambda[j] = b0 * (1 + b1 * cos(2 * pi() / 365.25 * (times_vec[i] - phi))) * sum(s_ij[j,]);
              inc[i][j] = lambda[j] * sigma_vect[j] * (Sp[j] + alpha_vect[j] * Ss[j]);
              }
        }

        return inc;
      }

}

data {
  int n_times;
  int n_steps;
  int n_times_wsteps;
  real t0;
  array[n_steps] int n_times_array;
  array[n_steps] int cumn_times_array;
  array[n_steps, n_times_wsteps] real times_array;
  int Ss_index;
  int Es_index;
  int Is_index;
  int Sp_index;
  int Ep_index;
  int Ip_index;
  int R_index;
  int nAges;
  vector[(nAges * 7)] init_conds;
  real total_population;
  real b0;
  real b1;
  real phi;
  real delta;
  real gammas;
  real gammap;
  real nu;
  vector[nAges] omega_vect;
  vector[nAges] sigma_vect;
  vector[nAges] alpha_vect;
  matrix[nAges, nAges] matrix_mean;
  vector[nAges] transition_rate;
  vector[nAges] rel_sizes;
}

parameters {

 }

 model {

 }

generated quantities{

  array[(n_times + 1)] vector[(nAges * 7)] gq_sim = cohort_ageing_stan(n_times,
                                                                       n_steps,
                                                                       t0,
                                                                       times_array,
                                                                       n_times_array,
                                                                       cumn_times_array,
                                                                       Sp_index,
                                                                       Ep_index,
                                                                       Ip_index,
                                                                       Ss_index,
                                                                       Es_index,
                                                                       Is_index,
                                                                       R_index,
                                                                       init_conds,
                                                                       nAges,
                                                                       total_population,
                                                                       b0,
                                                                       b1,
                                                                       phi,
                                                                       delta,
                                                                       gammas,
                                                                       gammap,
                                                                       nu,
                                                                       omega_vect,
                                                                       sigma_vect,
                                                                       alpha_vect,
                                                                       matrix_mean,
                                                                       transition_rate,
                                                                       rel_sizes);
}

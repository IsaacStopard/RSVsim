#' Vaccine efficacy against hospitalisation conditional on infection
#'
#' Calculates the vaccine efficacy against hospitalisation given infection using the
#' overall vaccine efficacy against hospitalisation and vaccine efficacy against infection.
#' If the number ages in VE_inf are greater than the number of ages in VE_hosp we just select the
#'
#' @param VE_hosp Vaccine efficacy against hospitalisation.
#' @param VE_inf Vaccine efficacy against infection.
#' @return Vaccine efficacy against hospitalisation given infection (\code{VE_hosp_given_inf}).
#' @export
#'
RSVsim_VE_hosp_given_inf <- function(VE_hosp, VE_inf){
  return(1 - (1 - VE_hosp)/(1-VE_inf))
}

#' Probability of hospitalisation given unvaccinated infection
#'
#' Calculates the age-specific probability of hospitalisation given infection, using observed hospitalisation rates per person among the population as a whole
#' and RSVsim to estimate number of infections per person.
#'
#' @param hosp_rate Vector of age-specific RSV hospitalisation rates (per person) in the status quo scenario.
#' @param age_limits_hosp_rate Vector of ages corresponding to the first age of each age category in (1) the hospitalisation rates (\code{hosp_rate}) and (2) the probability of hospitalisation given infection (\code{prob_hosp_given_inf}). These ages must be present in the age_limits used to run the simulations.
#' @param hosp_rate_min_time Minimum simulation time to select when summing the hospitalisations to calculate the hospitalisation rates.
#' @param hosp_rate_max_time Maximum simulation time to select when summing the hospitalisations to calculate the hospitalisation rates.
#' @param sim_status_quo \code{RSVsim_run_model} output used to estimate the number of infections when calculating the probability of hospitalisation given infection.
#' @param VE_hosp_given_inf Matrix of age-specific vaccine efficacy against hospitalisation given infection (\code{VE_hosp_given_inf}). Rows must correspond to the ages in \code{age_limits}, and columns must correspond to the vaccinated states used in \code{RSVsim_run_model}.
#' @param VE_hosp_vacc_states vector of integers indicating the vaccinated states used in \code{RSVsim_run_model}.
#' @return List of probability of hospitalisation given infection  and \code{sim_status_quo} model outputs with hospitalisations.
#' @export
#'
RSVsim_prob_hosp_given_inf <- function(hosp_rate,
                                       age_limits_hosp_rate,
                                       hosp_rate_min_time,
                                       hosp_rate_max_time,
                                       sim_status_quo,
                                       VE_hosp_given_inf,
                                       VE_hosp_vacc_states){

  age_limits <- unique(sim_status_quo$age)

  age_limits <- round(age_limits, digits = 5)

  age_limits_hosp_rate <- round(age_limits_hosp_rate, digits = 5)

  sim_status_quo[,"time"] <- round(sim_status_quo[, "time"], digits = 5)

  if(base::is.unsorted(age_limits_hosp_rate)){
    stop("RSVsim_prob_hosp_given_inf: age_limits_hosp_rate must be in ascending order and hosp_rate values must be in the same position as the corresponding age")
  }

  if(!all(age_limits_hosp_rate %in% age_limits)){
    stop("RSVsim_prob_hosp_given_inf: not all age_limits_hosp_rate are present in the sim_status_quo")
  }

  if(nrow(VE_hosp_given_inf) != length(age_limits) || ncol(VE_hosp_given_inf) != length(VE_hosp_vacc_states)){
    stop("RSVsim_prob_hosp_given_inf: VE_hosp_given_inf must be a matrix with the rows corresponding to the ages in age_limits
         and the columns corresponding to the VE_hosp_vacc_states")
  }

  age_limits_index <- data.table::data.table(age_index = base::findInterval(age_limits, age_limits_hosp_rate),
                                             age = age_limits)
  age_limits_index[, age_hosp := age_limits_hosp_rate[age_index]]

  sim_status_quo <- data.table::as.data.table(sim_status_quo)
  sim_status_quo[age_limits_index, on = .(age), age_hosp := i.age_hosp]

  # total population by ages in the hospitalisation rate data
  pop <- sim_status_quo[time == 0,
                        .(tot_pop = sum(Sp + Ss + Ep + Es + Ip + Is + R)),
                        by = age_hosp]

  # VE by ages and vacc_state
  VE_hosp_df <- data.table::as.data.table(VE_hosp_given_inf)
  data.table::setnames(VE_hosp_df, as.character(VE_hosp_vacc_states))
  VE_hosp_df[ , age := age_limits]
  VE_hosp_df <- data.table::melt(VE_hosp_df,
                                 id.vars = "age",
                                 variable.name = "vacc_state",
                                 value.name = "VE_hosp_given_inf")
  VE_hosp_df[, vacc_state := as.numeric(as.character(vacc_state))]

  sim_status_quo <- VE_hosp_df[sim_status_quo, on = .(vacc_state, age)]
  sim_status_quo[, Incidence_VE := Incidence * (1 - data.table::fcoalesce(VE_hosp_given_inf, 0))]

  inc_per_person_df <- sim_status_quo[time >= hosp_rate_min_time & time <= hosp_rate_max_time,
                                      .(tot_inc_VE = sum(Incidence_VE)),
                                      by = age_hosp][pop, on = .(age_hosp)]
  inc_per_person_df[order(age_hosp), ':='(hosp_rate = hosp_rate,
                                          inc_per_person = tot_inc_VE / tot_pop)][,prob_hosp_given_inf_uv := hosp_rate / inc_per_person]
  inc_per_person_df <- inc_per_person_df[, .(age_hosp, hosp_rate, prob_hosp_given_inf_uv)]

  if(any(inc_per_person_df[,"prob_hosp_given_inf_uv"] > 1)){
    stop("RSVsim_hospitalisations: the predicted incidence per person is less than the hospitalisations per person")
  }

  sim_status_quo[inc_per_person_df, on = .(age_hosp), prob_hosp_given_inf_uv := i.prob_hosp_given_inf_uv]
  sim_status_quo[VE_hosp_df, on = .(vacc_state, age), VE_hosp_given_inf := i.VE_hosp_given_inf]
  sim_status_quo[, ":="(Hospitalisations = Incidence * (1 - data.table::fcoalesce(VE_hosp_given_inf, 0)) * prob_hosp_given_inf_uv,
                        prob_hosp_given_inf_uv = NULL,
                        VE_hosp_given_inf = NULL,
                        Incidence_VE = NULL,
                        age_hosp = NULL)]

  return(list("prob_hosp_given_inf_uv" = as.vector(unlist(inc_per_person_df[, .(prob_hosp_given_inf_uv)])),
              "sim_status_quo" = as.data.frame(sim_status_quo))
         )
}

#' Calculate number of hospitalisations for \code{RSVsim_run_model} output
#'
#' Function to estimate the number of hospitalisations from the modelled incidence and probability of hospitalisation given unvaccinated infection.
#'
#' @inheritParams RSVsim_prob_hosp_given_inf
#' @param prob_hosp_given_inf Age-specific probability of hospitalisation given unvaccinated infection. Ages must correspond to those in \code{age_limits_hosp_rate}.
#' @param sim_scenario \code{RSVsim_run_model} output. Hospitalisations are calculated for this scenario using the estimated probability of hospitalisation given infection.
#' @return sim_scenario dataframe with hospitalisations column
RSVsim_hospitalisations <- function(prob_hosp_given_inf,
                                    age_limits_hosp_rate,
                                    sim_scenario,
                                    VE_hosp_given_inf,
                                    VE_hosp_vacc_states){

  age_limits <- unique(sim_scenario$age)

  age_limits <- round(age_limits, digits = 5)

  age_limits_hosp_rate <- round(age_limits_hosp_rate, digits = 5)

  sim_scenario[,"time"] <- round(sim_scenario[, "time"], digits = 5)

  if(base::is.unsorted(age_limits_hosp_rate)){
    stop("RSVsim_hospitalisations: age_limits_hosp_rate must be in ascending order and hosp_rate values must be in the same position as the corresponding age")
  }

  if(!all(age_limits_hosp_rate %in% age_limits)){
    stop("RSVsim_hospitalisations: not all age_limits_hosp_rate are present in sim_scenario")
  }

  if(nrow(VE_hosp_given_inf) != length(age_limits) || ncol(VE_hosp_given_inf) != length(VE_hosp_vacc_states)){
    stop("RSVsim_hospitalisations: VE_hosp_given_inf must be a matrix with the rows corresponding to the ages in age_limits
         and the columns corresponding to the VE_hosp_vacc_states")
  }

  if(length(prob_hosp_given_inf) != length(age_limits_hosp_rate)){
    stop("RSVsim_hospitalisations: prob_hosp_given_inf must be the same length as age_limits_hosp_rate")
  }

  age_limits_index <- data.table::data.table(age_index = base::findInterval(age_limits, age_limits_hosp_rate),
                                             age = age_limits)
  age_limits_index[, age_hosp := age_limits_hosp_rate[age_index]]

  sim_scenario <- data.table::as.data.table(sim_scenario)
  sim_scenario[age_limits_index, on = .(age), age_hosp := i.age_hosp]

  inc_per_person_df <- data.table::data.table(prob_hosp_given_inf_uv = prob_hosp_given_inf,
                                              age_hosp = age_limits_hosp_rate)

  VE_hosp_df <- data.table::as.data.table(VE_hosp_given_inf)
  data.table::setnames(VE_hosp_df, as.character(VE_hosp_vacc_states))
  VE_hosp_df[ , age := age_limits]
  VE_hosp_df <- data.table::melt(VE_hosp_df,
                                 id.vars = "age",
                                 variable.name = "vacc_state",
                                 value.name = "VE_hosp_given_inf")
  VE_hosp_df[, vacc_state := as.numeric(as.character(vacc_state))]

  sim_scenario[inc_per_person_df, on = .(age_hosp), prob_hosp_given_inf_uv := i.prob_hosp_given_inf_uv]
  sim_scenario[VE_hosp_df, on = .(vacc_state, age), VE_hosp_given_inf := i.VE_hosp_given_inf]
  sim_scenario[, ":="(Hospitalisations = Incidence * (1 - data.table::fcoalesce(VE_hosp_given_inf, 0)) * prob_hosp_given_inf_uv,
                      prob_hosp_given_inf_uv = NULL,
                      VE_hosp_given_inf = NULL,
                      age_hosp = NULL)]

  return(as.data.frame(sim_scenario))
}



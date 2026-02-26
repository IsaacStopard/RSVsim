#' Runs the RSV model using odin and dust
#'
#' Function to run the transmission model with cohort aging. Prevalence, incidence between the given time-steps and the incidence rate per day are also calculated.
#'
#' @param parameters List of parameters from \code{RSVsim_parameters} function.
#' @param times Simulation times. Default: 0 to 365.25 days with intervals of 0.25 days. The minimum interval size is 1E-5.
#' @param cohort_step_size Time steps to run the model over before adjusting the ages of all cohorts. Default: 1 month. If \code{is.numeric(cohort_step_size) == FALSE} then cohort ageing is not applied. Can have a maximum of 7 decimal places.
#' @param warm_up Length of time-points to exclude before calculating the likelihood. Default: \code{NULL}.
#' @return Simulation output (dataframe). In the dataframe, age refers to the lowest age in the age group.
#' @export
RSVsim_run_model <- function(parameters,
                             times = seq(0, 365.25*1, 0.25),
                             cohort_step_size = 1/12 * 365.25,
                             warm_up = NULL
                             ){

  ##########################################################

  max_t <- max(times)

  # runs the model with cohort ageing
  RSV_dust <- dust2::dust_system_create(generator = RSV_ODE,
                                        pars = parameters,
                                        n_particles = 1,
                                        n_groups = 1,
                                        time = 0,
                                        deterministic = TRUE
                                        ) |> invisible()


  dust2::dust_system_set_state_initial(RSV_dust)

  # do not change this
  # the order is determined by the order in the odin RSV_ODE.R file
  states <- dust2::dust_unpack_index(RSV_dust)

  states_order <- names(states)

  # creating the infection state, vaccination state and age corresponding to each position in dust
  states <- as.character(rep(states_order, lengths(states)))

  n_states <- length(states_order)

  # important for maintaining the same order as RSV_ODE
  ages <- rep(rep(parameters$age_limits, parameters$nVaccStates), n_states) |> base::round(digits = 5)
  vacc_states <- rep(rep(1:parameters$nVaccStates, each = parameters$nAges), n_states) |> base::round(digits = 5)
  age_chr = rep(rep(parameters$age_chr, parameters$nVaccStates), n_states)

  # no cohort aging in these variables
  output_variable_names <- states_order[!(states_order %in% c("Sp", "Ep", "Ip", "Ss", "Es", "Is", "R"))]

  dust_df <- data.frame("state" = as.character(states),
                        "age" = ages,
                        "age_chr" = age_chr,
                        "vacc_state" = vacc_states,
                        "index" = 1:length(states),
                        stringsAsFactors = FALSE
                        )

  if(is.numeric(cohort_step_size)){

    size_cohorts <- parameters$size_cohorts

    rel_sizes <- parameters$rel_sizes

    transition_rate <- cohort_step_size/size_cohorts

    # running for the time of the smallest cohort
    if(base::round(max(diff(times)), digits = 5) >= base::round(cohort_step_size, digits = 5)){
      stop("The maximum time difference is greater than or equal to the cohort step size")
    }

    if(base::round(min(size_cohorts), digits = 5) < base::round(cohort_step_size, digits = 5)){
      stop("The smallest cohort age size is smaller than the cohort_step_size: increase the differences in age limits or decrease the cohort_step_size")
    }

    # calculate indexing required for cohort ageing
    dust_df <- dust_df |>
      dplyr::arrange(index) |>
      dplyr::group_by(state, vacc_state) |>
      dplyr::mutate(
        cohort_ageing = dplyr::if_else(state %in% output_variable_names, 0, 1),
        lag_index = dplyr::if_else(cohort_ageing == 0, index, dplyr::lag(index, 1))) |>
      dplyr::ungroup() |>
      as.data.frame()

    n_dust <- nrow(dust_df)

    transition_rate_all <- transition_rate[match(dust_df[,"age"], round(parameters[["age_limits"]], digits = 5))]

    birth_index <- which(is.na(dust_df[,"lag_index"]) & dust_df[,"age"] == 0 & dust_df[,"state"] == "Sp" & dust_df[,"vacc_state"] == 1)
    n_births <- rel_sizes[1] * transition_rate[1] * parameters$total_population
    births_zero_index <- (birth_zero_index <- which(is.na(dust_df[,"lag_index"])))[birth_zero_index != birth_index]

    n_steps <- ceiling(max_t / cohort_step_size)

    # running the model with cohort aging (run for a single cohort, move cohort, change initial states, repeat)
    out_list <- vector(mode = "list", length = n_steps)

    times_all <- sort(unique(base::round(c(times, 1:n_steps * cohort_step_size), digits = 5)))

    times_in <- lapply(1:n_steps, FUN = function(i){
      c(times_all[times_all >= base::round(((i - 1) * cohort_step_size), digits = 5) & times_all < base::round((i * cohort_step_size), digits = 5)])
    })

    if(max(times_in[[n_steps]]) < max_t){
      times_in[[n_steps]] <- c(times_in[[n_steps]], max_t)
    }

    for(i in 1:n_steps){

      out <- dust2::dust_system_simulate(RSV_dust, times = times_in[[i]]) |> t()

      # if(ncol(out) != n_dust || nrow(out) != length(times_in[[i]])){
      #    stop("dust size error")
      #   }

      next_state <- out[base::length(times_in[[i]]),]

      out <- base::cbind(out, times_in[[i]])

      out_list[[i]] <- out

      # calculating the initial states with cohort ageing in the S, E, I, R compartments only
      # assumes that ages are sequential
      transition_ct <- next_state * transition_rate_all
      lag_transition_ct <- transition_ct[dust_df[,"lag_index"]]

      # filling in births to keep the population size constant
      lag_transition_ct[birth_index] <- n_births
      lag_transition_ct[births_zero_index] <- 0 # everyone is born unvaccinated

      next_value <- next_state - transition_ct + lag_transition_ct

      dust2::dust_system_set_state(sys = RSV_dust, state = next_value)

      # if(all(dust2::dust_system_state(RSV_dust) != next_state)){
      #   stop("RSVsim_run_model: states have not been updated")
      # }
    }

    out_checkout <- do.call(rbind, out_list) |> as.data.frame()

    rm(list = c("out_list"))

    } else{

      times <- base::round(times, digits = 5)

      out_checkout <- dust2::dust_system_simulate(RSV_dust, times = times) |> t() |> as.data.frame()

      out_checkout[,"time"] <- times

    }

  colnames(out_checkout) <- c(base::as.character(dust_df[,"index"]), "time")

  #if(any(abs(unique(round(rowSums(out_checkout[, which(!(states %in% output_variable_names))]), digits = 5)) - parameters$total_population) > 1E-1)){
  #  stop("RSVsim_run_model: matrix population not correct")
  #}

  # warm up
  if(!is.null(warm_up)){
    out_checkout <- out_checkout |> dplyr::filter(time >= times[min(which(times >= warm_up)) - 1])
  }

  # incidence calculation
  data.table::setDTthreads(1)
  data.table::setDT(out_checkout)
  data.table::setDT(dust_df)

  # pivot longer
  out_checkout <- data.table::melt(out_checkout,
                                   id.vars = "time",
                                   variable.name = "index",
                                   value.name = "value")

  # Combine with dust_df
  out_checkout[, index := as.integer(as.character(index))]
  out_checkout <- out_checkout[dust_df[, c("state", "age", "age_chr", "vacc_state", "index")], on = "index", nomatch = NA]

  out_checkout[, `:=`(age = round(age, 5), vacc_state = as.integer(vacc_state), time = round(time, 5))]

  # pivot wider
  out_checkout <- data.table::dcast(out_checkout, age + age_chr + vacc_state + time ~ state, value.var = "value", drop = TRUE)

  # calculating the population sizes by time
  pop_check <- base::as.data.frame(out_checkout[, .(total = sum(rowSums(.SD))), by = time, .SDcols = c("Sp", "Ep", "Ip", "Ss", "Es", "Is", "R")])

  # arrange
  data.table::setorder(out_checkout, vacc_state, age, time)
  data.table::setcolorder(out_checkout, c("age", "age_chr", "vacc_state", "time", states_order))

  out_checkout <- as.data.frame(out_checkout) |>
                                dplyr::group_by(age_chr, vacc_state) |>
                                dplyr::mutate(Incidence = tidyr::replace_na(Incidence - dplyr::lag(Incidence, 1), 0),
                                              DetIncidence = tidyr::replace_na(DetIncidence - dplyr::lag(DetIncidence, 1), 0),
                                              doses = tidyr::replace_na(doses - dplyr::lag(doses, 1), 0)) |>
                                dplyr::ungroup()

  # warm up
  if(!is.null(warm_up)){
      out_checkout <- out_checkout |> dplyr::filter(time >= warm_up) |> dplyr::mutate(time = time - warm_up)
  }

  # checking the total population is correct
  if(any(abs(pop_check[,"total"] - parameters$total_population) > 1E-1)){
    stop("RSVsim_run_model: population does not sum to the correct number")
  }

  return(out_checkout |> as.data.frame())
}


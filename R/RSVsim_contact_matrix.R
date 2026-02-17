#' Create contact matrix for RSVsim
#'
#' Uses the \code{socialmixr} package to access the POLYMOD data and calculate per-person mean contact matrix.
#' If age groupings smaller than 1-year are required, the function calculates the total contacts in the
#' population with the nearest integer age groupings (required for \code{socialmixr} package) and splits the total contacts
#' into smaller age groupings by dividing them uniformly.
#'
#' @param country Country for use in the \code{contact_matrix} function in the \code{socialmixr} package. Can be given as country name or 2 digit ISO code. United Kingdom default.
#' @param age_limits Lower limits of the age groups to run the simulation. Ages must be in years. Note the smallest age limits allowed are 1E-8.
#' @return List including the contact matrix (\code{matrix_mean} is the mean per person, \code{matrix_contacts} is the total contacts),
#' the age limits, the age distribution used to calculate the total contacts and the maximum age in the contact data.
#' @export
RSVsim_contact_matrix <- function(country = "United Kingdom",
                                  age_limits = c(seq(0,5,2/12), seq(10,70,5))){

  if(!is.numeric(age_limits)){
    stop("age_limits must be numeric")
  }

  if(!country %in% c("Italy", "Germany", "Luxembourg", "Netherlands", "Poland", "United Kingdom", "Finland", "Belgium")){
    stop("for the polymod data the country must be one of the following characters:
         Italy, Germany, Luxembourg, Netherlands, Poland, United Kingdom, Finland or Belgium")
  }


  if(max(age_limits) > 75){
    warning(paste0("polymod age groupings only go up to 75. Age limits above this have therefore been omitted."))
    age_limits <- age_limits[age_limits <= 75]
  }

  if(is.unsorted(age_limits)){
     age_limits <-  sort(age_limits)
     warning("age_limits has been sorted")
  }

  nAges <- length(age_limits)

  age_limits <- base::round(age_limits, digits = 10)

  if(all(age_limits %% 1 == 0)){
    age_limits_default <- age_limits

  } else{
    age_limits_default <- age_limits[age_limits %% 1 == 0]
    warning("Not all age_limits are integers.
            The contact matrix was calculated for the age groups given by the integers
            and divided uniformly within the smaller age groups.
            A warning from socialmixr about linearly estimating the integer age groups within the 5-year age bands will also appear.")

  }

  polymod_cm <- socialmixr::contact_matrix(survey = socialmixr::polymod,
                                           countries = country,
                                           age_limits = age_limits_default,
                                           symmetric = TRUE,
                                           per_capita = FALSE,
                                           counts = FALSE,
                                           missing_contact_age = "remove",
                                           missing_participant_age = "remove",
                                           return_demography = TRUE)

  D <- polymod_cm$demography$population
  M <- polymod_cm$matrix # total number of contacts matrix
  total_contacts_M <- D * M

  if(!isSymmetric(total_contacts_M, check.attributes = FALSE)){
    stop("total contacts matrix is not symmetric")
  }


  max_age <- socialmixr::polymod$participants$part_age |> subset(country == country) |> na.omit() |> max()

  out <- list("age_limits" = age_limits,
              "max_age" = max_age)

  # labels for the ages
  age_chr <- c()

  for(i in 1:(nAges - 1)){
    age_chr <- c(age_chr, c(paste0("[",round(age_limits[i], digits = 2),",", round(age_limits[i+1], digits = 2),")")))
  }

  age_chr <- c(age_chr, paste0("[",round(age_limits[nAges], digits = 2),",", round(max_age, digits = 2),")"))

  if(any(age_limits %% 1 != 0)){
    matrix_out <- matrix(NA, ncol = nAges, nrow = nAges)
    D_out <- rep(NA, length = nAges)

    max_age_default <- max(age_limits_default)

    for(i in 1:nAges){

      m_i <- findInterval(age_limits[i], age_limits_default)

      default_size_i <- if(age_limits_default[m_i] <  max_age_default){
        age_limits_default[m_i + 1] - age_limits_default[m_i]
      } else{
        max_age - age_limits_default[m_i]
      }

      actual_size_i <- if(age_limits[i] < max(age_limits)){
        age_limits[i + 1] - age_limits[i]
      } else{
        max_age - age_limits[i]
      }

      D_out[i] <- D[m_i] / (default_size_i / actual_size_i)

      for(j in 1:nAges){

        # default interval
        # i - rows
        # j - columns

        m_j <- findInterval(age_limits[j], age_limits_default)

        contacts <- total_contacts_M[m_i, m_j]

        # size interval of default ages
        default_size_j <- if(age_limits_default[m_j] <  max_age_default){
          age_limits_default[m_j + 1] - age_limits_default[m_j]
          } else{
          max_age - age_limits_default[m_j]
          }

        actual_size_j <- if(age_limits[j] < max(age_limits)){
        age_limits[j + 1] - age_limits[j]
      } else{
          max_age - age_limits[j]
      }

      div <- (default_size_j / actual_size_j) * (default_size_i / actual_size_i)

      matrix_out[i, j] <- contacts / div
    }
    }


  rownames(matrix_out) <- colnames(matrix_out) <- age_chr
  # rownames(M) <- colnames(M) <- age_limits_default

  per_person_matrix_out <- matrix_out / D_out

  if(abs(sum(total_contacts_M) - sum(matrix_out)) > 1E-5){
    stop("Total contacts do not add up correctly")
  }

  if(abs(sum(D_out) - sum(D)) > 1E-5){
    stop("Total population does not add up correctly")
  }


  if(!isSymmetric(matrix_out)){
    stop("Adjusted matrix should be symmetric")
  }

  out <- append(out, list("matrix_contacts" = matrix_out,
                          "matrix_per_person" = per_person_matrix_out,
                          "population" = D_out,
                          "age_distribution" = D_out/sum(D_out))
              )
  } else{
    out <- append(out, list("matrix_per_person" = M,
                          "matrix_contacts" = total_contacts_M,
                          "population" = D,
                          "age_distribution" = D/sum(D)))
  }

  if(abs(sum(out$age_distribution) - 1) > 1E-5){
    stop("age_distribution does not sum to 1")
  }



  # age differences in days
  size_cohorts <- c(diff(age_limits * 365.25), max_age*365.25 - age_limits[length(age_limits)]*365.25)

  rel_sizes <- size_cohorts/sum(size_cohorts)

  out <- append(out, list("age_chr" = age_chr,
                          "size_cohorts" = size_cohorts,
                          "rel_sizes" = rel_sizes))

  return(out)
}

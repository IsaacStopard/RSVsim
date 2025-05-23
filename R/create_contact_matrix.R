#' Create contact matrix for RSVsim
#'
#' Uses the \code{socialmixr} package to access the POLYMOD data and calculate per-person mean contact matrix.
#' If age groupings smaller than 1-year are required, the function calculates the total contacts in the
#' population with the nearest integer age groupings (required for \code{socialmixr} package) and splits the total contacts
#' into smaller age groupings by dividing them uniformly.
#'
#' @param country Country for use in the \code{contact_matrix} function in the \code{socialmixr} package. Can be given as country name or 2 digit ISO code. United Kingdom default.
#' @param age.limits Lower limits of the age groups to run the simulation. Ages must be in years.
#' @return Contact matrix.
#' @export
create_contact_matrix <- function(country = "United Kingdom",
                                  age.limits = c(seq(0,5,1/12), seq(10,70,5))){

  if(!is.numeric(age.limits)){
    stop("age.limits must be numeric")
  }

  if(!country %in% c("Italy", "Germany", "Luxembourg", "Netherlands", "Poland", "United Kingdom", "Finland", "Belgium")){
    stop("for the polymod data the country must be one of the following characters:
         Italy, Germany, Luxembourg, Netherlands, Poland, United Kingdom, Finland or Belgium")
  }


  if(max(age.limits) > 75){
    warning(paste0("polymod age groupings only go up to 75. Age limits above this have therefore been omitted."))
    age.limits <- age.limits[age.limits <= 75]
  }

  if(is.unsorted(age.limits)){
     age.limits <-  sort(age.limits)
     warning("age.limits has been sorted")
  }

  if(all(age.limits %% 1 == 0)){
    age.limits.default <- age.limits

  } else{
    age.limits.default <- age.limits[age.limits %% 1 == 0]
    warning("Not all age.limits are integers.
            The contact matrix was calculated for the age groups given by the integers
            and divided uniformly within the smaller age groups.
            A warning from socialmixr about linearly estimating the integer age groups within the 5-year age bands will also appear.")

  }

  polymod_cm <- socialmixr::contact_matrix(survey = socialmixr::polymod,
                                           countries = country,
                                           age.limits = age.limits.default,
                                           symmetric = TRUE,
                                           per.capita = FALSE,
                                           counts = FALSE,
                                           missing.contact.age = "remove",
                                           missing.participant.age = "remove",
                                           return.demography = TRUE)

  D <- polymod_cm$demography$population
  M <- polymod_cm$matrix # total number of contacts matrix
  total_contacts_M <- D * M

  if(!isSymmetric(total_contacts_M, check.attributes = FALSE)){
    stop("total contacts matrix is not symmetric")
  }

  out <- list()

  if(any(age.limits %% 1 != 0)){

    nAges <- length(age.limits)
    matrix_out <- matrix(NA, ncol = nAges, nrow = nAges)
    D_out <- rep(NA, length = nAges)

    max_age <- socialmixr::polymod$participants$part_age |> subset(country == country) |> na.omit() |> max()

    max_age_default <- max(age.limits.default)

    for(i in 1:nAges){

      m_i <- findInterval(age.limits[i], age.limits.default)

      default_size_i <- if(age.limits.default[m_i] <  max_age_default){
        age.limits.default[m_i + 1] - age.limits.default[m_i]
      } else{
        max_age - age.limits.default[m_i]
      }

      actual_size_i <- if(age.limits[i] < max(age.limits)){
        age.limits[i + 1] - age.limits[i]
      } else{
        max_age - age.limits[i]
      }

      D_out[i] <- D[m_i] / (default_size_i / actual_size_i)

      for(j in 1:nAges){

        # default interval
        # i - rows
        # j - columns

        m_j <- findInterval(age.limits[j], age.limits.default)

        contacts <- total_contacts_M[m_i, m_j]

        # size interval of default ages
        default_size_j <- if(age.limits.default[m_j] <  max_age_default){
          age.limits.default[m_j + 1] - age.limits.default[m_j]
          } else{
          max_age - age.limits.default[m_j]
          }

        actual_size_j <- if(age.limits[j] < max(age.limits)){
        age.limits[j + 1] - age.limits[j]
      } else{
          max_age - age.limits[j]
      }

      div <- (default_size_j / actual_size_j) * (default_size_i / actual_size_i)

      matrix_out[i, j] <- contacts / div
    }
    }

    rownames(matrix_out) <- colnames(matrix_out) <- age.limits
  rownames(M) <- colnames(M) <- age.limits.default

  mean_matrix_out <- matrix_out / D_out

  if(abs(sum(total_contacts_M) - sum(matrix_out)) > 1E-5){
    stop("Total contacts do not add up correctly")
  }

  if(abs(sum(D_out) - sum(D)) > 1E-5){
    stop("Total population does not add up correctly")
  }


  if(!isSymmetric(matrix_out)){
    stop("Adjusted matrix should be symmetric")
  }

  out <- append(out, list("adjusted_matrix_contacts" = matrix_out,
                          "adjusted_matrix_mean" = mean_matrix_out,
                          "adjusted_population" = D_out)
              )
  }

  out <- append(out, list("default_matrix_mean" = M,
                          "default_matrix_contacts" = total_contacts_M,
                          "default_population" = D,
                          "max_age" = max_age))

  return(out)
}

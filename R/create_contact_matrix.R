#' Calculate the contact matrix with unequal sized age groupings
#' Uses the socialmixr package to access the POLYMOD data and calculate a per-capita, symmetric contact matrix for the nearest integer age groupings
#' @param country country for use in the \code{contact_matrix} function in the socialmixr package. Can be given as country name or 2 digit ISO code. United Kingdom default.
#' @param age.limits lower limits of the age groups to run the simulation
#' @return contact matrix
#' @export
create_contact_matrix <- function(country = "United Kingdom",
                                  age.limits = c(seq(0,5,1/12), seq(10,70,5))){

  if(is.unsorted(age.limits)){
     age.limits <-  sort(age.limits)
     warning("age.limits has been sorted")
  }

  base <- 5
  max <- 70
  age.limits.default <- seq(0, max, base) # to align with Mossong age groups

  floor_base <- function(x, base){
    base * floor(x / base)
  }

  polymod_cm <- socialmixr::contact_matrix(survey = socialmixr::polymod,
                                           countries = country,
                                           age.limits = age.limits.default,
                                           symmetric = FALSE,
                                           per.capita = FALSE,
                                           counts = FALSE,
                                           missing.contact.age = "remove",
                                           missing.participant.age = "remove",
                                           return.demography = TRUE)

  D <- polymod_cm$demography$population

  M <- polymod_cm$matrix # contact matrix (average number of contacts)

  # making the matrix symmetric
  sM <- 0.5 / D * (D * M + D * t(M))

  nAges <- length(age.limits)

  matrix_out <- matrix(NA, ncol = nAges, nrow = nAges)

  # assumes the same number of mean contacts for the closest age grouping
  max_age <- socialmixr::polymod$participants$part_age |> subset(country == country) |> na.omit() |> max()

  for(i in 1:nAges){
    for(j in 1:nAges){

      m_i <- if(age.limits[i] > max){length(age.limits.default)} else{which(age.limits.default == floor_base(age.limits[i], base))}
      m_j <- if(age.limits[j] > max){length(age.limits.default)} else{which(age.limits.default == floor_base(age.limits[j], base))}

      contacts <- sM[m_i, m_j]

      default_size_j <- if(age.limits.default[m_j] < max(age.limits.default)){
        age.limits.default[m_j + 1] - age.limits.default[m_j]
      } else{
          max_age - age.limits.default[m_j]
        }

      actual_size_j <- if(age.limits[j] < max(age.limits)){
        age.limits[j + 1] - age.limits[j]
      } else{
          max_age - age.limits[j]
        }

      div_j <- default_size_j / actual_size_j

      matrix_out[i, j] <- contacts / div_j
    }
  }

  rownames(matrix_out) <- colnames(matrix_out) <- age.limits
  rownames(sM) <- colnames(sM) <- age.limits.default

  return(list("default_matrix" = sM,
              "adjusted_matrix" = matrix_out))
}

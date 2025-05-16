#' Calculate the contact matrix with unequal sized age groupings
#' Uses the socialmixr package to access the POLYMOD data and calculate a per-capita, symmetric contact matrix for the nearest integer age groupings
#' @param country country for use in the \code{contact_matrix} function in the socialmixr package. Can be given as country name or 2 digit ISO code. United Kingdom default.
#' @param age.limits lower limits of the age groups to run the simulation
#' @return contact matrix
#' @export
create_contact_matrix <- function(country = "United Kingdom",
                                  age.limits = c(seq(0,5,1/12), seq(10,75,5))){

  if(is.unsorted(age.limits)){
     age.limits <-  sort(age.limits)
     warning("age.limits has been sorted")
  }

  age.limits.default <- age.limits |> floor() |> as.integer() |> unique()

  mixing <- socialmixr::contact_matrix(survey = socialmixr::polymod,
                                       countries = country,
                                       age.limits = age.limits.default,
                                       symmetric = TRUE,
                                       per.capita = TRUE,
                                       missing.contact.age = "remove",
                                       missing.participant.age = "remove")

  nAges <- length(age.limits)

  matrix_out <- matrix(NA, ncol = nAges, nrow = nAges)

  # assumes the same number of per-capita contacts for the closest age grouping
  for(i in 1:nAges){
    for(j in 1:nAges){

      m_i <- which(age.limits.default == floor(age.limits[i]))
      m_j <- which(age.limits.default == floor(age.limits[j]))

      contacts <- mixing$matrix.per.capita[m_i, m_j]

      max_age <- socialmixr::polymod$participants$part_age |> subset(country == country) |> na.omit() |> max()

      actual_size_i <- if(i != nAges){age.limits[i + 1] - age.limits[i]} else{max_age - age.limits[i]}
      actual_size_j <- if(j != nAges){age.limits[j + 1] - age.limits[j]} else{max_age - age.limits[j]}

      div_i <- actual_size_j / actual_size_i

      matrix_out[i, j] <- contacts / div_i
    }
  }

  return(list("default_matrix_per_capita" = mixing,
              "adjusted_matrix_per_capita" = matrix_out))
}

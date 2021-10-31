require(deSolve)

npiRate = function(npi_rates, time) {
  t = time + 1
  if (t < length(npi_rates)) {
    npi_rates[t]
  } else {
    tail(npi_rates, 1)
  }
}

#' Numerically solve an SEIR model with parameters take from
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7445013/.
#' The model starts with 1 exposed individual.
#' 
#' @param population the number of individuals
#' @param r_0 the base reproductive rate.
#' @param end_day when the simulation ends.
#' @param npi_rates a sequence of multiplicative modifiers to |r_0|, so that on
#'     day n, the reproducive rate will be |r_0*npi_rates[n]|.
seirModel = function(population, r_0, end_day, npi_rates) {
# see https://www.freecodecamp.org/news/how-to-model-an-epidemic-with-r/ for code
  gamma_1 = 1 / 3.69
  gamma_2 = 1 / 3.48
  SEIR <- function(time, current_state, params) {
    npi_rate = npiRate(npi_rates, time)
    t = time - 1
    beta = r_0 * gamma_2 * npi_rate
    with(as.list(c(current_state, params)), {
      dS <- -(beta * S * I) / n
      dE <- (beta * S * I) / n - sigma * E
      dI <- sigma * E - gamma * I
      dR <- gamma * I
      
      return(list(c(dS, dE, dI, dR)))
    })
  }
  params <- c(sigma = gamma_1,
              gamma = gamma_2,
              n = POPULATION)
  
  initial_state <- c(S = POPULATION - 1,
                     E = 1,
                     I = 0,
                     R = 0)
  
  times <- 0:end_day
  
  ode(initial_state, times, SEIR, params)
}

#' Uses \code{seirModel} to find he first day that infections reaches 1.
findFirstInfection = function(population, r_0, end_day) {
  model = as.data.frame(seirModel(population, r_0, end_day, c(1)))
  infected = model$I[1:end_day]
  min(which(infected >= 1))
}

#' Generates a sequences on numbers, starting with a run of 1 and switching
#' to a run of \code{rate}s from \code{start} to \code{end}.
#' This represents an NPI that has a constant effect between the start and end.
NPI = function(rate, start, end) {
  c(rep(1, start - 1), rep(rate, end - start + 1))
}

#' Generates a dataframe of NPIs with sequential names.
NPIs = function(npis) {
  # This is so ugly.
  final_frame = data.frame("npi_1" = npis[[1]])
  
  i = 2
  while (i <= length(npis)) {
    npi = npis[i]
    frame = data.frame("npi" = npi)
    names(frame)[1] <- paste("npi_", i, sep="")
    final_frame = cbind(final_frame, frame)
    i = i + 1
  }
  final_frame
}



#' Returns a sequence of numbers produced by multiplicatively merging multiple
#' NPIs.
mergeNPIs = function(npis) {
  Reduce(npis, f = "*")
}
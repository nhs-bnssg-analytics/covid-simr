#' Covid-19 ICU bed-occupancy simulation
#'
#' @param tscenario (\code{string}), name to give the simulation scenario
#' @param cases (\code{data.frame}), hosptialisation data, one column named \code{date} (string), and one column named \code{hospitalisations} which is an integer recording new cases for that date.
#' @param galpha (\code{positive numeric}), median length of stay for admitted covid19 pateints, to be estimated by user from data or other source
#' @param gbeta (\code{positive numeric}), 0.95 quantile length of stay for admitted covid19 patients, to be estimated by user from data or other source
#' @param cap (\code{positive integer}), maximum number of patients who can be concurrently admitted to the hospital/unit
#' @param pfat (\code{numeric} in range 0 (no chance of death)`  to 1 (certain death)), probability that a patient who arrives but is rejected/cannot be admitted (because there are no beds free) dies
#' @param tol (\code{numeric} between 0 and 100, with 0 corresponding to full confidence inputs are correct, with 100 least level of confidence), subjective assessment of reliability of input hospitalisation estiamtes
#' @param nreps (\code{positive integer}), number of simulation replications to perform - larger number means better results, but longer time to compute
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
covid_simr <- function(tscenario,
                       cases,
                       galpha,
                       gbeta,
                       cap,
                       pfat,
                       tol = 25,
                       nreps = 100) {

  # pre-alloacte symbols to NULL to avoid NOTES from CMD check
  dates <- metric <- value <- fill <- time <- outcome <- type <- value_cum <-
  q01 <- q025 <- q05 <- q15 <- q30 <- q70 <- q85 <- q95 <- q975 <- q99 <-
  hospitalisations <- occ_q01 <- occ_q025 <- occ_q05 <- occ_q15 <- occ_q30 <-
  occ_q70 <- occ_q85 <- occ_q95 <- occ_q975 <- occ_q99 <- occ_median <-
  rejected_died_q01 <- rejected_died_q025 <- rejected_died_q05 <-
  rejected_died_q15 <- rejected_died_q30 <- rejected_died_q70 <-
  rejected_died_q85 <- rejected_died_q95 <- rejected_died_q975 <-
  rejected_died_q99 <- rejected_died_median <-
  admitted_q01 <- admitted_q025 <- admitted_q05 <- admitted_q15 <-
  admitted_q30 <- admitted_q70 <- admitted_q85 <- admitted_q95 <-
  admitted_q975 <- admitted_q99 <- admitted_median <-
  rejected_survived_q01 <- rejected_survived_q025 <- rejected_survived_q05 <-
  rejected_survived_q15 <- rejected_survived_q30 <- rejected_survived_q70 <-
  rejected_survived_q85 <- rejected_survived_q95 <- rejected_survived_q975 <-
  rejected_survived_q99 <- rejected_survived_median <- NULL

  # calculate lognormal distribution parameters from given information (los_mean and los_95)
  #erfinv<-function(x) qnorm((1+x)/2)/sqrt(2)
  #meanlog<-log(los_median)
  #sdlog<-((sqrt(2)*erfinv(0.9))^-1)*log(los_95/los_median)

  # determine whether epidemic curve uncertainty is to be considered
  if (tol > 0 & tol <= 100) {
    uncert <- TRUE
  } else {
    uncert <- FALSE
  }


  #This section Reads in given daily hospitalisations esimtate from a CSV (based on epidemic curve data, not caculated by this package).
  #    Generates confidence intervals by sampling noise from normal distrubtion, scaling that by a tolerance parameter (given as
  #    a percentage uncertainly as an integer in the range 0-100)  and scaling the orignial hospitalisations estimate using
  #    that (i.e. orignial value * 1+scaling, where the scaling can be +-, floored at -1). The output of this function is
  #    then used internally as the input dataframe by the run_sim_fun function in this package.


  #create three column dataframe:
  #date (date format)
  #name of metric type (e.g. hospitalisations) (character format)
  #value for date/type pair (numeric format)
  inputs <- cases
  inputs <- inputs %>%
    dplyr::mutate(dates = as.Date(dates, "%d/%m/%Y")) %>%
    tidyr::pivot_longer(cols = -dates,
                 names_to = "metric",
                 values_to = "value")

  #If a nonzero tolerance has been set (strictly between 0 and 100 - for percent)
  #then uncertaintly about the epidemic curve is assumed
  #and daily arrivals are modified by adding some normally distrubuted noise
  #(floored at 0, and increased by a maximum of the upper positive limit of a sample from
  #the normal distrubtion about zero with sd 1/3), modified by choice of tolerance parameter (to sim, not normal)
  #For each date, is is repeated reps number of times and the results use to derive quantile estimates for possible arrivals
  #to be used as confidence intervals on the original central estimate
  if (uncert == TRUE) {
    inputs_sim <- lapply(1:nreps, function(i) {
      set.seed(i)
      scaler <-
        tol / 100 * max(stats::rnorm(1, mean = 0, sd = 1 / 3),-1)
      data.frame(
        dates = inputs$dates,
        metric = paste0("rep", i),
        value = inputs$value * (1 + scaler)
      )
    })
    inputs_sim <- do.call("rbind", inputs_sim)
    inputs <- rbind(inputs, inputs_sim)
    #get CIs
    inputs_sim_ci <- inputs %>%
      dplyr::group_by(dates) %>%
      dplyr::filter(substr(metric, 1, 3) == "rep") %>%
      dplyr::summarise(
        q30 = stats::quantile(value, 0.3),
        q15 = stats::quantile(value, 0.15),
        q05 = stats::quantile(value, 0.05),
        q025 = stats::quantile(value, 0.025),
        q01 = stats::quantile(value, 0.01),
        q70 = stats::quantile(value, 0.7),
        q85 = stats::quantile(value, 0.85),
        q95 = stats::quantile(value, 0.95),
        q975 = stats::quantile(value, 0.975),
        q99 = stats::quantile(value, 0.99)
      ) %>%
      tidyr::pivot_longer(cols = -dates,
                   names_to = "metric",
                   values_to = "value")
    #adds original inputs datframe now has multiple rows for each date
    #one of orignial hospitalisations estimate,
    #for each rep, one with a modified hospitalisations estimate (scaled by the normal sample - itself scaled by the tolerance value)
    #and one for each of the confidence quantiles derived rmo the modified hospitalisations
    inputs <- rbind(inputs, inputs_sim_ci)
  }


  simfn <- function(rep) {
    #fix the random number stream for this replication
    #(so that if run again without changing any substantive parameters, the results will be the same)
    #uses rep number so that a different random number stream is initialised for each replication
    set.seed(rep)
    dmax <- length(unique(inputs$dates))
    #get num arrivals by day
    #then spread them out over time within days
    arr_times <-
      unlist(sapply(1:dmax, function(x) {
        #hospitalisation estimate for each day
        #either the original value from the epidemic curve
        #or, if a tolrenace value has been set, the modified value for the given replication
        narr <-
          inputs$value[which(inputs$metric == ifelse(uncert == TRUE, paste0("rep", rep), "hospitalisations"))][x]
        #For each date, a (re)sampled total hospitalisation number, using the original value as the rate parameter
        #to the Poisson distribution
        narr.pois <- stats::rpois(1, lambda = narr)
        #spread these arrivals out randomly throughout the given day, ordering first to last
        sort(stats::runif(narr.pois, 0, 1) + x - 1)
      }))
    #"cal" - the simulation schedule
    #used internally within the simulation loop to track events and determine admissions, rejections, service ends
    #dataframe listing the arrivals for the given day in order:
    #"id" - positive integer, index column
    #"time" - numeric, decimal time between first and last arrivals
    #"event" - character, describes type of event - intinally, here, that is just a patient arrival
    cal <-
      data.frame(id = 1:length(arr_times),
                 time = arr_times,
                 event = "arrival")
    #The simulation clock - used to update both the simulation schedule "cal" and (by rounding up) the outputs dataframe "res"
    tx <- 0
    #initialise outcome dataframe "res" for each day
    #number of rows is approx number of days plus 1/3
    #each row tracks the state of the hospital in terms of total occupancy at, and total admissions and rejections (separated into
    #survived and died) by that time index
    #dataframe cols:
    #"time" integer, tracks days from start to end + ~1/3
    #"occ", integer, current total numeber of beds occupied
    #"admitted", integer, cumulative number patients admitted since start
    #"rejected_died", integer, cumulative number of patients who were refused admission and as a result died, since start
    #"rejected_survived", integer, cumulative number of patients who were refused admission but survived, since start
    #CHANGE - initialise results with occupancy NA instead of occupancy zero - then handle occupancies at end ####
    res <-
      data.frame(
        time = 1:round(dmax * 1.33),
        occ = NA,
        admitted = 0,
        rejected_died = 0,
        rejected_survived = 0
      )
    #initialise hospital occupancy before starting simulation loop ####
    occ <- 0  #occupancy (number in unit)
    while (nrow(cal) > 0 &
           tx < max(res$time)) {
      #find all rows in the schedule which have not happened yet (time greater than simulation clock time tx)
      #and which are either arrivals or service completions
      #return as a vector of row indices to the schedule dataframe "cal"
      ind1 <- which(cal$time > tx &
                      cal$event %in% c("arrival", "endsrv"))
      #find the row of the schedule which corresponds to the earliest of those events and return its index
      ind <- ind1[which.min(cal$time[ind1])]
      #advance the simulation clock to that time
      tx <- cal$time[ind]
      #return the start time of the next full schedule day after the current simulation clock time
      tx_day <- ceiling(tx)
      #if the next day is beyond the time limit (the maximum arrivals schedule date plus approx 1/3), then stop the simulation loop
      #,and proceed to return the results dataframe as the function output
      if (tx_day > max(res$time))
        break
      #otherwise, process events 3 types (1.arrival-admitted, 2.arrival-rejected[died/survived], and 3.end of service)
      #Arrival
      if (cal$event[ind] == "arrival") {
        #1. Arrival-admission: if there is a free bed, they are admitted to the hospital
        if (occ < cap) {
          #increment the cumulative tally of admitted patients at the start of the next full calendar date by one in the results
          #dataframe. This is an output measure and will not be used to decide whether or not subsequent arrivals can be admitted
          res$admitted[tx_day] <-
            res$admitted[tx_day] + 1
          #add the start of this patient's length of stay as an event, to the end of the simulation schedule, at the decimal start time
          #(n.b. not the rounded up to start of next day time used to track occupancy)
          #POTENTIALLY SLOW, AND DONE LOTS OF TIMES ####
          cal <-
            rbind(cal,
                  data.frame(
                    id = cal$id[ind],
                    time = tx,
                    event = "startsrv"
                  ))
          #for the patient just admitted, sample a length of stay from the lognormal distribution using the function input parameters,
          #which scale it to expected median length of stay for all patients (dervived externally by user
          #from case data or other source)
          los <-
            stats::rgamma(1, shape = galpha, rate = gbeta)
          #add the end of this patient's length of stay as an event to the end of the simulation schedule, specifying the event type
          #and the decimal time (given by the patient's service start time plus the LOS just sampled)
          #POTENTIALLY SLOW, AND DONE LOTS OF TIMES ####
          cal <-
            rbind(cal,
                  data.frame(
                    id = cal$id[ind],
                    time = tx + los,
                    event = "endsrv"
                  ))
          #increment the current actual occupancy (at this point in time by the simulation clock tx) by one
          #this will be used to check whether the next arrival can be admitted
          occ <- occ + 1
        } else {
          ##2. Arrival-rejection: There is not a free bed, so the patient is rejected
          #remove the event row from the simulation schedule
          #POTENTIALLY SLOW, AND DONE LOTS OF TIMES ####
          cal <-
            cal[-which(cal$id == cal$id[ind]),]
          #Does the rejected patient die? Sample from the range [0,1] with each value being equally probable. If the sample is less
          #than the user-inputed estiamte of a rejected patient dying, they die
          if (stats::runif(1, 0, 1) < pfat) {
            #2.1 Arrival-rejection-died
            #increment the running tally of patients who have been rejected and died (at the start of the next full day, in
            #the results dataframe). This is an output, it is not used to decide future events within the simulation loop.
            res$rejected_died[tx_day] <-
              res$rejected_died[tx_day] + 1
          } else {
            #2.2 Arrival-rejection-survived
            #sample from unif was greater than the given chance of dying, so patient survives
            #increment the running tally fo rejected-survived in the outputs dataframe (res) a the row corresponding to the start of
            #the next day. Thisis an output, it is not used to decide future events within the simulation loop.
            res$rejected_survived[tx_day] <-
              res$rejected_survived[tx_day] + 1
          }
        }
        #3. End of service
        #Note that this may be survival or death - but not death due to capacity constraint
        #i.e. it is normal outcome in event of sufficient capacity for this patient
      } else if (cal$event[ind] == "endsrv") {
        #Remove row for this patient from the simulation schedule
        #POTENTIALLY SLOW, AND DONE LOTS OF TIMES ####
        cal <-
          cal[-which(cal$id == cal$id[ind]),]    # this is no longer needed since p measures recorded separately at each iteration
        #de-increment occuapncy at current simulation time
        occ <- occ - 1
      }
      #sort the simulation schedule into ascending time order, by the time column
      #THIS WILL BE SLOW - AND DONE LOTS OF TIMES ####
      cal <- cal[order(cal$time),]
      #save results, extract performance measures
      #Each row in the output dataframe "res" up to this point has measures for that calendar date - cumulative admissions,
      #cumulative rejections who survived, cumulative rejections who died, and current occupancy for that date
      #but the occupancy remains zero
      #this sets the occupancy to value it had during the last event which occured during that day.
      res$occ[tx_day] <- occ
    }
    #for time periods where there was no event, the occupancy field in the results will be NA
    #replace this to be (i) zero if it is the first day
    #                   (i) the same occupancy as the last non-NA day otherwise
    if (is.na(res$occ[1])) {
      res$occ[1] <- 0
    }
    res <- res %>% fill(occ)
    return(cbind(data.frame(rep = rep), res))
  }

  #Now run nreps replications of the simulation function to get results for nreps possible
  #actualisations of the parameters

  #Do the runs in parallel to speed up processing time
  #Set up the clusters (set to one fewer than the total number of logical cores in the user's comput)
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  #send necessary objects created in the global environment to the clusters
  parallel::clusterExport(
    cl = cl,
    varlist = c("uncert", "inputs", "galpha", "gbeta", "cap", "pfat"),
    envir = environment()
  )
  #make the libraries used within teh simfn function available to the clusters
  parallel::clusterEvalQ(cl = cl,
               c(library(tidyr),
                 library(dplyr)))
  #return results of the replications
  RES <- parallel::parLapply(cl, 1:nreps, simfn)
  parallel::stopCluster(cl)

  #combine the simulation outputs into a single dataframe
  outputs_sim <- do.call("rbind", RES)
  outputs <- outputs_sim %>%
    tidyr::pivot_longer(
      cols = -c(rep, time),
      names_to = "outcome",
      values_to = "value"
    ) %>%
    dplyr::group_by(time, outcome) %>%
    dplyr::summarise(
      mean = mean(value),
      median = stats::median(value),
      q30 = stats::quantile(value, 0.3),
      q25 = stats::quantile(value, 0.25),
      q15 = stats::quantile(value, 0.15),
      q05 = stats::quantile(value, 0.05),
      q025 = stats::quantile(value, 0.025),
      q01 = stats::quantile(value, 0.01),
      q70 = stats::quantile(value, 0.7),
      q75 = stats::quantile(value, 0.75),
      q85 = stats::quantile(value, 0.85),
      q95 = stats::quantile(value, 0.95),
      q975 = stats::quantile(value, 0.975),
      q99 = stats::quantile(value, 0.99)
    ) %>%
    tidyr::pivot_longer(
      cols = -c(time, outcome),
      names_to = "type",
      values_to = "value"
    ) %>%
    tidyr::pivot_wider(names_from = c(outcome, type),
                values_from = value) %>%
    dplyr::mutate(dates = seq(
      from = min(inputs$dates) + time - 1,
      by = 1,
      length.out = 1
    ))

  outputs_cum <- outputs_sim %>%
    tidyr::pivot_longer(
      cols = -c(rep, time),
      names_to = "outcome",
      values_to = "value"
    ) %>%
    dplyr::group_by(outcome, rep) %>%
    dplyr::mutate(value_cum = cumsum(value)) %>%
    dplyr::select(-value) %>%
    dplyr::group_by(time, outcome) %>%
    dplyr::summarise(
      mean = mean(value_cum),
      median = stats::median(value_cum),
      q30 = stats::quantile(value_cum, 0.3),
      q25 = stats::quantile(value_cum, 0.25),
      q15 = stats::quantile(value_cum, 0.15),
      q05 = stats::quantile(value_cum, 0.05),
      q025 = stats::quantile(value_cum, 0.025),
      q01 = stats::quantile(value_cum, 0.01),
      q70 = stats::quantile(value_cum, 0.7),
      q75 = stats::quantile(value_cum, 0.75),
      q85 = stats::quantile(value_cum, 0.85),
      q95 = stats::quantile(value_cum, 0.95),
      q975 = stats::quantile(value_cum, 0.975),
      q99 = stats::quantile(value_cum, 0.99)
    ) %>%
    tidyr::pivot_longer(
      cols = -c(time, outcome),
      names_to = "type",
      values_to = "value_cum"
    ) %>%
    tidyr::pivot_wider(names_from = c(outcome, type),
                values_from = value_cum) %>%
    dplyr::mutate(dates = seq(
      from = min(inputs$dates) + time - 1,
      by = 1,
      length.out = 1
    ))


  ####################################################################################################################################################################

  #plot the results

  #daily hospitalisations (including resampling and tolerance to give confidence ranges)
  plot1 <- inputs %>%
    tidyr::pivot_wider(names_from = "metric", values_from = "value") %>%
    ggplot2::ggplot(ggplot2::aes(x = dates)) +
    ggplot2::labs(title = "Inputs: Cases requiring hospitalisation (per day)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::scale_x_date(
      date_breaks = "months",
      date_labels = "%b-%y",
      limits = c(min(outputs$dates), max(outputs$dates))
    )
  if (uncert == TRUE) {
    plot1 <- plot1 +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q01, ymax = q025),
                  fill = "grey",
                  alpha = 0.6) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q025, ymax = q05),
                  fill = "dodgerblue1",
                  alpha = 0.6) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q05, ymax = q15),
                  fill = "dodgerblue2",
                  alpha = 0.6) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q15, ymax = q30),
                  fill = "dodgerblue3",
                  alpha = 0.6) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q30, ymax = q70),
                  fill = "dodgerblue4",
                  alpha = 0.6) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q70, ymax = q85),
                  fill = "dodgerblue3",
                  alpha = 0.6) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q85, ymax = q95),
                  fill = "dodgerblue2",
                  alpha = 0.6) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q95, ymax = q975),
                  fill = "dodgerblue1",
                  alpha = 0.6) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q975, ymax = q99),
                  fill = "grey",
                  alpha = 0.6)
  }
  plot1 <- plot1 +
    ggplot2::geom_line(ggplot2::aes(y = hospitalisations))

  #bed occupancy over simulation period
  plot2 <- outputs %>%
    ggplot2::ggplot(ggplot2::aes(x = dates)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = occ_q01, ymax = occ_q025),
                fill = "darkgray",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = occ_q025, ymax = occ_q05),
                fill = "darkorchid1",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = occ_q05, ymax = occ_q15),
                fill = "darkorchid2",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = occ_q15, ymax = occ_q30),
                fill = "darkorchid3",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = occ_q30, ymax = occ_q70),
                fill = "darkorchid4",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = occ_q70, ymax = occ_q85),
                fill = "darkorchid3",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = occ_q85, ymax = occ_q95),
                fill = "darkorchid2",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = occ_q95, ymax = occ_q975),
                fill = "darkorchid1",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = occ_q975, ymax = occ_q99),
                fill = "darkgray",
                alpha = 0.6) +
    ggplot2::geom_hline(yintercept = cap,
               linetype = "dashed",
               colour = "darkgray") +
    ggplot2::geom_line(ggplot2::aes(y = occ_median)) +
    ggplot2::labs(title = paste0("Outputs: Occupied beds (capacity inputted=", cap, ")")) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::scale_x_date(date_breaks = "months", date_labels = "%b-%y")

  #deaths resulting from insufficient capacity over simulation period
  plot3 <- outputs %>%
    ggplot2::ggplot(ggplot2::aes(x = dates)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q01, ymax = rejected_died_q025),
      fill = "yellow1",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q025, ymax = rejected_died_q05),
      fill = "orange2",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q05, ymax = rejected_died_q15),
      fill = "orangered3",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q15, ymax = rejected_died_q30),
      fill = "red3",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q30, ymax = rejected_died_q70),
      fill = "red2",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q70, ymax = rejected_died_q85),
      fill = "red3",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q85, ymax = rejected_died_q95),
      fill = "orangered3",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q95, ymax = rejected_died_q975),
      fill = "orange2",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q975, ymax = rejected_died_q99),
      fill = "yellow1",
      alpha = 0.6
    ) +
    ggplot2::geom_line(ggplot2::aes(y = rejected_died_median)) +
    ggplot2::labs(title = "Outputs: Capacity-dependent deaths (per day)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::scale_x_date(date_breaks = "months", date_labels = "%b-%y")

  #cumulative total of admitted patients over simulation period
  plot4 <- outputs_cum %>%
    ggplot2::ggplot(ggplot2::aes(x = dates)) +
    #admitted
    ggplot2::geom_ribbon(ggplot2::aes(ymin = admitted_q01, ymax = admitted_q025),
                fill = "darkorange",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = admitted_q025, ymax = admitted_q05),
                fill = "darkorange1",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = admitted_q05, ymax = admitted_q15),
                fill = "darkorange2",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = admitted_q15, ymax = admitted_q30),
                fill = "darkorange3",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = admitted_q30, ymax = admitted_q70),
                fill = "darkorange4",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = admitted_q70, ymax = admitted_q85),
                fill = "darkorange3",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = admitted_q85, ymax = admitted_q95),
                fill = "darkorange2",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = admitted_q95, ymax = admitted_q975),
                fill = "darkorange1",
                alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = admitted_q975, ymax = admitted_q99),
                fill = "darkorange",
                alpha = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = admitted_median)) +
    ggplot2::labs(title = "Outputs: Cumulative total for admitted patients") +
    ggplot2::scale_x_date(date_breaks = "months", date_labels = "%b-%y") +
    ggplot2::ylim(
      0,
      max(
        outputs_cum$admitted_q99,
        outputs_cum$rejected_died_q99,
        outputs_cum$rejected_survived_q99
      )
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  #cumulative total of patients who could not be admitted because of capacity constraints but did
  #surivive, over simulation period
  plot5 <- outputs_cum %>%
    ggplot2::ggplot(ggplot2::aes(x = dates)) +
    #rejected-survived
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_survived_q01, ymax = rejected_survived_q025),
      fill = "chartreuse",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_survived_q025, ymax = rejected_survived_q05),
      fill = "chartreuse1",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
    ggplot2::aes(ymin = rejected_survived_q05, ymax = rejected_survived_q15),
      fill = "chartreuse2",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_survived_q15, ymax = rejected_survived_q30),
      fill = "chartreuse3",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_survived_q30, ymax = rejected_survived_q70),
      fill = "chartreuse4",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_survived_q70, ymax = rejected_survived_q85),
      fill = "chartreuse3",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_survived_q85, ymax = rejected_survived_q95),
      fill = "chartreuse2",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_survived_q95, ymax = rejected_survived_q975),
      fill = "chartreuse1",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_survived_q975, ymax = rejected_survived_q99),
      fill = "chartreuse",
      alpha = 0.6
    ) +
    ggplot2::geom_line(ggplot2::aes(y = rejected_survived_median)) +
    ggplot2::labs(title = "Outputs: Cumulative total for rejected and survived") +
    ggplot2::scale_x_date(date_breaks = "months", date_labels = "%b-%y") +
    ggplot2::ylim(
      0,
      max(
        outputs_cum$admitted_q99,
        outputs_cum$rejected_died_q99,
        outputs_cum$rejected_survived_q99
      )
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  #cumulative total of patients who died as a result of insufficent capacity over simulation period
  plot6 <- outputs_cum %>%
    ggplot2::ggplot(ggplot2::aes(x = dates)) +
    #rejected-died
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q01, ymax = rejected_died_q025),
      fill = "yellow1",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q025, ymax = rejected_died_q05),
      fill = "orange2",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q05, ymax = rejected_died_q15),
      fill = "orangered3",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q15, ymax = rejected_died_q30),
      fill = "red3",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q30, ymax = rejected_died_q70),
      fill = "red2",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q70, ymax = rejected_died_q85),
      fill = "red3",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q85, ymax = rejected_died_q95),
      fill = "orangered3",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q95, ymax = rejected_died_q975),
      fill = "orange2",
      alpha = 0.6
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rejected_died_q975, ymax = rejected_died_q99),
      fill = "yellow1",
      alpha = 0.6
    ) +
    ggplot2::geom_line(ggplot2::aes(y = rejected_died_median)) +
    ggplot2::labs(title = "Outputs: Cumulative total for rejected and died") +
    ggplot2::scale_x_date(date_breaks = "months", date_labels = "%b-%y") +
    ggplot2::ylim(
      0,
      max(
        outputs_cum$admitted_q99,
        outputs_cum$rejected_died_q99,
        outputs_cum$rejected_survived_q99
      )
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )


  #string to use in filename to save outputs
  filename_ext <- paste0("_", tscenario)

  #save plots above as a grid in a pdf in same folder location as script is run from
  grDevices::pdf(paste0("output_plot", filename_ext, ".pdf"),
      height = 8,
      width = 12)
  gridExtra::grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, nrow = 2)
  grDevices::dev.off()
  grDevices::png(
    paste0("output_plot", filename_ext, ".png"),
    height = 8,
    width = 12,
    units = "in",
    res = 800
  )
  gridExtra::grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, nrow = 2)
  grDevices::dev.off()

  #save inputs (with tolerances/confinece range) and outputs as csvs,
  #in same location script was run from
  #write.csv(inputs,paste0("inputs",filename_ext,".csv"),row.names=FALSE)
  utils::write.csv(
    cbind(data.frame(cap = cap), as.data.frame(outputs)),
    paste0("outputs", filename_ext, ".csv"),
    row.names = FALSE
  )
  utils::write.csv(
    cbind(data.frame(cap = cap), as.data.frame(outputs_cum)),
    paste0("outputs_cum", filename_ext, ".csv"),
    row.names = FALSE
  )

  #return outputs dataframe
  return(outputs)

}

library(randtoolbox)
library(dplyr)
library(rsm)

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

iterations <- 10

search_space <- NULL
results <- NULL

sobol_dim <- 2
starting_sobol_n <- 20

sobol_n <- starting_sobol_n

level_min <- 0
level_max <- 1

supersim_dir <- "../../supersim/"

measure <- function(configuration) {
    start_time <- as.integer(format(Sys.time(), "%s"))

    cmd <- paste(supersim_dir,
                 "bin/supersim",
                 " ../settings/settings.json",
                 " workload.message_log.file=string=./output.mpf.gz",
                 " workload.applications[0].blast_terminal.request_injection_rate=float=",
                 configuration$Application_1,
                 " workload.applications[0].blast_terminal.request_injection_rate=float=",
                 configuration$Application_2,
                 sep = "")

    print(cmd)
    quiet(system(cmd))

    elapsed_time <- as.integer(format(Sys.time(), "%s")) - start_time

    cmd <- paste(supersim_dir,
                 "scripts/ssparse/bin/ssparse",
                 " -p application_1.csv -f +app=0 output.mpf.gz",
                 sep = "")

    print(cmd)
    system(cmd)

    application_1 <- read.csv("application_1.csv", header = FALSE)
    names(application_1) <- c("packet_start", "packet_finish",
                              "id_1", "id_2", "id_3")

    application_1$id <- "Application_1"

    application_1 <- application_1 %>%
        filter(row_number() == 1 | row_number() == n())

    cmd <- paste(supersim_dir,
                 "scripts/ssparse/bin/ssparse",
                 " -p application_2.csv -f +app=1 output.mpf.gz",
                 sep = "")

    print(cmd)
    system(cmd)

    application_2 <- read.csv("application_2.csv", header = FALSE)
    names(application_2) <- c("packet_start", "packet_finish",
                              "id_1", "id_2", "id_3")

    application_2$id <- "Application_2"
    application_2 <- application_2 %>%
        filter(row_number() == 1 | row_number() == n())

    system("rm output.mpf.gz application_1.csv application_2.csv")

    return(data.frame(Application_1 = c(application_1[2, "packet_finish"] -
                                        application_1[1, "packet_start"]),
                      Application_2 = c(application_2[2, "packet_finish"] -
                                        application_2[1, "packet_start"]),
                      injection_rate_1 = c(configuration$Application_1),
                      injection_rate_2 = c(configuration$Application_2),
                      duration = c(elapsed_time)))
}

mean_execution_time <- function(results) {
    return(results %>%
           mutate(performance_metric = (Application_1 + Application_2) / 2))
}

weighted_injection_time <- function(results) {
    execution_time_weight <- 1e-8 # Find a better weight
    return(results %>%
           mutate(performance_metric = (1 / (injection_rate_1 + injection_rate_2)) +
                      (execution_time_weight * ((Application_1 + Application_2) / 2))))
}

map_intervals <- function(x, interval_from, interval_to) {
    new_x <- x - interval_from[1]
    new_x <- new_x / (interval_from[2] - interval_from[1])
    new_x <- new_x * (interval_to[2] - interval_to[1])
    new_x <- new_x + interval_to[1]

    return(new_x)
}

for(i in 1:iterations){
    temp_sobol <- sobol(n = 2 * sobol_n,
                        dim = sobol_dim,
                        scrambling = 3,
                        seed = as.integer((99999 - 10000) * runif(1) + 10000),
                        init = TRUE)

    rm(temp_sobol)
    quiet(gc())

    design <- sobol(n = 2 * sobol_n,
                    dim = sobol_dim,
                    scrambling = 3,
                    seed = as.integer((99999 - 10000) * runif(1) + 10000),
                    init = FALSE)

    df_design <- data.frame(design)

    names(df_design) <- c("Application_1",
                          "Application_2")

    design <- design %>% filter(Application_1 + Application_2 < 1.0)

    print(df_design)

    # df_design <- map_intervals(df_design, c(0.0, 1.0), c(0.0, 0.5))

    # print(df_design)

    current_results <- bind_rows(lapply(1:nrow(df_design), function(row) { measure(df_design[row, ]) }))
    #current_results <- mean_execution_time(current_results)
    current_results <- weighted_injection_time(current_results)

    print(current_results)

    best_points <- filter(current_results, performance_metric == min(performance_metric))

    best_points$id <- i
    best_points$points <- sobol_n

    if(is.null(results)){
        results <- best_points
    } else{
        results <- bind_rows(results, best_points)
    }

    if(is.null(search_space)){
        search_space <- current_results
    } else{
        search_space <- bind_rows(search_space, current_results)
    }

    write.csv(results,
              paste("rs_",
                    starting_sobol_n,
                    "_samples_",
                    iterations,
                    "_iterations.csv",
                    sep = ""),
              row.names = FALSE)

    write.csv(search_space,
              paste("rs_",
                    starting_sobol_n,
                    "_samples_",
                    iterations,
                    "_iterations_search_space.csv",
                    sep = ""),
              row.names = FALSE)
}

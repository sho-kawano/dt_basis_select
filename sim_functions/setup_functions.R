library(tidycensus)
library(LaplacesDemon)
library(tidyverse)

# CHANGES (2025.06)
# * took out the county stuff since we won't use it 
# * folder name is different for this set of sims
# Changed inputs from (data, y_name, transform, ncomps) 
#                 to  (z_mean, z_var, ncomps)

setup_comp <- function(z_mean, z_var, ncomps, results_dir){
  # set up the sim_files folder containing all the results
  base_folder = file.path(getwd(), results_dir)
  dir.create(base_folder, showWarnings = FALSE)
  
  # main loop
  for (sim in 1:ncomps){
    # create simulation folder
    sim_folder = file.path(base_folder, sprintf("comparison_%03d", sim))
    dir.create(sim_folder, showWarnings = FALSE)

    # create dt_1fold, dt_5fold, emp_sim, & fit_on_z  within each simulation folder
    dir.create(file.path(sim_folder, "dt_1fold"), showWarnings = FALSE)
    dir.create(file.path(sim_folder, "dt_5fold"), showWarnings = FALSE)
    dir.create(file.path(sim_folder, "emp_sim"), showWarnings = FALSE)
    dir.create(file.path(sim_folder, "fit_on_z"), showWarnings = FALSE)

    # generates 'z' - observed direct estimate
    set.seed(sim)
    z = rmvn(mu=z_mean, Sigma=diag(z_var))

    # save results
    saveRDS(z, file=file.path(sim_folder, "z.RDS"))
  }
}


# NOTE: niters is the number of iteration per empirical simulation study
esim_helper <- function(comp_no, w_var, niters, results_dir){
  # set up the folder
  base_folder = file.path(getwd(), results_dir)
  sim_folder = file.path(base_folder, sprintf("comparison_%03d", comp_no))
  
  # load direct estimates / truth for a given emp.study
  z = readRDS(file.path(sim_folder, "z.RDS")) %>% as.numeric()
  
  # generate synthetic data for one simulation study
  for (iter in 1:niters){
    sub_folder = sprintf("%03d", iter)
    sub_path = file.path(sim_folder, "emp_sim", sub_folder)
    dir.create(sub_path, recursive = TRUE, showWarnings = FALSE)
    set.seed(iter*2)
    w = rmvn(mu=z, Sigma=diag(w_var))
    # save results
    saveRDS(w, file=file.path(sub_path, "w.RDS"))
  }
}

# iterates over each comparison
setup_esim <- function(w_var, ncomps, results_dir){
  lapply(1:ncomps, function(k){ esim_helper(comp_no=k, w_var=w_var, niters=100, results_dir=results_dir)})
  return("Setup complete!")
}

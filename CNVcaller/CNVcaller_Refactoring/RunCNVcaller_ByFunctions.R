library('futile.logger')
library(infercnv)

# Load functions
source('CNVcaller/inferCNV_ops_trimmed.R')
source('CNVcaller/inferCNV.R')

refactored_run <- function(infercnv_obj,
                
                # gene filtering settings
                #cutoff=1,
                min_cells_per_gene=3,
                
                #out_dir=NULL,
                
                ## smoothing params
                window_length=101,
                smooth_method=c('pyramidinal', 'runmeans', 'coordinates'),
                
                num_ref_groups=NULL,
                ref_subtract_use_mean_bounds=TRUE,
                
                # observation cell clustering settings
                #cluster_by_groups=FALSE,
                cluster_references=TRUE,
                k_obs_groups=1,
                
                hclust_method='ward.D2',
                
                max_centered_threshold=3, # or set to a specific value or "auto", or NA to turn off
                scale_data=FALSE,
                
                ## HMM opts
                #HMM=FALSE, # turn on to auto-run the HMM prediction of CNV levels
                ## tumor subclustering opts
                HMM_transition_prob=1e-6,
                HMM_report_by=c("subcluster","consensus","cell"),
                HMM_type=c('i6', 'i3'),
                HMM_i3_pval=0.05,
                HMM_i3_use_KS=TRUE,
                BayesMaxPNormal=0.5,
                
                ## some experimental params
                #sim_method=c('meanvar', 'simple', 'splatter'), ## only meanvar supported, others experimental
                sim_method='meanvar',
                sim_foreground=FALSE, ## experimental
                reassignCNVs=TRUE,
                
                
                ## tumor subclustering options
                analysis_mode=c('samples', 'subclusters', 'cells'), # for filtering and HMM
                tumor_subcluster_partition_method=c('random_trees', 'qnorm', 'pheight', 'qgamma', 'shc'),
                tumor_subcluster_pval=0.1,
                
                
                ## noise settings
                #denoise=FALSE,
                noise_filter=NA,
                sd_amplifier = 1.5,
                noise_logistic=FALSE, # if false, does complete 'noise' elimination.
                
                # outlier adjustment settings
                outlier_method_bound="average_bound",
                outlier_lower_bound=NA,
                outlier_upper_bound=NA,
                
                ## misc options
                final_scale_limits = NULL,
                final_center_val = NULL,
                debug=FALSE, #for debug level logging
                #num_threads = 4,
                plot_steps=FALSE,
                resume_mode=TRUE,
                png_res=300,
                plot_probabilities = TRUE,
                save_rds = TRUE,
                save_final_rds = TRUE,
                diagnostics = FALSE,
                
                ## experimental options
                remove_genes_at_chr_ends=FALSE,
                prune_outliers=FALSE,
                
                mask_nonDE_genes=FALSE,
                mask_nonDE_pval=0.05, # use permissive threshold
                test.use='wilcoxon',
                require_DE_all_normals="any",
                
                
                hspike_aggregate_normals = FALSE,
                
                no_plot = FALSE,
                no_prelim_plot = FALSE,
                output_format = "png",
                useRaster = TRUE,
                
                up_to_step=100
                
) {
  smooth_method = match.arg(smooth_method)
  HMM_type = match.arg(HMM_type)
  if (HMM && HMM_type == 'i6' && smooth_method == 'coordinates') {
    flog.error("Cannot use 'coordinate' smoothing method with i6 HMM model at this time.")
    stop("Incompatible HMM mode and smoothing method.")
  }
  if (smooth_method == 'coordinates' && window_length < 10000) {
    window_length=10000000
    flog.warn(paste0("smooth_method set to 'coordinates', but window_length ",
                     "is less than 10.000, setting it to 10.000.000. Please ",
                     "set a different value > 10.000 if this default does not seem suitable."))
  }
  
  HMM_report_by = match.arg(HMM_report_by)
  analysis_mode = match.arg(analysis_mode)
  if(analysis_mode == "cells" && HMM_report_by != "cell") {
    flog.warn(paste0("analysis_mode is \"cells\" but HMM_report_by is not, "),
              "changing HMM_report_by to \"cells\".")
    HMM_report_by = "cell"
  }
  tumor_subcluster_partition_method = match.arg(tumor_subcluster_partition_method)
  
  if (debug) {
    flog.threshold(DEBUG)
  } else {
    flog.threshold(INFO)
  }
  
  flog.info(paste("::process_data:Start", sep=""))
  
  #infercnv.env$GLOBAL_NUM_THREADS <- num_threads
  #if (is.null(out_dir)) {
    #flog.error("Error, out_dir is NULL, please provide a path.")
    #stop("out_dir is NULL")
  #}
  
  call_match = match.call()
  
  arg_names = names(call_match)
  arg_names = arg_names[-which(arg_names == "" | arg_names == "infercnv_obj")]
  
  for (n in arg_names) {
    call_match[[n]] = get(n)
  }
  
  current_args = as.list(call_match[-which(names(call_match) == "" | names(call_match) == "infercnv_obj")])
  
  if(out_dir != "." && !file.exists(out_dir)){
    flog.info(paste0("Creating output path ", out_dir))
    dir.create(out_dir)
  }
  
  infercnv_obj@options = c(infercnv_obj@options, current_args)
  
  #run_call <- match.call()
  call_match[[1]] <- as.symbol(".get_relevant_args_list")
  reload_info = eval(call_match)
  
  reload_info$relevant_args
  reload_info$expected_file_names
  
  skip_hmm = 0
  skip_mcmc = 0
  skip_past = 0
  
  if (resume_mode) {
    flog.info("Checking for saved results.")
    for (i in rev(seq_along(reload_info$expected_file_names))) {
      if (file.exists(reload_info$expected_file_names[[i]])) {
        flog.info(paste0("Trying to reload from step ", i))
        if ((i == 17 || i == 18 || i == 19) && skip_hmm == 0) {
          hmm.infercnv_obj = readRDS(reload_info$expected_file_names[[i]])
          if (!.compare_args(infercnv_obj@options, unlist(reload_info$relevant_args[1:i]), hmm.infercnv_obj@options)) {
            rm(hmm.infercnv_obj)
            invisible(gc())
          }
          else {
            hmm.infercnv_obj@options = infercnv_obj@options
            skip_hmm = i - 16
            flog.info(paste0("Using backup HMM from step ", i))
          }
        }
        else if (i != 17 && i != 18 && i != 19) {
          reloaded_infercnv_obj = readRDS(reload_info$expected_file_names[[i]])
          if (skip_past > i) { # in case denoise was found
            if (20 > i) { # if 21/20 already found and checked HMM too, stop
              break
            }
            else { # if 21 denoise found, don't check 20 maskDE
              next
            }
          }
          if (.compare_args(infercnv_obj@options, unlist(reload_info$relevant_args[1:i]), reloaded_infercnv_obj@options)) {
            options_backup = infercnv_obj@options
            infercnv_obj = reloaded_infercnv_obj # replace input infercnv_obj
            rm(reloaded_infercnv_obj) # remove first (temporary) reference so there's no duplication when they would diverge
            infercnv_obj@options = options_backup
            skip_past = i
            flog.info(paste0("Using backup from step ", i))
            if (i < 21) { # do not stop check right away if denoise was found to also allow to check for HMM
              break
            }
          }
          else {
            rm(reloaded_infercnv_obj)
            invisible(gc())
          }
        }
      }
    }
  }
  
  step_count = 0;
  reaching_step(up_to_step,step_count)
  # 1
  step_get_incoming_data(step_count,skip_past,save_rds,reload_info)
  reaching_step(up_to_step, step_count)
  # 2
  step_remove_low_expr_genes(step_count,cutoff,min_cells_per_gene,skip_past,reload_info)
  reaching_step(up_to_step, step_count)
  # 3
  step_normalize_by_seq_depth(step_count,HMM,HMM_type,sim_method,hspike_aggregate_normals,reload_info)
  reaching_step(up_to_step, step_count)
  
  #our step
  step_normalize_by_pathway()
  return()
  # 4
  step_log_transform(step_count,skip_past,reload_info,plot_steps,k_obs_groups,cluster_by_groups,
                     cluster_references,out_dir,output_format,png_res,useRaster)
  reaching_step(up_to_step, step_count)
  # 5
  #step_scale_expr_data(step_count,skip_past,reload_info,plot_steps,k_obs_groups,cluster_by_groups,
                        #cluster_references,out_dir,output_format,png_res,useRaster)
  #reaching_step(up_to_step, step_count)
  # 6
  #step_split_ref_data_into_groups(step_count,skip_past,reload_info, num_ref_groups,hclust_method)
  #reaching_step(up_to_step, step_count)
  # 7
  #step_compute_subclusters_random_trees(step_count,skip_past,reload_info,analysis_mode,tumor_subcluster_partition_method,
                                          #tumor_subcluster_pval,hclust_method,cluster_by_groups,plot_steps,k_obs_groups,
                                          #cluster_references,out_dir,output_format,png_res,useRaster)
 # reaching_step(up_to_step, step_count)
  # 8
  step_subtract_average_ref(step_count,skip_past,reload_info,ref_subtract_use_mean_bounds,plot_steps,
                            k_obs_groups,cluster_by_groups,cluster_references,out_dir,output_format,
                            png_res,useRaster)
  reaching_step(up_to_step, step_count)
  # 9
  step_apply_max_centered_expr_threshold(step_count,skip_past,reload_info,max_centered_threshold,threshold,plot_steps,
                                         k_obs_groups,cluster_by_groups,cluster_references,out_dir,output_format,
                                         png_res,useRaster)
  reaching_step(up_to_step, step_count)
  # 10
  step_smooth_by_chr(step_count,skip_past,reload_info,smooth_method,window_length,plot_steps,k_obs_groups,
                      cluster_by_groups,cluster_references,out_dir,output_format,png_res,useRaster)
  reaching_step(up_to_step, step_count)
  # 11
  step_re_centers_data(step_count,skip_past,reload_info,plot_steps,k_obs_groups,cluster_by_groups,
                        cluster_references,out_dir,output_format,png_res,useRaster)
  reaching_step(up_to_step, step_count)
  # 12
  step_remove_avg_after_smoothing(step_count,skip_past,reload_info,ref_subtract_use_mean_bounds,
                                  plot_steps,k_obs_groups,cluster_by_groups,cluster_references,
                                  out_dir,output_format,png_res,useRaster)
  reaching_step(up_to_step, step_count)
  # 13
  #step_remove_chr_ends(step_count,skip_past,reload_info,remove_genes_at_chr_ends,smooth_method,
                       #window_length,plot_steps,k_obs_groups,cluster_by_groups,
                       #cluster_references,out_dir,output_format,png_res,useRaster)
  #reaching_step(up_to_step, step_count)
  # 14
  step_invert_log_transform(step_count,skip_past,reload_info,plot_steps,k_obs_groups,cluster_by_groups,
                            cluster_references,out_dir,output_format,png_res,useRaster)
  reaching_step(up_to_step, step_count)
  # 15
  step_cluster_samples(step_count,skip_past,reload_info,analysis_mode,tumor_subcluster_partition_method,
                       tumor_subcluster_pval,hclust_method,cluster_by_groups,tumor_subcluster_partition_method,
                       plot_steps,k_obs_groups,cluster_by_groups,cluster_references,out_dir,output_format,
                       png_res,useRaster)
  reaching_step(up_to_step, step_count)
  # 16
  #step_remove_outliers(step_count,skip_past,reload_info,prune_outliers,outlier_method_bound,outlier_lower_bound,
                       #outlier_upper_bound,plot_steps,k_obs_groups,cluster_by_groups,cluster_references,out_dir,
                       #output_format,png_res,useRaster)
  #reaching_step(up_to_step, step_count)
  # 17
  step_run_hmm(step_count,skip_past,reload_info,analysis_mode,skip_hmm,HMM,HMM_type,HMM_transition_prob,
               HMM_i3_pval,HMM_i3_use_KS,tumor_subcluster_partition_method,hclust_method,hmm_resume_file_token,
               out_dir,hmm_center,HMM_report_by,no_plot,k_obs_groups,cluster_by_groups,cluster_references,
               output_format,hmm_state_range,png_res,useRaster)
  reaching_step(up_to_step, step_count)
  # 18
  step_run_bayesian_network(step_count,skip_past,reload_info,skip_hmm,HMM,BayesMaxPNormal,out_dir,no_plot,
                            hmm_resume_file_token,num_threads,plot_probabilities,diagnostics,HMM_type,
                            k_obs_groups,cluster_by_groups,reassignCNVs,cluster_references,output_format,
                            png_res,useRaster)
  reaching_step(up_to_step, step_count)
  # 19
  step_convert_states(step_count,skip_past,reload_info,skip_hmm,HMM,HMM_type,no_plot,k_obs_groups,
                      cluster_by_groups,cluster_references,out_dir,hmm_resume_file_token,
                      BayesMaxPNormal,output_format,png_res,useRaster)
  reaching_step(up_to_step, step_count)
  # 20
  #step_filter_DE_genes(step_count,skip_past,reload_info,mask_nonDE_genes,mask_nonDE_pval,test.use,
                       #require_DE_all_normals,plot_steps,k_obs_groups,cluster_by_groups,cluster_references,
                       #out_dir,output_format,png_res,useRaster)
  #reaching_step(up_to_step, step_count)
  # 21
  #step_denoise(step_count,skip_past,reload_info,denoise,noise_filter,noise_logistic,sd_amplifier,no_plot,
               #k_obs_groups,cluster_by_groups,cluster_references,out_dir,output_format,png_res,useRaster)
  
  return(infercnv_obj)
}

reaching_step <- function(){
  # Returns inferCNV obj at current step
  if (up_to_step == step_count) {
    flog.info("Reached up_to_step")
    return(infercnv_obj)
  }
}

get_incoming_data <- function(){
  #
  step_count = step_count + 1 # 1
  flog.info(sprintf("\n\n\tSTEP %d: incoming data\n", step_count))
  
  if (skip_past < step_count && save_rds) {
    saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
  }
}

remove_lowly_expr_genes <- function(){
  # Removes insufficiently expressed genes
  step_count = step_count + 1 # 2
  flog.info(sprintf("\n\n\tSTEP %02d: Removing lowly expressed genes\n", step_count))
  
  ## Remove genes that aren't sufficiently expressed, according to min mean count cutoff.
  ## Examines the original (non-log-transformed) data, gets mean for each gene, and removes genes
  ##  with mean values below cutoff.
  
  if (skip_past < step_count) {
    # }
    infercnv_obj <- require_above_min_mean_expr_cutoff(infercnv_obj, cutoff)
    
    ## require each gene to be present in a min number of cells for ref sets
    
    infercnv_obj <- require_above_min_cells_ref(infercnv_obj, min_cells_per_gene=min_cells_per_gene)
    
    if (save_rds) {
      saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
    }
  }
}

normalize_by_seq_depth <- function(){
  # Normalizes by sequencing depth
  step_count = step_count + 1 # 3
  flog.info(sprintf("\n\n\tSTEP %02d: normalization by sequencing depth\n", step_count))
  
  
  resume_file_token = ifelse( (HMM), paste0("HMM",HMM_type), "")
  
  if (skip_past < step_count) {
    infercnv_obj <- normalize_counts_by_seq_depth(infercnv_obj)
    
    if (HMM && HMM_type == 'i6') {
      ## add in the hidden spike needed by the HMM
      infercnv_obj <- .build_and_add_hspike(infercnv_obj, sim_method=sim_method, aggregate_normals=hspike_aggregate_normals)
      
      if (sim_foreground) {
        infercnv_obj <- .sim_foreground(infercnv_obj, sim_method=sim_method)
      }
    }
    
    if (save_rds) {
      saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
    }
  }
}

log_transform <- function(){
  # Transforms to log scale
  step_count = step_count + 1 # 4
  flog.info(sprintf("\n\n\tSTEP %02d: log transformation of data\n", step_count))
  
  if (skip_past < step_count) {
    infercnv_obj <- log2xplus1(infercnv_obj)
    
    if (save_rds) {
      saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
    }
    invisible(gc())
    
    ## Plot incremental steps.
    if (plot_steps){
      plot_cnv(infercnv_obj=infercnv_obj,
               k_obs_groups=k_obs_groups,
               cluster_by_groups=cluster_by_groups,
               cluster_references=cluster_references,
               out_dir=out_dir,
               title=sprintf("%02d_log_transformed_data",step_count),
               output_filename=sprintf("infercnv.%02d_log_transformed",step_count),
               output_format=output_format,
               write_expr_matrix=TRUE,
               png_res=png_res,
               useRaster=useRaster
      )
    }
  }
}

scale_expr_data <- function(){
  # Scales all expression data
  step_count = step_count + 1 # 5
  if (scale_data) {
    
    flog.info(sprintf("\n\n\tSTEP %02d: scaling all expression data\n", step_count))
    
    if (skip_past < step_count) {
      infercnv_obj <- scale_infercnv_expr(infercnv_obj)
      
      if (save_rds) {
        saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
      }
      invisible(gc())
      
      ## Plot incremental steps.
      if (plot_steps){
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 cluster_references=cluster_references,
                 out_dir=out_dir,
                 title=sprintf("%02d_scaled",step_count),
                 output_filename=sprintf("infercnv.%02d_scaled",step_count),
                 output_format=output_format,
                 write_expr_matrix=TRUE,
                 png_res=png_res,
                 useRaster=useRaster)
        
      }
    }
  }
}

split_ref_data_into_groups <- function(){
  # Split the reference data into groups if requested
  if (!is.null(num_ref_groups)) {
    
    if (! has_reference_cells(infercnv_obj)) {
      stop("Error, no reference cells defined. Cannot split them into groups as requested")
    }
    
    flog.info(sprintf("\n\n\tSTEP %02d: splitting reference data into %d clusters\n", step_count, num_ref_groups))
    
    if (skip_past < step_count) {
      infercnv_obj <- split_references(infercnv_obj,
                                       num_groups=num_ref_groups,
                                       hclust_method=hclust_method)
      if (save_rds) {
        saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
      }
    }
  }
}

compute_subclusters_random_trees <- function(){
  # Computes tumor subclusters via random trees
  # Only if splitting the reference data into groups is requested
  step_count = step_count + 1 # 7
  if (analysis_mode == 'subclusters' & tumor_subcluster_partition_method == 'random_trees') {
    
    flog.info(sprintf("\n\n\tSTEP %02d: computing tumor subclusters via %s\n", step_count, tumor_subcluster_partition_method))
    
    resume_file_token = paste0(resume_file_token, ".rand_trees")
    
    if (skip_past < step_count) {
      infercnv_obj <- define_signif_tumor_subclusters_via_random_smooothed_trees(infercnv_obj,
                                                                                 p_val=tumor_subcluster_pval,
                                                                                 hclust_method=hclust_method,
                                                                                 cluster_by_groups=cluster_by_groups)
      if (save_rds) {
        saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
      }
      invisible(gc())
      
      if (plot_steps) {
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 cluster_references=cluster_references,
                 out_dir=out_dir,
                 title=sprintf("%02d_tumor_subclusters.%s", step_count, tumor_subcluster_partition_method),
                 output_filename=sprintf("infercnv.%02d_tumor_subclusters.%s", step_count, tumor_subcluster_partition_method),
                 output_format=output_format,
                 write_expr_matrix=TRUE,
                 png_res=png_res,
                 useRaster=useRaster)
      }
    }
  }
}

subtract_average_ref <- function(){
  # Subtract average of reference data (before smoothing)
  # Since we're in log space, this now becomes log(fold_change)
  step_count = step_count + 1 # 8
  flog.info(sprintf("\n\n\tSTEP %02d: removing average of reference data (before smoothing)\n", step_count))
  
  if (skip_past < step_count) {
    infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj, inv_log=FALSE, use_bounds=ref_subtract_use_mean_bounds)
    
    if (save_rds) {
      saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
    }
    invisible(gc())
    
    if (plot_steps) {
      plot_cnv(infercnv_obj,
               k_obs_groups=k_obs_groups,
               cluster_by_groups=cluster_by_groups,
               cluster_references=cluster_references,
               out_dir=out_dir,
               title=sprintf("%02d_remove_average",step_count),
               output_filename=sprintf("infercnv.%02d_remove_average", step_count),
               output_format=output_format,
               write_expr_matrix=TRUE,
               png_res=png_res,
               useRaster=useRaster)
    }
  }
}

apply_max_centered_expr_threshold <- function(){
  # Applies maximum centered expression thresholds to data
  # Caps values between threshold and -threshold, retaining earlier center
  step_count = step_count + 1 # 9
  if (! is.na(max_centered_threshold)) {
    flog.info(sprintf("\n\n\tSTEP %02d: apply max centered expression threshold: %s\n", step_count, max_centered_threshold))
    
    if (skip_past < step_count) {
      threshold = max_centered_threshold
      if (is.character(max_centered_threshold) && max_centered_threshold == "auto") {
        threshold = mean(abs(get_average_bounds(infercnv_obj)))
        flog.info(sprintf("Setting max centered threshoolds via auto to: +- %g", threshold))
      }
      
      infercnv_obj <- apply_max_threshold_bounds(infercnv_obj, threshold=threshold)
      
      if (save_rds) {
        saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
      }
      invisible(gc())
      
      ## Plot incremental steps.
      if (plot_steps){
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 cluster_references=cluster_references,
                 out_dir=out_dir,
                 title=sprintf("%02d_apply_max_centered_expr_threshold",step_count),
                 output_filename=sprintf("infercnv.%02d_apply_max_centred_expr_threshold",step_count),
                 output_format=output_format,
                 write_expr_matrix=TRUE,
                 png_res=png_res,
                 useRaster=useRaster)
        
      }
    }
  }
}

smooth_by_chr <- function(){
  # For each cell, smoothes the data along chromosome with gene windows
  step_count = step_count + 1 # 10
  flog.info(sprintf("\n\n\tSTEP %02d: Smoothing data per cell by chromosome\n", step_count))
  
  if (skip_past < step_count) {
    
    if (smooth_method == 'runmeans') {
      
      infercnv_obj <- smooth_by_chromosome_runmeans(infercnv_obj, window_length)
    } else if (smooth_method == 'pyramidinal') {
      
      infercnv_obj <- smooth_by_chromosome(infercnv_obj, window_length=window_length, smooth_ends=TRUE)
    } else if (smooth_method == 'coordinates') {
      infercnv_obj <- smooth_by_chromosome_coordinates(infercnv_obj, window_length=window_length)
    } else {
      stop(sprintf("Error, don't recognize smoothing method: %s", smooth_method))
    }
    
    if (save_rds) {
      saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
    }
    invisible(gc())
    
    ## Plot incremental steps.
    if (plot_steps){
      
      plot_cnv(infercnv_obj,
               k_obs_groups=k_obs_groups,
               cluster_by_groups=cluster_by_groups,
               cluster_references=cluster_references,
               out_dir=out_dir,
               title=sprintf("%02d_smoothed_by_chr",step_count),
               output_filename=sprintf("infercnv.%02d_smoothed_by_chr", step_count),
               output_format=output_format,
               write_expr_matrix=TRUE,
               png_res=png_res,
               useRaster=useRaster)
    }
  }
}

re_centers_data <- function(){
  # Centers cells/observations after smoothing. 
  # Helps reduce the effect of complexity.
  step_count = step_count + 1 # 11
  flog.info(sprintf("\n\n\tSTEP %02d: re-centering data across chromosome after smoothing\n", step_count))
  
  if (skip_past < step_count) {
    infercnv_obj <- center_cell_expr_across_chromosome(infercnv_obj, method="median")
    
    if (save_rds) {
      saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
    }
    invisible(gc())
    
    ## Plot incremental steps.
    if (plot_steps) {
      
      plot_cnv(infercnv_obj,
               k_obs_groups=k_obs_groups,
               cluster_by_groups=cluster_by_groups,
               cluster_references=cluster_references,
               out_dir=out_dir,
               title=sprintf("%02d_centering_of_smoothed",step_count),
               output_filename=sprintf("infercnv.%02d_centering_of_smoothed", step_count),
               output_format=output_format,
               write_expr_matrix=TRUE,
               png_res=png_res,
               useRaster=useRaster)
      
    }
  }
}

remove_avg_after_smoothing <- function(){
  # Subtracts average reference (adjustment after smoothing)
  step_count = step_count + 1 # 12
  flog.info(sprintf("\n\n\tSTEP %02d: removing average of reference data (after smoothing)\n", step_count))
  
  if (skip_past < step_count) {
    infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj, inv_log=FALSE, use_bounds=ref_subtract_use_mean_bounds)
    
    if (save_rds) {
      saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
    }
    invisible(gc())
    
    if (plot_steps) {
      plot_cnv(infercnv_obj,
               k_obs_groups=k_obs_groups,
               cluster_by_groups=cluster_by_groups,
               cluster_references=cluster_references,
               out_dir=out_dir,
               title=sprintf("%02d_remove_average",step_count),
               output_filename=sprintf("infercnv.%02d_remove_average", step_count),
               output_format=output_format,
               write_expr_matrix=TRUE,
               png_res=png_res,
               useRaster=useRaster)
    }
  }
}

remove_ends <- function(){
  # If requested, removes genes at chromosome ends
  step_count = step_count + 1 # 13
  if (remove_genes_at_chr_ends == TRUE && smooth_method != 'coordinates') {
    
    flog.info(sprintf("\n\n\tSTEP %02d: removing genes at chr ends\n", step_count))
    
    if (skip_past < step_count) {
      infercnv_obj <- remove_genes_at_ends_of_chromosomes(infercnv_obj, window_length)
      
      if (save_rds) {
        saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
      }
      invisible(gc())
      
      ## Plot incremental steps.
      if (plot_steps){
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 cluster_references=cluster_references,
                 out_dir=out_dir,
                 title=sprintf("%02d_remove_genes_at_chr_ends",step_count),
                 output_filename=sprintf("infercnv.%02d_remove_genes_at_chr_ends",step_count),
                 output_format=output_format,
                 write_expr_matrix=TRUE,
                 png_res=png_res)
        
      }
    }
  }
}

invert_log_transform <- function(){
  # Inverts log transform  (converts from log(FC) to FC)
  step_count = step_count + 1 # 14
  flog.info(sprintf("\n\n\tSTEP %02d: invert log2(FC) to FC\n", step_count))
  
  if (skip_past < step_count) {
    
    infercnv_obj <- invert_log2(infercnv_obj)
    
    if (save_rds) {
      saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
    }
    invisible(gc())
    
    if (plot_steps) {
      plot_cnv(infercnv_obj,
               k_obs_groups=k_obs_groups,
               cluster_by_groups=cluster_by_groups,
               cluster_references=cluster_references,
               out_dir=out_dir,
               title=sprintf("%02d_invert_log_transform log(FC)->FC",step_count),
               output_filename=sprintf("infercnv.%02d_invert_log_FC",step_count),
               output_format=output_format,
               write_expr_matrix=TRUE,
               png_res=png_res,
               useRaster=useRaster)
      
    }
  }
}

cluster_samples <- function(){
  # Clusters samples if tumor subclusters not requested
  # Else computes tumor subclusters via a partition method other than random trees
  #
  if (analysis_mode == 'subclusters' & tumor_subcluster_partition_method != 'random_trees') {
    
    resume_file_token = paste0(resume_file_token, '.', tumor_subcluster_partition_method)
    
    flog.info(sprintf("\n\n\tSTEP %02d: computing tumor subclusters via %s\n", step_count, tumor_subcluster_partition_method))
    
    if (skip_past < step_count) {
      infercnv_obj <- define_signif_tumor_subclusters(infercnv_obj,
                                                      p_val=tumor_subcluster_pval,
                                                      hclust_method=hclust_method,
                                                      cluster_by_groups=cluster_by_groups,
                                                      partition_method=tumor_subcluster_partition_method)
      
      if (save_rds) {
        saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
      }
      invisible(gc())
      
      if (plot_steps) {
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 cluster_references=cluster_references,
                 out_dir=out_dir,
                 title=sprintf("%02d_tumor_subclusters",step_count),
                 output_filename=sprintf("infercnv.%02d_tumor_subclusters",step_count),
                 output_format=output_format,
                 write_expr_matrix=TRUE,
                 png_res=png_res,
                 useRaster=useRaster)
      }
      
    }
  }
  else if (analysis_mode != 'subclusters') {
    
    flog.info(sprintf("\n\n\tSTEP %02d: Clustering samples (not defining tumor subclusters)\n", step_count))
    
    if (skip_past < step_count) {
      
      
      infercnv_obj <- define_signif_tumor_subclusters(infercnv_obj,
                                                      p_val=tumor_subcluster_pval,
                                                      hclust_method=hclust_method,
                                                      cluster_by_groups=cluster_by_groups,
                                                      partition_method='none')
      
      if (save_rds) {
        saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
      }
    }
    
  }
  
  ## This is a milestone step and results should always be examined here.
  if (skip_past < step_count) {
    if (save_rds) {
      saveRDS(infercnv_obj, file=file.path(out_dir, "preliminary.infercnv_obj"))
    }
    
    invisible(gc())
    
    if (! (no_prelim_plot | no_plot) ) {
      
      prelim_heatmap_png = "infercnv.preliminary.png"
      
      if (! file.exists(file.path(out_dir, prelim_heatmap_png))) {
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 cluster_references=cluster_references,
                 out_dir=out_dir,
                 title="Preliminary infercnv (pre-noise filtering)",
                 output_filename="infercnv.preliminary", # png ext auto added
                 output_format=output_format,
                 write_expr_matrix=TRUE,
                 png_res=png_res,
                 useRaster=useRaster)
      }
    }
  }
}

remove_outliers <- function(){
  # Removes outliers for visualization (optional)
  step_count = step_count + 1 # 16
  if (prune_outliers) {
    
    ## ################################
    ## STEP: Remove outliers for viz
    
    flog.info(sprintf("\n\n\tSTEP %02d: Removing outliers\n", step_count))
    
    if (skip_past < step_count) {
      
      infercnv_obj = remove_outliers_norm(infercnv_obj,
                                          out_method=outlier_method_bound,
                                          lower_bound=outlier_lower_bound,
                                          upper_bound=outlier_upper_bound)
      
      if (save_rds) {
        saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
      }
      invisible(gc())
      
      ## Plot incremental steps.
      if (plot_steps) {
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 cluster_references=cluster_references,
                 out_dir=out_dir,
                 title=sprintf("%02d_removed_outliers",step_count),
                 output_filename=sprintf("infercnv.%02d_removed_outliers", step_count),
                 output_format=output_format,
                 write_expr_matrix=TRUE,
                 png_res=png_res,
                 useRaster=useRaster)
      }
    }
  }
}

run_hmm <- function(){
  # Computes HMM-based CNV prediction
  step_count = step_count + 1 # 17
  hmm_resume_file_token = paste0(resume_file_token, ".hmm_mode-", analysis_mode)
  if (skip_hmm < 1) {
    if (HMM) {
      flog.info(sprintf("\n\n\tSTEP %02d: HMM-based CNV prediction\n", step_count))
      
      if (HMM_type == 'i6') {
        hmm_center = 3
        hmm_state_range = c(0,6)
      } else {
        ## i3
        hmm_center = 2
        hmm_state_range = c(1,3)
      }
      
      if (analysis_mode == 'subclusters') {
        
        if (HMM_type == 'i6') {
          hmm.infercnv_obj <- predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj,
                                                                       t=HMM_transition_prob)
        } else if (HMM_type == 'i3') {
          hmm.infercnv_obj <- i3HMM_predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj,
                                                                             i3_p_val=HMM_i3_pval,
                                                                             t=HMM_transition_prob,
                                                                             use_KS=HMM_i3_use_KS)
        } else {
          stop("Error, not recognizing HMM_type")
        }
        
        if (tumor_subcluster_partition_method == 'random_trees') {
          ## need to redo the hierarchicial clustering, since the subcluster assignments dont always perfectly line up with the top-level dendrogram.
          hmm.infercnv_obj <- .redo_hierarchical_clustering(hmm.infercnv_obj, hclust_method=hclust_method)
        }
        
      } else if (analysis_mode == 'cells') {
        
        if (HMM_type == 'i6') {
          hmm.infercnv_obj <- predict_CNV_via_HMM_on_indiv_cells(infercnv_obj, t=HMM_transition_prob)
        } else if (HMM_type == 'i3') {
          hmm.infercnv_obj <- i3HMM_predict_CNV_via_HMM_on_indiv_cells(infercnv_obj,
                                                                       i3_p_val=HMM_i3_pval,
                                                                       t=HMM_transition_prob,
                                                                       use_KS=HMM_i3_use_KS)
        } else {
          stop("Error, not recognizing HMM_type")
        }
        
        
      } else {
        ## samples mode
        
        if (HMM_type == 'i6') {
          hmm.infercnv_obj <- predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_obj, t=HMM_transition_prob)
        } else if (HMM_type == 'i3') {
          hmm.infercnv_obj <- i3HMM_predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj,
                                                                             i3_p_val=HMM_i3_pval,
                                                                             t=HMM_transition_prob,
                                                                             use_KS=HMM_i3_use_KS
          )
        } else {
          stop("Error, not recognizing HMM_type")
        }
        
      }
      
      ## ##################################
      ## Note, HMM invercnv object is only leveraged here, but stored as file for future use:
      ## ##################################
      
      
      if (save_rds) {
        saveRDS(hmm.infercnv_obj, reload_info$expected_file_names[[step_count]])
      }
      
      ## report predicted cnv regions:
      generate_cnv_region_reports(hmm.infercnv_obj,
                                  output_filename_prefix=sprintf("%02d_HMM_pred%s", step_count, hmm_resume_file_token),
                                  out_dir=out_dir,
                                  ignore_neutral_state=hmm_center,
                                  by=HMM_report_by)
      
      invisible(gc())
      
      if (! no_plot) {
        
        ## Plot HMM pred img
        plot_cnv(infercnv_obj=hmm.infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 cluster_references=cluster_references,
                 out_dir=out_dir,
                 title=sprintf("%02d_HMM_preds",step_count),
                 output_filename=sprintf("infercnv.%02d_HMM_pred%s", step_count, hmm_resume_file_token),
                 output_format=output_format,
                 write_expr_matrix=TRUE,
                 x.center=hmm_center,
                 x.range=hmm_state_range,
                 png_res=png_res,
                 useRaster=useRaster
        )
      }
    }
  }
  
}

run_bayesian_network <- function(){
  # Runs Bayesian Network Model on HMM predicted CNV's
  step_count = step_count + 1 # 18
  if (skip_hmm < 2) {
    if (HMM == TRUE && BayesMaxPNormal > 0 && length(unique(apply(hmm.infercnv_obj@expr.data,2,unique))) != 1 ) {
      flog.info(sprintf("\n\n\tSTEP %02d: Run Bayesian Network Model on HMM predicted CNV's\n", step_count))
      
      ## the MCMC  object
      
      mcmc_obj <- infercnv::inferCNVBayesNet( infercnv_obj     = infercnv_obj,
                                              HMM_states        = hmm.infercnv_obj@expr.data,
                                              file_dir          = out_dir,
                                              no_plot           = no_plot,
                                              postMcmcMethod    = "removeCNV",
                                              out_dir           = file.path(out_dir, sprintf("BayesNetOutput.%s", hmm_resume_file_token)),
                                              resume_file_token = hmm_resume_file_token,
                                              quietly           = TRUE,
                                              CORES             = num_threads,
                                              plotingProbs      = plot_probabilities,
                                              diagnostics       = diagnostics,
                                              HMM_type          = HMM_type,
                                              k_obs_groups      = k_obs_groups,
                                              cluster_by_groups = cluster_by_groups,
                                              reassignCNVs      = reassignCNVs)
      
      mcmc_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj",
                                                 step_count, hmm_resume_file_token))
      
      if (save_rds) {
        saveRDS(mcmc_obj, file=mcmc_obj_file)
      }
      
      ## Filter CNV's by posterior Probabilities
      mcmc_obj_hmm_states_list <- infercnv::filterHighPNormals( MCMC_inferCNV_obj = mcmc_obj,
                                                                HMM_states = hmm.infercnv_obj@expr.data, 
                                                                BayesMaxPNormal   = BayesMaxPNormal)
      
      hmm_states_highPnormCNVsRemoved.matrix <- mcmc_obj_hmm_states_list[[2]]
      
      # replace states
      hmm.infercnv_obj@expr.data <- hmm_states_highPnormCNVsRemoved.matrix
      
      ## Save the MCMC inferCNV object
      if (save_rds) {
        saveRDS(hmm.infercnv_obj, reload_info$expected_file_names[[step_count]])
      }
      invisible(gc())
      
      if (! no_plot) {
        ## Plot HMM pred img after cnv removal
        plot_cnv(infercnv_obj=hmm.infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 cluster_references=cluster_references,
                 out_dir=out_dir,
                 title=sprintf("%02d_HMM_preds_Bayes_Net",step_count),
                 output_filename=sprintf("infercnv.%02d_HMM_pred.Bayes_Net.Pnorm_%g",step_count, BayesMaxPNormal),
                 output_format=output_format,
                 write_expr_matrix=TRUE,
                 x.center=3,
                 x.range=c(0,6),
                 png_res=png_res,
                 useRaster=useRaster
        )
      }
      ## write the adjusted CNV report files
      ## report predicted cnv regions:
      adjust_genes_regions_report(mcmc_obj_hmm_states_list[[1]],
                                  # input_filename_prefix=sprintf("%02d_HMM_preds", (step_count-1)),
                                  input_filename_prefix=sprintf("%02d_HMM_pred%s", (step_count-1), hmm_resume_file_token),
                                  output_filename_prefix=sprintf("HMM_CNV_predictions.%s.Pnorm_%g", hmm_resume_file_token, BayesMaxPNormal),
                                  out_dir=out_dir)
    }
  }
}

convert_states <- function(){
  # Converts HMM-based CNV states to representative intensity values
  step_count = step_count + 1 # 19
  if (skip_hmm < 3) {
    if (HMM) {
      flog.info(sprintf("\n\n\tSTEP %02d: Converting HMM-based CNV states to repr expr vals\n", step_count))
      
      if (HMM_type == 'i6') {
        hmm.infercnv_obj <- assign_HMM_states_to_proxy_expr_vals(hmm.infercnv_obj)
      } else if (HMM_type == 'i3') {
        hmm.infercnv_obj <- i3HMM_assign_HMM_states_to_proxy_expr_vals(hmm.infercnv_obj)
      }
      
      if (save_rds) {
        saveRDS(hmm.infercnv_obj, reload_info$expected_file_names[[step_count]])
      }
      invisible(gc())
      
      ## Plot HMM pred img
      if (! no_plot) {
        plot_cnv(infercnv_obj=hmm.infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 cluster_references=cluster_references,
                 out_dir=out_dir,
                 title=sprintf("%02d_HMM_preds.repr_intensities",step_count),
                 output_filename=sprintf("infercnv.%02d_HMM_pred%s.Pnorm_%g.repr_intensities", step_count, hmm_resume_file_token, BayesMaxPNormal),
                 output_format=output_format,
                 write_expr_matrix=TRUE,
                 x.center=1,
                 x.range=c(-1,3),
                 png_res=png_res,
                 useRaster=useRaster
        )
      }
    }
  }
}


library(optimx)         # You need to have some version of optimx available
# as it is a BioGeoBEARS dependency; however, if you
# don't want to use optimx, and use optim() (from R core) 
# you can set:
# BioGeoBEARS_run_object$use_optimx = FALSE
# ...everything should work either way -- NJM 2014-01-08
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)
source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
# slight speedup hopefully


######this script runs all 6 models in BioGeoBEARS and outputs a table with AICs sortes
#run_BioGeoBEARS_models('./biogeo_test/Otomys.tree','./biogeo_test/Otomys_geogdata.txt','./biogeo_test/','Otomys')
run_BioGeoBEARS_models<-function(treefile,geographyfile,path,name){
  
  #extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
  #list.files(extdata_dir)
  setwd(path)
  tr = read.tree(treefile)
  geogfn<-geographyfile
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  max_range_size = 2
  #######################################################
  #######################################################
  # DEC AND DEC+J ANALYSIS
  #######################################################
  #######################################################
  # NOTE: The BioGeoBEARS "DEC" model is identical with 
  # the Lagrange DEC model, and should return identical
  # ML estimates of parameters, and the same 
  # log-likelihoods, for the same datasets.
  #
  # Ancestral state probabilities at nodes will be slightly 
  # different, since BioGeoBEARS is reporting the 
  # ancestral state probabilities under the global ML
  # model, and Lagrange is reporting ancestral state
  # probabilities after re-optimizing the likelihood
  # after fixing the state at each node. These will 
  # be similar, but not identical. See Matzke (2014),
  # Systematic Biology, for discussion.
  #
  # Also see Matzke (2014) for presentation of the 
  # DEC+J model.
  #######################################################
  #######################################################
  
  #######################################################
  #######################################################
  
  #######################################################
  # Run DEC
  #######################################################
  
  # Intitialize a default model (DEC model)
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  
  # Give BioGeoBEARS the location of the phylogeny Newick file
  BioGeoBEARS_run_object$trfn = treefile
  
  # Give BioGeoBEARS the location of the geography text file
  BioGeoBEARS_run_object$geogfn = geographyfile
  
  # Input the maximum range size
  BioGeoBEARS_run_object$max_range_size = max_range_size
  
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  # 1. Here, un-comment ONLY the files you want to use.
  # 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
  # 3. For example files see (a) extdata_dir, 
  #  or (b) http://phylo.wikidot.com/biogeobears#files
  #  and BioGeoBEARS Google Group posts for further hints)
  #
  # Uncomment files you wish to use in time-stratified analyses:
  #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
  BioGeoBEARS_run_object$num_cores_to_use = 1
  # (use more cores to speed it up; this requires
  # library(parallel) and/or library(snow). The package "parallel" 
  # is now default on Macs in R 3.0+, but apparently still 
  # has to be typed on some Windows machines. Note: apparently 
  # parallel works on Mac command-line R, but not R.app.
  # BioGeoBEARS checks for this and resets to 1
  # core with R.app)
  
  # Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
  # I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
  # but the results are imprecise and so I haven't explored it further.
  # In a Bayesian analysis, it might work OK, but the ML point estimates are
  # not identical.
  # Also, I have not implemented all functions to work with force_sparse=TRUE.
  # Volunteers are welcome to work on it!!
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up DEC model
  # (nothing to do; defaults)
  
  # Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
  BioGeoBEARS_run_object
  
  # This contains the model object
  BioGeoBEARS_run_object$BioGeoBEARS_model_object
  
  # This table contains the parameters of the model 
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
  
  # Run this to check inputs. Read the error messages if you get them!
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # For a slow analysis, run once, then set runslow=FALSE to just 
  # load the saved result.
  runslow = TRUE
  resfn = paste(name,"_DEC.Rdata",sep='')
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    resDEC = res
  } else {
    # Loads to "res"
    load(resfn)
    resDEC = res
  }
  
  #######################################################
  # Run DEC+J
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = treefile
  BioGeoBEARS_run_object$geogfn = geographyfile
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
  BioGeoBEARS_run_object$num_cores_to_use = 1
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up DEC+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDEC$outputs@params_table["d","est"]
  estart = resDEC$outputs@params_table["e","est"]
  jstart = 0.0001
  
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  # Add j as a free parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  resfn = paste(name,"_DEC+J.Rdata",sep='')
  runslow = TRUE
  if (runslow)
  {
    #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
    
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    
    resDECj = res
  } else {
    # Loads to "res"
    load(resfn)
    resDECj = res
  }
  
  ##########################################################
  #### PDF plots
  ##########################################################
  ###pdffn = "Otomys_DEC_vs_DEC+J_M0_unconstrained_v1.pdf"
  ###pdf(pdffn, width=6, height=6)
  ###
  ##########################################################
  #### Plot ancestral states - DEC
  ##########################################################
  ###analysis_titletxt ="BioGeoBEARS DEC on Otomys M0_unconstrained"
  ###
  #### Setup
  ###results_object = resDEC
  ###scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  ###
  #### States
  ###res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  ###
  #### Pie chart
  ###plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  ###
  ##########################################################
  #### Plot ancestral states - DECJ
  ##########################################################
  ###analysis_titletxt ="BioGeoBEARS DEC+J on Otomys M0_unconstrained"
  ###
  #### Setup
  ###results_object = resDECj
  ###scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  ###
  #### States
  ###res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  ###
  #### Pie chart
  ###plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  ###
  ###dev.off()  # Turn off PDF
  ###cmdstr = paste("open ", pdffn, sep="")
  ###system(cmdstr) # Plot it
  
  #######################################################
  #######################################################
  # DIVALIKE AND DIVALIKE+J ANALYSIS
  #######################################################
  #######################################################
  # NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with 
  # Ronquist (1997)'s parsimony DIVA. It is a likelihood
  # interpretation of DIVA, constructed by modelling DIVA's
  # processes the way DEC does, but only allowing the 
  # processes DIVA allows (widespread vicariance: yes; subset
  # sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
  #
  # DIVALIKE is a likelihood interpretation of parsimony
  # DIVA, and it is "like DIVA" -- similar to, but not
  # identical to, parsimony DIVA.
  #
  # I thus now call the model "DIVALIKE", and you should also. ;-)
  #######################################################
  #######################################################
  
  #######################################################
  # Run DIVALIKE
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = treefile
  BioGeoBEARS_run_object$geogfn = geographyfile
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
  BioGeoBEARS_run_object$num_cores_to_use = 1
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up DIVALIKE model
  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
  
  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
  
  # No jump dispersal/founder-event speciation
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  runslow = TRUE
  resfn = paste(name,"_DIVALIKE_M0_unconstrained_v1.Rdata",sep='')
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    resDIVALIKE = res
  } else {
    # Loads to "res"
    load(resfn)
    resDIVALIKE = res
  }
  
  #######################################################
  # Run DIVALIKE+J
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = treefile
  BioGeoBEARS_run_object$geogfn = geographyfile
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
  BioGeoBEARS_run_object$num_cores_to_use = 1
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up DIVALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDIVALIKE$outputs@params_table["d","est"]
  estart = resDIVALIKE$outputs@params_table["e","est"]
  jstart = 0.0001
  
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
  
  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
  
  # Add jump dispersal/founder-event speciation
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  resfn = paste(name,"_DIVALIKE+J_M0_unconstrained_v1.Rdata",sep='')
  runslow = TRUE
  if (runslow)
  {
    #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
    
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    
    resDIVALIKEj = res
  } else {
    # Loads to "res"
    load(resfn)
    resDIVALIKEj = res
  }
  
  ###pdffn = "Otomys_DIVALIKE_vs_DIVALIKE+J_M0_unconstrained_v1.pdf"
  ###pdf(pdffn, width=6, height=6)
  ###
  ##########################################################
  #### Plot ancestral states - DIVALIKE
  ##########################################################
  ###analysis_titletxt ="BioGeoBEARS DIVALIKE on Otomys M0_unconstrained"
  ###
  #### Setup
  ###results_object = resDIVALIKE
  ###scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  ###
  #### States
  ###res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  ###
  #### Pie chart
  ###plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  ###
  ##########################################################
  #### Plot ancestral states - DIVALIKE+J
  ##########################################################
  ###analysis_titletxt ="BioGeoBEARS DIVALIKE+J on Otomys M0_unconstrained"
  ###
  #### Setup
  ###results_object = resDIVALIKEj
  ###scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  ###
  #### States
  ###res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  ###
  #### Pie chart
  ###plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  ###
  ###dev.off()
  ###cmdstr = paste("open ", pdffn, sep="")
  ###system(cmdstr)
  
  #######################################################
  #######################################################
  # BAYAREALIKE AND BAYAREALIKE+J ANALYSIS
  #######################################################
  #######################################################
  # NOTE: As with DIVA, the BioGeoBEARS BayArea-like model is 
  # not identical with the full Bayesian model implemented 
  # in the "BayArea" program of Landis et al. (2013). 
  #
  # Instead, this is a simplified likelihood interpretation
  # of the model.  Basically, in BayArea and BioGeoBEARS-BAYAREALIKE, 
  # "d" and "e" work like they do in the DEC model of Lagrange 
  # (and BioGeoBEARS), and then BayArea's cladogenesis assumption
  # (which is that nothing in particular happens at cladogenesis) is 
  # replicated by BioGeoBEARS.
  #
  # This leaves out 3 important things that are in BayArea:
  # 1. Distance dependence (you can add this with a distances 
  #    matrix + the "x" parameter in BioGeoBEARS, however)
  # 2. A correction for disallowing "e" events that drive
  #    a species extinct (a null geographic range)
  # 3. The neat Bayesian sampling of histories, which allows
  #    analyses on large numbers of areas.
  #
  # The main purpose of having a "BAYAREALIKE" model is 
  # to test the importance of the cladogenesis model on 
  # particular datasets. Does it help or hurt the data 
  # likelihood if there is no special cladogenesis process?
  # 
  # BAYAREALIKE is a likelihood interpretation of BayArea,
  # and it is "like BayArea" -- similar to, but not
  # identical to, Bayesian BayArea.
  # I thus now call the model "BAYAREALIKE", and you should also. ;-)
  #######################################################
  #######################################################
  
  #######################################################
  # Run BAYAREALIKE
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = treefile
  BioGeoBEARS_run_object$geogfn = geographyfile
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
  BioGeoBEARS_run_object$num_cores_to_use = 1
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up BAYAREALIKE model
  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  # No jump dispersal/founder-event speciation
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
  
  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
  
  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  
  # Check the inputs
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  runslow = TRUE
  resfn = paste(name,"_BAYAREALIKE_M0_unconstrained_v1.Rdata",sep='')
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    resBAYAREALIKE = res
  } else {
    # Loads to "res"
    load(resfn)
    resBAYAREALIKE = res
  }
  
  #######################################################
  # Run BAYAREALIKE+J
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = treefile
  BioGeoBEARS_run_object$geogfn = geographyfile
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE
  BioGeoBEARS_run_object$num_cores_to_use = 1
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up BAYAREALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resBAYAREALIKE$outputs@params_table["d","est"]
  estart = resBAYAREALIKE$outputs@params_table["e","est"]
  jstart = 0.0001
  
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
  
  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
  
  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  
  # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
  # machines. I can't replicate this on my Mac machines, but it is almost certainly
  # just some precision under-run issue, when optim/optimx tries some parameter value 
  # just below zero.  The "min" and "max" options on each parameter are supposed to
  # prevent this, but apparently optim/optimx sometimes go slightly beyond 
  # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
  # slightly for each parameter:
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  resfn = paste(name,"_BAYAREALIKE+J_M0_unconstrained_v1.Rdata",sep='')
  runslow = TRUE
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    
    resBAYAREALIKEj = res
  } else {
    # Loads to "res"
    load(resfn)
    resBAYAREALIKEj = res
  }
  
  ###pdffn = "Otomys_BAYAREALIKE_vs_BAYAREALIKE+J_M0_unconstrained_v1.pdf"
  ###pdf(pdffn, width=6, height=6)
  ###
  ##########################################################
  #### Plot ancestral states - BAYAREALIKE
  ##########################################################
  ###analysis_titletxt ="BioGeoBEARS BAYAREALIKE on Otomys M0_unconstrained"
  ###
  #### Setup
  ###results_object = resBAYAREALIKE
  ###scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  ###
  #### States
  ###res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  ###
  #### Pie chart
  ###plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  ###
  ##########################################################
  #### Plot ancestral states - BAYAREALIKE+J
  ##########################################################
  ###analysis_titletxt ="BioGeoBEARS BAYAREALIKE+J on Otomys M0_unconstrained"
  ###
  #### Setup
  ###results_object = resBAYAREALIKEj
  ###scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  ###
  #### States
  ###res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  ###
  #### Pie chart
  ###plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  ###
  ###dev.off()
  ###cmdstr = paste("open ", pdffn, sep="")
  ###system(cmdstr)
  ###
  ############################################################################
  ############################################################################
  #########################################################################
  #########################################################################
  # 
  # CALCULATE SUMMARY STATISTICS TO COMPARE
  # DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
  # 
  #########################################################################
  #########################################################################
  #########################################################################
  #########################################################################
  
  #########################################################################
  #########################################################################
  # REQUIRED READING:
  #
  # Practical advice / notes / basic principles on statistical model 
  #    comparison in general, and in BioGeoBEARS:
  # http://phylo.wikidot.com/advice-on-statistical-model-comparison-in-biogeobears
  #########################################################################
  #########################################################################
  
  # Set up empty tables to hold the statistical results
  restable = NULL
  teststable = NULL
  
  #######################################################
  # Statistics -- DEC vs. DEC+J
  #######################################################
  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)
  
  numparams1 = 3
  numparams2 = 2
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  stats
  
  # DEC, null model for Likelihood Ratio Test (LRT)
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # DEC+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  # The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
  # confer the same likelihood on the data. See: Brian O'Meara's webpage:
  # http://www.brianomeara.info/tutorials/aic
  # ...for an intro to LRT, AIC, and AICc
  
  rbind(res2, res1)
  tmp_tests = conditional_format_table(stats)
  
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
  
  #######################################################
  # Statistics -- DIVALIKE vs. DIVALIKE+J
  #######################################################
  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)
  
  numparams1 = 3
  numparams2 = 2
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  stats
  
  # DIVALIKE, null model for Likelihood Ratio Test (LRT)
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  rbind(res2, res1)
  conditional_format_table(stats)
  
  tmp_tests = conditional_format_table(stats)
  
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
  
  #######################################################
  # Statistics -- BAYAREALIKE vs. BAYAREALIKE+J
  #######################################################
  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)
  
  numparams1 = 3
  numparams2 = 2
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  stats
  
  # BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  rbind(res2, res1)
  conditional_format_table(stats)
  
  tmp_tests = conditional_format_table(stats)
  
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
  
  #########################################################################
  # RESULTS: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
  #########################################################################
  teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
  teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
  row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
  
  # Look at the results!!
  restable
  teststable
  
  #######################################################
  # Save the results tables for later -- check for e.g.
  # convergence issues
  #######################################################
  
  # Loads to "restable"
  save(restable, file=paste(name,"_restable_v1.Rdata",sep=''))
  load(file=paste(name,"_restable_v1.Rdata",sep=''))
  
  # Loads to "teststable"
  save(teststable, file=paste(name,"_restable_v1.Rdata",sep=''))
  load(file=paste(name,"_restable_v1.Rdata",sep=''))
  
  # Also save to text files
  write.table(restable, file=paste(name,"_restable.txt",sep=''), quote=FALSE, sep="\t")
  write.table(unlist_df(teststable), file=paste(name,"_teststable.txt",sep=''), quote=FALSE, sep="\t")
  
  #######################################################
  # Model weights of all six models
  #######################################################
  restable2 = restable
  
  # With AICs:
  AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
  restable = cbind(restable, AICtable)
  restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
  restable_AIC_rellike
  
  # With AICcs -- factors in sample size
  samplesize = length(tr$tip.label)
  AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
  restable2 = cbind(restable2, AICtable)
  restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AIC")
  restable_AICc_rellike
  free_params = row.names(resDECj$output@params_table[resDECj$output@params_table$type=="free",])
  names(restable_AICc_rellike) = c("LnL", "numparams", free_params, "AICc", "AICc_wt")
  
  # Also save to text files
  write.table(restable_AIC_rellike, file=paste(name,"_restable_AIC_rellike.txt",sep=''), quote=FALSE, sep="\t")
  write.table(restable_AICc_rellike, file=paste(name,"_restable_AICc_rellike.txt",sep=''), quote=FALSE, sep="\t")
  
  # Save with nice conditional formatting
  write.table(conditional_format_table(restable_AIC_rellike), file=paste(name,"_restable_AIC_rellike_formatted.txt",sep=''), quote=FALSE, sep="\t")
  write.table(conditional_format_table(restable_AICc_rellike), file=paste(name,"_restable_AICc_rellike_formatted.txt",sep=''), quote=FALSE, sep="\t")
  
}
##setwd('../')
##resDEC<-run_BioGeoBEARS_selectedmodel_BSM('./Otomys.tree','./Otomys_geogdata.txt','./biogeo_test/','Otomys','DEC')
##setwd('../')
##resDECJ<-run_BioGeoBEARS_selectedmodel_BSM('./Otomys.tree','./Otomys_geogdata.txt','./biogeo_test/','Otomys','DEC+J')
##setwd('../')
##resDIVALIKE<-run_BioGeoBEARS_selectedmodel_BSM('./Otomys.tree','./Otomys_geogdata.txt','./biogeo_test/','Otomys','DIVALIKE')
##setwd('../')
##resDIVALIKEJ<-run_BioGeoBEARS_selectedmodel_BSM('./Otomys.tree','./Otomys_geogdata.txt','./biogeo_test/','Otomys','DIVALIKE+J')
##setwd('../')
##resBAYAREALIKE<-run_BioGeoBEARS_selectedmodel_BSM('./Otomys.tree','./Otomys_geogdata.txt','./biogeo_test/','Otomys','BAYAREALIKE')
##setwd('../')
##resBAYAREALIKEJ<-run_BioGeoBEARS_selectedmodel_BSM('./Otomys.tree','./Otomys_geogdata.txt','./biogeo_test/','Otomys','BAYAREALIKE+J')

########

#ptm <- proc.time()
#afro.BSM<-run_BioGeoBEARS_selectedmodel_BSM('./mammals_Rolland_terrestrial_IUCN_afrotropical.tree','./afrotropical_geographyfile.txt','./biogeo_test/','afrotropical','BAYAREALIKE+J')
#save(afro.BSM,file='./afro.BSM.Rsave')
#cat(proc.time() - ptm,'\n',file='/Users/javier/Dropbox/Work_in_progress/hotspots/hotspots_vertebrates/timer_file.txt')
#setwd('../')
#ptm <- proc.time()
#austral.BSM<-run_BioGeoBEARS_selectedmodel_BSM('./mammals_Rolland_terrestrial_IUCN_austral.tree','./austral_geographyfile.txt','./biogeo_test/','austral','DEC+J')
#save(austral.BSM,file='./austral.BSM.Rsave')
#cat(proc.time() - ptm,'\n',file='/Users/javier/Dropbox/Work_in_progress/hotspots/hotspots_vertebrates/timer_file.txt')
#setwd('../')
#ptm <- proc.time()
#indo.BSM<-run_BioGeoBEARS_selectedmodel_BSM('./mammals_Rolland_terrestrial_IUCN_indo.tree','./indo_geographyfile.txt','./biogeo_test/','indo','BAYAREALIKE+J')
#save(indo.BSM,file='./indo.BSM.Rsave')
#cat(proc.time() - ptm,'\n',file='/Users/javier/Dropbox/Work_in_progress/hotspots/hotspots_vertebrates/timer_file.txt')
#setwd('../')
#ptm <- proc.time()
#nearctic.BSM<-run_BioGeoBEARS_selectedmodel_BSM('./mammals_Rolland_terrestrial_IUCN_nearctic.tree','./nearctic_geographyfile.txt','./biogeo_test/','nearctic','DEC+J')
#save(nearctic.BSM,file='./nearctic.BSM.Rsave')
#cat(proc.time() - ptm,'\n',file='/Users/javier/Dropbox/Work_in_progress/hotspots/hotspots_vertebrates/timer_file.txt')


#afro.BSM$all_dispersals_counts_fromto_means
#afro.BSM$all_dispersals_counts_fromto_sds

#austral.BSM$all_dispersals_counts_fromto_means
#austral.BSM$all_dispersals_counts_fromto_sds

#indo.BSM$all_dispersals_counts_fromto_means
#indo.BSM$all_dispersals_counts_fromto_sds

#nearctic.BSM$all_dispersals_counts_fromto_means
#nearctic.BSM$all_dispersals_counts_fromto_sds


run_BioGeoBEARS_selectedmodel_BSM<-function(treefile,geographyfile,path,name,model_name){
  #extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
  #list.files(extdata_dir)
  setwd(path)
  tr = read.tree(treefile)
  geogfn<-geographyfile
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  max_range_size = 2
  #####
  if (model_name=='DEC'){
    #######################################################
    #######################################################
    # DEC AND DEC+J ANALYSIS
    #######################################################
    #######################################################
    # NOTE: The BioGeoBEARS "DEC" model is identical with 
    # the Lagrange DEC model, and should return identical
    # ML estimates of parameters, and the same 
    # log-likelihoods, for the same datasets.
    #
    # Ancestral state probabilities at nodes will be slightly 
    # different, since BioGeoBEARS is reporting the 
    # ancestral state probabilities under the global ML
    # model, and Lagrange is reporting ancestral state
    # probabilities after re-optimizing the likelihood
    # after fixing the state at each node. These will 
    # be similar, but not identical. See Matzke (2014),
    # Systematic Biology, for discussion.
    #
    # Also see Matzke (2014) for presentation of the 
    # DEC+J model.
    #######################################################
    #######################################################
    
    #######################################################
    #######################################################
    
    #######################################################
    # Run DEC
    #######################################################
    
    # Intitialize a default model (DEC model)
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    
    # Give BioGeoBEARS the location of the phylogeny Newick file
    BioGeoBEARS_run_object$trfn = treefile
    
    # Give BioGeoBEARS the location of the geography text file
    BioGeoBEARS_run_object$geogfn = geographyfile
    
    # Input the maximum range size
    BioGeoBEARS_run_object$max_range_size = max_range_size
    
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    # 1. Here, un-comment ONLY the files you want to use.
    # 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
    # 3. For example files see (a) extdata_dir, 
    #  or (b) http://phylo.wikidot.com/biogeobears#files
    #  and BioGeoBEARS Google Group posts for further hints)
    #
    # Uncomment files you wish to use in time-stratified analyses:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    # (use more cores to speed it up; this requires
    # library(parallel) and/or library(snow). The package "parallel" 
    # is now default on Macs in R 3.0+, but apparently still 
    # has to be typed on some Windows machines. Note: apparently 
    # parallel works on Mac command-line R, but not R.app.
    # BioGeoBEARS checks for this and resets to 1
    # core with R.app)
    
    # Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
    # I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
    # but the results are imprecise and so I haven't explored it further.
    # In a Bayesian analysis, it might work OK, but the ML point estimates are
    # not identical.
    # Also, I have not implemented all functions to work with force_sparse=TRUE.
    # Volunteers are welcome to work on it!!
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DEC model
    # (nothing to do; defaults)
    
    # Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
    BioGeoBEARS_run_object
    
    # This contains the model object
    BioGeoBEARS_run_object$BioGeoBEARS_model_object
    
    # This table contains the parameters of the model 
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
    
    # Run this to check inputs. Read the error messages if you get them!
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # For a slow analysis, run once, then set runslow=FALSE to just 
    # load the saved result.
    runslow = TRUE
    resfn = paste(name,"_DEC.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resDEC = res
    } else {
      # Loads to "res"
      load(resfn)
      resDEC = res
    }
    
  } else if (model_name=='DEC+J'){
    #######################################################
    # Run DEC+J
    #######################################################
    # Intitialize a default model (DEC model)
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    
    # Give BioGeoBEARS the location of the phylogeny Newick file
    BioGeoBEARS_run_object$trfn = treefile
    
    # Give BioGeoBEARS the location of the geography text file
    BioGeoBEARS_run_object$geogfn = geographyfile
    
    # Input the maximum range size
    BioGeoBEARS_run_object$max_range_size = max_range_size
    
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    # 1. Here, un-comment ONLY the files you want to use.
    # 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
    # 3. For example files see (a) extdata_dir, 
    #  or (b) http://phylo.wikidot.com/biogeobears#files
    #  and BioGeoBEARS Google Group posts for further hints)
    #
    # Uncomment files you wish to use in time-stratified analyses:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    # (use more cores to speed it up; this requires
    # library(parallel) and/or library(snow). The package "parallel" 
    # is now default on Macs in R 3.0+, but apparently still 
    # has to be typed on some Windows machines. Note: apparently 
    # parallel works on Mac command-line R, but not R.app.
    # BioGeoBEARS checks for this and resets to 1
    # core with R.app)
    
    # Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
    # I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
    # but the results are imprecise and so I haven't explored it further.
    # In a Bayesian analysis, it might work OK, but the ML point estimates are
    # not identical.
    # Also, I have not implemented all functions to work with force_sparse=TRUE.
    # Volunteers are welcome to work on it!!
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DEC model
    # (nothing to do; defaults)
    
    # Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
    BioGeoBEARS_run_object
    
    # This contains the model object
    BioGeoBEARS_run_object$BioGeoBEARS_model_object
    
    # This table contains the parameters of the model 
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
    
    # Run this to check inputs. Read the error messages if you get them!
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # For a slow analysis, run once, then set runslow=FALSE to just 
    # load the saved result.
    runslow = TRUE
    resfn = paste(name,"_DEC.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resDEC = res
    } else {
      # Loads to "res"
      load(resfn)
      resDEC = res
    }
    
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DEC+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resDEC$outputs@params_table["d","est"]
    estart = resDEC$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    
    # Add j as a free parameter
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    resfn = paste(name,"_DEC+J.Rdata",sep='')
    runslow = TRUE
    if (runslow)
    {
      #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
      
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resDECj = res
    } else {
      # Loads to "res"
      load(resfn)
      resDECj = res
    }
    
  } else if (model_name=='DIVALIKE'){
    
    #######################################################
    #######################################################
    # DIVALIKE AND DIVALIKE+J ANALYSIS
    #######################################################
    #######################################################
    # NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with 
    # Ronquist (1997)'s parsimony DIVA. It is a likelihood
    # interpretation of DIVA, constructed by modelling DIVA's
    # processes the way DEC does, but only allowing the 
    # processes DIVA allows (widespread vicariance: yes; subset
    # sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
    #
    # DIVALIKE is a likelihood interpretation of parsimony
    # DIVA, and it is "like DIVA" -- similar to, but not
    # identical to, parsimony DIVA.
    #
    # I thus now call the model "DIVALIKE", and you should also. ;-)
    #######################################################
    #######################################################
    
    #######################################################
    # Run DIVALIKE
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DIVALIKE model
    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    # No jump dispersal/founder-event speciation
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = paste(name,"_DIVALIKE_M0_unconstrained_v1.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resDIVALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKE = res
    }
    
  } else if(model_name=='DIVALIKE+J'){
    #######################################################
    # Run DIVALIKE
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DIVALIKE model
    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    # No jump dispersal/founder-event speciation
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = paste(name,"_DIVALIKE_M0_unconstrained_v1.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resDIVALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKE = res
    }
      
    
    #######################################################
    # Run DIVALIKE+J
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DIVALIKE+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resDIVALIKE$outputs@params_table["d","est"]
    estart = resDIVALIKE$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    
    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    # Add jump dispersal/founder-event speciation
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    resfn = paste(name,"_DIVALIKE+J_M0_unconstrained_v1.Rdata",sep='')
    runslow = TRUE
    if (runslow)
    {
      #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
      
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resDIVALIKEj = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKEj = res
    }
    
  } else if (model_name=='BAYAREALIKE'){
    #######################################################
    # Run BAYAREALIKE
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up BAYAREALIKE model
    # No subset sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    # No vicariance
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
    
    # No jump dispersal/founder-event speciation
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
    
    # Adjust linkage between parameters
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
    
    # Only sympatric/range-copying (y) events allowed, and with 
    # exact copying (both descendants always the same size as the ancestor)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
    
    # Check the inputs
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = paste(name,"_BAYAREALIKE_M0_unconstrained_v1.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resBAYAREALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resBAYAREALIKE = res
    }
    
  } else if (model_name=='BAYAREALIKE+J'){
  #######################################################
  # Run BAYAREALIKE
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = treefile
  BioGeoBEARS_run_object$geogfn = geographyfile
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
  BioGeoBEARS_run_object$num_cores_to_use = 1
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up BAYAREALIKE model
  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  # No jump dispersal/founder-event speciation
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
  
  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
  
  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  
  # Check the inputs
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  runslow = TRUE
  resfn = paste(name,"_BAYAREALIKE_M0_unconstrained_v1.Rdata",sep='')
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    resBAYAREALIKE = res
  } else {
    # Loads to "res"
    load(resfn)
    resBAYAREALIKE = res
  }
  
  #######################################################
  # Run BAYAREALIKE+J
  #######################################################
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = treefile
  BioGeoBEARS_run_object$geogfn = geographyfile
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
  #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
  #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
  #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
  # Also: search script on "include_null_range" for other places to change
  
  # Set up a time-stratified analysis:
  #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
  #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
  #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
  #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
  #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
  # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE
  BioGeoBEARS_run_object$num_cores_to_use = 1
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
  # The stratified tree is described in this table:
  #BioGeoBEARS_run_object$master_table
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up BAYAREALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resBAYAREALIKE$outputs@params_table["d","est"]
  estart = resBAYAREALIKE$outputs@params_table["e","est"]
  jstart = 0.0001
  
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
  
  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
  
  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  
  # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
  # machines. I can't replicate this on my Mac machines, but it is almost certainly
  # just some precision under-run issue, when optim/optimx tries some parameter value 
  # just below zero.  The "min" and "max" options on each parameter are supposed to
  # prevent this, but apparently optim/optimx sometimes go slightly beyond 
  # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
  # slightly for each parameter:
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  resfn = paste(name,"_BAYAREALIKE+J_M0_unconstrained_v1.Rdata",sep='')
  runslow = TRUE
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    
    save(res, file=resfn)
    
    resBAYAREALIKEj = res
  } else {
    # Loads to "res"
    load(resfn)
    resBAYAREALIKEj = res
  }
  
  }

  #return(res)
  #######################################################
  # Stochastic mapping on best model
  #######################################################
  clado_events_tables = NULL
  ana_events_tables = NULL
  lnum = 0
  
  #######################################################
  # Get the inputs for Biogeographical Stochastic Mapping
  # Note: this can be slow for large state spaces and trees, since 
  # the independent likelihoods for each branch are being pre-calculated
  # E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
  # for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
  # for storage of "BSM_inputs_file.Rdata".
  # Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
  # the same settings will be used for get_inputs_for_stochastic_mapping().
  #######################################################
  BSM_inputs_fn = "BSM_inputs_file.Rdata"
  runInputsSlow = TRUE
  if (runInputsSlow)
  {
    stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
    save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
  } else {
    # Loads to "stochastic_mapping_inputs_list"
    load(BSM_inputs_fn)
  } # END if (runInputsSlow)
  
  # Check inputs (doesn't work the same on unconstr)
  names(stochastic_mapping_inputs_list)
  stochastic_mapping_inputs_list$phy2
  stochastic_mapping_inputs_list$COO_weights_columnar
  stochastic_mapping_inputs_list$unconstr
  set.seed(seed=as.numeric(Sys.time()))
  
  runBSMslow = TRUE
  if (runBSMslow == TRUE)
  {
    # Saves to: RES_clado_events_tables.Rdata
    # Saves to: RES_ana_events_tables.Rdata
    BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=1000, nummaps_goal=50, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
    
    RES_clado_events_tables = BSM_output$RES_clado_events_tables
    RES_ana_events_tables = BSM_output$RES_ana_events_tables
  } else {
    # Load previously saved...
    
    # Loads to: RES_clado_events_tables
    load(file="RES_clado_events_tables.Rdata")
    # Loads to: RES_ana_events_tables
    load(file="RES_ana_events_tables.Rdata")
    BSM_output = NULL
    BSM_output$RES_clado_events_tables = RES_clado_events_tables
    BSM_output$RES_ana_events_tables = RES_ana_events_tables
  } # END if (runBSMslow == TRUE)
  
  # Extract BSM output
  clado_events_tables = BSM_output$RES_clado_events_tables
  ana_events_tables = BSM_output$RES_ana_events_tables
  head(clado_events_tables[[1]])
  head(ana_events_tables[[1]])
  length(clado_events_tables)
  length(ana_events_tables)
  
  #######################################################
  # Plot one stochastic map, manual method
  #######################################################
  # (we have to convert the stochastic maps into event
  #  maps for plotting)
  
  ######################
  # Get the color scheme
  ######################
  include_null_range = TRUE
  areanames = names(tipranges@df)
  areas = areanames
  max_range_size = 2
  
  # Note: If you did something to change the states_list from the default given the number of areas, you would
  # have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
  states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
  
  colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
  
  ############################################
  # Setup for painting a single stochastic map
  ############################################
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  stratified=FALSE
  clado_events_table = clado_events_tables[[1]]
  ana_events_table = ana_events_tables[[1]]
  
  # cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
  # colnums = match(cols_to_get, names(ana_events_table))
  # ana_events_table_cols_to_add = ana_events_table[,colnums]
  # anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
  # ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
  # rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
  # master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)
  
  ##############################################
  ### Open a PDF
  ##############################################
  ##pdffn = paste0(model_name, "_single_stochastic_map_n1.pdf")
  ##pdf(file=pdffn, width=6, height=6)
  ##
  ### Convert the BSM into a modified res object
  ##master_table_cladogenetic_events = clado_events_tables[[1]]
  ##resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)
  ##
  ##plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)
  ##
  ### Paint on the branch states
  ##paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=stratified)
  ##
  ##plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)
  ##
  ##############################################
  ### Close PDF
  ##############################################
  ##dev.off()
  ##cmdstr = paste("open ", pdffn, sep="")
  ##system(cmdstr)
  ##
  #########################################################
  ### Plot all 50 stochastic maps to PDF
  #########################################################
  ### Setup
  ##include_null_range = include_null_range
  ##areanames = areanames
  ##areas = areanames
  ##max_range_size = max_range_size
  ##states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
  ##colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
  ##scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  ##stratified = stratified
  ##
  ### Loop through the maps and plot to PDF
  ##pdffn = paste0(model_name, "_", length(clado_events_tables), "BSMs_v1.pdf")
  ##pdf(file=pdffn, width=6, height=6)
  ##
  ##nummaps_goal = 50
  ##for (i in 1:nummaps_goal)
  ##{
  ##  clado_events_table = clado_events_tables[[i]]
  ##  analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
  ##  plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range)
  ##} # END for (i in 1:nummaps_goal)
  ##
  ##dev.off()
  ##cmdstr = paste("open ", pdffn, sep="")
  ##system(cmdstr)
  ##
  #######################################################
  # Summarize stochastic map tables
  #######################################################
  length(clado_events_tables)
  length(ana_events_tables)
  
  head(clado_events_tables[[1]][,-20])
  tail(clado_events_tables[[1]][,-20])
  
  head(ana_events_tables[[1]])
  tail(ana_events_tables[[1]])
  
  areanames = names(tipranges@df)
  actual_names = areanames
  actual_names
  
  # Get the dmat and times (if any)
  dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
  dmat_times
  
  # Extract BSM output
  clado_events_tables = BSM_output$RES_clado_events_tables
  ana_events_tables = BSM_output$RES_ana_events_tables
  
  # Simulate the source areas
  BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
  clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
  ana_events_tables = BSMs_w_sourceAreas$ana_events_tables
  
  # Count all anagenetic and cladogenetic events
  counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)
  
  summary_counts_BSMs = counts_list$summary_counts_BSMs
  print(conditional_format_table(summary_counts_BSMs))
  
  # Histogram of event counts
  hist_event_counts(counts_list, pdffn=paste0(model_name, "_histograms_of_event_counts.pdf"))
  
  #######################################################
  # Check that ML ancestral state/range probabilities and
  # the mean of the BSMs approximately line up
  #######################################################
  #library(MultinomialCI)    # For 95% CIs on BSM counts
  #check_ML_vs_BSM(res, clado_events_tables, model_name, tr=NULL, plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE)
  return(counts_list)
  
}
run_BioGeoBEARS_selectedmodel_BSM_object<-function(treefile,geographyfile,path,name,model_name,nreplicates){
  #extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
  #list.files(extdata_dir)
  setwd(path)
  tr = read.tree(treefile)
  geogfn<-geographyfile
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  max_range_size = 2
  #####
  if (model_name=='DEC'){
    #######################################################
    #######################################################
    # DEC AND DEC+J ANALYSIS
    #######################################################
    #######################################################
    # NOTE: The BioGeoBEARS "DEC" model is identical with 
    # the Lagrange DEC model, and should return identical
    # ML estimates of parameters, and the same 
    # log-likelihoods, for the same datasets.
    #
    # Ancestral state probabilities at nodes will be slightly 
    # different, since BioGeoBEARS is reporting the 
    # ancestral state probabilities under the global ML
    # model, and Lagrange is reporting ancestral state
    # probabilities after re-optimizing the likelihood
    # after fixing the state at each node. These will 
    # be similar, but not identical. See Matzke (2014),
    # Systematic Biology, for discussion.
    #
    # Also see Matzke (2014) for presentation of the 
    # DEC+J model.
    #######################################################
    #######################################################
    
    #######################################################
    #######################################################
    
    #######################################################
    # Run DEC
    #######################################################
    
    # Intitialize a default model (DEC model)
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    
    # Give BioGeoBEARS the location of the phylogeny Newick file
    BioGeoBEARS_run_object$trfn = treefile
    
    # Give BioGeoBEARS the location of the geography text file
    BioGeoBEARS_run_object$geogfn = geographyfile
    
    # Input the maximum range size
    BioGeoBEARS_run_object$max_range_size = max_range_size
    
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    # 1. Here, un-comment ONLY the files you want to use.
    # 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
    # 3. For example files see (a) extdata_dir, 
    #  or (b) http://phylo.wikidot.com/biogeobears#files
    #  and BioGeoBEARS Google Group posts for further hints)
    #
    # Uncomment files you wish to use in time-stratified analyses:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    # (use more cores to speed it up; this requires
    # library(parallel) and/or library(snow). The package "parallel" 
    # is now default on Macs in R 3.0+, but apparently still 
    # has to be typed on some Windows machines. Note: apparently 
    # parallel works on Mac command-line R, but not R.app.
    # BioGeoBEARS checks for this and resets to 1
    # core with R.app)
    
    # Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
    # I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
    # but the results are imprecise and so I haven't explored it further.
    # In a Bayesian analysis, it might work OK, but the ML point estimates are
    # not identical.
    # Also, I have not implemented all functions to work with force_sparse=TRUE.
    # Volunteers are welcome to work on it!!
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DEC model
    # (nothing to do; defaults)
    
    # Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
    BioGeoBEARS_run_object
    
    # This contains the model object
    BioGeoBEARS_run_object$BioGeoBEARS_model_object
    
    # This table contains the parameters of the model 
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
    
    # Run this to check inputs. Read the error messages if you get them!
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # For a slow analysis, run once, then set runslow=FALSE to just 
    # load the saved result.
    runslow = TRUE
    resfn = paste(name,"_DEC.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resDEC = res
    } else {
      # Loads to "res"
      load(resfn)
      resDEC = res
    }
    
  } else if (model_name=='DEC+J'){
    #######################################################
    # Run DEC+J
    #######################################################
    # Intitialize a default model (DEC model)
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    
    # Give BioGeoBEARS the location of the phylogeny Newick file
    BioGeoBEARS_run_object$trfn = treefile
    
    # Give BioGeoBEARS the location of the geography text file
    BioGeoBEARS_run_object$geogfn = geographyfile
    
    # Input the maximum range size
    BioGeoBEARS_run_object$max_range_size = max_range_size
    
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    # 1. Here, un-comment ONLY the files you want to use.
    # 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
    # 3. For example files see (a) extdata_dir, 
    #  or (b) http://phylo.wikidot.com/biogeobears#files
    #  and BioGeoBEARS Google Group posts for further hints)
    #
    # Uncomment files you wish to use in time-stratified analyses:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    # (use more cores to speed it up; this requires
    # library(parallel) and/or library(snow). The package "parallel" 
    # is now default on Macs in R 3.0+, but apparently still 
    # has to be typed on some Windows machines. Note: apparently 
    # parallel works on Mac command-line R, but not R.app.
    # BioGeoBEARS checks for this and resets to 1
    # core with R.app)
    
    # Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
    # I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
    # but the results are imprecise and so I haven't explored it further.
    # In a Bayesian analysis, it might work OK, but the ML point estimates are
    # not identical.
    # Also, I have not implemented all functions to work with force_sparse=TRUE.
    # Volunteers are welcome to work on it!!
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DEC model
    # (nothing to do; defaults)
    
    # Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
    BioGeoBEARS_run_object
    
    # This contains the model object
    BioGeoBEARS_run_object$BioGeoBEARS_model_object
    
    # This table contains the parameters of the model 
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
    
    # Run this to check inputs. Read the error messages if you get them!
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # For a slow analysis, run once, then set runslow=FALSE to just 
    # load the saved result.
    runslow = TRUE
    resfn = paste(name,"_DEC.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resDEC = res
    } else {
      # Loads to "res"
      load(resfn)
      resDEC = res
    }
    
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DEC+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resDEC$outputs@params_table["d","est"]
    estart = resDEC$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    
    # Add j as a free parameter
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    resfn = paste(name,"_DEC+J.Rdata",sep='')
    runslow = TRUE
    if (runslow)
    {
      #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
      
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resDECj = res
    } else {
      # Loads to "res"
      load(resfn)
      resDECj = res
    }
    
  } else if (model_name=='DIVALIKE'){
    
    #######################################################
    #######################################################
    # DIVALIKE AND DIVALIKE+J ANALYSIS
    #######################################################
    #######################################################
    # NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with 
    # Ronquist (1997)'s parsimony DIVA. It is a likelihood
    # interpretation of DIVA, constructed by modelling DIVA's
    # processes the way DEC does, but only allowing the 
    # processes DIVA allows (widespread vicariance: yes; subset
    # sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
    #
    # DIVALIKE is a likelihood interpretation of parsimony
    # DIVA, and it is "like DIVA" -- similar to, but not
    # identical to, parsimony DIVA.
    #
    # I thus now call the model "DIVALIKE", and you should also. ;-)
    #######################################################
    #######################################################
    
    #######################################################
    # Run DIVALIKE
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DIVALIKE model
    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    # No jump dispersal/founder-event speciation
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = paste(name,"_DIVALIKE_M0_unconstrained_v1.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resDIVALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKE = res
    }
    
  } else if(model_name=='DIVALIKE+J'){
    #######################################################
    # Run DIVALIKE
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DIVALIKE model
    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    # No jump dispersal/founder-event speciation
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = paste(name,"_DIVALIKE_M0_unconstrained_v1.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resDIVALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKE = res
    }
    
    
    #######################################################
    # Run DIVALIKE+J
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DIVALIKE+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resDIVALIKE$outputs@params_table["d","est"]
    estart = resDIVALIKE$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    
    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    # Add jump dispersal/founder-event speciation
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    resfn = paste(name,"_DIVALIKE+J_M0_unconstrained_v1.Rdata",sep='')
    runslow = TRUE
    if (runslow)
    {
      #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
      
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resDIVALIKEj = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKEj = res
    }
    
  } else if (model_name=='BAYAREALIKE'){
    #######################################################
    # Run BAYAREALIKE
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up BAYAREALIKE model
    # No subset sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    # No vicariance
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
    
    # No jump dispersal/founder-event speciation
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
    
    # Adjust linkage between parameters
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
    
    # Only sympatric/range-copying (y) events allowed, and with 
    # exact copying (both descendants always the same size as the ancestor)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
    
    # Check the inputs
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = paste(name,"_BAYAREALIKE_M0_unconstrained_v1.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resBAYAREALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resBAYAREALIKE = res
    }
    
  } else if (model_name=='BAYAREALIKE+J'){
    #######################################################
    # Run BAYAREALIKE
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up BAYAREALIKE model
    # No subset sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    # No vicariance
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
    
    # No jump dispersal/founder-event speciation
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
    
    # Adjust linkage between parameters
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
    
    # Only sympatric/range-copying (y) events allowed, and with 
    # exact copying (both descendants always the same size as the ancestor)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
    
    # Check the inputs
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = paste(name,"_BAYAREALIKE_M0_unconstrained_v1.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resBAYAREALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resBAYAREALIKE = res
    }
    
    #######################################################
    # Run BAYAREALIKE+J
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up BAYAREALIKE+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resBAYAREALIKE$outputs@params_table["d","est"]
    estart = resBAYAREALIKE$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    
    # No subset sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    # No vicariance
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
    
    # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
    
    # Adjust linkage between parameters
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
    
    # Only sympatric/range-copying (y) events allowed, and with 
    # exact copying (both descendants always the same size as the ancestor)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
    
    # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
    # machines. I can't replicate this on my Mac machines, but it is almost certainly
    # just some precision under-run issue, when optim/optimx tries some parameter value 
    # just below zero.  The "min" and "max" options on each parameter are supposed to
    # prevent this, but apparently optim/optimx sometimes go slightly beyond 
    # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
    # slightly for each parameter:
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    resfn = paste(name,"_BAYAREALIKE+J_M0_unconstrained_v1.Rdata",sep='')
    runslow = TRUE
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resBAYAREALIKEj = res
    } else {
      # Loads to "res"
      load(resfn)
      resBAYAREALIKEj = res
    }
    
  }
  
  #return(res)
  #######################################################
  # Stochastic mapping on best model
  #######################################################
  clado_events_tables = NULL
  ana_events_tables = NULL
  lnum = 0
  
  #######################################################
  # Get the inputs for Biogeographical Stochastic Mapping
  # Note: this can be slow for large state spaces and trees, since 
  # the independent likelihoods for each branch are being pre-calculated
  # E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
  # for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
  # for storage of "BSM_inputs_file.Rdata".
  # Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
  # the same settings will be used for get_inputs_for_stochastic_mapping().
  #######################################################
  BSM_inputs_fn = "BSM_inputs_file.Rdata"
  runInputsSlow = TRUE
  if (runInputsSlow)
  {
    stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
    save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
  } else {
    # Loads to "stochastic_mapping_inputs_list"
    load(BSM_inputs_fn)
  } # END if (runInputsSlow)
  
  # Check inputs (doesn't work the same on unconstr)
  names(stochastic_mapping_inputs_list)
  stochastic_mapping_inputs_list$phy2
  stochastic_mapping_inputs_list$COO_weights_columnar
  stochastic_mapping_inputs_list$unconstr
  set.seed(seed=as.numeric(Sys.time()))
  
  runBSMslow = TRUE
  if (runBSMslow == TRUE)
  {
    # Saves to: RES_clado_events_tables.Rdata
    # Saves to: RES_ana_events_tables.Rdata
    BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=1000, nummaps_goal=nreplicates, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
    
    RES_clado_events_tables = BSM_output$RES_clado_events_tables
    RES_ana_events_tables = BSM_output$RES_ana_events_tables
  } else {
    # Load previously saved...
    
    # Loads to: RES_clado_events_tables
    load(file="RES_clado_events_tables.Rdata")
    # Loads to: RES_ana_events_tables
    load(file="RES_ana_events_tables.Rdata")
    BSM_output = NULL
    BSM_output$RES_clado_events_tables = RES_clado_events_tables
    BSM_output$RES_ana_events_tables = RES_ana_events_tables
  } # END if (runBSMslow == TRUE)
  
  
  return(list(res,BSM_output))
  
}

run_BioGeoBEARS_selectedmodel_BSM_object_old<-function(treefile,geographyfile,path,name,model_name){
  #extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
  #list.files(extdata_dir)
  setwd(path)
  tr = read.tree(treefile)
  geogfn<-geographyfile
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  max_range_size = 2
  #####
  if (model_name=='DEC'){
    #######################################################
    #######################################################
    # DEC AND DEC+J ANALYSIS
    #######################################################
    #######################################################
    # NOTE: The BioGeoBEARS "DEC" model is identical with 
    # the Lagrange DEC model, and should return identical
    # ML estimates of parameters, and the same 
    # log-likelihoods, for the same datasets.
    #
    # Ancestral state probabilities at nodes will be slightly 
    # different, since BioGeoBEARS is reporting the 
    # ancestral state probabilities under the global ML
    # model, and Lagrange is reporting ancestral state
    # probabilities after re-optimizing the likelihood
    # after fixing the state at each node. These will 
    # be similar, but not identical. See Matzke (2014),
    # Systematic Biology, for discussion.
    #
    # Also see Matzke (2014) for presentation of the 
    # DEC+J model.
    #######################################################
    #######################################################
    
    #######################################################
    #######################################################
    
    #######################################################
    # Run DEC
    #######################################################
    
    # Intitialize a default model (DEC model)
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    
    # Give BioGeoBEARS the location of the phylogeny Newick file
    BioGeoBEARS_run_object$trfn = treefile
    
    # Give BioGeoBEARS the location of the geography text file
    BioGeoBEARS_run_object$geogfn = geographyfile
    
    # Input the maximum range size
    BioGeoBEARS_run_object$max_range_size = max_range_size
    
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    # 1. Here, un-comment ONLY the files you want to use.
    # 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
    # 3. For example files see (a) extdata_dir, 
    #  or (b) http://phylo.wikidot.com/biogeobears#files
    #  and BioGeoBEARS Google Group posts for further hints)
    #
    # Uncomment files you wish to use in time-stratified analyses:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    # (use more cores to speed it up; this requires
    # library(parallel) and/or library(snow). The package "parallel" 
    # is now default on Macs in R 3.0+, but apparently still 
    # has to be typed on some Windows machines. Note: apparently 
    # parallel works on Mac command-line R, but not R.app.
    # BioGeoBEARS checks for this and resets to 1
    # core with R.app)
    
    # Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
    # I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
    # but the results are imprecise and so I haven't explored it further.
    # In a Bayesian analysis, it might work OK, but the ML point estimates are
    # not identical.
    # Also, I have not implemented all functions to work with force_sparse=TRUE.
    # Volunteers are welcome to work on it!!
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DEC model
    # (nothing to do; defaults)
    
    # Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
    BioGeoBEARS_run_object
    
    # This contains the model object
    BioGeoBEARS_run_object$BioGeoBEARS_model_object
    
    # This table contains the parameters of the model 
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
    
    # Run this to check inputs. Read the error messages if you get them!
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # For a slow analysis, run once, then set runslow=FALSE to just 
    # load the saved result.
    runslow = TRUE
    resfn = paste(name,"_DEC.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resDEC = res
    } else {
      # Loads to "res"
      load(resfn)
      resDEC = res
    }
    
  } else if (model_name=='DEC+J'){
    #######################################################
    # Run DEC+J
    #######################################################
    # Intitialize a default model (DEC model)
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    
    # Give BioGeoBEARS the location of the phylogeny Newick file
    BioGeoBEARS_run_object$trfn = treefile
    
    # Give BioGeoBEARS the location of the geography text file
    BioGeoBEARS_run_object$geogfn = geographyfile
    
    # Input the maximum range size
    BioGeoBEARS_run_object$max_range_size = max_range_size
    
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    # 1. Here, un-comment ONLY the files you want to use.
    # 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
    # 3. For example files see (a) extdata_dir, 
    #  or (b) http://phylo.wikidot.com/biogeobears#files
    #  and BioGeoBEARS Google Group posts for further hints)
    #
    # Uncomment files you wish to use in time-stratified analyses:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    # (use more cores to speed it up; this requires
    # library(parallel) and/or library(snow). The package "parallel" 
    # is now default on Macs in R 3.0+, but apparently still 
    # has to be typed on some Windows machines. Note: apparently 
    # parallel works on Mac command-line R, but not R.app.
    # BioGeoBEARS checks for this and resets to 1
    # core with R.app)
    
    # Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
    # I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
    # but the results are imprecise and so I haven't explored it further.
    # In a Bayesian analysis, it might work OK, but the ML point estimates are
    # not identical.
    # Also, I have not implemented all functions to work with force_sparse=TRUE.
    # Volunteers are welcome to work on it!!
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DEC model
    # (nothing to do; defaults)
    
    # Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
    BioGeoBEARS_run_object
    
    # This contains the model object
    BioGeoBEARS_run_object$BioGeoBEARS_model_object
    
    # This table contains the parameters of the model 
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
    
    # Run this to check inputs. Read the error messages if you get them!
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # For a slow analysis, run once, then set runslow=FALSE to just 
    # load the saved result.
    runslow = TRUE
    resfn = paste(name,"_DEC.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resDEC = res
    } else {
      # Loads to "res"
      load(resfn)
      resDEC = res
    }
    
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DEC+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resDEC$outputs@params_table["d","est"]
    estart = resDEC$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    
    # Add j as a free parameter
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    resfn = paste(name,"_DEC+J.Rdata",sep='')
    runslow = TRUE
    if (runslow)
    {
      #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
      
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resDECj = res
    } else {
      # Loads to "res"
      load(resfn)
      resDECj = res
    }
    
  } else if (model_name=='DIVALIKE'){
    
    #######################################################
    #######################################################
    # DIVALIKE AND DIVALIKE+J ANALYSIS
    #######################################################
    #######################################################
    # NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with 
    # Ronquist (1997)'s parsimony DIVA. It is a likelihood
    # interpretation of DIVA, constructed by modelling DIVA's
    # processes the way DEC does, but only allowing the 
    # processes DIVA allows (widespread vicariance: yes; subset
    # sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
    #
    # DIVALIKE is a likelihood interpretation of parsimony
    # DIVA, and it is "like DIVA" -- similar to, but not
    # identical to, parsimony DIVA.
    #
    # I thus now call the model "DIVALIKE", and you should also. ;-)
    #######################################################
    #######################################################
    
    #######################################################
    # Run DIVALIKE
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DIVALIKE model
    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    # No jump dispersal/founder-event speciation
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = paste(name,"_DIVALIKE_M0_unconstrained_v1.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resDIVALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKE = res
    }
    
  } else if(model_name=='DIVALIKE+J'){
    #######################################################
    # Run DIVALIKE
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DIVALIKE model
    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    # No jump dispersal/founder-event speciation
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = paste(name,"_DIVALIKE_M0_unconstrained_v1.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resDIVALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKE = res
    }
    
    
    #######################################################
    # Run DIVALIKE+J
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up DIVALIKE+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resDIVALIKE$outputs@params_table["d","est"]
    estart = resDIVALIKE$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    
    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
    
    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
    
    # Add jump dispersal/founder-event speciation
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    resfn = paste(name,"_DIVALIKE+J_M0_unconstrained_v1.Rdata",sep='')
    runslow = TRUE
    if (runslow)
    {
      #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
      
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resDIVALIKEj = res
    } else {
      # Loads to "res"
      load(resfn)
      resDIVALIKEj = res
    }
    
  } else if (model_name=='BAYAREALIKE'){
    #######################################################
    # Run BAYAREALIKE
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up BAYAREALIKE model
    # No subset sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    # No vicariance
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
    
    # No jump dispersal/founder-event speciation
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
    
    # Adjust linkage between parameters
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
    
    # Only sympatric/range-copying (y) events allowed, and with 
    # exact copying (both descendants always the same size as the ancestor)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
    
    # Check the inputs
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = paste(name,"_BAYAREALIKE_M0_unconstrained_v1.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resBAYAREALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resBAYAREALIKE = res
    }
    
  } else if (model_name=='BAYAREALIKE+J'){
    #######################################################
    # Run BAYAREALIKE
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE     # if FALSE, use optim() instead of optimx()
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up BAYAREALIKE model
    # No subset sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    # No vicariance
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
    
    # No jump dispersal/founder-event speciation
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
    # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
    
    # Adjust linkage between parameters
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
    
    # Only sympatric/range-copying (y) events allowed, and with 
    # exact copying (both descendants always the same size as the ancestor)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
    
    # Check the inputs
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    runslow = TRUE
    resfn = paste(name,"_BAYAREALIKE_M0_unconstrained_v1.Rdata",sep='')
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      resBAYAREALIKE = res
    } else {
      # Loads to "res"
      load(resfn)
      resBAYAREALIKE = res
    }
    
    #######################################################
    # Run BAYAREALIKE+J
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = treefile
    BioGeoBEARS_run_object$geogfn = geographyfile
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change
    
    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
    
    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
    
    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
    
    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
    
    # Set up BAYAREALIKE+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resBAYAREALIKE$outputs@params_table["d","est"]
    estart = resBAYAREALIKE$outputs@params_table["e","est"]
    jstart = 0.0001
    
    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
    
    # No subset sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
    
    # No vicariance
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
    
    # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
    
    # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
    
    # Adjust linkage between parameters
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
    
    # Only sympatric/range-copying (y) events allowed, and with 
    # exact copying (both descendants always the same size as the ancestor)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
    
    # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
    # machines. I can't replicate this on my Mac machines, but it is almost certainly
    # just some precision under-run issue, when optim/optimx tries some parameter value 
    # just below zero.  The "min" and "max" options on each parameter are supposed to
    # prevent this, but apparently optim/optimx sometimes go slightly beyond 
    # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
    # slightly for each parameter:
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999
    
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
    
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)
    
    resfn = paste(name,"_BAYAREALIKE+J_M0_unconstrained_v1.Rdata",sep='')
    runslow = TRUE
    if (runslow)
    {
      res = bears_optim_run(BioGeoBEARS_run_object)
      res    
      
      save(res, file=resfn)
      
      resBAYAREALIKEj = res
    } else {
      # Loads to "res"
      load(resfn)
      resBAYAREALIKEj = res
    }
    
  }
  
  #return(res)
  #######################################################
  # Stochastic mapping on best model
  #######################################################
  clado_events_tables = NULL
  ana_events_tables = NULL
  lnum = 0
  
  #######################################################
  # Get the inputs for Biogeographical Stochastic Mapping
  # Note: this can be slow for large state spaces and trees, since 
  # the independent likelihoods for each branch are being pre-calculated
  # E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
  # for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
  # for storage of "BSM_inputs_file.Rdata".
  # Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
  # the same settings will be used for get_inputs_for_stochastic_mapping().
  #######################################################
  BSM_inputs_fn = "BSM_inputs_file.Rdata"
  runInputsSlow = TRUE
  if (runInputsSlow)
  {
    stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
    save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
  } else {
    # Loads to "stochastic_mapping_inputs_list"
    load(BSM_inputs_fn)
  } # END if (runInputsSlow)
  
  # Check inputs (doesn't work the same on unconstr)
  names(stochastic_mapping_inputs_list)
  stochastic_mapping_inputs_list$phy2
  stochastic_mapping_inputs_list$COO_weights_columnar
  stochastic_mapping_inputs_list$unconstr
  set.seed(seed=as.numeric(Sys.time()))
  
  runBSMslow = TRUE
  if (runBSMslow == TRUE)
  {
    # Saves to: RES_clado_events_tables.Rdata
    # Saves to: RES_ana_events_tables.Rdata
    BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=1000, nummaps_goal=50, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
    
    RES_clado_events_tables = BSM_output$RES_clado_events_tables
    RES_ana_events_tables = BSM_output$RES_ana_events_tables
  } else {
    # Load previously saved...
    
    # Loads to: RES_clado_events_tables
    load(file="RES_clado_events_tables.Rdata")
    # Loads to: RES_ana_events_tables
    load(file="RES_ana_events_tables.Rdata")
    BSM_output = NULL
    BSM_output$RES_clado_events_tables = RES_clado_events_tables
    BSM_output$RES_ana_events_tables = RES_ana_events_tables
  } # END if (runBSMslow == TRUE)
  
  
  return(BSM_output)
  
}

summarise_BSMobject<-function(geographyfile,path,name,BSM.objectfile){
  BSM.object<-readRDS(BSM.objectfile)
  setwd(path)
  geogfn<-geographyfile
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  res<-BSM.object[[1]]
  BSM_output<-BSM.object[[2]]
  # Extract BSM output
  clado_events_tables = BSM_output$RES_clado_events_tables
  ana_events_tables = BSM_output$RES_ana_events_tables
  head(clado_events_tables[[1]])
  head(ana_events_tables[[1]])
  length(clado_events_tables)
  length(ana_events_tables)
  
  clado_events_table = clado_events_tables[[1]]
  ana_events_table = ana_events_tables[[1]]
  
  # cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
  # colnums = match(cols_to_get, names(ana_events_table))
  # ana_events_table_cols_to_add = ana_events_table[,colnums]
  # anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
  # ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
  # rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
  # master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)
  
  length(clado_events_tables)
  length(ana_events_tables)
  
  head(clado_events_tables[[1]][,-20])
  tail(clado_events_tables[[1]][,-20])
  
  head(ana_events_tables[[1]])
  tail(ana_events_tables[[1]])
  
  areanames = names(tipranges@df)
  actual_names = areanames
  actual_names
  
  # Get the dmat and times (if any)
  dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
  dmat_times
  
  # Extract BSM output
  clado_events_tables = BSM_output$RES_clado_events_tables
  ana_events_tables = BSM_output$RES_ana_events_tables
  
  # Simulate the source areas
  BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
  clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
  ana_events_tables = BSMs_w_sourceAreas$ana_events_tables
  
  # Count all anagenetic and cladogenetic events
  counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)
  
  summary_counts_BSMs = counts_list$summary_counts_BSMs
  return(counts_list)
}
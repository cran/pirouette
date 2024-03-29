#' This function does nothing. It is intended to inherit is parameters'
#' documentation.
#' @param alignment a DNA alignment, of class \link[ape]{DNAbin}
#' @param alignment_params parameters to simulate an alignment,
#'   as can be created by \link{create_alignment_params}
#' @param alignment_rng_seed The random number generator seed used
#'   to generate an alignment
#' @param base_frequencies the four base frequencies (a, c, g, t) to be
#'   specified to create the rate matrix (i.e. Q matrix)
#'   used to simulate alignments
#' @param bd_mutation_rate the mutation rate when creating an alignment
#'   from a BD tree
#' @param bd_tree a phylogent of class \link[ape]{phylo},
#'   created by a Birth Death process
#' @param bd_tree_filename name of the file that stores a BD twin tree
#' @param beast2_bin_path path to BEAST2 binary file. The use of the
#'   binary BEAST2 file is required for estimation of the evidence (aka
#'   marginal likelihood). The default BEAST2 binary path can be
#'   obtained using \link[beastier]{get_default_beast2_bin_path}
#' @param beast2_input_filename path of the BEAST2 configuration file.
#'   By default, this file is put in a temporary folder with a random filename,
#'   as the user needs not read it: it is used as input of BEAST2.
#'   Specifying a \code{beast2_input_filename} allows
#'   to store that file in a more permanently stored location.
#' @param beast2_options BEAST2 options,
#'   as can be created by \link[beastier]{create_beast2_options}
#' @param beast2_optionses list of one or more BEAST2 options,
#'   as can be created by \link[beastier]{create_beast2_options}
#' @param beast2_options_inference BEAST2 options,
#'   as can be created by \link[beastier]{create_beast2_options}.
#'   The MCMC must be a normal MCMC,
#'   as can be created by \link[beautier]{create_mcmc}.
#' @param beast2_options_est_evidence BEAST2 options to estimate
#'   the evidence (aka marginal likelihood),
#'   as can be created by \link[beastier]{create_beast2_options}.
#'   The MCMC must be a Nested Sampling MCMC,
#'   as can be created by \link[beautier]{create_ns_mcmc}.
#' @param beast2_output_log_filename name of the log file created by BEAST2,
#'   containing the parameter estimates in time.
#'   By default, this file is put a temporary folder with a random filename,
#'   as the user needs not read it.
#'   Specifying a beast2_output_log_filename allows to store that file
#'   in a more permanently stored location.
#' @param beast2_output_state_filename name of the final state file
#'   created by BEAST2, containing the operator acceptances.
#'   By default, this file is put a temporary folder with a random filename,
#'   as the user needs not read it.
#'   Specifying a beast2_output_state_filename allows to store
#'   that file in a more permanently stored location.
#' @param beast2_output_trees_filename name of a trees files
#'   created by BEAST2.
#'   By default, this file is put a temporary folder with a random filename,
#'   as the user needs not read it: its content is parsed and
#'   compared to a true phylogeny to obtain the inference errors.
#'   Specifying \code{beast2_output_trees_filename} allows to store
#'   this file in a more permanently stored location.
#' @param beast2_output_trees_filenames	name of the one or more trees files
#'   created by BEAST2, one per alignment.
#'   By default, these files are put a temporary folder with a random filename,
#'   as the user needs not read it: its content is parsed and
#'   compared to a true phylogeny to obtain the inference errors.
#'   Specifying \code{beast2_output_trees_filenames} allows to store
#'   these one or more files in a more permanently stored location.
#' @param beast2_path Path to the
#'   BEAST2 jar file (\code{beast.jar})
#'   or BEAST2 binary file '(\code{beast})'.
#'   Use \link[beastier]{get_default_beast2_jar_path} for the default
#'   BEAST2 jar file path.
#'   Use \link[beastier]{get_default_beast2_bin_path} for the default
#'   BEAST2 binary file path.
#' @param beast2_rng_seed The random number generator seed used by BEAST2
#' @param branch_mutation_rate mutation rate along the branch.
#' See, among others, \link[nodeSub]{sim_unlinked} for more details
#' @param branch_subst_matrix substitution matrix along the branches.
#' See, among others, \link[nodeSub]{sim_unlinked} for more details
#' @param brts numeric vector of (all postive) branching times,
#'   in time units before the present. Assuming no stem, the heighest
#'   value equals the crown age.
#' @param burn_in_fraction the fraction of the posterior trees (starting
#'   from the ones generated first)
#'   that will be discarded,
#'   must be a value from 0.0 (keep all), to 1.0 (discard all).
#' @param chain_length something
#' @param check_input boolean to indicate if the input is checked.
#'   If set to \link{TRUE}, input is checked, resulting in a proper
#'   error message. Else, input is left unchecked, possibly resulting
#'   in unhelpful error messages.
#' @param clock_model a clock model,
#'   as created by \link[beautier]{create_clock_model}
#' @param clock_models a list of one or more clock models,
#'   as created by \link[beautier]{create_clock_model}
#' @param clock_model_name name of a clock model
#' @param consensus the order of which the taxon labels are plotted
#' @param crown_age the fixed crown age of the posterior. Set to NA
#'   to let it be estimated
#' @param df_long the output created by \code{\link{pir_run}} in the long form
#' @param do_measure_evidence boolean to indicate if the
#'   evidence (aka marginal likelihood) of an experiment must be
#'   measured
#' @param epsilon	measure of relative accuracy when estimating a model's
#'   evidence (also known as marginal likelihood).
#'   Smaller values result in more precise estimations, that take
#'   longer to compute
#' @param error_fun function that determines the error between
#'   a given phylogeny and a the trees in a Bayesian posterior.
#'   The function must have two arguments:
#'   \itemize{
#'     \item the one given phylogeny, of class \link[ape]{phylo}
#'     \item one or more posterior trees, of class \link[ape]{multiphylo}
#'   }
#'   The function must return as many errors as there are posterior
#'   trees given. The error must be lowest between identical trees.
#'   Example functions are:
#'   \itemize{
#'     \item \link{get_gamma_error_fun}: use the absolute difference
#'       in gamma statistic
#'     \item \link{get_nltt_error_fun}: use the nLTT statistic
#'   }
#' @param error_measure_params parameter set to specify how the
#'   error between the given phylogeny and the Bayesian
#'   posterior is determined.
#'   Use \link{create_error_measure_params} to create such
#'   a parameter set
#' @param errors a numeric vector of (positive) Bayesian inference errors.
#'   Use \link{NA} if these are not measured (yet)
#' @param errors_filename baseline name for errors filenames,
#' as created by \link{get_temp_errors_filename}
#' @param est_evidence_mcmc MCMC used in the estimation of
#'   the evidence (aka marginal likelihood).
#'   The MCMC must be a Nested Sampling MCMC,
#'   as can be created by \link[beautier]{create_ns_mcmc}.
#' @param evidence_epsilon relative error in estimating the
#'   evidence (aka marginal likelihood).
#' @param evidence_filename filename to store the estimated
#' evidences (aka marginal likelihoods),
#' as can be created by \link{get_temp_evidence_filename}.
#' Must be \link{NA} if there is evidence
#' estimation (as determined by \link{will_measure_evidence}).
#' @param exclude_model an inference model that has to be excluded, as can be
#'   created by \link[beautier]{create_inference_model}
#' @param experiment a \link{pirouette} experiment,
#'   as can be created by \link{create_experiment}
#' @param experiments a list of one or more \link{pirouette} experiments,
#'   as can be created by \link{create_experiment}. If more than one experiment
#'   is provided and a "generative" experiment is part of them, the "generative"
#'   one has to be the first in the list. See also:
#'   \itemize{
#'     \item Use \link{check_experiments} to check the list of
#'       experiments for validity
#'     \item Use \link{create_all_experiments} to create experiments with
#'       all combinations of tree model, clock model and tree priors
#'     \item Use \link{create_all_bd_experiments} to create experiments
#'       with all combinations of tree model, clock model and tree priors,
#'       except for only using birth-death tree priors
#'     \item Use \link{create_all_coal_experiments} to create all experiments
#'       with all combinations of tree model, clock model and tree priors,
#'       except for only coalescent tree priors
#'     \item Use \link{shorten_experiments} to shorten the run time
#'       of the list of experiments
#'   }
#' @param extinction_rate per-species extinction rate
#' @param fasta_filename name of a FASTA file.
#'   Use \link{get_alignment_id} to get the ID of the alignment
#' @param filename the file's name, without the path
#' @param folder_name name of the main folder
#' @param folder_names one or more folder names
#' @param inference_conditions conditions under which the inference model
#'   is used in the inference
#' @param inference_model an inference model, which is a combination
#'   of site model, clock model, tree prior and BEAST2 input and
#'   input filenames.
#' @param init_speciation_rate a speciation rate
#' @param init_extinction_rate an extinction rate
#' @param lambda per-lineage speciation rate
#' @param log_evidence the natural logarithm of the evidence (aka marginal
#'   likelihood). Can be NA if this is not measured
#' @param marg_lik_filename name of the file the marginal
#'   likelihoods (also known as 'evidences') are saved to
#' @param marg_liks a data frame with marginal likelihoods/evidences.
#'   A test data frame can be created by \link{create_test_marg_liks}
#' @param max_evidence_epsilon set the maximum acceptable threshold for the
#'   parameter \code{evidence_epsilon}
#' @param max_n_tries maximum number of tries before giving up
#' @param mcmc MCMC options, as created by \link[beautier]{create_mcmc}
#' @param mbd_l_matrix the L matrix of an MBD tree
#' @param mbd_mutation_rate the mutation rate when creating an alignment
#'   from a MBD tree
#' @param mbd_tree an MBD tree
#' @param method determines how to create the twin tree
#' \itemize{
#'     \item 'random_tree' just produces a random tree;
#'     \item 'max_clade_cred' simulates \code{n_replicates} trees and
#'       uses \link[phangorn]{maxCladeCred} to create a consensus tree;
#'     \item 'max_likelihood' simulates \code{n_replicates} trees
#'      and selects the most likely;
#'   }
#' @param model_selection one ways to select the models used in
#'   inference, for example, \code{generative} picks the generative
#'   model, where \code{most_evidence} picks the model with most
#'   evidence. See \link{get_model_selections} for a list of
#' @param model_type type of inference model supplied for an experiment.
#'   Possible values:
#'   \itemize{
#'     \item \code{generative}: the inference model is (or is assumed to be)
#'       the inference model underlying the phylogeny
#'     \item \code{candidate}: the inference model is a candidate model,
#'       that competes with other models for having the most
#'       evidence (aka highest marginal likelihood)
#'   }
#' @param mrca_prior an MRCA prior,
#'   as created by \link[beautier]{create_mrca_prior}
#' @param mu per-species extinction rate
#' @param mutation_rate the mutation rate per base pair per time unit.
#'   Use \link{check_mutation_rate} to check if a mutation rate is valid.
#' @param n_0 number of starting species
#' @param n_mutations costrained number of mutations
#' @param n_taxa number of tree tips
#' @param n_replicates number of replicas to evaluate in order to create the
#'   twin tree
#' @param node_mutation_rate mutation rate on the node.
#' See, among others, \link[nodeSub]{sim_unlinked} for more details
#' @param node_subst_matrix substitution matrix on the nodes.
#' See, among others, \link[nodeSub]{sim_unlinked} for more details
#' @param node_time amount of time spent at the nodes.
#' See, among others, \link[nodeSub]{sim_unlinked} for more details
#' @param nu the rate at which a multiple-birth specation is triggered
#' @param nu_events the number of nu-triggered events that have to be
#'  present in the simulated tree
#' @param os name of the operating system, can be \code{mac}, \code{unix}
#'   or \code{win}. Use \link[beastier]{check_os} if the operating system
#'   is valid.
#' @param parameter_filename full path to a 'parameters.csv' file
#' @param parameters_filename full path to a 'parameters.csv' file
#' @param phylo a phylogeny of class \link[ape]{phylo}
#' @param phylogenies a list of phylogenies,
#'   each phylogeny being of class \link[ape]{phylo}
#' @param phylogeny a phylogeny of class \link[ape]{phylo}
#' @param pir_params the parameters of \link[pirouette]{pirouette}.
#'   They are created by \link{create_pir_params}.
#' @param pir_paramses a list of \link[pirouette]{pirouette} parameters,
#'   each element created by \link{create_pir_params}.
#' @param pir_out the output of \link{pir_run}
#' @param pir_outs the output of \link{pir_runs}
#' @param posterior_trees phylogenetic trees in a BEAST2 posterior,
#'   of class \code{multiphylo}
#' @param precision define the precision of the approximation.
#' @param project_folder_name project folder name
#' @param rename_fun a function to rename a filename,
#' as can be checked by \link{check_rename_fun}. This function should
#' have one argument, which will be a filename or \link{NA}. The
#' function should \link{return} one filename (when passed one filename) or
#' one \link{NA} (when passed one \link{NA}).
#' Example rename functions are:
#' \itemize{
#'   \item \link{get_remove_dir_fun} function that removes the directory
#'     paths from the filenames, in effect turning these into local files
#'   \item \link{get_replace_dir_fun} function that replaces the directory
#'     paths from the filenames
#' }
#' @param result results from measurements. These are:
#'   \itemize{
#'     \item log_evidence the natural logarithm of the evidence (aka marginal
#'       likelihood). Can be NA if this is not measured
#'     \item weight the weight of the model, compared to other (candidate)
#'       models. This weight will be between 0.0 (there is no evidence for
#'       this model) to 1.0 (all evidence indicates this is the best model).
#'       A weight of NA denotes that the weight is not measured
#'     \item errors a numeric vector of (positive) Bayesian inference errors.
#'       Will be NA if these are not measured.
#'   }
#' @param rng_seed a random number generator seed
#' @param rng_seeds a vector of random number generator seeds
#' @param rng_seed_twin_alignment the random number generator seed
#'   as used in the simulation of a twin alignment
#' @param rng_seed_twin_tree the random number generator seed as used in the
#'   simulation of a twin tree
#' @param root_sequence the DNA sequence at the root of the phylogeny.
#'   By default, this will consist out of an equal amount of each letter
#'   Use \link{check_root_sequence} to check if a root sequence is valid.
#' @param run_if the condition for an experiment's inference model to be run.
#'   Possible values:
#'   \itemize{
#'     \item \code{always}: always
#'     \item \code{best_candidate}: if the inference model is the
#'       candidate model with the most evidence (aka highest marginal
#'       likelihood)
#'   }
#' @param run_experiment one \link{pirouette} run experiment.
#'   A run experiment has these attributes:
#'   \itemize{
#'     \item experiment the (original) experiment
#'     \item true_result the result of running the original experiment on
#'       the true phylogeny
#'     \item twin_result the result of running the original experiment on
#'       the twin phylogeny
#'   }
#' @param run_experiments a list of one or more \link{pirouette} run experiments
#' @param sample_interval the interval at which the MCMC algorithm
#'   makes a measurement
#' @param sequence_length the length of each DNA sequence in an alignment
#' @param seed a random number generator seed
#' @param sim_phylo_fun function that, each time when called,
#' simulates one random tree.
#' @param sim_tral_fun function to simulate a
#' true alignment with.
#' This function must have two arguments,
#' called \code{true_phylogeny} (which will hold the true phylogeny)
#' and \code{root_sequence} (which holds the DNA root sequence).
#' The return type must be \link[ape]{DNAbin}.
#'
#' Use \link{check_sim_tral_fun} to verify if the function
#' has the right signature and output.
#'
#' Some standard functions:\cr
#' \itemize{
#'   \item Use \link{get_sim_tral_with_std_nsm_fun}
#'   to get a function (\link{sim_tral_with_std_nsm})
#'   the use a standard site model.
#'   \item Use
#'   \link{get_sim_tral_with_lns_nsm_fun}
#'   to get a function
#'   (\link{sim_tral_with_lns_nsm})
#'   the use a linked node substitution site model.
#'   \item Use
#'   \link{get_sim_tral_with_uns_nsm_fun}
#'   to get a function
#'   (\link{sim_tral_with_uns_nsm})
#'   the use an unlinked node substitution site model.
#' }
#' @param sim_twal_fun function to simulate a
#' twin alignment with.
#' This function must have two arguments called \code{twin_phylogeny} (which
#' will hold the twin phylogeny) and \code{true_alignment} (which will
#' hold the alignment simulated from the true phylogeny). The
#' return type must be \link[ape]{DNAbin}.
#'
#' Use \link{check_sim_twal_fun} to verify if the function
#' has the right signature and output.
#'
#' Some standard functions:\cr
#' \itemize{
#'   \item Use \link{get_copy_tral_fun}
#'     to get a function
#'     (\link{copy_true_alignment})
#'     that copies a true to alignment to create a twin alignment
#'   \item Use \link{get_sim_twal_with_std_nsm_fun}
#'     to get a function
#'     (\link{sim_twal_with_std_nsm})
#'     that simulates a twin alignment using a standard site model
#'   \item Use \link{get_sim_twal_same_n_muts_fun}
#'     to get a function
#'     (\link{sim_twal_with_same_n_mutation})
#'     that simulates -using a standard model- a twin alignment with as much
#'     mutations compared to the root sequence as the true alignment has
#'   \item Use \link{sim_twal_with_lns_nsm}
#'     that simulates a twin alignment using a linked node substitution
#'     model
#'   \item Use \link{sim_twal_with_uns_nsm}
#'     that simulates a twin alignment using an unlinked node substitution
#'     model
#' }
#' @param sim_twin_tree_fun function to simulate a twin tree with.
#' This function must have one argument called \code{phylogeny}
#' of type \link[ape]{phylo} and have a return type of type \link[ape]{phylo}
#' as well.
#'
#' Some standard functions:\cr
#' \itemize{
#'   \item Use \link{create_sim_yule_twin_tree_fun} to use a
#'     Yule (aka Pure Birth) process
#'   \item Use \link{create_copy_twtr_from_true_fun} to for a function
#'     that copies the true tree
#'   \item Use \link{get_sim_bd_twin_tree_fun} to use a
#'     Birth-Death process
#' }
#' @param site_model a nucleotide substitution model,
#'   which can be:
#'   \itemize{
#'     \item{
#'       A standard nucloetide substitution model,
#'       as created by \link[beautier]{create_site_model}
#'     }
#'     \item{
#'       \code{lns}: a linked node-substitution model
#'     }
#'     \item{
#'       \code{uns}: an unlinked node-substitution model
#'     }
#'   }
#' @param site_models a list of one or more site models,
#'   as created by \link[beautier]{create_site_model}
#' @param site_model_name name of a site model
#' @param subst_matrix nucleotide substitution matrix
#' @param sub_chain_length length of the sub-chain used by the Nested Sampling
#'   algorithm to estimate the marginal likelihood
#' @param tree an ultrametric phylogenetic tree of class \link[ape]{phylo}
#' @param tree_and_model one combination of a tree and model,
#'   as created by \link{get_tree_and_model_values}
#' @param tree_and_models one or more combination of a tree and model,
#'   as created by \link{get_tree_and_model_values}
#' @param tree_and_model_descriptions tabular data that maps
#'   a \code{tree_and_model} (e.g. \code{generative_true}) to
#'   a description (e.g. "Generative, true"),
#'   as created by \link{get_tree_and_model_descriptions}
#' @param tree_and_model_errors a tibble of a \code{tree_and_model}
#'   and errors, which passes \link{check_tree_and_model_errors}
#' @param treelog_filename name of the MCMC's treelog file,
#'   which is \code{$(tree).trees} by default.
#'   Use \link{complete_treelog_filename} to obtain the complete path to
#'   the MCMC's treelog file.
#' @param tree_model model used to simulate the tree
#' @param tree_prior a tree prior,
#'   as created by \link[beautier]{create_tree_prior}
#' @param tree_priors a list of one or more tree priors,
#'   as created by \link[beautier]{create_tree_prior}
#' @param tree_prior_name name of a tree prior
#' @param tree_type type of tree, can be \code{true} for the true
#'   phylogeny, and \code{twin} for its twin tree
#' @param tree_types types of tree, a vector of \code{true} for a true
#'   phylogeny, and \code{twin} for a twin tree
#' @param tree_filename name of the phylogeny file
#' @param true_alignment a DNA alignment, of class \link[ape]{DNAbin}
#' @param true_phylogeny the true phylogeny; the actual evolutionary
#' history of the species, of class \link[ape]{phylo}
#' @param true_result result obtained from using the true tree
#' @param twin_alignment a DNA alignment, of class \link[ape]{DNAbin}
#' @param twin_alignment_filename name of the FASTA file the twin
#'   alignment will be saved to
#' @param twin_evidence_filename filename to store the estimated
#'   evidences (aka marginal likelihoods) of the twin tree
#' @param twin_model the model you want to use to generate the twin tree:
#'   \itemize{
#'     \item \code{birth_death}: birth death
#'     \item \code{yule}: Yule or pure-birth
#'     \item \code{copy_true}: use a copy of the true tree in the twinning
#'      pipeline
#'   }
#'   See \link{get_twin_models} to see all possible
#'   values of \code{twin_model}
#' @param twin_phylogeny a phylogeny of class \link[ape]{phylo}
#' @param twin_result result obtained from using the twin tree
#' @param twin_tree_filename  name of the (\code{.newick}) file the twin
#'   tree will be saved to
#' @param twinning_params can be \code{NA} if no twinning is desired,
#'   or can be the twinning parameters,
#'   as can be created by \link{create_twinning_params}
#' @param type one or more ways to select the models used in inference:
#'   \itemize{
#'     \item \code{"generative"}: pick the generative model
#'     \item \code{most_evidence} picks the model with most evidence
#'   }
#'   See \link{get_model_selections} for a list.
#' @param verbose if TRUE, show more output
#' @param weight the weight of the model, compared to other (candidate)
#'   models. This weight will be between 0.0 (there is no evidence for
#'   this model) to 1.0 (all evidence indicates this is the best model).
#'   A weight of NA denotes that the weight is not measured
#' @author Documentation by Giovanni Laudanno,
#'   use of this function by Richèl J.C. Bilderbeek
#' @note This is an internal function, so it should be marked with
#'   \code{@noRd}. This is not done, as this will disallow all
#'   functions to find the documentation parameters
default_params_doc <- function(
  alignment,
  alignment_params,
  alignment_rng_seed,
  base_frequencies,
  bd_mutation_rate,
  bd_tree,
  bd_tree_filename,
  beast2_bin_path,
  beast2_input_filename,
  beast2_options,
  beast2_optionses,
  beast2_options_inference,
  beast2_options_est_evidence,
  beast2_output_log_filename,
  beast2_output_state_filename,
  beast2_output_trees_filename,
  beast2_output_trees_filenames,
  beast2_path,
  beast2_rng_seed,
  branch_mutation_rate,
  branch_subst_matrix,
  brts,
  burn_in_fraction,
  chain_length,
  check_input,
  clock_model, clock_models,
  clock_model_name,
  consensus,
  crown_age,
  df_long,
  do_measure_evidence,
  epsilon,
  error_fun,
  error_measure_params,
  errors,
  errors_filename,
  est_evidence_mcmc,
  evidence_epsilon,
  evidence_filename,
  exclude_model,
  experiment, experiments,
  extinction_rate,
  fasta_filename,
  filename,
  folder_name,
  folder_names,
  inference_model,
  inference_conditions,
  init_speciation_rate,
  init_extinction_rate,
  lambda,
  log_evidence,
  marg_lik_filename,
  marg_liks,
  max_evidence_epsilon,
  max_n_tries,
  mbd_l_matrix,
  mbd_mutation_rate,
  mbd_tree,
  mcmc,
  method,
  model_selection,
  model_type,
  mrca_prior,
  mu,
  mutation_rate,
  n_0,
  n_mutations,
  n_taxa,
  n_replicates,
  node_mutation_rate,
  node_subst_matrix,
  node_time,
  nu,
  nu_events,
  os,
  parameter_filename,
  parameters_filename,
  phylo,
  phylogenies,
  phylogeny,
  pir_params,
  pir_paramses,
  pir_out,
  pir_outs,
  posterior_trees,
  precision,
  project_folder_name,
  rename_fun,
  result,
  rng_seed, rng_seeds,
  rng_seed_twin_alignment,
  rng_seed_twin_tree,
  root_sequence,
  run_experiment,
  run_experiments,
  run_if,
  sample_interval,
  seed,
  sequence_length,
  sim_phylo_fun,
  sim_tral_fun,
  sim_twal_fun,
  sim_twin_tree_fun,
  site_model,
  site_models,
  site_model_name,
  sub_chain_length,
  subst_matrix,
  tree,
  tree_and_model, tree_and_models,
  tree_and_model_descriptions,
  tree_and_model_errors,
  treelog_filename,
  tree_filename,
  tree_model,
  tree_prior, tree_priors,
  tree_prior_name,
  tree_type,
  tree_types,
  true_alignment,
  true_phylogeny,
  true_result,
  twin_alignment,
  twin_alignment_filename,
  twin_evidence_filename,
  twin_phylogeny,
  twin_model,
  twin_result,
  twin_tree_filename,
  twinning_params,
  type,
  verbose,
  weight
) {
  # Nothing
}

## ---- message = FALSE---------------------------------------------------------
library(pirouette)

## ---- message = FALSE---------------------------------------------------------
library(pryr)

## -----------------------------------------------------------------------------
true_phylogeny <- ape::read.tree(text = "((A:1, B:1):1, C:2);")
ape::plot.phylo(true_phylogeny, main = "True phylogeny")

## -----------------------------------------------------------------------------
twin_tree <- pirouette::create_twin_tree(true_phylogeny)
ape::plot.phylo(twin_tree, main = "Twin phylogeny")

## -----------------------------------------------------------------------------
twinning_params <- pirouette::create_twinning_params()

## -----------------------------------------------------------------------------
twin_tree <- pirouette::create_twin_tree(
  phylogeny = true_phylogeny,
  twinning_params = twinning_params
)
ape::plot.phylo(twin_tree, main = "Twin phylogeny")

## -----------------------------------------------------------------------------
twinning_params <- pirouette::create_twinning_params(
  sim_twin_tree_fun = pirouette::get_sim_bd_twin_tree_fun()
)

## -----------------------------------------------------------------------------
twin_tree <- pirouette::create_twin_tree(
  phylogeny = true_phylogeny,
  twinning_params = twinning_params
)
ape::plot.phylo(twin_tree, main = "Twin phylogeny")

## -----------------------------------------------------------------------------
# Define my function
my_fun <- function(true_phylogeny) {
  new_phylo <- ape::rcoal(n = ape::Ntip(true_phylogeny))
  new_phylo$tip.label <- true_phylogeny$tip.label # nolint ape style, not mine
  new_phylo
}
# Put my function in the twinning_params
twinning_params <- pirouette::create_twinning_params(
  sim_twin_tree_fun = my_fun
)
# Create a twin tree using my function
twin_tree <- pirouette::create_twin_tree(
  phylogeny = true_phylogeny,
  twinning_params = twinning_params
)
# Show the twin tree created by my function
ape::plot.phylo(twin_tree, main = "Twin phylogeny")

## -----------------------------------------------------------------------------
head(sim_bd_twin_tree)

## -----------------------------------------------------------------------------
tryCatch(
  pirouette::create_twinning_params(
    sim_twin_tree_fun = sim_bd_twin_tree
  ),
  error = function(e) {
    cat(e$message)
  }
)

## -----------------------------------------------------------------------------
# Create my partially evaluated function
sim_twin_tree_fun <- pryr::partial(
  sim_bd_twin_tree,
  method = "random_tree",
  n_replicates = 1
)
# Create twinning_params with my function
twinning_params <- pirouette::create_twinning_params(
  sim_twin_tree_fun = sim_twin_tree_fun
)
# Create a twin tree using my function
twin_tree <- pirouette::create_twin_tree(
  phylogeny = true_phylogeny,
  twinning_params = twinning_params
)
# Show the twin tree
ape::plot.phylo(twin_tree, main = "Twin phylogeny")

## -----------------------------------------------------------------------------
# Create twinning_params with my function
twinning_params <- pirouette::create_twinning_params(
  sim_twin_tree_fun = pirouette::get_sim_bd_twin_tree_fun(
    method = "random_tree",
    n_replicates = 1
  )
)
# Create a twin tree using my function
twin_tree <- pirouette::create_twin_tree(
  phylogeny = true_phylogeny,
  twinning_params = twinning_params
)
# Show the twin tree
ape::plot.phylo(twin_tree, main = "Twin phylogeny")

## -----------------------------------------------------------------------------
# Create twinning_params with my function
twinning_params <- pirouette::create_twinning_params(
  sim_twin_tree_fun = pirouette::create_sim_yule_twin_tree_fun(
    method = "random_tree",
    n_replicates = 1
  )
)
# Create a twin tree using my function
twin_tree <- create_twin_tree(
  phylogeny = true_phylogeny,
  twinning_params = twinning_params
)
# Show the twin tree
ape::plot.phylo(twin_tree, main = "Twin phylogeny")

## -----------------------------------------------------------------------------
# Create twinning_params
twinning_params <- pirouette::create_twinning_params(
  sim_twin_tree_fun = pirouette::create_copy_twtr_from_true_fun()
)
# Create a twin tree using my function
twin_tree <- create_twin_tree(
  phylogeny = true_phylogeny,
  twinning_params = twinning_params
)
# Show the twin tree
ape::plot.phylo(twin_tree, main = "Twin phylogeny")

## -----------------------------------------------------------------------------
true_phylogeny <- ape::read.tree(text = "((A:1, B:1):1, C:2);")
ape::plot.phylo(true_phylogeny, main = "True phylogeny")

## -----------------------------------------------------------------------------
root_sequence <- pirouette::create_blocked_dna(length = 16)

## -----------------------------------------------------------------------------
alignment_params <- pirouette::create_alignment_params(
  root_sequence = root_sequence
)
true_alignment <- pirouette::sim_true_alignment(
  true_phylogeny = true_phylogeny,
  alignment_params = alignment_params
)
ape::image.DNAbin(true_alignment, main = "True alignment", legend = FALSE)

## -----------------------------------------------------------------------------
alignment_params <- pirouette::create_alignment_params(
  sim_tral_fun = sim_tral_with_std_nsm,
  root_sequence = root_sequence
)
true_alignment <- pirouette::sim_true_alignment(
  true_phylogeny = true_phylogeny,
  alignment_params = alignment_params
)
ape::image.DNAbin(true_alignment, main = "True alignment", legend = FALSE)

## -----------------------------------------------------------------------------
# My function
my_fun <- function(
  true_phylogeny,
  root_sequence
) {
  sequences <- list()
  for (i in seq_len(ape::Ntip(true_phylogeny))) {
    sequences[[i]] <- rep(
      sample(c("a", "c", "g", "t"), size = 1),
      nchar(root_sequence)
    )
  }
  ape::as.DNAbin(sequences)
}
# Putting my function in the alignment_params
alignment_params <- pirouette::create_alignment_params(
  sim_tral_fun = my_fun,
  root_sequence = root_sequence
)
# Simulate the true alignment using my function
true_alignment <- pirouette::sim_true_alignment(
  true_phylogeny = true_phylogeny,
  alignment_params = alignment_params
)
# Show the true alignment
ape::image.DNAbin(true_alignment, main = "True alignment", legend = FALSE)

## -----------------------------------------------------------------------------
head(sim_tral_with_std_nsm)

## -----------------------------------------------------------------------------
sim_tral_fun <- pryr::partial(
  sim_tral_with_std_nsm,
  mutation_rate = 0.5,
  site_model = beautier::create_hky_site_model()
)
head(sim_tral_fun)

## -----------------------------------------------------------------------------
alignment_params <- pirouette::create_alignment_params(
  sim_tral_fun = sim_tral_fun,
  root_sequence = root_sequence
)
true_alignment <- pirouette::sim_true_alignment(
  true_phylogeny = true_phylogeny,
  alignment_params = alignment_params
)
ape::image.DNAbin(true_alignment, main = "True alignment", legend = FALSE)

## -----------------------------------------------------------------------------
alignment_params <- pirouette::create_alignment_params(
  sim_tral_fun =
    pirouette::get_sim_tral_with_std_nsm_fun(
    mutation_rate = 0.5,
    site_model = beautier::create_hky_site_model()
  ),
  root_sequence = root_sequence
)
true_alignment <- pirouette::sim_true_alignment(
  true_phylogeny = true_phylogeny,
  alignment_params = alignment_params
)
ape::image.DNAbin(true_alignment, main = "True alignment", legend = FALSE)

## -----------------------------------------------------------------------------
if (1 == 2) { # nodeSub is not yet on CRAN
  alignment_params <- pirouette::create_alignment_params(
    sim_tral_fun =
      get_sim_tral_with_lns_nsm_fun(
      branch_mutation_rate = 0.1,
      node_mutation_rate = 0.2
    ),
    root_sequence = root_sequence
  )
  true_alignment <- pirouette::sim_true_alignment(
    true_phylogeny = true_phylogeny,
    alignment_params = alignment_params
  )
  ape::image.DNAbin(true_alignment, main = "True alignment", legend = FALSE)
}

## -----------------------------------------------------------------------------
if (1 == 2) { # nodeSub is not yet on CRAN
  alignment_params <- pirouette::create_alignment_params(
    sim_tral_fun =
      pirouette::get_sim_tral_with_uns_nsm_fun(
        branch_mutation_rate = 1.0,
        node_mutation_rate = 2.0,
        node_time = 0.1
      ),
    root_sequence = root_sequence
  )
  true_alignment <- pirouette::sim_true_alignment(
    true_phylogeny = true_phylogeny,
    alignment_params = alignment_params
  )
  ape::image.DNAbin(true_alignment, main = "True alignment", legend = FALSE)
}

## -----------------------------------------------------------------------------
twin_phylogeny <- ape::read.tree(text = "((A:2, B:2):1, C:3);")
ape::plot.phylo(twin_phylogeny, main = "Twin phylogeny")

## -----------------------------------------------------------------------------
root_sequence <- pirouette::create_blocked_dna(length = 20)

## -----------------------------------------------------------------------------
true_alignment <- pirouette::get_test_alignment(
  n_taxa = ape::Ntip(true_phylogeny),
  sequence_length = nchar(root_sequence)
)
ape::image.DNAbin(true_alignment, main = "True alignment", legend = FALSE)

## -----------------------------------------------------------------------------
my_fun <- function(
  twin_phylogeny = "irrelevant",
  true_alignment,
  root_sequence = "irrelevant"
) {
  true_alignment
}

## -----------------------------------------------------------------------------
twinning_params <- pirouette::create_twinning_params(
  sim_twal_fun = my_fun
)

## -----------------------------------------------------------------------------
twin_alignment <- pirouette::sim_twin_alignment(
  twin_phylogeny = twin_phylogeny,
  true_alignment = true_alignment,
  alignment_params = pirouette::create_test_alignment_params(
    root_sequence = root_sequence
  ),
  twinning_params = twinning_params
)
ape::image.DNAbin(twin_alignment, main = "Twin alignment")

## -----------------------------------------------------------------------------
head(sim_twal_with_std_nsm)

## -----------------------------------------------------------------------------
sim_twin_align_fun <- pryr::partial(
  sim_twal_with_std_nsm,
  mutation_rate = 0.1
)
head(sim_twin_align_fun)
pirouette::check_sim_twal_fun(sim_twin_align_fun)

## -----------------------------------------------------------------------------
twin_alignment <- pirouette::sim_twin_alignment(
  twin_phylogeny = twin_phylogeny,
  true_alignment = true_alignment,
  alignment_params = pirouette::create_test_alignment_params(),
  twinning_params = pirouette::create_twinning_params(
    sim_twal_fun = sim_twin_align_fun
  )
)
ape::image.DNAbin(twin_alignment, main = "Twin alignment", legend = FALSE)

## -----------------------------------------------------------------------------
twin_alignment <- pirouette::sim_twin_alignment(
  twin_phylogeny = twin_phylogeny,
  true_alignment = true_alignment,
  alignment_params = pirouette::create_test_alignment_params(),
  twinning_params = pirouette::create_twinning_params(
    sim_twal_fun =
      pirouette::get_sim_twal_with_std_nsm_fun(
        mutation_rate = 0.1
      )
  )
)
ape::image.DNAbin(twin_alignment, main = "Twin alignment", legend = FALSE)

## -----------------------------------------------------------------------------
head(sim_twal_with_same_n_mutation)

## -----------------------------------------------------------------------------
sim_twin_align_fun <- pryr::partial(
  sim_twal_with_same_n_mutation,
  mutation_rate = 0.5
)
head(sim_twin_align_fun)
pirouette::check_sim_twal_fun(sim_twin_align_fun)

## -----------------------------------------------------------------------------
twin_alignment <- pirouette::sim_twin_alignment(
  twin_phylogeny = twin_phylogeny,
  true_alignment = true_alignment,
  alignment_params = pirouette::create_test_alignment_params(
    root_sequence = root_sequence
  ),
  twinning_params = pirouette::create_twinning_params(
    sim_twal_fun = sim_twin_align_fun
  )
)
ape::image.DNAbin(twin_alignment, main = "Twin alignment", legend = FALSE)

## -----------------------------------------------------------------------------
twin_alignment <- pirouette::sim_twin_alignment(
  twin_phylogeny = twin_phylogeny,
  true_alignment = true_alignment,
  alignment_params = pirouette::create_test_alignment_params(),
  twinning_params = pirouette::create_twinning_params(
    sim_twal_fun =
      pirouette::get_sim_twal_with_std_nsm_fun(
        mutation_rate = 0.1
      )
  )
)
ape::image.DNAbin(twin_alignment, main = "Twin alignment", legend = FALSE)


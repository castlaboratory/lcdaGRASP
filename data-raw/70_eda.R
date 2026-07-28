# data-raw/70_eda.R
#
# Tukey EDA fuel: replicate Q and H over (graph x algorithm x centrality x
# similarity) so the vignette can draw the boxplots / ECDFs / (Q,H) scatter
# the paper omits (it reports only means + significance letters).
#
# -> inst/extdata/eda_replicates.rds

if (!exists("save_dataset")) {
  hp <- Filter(file.exists, c("00_helpers.R", "data-raw/00_helpers.R",
                              "lcdaGRASP/data-raw/00_helpers.R"))
  source(hp[1])
}
suppressPackageStartupMessages(library(purrr))

set.seed(1)
graphs <- list(
  Karate   = igraph::make_graph("Zachary"),
  Dolphins = load_benchmark("dolphins"),
  Football = load_benchmark("football"),
  Books    = load_benchmark("polbooks"))

REPS <- 30L
grid <- tidyr::expand_grid(
  graph      = names(graphs),
  algorithm  = c("LCDA-GRASP 1", "LCDA-GR 1"),
  centrality = c("eigen", "closeness"),
  similarity = c("hpi", "dice", "jaccard"),
  rep        = seq_len(REPS))
cli::cli_h1("EDA replicates: {nrow(grid)} runs")

results <- grid |>
  dplyr::mutate(.row = dplyr::row_number()) |>
  purrr::pmap_dfr(function(graph, algorithm, centrality, similarity, rep, .row) {
    if (.row %% 100 == 0) cli::cli_alert_info("eda {(.row)}/{nrow(grid)}")
    g <- graphs[[graph]]
    t0 <- Sys.time()
    res <- if (algorithm == "LCDA-GRASP 1")
      lcda_grasp(g, alpha_c = 0.1, alpha_s = 0.3, B = 50,
                 centrality = centrality, similarity = similarity, seed = rep)
    else
      lcda_gr(g, variant = 1, B = 150, centrality = centrality, similarity = similarity, seed = rep)
    tibble::tibble(graph, algorithm, centrality, similarity, rep,
                   Q = res$best$Q, H = res$best$H,
                   secs = as.numeric(Sys.time() - t0, units = "secs"))
  })

meta <- make_meta(
  title = "Replicated (Q, H) for EDA (Tukey lens)",
  description = "30 replications per (graph x algorithm x centrality x similarity) of best Q and H, with timing. Real benchmark networks (Karate/Dolphins/Football/Books).",
  network_source = src_newman(c("karate", "dolphins", "football", "polbooks")),
  seed = 1, reps = REPS,
  limitations = c("PolBlogs omitted here for runtime; covered in the benchmark-reproduction dataset.",
                  "alpha fixed at the paper defaults (0.1, 0.3) for the fixed-GRASP arm."))
save_dataset("eda_replicates", results, meta)

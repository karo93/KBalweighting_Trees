# Select simulation setting and replication

dataset_names <- list.files(here::here("simulation", "data"))
dataset_name <- dataset_names[8] #j=1
dat1 <- readRDS(here::here("simulation", "data", dataset_name) )

# Prepare training sample in the same way as in the simulation analysis
Xcol_select     <- floor(((ncol(dat1) + 3) / 2))
Xtestcol_select <- Xcol_select + 1
datTrain <- dat1 %>%filter(sim_id == 1) %>%select(2:Xcol_select) 

X <- datTrain %>% select(X1:last_col()) %>% as.data.frame()
Z <- datTrain$T_ind

# IPW
w_ipw <- ifelse(Z == 1, 1/datTrain$ps, 1/(1 - datTrain$ps))

# kbal
kbalout_treated <- kbal::kbal(
  allx = X,
  sampled = Z,
  sampledinpop = TRUE,
  ebal.tol = 1e-6,
  b = ncol(X),
  printprogress = TRUE,
  linkernel = TRUE
)

kbalout_control <- kbal::kbal(
  allx = X,
  sampled = 1 - Z,
  sampledinpop = TRUE,
  ebal.tol = 1e-6,
  b = ncol(X),
  printprogress = TRUE,
  linkernel = TRUE
)

# KBAL (combine treated + control)
w_kbal <- rep(NA_real_, length(Z))
w_kbal[Z == 1] <- kbalout_treated$w[Z == 1]
w_kbal[Z == 0] <- kbalout_control$w[Z == 0]

bt_all <- cobalt::bal.tab(
  X,
  treat = Z,
  weights = list(IPW = w_ipw, KBAL = w_kbal), 
  method = "weighting",
  s.d.denom = "pooled",
  un = TRUE,
  quick = FALSE
)

p <- cobalt::love.plot(
  bt_all,
  abs = TRUE,
  thresholds = c(m = 0.1),
  var.order = "unadjusted",
  stars = "raw",
  ggplot = TRUE
)
p <- p +
  labs(
    x = "Absolute standardized mean difference",
    y = "Covariates",
    title=NULL
  )  +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
plot_name <- paste0("covariate_balance_",  sub("^sim_", "", tools::file_path_sans_ext(dataset_name)),".png")

ggsave(
  filename = here::here(
    "output",
    "simulation_plots",
    plot_name
    ),  
  plot     = p,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 300
)

#Plot Propensities

data1<- sim_s3 %>% filter(sim_id==1)
data2<- sim_s3_moderate_overlap %>% filter(sim_id==1)
data3<- sim_s3_no_overlap %>% filter(sim_id==1)

plot_ps <- function(df, title) {
  ggplot(df, aes(x = ps, fill = factor(T_ind))) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = c("#1f77b4", "#ff7f0e"),
                      name = "Treatment Group",
                      labels = c("Control", "Treated")) +
    labs(x = "True propensity score", y = "Density", title = title) +
    theme_minimal(base_size = 14)
}
# Create plots
p1 <- plot_ps(data1, "Scenario 1: Low Overlap")
p2 <- plot_ps(data2, "Scenario 2: Low Overlap & Unbalanced")
p3 <- plot_ps(data3, "Scenario 3: Low Overlap & Balanced")
df_all <- bind_rows(
  mutate(data1, Scenario = "Good overlap"),
  mutate(data2, Scenario = "Moderate overlap"),
  mutate(data3, Scenario = "Poor overlap")
) %>%
  mutate(
    Z = factor(T_ind, levels = c(0,1), labels = c("Control","Treatment")) 
  )

ggplot(df_all, aes(x = ps, fill = Z, colour = Z)) +
  geom_density(alpha = 0.4, linewidth = 0.8, adjust = 1) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +    
  scale_colour_manual(values = c("#1f77b4", "#ff7f0e")) +
  facet_wrap(~ Scenario, nrow = 1, scales = "fixed") +    
  labs(
    x = "PS",
    y = "density",
    fill = "Treatment Group",      # legend title
    colour = "Treatment Group",    # legend title
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",             # legend on right
    strip.text = element_blank(),           # no facet headlines
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0)  # caption below, left-aligned
  )
ggsave(here::here("output" , "simulation_plots" , "ps_overlap_s3_plot.png"), width=12, height=4, dpi=300)

# Now Application 
train_df <- IHDP_A_train[[2]]
test_df  <- IHDP_A_test[[2]]

Data_A_all <- dplyr::bind_rows(
  mutate(train_df, set = "train"),
  mutate(test_df,  set = "test")
)

Data_B_all <- IHDP_B_no_overlap %>% filter(iteration == 1)

X_all_B <- Data_B_all[, 2:26]
t_all_B <- Data_B_all$z
cvfit_B <- glmnet::cv.glmnet(as.matrix(X_all_B), t_all_B, family = "binomial")
best_lambda <- cvfit_B$lambda.min

# Predict propensity scores (probabilities of treatment = 1)
Propensity_Logit_B <- predict(cvfit_B,
                              newx = as.matrix(X_all_B),
                              s = best_lambda,
                              type = "response")

# Make it a simple numeric vector
Propensity_Logit_B <- as.vector(Propensity_Logit_B)
X_all_A <- Data_A_all[, 6:30]
t_all_A <- Data_A_all$treatment
cvfit_A <- glmnet::cv.glmnet(as.matrix(X_all_A), t_all_A, family = "binomial")
best_lambda <- cvfit_A$lambda.min

# Predict propensity scores (probabilities of treatment = 1)
Propensity_Logit_A <- predict(cvfit_A,
                              newx = as.matrix(X_all_A),
                              s = best_lambda,
                              type = "response")

# Make it a simple numeric vector
Propensity_Logit_A <- as.vector(Propensity_Logit_A)
df_all <- bind_rows(
  tibble(
    Propensity = Propensity_Logit_A,
    Treatment  = factor(t_all_A, levels = c(0,1), labels = c("Control","Treatment")),
    Scenario   = "Left: Dataset A"
  ),
  tibble(
    Propensity = Propensity_Logit_B,
    Treatment  = factor(t_all_B, levels = c(0,1), labels = c("Control","Treatment")),
    Scenario   = "Right: Dataset B"
  )
)

ggplot(df_all, aes(x = Propensity, fill = Treatment, colour = Treatment)) +
  geom_density(alpha = 0.4, linewidth = 0.8, adjust = 1) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e"))+
  scale_colour_manual(values = c("#1f77b4", "#ff7f0e")) +
  facet_wrap(~ Scenario, nrow = 1, scales = "fixed") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    x = "PS",
    y = "density",
    fill = "Treatment",
    colour = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    strip.text = element_blank(),
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0)
  )
ggsave(
  filename = here::here(
    "output",
    "application_plots",
    "ps_overlap_plot_application.png"
  ),
  width = 12,
  height = 4,
  dpi = 300
)

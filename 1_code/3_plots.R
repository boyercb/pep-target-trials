

# main plot ---------------------------------------------------------------

sims <- bind_rows(
  sim0,
  sim40,
  .id = 'scenario'
) |> 
  mutate(
    rr = if_else(str_detect(type, "_rr"), 1, 0),
    type = factor(
      type,
      levels = c(
        "leave",
        "move",
        "cox",
        "target_trial",
        "leave_rr",
        "move_rr",
        "target_trial_rr"
      ), 
      labels = c(
        "naive, leave",
        "naive, move",
        "time-varying cox",
        "target trial",
        "naive, leave",
        "naive, move",
        "target trial"
      )
    ),
    scenario = case_when(
      scenario == 1 ~ "Scenario 1: VE = 0%",
      scenario == 2 & rr == 1 ~ "Scenario 2: VE = 31.6%",
      scenario == 2 & rr == 0 ~ "Scenario 2: VE = 40%",
      scenario == 3 ~ "Scenario 3: time-dependent VE",
      scenario == 4 ~ "Scenario 4: confounding"
    ),
  )

p1 <- ggplot(
  filter(sims, rr == 0), 
  aes(x = type, y = VE, fill = type)
) +
  facet_grid(~scenario) +
  geom_hline(aes(yintercept = truth), linetype = 'dashed') +
  geom_boxplot() +
  scale_color_manual(values = brewer_pal(palette = "Blues")(5)[1:4]) +
  scale_fill_manual(values = brewer_pal(palette = "Blues")(5)[1:4]) +
  scale_y_continuous(labels = label_percent(), limits = c(-0.35, 0.65)) +
  scale_x_discrete(labels = scales::label_wrap(10)) +
  labs(
    x = NULL,
    y = "Estimated vaccine efficacy"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p2 <- ggplot(
  filter(sims, rr == 1), 
  aes(x = type, y = VE, fill = type)
) +
  facet_grid(~scenario) +
  geom_hline(aes(yintercept = truth), linetype = 'dashed') +
  geom_boxplot() +
  scale_color_manual(values = brewer_pal(palette = "Blues")(5)[c(1,2,4)]) +
  scale_fill_manual(values = brewer_pal(palette = "Blues")(5)[c(1,2,4)]) +
  scale_y_continuous(labels = label_percent(), limits = c(-0.35, 0.65)) +
  scale_x_discrete(labels = scales::label_wrap(10)) +
  labs(
    x = NULL,
    y = "Estimated vaccine efficacy"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(
  filename = "3_figures/sim_hr.pdf",
  plot = p1,
  device = "pdf",
  width = 6.5,
  height = 4
)

ggsave(
  filename = "3_figures/sim_rr.pdf",
  plot = p2,
  device = "pdf",
  width = 6.5,
  height = 4
)


# heterogeneity -----------------------------------------------------------

simV_df <- 
  simV |>
  mutate(
    rr = if_else(str_detect(type, "_rr"), 1, 0),
    type = factor(
      type,
      levels = c(
        "leave",
        "move",
        "cox",
        "target_trial",
        "leave_rr",
        "move_rr",
        "target_trial_rr"
      ), 
      labels = c(
        "naive, leave",
        "naive, move",
        "cox",
        "target trial",
        "naive, leave",
        "naive, move",
        "target trial"
      )
    ),
  )

p <- ggplot(filter(simV_df, rr == 0), aes(x = factor(X), y = VE, fill = type)) +
  facet_wrap(~type) +
  geom_boxplot() +
  geom_point(aes(y = truth), color = "red", fill = "red", shape = 23) +
  # geom_lin(aes(y = truth)) +
  scale_color_brewer(palette = "Blues") +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(labels = label_percent(), limits = c(-0.25, 1)) +
  labs(
    x = "Day",
    y = "Estimated vaccine efficacy"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")


ggsave(
  filename = "3_figures/sim_hetx.pdf",
  plot = p,
  device = "pdf",
  width = 6.5,
  height = 5
)


# multiple designs plot ---------------------------------------------------

simTT_df <- 
  simTT |>
  mutate(
    type = factor(
      type,
      levels = c("fixed", "grace", rep("fixed_tv", 7), rep("dayzero", 7)), 
      labels = c("fixed enrollment (pooled)", "grace period", rep("fixed enrollment", 7), rep("day-zero", 7))
    )
  )

p <- ggplot(simTT_df, aes(x = factor(X), y = VE, fill = type)) +
  geom_boxplot() +
  geom_point(aes(y = truth), color = "red", fill = "red", shape = 23) +
  # geom_lin(aes(y = truth)) +
  coord_flip() +
  scale_color_brewer(palette = "Blues") +
  facet_wrap(~ type, drop = TRUE, shrink = TRUE, scales = "free") +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(labels = label_percent(), limits = c(-0.25, 1)) +
  labs(
    x = "Day",
    y = "Estimated vaccine efficacy"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")


ggsave(
  filename = "3_figures/sim_TT.pdf",
  plot = p,
  device = "pdf",
  width = 6.5,
  height = 5
)

# overlap plot ------------------------------------------------------------

# sim_overlap$type <- rep(
#   c(
#     "leave",
#     "move",
#     "cox",
#     "target_trial",
#     "leave_rr",
#     "move_rr",
#     "target_trial_rr"
#   ),
#   3000
# )

sim_overlap_df <- 
  sim_overlap |> 
  mutate(
    rr = if_else(str_detect(type, "_rr"), 1, 0),
    type = factor(
      type,
      levels = c(
        "leave",
        "move",
        "cox",
        "target_trial",
        "leave_rr",
        "move_rr",
        "target_trial_rr"
      ), 
      labels = c(
        "naive, leave",
        "naive, move",
        "cox",
        "target trial",
        "naive, leave",
        "naive, move",
        "target trial"
      )
    ),
    # meanlog = exp(meanlog),
    flab = factor(
      paste0("mu == ", meanlog),
      levels = c(
        "mu == 3",
        "mu == 9",
        "mu == 15"
      )
    )
  )

p1 <- ggplot(
  filter(sim_overlap_df, rr == 1), 
  aes(x = type, y = VE, fill = type)
) +
  facet_wrap(~flab, labeller = label_parsed) +
  geom_hline(aes(yintercept = truth), linetype = 'dashed') +
  geom_boxplot() +
  scale_color_manual(name = "", guide = "none", values = brewer_pal(palette = "Blues")(5)[1:4]) +
  scale_fill_manual(name = "", guide = "none", values = brewer_pal(palette = "Blues")(5)[1:4]) +
  scale_y_continuous(labels = label_percent(), limits = c(-0.55, 0.55)) +
  scale_x_discrete(labels = scales::label_wrap(10)) +
  labs(
    x = NULL,
    y = "Estimated vaccine efficacy"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")


df <- tibble(
  flab = rep(c(
    "mu == 3",
    "mu == 9",
    "mu == 15"
  ), each = 50000),
  X = rpois(50000 * 3, 5),
  T = floor(
    rlnorm(50000 * 3, 
           rep(c(log(3), log(9), log(15)), each = 50000),
           0.59)
    )
)

df$flab <- 
  factor(
    df$flab,
    levels = c(
      "mu == 3",
      "mu == 9",
      "mu == 15"
    )
  )

p2 <- 
  ggplot(pivot_longer(df, -flab), aes(fill = name, x = value)) +
  facet_wrap(~flab, labeller = label_parsed) +
  geom_bar(
    aes(y = ..prop..),
    alpha = 0.45,
    position = position_identity(),
    color = "black",
    width = 1
  ) +
  xlim(c(-1, 22)) +
  scale_fill_manual(
    name = "",
    values = c("blue", "orange"),
    labels = c(
      "Symptom onset time (T*)", 
      "Vaccine administration time (X*)"
    )
  ) +
  theme_bw(base_size = 12) +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme(
    legend.position = "bottom", 
    legend.background = element_blank(),
    legend.spacing.y = unit(0.1, 'cm'),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Day",
    y = "Probability"
  )

ggsave(
  filename = "3_figures/sim_overlap.pdf",
  plot = p2 / p1 + plot_layout(guides = 'collect') &
    theme(legend.position='top'),
  device = "pdf",
  width = 6.5,
  height = 6.5
)
# 
# p3 <- 
#   ggplot(pivot_longer(df, -flab) |> filter(name == "T" & flab == "mu == 9"), aes(fill = name, x = value)) +
#   #facet_wrap(~flab, labeller = label_parsed) +
#   geom_bar(
#     aes(y = ..prop..),
#     alpha = 0.45,
#     position = position_identity(),
#     color = "black",
#     width = 1
#   ) +
#   xlim(c(-1, 22)) +
#   scale_fill_manual(
#     name = "",
#     values = c("blue", "orange"),
#     labels = c(
#       "Symptom onset time (T*)", 
#       "Vaccine administration time (X*)"
#     )
#   ) +
#   theme_bw(base_size = 12) +
#   guides(fill = guide_legend(byrow = TRUE)) +
#   theme(
#     legend.position = "none", 
#     legend.background = element_blank(),
#     legend.spacing.y = unit(0.1, 'cm'),
#     panel.grid = element_blank()
#   ) +
#   ylim(c(0, 0.17)) +
#   labs(
#     x = "day",
#     y = "probability"
#   )
# 
# ggsave(
#   filename = "3_figures/dist_2.pdf",
#   plot = p3,
#   device = "pdf",
#   width = 5,
#   height = 3.5
# )




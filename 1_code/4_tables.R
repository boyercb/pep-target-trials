files <- paste0("0_data/sim", 1:5, ".rds")

sims_list <- lapply(files, function(file) read_rds(file))

sims <- rbindlist(sims_list, idcol = "run")

sims <- sims[, sim := sim + 100 * (run - 1)]
sims <- 
  sims[
    str_detect(type, "_rr") & truth == log(1 - 0.4) & mu == log(3), truth := log(1 - calc_true_rr_cond(0.4, log(3), 0.59))
    ][
      str_detect(type, "_rr") & truth == log(1 - 0.4) & mu == log(9), truth := log(1 - calc_true_rr_cond(0.4, log(9), 0.59))
    ][
      str_detect(type, "_rr") & truth == log(1 - 0.4) & mu == log(15), truth := log(1 - calc_true_rr_cond(0.4, log(15), 0.59))
    ][
      str_detect(type, "_rr") & truth == log(1 - 0.8) & mu == log(3), truth := log(1 - calc_true_rr_cond(0.8, log(3), 0.59))
    ][
      str_detect(type, "_rr") & truth == log(1 - 0.8) & mu == log(9), truth := log(1 - calc_true_rr_cond(0.8, log(9), 0.59))
    ][
      str_detect(type, "_rr") & truth == log(1 - 0.8) & mu == log(15), truth := log(1 - calc_true_rr_cond(0.8, log(15), 0.59))
    ][]

sim_results <-
  sims[, .(est = mean(est),
           VE = mean(1 - exp(est)),
           bias = mean(est - mean(truth)),
           ese = mean(se),
           sd = sd(est),
           rmse = sqrt(mean((est - mean(truth))^2)),
           coverage = mean(
             (est + qnorm(0.025) * se) < truth & truth < (est + qnorm(0.975) * se)
           )
  ),
  by = list(mu, truth, type)]

tab_rr <- sim_results[str_detect(sim_results$type, "_rr"), ]

tab <- sim_results[!str_detect(sim_results$type, "_rr"), ]

tab |>
  arrange(desc(truth), mu) |>
  select(-c(mu, truth, rmse, VE)) |>
  mutate(
    type = case_when(
      type == "leave" ~ "naive, leave",
      type == "move" ~ "naive, move",
      type == "cox" ~ "time-varying cox",
      type == "target_trial" ~ "target trial"
    )
  ) |>
  kable(
    format = "latex",
    digits = 3,
    booktabs = TRUE,
    align = "lccccc",
    col.names = c(
      "Estimator",
      "Mean",
      "Bias",
      "ESE",
      "SE",
      "Coverage"
    ),
    linesep = ""
  ) |>
  kable_styling() |>
  group_rows(
    group_label = "VE = 0\\% and $\\\\mu = 3$",
    start_row = 1,
    end_row = 4,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 0\\% and $\\\\mu = 9$",
    start_row = 5,
    end_row = 8,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 0\\% and $\\\\mu = 15$",
    start_row = 9,
    end_row = 12,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 40\\% and $\\\\mu = 3$",
    start_row = 13,
    end_row = 16,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 40\\% and $\\\\mu = 9$",
    start_row = 17,
    end_row = 20,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 40\\% and $\\\\mu = 15$",
    start_row = 21,
    end_row = 24,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |> 
  group_rows(
    group_label = "VE = 80\\% and $\\\\mu = 3$",
    start_row = 25,
    end_row = 28,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 80\\% and $\\\\mu = 9$",
    start_row = 29,
    end_row = 32,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 80\\% and $\\\\mu = 15$",
    start_row = 33,
    end_row = 36,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  )


tab_rr |>
  arrange(desc(truth), mu) |>
  select(-c(mu, truth, rmse, VE)) |>
  mutate(
    type = case_when(
      type == "leave_rr" ~ "naive, leave",
      type == "move_rr" ~ "naive, move",
      type == "target_trial_rr" ~ "target trial"
    )
  ) |>
  kable(
    format = "latex",
    digits = 3,
    booktabs = TRUE,
    align = "lccccc",
    col.names = c(
      "Estimator",
      "Mean",
      "Bias",
      "ESE",
      "SE",
      "Coverage"
    ),
    linesep = ""
  ) |>
  kable_styling() |>
  group_rows(
    group_label = "VE = 0\\% and $\\\\mu = 3$",
    start_row = 1,
    end_row = 3,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 0\\% and $\\\\mu = 9$",
    start_row = 4,
    end_row = 6,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 0\\% and $\\\\mu = 15$",
    start_row = 7,
    end_row = 9,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 40\\% and $\\\\mu = 3$",
    start_row = 10,
    end_row = 12,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 40\\% and $\\\\mu = 9$",
    start_row = 13,
    end_row = 15,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 40\\% and $\\\\mu = 15$",
    start_row = 16,
    end_row = 18,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |> 
  group_rows(
    group_label = "VE = 80\\% and $\\\\mu = 3$",
    start_row = 19,
    end_row = 21,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 80\\% and $\\\\mu = 9$",
    start_row = 22,
    end_row = 24,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  group_rows(
    group_label = "VE = 80\\% and $\\\\mu = 15$",
    start_row = 25,
    end_row = 27,
    escape = FALSE,
    bold = FALSE,
    italic = TRUE
  ) |>
  footnote(
    general = "Test",
    threeparttable = TRUE)



sim_overlap_df <- sims |>
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
      paste0("mu == ", exp(mu)),
      levels = c(
        "mu == 3",
        "mu == 9",
        "mu == 15"
      )
    ),
    truth.fct = factor(truth, labels = c("80", "40", "0")),
    truth = pmap_dbl(list(truth, mu), function(truth, mu) calc_true_rr_cond(1 - exp(truth), mu, 0.59)),
    VE = 1 - exp(est)
  )
  

p1 <- ggplot(
  filter(sim_overlap_df, rr == 1 & truth.fct == "0"), 
  aes(x = type, y = VE, fill = type)
) +
  facet_wrap(~flab, labeller = label_parsed) +
  geom_hline(aes(yintercept = truth), linetype = 'dashed') +
  geom_boxplot() +
  scale_color_manual(name = "", guide = "none", values = brewer_pal(palette = "Blues")(5)[c(1, 2, 4)]) +
  scale_fill_manual(name = "", guide = "none", values = brewer_pal(palette = "Blues")(5)[c(1, 2, 4)]) +
  scale_y_continuous(labels = label_percent(), limits = c(-0.4, 0.4)) +
  scale_x_discrete(labels = scales::label_wrap(10)) +
  labs(
    x = NULL,
    y = "Estimated vaccine efficacy"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
p1


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
    aes(y = after_stat(prop)),
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
  filename = "3_figures/sim_overlap_rr.pdf",
  plot = p2 / p1 + plot_layout(guides = 'collect') &
    theme(legend.position='top'),
  device = "pdf",
  width = 6.5,
  height = 6.5
)
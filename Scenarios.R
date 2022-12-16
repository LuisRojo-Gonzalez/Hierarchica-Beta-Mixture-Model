# This script constructs the bounds for the model
# based on the literature review

# Working directory
setwd("~/Desktop/ObesityChile/BetaMixture/AlioEuro")

# libraries
library(tikzDevice)
options(mc.cores = parallel::detectCores(),
        tz="CA")
library(haven)
library(tidyverse)
library(xtable)
library(xlsx)
library(readxl)
library(ggpubr)
library(distr)
# library(ggtern)
# library(ggalt)
library(maxLik)
library(sirt) # dirichlet.mle

# citation(package = "sirt")

# Transition probabilities, feasible bounds

# Men
LitBoundsMen = read_excel("LitBounds.xlsx", 
                          range = "D3:F24")
LitBoundsMen$NS = c(rep("Underweight", 3),
                    rep("Normal", 5),
                    rep("Overweight", 7),
                    rep("Obese", 6))
LitBoundsMen$Sex = "Men"

# Women
LitBoundsWomen = read_excel("LitBounds.xlsx", 
                            range = "G3:I24")
LitBoundsWomen$NS = c(rep("Underweight", 3),
                      rep("Normal", 5),
                      rep("Overweight", 7),
                      rep("Obese", 6))
LitBoundsWomen$Sex = "Women"

# Both
LitBounds = rbind(LitBoundsMen, LitBoundsWomen)
rm(LitBoundsMen, LitBoundsWomen)

LitBounds = LitBounds %>%
  mutate(NS = factor(NS, levels = c("Underweight", "Normal",
                                    "Overweight", "Obese")),
         Sex = factor(Sex, levels = c("Men", "Women")))

# plot density curves
LitBounds %>%
  gather(key = "Transition", value = "value",
         -NS, -Sex, factor_key = TRUE) %>%
  mutate(Transition = factor(Transition,
                             levels = c("Decrease", "Remain", "Increase"),
                             labels = c(expression("To decrease ("*beta*")"),
                                        expression("To remain ("*phi*")"),
                                        expression("To increase ("*alpha*")"))),
         NS = factor(NS, levels = c("Underweight", "Normal",
                                    "Overweight", "Obese"),
                     labels = c("Underweight", "Normal",
                                "Overweight", "Obese")),
         value = as.numeric(value)) %>%
  ggplot() +
  geom_density(aes(x = value, fill = NS,
                   colour = NS),
               alpha = 0.3) +
  facet_grid(Sex ~ Transition,
             scales = "free_y",
             space = "fixed",
             switch = "y",
             labeller = label_parsed) +
  theme_bw() +
  theme(plot.title = element_text(colour = "black",
                                  size = 10),
        axis.text.x = element_text(angle = 0, # whether to show the lines on
                                   hjust = 0.5, # the axis
                                   vjust = 0.5,
                                   colour = "black",
                                   size = 10),
        axis.text.y = element_text(angle = 0, # whether to show the lines on
                                   hjust = 0.5, # the axis
                                   vjust = 0.5,
                                   colour = "black",
                                   size = 10),
        # axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank(), # to control the appearance
        panel.spacing = unit(1, "lines"), # space around the panel in grid
        legend.margin = margin(2, 2, 2, 2), # space around the legend box
        legend.position = "bottom",     # of the plot
        # legend.key.width = unit(3.9, "cm"), # length of the legend
        legend.box = "vertical",
        legend.spacing = unit(0.2, "lines"),
        legend.box.margin = margin(-8, 0, 0, 0), #top, right, bottom, left
        legend.background = element_rect(fill = "white",
                                         size = 0.5,
                                         linetype = "solid",
                                         colour = "black"),
        legend.text = element_text(colour = "black",
                                   size = 10),
        legend.title = element_text(colour = "black",
                                    size = 10),
        strip.text.x = element_text(colour = "black",
                                    size = 10),
        strip.text.y = element_text(colour = "black",
                                    size = 10),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + #top, right, bottom, left
  labs(x = "Transition probability",
       y = "density",
       fill = "Nutritional status",
       colour = "Nutritional status") +
  scale_x_continuous(labels = scales::percent)

# summarizes the probabilities by combination of factor
resBounds = LitBounds  %>%
  # replace(is.na(.), 0) %>%
  gather(key = "Transition", value = "value",
         -NS, -Sex, factor_key = TRUE) %>%
  unite(temp, Transition, Sex) %>%
  group_by(NS, temp) %>%
  summarise(Min = round(min(value, na.rm = TRUE), 2),
            Q1 = round(quantile(value, probs = 0.25, na.rm = TRUE), 2),
            Q2 = round(quantile(value, probs = 0.5, na.rm = TRUE), 2),
            Mean = round(mean(value, na.rm = TRUE), 2),
            Q3 = round(quantile(value, probs = 0.75, na.rm = TRUE), 2),
            Max = round(max(value, na.rm = TRUE), 2),
            SD = round(sd(value, na.rm = TRUE), 2),
            Zeros = sum(is.na(value))) %>%
  ungroup() %>%
  mutate(temp = factor(temp, levels = c("Decrease_Men", "Remain_Men",
                                        "Increase_Men", "Decrease_Women",
                                        "Remain_Women", "Increase_Women"))) %>%
  gather(key = "Variable", value = "value", -NS, -temp) %>%
  mutate(Variable = factor(Variable, levels = c("Min", "Q1", "Q2",
                                                "Mean", "Q3", "Max", "SD", "Zeros"))) %>%
  filter(!is.infinite(value)) %>%
  spread(temp, value, fill = "-")

# print table in latex format
print(xtable(resBounds,
             digits = 2, label = "tab:boundssummary",
             caption = "Descriptive statistics for transition probabilities between nutritional statuses found in the leterature."),
      caption.placement = "top",
      comment = FALSE,
      include.rownames = FALSE,
      include.colnames = TRUE,
      hline.after = c(-1, 0, 7, 14, 21, 28))

# fit beta distribution from the bounds for each type of transition
# and combination of factors
beta_margin_params = data.frame(Sex = c("Men", "Women"))

for (i in 1:nrow(beta_margin_params)) {
  x = LitBounds %>%
    replace_na(list(Decrease = 0,
                    Increase = 0)) %>%
    filter(Sex == beta_margin_params$Sex[i]) %>%
    dplyr::select(-c("NS", "Sex")) %>%
    as.matrix()
  # count = apply(aux, 2, function(x) sum(x != 0)) %>%
  #   as.vector()
  beta_margin_params[i, 2:14] = x %>%
    dirichlet.mle(x = ., progress = FALSE, maxit = 5000) %>%
    unlist() %>%
    t() %>%
    as_tibble() %>%
    mutate(shape1_Decrease = exp(.[[1]]),
           shape1_Remain = exp(.[[2]]),
           shape1_Increase = exp(.[[3]]),
           shape2_Decrease = exp(.[[2]]) * exp(.[[3]]),
           shape2_Remain = exp(.[[1]]) * exp(.[[3]]),
           shape2_Increase = exp(.[[1]]) * exp(.[[2]]))
  # x ~ beta(alpha, alpha_0 - alpha)
}

beta_margin_params = beta_margin_params[, c(1, 9:14)] %>%
  gather("Parameter", "value", -Sex) %>%
  separate(., col = "Parameter", into = c("shape", "Transition")) %>%
  spread(shape, value)

# simulating dirichlet random numbers
set.seed(12345)
n_sim = 1000
sample_dir = data.frame()
for (i in 1:nrow(beta_margin_params)) {
  sample_dir = rbind(sample_dir,
                     data.frame(Sex = beta_margin_params$Sex[i],
                                Transition = beta_margin_params$Transition[i],
                                prob = rbeta(n_sim, shape1 = beta_margin_params$shape1[i],
                                             shape2 = beta_margin_params$shape2[i])))
}

sample_dir %>%
  ggplot() +
  geom_density(aes(x = prob)) +
  facet_grid(Sex ~ Transition)

# EM algorithm initialization

# step 1: create the grid search
control = resBounds %>%
  filter(Variable %in% c("Mean", "SD")) %>%
  gather("combination", "value",
         -NS, -Variable) %>%
  separate(., combination, into = c("Transition", "Sex")) %>%
  filter(value != "-") %>%
  spread(Variable, value) %>%
  # select(c("NS", "Transition", "Sex")) %>%
  data.frame()

weights = resBounds %>%
  filter(Variable == "Zeros") %>%
  dplyr::select(-c("Variable")) %>%
  gather(key = "temp", value, -NS) %>%
  separate(temp, into = c("Transition", "Sex")) %>%
  mutate(value = as.numeric(value)) %>%
  group_by(Transition, Sex) %>%
  mutate(value = value + 1,
         weight_aux = (sum(value)-value)/sum(sum(value)-value),
         weight = ifelse(is.na(weight_aux), 0.25, weight_aux)) %>%
  ungroup() %>%
  dplyr::select(c("Sex", "Transition", "weight"))

grid = expand_grid(sex = c("Men", "Women"),
                   tran = c("Decrease", "Remain", "Increase"),
                   ns = c("Underweight", "Normal",
                          "Overweight", "Obese"))

grid_2 = expand_grid(sex = c("Men", "Women"),
                     tran = c("Decrease", "Remain", "Increase"),
                     algo = c("Quasi", "NM", "SANN"),
                     w_inverse = c("No", "Yes"))

# saves the file required by the EM algorithm
saveRDS(sample_dir, file = "sample_dir.rds")
saveRDS(control, file = "control.rds")
saveRDS(grid, file = "grid.rds")
saveRDS(grid_2, file = "grid_2.rds")
saveRDS(beta_margin_params, file = "beta_margin.rds")
saveRDS(weights, file = "weights.rds")

# reading the EM algorithm output files from NLHPC
# em = list.files(path = "~/Desktop/ObesityChile/BetaMixture/",
#                 pattern = "Par_") %>%
#   map_dfr(readRDS)
# 
# res = list.files(path = "/BetaMixture/",
#                  pattern = "Res_") %>%
#   map_dfr(readRDS)

em_params_aux = data.frame()
res = data.frame()
for (i in 1:nrow(grid_2)) {
  em_params_aux = rbind(em_params_aux,
                        readRDS(paste0("Experiments/Par_",
                                       grid_2$sex[i], "_",
                                       grid_2$tran[i], "_",
                                       as.numeric(i),".rds")))
  
  res = rbind(res, readRDS(paste0("Experiments/Res_",
                                  grid_2$sex[i], "_",
                                  grid_2$tran[i], "_",
                                  as.numeric(i), ".rds")) %>%
                mutate(sex = grid_2$sex[i],
                       tran = grid_2$tran[i]))
}

# complete serie to assess the evolution of the algorithm
a_aux_1_serie = res %>%
  left_join(., beta_margin_params %>%
              rename(sex = Sex, tran = Transition),
            by = c("sex", "tran"))

for (i in c("Men", "Women")) {
  for (j in c("Decrease", "Remain", "Increase")) {
   assign(paste0("plot_", i, "_", j),
          a_aux_1_serie %>%
            dplyr::select(c("id", "Iteration", "Loglik",
                            "area", "Method", "Time", "w_inverse",
                            "sex", "tran")) %>%
            filter(Iteration != 0,
                   sex == i,
                   tran == j) %>%
            ggplot(aes(x = Time,
                       y = Loglik,
                       group = id)) +
            geom_point(aes(colour = Method,
                           shape = w_inverse)) +
            geom_line(aes(colour = Method,
                          linetype = w_inverse)) +
            geom_point(data = a_aux_1_serie %>%
                         dplyr::select(c("id", "Iteration", "Loglik",
                                         "area", "Method", "Time", "w_inverse",
                                         "sex", "tran")) %>%
                         filter(Iteration != 0,
                                sex == i,
                                tran == j) %>%
                         group_by(id) %>%
                         filter(Iteration == max(Iteration)) %>%
                         ungroup() %>%
                         filter(area == max(area)),
                       aes(x = Time, y = Loglik),
                       shape = 1,
                       size = 5,
                       show.legend = FALSE) +
            theme(plot.title = element_text(#family = "Times",
              #size = 10,
              colour = "black"),
              axis.text.x = element_text(angle = 0, # whether to show the lines on
                                         hjust = 0.5, # the axis
                                         # family = "Times",
                                         # size = 10,
                                         vjust = 0.5,
                                         colour = "black"),
              axis.text.y = element_text(angle = 0, # whether to show the lines on
                                         hjust = 0.5, # the axis
                                         # family = "Times",
                                         # size = 10,
                                         vjust = 0.5,
                                         colour = "black"),
              axis.title = element_text(angle = 0, # whether to show the lines on
                                        hjust = 0.5, # the axis
                                        # family = "Times",
                                        # size = 10,
                                        vjust = 0.5,
                                        colour = "black"),
              # axis.ticks.x = element_blank(),
              # axis.ticks.y = element_blank(), # to control the appearance
              panel.spacing = unit(1.5, "in"), # space around the panel in grid
              legend.margin = margin(2, 2, 2, 2), # space around the legend box
              legend.position = "bottom",     # of the plot
              # legend.key.width = unit(3.9, "cm"), # length of the legend
              legend.box = "vertical",
              legend.spacing = unit(0.2, "lines"),
              legend.box.margin = margin(-8, 0, 0, 0), #top, right, bottom, left
              legend.background = element_rect(fill = "white",
                                               size = 0.5,
                                               linetype = "solid",
                                               colour = "black"),
              legend.text = element_text(#family = "Times",
                #size = 10,
                colour = "black"),
              legend.title = element_text(#family = "Times",
                #size = 10,
                colour = "black"),
              strip.text = element_text(#family = "Times",
                #size = 10,
                colour = "black"),
              plot.margin = unit(c(5, 0.5, 0.5, 0.5), "in")) + #top, right, bottom, left
            scale_y_continuous(breaks = round(seq.int(from = ceiling(a_aux_1_serie %>%
                                                                       filter(sex == i,
                                                                              tran == j) %>%
                                                                       dplyr::select(c("Loglik",
                                                                                       "sex", "tran")) %>%
                                                                       dplyr::select(c("Loglik")) %>%
                                                                       filter(Loglik == min( Loglik[Loglik!= min(Loglik)]),
                                                                              !duplicated(Loglik)) %>%
                                                                       as.numeric() - 1),
                                                      to = ceiling(a_aux_1_serie %>%
                                                                     filter(sex == i,
                                                                            tran == j) %>%
                                                                     dplyr::select(c("Loglik",
                                                                                     "sex", "tran")) %>%
                                                                     dplyr::select(c("Loglik")) %>%
                                                                     filter(Loglik == max(Loglik),
                                                                            !duplicated(Loglik)) %>%
                                                                     as.numeric()),
                                                      length.out = 2 + 3), 1),
                               limits = c(ceiling(a_aux_1_serie %>%
                                                    filter(sex == i,
                                                           tran == j) %>%
                                                    dplyr::select(c("Loglik",
                                                                    "sex", "tran")) %>%
                                                    dplyr::select(c("Loglik")) %>%
                                                    filter(Loglik == min( Loglik[Loglik!= min(Loglik)]),
                                                           !duplicated(Loglik)) %>%
                                                    as.numeric() - 1), ceiling(a_aux_1_serie %>%
                                                       filter(sex == i,
                                                              tran == j) %>%
                                                       dplyr::select(c("Loglik",
                                                                       "sex", "tran")) %>%
                                                       dplyr::select(c("Loglik")) %>%
                                                       filter(Loglik == max(Loglik),
                                                              !duplicated(Loglik)) %>%
                                                       as.numeric()))) +
            scale_x_continuous(breaks = round(seq.int(from = 0,
                                                      to = ceiling(a_aux_1_serie %>%
                                                                     filter(sex == i,
                                                                            tran == j) %>%
                                                                     dplyr::select(c("Time",
                                                                                     "sex", "tran")) %>%
                                                                     dplyr::select(c("Time")) %>%
                                                                     filter(Time == max(Time),
                                                                            !duplicated(Time)) %>%
                                                                     as.numeric()),
                                                      length.out = 2 + 3), 0),
                               limits = c(0, ceiling(a_aux_1_serie %>%
                                                       filter(sex == i,
                                                              tran == j) %>%
                                                       dplyr::select(c("Time",
                                                                       "sex", "tran")) %>%
                                                       dplyr::select(c("Time")) %>%
                                                       filter(Time == max(Time),
                                                              !duplicated(Time)) %>%
                                                       as.numeric()))) +
            # guides(colour = guide_legend(override.aes = list(size = 2.5)),
            #        shape = guide_legend(override.aes = list(size = 2.5))) +
            labs(x = "Time (s)") +
            # scale_colour_manual(values = c("purple", "grey")) +
            theme_bw())
  }
}

tikz(file = "trace.tex",
     width = 6.5,
     height = 4.5,
     standAlone = FALSE,
     onefile = FALSE,
     engine = "pdftex",
     sanitize = TRUE)

ggarrange(plot_Men_Decrease, plot_Men_Remain, plot_Men_Increase,
          plot_Women_Decrease, plot_Women_Remain, plot_Women_Increase,
          common.legend = TRUE,
          legend = "bottom",
          ncol = 3,
          nrow = 2,
          align = "hv",
          hjust = 0,
          vjust = 0,
          font.label = list(size = 10),
          labels = c("a) Men, to decrease", "b) Men, to remain", "c) Men, to increase",
                     "d) Women, to decrease", "e) Women, to remain", "f) Women, to increase")) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0.5, "cm")) #top, right, bottom, left

dev.off()

# get the last iteration for each experiment
a_aux_1 = res %>%
  left_join(., beta_margin_params %>%
              rename(sex = Sex, tran = Transition),
            by = c("sex", "tran")) %>%
  group_by(id) %>%
  filter(Iteration == max(Iteration)) %>%
  ungroup()

# get the parameters at that iterations
em_params = em_params_aux %>%
  filter(id %in% (a_aux_1 %>%
                    dplyr::select(c("id", "area", "tran", "sex")) %>%
                    group_by(sex, tran) %>%
                    filter(area == max(area),
                           !duplicated(area)) %>%
                    ungroup() %>%
                    dplyr::select(c("id")) %>%
                    pull()))

handle = function(ns, sex, tran) {
  if(ns == "ALL") {
    em_params %>%
      filter(Sex == sex,
             Transition == tran)
  } else {
    em_params %>%
      filter(NS == ns,
             Sex == sex,
             Transition == tran)
  }
}

dens_mix = function(sex, tran) {
  ggplot() + theme_bw() +
    geom_function(fun = dbeta,
                  args = list(shape1 = beta_margin_params %>%
                                filter(Sex == sex, Transition == tran) %>%
                                dplyr::select(c("shape1")) %>% as.numeric(),
                              shape2 = beta_margin_params %>%
                                filter(Sex == sex, Transition == tran) %>%
                                dplyr::select(c("shape2")) %>% as.numeric()),
                  aes(colour = "Marginal"),
                  linetype = "dashed",
                  show.legend = FALSE,
                  col = "black") +
    # geom_density(data = tibble(value = LitBounds %>%
    #                              replace_na(list(Decrease = 0,
    #                                              Increase = 0)) %>%
    #                              filter(Sex == sex) %>%
    #                              dplyr::select(c(tran)) %>%
    #                              pull() %>%
    #                              as.numeric()),
    #              aes(x = value,
    #                  colour = "Literature"),
    #              kernel = "epanechnikov",
    #              col = "black",
    #              show.legend = FALSE) +
    geom_function(data = handle("Underweight", sex, tran),
                  fun = function(x, shape1, shape2, weight) {
                    dbeta(x, shape1, shape2)*weight },
                  args = list(weight = handle("Underweight", sex, tran)$weight,
                              shape1 = handle("Underweight", sex, tran)$alpha,
                              shape2 = handle("Underweight", sex, tran)$beta),
                  aes(colour = "Underweight")) +
    geom_function(data = handle("Normal", sex, tran),
                  fun = function(x, shape1, shape2, weight) {
                    dbeta(x, shape1, shape2)*weight },
                  args = list(weight = handle("Normal", sex, tran)$weight,
                              shape1 = handle("Normal", sex, tran)$alpha,
                              shape2 = handle("Normal", sex, tran)$beta),
                  aes(colour = "Normal")) +
    geom_function(data = handle("Overweight", sex, tran),
                  fun = function(x, shape1, shape2, weight) {
                    dbeta(x, shape1, shape2)*weight },
                  args = list(weight = handle("Overweight", sex, tran)$weight,
                              shape1 = handle("Overweight", sex, tran)$alpha,
                              shape2 = handle("Overweight", sex, tran)$beta),
                  aes(colour = "Overweight")) +
    geom_function(data = handle("Obese", sex, tran),
                  fun = function(x, shape1, shape2, weight) {
                    dbeta(x, shape1, shape2)*weight },
                  args = list(weight = handle("Obese", sex, tran)$weight,
                              shape1 = handle("Obese", sex, tran)$alpha,
                              shape2 = handle("Obese", sex, tran)$beta),
                  aes(colour = "Obese")) +
    # geom_function(data = handle("NA", sex, tran),
    #               fun = function(x, shape1, shape2, weight) {
    #                 dbeta(x, shape1, shape2)*weight },
    #               args = list(weight = handle("NA", sex, tran)$weight,
    #                           shape1 = handle("NA", sex, tran)$alpha,
    #                           shape2 = handle("NA", sex, tran)$beta),
    #               aes(colour = "Non-observed")) +
    theme(plot.title = element_text(#family = "Times",
                                    #size = 10,
                                    colour = "black"),
          axis.text.x = element_text(angle = 0, # whether to show the lines on
                                     hjust = 0.5, # the axis
                                     # family = "Times",
                                     # size = 10,
                                     vjust = 0.5,
                                     colour = "black"),
          axis.text.y = element_text(angle = 0, # whether to show the lines on
                                     hjust = 0.5, # the axis
                                     # family = "Times",
                                     # size = 10,
                                     vjust = 0.5,
                                     colour = "black"),
          axis.title = element_text(angle = 0, # whether to show the lines on
                                    hjust = 0.5, # the axis
                                    # family = "Times",
                                    # size = 10,
                                    vjust = 0.5,
                                    colour = "black"),
          # axis.ticks.x = element_blank(),
          # axis.ticks.y = element_blank(), # to control the appearance
          panel.spacing = unit(1, "lines"), # space around the panel in grid
          legend.margin = margin(2, 2, 2, 2), # space around the legend box
          legend.position = "bottom",     # of the plot
          # legend.key.width = unit(3.9, "cm"), # length of the legend
          legend.box = "vertical",
          legend.spacing = unit(0.2, "lines"),
          legend.box.margin = margin(-8, 0, 0, 0), #top, right, bottom, left
          # legend.background = element_rect(fill = "white",
          #                                  size = 0.5,
          #                                  linetype = "solid",
          #                                  colour = "black"),
          legend.text = element_text(#family = "Times",
                                     #size = 10,
                                     colour = "black"),
          legend.title = element_text(#family = "Times",
                                      #size = 10,
                                      colour = "black"),
          strip.text = element_text(#family = "Times",
                                    #size = 10,
                                    colour = "black"),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + #top, right, bottom, left
    guides(col = guide_legend(override.aes = list(size = 2))) + #increase legend figures
    labs(x = "Probability",
         y = "density",
         colour = "Nutritional status") +
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                       limits = c(0, 1)) +
    ylim(c(0, 2.5)) +
    scale_colour_discrete(breaks = c("Marginal", "Literature",
                                     "Underweight", "Normal",
                                     "Overweight", "Obese"))#,
                                     # "Non-observed"))
}

for (i in c("Men", "Women")) {
  for (j in c("Decrease", "Remain", "Increase")) {
    assign(paste0("plot_mix_", i, "_", j),
           dens_mix(i, j) +
             geom_function(data = handle("ALL", i, j),
                           fun = function(x, shape1, shape2, weight) {
                             dbeta(x, shape1[4], shape2[4])*weight[4] + # underweight
                               dbeta(x, shape1[1], shape2[1])*weight[1] + # normal
                               dbeta(x, shape1[3], shape2[3])*weight[3] + # overweight
                               dbeta(x, shape1[2], shape2[2])*weight[2]}, # obese
                           args = list(weight = handle("ALL", i, j)$weight,
                                       shape1 = handle("ALL", i, j)$alpha,
                                       shape2 = handle("ALL", i, j)$beta),
                           # aes(colour = "Mixture"),
                           linetype = "dotted",
                           show.legend = FALSE) +
             geom_area(data = NULL, aes(c(0, 1)),
                       stat = "function",
                       fun = function(x, shape1, shape2, weight) {
                         pmin(dbeta(x, shape1[5], shape2 = shape2[5]),
                              dbeta(x, shape1[4], shape2[4])*weight[4] + # underweight
                                dbeta(x, shape1[1], shape2[1])*weight[1] + # normal
                                dbeta(x, shape1[3], shape2[3])*weight[3] + # overweight
                                dbeta(x, shape1[2], shape2[2])*weight[2])}, # obese
                       fill = "grey",
                       xlim = c(0, 1),
                    alpha = 0.5,
                    args = list(weight = add_row(handle("ALL", i, j) %>%
                                                   dplyr::select(-c("NS", "id")),
                                                 beta_margin_params %>%
                                                   filter(Sex == i,
                                                          Transition == j) %>%
                                                   rename(alpha = shape1,
                                                          beta = shape2) %>%
                                                   mutate(weight = 1))$weight,
                                shape1 = add_row(handle("ALL", i, j) %>%
                                                   dplyr::select(-c("NS", "id")),
                                                 beta_margin_params %>%
                                                   filter(Sex == i,
                                                          Transition == j) %>%
                                                   rename(alpha = shape1,
                                                          beta = shape2) %>%
                                                   mutate(weight = 1))$alpha,
                                shape2 = add_row(handle("ALL", i, j) %>%
                                                   dplyr::select(-c("NS", "id")),
                                                 beta_margin_params %>%
                                                   filter(Sex == i, Transition == j) %>%
                                                   rename(alpha = shape1,
                                                          beta = shape2) %>%
                                                   mutate(weight = 1))$beta),
                    # aes(colour = "Mixture"),
                    linetype = "longdash",
                    show.legend = FALSE))
  }
}

tikz(file = "mixplot.tex",
     width = 6.5,
     height = 4.5,
     standAlone = FALSE,
     onefile = FALSE,
     engine = "pdftex",
     sanitize = TRUE)

ggarrange(plot_mix_Men_Decrease, plot_mix_Men_Remain, plot_mix_Men_Increase,
          plot_mix_Women_Decrease, plot_mix_Women_Remain, plot_mix_Women_Increase,
          common.legend = TRUE,
          legend = "bottom",
          align = "hv",
          hjust = 0,
          vjust = 1.1,
          font.label = list(size = 10),
          labels = c("a) Men, to decrease", "b) Men, to remain", "c) Men, to increase",
                     "d) Women, to decrease", "e) Women, to remain", "f) Women, to increase"))

dev.off()

ggsave("mix_density.png",
       width = 6.5, height = 4.5, units = "in")

# print table in latex format
print(xtable(em_params %>%
               dplyr::select(-c("id")) %>%
               add_row(beta_margin_params %>%
                         mutate(NS = "Marginal") %>%
                         rename(alpha = shape1, beta = shape2)) %>%
               mutate(NS = factor(NS, levels = c("Underweight", "Normal",
                                                 "Overweight", "Obese", #"NA",
                                                 "Marginal"),
                                  labels = c("Underweight", "Normal",
                                             "Overweight", "Obese", #"Non-observed",
                                             "Marginal"))) %>%
               gather(key = "Parameter", value = "value",
                      -Sex, -Transition, -NS) %>%
               mutate(Parameter = factor(Parameter,
                                         levels = c("alpha", "beta", "weight"))) %>%
               unite(temp, NS, Parameter) %>%
               filter(!is.na(value)) %>%
               mutate(temp = factor(temp, levels = c(paste("Underweight",
                                                           c("alpha", "beta", "weight"),
                                                           sep = "_"),
                                                     paste("Normal",
                                                           c("alpha", "beta", "weight"),
                                                           sep = "_"),
                                                     paste("Overweight",
                                                           c("alpha", "beta", "weight"),
                                                           sep = "_"),
                                                     paste("Obese",
                                                           c("alpha", "beta", "weight"),
                                                           sep = "_"),
                                                     paste("Marginal",
                                                           c("alpha", "beta"),
                                                           sep = "_"))),
                      Transition = factor(Transition, levels = c("Decrease",
                                                                 "Remain",
                                                                 "Increase"),
                                          labels = c("To decrease",
                                                     "To remain",
                                                     "To increase"))) %>%
               # select(-c("weight")) %>%
               spread(temp, value),
             digits = 2,
             label = "tab:mixture",
             caption = "Parameters of the mixture density for combination of nutritional status, sex and kind of transition."),
      caption.placement = "top",
      comment = FALSE,
      include.rownames = FALSE,
      include.colnames = TRUE,
      type = "latex", sanitize.text.function = function(x){x},
      hline.after = c(-1, 0, 5))

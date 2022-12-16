beta_margin_params = readRDS("beta_margin.rds")

beta_margin_plot = beta_margin_params %>%
  filter(Sex == "Men",
         Transition == "Remain")

grid_2 = readRDS("grid_2.rds")

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

plot_beta = em_params_aux %>%
  filter(Sex == "Men",
         Transition == "Remain")

# 325 - 648
# plot_id = c(330, 342, 355)
plot_id = c(330)

plot_beta %>%
  filter(id %in% plot_id)

for (i in 1:length(plot_id)) {
  assign(paste0("plot_beta_", i),
         ggplot() +
           theme_bw() +
           # marginal
           geom_function(data = beta_margin_plot,
                         fun = dbeta,
                         args = list(shape1 = beta_margin_plot$shape1,
                                     shape2 = beta_margin_plot$shape2),
                         aes(colour = "Marginal"),
                         linetype = "dashed",
                         show.legend = FALSE,
                         col = "black") +
           # underweight
           geom_function(data = plot_beta,
                         fun = function(x, shape1, shape2, weight) {
                           dbeta(x, shape1, shape2)*weight },
                         args = list(weight = plot_beta %>%
                                       filter(id == plot_id[i],
                                              NS == "Underweight") %>%
                                       dplyr::select(c("weight")) %>%
                                       as.numeric(),
                                     shape1 = plot_beta %>%
                                       filter(id == plot_id[i],
                                              NS == "Underweight") %>%
                                       dplyr::select(c("alpha")) %>%
                                       as.numeric(),
                                     shape2 = plot_beta %>%
                                       filter(id == plot_id[i],
                                              NS == "Underweight") %>%
                                       dplyr::select(c("beta")) %>%
                                       as.numeric()),
                         aes(colour = "Underweight")) +
           # normal
           geom_function(data = plot_beta,
                         fun = function(x, shape1, shape2, weight) {
                           dbeta(x, shape1, shape2)*weight },
                         args = list(weight = plot_beta %>%
                                       filter(id == plot_id[i],
                                              NS == "Normal") %>%
                                       dplyr::select(c("weight")) %>%
                                       as.numeric(),
                                     shape1 = plot_beta %>%
                                       filter(id == plot_id[i],
                                              NS == "Normal") %>%
                                       dplyr::select(c("alpha")) %>%
                                       as.numeric(),
                                     shape2 = plot_beta %>%
                                       filter(id == plot_id[i],
                                              NS == "Normal") %>%
                                       dplyr::select(c("beta")) %>%
                                       as.numeric()),
                         aes(colour = "Normal")) +
           # overweight
           geom_function(data = plot_beta,
                         fun = function(x, shape1, shape2, weight) {
                           dbeta(x, shape1, shape2)*weight },
                         args = list(weight = plot_beta %>%
                                       filter(id == plot_id[i],
                                              NS == "Overweight") %>%
                                       dplyr::select(c("weight")) %>%
                                       as.numeric(),
                                     shape1 = plot_beta %>%
                                       filter(id == plot_id[i],
                                              NS == "Overweight") %>%
                                       dplyr::select(c("alpha")) %>%
                                       as.numeric(),
                                     shape2 = plot_beta %>%
                                       filter(id == plot_id[i],
                                              NS == "Overweight") %>%
                                       dplyr::select(c("beta")) %>%
                                       as.numeric()),
                         aes(colour = "Overweight")) +
           # obese
           geom_function(data = plot_beta,
                         fun = function(x, shape1, shape2, weight) {
                           dbeta(x, shape1, shape2)*weight },
                         args = list(weight = plot_beta %>%
                                       filter(id == plot_id[i],
                                              NS == "Obese") %>%
                                       dplyr::select(c("weight")) %>%
                                       as.numeric(),
                                     shape1 = plot_beta %>%
                                       filter(id == plot_id[i],
                                              NS == "Obese") %>%
                                       dplyr::select(c("alpha")) %>%
                                       as.numeric(),
                                     shape2 = plot_beta %>%
                                       filter(id == plot_id[i],
                                              NS == "Obese") %>%
                                       dplyr::select(c("beta")) %>%
                                       as.numeric()),
                         aes(colour = "Obese")) +
           # mixture
           geom_function(data = plot_beta,
                         fun = function(x, shape1, shape2, weight) {
                           dbeta(x, shape1[4], shape2[4])*weight[4] + # underweight
                             dbeta(x, shape1[1], shape2[1])*weight[1] + # normal
                             dbeta(x, shape1[3], shape2[3])*weight[3] + # overweight
                             dbeta(x, shape1[2], shape2[2])*weight[2]}, # obese
                         args = list(weight = plot_beta %>%
                                       filter(id == plot_id[i]) %>%
                                       dplyr::select(c("weight")) %>%
                                       pull() %>%
                                       as.numeric(),
                                     shape1 = plot_beta %>%
                                       filter(id == plot_id[i]) %>%
                                       dplyr::select(c("alpha")) %>%
                                       pull() %>%
                                       as.numeric(),
                                     shape2 = plot_beta %>%
                                       filter(id == plot_id[i]) %>%
                                       dplyr::select(c("beta")) %>%
                                       pull() %>%
                                       as.numeric()),
                         aes(colour = "Mixture"),
                         linetype = "dotted",
                         show.legend = FALSE) +
           # area
           geom_area(data = NULL,
                     aes(c(0, 1)),
                     stat = "function",
                     fun = function(x, shape1, shape2, weight) {
                       pmin(dbeta(x, shape1[5], shape2 = shape2[5]), # marginal
                            dbeta(x, shape1[4], shape2[4])*weight[4] + # underweight
                              dbeta(x, shape1[1], shape2[1])*weight[1] + # normal
                              dbeta(x, shape1[3], shape2[3])*weight[3] + # overweight
                              dbeta(x, shape1[2], shape2[2])*weight[2]) }, # obese
                     fill = "grey",
                     xlim = c(0, 1),
                     alpha = 0.5,
                     args = list(weight = plot_beta %>%
                                   filter(id == plot_id[i]) %>%
                                   dplyr::select(c("weight")) %>%
                                   pull() %>%
                                   as.numeric(),
                                 shape1 = c(plot_beta %>%
                                              filter(id == plot_id[i]) %>%
                                              dplyr::select(c("alpha")) %>%
                                              pull() %>%
                                              as.numeric(),
                                            beta_margin_plot$shape1),
                                 shape2 = c(plot_beta %>%
                                              filter(id == plot_id[i]) %>%
                                              dplyr::select(c("beta")) %>%
                                              pull() %>%
                                              as.numeric(),
                                            beta_margin_plot$shape2)),
                     # aes(colour = "Mixture"),
                     linetype = "longdash",
                     show.legend = FALSE) +
           # add the vectical lines indicating the intervals R
           geom_vline(xintercept = c(0, 1),
                      linetype = "dotdash",
                      col = "brown",
                      alpha = 0.25,
                      size = 1.5) +
           geom_vline(xintercept = optim(par = 0.5,
                                         function(x, shape1, shape2, weight) {
                                           shape1 = c(plot_beta %>%
                                                        filter(id == plot_id[i]) %>%
                                                        dplyr::select(c("alpha")) %>%
                                                        pull() %>%
                                                        as.numeric(),
                                                      beta_margin_plot$shape1)
                                           shape2 = c(plot_beta %>%
                                                        filter(id == plot_id[i]) %>%
                                                        dplyr::select(c("beta")) %>%
                                                        pull() %>%
                                                        as.numeric(),
                                                      beta_margin_plot$shape2)
                                           weight = plot_beta %>%
                                             filter(id == plot_id[i]) %>%
                                             dplyr::select(c("weight")) %>%
                                             pull() %>%
                                             as.numeric()
                                           -pmin(dbeta(x, shape1[5], shape2 = shape2[5]), # marginal
                                                dbeta(x, shape1[4], shape2[4])*weight[4] + # underweight
                                                  dbeta(x, shape1[1], shape2[1])*weight[1] + # normal
                                                  dbeta(x, shape1[3], shape2[3])*weight[3] + # overweight
                                                  dbeta(x, shape1[2], shape2[2])*weight[2]) }, # obese,
                                         method = "BFGS")$par,
                      linetype = "dotdash",
                      col = "brown",
                      alpha = 0.25,
                      size = 1.5) +
           geom_text(aes(label = "M1",
                         x = optim(par = 0.5,
                                   function(x, shape1, shape2, weight) {
                                     shape1 = c(plot_beta %>%
                                                  filter(id == plot_id[i]) %>%
                                                  dplyr::select(c("alpha")) %>%
                                                  pull() %>%
                                                  as.numeric(),
                                                beta_margin_plot$shape1)
                                     shape2 = c(plot_beta %>%
                                                  filter(id == plot_id[i]) %>%
                                                  dplyr::select(c("beta")) %>%
                                                  pull() %>%
                                                  as.numeric(),
                                                beta_margin_plot$shape2)
                                     weight = plot_beta %>%
                                       filter(id == plot_id[i]) %>%
                                       dplyr::select(c("weight")) %>%
                                       pull() %>%
                                       as.numeric()
                                     -pmin(dbeta(x, shape1[5], shape2 = shape2[5]), # marginal
                                           dbeta(x, shape1[4], shape2[4])*weight[4] + # underweight
                                             dbeta(x, shape1[1], shape2[1])*weight[1] + # normal
                                             dbeta(x, shape1[3], shape2[3])*weight[3] + # overweight
                                             dbeta(x, shape1[2], shape2[2])*weight[2]) }, # obese,
                                   method = "BFGS")$par/2,
                         y = 2.5),
                     parse = TRUE,
                     col = "black") +
           geom_text(aes(label = "M2",
                         x = (1+optim(par = 0.5,
                                      function(x, shape1, shape2, weight) {
                                        shape1 = c(plot_beta %>%
                                                     filter(id == plot_id[i]) %>%
                                                     dplyr::select(c("alpha")) %>%
                                                     pull() %>%
                                                     as.numeric(),
                                                   beta_margin_plot$shape1)
                                        shape2 = c(plot_beta %>%
                                                     filter(id == plot_id[i]) %>%
                                                     dplyr::select(c("beta")) %>%
                                                     pull() %>%
                                                     as.numeric(),
                                                   beta_margin_plot$shape2)
                                        weight = plot_beta %>%
                                          filter(id == plot_id[i]) %>%
                                          dplyr::select(c("weight")) %>%
                                          pull() %>%
                                          as.numeric()
                                        -pmin(dbeta(x, shape1[5], shape2 = shape2[5]), # marginal
                                              dbeta(x, shape1[4], shape2[4])*weight[4] + # underweight
                                                dbeta(x, shape1[1], shape2[1])*weight[1] + # normal
                                                dbeta(x, shape1[3], shape2[3])*weight[3] + # overweight
                                                dbeta(x, shape1[2], shape2[2])*weight[2]) }, # obese,
                                      method = "BFGS")$par)/2,
                         y = 2.5),
                     parse = TRUE,
                     col = "black") +
           # raices
           geom_label(aes(label = c("R1jk", "R2jk", "R3jk"),
                          x = c(0, optim(par = 0.5,
                                         function(x, shape1, shape2, weight) {
                                           shape1 = c(plot_beta %>%
                                                        filter(id == plot_id[i]) %>%
                                                        dplyr::select(c("alpha")) %>%
                                                        pull() %>%
                                                        as.numeric(),
                                                      beta_margin_plot$shape1)
                                           shape2 = c(plot_beta %>%
                                                        filter(id == plot_id[i]) %>%
                                                        dplyr::select(c("beta")) %>%
                                                        pull() %>%
                                                        as.numeric(),
                                                      beta_margin_plot$shape2)
                                           weight = plot_beta %>%
                                             filter(id == plot_id[i]) %>%
                                             dplyr::select(c("weight")) %>%
                                             pull() %>%
                                             as.numeric()
                                           -pmin(dbeta(x, shape1[5], shape2 = shape2[5]), # marginal
                                                 dbeta(x, shape1[4], shape2[4])*weight[4] + # underweight
                                                   dbeta(x, shape1[1], shape2[1])*weight[1] + # normal
                                                   dbeta(x, shape1[3], shape2[3])*weight[3] + # overweight
                                                   dbeta(x, shape1[2], shape2[2])*weight[2]) }, # obese,
                                         method = "BFGS")$par, 1),
                          y = c(2, 2, 2),
                          angle = 90,
                          label.size = NULL),
                      col = "black") +
           theme(plot.title = element_text(colour = "black"),
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
           guides(col = guide_legend(override.aes = list(size = 2),
                                     nrow = 2, byrow = TRUE)) + #increase legend figures
           labs(x = "Transition probability",
                y = "density",
                colour = "") +
           scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2),
                              limits = c(-0.15, 1.15)) +
           ylim(c(0, 2.5)) +
           scale_colour_discrete(breaks = c("Underweight",
                                            "Normal",
                                            "Overweight",
                                            "Obese")))
}

tikz(file = "Document/Figures/BetaIllustration.tex",
     width = 3.5,
     height = 3,
     standAlone = FALSE,
     onefile = FALSE,
     engine = "pdftex",
     sanitize = TRUE)

ggarrange(plot_beta_1,
          #plot_beta_2,
          #plot_beta_3,
          ncol = 1,
          common.legend = FALSE,
          legend = "none",
          align = "hv",
          hjust = -1,
          vjust = 1,
          font.label = list(size = 10))

dev.off()

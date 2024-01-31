
library(ggalt)
library(ggtext)
library(ggrepel)
library(ggpubr)
#  trace plot and density:
# log-likelihood
trace_sp_lkhd <-  ggplot(data = spprDraws) + 
  geom_line(aes(x = index, y = lg_lkhd)) +
  labs(title = expression(paste("Trace of Log-likelihood from point-level spatial model")),
       x = "Iteration Number",
       y = "Log likelihood") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))

den_sp_lkhd <- ggplot(data = spprDraws) +
  geom_histogram(aes(x = lg_lkhd, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  # geom_density(aes( x = rnorm(5000, 20, 10000)), col = "seagreen", lwd = 1) +
  labs(title = expression(paste("Density of Log-likelihood from point-level spatial model")),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# log likelihood from MVR
trace_mv_lkhd <-  ggplot(data = mvrDraws) + 
  geom_line(aes(x = index, y = lg_lkhd)) +
  labs(title = expression(paste("Trace of Log-likelihood from non-spatial regression")),
       x = "Iteration Number",
       y = "Log likelihood") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))

den_mv_lkhd <- ggplot(data = mvrDraws) +
  geom_histogram(aes(x = lg_lkhd, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") +
  # geom_density(aes( x = lg_lkhd), col = "seagreen", linewidth = 1) +
  labs(title = expression(paste("Density of Log-likelihood from non-spatial regression")),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
    
# theta bar
trace_sp_thbar <-  ggplot(data = spprDraws) + 
  geom_line(aes(x = index, y = th_mean)) +
  labs(title = expression(paste("Trace of ", bar(theta))),
       x = "Iteration Number",
       y = expression(bar(theta))) +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))
  # scale_color_manual(values = col_pal)

den_sp_thbar <- ggplot(data = spprDraws) +
  geom_histogram(aes(x = th_mean, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  # geom_density(aes( x = rnorm(5000, 20, 10000)), col = "seagreen", lwd = 1) +
  labs(title = expression(paste("Density of ", bar(theta))),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# beta_0
trace_sp_b0 <-  ggplot(data = spprDraws) + 
  geom_line(aes(x = index, y = beta0)) +
  labs(title = expression(paste("Trace of ", beta[0])),
       x = "Iteration Number",
       y = expression(beta[0])) +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))
# scale_color_manual(values = col_pal)

den_sp_b0 <- ggplot(data = spprDraws) +
  geom_histogram(aes(x = beta0, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  # geom_density(aes( x = rnorm(5000, 20, 10000)), col = "seagreen", lwd = 1) +
  labs(title = expression(paste("Density of ", beta[0])),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# beta_1
trace_sp_b1 <-  ggplot(data = spprDraws) + 
  geom_line(aes(x = index, y = beta1)) +
  labs(title = expression(paste("Trace of ", beta[1])),
       x = "Iteration Number",
       y = expression(beta[1])) +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))
# scale_color_manual(values = col_pal)

den_sp_b1 <- ggplot(data = spprDraws) +
  geom_histogram(aes(x = beta1, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  # geom_density(aes( x = rnorm(5000, 20, 10000)), col = "seagreen", lwd = 1) +
  labs(title = expression(paste("Density of ", beta[1])),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
# beta_2
trace_sp_b2 <-  ggplot(data = spprDraws) + 
  geom_line(aes(x = index, y = beta2)) +
  labs(title = expression(paste("Trace of ", beta[2])),
       x = "Iteration Number",
       y = expression(beta[2])) +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))
# scale_color_manual(values = col_pal)

den_sp_b2 <- ggplot(data = spprDraws) +
  geom_histogram(aes(x = beta2, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  # geom_density(aes( x = rnorm(5000, 20, 10000)), col = "seagreen", lwd = 1) +
  labs(title = expression(paste("Density of ", beta[2])),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
# sigma2
trace_sp_s2 <-  ggplot(data = spprDraws) + 
  geom_line(aes(x = index, y = sigma2)) +
  labs(title = expression(paste("Trace of ", sigma^2)),
       x = "Iteration Number",
       y = expression(sigma^2)) +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))
# scale_color_manual(values = col_pal)

den_sp_s2 <- ggplot(data = spprDraws) +
  geom_histogram(aes(x = sigma2, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  # geom_density(aes( x = rnorm(5000, 20, 10000)), col = "seagreen", lwd = 1) +
  labs(title = expression(paste("Density of ", sigma^2)),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# tau_2
trace_sp_tu2 <-  ggplot(data = spprDraws) + 
  geom_line(aes(x = index, y = tau2)) +
  labs(title = expression(paste("Trace of ", tau^2)),
       x = "Iteration Number",
       y = expression(tau^2)) +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))
# scale_color_manual(values = col_pal)

den_sp_tu2 <- ggplot(data = spprDraws) +
  geom_histogram(aes(x = tau2, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  # geom_density(aes( x = rnorm(5000, 20, 10000)), col = "seagreen", lwd = 1) +
  labs(title = expression(paste("Density of ", tau^2)),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# lambda
trace_sp_lambda <-  ggplot(data = spprDraws) + 
  geom_line(aes(x = index, y = lambda)) +
  labs(title = expression(paste("Trace of ", lambda)),
       x = "Iteration Number",
       y = expression(lambda)) +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))
# scale_color_manual(values = col_pal)

den_sp_lambda <- ggplot(data = spprDraws) +
  geom_histogram(aes(x = lambda, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  # geom_density(aes( x = rnorm(5000, 20, 10000)), col = "seagreen", lwd = 1) +
  labs(title = expression(paste("Density of ", lambda)),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# for non-spatial

trace_mv_b0 <-  ggplot(data = mvrDraws) + 
  geom_line(aes(x = index, y = beta0)) +
  labs(title = expression(paste("Trace of Non-spatial ", beta[0])),
       x = "Iteration Number",
       y = expression(beta[0])) +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))
# scale_color_manual(values = col_pal)

den_mv_b0 <- ggplot(data = mvrDraws) +
  geom_histogram(aes(x = beta0, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  # geom_density(aes( x = rnorm(5000, 20, 10000)), col = "seagreen", lwd = 1) +
  labs(title = expression(paste("Density of ", beta[0])),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# beta1
trace_mv_b1 <-  ggplot(data = mvrDraws) + 
  geom_line(aes(x = index, y = beta1)) +
  labs(title = expression(paste("Trace of Non-spatial ", beta[1])),
       x = "Iteration Number",
       y = expression(beta[1])) +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))
# scale_color_manual(values = col_pal)

den_mv_b1 <- ggplot(data = mvrDraws) +
  geom_histogram(aes(x = beta1, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  # geom_density(aes( x = rnorm(5000, 20, 10000)), col = "seagreen", lwd = 1) +
  labs(title = expression(paste("Density of ", beta[1])),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
# beta2
trace_mv_b2 <-  ggplot(data = mvrDraws) + 
  geom_line(aes(x = index, y = beta2)) +
  labs(title = expression(paste("Trace of Non-spatial ", beta[2])),
       x = "Iteration Number",
       y = expression(beta[2])) +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))
# scale_color_manual(values = col_pal)

den_mv_b2 <- ggplot(data = mvrDraws) +
  geom_histogram(aes(x = beta2, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  # geom_density(aes( x = rnorm(5000, 20, 10000)), col = "seagreen", lwd = 1) +
  labs(title = expression(paste("Density of ", beta[2])),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
# sigma2
trace_mv_s2 <-  ggplot(data = mvrDraws) + 
  geom_line(aes(x = index, y = sigma2)) +
  labs(title = expression(paste("Trace of Non-spatial ", sigma^2)),
       x = "Iteration Number",
       y = expression(sigma^2)) +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, size = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))
# scale_color_manual(values = col_pal)

den_mv_s2 <- ggplot(data = mvrDraws) +
  geom_histogram(aes(x = sigma2, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  # geom_density(aes( x = rnorm(5000, 20, 10000)), col = "seagreen", lwd = 1) +
  labs(title = expression(paste("Density of ", sigma^2)),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# prediction plot
legendColor <- c("True response" = "red3", "Spatial regression" = "lightseagreen", "Non-spatial regression" = "indianred1")
pred_sp_y1 <- ggplot() +
  # geom_histogram(aes(x = lambda, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  geom_density(data = sp_predictiondf, aes(x = y_1, col = "Spatial regression"), lwd = 1) + 
  geom_density(data = mv_predictiondf, aes(x = y_1,col = "Non-spatial regression"), lwd = 1) +
  geom_vline(xintercept = c(true_response[1], sp_pred_medians["y_1"], mv_pred_medians["y_1"], sp_pred_quantiles[, "y_1"], mv_pred_quantiles[, "y_1"]),
             col = c("red3", "lightseagreen", "indianred1", "lightseagreen", "lightseagreen", "indianred1", "indianred1"),
             linetype = c("solid", "F1", "dashed", "F1", "F1", "dashed", "dashed"), linewidth = 1) +

  labs(title = expression(paste("Prediction of ", y[1])),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  scale_color_manual(values = legendColor)
# print(pred_sp_y1)

pred_sp_y2 <- ggplot() +
  # geom_histogram(aes(x = lambda, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  geom_density(data = sp_predictiondf, aes(x = y_2), col = "lightseagreen", lwd = 1) + 
  geom_density(data = mv_predictiondf, aes(x = y_2), col = "indianred1", lwd = 1) +
  geom_vline(xintercept = c(true_response[2], sp_pred_medians["y_2"], mv_pred_medians["y_2"], sp_pred_quantiles[, "y_2"], mv_pred_quantiles[, "y_2"]),
             col = c("red4", "lightseagreen", "indianred1", "lightseagreen", "lightseagreen", "indianred1", "indianred1"),
             linetype = c("solid", "F1", "dashed", "F1", "F1", "dashed", "dashed"), linewidth = 1) +
  
  labs(title = expression(paste("Prediction of ", y[2])),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

pred_sp_y3 <- ggplot() +
  # geom_histogram(aes(x = lambda, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  geom_density(data = sp_predictiondf, aes(x = y_3), col = "lightseagreen", lwd = 1) + 
  geom_density(data = mv_predictiondf, aes(x = y_3), col = "indianred1", lwd = 1) +
  geom_vline(xintercept = c(true_response[3], sp_pred_medians["y_3"], mv_pred_medians["y_3"], sp_pred_quantiles[, "y_3"], mv_pred_quantiles[, "y_3"]),
             col = c("red4", "lightseagreen", "indianred1", "lightseagreen", "lightseagreen", "indianred1", "indianred1"),
             linetype = c("solid", "F1", "dashed", "F1", "F1", "dashed", "dashed"), linewidth = 1) +
  
  labs(title = expression(paste("Prediction of ", y[3])),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

pred_sp_y4 <- ggplot() +
  # geom_histogram(aes(x = lambda, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  geom_density(data = sp_predictiondf, aes(x = y_4), col = "lightseagreen", lwd = 1) + 
  geom_density(data = mv_predictiondf, aes(x = y_4), col = "indianred1", lwd = 1) +
  geom_vline(xintercept = c(true_response[4], sp_pred_medians["y_4"], mv_pred_medians["y_4"], sp_pred_quantiles[, "y_4"], mv_pred_quantiles[, "y_4"]),
             col = c("red4", "lightseagreen", "indianred1", "lightseagreen", "lightseagreen", "indianred1", "indianred1"),
             linetype = c("solid", "F1", "dashed", "F1", "F1", "dashed", "dashed"), linewidth = 1) +
  
  labs(title = expression(paste("Prediction of ", y[4])),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

pred_sp_y5 <- ggplot() +
  # geom_histogram(aes(x = lambda, y = after_stat(density)), bins = 30, fill = "white", col = "#003153") + 
  geom_density(data = sp_predictiondf, aes(x = y_5), col = "lightseagreen", lwd = 1) + 
  geom_density(data = mv_predictiondf, aes(x = y_5), col = "indianred1", lwd = 1) +
  geom_vline(xintercept = c(true_response[5], sp_pred_medians["y_5"], mv_pred_medians["y_5"], sp_pred_quantiles[, "y_5"], mv_pred_quantiles[, "y_5"]),
             col = c("red4", "lightseagreen", "indianred1", "lightseagreen", "lightseagreen", "indianred1", "indianred1"),
             linetype = c("solid", "F1", "dashed", "F1", "F1", "dashed", "dashed"), linewidth = 1) +
  
  labs(title = expression(paste("Prediction of ", y[5])),
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 12, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


prdtemp <- ggarrange(pred_sp_y1, pred_sp_y2, pred_sp_y3, pred_sp_y4, pred_sp_y5,  
                               common.legend = T, nrow = 2, ncol = 3)

predplot <- annotate_figure(prdtemp, top = text_grob(paste0("Distribtuion of prediction(with true response) using spatial and non-spatial regression (median and 90% credible interval)")
                        , face = "bold", size = 12))
# print(predplot)

trace_sp_plot1 <-  annotate_figure(ggarrange(trace_sp_b0, trace_sp_b1, trace_sp_b2, den_sp_b0, den_sp_b1, den_sp_b2,  
                        common.legend = T, nrow = 2, ncol = 3), top = text_grob(paste0("Trace plot of parameters from spatial point level regression")
                                                                                , face = "bold", size = 12)) 
# print(trace_sp_plot1)
trace_sp_plot2 <-  annotate_figure(ggarrange(trace_sp_thbar, trace_sp_tu2, trace_sp_s2, trace_sp_lambda, den_sp_thbar, den_sp_tu2, den_sp_s2, den_sp_lambda, 
                                             common.legend = T, nrow = 2, ncol = 4), 
                                   top = text_grob(paste0("Trace plot of parameters from spatial point level regression"), face = "bold", size = 12)) 

trace_mv_plot1 <-  annotate_figure(ggarrange(trace_mv_b0, trace_mv_b1, trace_mv_b2, trace_mv_s2, den_mv_b0, den_mv_b1, den_mv_b2, den_mv_s2,  
                                             common.legend = T, nrow = 2, ncol = 4), 
                                   top = text_grob(paste0("Trace plot of parameters from non-spatial regression"), face = "bold", size = 12)) 

lkhdplot <- annotate_figure(ggarrange(trace_sp_lkhd, trace_mv_lkhd,  
                      common.legend = T, nrow = 2, ncol = 1), top = text_grob(paste0("Log-likelihood plot spatial and non-spatial regression")
                                                                              , face = "bold", size = 12))

filename<-  paste0("summary output 1 acc 25_48.pdf")
pdf(filename, width = 14, height = 8.5, onefile = TRUE)
print(lkhdplot)
print(trace_sp_plot1)
print(trace_sp_plot2)
print(trace_mv_plot1)
print(predplot)
dev.off()
# graphics.off()
print(den_sp_lambda)

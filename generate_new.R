library(tidyverse)
library(brms)
library(tidybayes)
library(modelr)

# this script requires loading one of the brms model objects from
# photosynthesis_respiration_log_corr.R or
# photosynthesis_respiration_log_nocorr.R
# (which are too large for github, so you have to refit the model — sorry)

# The purpose here is to draw parameter samples from the joint posterior, which
# will represent ecologically plausible reaction norms for hypothetical species.
# They will be used in simulations of multi-species competitive communities whose
# dynamics play out under changing temperature conditions. The benefit of fitting
# a model with group-level correlations, is that if, for example, species with
# high lnc for photosynthesis also tend to have high lnc for respitation (they do),
# this correlation will constrain the samples we draw and will be represented in
# them




get_variables(m) # to see the names used by brms

# ignoring replicate variation:                                             ####
dr = m %>%
  spread_draws( b_lnGP_lnc_Intercept, r_species__lnGP_lnc[s,Intercept],
                b_lnGP_Ea_Intercept, r_species__lnGP_Ea[s,Intercept],
                b_lnGP_Eh_Intercept, r_species__lnGP_Eh[s,Intercept],
                b_lnGP_Th_Intercept, r_species__lnGP_Th[s,Intercept],

                b_lnR_lnc_Intercept, r_species__lnR_lnc[s,Intercept],
                b_lnR_Ea_Intercept, r_species__lnR_Ea[s,Intercept],
                b_lnR_Eh_Intercept, r_species__lnR_Eh[s,Intercept],
                b_lnR_Th_Intercept, r_species__lnR_Th[s,Intercept],
                ndraws = 100
                #,seed = 123
  ) %>%
  mutate(.keep = "unused",
         # we exponentiate to get back to the original scale
         gp_c =  exp(b_lnGP_lnc_Intercept + r_species__lnGP_lnc),
         gp_Ea = b_lnGP_Ea_Intercept + r_species__lnGP_Ea,
         gp_Eh = b_lnGP_Eh_Intercept + r_species__lnGP_Eh,
         gp_Th = b_lnGP_Th_Intercept + r_species__lnGP_Th,

         # we exponentiate to get back to the original scale
         r_c =   exp(b_lnR_lnc_Intercept + r_species__lnR_lnc),
         r_Ea =  b_lnR_Ea_Intercept + r_species__lnR_Ea,
         r_Eh =  b_lnR_Eh_Intercept + r_species__lnR_Eh,
         r_Th =  b_lnR_Th_Intercept + r_species__lnR_Th) %>%
  #filter(gp_c > 0 & r_c > 0) %>%
  filter(s != "sc") %>%
  ungroup() %>%
  slice_sample(n = 16)

# including replicate variation:                                            ####
dr = m %>%
  spread_draws(
    b_lnGP_lnc_Intercept,  r_species__lnGP_lnc[s,Intercept], r_rep__lnGP_lnc[r,Intercept],
    b_lnGP_Ea_Intercept, r_species__lnGP_Ea[s,Intercept], r_rep__lnGP_Ea[r,Intercept],
    b_lnGP_Eh_Intercept, r_species__lnGP_Eh[s,Intercept], r_rep__lnGP_Eh[r,Intercept],
    b_lnGP_Th_Intercept, r_species__lnGP_Th[s,Intercept], r_rep__lnGP_Th[r,Intercept],

    b_lnR_lnc_Intercept, r_species__lnR_lnc[s,Intercept], r_rep__lnR_lnc[r,Intercept],
    b_lnR_Ea_Intercept, r_species__lnR_Ea[s,Intercept], r_rep__lnR_Ea[r,Intercept],
    b_lnR_Eh_Intercept, r_species__lnR_Eh[s,Intercept], r_rep__lnR_Eh[r,Intercept],
    b_lnR_Th_Intercept, r_species__lnR_Th[s,Intercept], r_rep__lnR_Th[r,Intercept],
    ndraws = 100
    ,seed = 123
  ) %>%
  mutate(.keep = "unused",
         # we exponentiate to get back to the original scale
         gp_c =  exp(b_lnGP_lnc_Intercept + r_species__lnGP_lnc  + r_rep__lnGP_lnc),
         gp_Ea = b_lnGP_Ea_Intercept + r_species__lnGP_Ea + r_rep__lnGP_Ea,
         gp_Eh = b_lnGP_Eh_Intercept + r_species__lnGP_Eh + r_rep__lnGP_Eh,
         gp_Th = b_lnGP_Th_Intercept + r_species__lnGP_Th + r_rep__lnGP_Th,

         # we exponentiate to get back to the original scale
         r_c =   exp(b_lnR_lnc_Intercept + r_species__lnR_lnc  + r_rep__lnR_lnc) ,
         r_Ea =  b_lnR_Ea_Intercept + r_species__lnR_Ea + r_rep__lnR_Ea,
         r_Eh =  b_lnR_Eh_Intercept + r_species__lnR_Eh + r_rep__lnR_Eh,
         r_Th =  b_lnR_Th_Intercept + r_species__lnR_Th + r_rep__lnR_Th) %>%
  #filter(gp_c > 0 & r_c > 0) %>%
  filter(s != "sc") %>%
  ungroup() %>%
  slice_sample(n = 16)








temps <- tibble(K = seq(280.15, 320.15, length.out = 100))  # 0°C to 40°C

# Define Boltzmann constant
k = 8.62e-5

# Join temperature data with posterior draws
# draw_curves <- dr %>%
#   mutate(draw_id = row_number()) %>%
#   crossing(temps) %>%
#   mutate(
#     gp_rate = gp_c * exp(gp_Ea / k * (1 / 293.15 - 1 / K)) /
#       (1 + exp(gp_Eh / k * (1 / gp_Th - 1 / K))),
#
#     r_rate = r_c * exp(r_Ea / k * (1 / 293.15 - 1 / K)) /
#       (1 + exp(r_Eh / k * (1 / r_Th - 1 / K)))
#   )


draw_curves = dr %>%
  mutate(draw_id = row_number()) %>%
  crossing(temps) %>%
  mutate(
    gp_rate = gp_c * exp(gp_Ea / k * (1 / 293.15 - 1 / K)) /
      (1 + exp(gp_Eh / k * (1 / gp_Th - 1 / K))),

    r_rate = r_c * exp(r_Ea / k * (1 / 293.15 - 1 / K)) /
      (1 + exp(r_Eh / k * (1 / r_Th - 1 / K))),

    net_rate = gp_rate - r_rate)

draw_curves %>%
  pivot_longer(cols = c(gp_rate, r_rate, net_rate),
               names_to = "rate_type", values_to = "rate") %>%
  mutate(
    rate_type = recode(rate_type,
                       gp_rate = "Gross Photosynthesis",
                       r_rate = "Respiration",
                       net_rate = "Net Photosynthesis"),
    rate_type = factor(rate_type, levels = c(
      "Gross Photosynthesis",
      "Respiration",
      "Net Photosynthesis"
    ))
  ) %>%
  ggplot(aes(x = K - 273.15, y = rate, group = draw_id, color = factor(draw_id))) +
  geom_line(alpha = 0.8, linewidth = .75) +
  facet_wrap(~rate_type, #scales = "free_y",
             ncol = 3) +
  labs(x = "Temperature (°C)", y = expression("Rate day"^{-1}), color = "Draw") +
  theme_minimal() +
  theme(legend.position = "none")

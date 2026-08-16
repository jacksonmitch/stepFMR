scenarios <- list(
  two_group_effects = list(
    betas = rbind(
      g1 = c("(Intercept)" = -1.5, x1 = -3, x2 = 2, x3 = 0.5),
      g2 = c("(Intercept)" = 1.5, x1 = 3, x2 = 2, x3 = 0.5)
    ),
    pi = c(0.4, 0.6),
    sigma = c(0.5, 0.5)
  ),
  two_group_effects_poisson = list(
    betas = rbind(
      g1 = c("(Intercept)" = -0.5, x1 = -1.0, x2 = 0.3, x3 = 0.1),
      g2 = c("(Intercept)" = 0.5, x1 = 1.0, x2 = 0.3, x3 = 0.1)
    ),
    pi = c(0.4, 0.6),
    family = "poisson"
  ),
  two_group_effects_binomial = list(
    betas = rbind(
      g1 = c("(Intercept)" = -1.0, x1 = -1.5, x2 = 0.5, x3 = 0.2),
      g2 = c("(Intercept)" = 1.0, x1 = 1.5, x2 = 0.5, x3 = 0.2)
    ),
    pi = c(0.5, 0.5),
    family = "binomial",
    size = 25
  ),
  two_group_four_variables = list(
    betas = rbind(
      g1 = c("(Intercept)" = -1.0, x1 = -1.5, x2 = 0.5, x3 = 0.3, x4 = 0),
      g2 = c("(Intercept)" = 1.0, x1 = 1.5, x2 = 0.5, x3 = 0.3, x4 = 0)
    ),
    pi = c(0.5, 0.5),
    sigma = c(1, 1)
  ),
  three_group_four_variables = list(
    betas = rbind(
      g1 = c("(Intercept)" = -3, x1 = 0.3, x2 = -0.4, x3 = 0.25, x4 = 0),
      g2 = c("(Intercept)" = 0, x1 = 0.6, x2 = -0.4, x3 = 0.25, x4 = 0),
      g3 = c("(Intercept)" = 3, x1 = 0.8, x2 = -0.4, x3 = 0.25, x4 = 0)
    ),
    pi = c(0.2, 0.3, 0.5),
    sigma = c(0.5, 0.5, 0.5)
  ),
  three_group_twelve_variables_gaussian = list(
    betas = rbind(
      g1 = c(
        "(Intercept)" = -2, x1 = 0.4, x2 = 0.2, x3 = -0.5,
        x4 = -1.1, x5 = 1.0, x6 = 0.2, x7 = -0.3,
        x8 = 0.3, x9 = -0.4, x10 = 0.2, x11 = 0,
        x12 = 0
      ),
      g2 = c(
        "(Intercept)" = 0, x1 = 0.8, x2 = 0.5, x3 = 0,
        x4 = -0.5, x5 = 1.4, x6 = 0.3, x7 = -0.2,
        x8 = 0.3, x9 = -0.4, x10 = 0.2, x11 = 0,
        x12 = 0
      ),
      g3 = c(
        "(Intercept)" = 2, x1 = 1.1, x2 = 0.6, x3 = 0.4,
        x4 = 0.1, x5 = 1.8, x6 = 0.4, x7 = -0.1,
        x8 = 0.3, x9 = -0.4, x10 = 0.2, x11 = 0,
        x12 = 0
      )
    ),
    pi = c(0.2, 0.5, 0.3),
    sigma = c(0.5, 0.5, 0.5)
  ),
  four_group_twelve_variables = list(
    betas = rbind(
      g1 = c(
        "(Intercept)" = -2, het1 = 0.4, het2 = 0.2,
        het3 = -0.5, het4 = -1.1, het5 = 1.0,
        hom1 = 0.6, hom2 = -0.3, hom3 = 0.2,
        null1 = 0, null2 = 0, null3 = 0, null4 = 0
      ),
      g2 = c(
        "(Intercept)" = 0, het1 = 0.8, het2 = 0.5, 
        het3 = 0, het4 = -0.5, het5 = 1.4,
        hom1 = 0.6, hom2 = -0.3, hom3 = 0.2,
        null1 = 0, null2 = 0, null3 = 0, null4 = 0
      ),
      g3 = c(
        "(Intercept)" = 2, het1 = 1.1, het2 = 0.6, 
        het3 = 0.4, het4 = 0.1, het5 = 1.8,
        hom1 = 0.6, hom2 = -0.3, hom3 = 0.2,
        null1 = 0, null2 = 0, null3 = 0, null4 = 0
      ),
      g4 = c(
        "(Intercept)" = 1, het1 = 0.9, het2 = 0.3, 
        het3 = 0.7, het4 = -0.3, het5 = -1.6,
        hom1 = 0.6, hom2 = -0.3, hom3 = 0.2,
        null1 = 0, null2 = 0, null3 = 0, null4 = 0
      )
    ),
    pi = c(0.2, 0.3, 0.25, 0.25),
    sigma = c(0.6, 0.5, 0.4, 0.6)
  ),
  four_group_twelve_variables_doubled = list(
    betas = rbind(
      g1 = c(
        "(Intercept)" = -4, het1 = 0.8, het2 = 0.4,
        het3 = -1, het4 = -2.2, het5 = 2,
        hom1 = 1.2, hom2 = -0.6, hom3 = 0.4,
        null1 = 0, null2 = 0, null3 = 0, null4 = 0
      ),
      g2 = c(
        "(Intercept)" = 0, het1 = 1.6, het2 = 1, 
        het3 = 0, het4 = -1, het5 = 2.8,
        hom1 = 1.2, hom2 = -0.6, hom3 = 0.4,
        null1 = 0, null2 = 0, null3 = 0, null4 = 0
      ),
      g3 = c(
        "(Intercept)" = 4, het1 = 2.2, het2 = 1.2, 
        het3 = 0.8, het4 = 0.2, het5 = 3.6,
        hom1 = 1.2, hom2 = -0.6, hom3 = 0.4,
        null1 = 0, null2 = 0, null3 = 0, null4 = 0
      ),
      g4 = c(
        "(Intercept)" = 2, het1 = 1.8, het2 = 0.6, 
        het3 = 1.4, het4 = -0.6, het5 = -3.2,
        hom1 = 1.2, hom2 = -0.6, hom3 = 0.4,
        null1 = 0, null2 = 0, null3 = 0, null4 = 0
      )
    ),
    pi = c(0.2, 0.3, 0.25, 0.25),
    sigma = c(0.6, 0.5, 0.4, 0.6)
  )
)

library(glue)

# Load assortativity results for E9.5 and E13.5 from assortnet::assortment.discrete()
r1  <- 0.104899 # assortativity coefficient for E9.5
r2  <- 0.1103596 # assortativity coefficient for E13.5
se1 <- 0.05149285 # standard error for E9.5 assortativity coefficient (Newman 2003 jackknife)
se2 <- 0.05149693 # standard error for E13.5 assortativity coefficient (Newman 2003 jackknife)

z_stat  <- (r2 - r1) / sqrt(se1^2 + se2^2)
p_value <- 2 * pnorm(-abs(z_stat))

# Print results
cat(glue("\n{'='*60}\n"))
cat("  Cross-stage assortativity z-test (Newman 2003 jackknife)\n")
cat(glue("{'='*60}\n"))
cat(glue("  E9.5  : r = {round(r1, 6)},  SE = {round(se1, 6)},  M = {results[['E9.5']]$M}\n"))
cat(glue("  E13.5 : r = {round(r2, 6)},  SE = {round(se2, 6)},  M = {results[['E13.5']]$M}\n"))
cat(glue("  Delta r (E13.5 - E9.5) = {round(r2 - r1, 6)}\n"))
cat(glue("  z-statistic = {round(z_stat, 4)}\n"))
cat(glue("  p-value (two-sided) = {signif(p_value, 4)}\n"))
cat(glue("{'='*60}\n\n"))
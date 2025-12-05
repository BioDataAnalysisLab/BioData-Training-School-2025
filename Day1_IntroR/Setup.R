###############################################
######### R Environment Test Script ###########
###############################################

required <- c("ggplot2","dplyr","readr","tidyr","knitr","rmarkdown")
results <- list()

cat("Checking required packages...\n")

# Install + load packages
for (pkg in required) {
  results[[pkg]] <- FALSE
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg)
  }
  
  # Try loading
  ok <- tryCatch({
    library(pkg, character.only = TRUE)
    TRUE
  }, error = function(e) FALSE)
  
  results[[pkg]] <- ok
}

# Simple R functionality test
results[["basic_math"]] <- (2 + 2 == 4)

# Data frame test
results[["data_frame"]] <- tryCatch({
  df <- data.frame(a = 1:3, b = 3:1)
  stopifnot(nrow(df) == 3)
  TRUE
}, error = function(e) FALSE)

# Plot test
results[["plotting"]] <- tryCatch({
  ggplot(data.frame(x = 1:5, y = 1:5), aes(x, y)) + geom_point()
  TRUE
}, error = function(e) FALSE)

###############################################
# PASS / FAIL SUMMARY
###############################################

cat("\n------ SUMMARY ------\n")
summary_lines <- c()
for (name in names(results)) {
  status <- ifelse(results[[name]], "PASS", "FAIL")
  line <- sprintf("%-12s : %s", name, status)
  summary_lines <- c(summary_lines, line)
  cat(line, "\n")
}

overall_pass <- all(unlist(results))

if (overall_pass) {
  cat("\n🎉 All tests passed! You're ready for the workshop.\n")
} else {
  cat("\n⚠️ Some tests failed. Please ask an instuctor for help.\n")
}

###############################################
# WRITE SUMMARY TO FILE
###############################################

outfile <- "R_environment_test_summary.txt"
writeLines(summary_lines, outfile)

cat("\nSummary written to:", normalizePath(outfile), "\n")
###############################################

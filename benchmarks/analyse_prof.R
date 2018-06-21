#!/usr/bin/env Rscript

# print some stats from Rprof output

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
    stop("Usage analyse_prof.R <PROFILE_FILE>")
}
prof_file <- args[1]
cat("Loading profile file:", prof_file, "\n")

library(proftools)

pd <- readProfileData(filename = prof_file)

cat("\nFunction summary (top 20)...\n")
head(funSummary(pd), 20)

cat("\nCall summary (top 20)...\n")
head(callSummary(pd), 20)

cat("\nHot paths (>10%)...\n")
hotPaths(pd, total.pct=10.0)

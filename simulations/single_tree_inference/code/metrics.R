# Metrics Functions for Single-Tree Inference Study
# Created: 2026-04-29
#
# Purpose: Compute similarity and coverage metrics
#
# Study 1 (Similarity) metrics:
# - structure_match(): Check if two tree structures are identical
# - leaf_rmse(): Compute RMSE between full-sample and cross-fitted leaf values
# - max_absolute_diff(): Max |μ_full - μ_crossfit| across leaves
#
# Study 2 (Coverage) metrics:
# - compute_bias_estimate(): B_hat = theta_full - theta_crossfit
# - construct_ci(): Build CI with specified adjustment method
# - check_coverage(): Does CI contain true theta?

# TODO: Implement metric functions here

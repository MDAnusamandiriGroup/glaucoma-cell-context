"""Numerically stable, unsigned Wakefield ABF colocalization.

Equations follow Giambartolomei et al. 2014 and the uniform-prior core of
the official coloc implementation. This is a Python implementation, not
an execution of the R coloc package. No network calls occur on import.
"""
from __future__ import annotations

import numpy as np
from scipy.special import logsumexp
from scipy.stats import t as t_distribution


def derived_standard_error(beta, nominal_p, nominal_df):
    beta = np.asarray(beta, dtype=float)
    nominal_p = np.asarray(nominal_p, dtype=float)
    nominal_df = np.asarray(nominal_df, dtype=float)
    if np.any(~np.isfinite(beta)):
        raise ValueError("Non-finite beta")
    if np.any((nominal_p <= 0) | (nominal_p >= 1) | ~np.isfinite(nominal_p)):
        raise ValueError("Require 0 < nominal p < 1; never replace truncated p-values")
    if np.any((nominal_df <= 0) | ~np.isfinite(nominal_df)):
        raise ValueError("Invalid nominal-test degrees of freedom")
    absolute_t = t_distribution.isf(nominal_p / 2, nominal_df)
    result = np.abs(beta) / absolute_t
    if np.any((result <= 0) | ~np.isfinite(result)):
        raise ValueError("Beta/p pair cannot yield a valid standard error")
    return result


def log_abf(beta, standard_error, prior_sd):
    beta = np.asarray(beta, dtype=float)
    standard_error = np.asarray(standard_error, dtype=float)
    if beta.shape != standard_error.shape or beta.ndim != 1 or len(beta) == 0:
        raise ValueError("Require nonempty aligned 1D beta and SE")
    if np.any(~np.isfinite(beta)) or np.any(~np.isfinite(standard_error)) or np.any(standard_error <= 0):
        raise ValueError("Invalid beta/SE")
    if not np.isfinite(prior_sd) or prior_sd <= 0:
        raise ValueError("Invalid prior SD")
    variance = standard_error**2
    prior_variance = float(prior_sd)**2
    shrinkage = prior_variance / (variance + prior_variance)
    # Algebraically log(1-r), expressed without cancellation at r ~= 1.
    result = 0.5 * (np.log(variance) - np.log(variance + prior_variance) + shrinkage * (beta / standard_error)**2)
    if np.any(~np.isfinite(result)):
        raise ValueError("Non-finite log Bayes factor")
    return result


def combine_log_abf(l1, l2, p1=1e-4, p2=1e-4, p12=1e-5):
    l1 = np.asarray(l1, dtype=float)
    l2 = np.asarray(l2, dtype=float)
    if l1.shape != l2.shape or l1.ndim != 1 or len(l1) < 2:
        raise ValueError("Require at least two aligned regional variants")
    if np.any(~np.isfinite(l1)) or np.any(~np.isfinite(l2)):
        raise ValueError("Non-finite log Bayes factor")
    if not (0 < p12 <= min(p1, p2) < 1):
        raise ValueError("Invalid per-variant priors")
    if len(l1) * max(p1, p2, p12) >= 1:
        raise ValueError("Region too large for frozen per-variant priors")
    # H3 sums pairs with different SNPs. Prefix/suffix log sums avoid loss
    # of precision from subtracting two nearly equal huge Bayes factors.
    prefix = np.logaddexp.accumulate(l2)
    suffix = np.logaddexp.accumulate(l2[::-1])[::-1]
    left = np.concatenate(([-np.inf], prefix[:-1]))
    right = np.concatenate((suffix[1:], [-np.inf]))
    excluding_same = np.logaddexp(left, right)
    log_h3_sum = logsumexp(l1 + excluding_same)
    log_weights = np.array([
        0.0,
        np.log(p1) + logsumexp(l1),
        np.log(p2) + logsumexp(l2),
        np.log(p1) + np.log(p2) + log_h3_sum,
        np.log(p12) + logsumexp(l1 + l2),
    ])
    posterior = np.exp(log_weights - logsumexp(log_weights))
    conditional_snp_h4 = np.exp(l1 + l2 - logsumexp(l1 + l2))
    if not np.isclose(posterior.sum(), 1, atol=1e-12):
        raise RuntimeError("Posterior normalization failed")
    return posterior, conditional_snp_h4

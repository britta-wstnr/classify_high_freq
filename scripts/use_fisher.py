import scipy.stats as stats

# hits and misses split as follows:
# RF: 2835 hits, 1435 misses
# SVM: 2922 hits, 1348 misses

rf_h = 2835
rf_m = 1435
svm_h_rbf = 2922  # svm rbf kernel
svm_m_rbf = 1348  # svm rbf kernel
svm_h = 2697      # svm linear kernel
svm_m = 1573      # svm linear kernel

# cont = ([[rf_h, svm_h], [rf_m, svm_m]])
cont = ([[svm_h, rf_h], [svm_m, rf_m]])

odds_rat, p_val = stats.fisher_exact(cont)

print("P-value of %.3f with an odds ratio of %.2f %%"
      % (p_val, odds_rat * 100))

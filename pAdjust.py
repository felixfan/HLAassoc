#!/usr/bin/env python

def adjustP(pvalues, method = "Benjamini-Hochberg"):                
	"""                                                                                                   
	correct p-values for multiple testing
	methods: Bonferroni, Bonferroni-Holm or Holm, Benjamini-Hochberg or FDR
	"""
	n = len(pvalues)
	cp = [1]*n
	if method == "Bonferroni":
		cp = map(lambda x:min(x*n,1.0), pvalues)
	elif method == "Bonferroni-Holm" or method == "Holm":
		values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
		values = sorted(values)
		for rank, vals in enumerate(values):
			pvalue, i = vals
			cp[i] = (n-rank) * pvalue
		for rank, vals in enumerate(values):
			pvalue, i = vals                                                      
			if rank > 0:
				cp[i] = min(1.0, max(cp[i], cp[j]))
			else:
				cp[i] = min(1.0, cp[i])
			j = i
	elif method == "Benjamini-Hochberg" or method == "FDR":
		values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
		values = sorted(values,reverse=True)
		for rank, vals in enumerate(values):
			pvalue, i = vals
			cp[i] = n * pvalue / (n-rank)
		for rank, vals in enumerate(values):
			pvalue, i = vals
			if rank > 0:
				cp[i] = min(1.0, min(cp[i], cp[j]))
			else:
				cp[i] = min(1.0, cp[i])
			j = i
	elif method == "Benjamini-Yekutieli" or method == "FDR_BY":
		q = 0
		for i in range(1,n+1):
			q += 1.0 / i
		values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
		values = sorted(values,reverse=True)
		for rank, vals in enumerate(values):
			pvalue, i = vals
			cp[i] = q * pvalue * n/(n-rank)
		for rank, vals in enumerate(values):
			pvalue, i = vals
			if rank > 0:
				cp[i] = min(1.0, min(cp[i], cp[j]))
			else:
				cp[i] = min(1.0, cp[i])
			j = i
	return cp


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
			j = i
    return cp


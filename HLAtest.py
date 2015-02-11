from __future__ import division
import scipy.stats
import math
import HLAcount

def assocADRChiFisher(caseAlleles, ctrlAlleles, np, nc, allele, freq, test):
	'''
	Association Analysis for Allelic, Dominant or Recessive Model
	Pearson's Chi-squared test or Fisher exact test
	return allele counts, frequency, [chi-square, df,] p, and OR
	'''
	assoc = {}
	case = {}
	ctrl = {}
	### counts
	for a in caseAlleles:
		if a.startswith(allele):
			case[a] = caseAlleles[a]
	for a in ctrlAlleles:
		if a.startswith(allele):
			ctrl[a] = ctrlAlleles[a]
	### freq
	freqCase = {}
	freqCtrl = {}
	for a in case:
		freqCase[a] = 1.0 * case[a] / np[allele]
	for a in ctrl:
		freqCtrl[a] = 1.0 * ctrl[a] / nc[allele]
	### assoc
	for a in case:
		if a in ctrl:
			if freqCase[a] > freq or freqCtrl[a] > freq:
				n1 = case[a]
				n2 = np[allele] - case[a]
				n3 = ctrl[a]
				n4 = nc[allele] - ctrl[a]
				data = [[n1, n2], [n3, n4]]

				if test == "chisq":
					chi2, p, dof, expected = scipy.stats.chi2_contingency(data)

				OR, pvalue = scipy.stats.fisher_exact(data)

				se = math.sqrt(1.0/n1  + 1.0/n2 +  1.0/n3 + 1.0/n4)
				l95 = math.exp(math.log(OR) - 1.96 * se)
				u95 = math.exp(math.log(OR) + 1.96 * se)

				s1 = a + '\t' + str(n1) + '\t' + str(n2) + '\t' + str(n3) + '\t' + str(n4)
				s2 = '\t' + str(round(freqCase[a],4)) + '\t' + str(round(freqCtrl[a],4))
				if test == "chisq":
					s3 = '\t' + str(round(chi2,4)) + '\t' + str(dof) + '\t' + str(round(p,6))
				elif test == "fisher":
					s3 = '\t' + str(round(pvalue,6))
				s4 = '\t' + str(round(OR,4)) + '\t' + str(round(l95,4)) + '\t' + str(round(u95,4))
				assoc[a] = s1 + s2 + s3 + s4
	return assoc

def runAssoc(infile, digit, freq, model, test):
	geno = {}  # get all genes name
	f = open(infile)
	for line in f:
		alleles = line.split()
		for allele in alleles[2:]:
			gene = allele.split('*')
			geno[gene[0]] = 1
		break
	f.close()

	if model == 'allelic':
		caseAlleles, ctrlAlleles, np, nc = HLAcount.allelicCount(infile,digit)
	elif model == 'dom':
		caseAlleles, ctrlAlleles, np, nc = HLAcount.domCount(infile,digit)
	elif model == 'rec':
		caseAlleles, ctrlAlleles, np, nc = HLAcount.recCount(infile,digit)
	

	result = []
	gs = geno.keys()
	gs = sorted(gs)
	for g in gs:
		assocA = assocADRChiFisher(caseAlleles, ctrlAlleles, np, nc, g, freq, test)
		result.append(assocA)
			
	return result

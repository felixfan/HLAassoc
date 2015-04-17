#!/usr/bin/env python

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
	freqAll = {}
	for a in case:
		freqCase[a] = 1.0 * case[a] / np[allele]
	for a in ctrl:
		freqCtrl[a] = 1.0 * ctrl[a] / nc[allele]
		if a in case:
			freqAll[a] = 1.0 * (case[a] + ctrl[a]) / (np[allele] + nc[allele])
	### assoc
	for a in case:
		if a in ctrl:
			# if freqCase[a] > freq or freqCtrl[a] > freq:
			if freqAll[a] > freq:
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
				ss = []
				ss.append(a)
				ss.append(n1)
				ss.append(n2)
				ss.append(n3)
				ss.append(n4)
				ss.append(freqCase[a])
				ss.append(freqCtrl[a])
				ss.append(freqAll[a])
				if test == "chisq":
					ss.append(chi2)
					ss.append(dof)
					ss.append(p)
				elif test == "fisher":
					ss.append(pvalue)
				ss.append(OR)
				ss.append(l95)
				ss.append(u95)
				assoc[a] = ss
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

def assocRaw(caseAlleles, ctrlAlleles, np, nc, freq, test):
	'''
	Association Analysis (2 x m)
	Pearson's Chi-squared test
	return chi-square, df, p
	'''
	assoc = {}
	### genes
	gene = {}  # get all genes name
	for a in caseAlleles:
		temp = a.split('*')
		gene[temp[0]] = 1
	for g in gene:
		### counts
		case = {}
		ctrl = {}
		for a in caseAlleles:
			if a.startswith(g):
				case[a] = caseAlleles[a]
		for a in ctrlAlleles:
			if a.startswith(g):
				ctrl[a] = ctrlAlleles[a]
		### freq
		freqCase = {}
		freqCtrl = {}
		freqAll = {}
		for a in case:
			freqCase[a] = 1.0 * case[a] / np[g]
		for a in ctrl:
			freqCtrl[a] = 1.0 * ctrl[a] / nc[g]
			if a in case:
				freqAll[a] = 1.0 * (case[a] + ctrl[a]) / (np[g] + nc[g])
		### assoc
		n1 = []
		n2 = []
		for a in case:
			if a in ctrl:
				# if freqCase[a] > freq or freqCtrl[a] > freq:
				if freqAll[a] > freq:
					n1.append(case[a])
					n2.append(ctrl[a])
		data = [n1, n2]
		if test == "raw":
			chi2, p, dof, expected = scipy.stats.chi2_contingency(data)
		ss = []
		if not isinstance(chi2, float):
			ss.append('NA')
		else:
			ss.append(chi2)
		if not isinstance(dof, int):
			ss.append('NA')
		else:
			ss.append(dof)
		if not isinstance(p, float):
			ss.append('NA')
		else:
			ss.append(p)
		assoc[g] = ss
	return assoc
def assocScoreU(caseAlleles, ctrlAlleles, np, nc, freq, test):
	'''
	Association Analysis (2 x m)
	Score test
	return score test U
	'''
	assoc = {}
	### genes
	gene = {}  # get all genes name
	for a in caseAlleles:
		temp = a.split('*')
		gene[temp[0]] = 1
	for g in gene:
		### counts
		case = {}
		ctrl = {}
		for a in caseAlleles:
			if a.startswith(g):
				case[a] = caseAlleles[a]
		for a in ctrlAlleles:
			if a.startswith(g):
				ctrl[a] = ctrlAlleles[a]
		### freq
		freqCase = {}
		freqCtrl = {}
		freqAll = {}
		n1 = 0
		n2 = 0
		for a in case:
			freqCase[a] = 1.0 * case[a] / np[g]
		for a in ctrl:
			freqCtrl[a] = 1.0 * ctrl[a] / nc[g]
			if a in case:
				freqAll[a] = 1.0 * (case[a] + ctrl[a]) / (np[g] + nc[g])
		### score test U
		n1 = np[g]
		u = 0
		for a in freqAll:
			# if freqCase[a] > freq or freqCtrl[a] > freq:
			if freqAll[a] > freq:
				u = u + (case[a] - n1 * freqAll[a]) ** 2 / freqAll[a] - (case[a] - n1 * freqAll[a]) / freqAll[a]
		if not isinstance(u, float):
			u = 'NA'
		assoc[g] = u
	return assoc
#!/usr/bin/env python

import scipy.stats
import random
import pandas as pd
import statsmodels.formula.api as smf
import os
import sys
import HLAcountPerm
import HLAcount
import HLArecode

def chisqFisherPerm(infile, digit, model, test, perm, seed):
	'''
	permutation
	input: genotype, digit to test, genetic model, chisq or Fisher test, number of permutation to run.
	output: dictionary, key:allele, value: p-value
	'''
	# original test
	origP = {}
	if model == 'allelic':
		case, ctrl, np, nc = HLAcount.allelicCount(infile,digit)
	elif model == 'dom':
		case, ctrl, np, nc = HLAcount.domCount(infile,digit)
	elif model == 'rec':
		case, ctrl, np, nc = HLAcount.recCount(infile,digit)
	for a in case:
		if a in ctrl:
			n1 = case[a]
			n2 = np[a.split('*')[0]] - n1
			n3 = ctrl[a]
			n4 = nc[a.split('*')[0]] - n3
			data = [[n1, n2], [n3, n4]]

			if test == "chisq":
				chi2, p, dof, expected = scipy.stats.chi2_contingency(data)
			elif test == 'fisher':
				OR, p = scipy.stats.fisher_exact(data)
			origP[a] = p
	# premutation
	random.seed(seed)
	permP = {}
	pf = perm / 10
	for i in range(perm):
		if i % pf == 1:
			print 'permutation {}/{} ...'.format(i, perm)
		if model == 'allelic':
			case, ctrl, np, nc = HLAcountPerm.allelicCount(infile,digit)
		elif model == 'dom':
			case, ctrl, np, nc = HLAcountPerm.domCount(infile,digit)
		elif model == 'rec':
			case, ctrl, np, nc = HLAcountPerm.recCount(infile,digit)
		for a in case:
			if a in ctrl:
				n1 = case[a]
				n2 = np[a.split('*')[0]] - n1
				n3 = ctrl[a]
				n4 = nc[a.split('*')[0]] - n3
				data = [[n1, n2], [n3, n4]]

				if test == "chisq":
					chi2, p, dof, expected = scipy.stats.chi2_contingency(data)
				elif test == 'fisher':
					OR, p = scipy.stats.fisher_exact(data)
				if a in origP:
					if isinstance(p, float) and p < origP[a]:
						if a in permP:
							permP[a] += 1
						else:
							permP[a] = 1
	for a in origP:
		if a in permP:
			permP[a] = 1.0 * (permP[a] + 1) / (perm + 1)
		else:
			permP[a] = 'NA'
	return permP

def regressionPerm(infile, digits, method,perm,seed):
	'''
	linear regression or logistitic regression
	output: dictionary, key: allele, value: p-value
	'''
	tfile = HLArecode.writeRecode(infile, digits, method)
	geno = pd.read_csv(tfile,delim_whitespace= True, header = 0)
	os.remove(tfile)
	alleles = list(geno.columns.values)[2:]
	# original test
	assoc = {}
	for allele in alleles:
		myformula = 'PHT ~ ' + allele
		if method == 'logistic':
			try:
				lr = smf.logit(formula = myformula, data = geno).fit(maxiter=100, disp=False)
				p = lr.pvalues[1]	
			except:
				p = 'NA'
		elif method == 'linear':
			try:
				lr = smf.ols(formula = myformula, data = geno).fit(maxiter=100, disp=False)
				p = lr.pvalues[1]
			except:
				p = 'NA'
		aname = allele.split('_')
		nname = aname[0] + '*' + aname[1]
		if digits == 4:
			nname = nname + ':' + aname[2]
		elif digits == 6:
			nname = nname + ':' + aname[2] + ':' + aname[3]
		assoc[nname] = p
	# premutation
	random.seed(seed)
	permP = {}
	pf = perm / 10
	for i in range(perm):
		if i % pf == 1:
			print 'permutation {}/{} ...'.format(i, perm)
		random.shuffle(geno['PHT'])
		for allele in alleles:
			aname = allele.split('_')
			nname = aname[0] + '*' + aname[1]
			if digits == 4:
				nname = nname + ':' + aname[2]
			elif digits == 6:
				nname = nname + ':' + aname[2] + ':' + aname[3]
			if nname in assoc:
				myformula = 'PHT ~ ' + allele
				if method == 'logistic':
					try:
						lr = smf.logit(formula = myformula, data = geno).fit(maxiter=100, disp=False)
						p = lr.pvalues[1]	
					except:
						p = 'NA'
				elif method == 'linear':
					try:
						lr = smf.ols(formula = myformula, data = geno).fit(maxiter=100, disp=False)
						p = lr.pvalues[1]
					except:
						p = 'NA'
				if p != 'NA' and assoc[nname] != 'NA' and p < assoc[nname]:
					if nname in permP:
						permP[nname] += 1
					else:
						permP[nname] = 1
	for a in assoc:
		if a in permP:
			permP[a] = 1.0 * (permP[a] + 1) / (perm + 1)
		else:
			permP[a] = 'NA'
	return permP

def regressionCovPerm(infile, digits, method, perm, covfile, covname,seed):
	'''
	linear regression or logistitic regression with covariants
	output: dictionary, key: allele, value: p-value
	'''
	tfile = HLArecode.writeRecode(infile, digits,method)
	geno = pd.read_csv(tfile,delim_whitespace= True, header = 0)
	cov = pd.read_csv(covfile,delim_whitespace= True, header = 0)
	os.remove(tfile)
	alleles = list(geno.columns.values)[2:]
	covindex = list(cov.columns.values)[1:]
	if not covname:                           # default: use all covariants
		covname = covindex
	else:
		covname = covname.split(',')
	# original test
	assoc = {}
	for allele in alleles:
			geno9 = geno.ix[:, ['IID', 'PHT', allele]]
			mydata = pd.merge(geno9, cov, on='IID', how='inner')
			myformula = 'PHT ~ ' + allele
			for name in covname:
				if name in covindex:
					myformula = myformula + ' + ' + name
				else:
					print 'can not find covariant name ' + '"' + name + '" in covariant file'
					sys.exit()
			if method == 'logistic':
				try:
					lr = smf.logit(formula = myformula, data = mydata).fit(maxiter=100, disp=False)
					p = lr.pvalues[1]
				except:
					p = 'NA'
			elif method == 'linear':
				try:
					lr = smf.ols(formula = myformula, data = mydata).fit(maxiter=100, disp=False)
					p = lr.pvalues[1]
				except:
					p = 'NA'
			aname = allele.split('_')
			nname = aname[0] + '*' + aname[1]
			if digits == 4:
				nname = nname + ':' + aname[2]
			elif digits == 6:
				nname = nname + ':' + aname[2] + ':' + aname[3]
			assoc[nname] = p
	# permutation test
	random.seed(seed)
	permP = {}
	pf = perm / 10
	for i in range(perm):
		if i % pf == 1:
			print 'permutation {}/{} ...'.format(i, perm)
		random.shuffle(geno['PHT'])
		for allele in alleles:
			aname = allele.split('_')
			nname = aname[0] + '*' + aname[1]
			if digits == 4:
				nname = nname + ':' + aname[2]
			elif digits == 6:
				nname = nname + ':' + aname[2] + ':' + aname[3]
			if nname in assoc:
				geno9 = geno.ix[:, ['IID', 'PHT', allele]]
				mydata = pd.merge(geno9, cov, on='IID', how='inner')
				myformula = 'PHT ~ ' + allele
				for name in covname:
					if name in covindex:
						myformula = myformula + ' + ' + name
					else:
						print 'can not find covariant name ' + '"' + name + '" in covariant file'
						sys.exit()
				if method == 'logistic':
					try:
						lr = smf.logit(formula = myformula, data = mydata).fit(maxiter=100, disp=False)
						p = lr.pvalues[1]
					except:
						p = 'NA'
				elif method == 'linear':
					try:
						lr = smf.ols(formula = myformula, data = mydata).fit(maxiter=100, disp=False)
						p = lr.pvalues[1]
					except:
						p = 'NA'
				if p != 'NA' and assoc[nname] != 'NA' and p < assoc[nname]:
					if nname in permP:
						permP[nname] += 1
					else:
						permP[nname] = 1
	for a in assoc:
		if a in permP:
			permP[a] = 1.0 * (permP[a] + 1) / (perm + 1)
		else:
			permP[a] = 'NA'
	return permP
def rawPerm(infile, digit, perm, seed, freq):
	'''
	permutation
	input: genotype, digit to test, number of permutation to run.
	output: dictionary, key:gene, value: p-value
	'''
	# original test
	caseAlleles, ctrlAlleles, np, nc = HLAcount.allelicCount(infile,digit)
	origP = {}
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
		### counts
		n1 = []
		n2 = []
		for a in case:
			if a in ctrl:
				if freqCase[a] > freq or freqCtrl[a] > freq:
					n1.append(case[a])
					n2.append(ctrl[a])
		data = [n1, n2]
		chi2, p, dof, expected = scipy.stats.chi2_contingency(data)
		if not isinstance(p, float):
			origP[g] = 'NA'
		else:
			origP[g] = p
	# premutation
	random.seed(seed)
	permP = {}
	pf = perm / 10
	for i in range(perm):
		if i % pf == 1:
			print 'permutation {}/{} ...'.format(i, perm)
		caseAlleles, ctrlAlleles, np, nc = HLAcountPerm.allelicCount(infile,digit)
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
			### counts
			n1 = []
			n2 = []
			for a in case:
				if a in ctrl:
					if freqCase[a] > freq or freqCtrl[a] > freq:
						n1.append(case[a])
						n2.append(ctrl[a])
			data = [n1, n2]
			chi2, p, dof, expected = scipy.stats.chi2_contingency(data)

			if g in origP:
				if isinstance(p, float) and origP[g] != 'NA' and p < origP[g]:
					if g in permP:
						permP[g] += 1
					else:
						permP[g] = 1
	for a in origP:
		if a in permP:
			permP[a] = 1.0 * (permP[a] + 1) / (perm + 1)
		else:
			permP[a] = 'NA'
	return permP
def scorePerm(infile, digit, perm, seed, freq):
	'''
	permutation
	input: genotype, digit to test, number of permutation to run.
	output: dictionary, key:gene, value: p-value
	'''
	# original test
	caseAlleles, ctrlAlleles, np, nc = HLAcount.allelicCount(infile,digit)
	origP = {}
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
				n1 += case[a]
				n2 += ctrl[a]
		### score test U
		u = 0
		for a in freqAll:
			if freqCase[a] > freq or freqCtrl[a] > freq:
				u = u + (case[a] - n1 * freqAll[a]) ** 2 / freqAll[a] - (case[a] - n1 * freqAll[a]) / freqAll[a]
		if not isinstance(u, float):
			origP[g] = 'NA'
		else:
			origP[g] = u
	# premutation
	random.seed(seed)
	permP = {}
	pf = perm / 10
	for i in range(perm):
		if i % pf == 1:
			print 'permutation {}/{} ...'.format(i, perm)
		caseAlleles, ctrlAlleles, np, nc = HLAcountPerm.allelicCount(infile,digit)
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
					n1 += case[a]
					n2 += ctrl[a]
			### score test U
			u = 0
			for a in freqAll:
				if freqCase[a] > freq or freqCtrl[a] > freq:
					u = u + (case[a] - n1 * freqAll[a]) ** 2 / freqAll[a] - (case[a] - n1 * freqAll[a]) / freqAll[a]
			
			if g in origP:
				if isinstance(u, float) and origP[g] != 'NA' and u > origP[g]:
					if g in permP:
						permP[g] += 1
					else:
						permP[g] = 1
	for a in origP:
		if a in permP:
			permP[a] = 1.0 * (permP[a] + 1) / (perm + 1)
		else:
			permP[a] = 'NA'
	return permP
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

def chisqFisherPerm(infile, digit, model, test, perm, seed, freq):
	'''
	permutation
	input: genotype, digit to test, genetic model, chisq or Fisher test, number of permutation to run.
	output: dictionary, key:allele, value: p-value
	'''
	# original test
	origP = {}
	case = {}
	ctrl = {}
	np = {}
	nc = {}
	comalleles = []
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
			if 1.0 * (n1 + n3) / (n1 + n2 + n3 + n4) > freq:
            			comalleles.append(a)
				if test == "chisq":
					chi2, p, dof, expected = scipy.stats.chi2_contingency(data)
				elif test == 'fisher':
					OR, p = scipy.stats.fisher_exact(data)
				origP[a] = p
	# premutation
	random.seed(seed)
	permP = {}
	permN = {}
	permNL = {}
	permNA = {}
	for a in origP:
		permNA[a] = 0
		permNL[a] = 0
		permN[a] = 0
	if perm >= 10:
		pf = perm / 10
	else:
		pf = 2
	pn = 0
	while 1:
		if model == 'allelic':
			case9, ctrl9, np9, nc9 = HLAcountPerm.allelicCount(infile,digit)
		elif model == 'dom':
			case9, ctrl9, np9, nc9 = HLAcountPerm.domCount(infile,digit)
		elif model == 'rec':
			case9, ctrl9, np9, nc9 = HLAcountPerm.recCount(infile,digit)
		ca = []
		for a in case9:
			if a in ctrl9:
				ca.append(a)
		if set(comalleles) <= set(ca):
        		for a in comalleles:
				n1 = case9[a]
				n2 = np9[a.split('*')[0]] - n1
				n3 = ctrl9[a]
				n4 = nc9[a.split('*')[0]] - n3
				data = [[n1, n2], [n3, n4]]
				if test == "chisq":
					chi2, p, dof, expected = scipy.stats.chi2_contingency(data)
				elif test == 'fisher':
					OR, p = scipy.stats.fisher_exact(data)
				if not isinstance(p, float):
					permNA[a] += 1
				else:
					if origP[a] == 'NA':
						permNA[a] += 1
					else:
						if p < origP[a]:
							permN[a] += 1
						else:
							permNL[a] += 1
			pn += 1
			if pn % pf == 1:
				print 'permutation {}/{} ...'.format(pn, perm)
		if pn == perm:
			break
	for a in origP:
		if origP[a] == 'NA':
			permP[a] = 'NA'
		else:
			if permNA[a] == perm:
				permP[a] = 'NA'
			else:
				permP[a] = 1.0 * (permN[a] + 1) / (perm + 1 - permNA[a])
				# permP[a] = 1.0 * (permN[a] + 1) / (permN[a] + permNL[a] + 1)
	return permP, permN, permNA, permNL

def regressionPerm(infile, digits, method,perm,seed, freq):
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
		n1 = int(geno[allele].sum(axis=0)) # number of allele
		n2 = int(geno[allele].count()) * 2 # total allele
		f12 = 1.0 * n1 / n2
		if f12 > freq:
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
	permN = {}
	permNL = {}
	permNA = {}
	for a in assoc:
		permNA[a] = 0
		permNL[a] = 0
		permN[a] = 0
	if perm >= 10:
		pf = perm / 10
	else:
		pf = 2
	for i in range(perm):
		if i % pf == 1:
			print 'permutation {}/{} ...'.format(i, perm)
		# random.shuffle(geno['PHT'])
		xxx = geno['PHT'].values.flatten()
		random.shuffle(xxx)
		geno['PHT'] = xxx
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
						# lr = smf.ols(formula = myformula, data = geno).fit(maxiter=100, disp=False)
						lr = smf.ols(formula = myformula, data = geno).fit()
						p = lr.pvalues[1]
					except:
						p = 'NA'
				if p == 'NA':
					permNA[nname] += 1
				else:
					if assoc[nname] == 'NA':
						permNA[nname] += 1
					else:
						if p < assoc[nname]:
							permN[nname] += 1
						else:
							permNL[nname] += 1
	for a in assoc:
		if assoc[a] == 'NA':
			permP[a] = 'NA'
		else:
			if permNA[a] == perm:
				permP[a] = 'NA'
			else:
				permP[a] = 1.0 * (permN[a] + 1) / (perm + 1 - permNA[a])
	return permP, permN, permNA, permNL

def regressionCovPerm(infile, digits, method, perm, covfile, covname,seed, freq):
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
		n1 = int(geno[allele].sum(axis=0)) # number of allele
		n2 = int(geno[allele].count()) * 2 # total allele
		f12 = 1.0 * n1 / n2
		if f12 > freq:
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
	permN = {}
	permNL = {}
	permNA = {}
	for a in assoc:
		permNA[a] = 0
		permNL[a] = 0
		permN[a] = 0
	if perm >= 10:
		pf = perm / 10
	else:
		pf = 2
	for i in range(perm):
		if i % pf == 1:
			print 'permutation {}/{} ...'.format(i, perm)
		# random.shuffle(geno['PHT'])
		xxx = geno['PHT'].values.flatten()
		random.shuffle(xxx)
		geno['PHT'] = xxx
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
				if p == 'NA':
					permNA[nname] += 1
				else:
					if assoc[nname] == 'NA':
						permNA[nname] += 1
					else:
						if p < assoc[nname]:
							permN[nname] += 1
						else:
							permNL[nname] += 1
	for a in assoc:
		if assoc[a] == 'NA':
			permP[a] = 'NA'
		else:
			if permNA[a] == perm:
				permP[a] = 'NA'
			else:
				permP[a] = 1.0 * (permN[a] + 1) / (perm + 1 - permNA[a])
	return permP, permN, permNA, permNL

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
	commalleles = []
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
		for a in freqAll:
			if freqAll[a] > freq:
				commalleles.append(a)
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
	permN = {}
	permNL = {}
	permNA = {}
	for a in origP:
		permNA[a] = 0
		permNL[a] = 0
		permN[a] = 0
	if perm >= 10:
		pf = perm / 10
	else:
		pf = 2
	pn = 0
	while True:
		ca = []
		caseAlleles9, ctrlAlleles9, np9, nc9 = HLAcountPerm.allelicCount(infile,digit)
		for a in caseAlleles9:
			if a in ctrlAlleles9:
				ca.append(a)
		if set(commalleles) <= set(ca):
			for g in gene:
				### counts
				case = {}
				ctrl = {}
				for a in caseAlleles9:
					if a.startswith(g):
						case[a] = caseAlleles9[a]
				for a in ctrlAlleles9:
					if a.startswith(g):
						ctrl[a] = ctrlAlleles9[a]
				n1 = []
				n2 = []
				for a in commalleles:
					if a.startswith(g):
						n1.append(case[a])
						n2.append(ctrl[a])
				data = [n1, n2]
				if not isinstance(p, float):
					permNA[g] += 1
				else:
					if origP[g] == 'NA':
						permNA[g] += 1
					else:
						if p < origP[g]:
							permN[g] += 1
						else:
							permNL[g] += 1
			pn += 1
			if pn % pf == 1:
				print 'permutation {}/{} ...'.format(pn, perm)
			if pn == perm:
				break
	for a in origP:
		if origP[a] == 'NA':
			permP[a] = 'NA'
		else:
			if permNA[a] == perm:
				permP[a] = 'NA'
			else:
				permP[a] = 1.0 * (permN[a] + 1) / (perm + 1 - permNA[a])
	return permP, permN, permNA, permNL

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
	commalleles = []
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
			if freqAll[a] > freq:
				commalleles.append(a)
				u = u + (case[a] - n1 * freqAll[a]) ** 2 / freqAll[a] - (case[a] - n1 * freqAll[a]) / freqAll[a]
		if not isinstance(u, float):
			origP[g] = 'NA'
		else:
			origP[g] = u
	# premutation
	random.seed(seed)
	permP = {}
	permN = {}
	permNL = {}
	permNA = {}
	for a in origP:
		permNA[a] = 0
		permNL[a] = 0
		permN[a] = 0
	if perm >= 10:
		pf = perm / 10
	else:
		pf = 2
	pn = 0
	while True:
		ca = []
		caseAlleles, ctrlAlleles, np, nc = HLAcountPerm.allelicCount(infile,digit)
		for a in caseAlleles:
			if a in ctrlAlleles:
				ca.append(a)
		if set(commalleles) <= set(ca):
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
				freqAll = {}
				n1 = 0
				n2 = 0
				for a in commalleles:
					if a.startswith(g):
						freqAll[a] = 1.0 * (case[a] + ctrl[a]) / (np[g] + nc[g])
						n1 += case[a]
						n2 += ctrl[a]
				### score test U
				u = 0
				for a in commalleles:
					if a.startswith(g):
						u = u + (case[a] - n1 * freqAll[a]) ** 2 / freqAll[a] - (case[a] - n1 * freqAll[a]) / freqAll[a]				
				if not isinstance(u, float):
					permNA[g] = 'NA'
				else:
					if origP[g] == 'NA':
						permNA[g] = 'NA'
					else:
						if u > origP[g]:
							permN[g] += 1
						else:
							permNL[g] += 1
			pn += 1
			if pn % pf == 1:
				print 'permutation {}/{} ...'.format(pn, perm)
			if pn == perm:
				break
	for a in origP:
		if origP[a] == 'NA':
			permP[a]  = 'NA'
		else:
			if permNA[a] == perm:
				permP[a]  = 'NA'
			else:
				permP[a] = 1.0 * (permN[a] + 1) / (perm + 1 - permNA[a])
	return permP, permN, permNA, permNL

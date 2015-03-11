#!/usr/bin/env python

import pandas as pd
import statsmodels.formula.api as smf
import math
import os
import sys
import HLArecode

def regressionLogistic(infile, digits, freq, method):
	'''
	logistitic regression
	output: dictionary, key: allele, value: statistic includes count, freq, p and OR
	'''
	tfile = HLArecode.writeRecode(infile, digits, method)
	geno = pd.read_csv(tfile,delim_whitespace= True, header = 0)
	os.remove(tfile)
	alleles = list(geno.columns.values)[2:]
	assoc = {}
	for allele in alleles:
		n1 = int(geno[allele].where(geno['PHT'] == 1).sum(axis=0)) # case
		n2 = int(geno[allele].where(geno['PHT'] == 1).count()) * 2 # case
		n3 = int(geno[allele].where(geno['PHT'] == 0).sum(axis=0)) # control
		n4 = int(geno[allele].where(geno['PHT'] == 0).count()) * 2 # control
		f1 = 1.0 * n1 / n2
		f2 = 1.0 * n3 / n4
		f12 = 1.0 * (n1 + n3) / (n2 + n4)
		if f1 > freq or f2 > freq:
			n2 -= n1
			n4 -= n3
			myformula = 'PHT ~ ' + allele
			if method == 'logistic':
				try:
					lr = smf.logit(formula = myformula, data = geno).fit(maxiter=100, disp=False)
					p = lr.pvalues[1]
					try:
						OR = math.exp(lr.params[1])
					except:
						OR = 'NA'
					try:
						L95 = math.exp(lr.conf_int()[0][1])
					except:
						L95 = 'NA'
					try:
						U95 = math.exp(lr.conf_int()[1][1])
					except:
						U95 = 'NA'
				except:
					p = 'NA'
					OR = 'NA'
					L95 = 'NA'
					U95 = 'NA'
			aname = allele.split('_')
			nname = aname[0] + '*' + aname[1]
			if digits == 4:
				nname = nname + ':' + aname[2]
			elif digits == 6:
				nname = nname + ':' + aname[2] + ':' + aname[3]
			s1 = nname + '\t' + str(n1) + '\t' + str(n2) + '\t' + str(n3) + '\t' + str(n4)
			s2 = '\t' + str(round(f1,4)) + '\t' + str(round(f2,4)) + '\t' + str(round(f12,4))

			if p != 'NA':
				s3 = '\t' + str(round(p,6))
			else:
				s3 = '\t' + 'NA'

			if method == "logistic":	
				if OR != 'NA':
					s3 = s3 + '\t' + str(round(OR,4))
				else:
					s3 = s3 + '\t' + 'NA'
			if L95 != 'NA':
				s4 = '\t' + str(round(L95,4))
			else:
				s4 = '\t' + 'NA'
			if U95 != 'NA':
				s4 = s4 + '\t' + str(round(U95,4))
			else:
				s4 = s4 + '\t' + 'NA'
			assoc[nname] = s1 + s2 + s3 + s4
	return assoc
def regressionLinear(infile, digits, freq, method):
	'''
	linear regression
	output: dictionary, key: allele, value: statistic includes p and beta
	'''
	tfile = HLArecode.writeRecode(infile, digits, method)
	geno = pd.read_csv(tfile,delim_whitespace= True, header = 0)
	os.remove(tfile)
	alleles = list(geno.columns.values)[2:]
	assoc = {}
	for allele in alleles:
		n1 = int(geno[allele].sum(axis=0)) # number of allele
		n2 = int(geno[allele].count()) * 2 # total allele
		f12 = 1.0 * n1 / n2
		if f12 > freq:
			myformula = 'PHT ~ ' + allele
			if method == 'linear':
				try:
					lr = smf.ols(formula = myformula, data = geno).fit(maxiter=100, disp=False)
					p = lr.pvalues[1]
					beta = lr.params[1]
					L95 = lr.conf_int()[0][1]
					U95 = lr.conf_int()[1][1]
				except:
					p = 'NA'
					beta = 'NA'
					L95 = 'NA'
					U95 = 'NA'
			aname = allele.split('_')
			nname = aname[0] + '*' + aname[1]
			if digits == 4:
				nname = nname + ':' + aname[2]
			elif digits == 6:
					nname = nname + ':' + aname[2] + ':' + aname[3]
			s1 = nname + '\t' + str(round(f12,4))
			if p != 'NA':
				s3 = '\t' + str(round(p,6))
			else:
				s3 = '\t' + 'NA'

			if method == "linear":
				if beta != 'NA':
					s3 = s3 + '\t' + str(round(beta,4))
				else:
					s3 = s3 + '\t' + 'NA'
			if L95 != 'NA':
				s4 = '\t' + str(round(L95,4))
			else:
				s4 = '\t' + 'NA'
			if U95 != 'NA':
				s4 = s4 + '\t' + str(round(U95,4))
			else:
				s4 = s4 + '\t' + 'NA'
			assoc[nname] = s1 + s3 + s4
	return assoc
def logisticCov(infile, digits, freq, method, covfile, covname):
	'''
	logistitic regression with covariants
	output: dictionary, key: allele, value: statistic includes count, freq, p and OR
	'''
	tfile = HLArecode.writeRecode(infile, digits, method)
	geno = pd.read_csv(tfile,delim_whitespace= True, header = 0)
	cov = pd.read_csv(covfile,delim_whitespace= True, header = 0)
	os.remove(tfile)
	alleles = list(geno.columns.values)[2:]
	covindex = list(cov.columns.values)[1:]
	if not covname:                           # default: use all covariants
		covname = covindex
	else:
		covname = covname.split(',')
	assoc = {}
	for allele in alleles:
		n1 = int(geno[allele].where(geno['PHT'] == 1).sum(axis=0)) # case
		n2 = int(geno[allele].where(geno['PHT'] == 1).count()) * 2 # case
		n3 = int(geno[allele].where(geno['PHT'] == 0).sum(axis=0)) # control
		n4 = int(geno[allele].where(geno['PHT'] == 0).count()) * 2 # control
		f1 = 1.0 * n1 / n2
		f2 = 1.0 * n3 / n4
		f12 = 1.0 * (n1 + n3) / (n2 + n4)
		if f1 > freq or f2 > freq:
			n2 -= n1
			n4 -= n3
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
					try:
						OR = math.exp(lr.params[1])
					except:
						OR = 'NA'
					try:
						L95 = math.exp(lr.conf_int()[0][1])
					except:
						L95 = 'NA'
					try:
						U95 = math.exp(lr.conf_int()[1][1])
					except:
						U95 = 'NA'
				except:
					p = 'NA'
					OR = 'NA'
					L95 = 'NA'
					U95 = 'NA'
			aname = allele.split('_')
			nname = aname[0] + '*' + aname[1]
			if digits == 4:
				nname = nname + ':' + aname[2]
			elif digits == 6:
				nname = nname + ':' + aname[2] + ':' + aname[3]
			s1 = nname + '\t' + str(n1) + '\t' + str(n2) + '\t' + str(n3) + '\t' + str(n4)
			s2 = '\t' + str(round(f1,4)) + '\t' + str(round(f2,4)) + '\t' + str(round(f12,4))
			if p != 'NA':
				s3 = '\t' + str(round(p,6))
			else:
				s3 = '\t' + 'NA'

			if method == "logistic":	
				if OR != 'NA':
					s3 = s3 + '\t' + str(round(OR,4))
				else:
					s3 = s3 + '\t' + 'NA'
			if L95 != 'NA':
				s4 = '\t' + str(round(L95,4))
			else:
				s4 = '\t' + 'NA'
			if U95 != 'NA':
				s4 = s4 + '\t' + str(round(U95,4))
			else:
				s4 = s4 + '\t' + 'NA'
			assoc[nname] = s1 + s2 + s3 + s4
	return assoc
def linearCov(infile, digits, freq, method, covfile, covname):
	'''
	linear regression with covariants
	output: dictionary, key: allele, value: statistic includes freq, p and beta
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
			if method == 'linear':
				try:
					lr = smf.ols(formula = myformula, data = mydata).fit(maxiter=100, disp=False)
					p = lr.pvalues[1]
					beta = lr.params[1]
					L95 = lr.conf_int()[0][1]
					U95 = lr.conf_int()[1][1]
				except:
					p = 'NA'
					beta = 'NA'
					L95 = 'NA'
					U95 = 'NA'
			aname = allele.split('_')
			nname = aname[0] + '*' + aname[1]
			if digits == 4:
				nname = nname + ':' + aname[2]
			elif digits == 6:
				nname = nname + ':' + aname[2] + ':' + aname[3]
			s1 = nname + '\t' + str(round(f12,4))
			if p != 'NA':
				s3 = '\t' + str(round(p,6))
			else:
				s3 = '\t' + 'NA'

			if method == "linear":
				if beta != 'NA':
					s3 = s3 + '\t' + str(round(beta,4))
				else:
					s3 = s3 + '\t' + 'NA'
			if L95 != 'NA':
				s4 = '\t' + str(round(L95,4))
			else:
				s4 = '\t' + 'NA'
			if U95 != 'NA':
				s4 = s4 + '\t' + str(round(U95,4))
			else:
				s4 = s4 + '\t' + 'NA'
			assoc[nname] = s1 + s3 + s4
	return assoc

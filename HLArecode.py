#!/usr/bin/env python

import sys
import time
import string

def getAlleles(infile, digits):
	'''
	get all alleles for each gene
	return a dictionary, keys are the gene names, values are the alleles
	return second dictionary, keys are the gene names, values are the start column
	'''
	geneAlleles = {}
	geneCol = {}
	f = open(infile)
	for line in f:
		line = line.rstrip()
		alleles = line.split()
		for i in range(2,len(alleles)):
			if alleles[i] != 'NA':
				names = alleles[i].split(":")
				gene = names[0].split("*")
				if i not in geneCol:
					geneCol[i] = gene[0]
				if gene[0] not in geneAlleles:
					geneAlleles[gene[0]] = []
				if digits == 4:
					temp = names[0] + ":" + names[1]
					if temp not in geneAlleles[gene[0]]:
						geneAlleles[gene[0]].append(temp)
				elif digits == 2:
					temp = names[0]
					if temp not in geneAlleles[gene[0]]:
						geneAlleles[gene[0]].append(temp)
	f.close()
	return geneAlleles, geneCol

def allelicRecode(infile, digits):
	'''
	allele dosage coding
	Assume A*01:01 is the test allele, then A*01:01 A*01:01 is code as 2
	A*01:01 A*01:02 is code as 1, and A*01:02 A*01:03 is code as 0
	return a list contains the alleles
	return a dictionary contains the coding for alleles
	'''

	geneAlleles, geneCol = getAlleles(infile, digits)

	header = ['IID','PHT']
	ans = {}
	f = open(infile)
	for line in f:
		line = line.rstrip()
		alleles = line.split()
		ans[alleles[0]] = [alleles[0],int(alleles[1])-1]
		for i in range(2,len(alleles),2):
			j = i + 1
			if alleles[i] != 'NA' and alleles[j] != 'NA':
				gene1 = alleles[i].split('*')[0]
				gene2 = alleles[j].split('*')[0]
				if gene1 == gene2:
					if digits == 4:
						allele1 = alleles[i].split(':')[0] + ':' + alleles[i].split(':')[1]
						allele2 = alleles[j].split(':')[0] + ':' + alleles[j].split(':')[1]
					else:
						allele1 = alleles[i].split(':')[0]
						allele2 = alleles[j].split(':')[0]

					gAlleles = sorted(geneAlleles[gene1])
					for ga in gAlleles:
						if ga not in header:
							header.append(ga)
						if allele1 == ga and allele2 == ga:
							ans[alleles[0]].append(2)
						elif allele1 == ga or allele2 == ga:
							ans[alleles[0]].append(1)
						else:
							ans[alleles[0]].append(0)
				else:
					sys.exit("input format is wrong!")
			else:
				for gg in geneAlleles[geneCol[i]]:
					ans[alleles[0]].append('NA')
	return ans,header

def writeRecode(infile, digits):
	'''
	write coding to a temp file
	return the temp file name
	'''
	
	tmp = time.strftime("%H%M%S%d%b%Y")
	tmp = tmp + '.txt'
	f = open(tmp,'w')
	ans, header = allelicRecode(infile, digits)
	for i in header:
		i = string.replace(i, '*', '_') # chang A*01:01 to A_01_01
		i = string.replace(i, ':', '_')
		f.write("%12s" % i,)
	f.write('\n')

	for i in ans:
		for j in ans[i]:
			f.write("%12s" % j,)
		f.write('\n')
	f.close()
	return tmp

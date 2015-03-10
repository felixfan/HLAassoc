#!/usr/bin/env python
import sys

def allelicCount(infile, digits):
	'''
	count all alleles
	return the count of each allele in case and sontrol
	'''
	case = {}      # counts for each allele
	control = {}
	np = {}         # total non-NA alleles for each gene
	nc = {}
	f = open(infile)
	for line in f:
		line = line.rstrip()
		alleles = line.split()
		for i in range(2,len(alleles),2):
			j = i + 1
			if alleles[i] != 'NA' and alleles[j] != 'NA':
				names1 = alleles[i].split(":")
				names2 = alleles[j].split(":")
				if digits == 6:
					if len(names1) < 3 or len(names2) < 3:
						sys.exit("--digits 6 requires at least 6 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1] + ":" + names1[2]
					a4d2 = names2[0] + ":" + names2[1] + ":" + names2[2]
				if digits == 4:
					if len(names1) < 2 or len(names2) < 2:
						sys.exit("--digits 4 requires at least 4 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1]
					a4d2 = names2[0] + ":" + names2[1]
				elif digits == 2:
					if len(names1) < 1 or len(names2) < 1:
						sys.exit("--digits 2 requires at least 2 digits resolution genotype!")
					a4d1 = names1[0]
					a4d2 = names2[0]
				
				if alleles[1] == "2":
					if a4d1.split("*")[0] in np:
						np[a4d1.split("*")[0]] += 1
					else:
						np[a4d1.split("*")[0]] = 1

					if a4d1 in case:
						case[a4d1] += 1
					else:
						case[a4d1] = 1
					##
					if a4d2.split("*")[0] in np:
						np[a4d2.split("*")[0]] += 1
					else:
						np[a4d2.split("*")[0]] = 1

					if a4d2 in case:
						case[a4d2] += 1
					else:
						case[a4d2] = 1
				elif alleles[1] == "1":
					if a4d1.split("*")[0] in nc:
						nc[a4d1.split("*")[0]] += 1
					else:
						nc[a4d1.split("*")[0]] = 1

					if a4d1 in control:
						control[a4d1] += 1
					else:
						control[a4d1] = 1
					##
					if a4d2.split("*")[0] in nc:
						nc[a4d2.split("*")[0]] += 1
					else:
						nc[a4d2.split("*")[0]] = 1

					if a4d2 in control:
						control[a4d2] += 1
					else:
						control[a4d2] = 1
	f.close()
	return case, control, np, nc

def domCount(infile, digits):
	'''
	count all alleles
	return the count of each allele in case and sontrol
	'''
	case = {}     # counts for each allele
	control = {}
	nc = {}
	np = {}
	f = open(infile)
	for line in f:
		line = line.rstrip()
		alleles = line.split()
		for i in range(2,len(alleles),2):
			j = i + 1
			if alleles[i] != 'NA' and alleles[j] != 'NA':
				names1 = alleles[i].split(":")
				names2 = alleles[j].split(":")
				if digits == 6:
					if len(names1) < 3 or len(names2) < 3:
						sys.exit("--digits 6 requires at least 6 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1] + ":" + names1[2]
					a4d2 = names2[0] + ":" + names2[1] + ":" + names2[2]
				if digits == 4:
					if len(names1) < 2 or len(names2) < 2:
						sys.exit("--digits 4 requires at least 4 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1]
					a4d2 = names2[0] + ":" + names2[1]
				elif digits == 2:
					if len(names1) < 1 or len(names2) < 1:
						sys.exit("--digits 2 requires at least 2 digits resolution genotype!")
					a4d1 = names1[0]
					a4d2 = names2[0]

				if alleles[1] == "2":
					if a4d1.split('*')[0] in np:
						np[a4d1.split('*')[0]] += 1
					else:
						np[a4d1.split('*')[0]] = 1

					if a4d1 in case:
						case[a4d1] += 1
					else:
						case[a4d1] = 1

					if a4d1 != a4d2:
						if a4d2 in case:
							case[a4d2] += 1
						else:
							case[a4d2] = 1
						

				elif alleles[1] == "1":
					if a4d1.split('*')[0] in nc:
						nc[a4d1.split('*')[0]] += 1
					else:
						nc[a4d1.split('*')[0]] = 1

					if a4d1 in control:
						control[a4d1] += 1
					else:
						control[a4d1] = 1

					if a4d1 != a4d2:
						if a4d2 in control:
							control[a4d2] += 1
						else:
							control[a4d2] = 1	
	f.close()
	return case, control, np, nc

def recCount(infile, digits):
	'''
	count all alleles
	return the count of each allele in case and sontrol
	'''
	case = {}     # counts for each allele
	control = {}
	nc = {}
	np = {}
	f = open(infile)
	for line in f:
		line = line.rstrip()
		alleles = line.split()
		for i in range(2,len(alleles),2):
			j = i + 1
			if alleles[i] != 'NA' and alleles[j] != 'NA':
				names1 = alleles[i].split(":")
				names2 = alleles[j].split(":")
				if digits == 6:
					if len(names1) < 3 or len(names2) < 3:
						sys.exit("--digits 6 requires at least 6 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1] + ":" + names1[2]
					a4d2 = names2[0] + ":" + names2[1] + ":" + names2[2]
				if digits == 4:
					if len(names1) < 2 or len(names2) < 2:
						sys.exit("--digits 4 requires at least 4 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1]
					a4d2 = names2[0] + ":" + names2[1]
				elif digits == 2:
					if len(names1) < 1 or len(names2) < 1:
						sys.exit("--digits 2 requires at least 2 digits resolution genotype!")
					a4d1 = names1[0]
					a4d2 = names2[0]

				if alleles[1] == "2":
					if a4d1.split('*')[0] in np:
						np[a4d1.split('*')[0]] += 1
					else:
						np[a4d1.split('*')[0]] = 1

					if a4d1 == a4d2:
						if a4d2 in case:
							case[a4d2] += 1
						else:
							case[a4d2] = 1

				elif alleles[1] == "1":
					if a4d1.split('*')[0] in nc:
						nc[a4d1.split('*')[0]] += 1
					else:
						nc[a4d1.split('*')[0]] = 1

					if a4d1 == a4d2:
						if a4d2 in control:
							control[a4d2] += 1
						else:
							control[a4d2] = 1
	f.close()
	return case, control, np, nc

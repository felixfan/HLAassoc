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
		for allele in alleles[2:]:
			if allele != 'NA':
				names = allele.split(":")
				if digits == 4:
					a4d = names[0] + ":" + names[1]
				elif digits == 2:
					a4d = names[0]
				
				if alleles[1] == "2":
					if a4d.split("*")[0] in np:
						np[a4d.split("*")[0]] += 1
					else:
						np[a4d.split("*")[0]] = 1

					if a4d in case:
						case[a4d] += 1
					else:
						case[a4d] = 1
				elif alleles[1] == "1":
					if a4d.split("*")[0] in nc:
						nc[a4d.split("*")[0]] += 1
					else:
						nc[a4d.split("*")[0]] = 1

					if a4d in control:
						control[a4d] += 1
					else:
						control[a4d] = 1
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
				if digits == 4:
					a4d1 = names1[0] + ":" + names1[1]
					a4d2 = names2[0] + ":" + names2[1]
				elif digits == 2:
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
				if digits == 4:
					a4d1 = names1[0] + ":" + names1[1]
					a4d2 = names2[0] + ":" + names2[1]
				elif digits == 2:
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

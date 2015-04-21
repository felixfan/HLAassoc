#!/usr/bin/env python

import argparse
import HLAtest
import HLAregression
import pAdjust
import HLAperm
import HLAcount

parser = argparse.ArgumentParser(description='HLA Association Analysis', prog="HLAassoc.py")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.7.2')
parser.add_argument('-i', '--file', help='input file', required=True, type=str)
parser.add_argument('-d', '--digits', help='digits to test, default 4', default=4, type=int, choices=[2,4,6])
parser.add_argument('-m', '--model', help='genetic model, default allelic', default='allelic', type=str, choices=['allelic','dom','rec'])
parser.add_argument('-t', '--test', help='statistical test method, default chisq', default='chisq', type=str, choices=['chisq','fisher','logistic','linear','raw','score'])
parser.add_argument('-c', '--covar', help='covariants file', type=str)
parser.add_argument('-n', '--covarname', help='select a particular subset of covariates', type=str)
parser.add_argument('-f', '--freq', help='minimal frequency, default 0.05', default=0.05, type=float)
parser.add_argument('-a', '--adjust', help='p value correction, default FDR', default='FDR', type=str,choices=['FDR','FDR_BY','Bonferroni','Holm'])
parser.add_argument('-o', '--out', help='output file', default='hlaassoc.txt')
parser.add_argument('-V', '--print', help='print output to screen', type=str, default='False', choices=['False', 'True'])
parser.add_argument('-p', '--perm', help='number of permutation', type=int)
parser.add_argument('-s', '--seed', help='random seed', type=int)

args = vars(parser.parse_args())

INFILE = args['file']
OUTFILE = args['out']
FREQ = args['freq']
DIGIT = args['digits']
MODEL = args['model']
ADJUST = args['adjust']
TEST = args['test']
PRINT = (args['print'] == 'True')
SEED = args['seed']

COVFILE = False
COVNAME = False
if 'covar' in args:
	COVFILE = args['covar']
if 'covarname' in args:
	COVNAME = args['covarname']

PERM = 0
if 'perm' in args:
	PERM = args['perm']

#####################################################################
print "@-------------------------------------------------------------@"
print "|       HLAassoc       |     v1.7.2    |     15 Apr 2015      |"
print "|-------------------------------------------------------------|"
print "|  (C) 2015 Felix Yanhui Fan, GNU General Public License, v2  |"
print "|-------------------------------------------------------------|"
print "|    For documentation, citation & bug-report instructions:   |"
print "|        http://felixfan.github.io/HLAassoc                   |"
print "@-------------------------------------------------------------@"
print "\n\tOptions in effect:"
print "\t--file", INFILE
print "\t--digits", DIGIT
print "\t--test", TEST
if TEST == 'chisq' or TEST == 'fisher':
	print "\t--model", MODEL
elif TEST == 'logistic' or 'linear':
	if COVFILE:
		print "\t--covar", COVFILE
		if COVNAME:
			print "\t--covarname", COVNAME
print "\t--freq", FREQ
if TEST != 'raw' and TEST != 'score':
	print "\t--adjust", ADJUST
print "\t--print", PRINT
if PERM:
	print "\t--perm", PERM
	if SEED:
		print "\t--seed", SEED
print "\t--out", OUTFILE
print
#####################################################################
### permutation
if PERM:
	if TEST == 'chisq' or TEST == 'fisher':
		permp, permn, permna, permnl = HLAperm.chisqFisherPerm(INFILE, DIGIT, MODEL, TEST, PERM, SEED, FREQ)
	elif TEST == 'raw':
		permp, permn, permna, permnl = HLAperm.rawPerm(INFILE, DIGIT, PERM, SEED, FREQ)
	elif TEST == 'score':
		permp, permn, permna, permnl = HLAperm.scorePerm(INFILE, DIGIT, PERM, SEED, FREQ)
	elif TEST == 'logistic' or TEST == 'linear':
		if COVFILE:
			permp, permn, permna, permnl = HLAperm.regressionCovPerm(INFILE, DIGIT, TEST, PERM, COVFILE, COVNAME, SEED,FREQ)
		else:
			permp, permn, permna, permnl = HLAperm.regressionPerm(INFILE, DIGIT, TEST, PERM,SEED, FREQ)
### output header
f = open(OUTFILE,"w")
if TEST == 'chisq':
	header = ["Allele","A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","Freq","Chisq","DF","P_chisq","OR","L95","U95","P_adj"]
elif TEST == 'fisher':
	header = ["Allele","A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","Freq","P_Fisher","OR","L95","U95","P_adj"]
elif TEST == 'logistic':
	header = ["Allele","A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","Freq","P_logit","OR","L95","U95","P_adj"]
elif TEST == 'linear':
	header = ["Allele","Freq","P_linear","beta","L95","U95","P_adj"]
elif TEST == 'raw':
	header = ["Gene","Chisq","DF","P_chisq"]
elif TEST == 'score':
	header = ["Gene","U"]
if PERM:
	header.append('P_perm')
for h in header:
	if h == 'Allele' or h == 'Gene':
		if PRINT:
			print "%12s" % h,
		f.write("%12s" % h,)
	else:
		if PRINT:
			print "%8s" % h,
		f.write("%8s" % h,)
if PRINT:		
	print
f.write('\n')
############################################################
### write perm output
if PERM:
	fp = open(OUTFILE + '.perm', 'w')
	fp.write("%12s%12s%12s%12s%12s%12s\n" % ('ID', 'L', 'S', 'NA', 'Perm', 'Perm_P'))
	for a in sorted(permp.keys()):
		fp.write("%12s" % a,)
		fp.write("%12s" % permn[a],)
		fp.write("%12s" % permnl[a],)
		fp.write("%12s" % permna[a],)
		fp.write("%12s" % PERM,)
		if permp[a] != 'NA':
			fp.write("%12s\n" % str(round(permp[a],6)))
		else:
			fp.write("%12s\n" % 'NA')
	fp.close()
############################################################
### run association analysis
rs = []
rs2 ={}
if TEST == 'chisq' or TEST == 'fisher':
	rs = HLAtest.runAssoc(INFILE, DIGIT, FREQ, MODEL, TEST)
elif TEST == 'raw':
	caseAlleles, ctrlAlleles, np, nc = HLAcount.allelicCount(INFILE, DIGIT)
	rs2 = HLAtest.assocRaw(caseAlleles, ctrlAlleles, np, nc, FREQ, TEST)
elif TEST == 'score':
	caseAlleles, ctrlAlleles, np, nc = HLAcount.allelicCount(INFILE, DIGIT)
	rs2 = HLAtest.assocScoreU(caseAlleles, ctrlAlleles, np, nc, FREQ, TEST)
elif TEST == 'logistic' or TEST == 'linear':
	if COVFILE:
		if TEST == 'logistic':
			ans = HLAregression.logisticCov(INFILE, DIGIT, FREQ, TEST, COVFILE, COVNAME)
		else:
			ans = HLAregression.linearCov(INFILE, DIGIT, FREQ, TEST, COVFILE, COVNAME)
	else:
		if TEST == 'logistic':
			ans = HLAregression.regressionLogistic(INFILE, DIGIT, FREQ, TEST)
		else:
			ans = HLAregression.regressionLinear(INFILE, DIGIT, FREQ, TEST)
	alleles =sorted(ans.keys())
	genes = {}
	for allele in alleles:
		genes[allele.split('*')[0]] = 1
	sortedGene = sorted(genes.keys())
	for g in sortedGene:
		tmp = {}
		for allele in alleles:
			if allele.split('*')[0] == g:
				tmp[allele] = ans[allele]
		rs.append(tmp)

############################################################
### adjust p values and write output
### position of p value
if TEST == 'chisq':
	pp = 10
elif TEST == 'linear':
	pp = 2
else:                   # logistic and fisher
	pp = 8
if TEST != 'raw' and TEST != 'score':
	for r in rs:                        ### GENE BY GENE
		keys = r.keys()
		keys = sorted(keys)
		ps = []
		### get p values for a group
		for key in keys:
			words = r[key]
			if words[pp] != 'NA':
				ps.append(words[pp])
		### adjust p
		cp = pAdjust.adjustP(ps,ADJUST)
		### output
		for key in keys:
				words = r[key]
				for word in words: # test results
					if word == words[0]:
						if PRINT:
							print "%12s" % word,
						f.write("%12s" % word,)
					elif word == words[pp]:
						if PRINT:
							if word == 'NA':
								print "%8s" % word,
							elif word > 0.001:
								print "%8.4f" % word,
							else:
								print "%8.2e" % word,
						if word == 'NA':
							f.write("%8s" % word,)
						elif word > 0.001:
							f.write("%8.4f" % word,)
						else:
							f.write("%8.2e" % word,)
					else:
						if PRINT:
							if isinstance(word, int):
								print "%8d" % word,
							elif isinstance(word, float):
								print "%8.4f" % word,
							elif isinstance(word, str):
								print "%8s" % word,
						if isinstance(word, int):
							f.write("%8d" % word,)
						elif isinstance(word, float):
							f.write("%8.4f" % word,)
						elif isinstance(word, str):
							f.write("%8s" % word,)
				### adjust
				if words[pp] != 'NA':
					tmp = cp.pop(0)
					if PRINT:
						if tmp > 0.001:
							print "%8.4f" % tmp,
						else:
							print "%8.2e" % tmp,
					if tmp > 0.001:
						f.write('%8.4f' % tmp,)
					else:
						f.write('%8.2e' % tmp,)
				else:
					if PRINT:
						print "%8s" % 'NA',
					f.write('%8s' % 'NA',)
				### perm
				if PERM:
					if key in permp:
						if PRINT:
							if permp[key] != 'NA':
								if  permp[key] > 0.001:
									print "%8.4f" % permp[key],
								else:
									print "%8.2e" % permp[key],
							else:
								print "%8s" % 'NA',
						if permp[key] != 'NA':
							if  permp[key] > 0.001:
								f.write('%8.4f' % permp[key],)
							else:
								f.write('%8.2e' % permp[key],)
						else:
							f.write('%8s' % 'NA',)
					else:
						if PRINT:
							print "%8s" % 'NA',
						f.write('%8s' % 'NA',)
				# new line
				if PRINT:
					print
				f.write('\n')
	f.close()
elif TEST == 'raw': # raw test
	keys = sorted(rs2.keys())
	for key in keys:
		if PRINT:
			print "%12s" % key,
		f.write('%12s' % key)

		fs = rs2[key]
		for ff in fs:
			if PRINT:
				if ff == fs[2]:           # p-value column
					if ff == 'NA':
						print "%8s" % ff,
					elif ff > 0.001:
						print "%8.4f" % ff,
					else:
						print "%8.2e" % ff,
				else:
					if isinstance(ff, int):
						print "%8d" % ff,
					elif isinstance(ff, float):
						print "%8.4f" % ff,
					elif isinstance(ff, str):
						print "%8s" % word,
			if ff == fs[2]:           # p-value column
				if ff == 'NA':
					f.write("%8s" % ff,)
				elif ff > 0.001:
					f.write("%8.4f" % ff,)
				else:
					f.write("%8.2e" % ff,)
			else:
				if isinstance(ff, int):
					f.write("%8d" % ff,)
				elif isinstance(ff, float):
					f.write("%8.4f" % ff,)
				elif isinstance(ff, str):
					f.write("%8s" % word,)
		### perm
		if PERM:
			if key in permp:
				if PRINT:
					if permp[key] != 'NA':
						if  permp[key] > 0.001:
							print "%8.4f" % permp[key],
						else:
							print "%8.2e" % permp[key],
					else:
						print "%8s" % 'NA',
				if permp[key] != 'NA':
					if  permp[key] > 0.001:
						f.write('%8.4f' % permp[key],)
					else:
						f.write('%8.2e' % permp[key],)
				else:
					f.write('%8s' % 'NA',)
			else:
				if PRINT:
					print "%8s" % 'NA',
				f.write('%8s' % 'NA',)
		# new line
		if PRINT:
			print
		f.write('\n')
	f.close()
else: # score test U
	keys = sorted(rs2.keys())
	for key in keys:
		if PRINT:
			print "%12s" % key,
			print "%8.4f" % rs2[key],
		f.write('%12s' % key,)
		f.write('%8.4f' % rs2[key],)
		### perm
		if PERM:
			if key in permp:
				if PRINT:
					if permp[key] != 'NA':
						if  permp[key] > 0.001:
							print "%8.4f" % permp[key],
						else:
							print "%8.2e" % permp[key],
					else:
						print "%8s" % 'NA',
				if permp[key] != 'NA':
					if  permp[key] > 0.001:
						f.write('%8.4f' % permp[key],)
					else:
						f.write('%8.2e' % permp[key],)
				else:
					f.write('%8s' % 'NA',)
			else:
				if PRINT:
					print "%8s" % 'NA',
				f.write('%8s' % 'NA',)
		# new line
		if PRINT:
			print
		f.write('\n')
	f.close()
##############################################################
import time
print "Finished at ",
print time.strftime("%H:%M:%S %d %b %Y")
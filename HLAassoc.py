#!/usr/bin/env python

import argparse
import HLAtest
import HLAregression
import pAdjust
import HLAperm
import HLAcount

parser = argparse.ArgumentParser(description='HLA Association Analysis', prog="HLAassoc.py")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.6')
parser.add_argument('-i', '--file', help='input file', required=True, type=str)
parser.add_argument('-d', '--digits', help='digits to test, default 4', default=4, type=int, choices=[2,4,6])
parser.add_argument('-m', '--model', help='genetic model, default allelic', default='allelic', type=str, choices=['allelic','dom','rec'])
parser.add_argument('-t', '--test', help='statistical test method, default chisq', default='chisq', type=str, choices=['chisq','fisher','logistic','linear','raw','score'])
parser.add_argument('-c', '--covar', help='covariants file', type=str)
parser.add_argument('-n', '--covarname', help='select a particular subset of covariates', type=str)
parser.add_argument('-f', '--freq', help='minimal frequency, default 0.05', default=0.05, type=float)
parser.add_argument('-a', '--adjust', help='p value correction, default FDR', default='FDR', type=str,choices=['FDR','Bonferroni','Holm'])
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
print "|       HLAassoc       |     v 1.7     |     15 Apr 2015      |"
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
	header = ["Allele","A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","Freq","P_logistic","OR","L95","U95","P_adj"]
elif TEST == 'linear':
	header = ["Allele","Freq","P_linear","beta","L95","U95","P_adj"]
elif TEST == 'raw':
	header = ["Gene","Chisq","DF","P_chisq"]
elif TEST == 'score':
	header = ["Gene","U"]
if PERM:
	header.append('P_perm')
for h in header:
	if PRINT:
		print "%12s" % h,
	f.write("%12s" % h,)
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
else:
	pp = 8
if TEST != 'raw' and TEST != 'score':
	for r in rs:                        ### GENE BY GENE
		keys = r.keys()
		keys = sorted(keys)
		ps = []
		### get p values for a group
		for key in keys:
			words = r[key].split()
			if words[pp] != 'NA':
				ps.append(float(words[pp]))
		### adjust p
		cp = pAdjust.adjustP(ps,ADJUST)
		### output
		for key in keys:
				words = r[key].split()
				for word in words: # test results
					if PRINT:
						print "%12s" % word,
					f.write("%12s" % word,)
				### adjust
				if words[pp] != 'NA':
					tmp = cp.pop(0)
					if PRINT:
						print "%12s" % str(round(tmp,6)),
					f.write('%12s' % str(round(tmp,6)))
				else:
					if PRINT:
						print "%12s" % ('NA'),
					f.write('%12s' % ('NA'))
				### perm
				if PERM:
					if key in permp:
						if PRINT:
							if permp[key] != 'NA':
								print "%12s" % str(round(permp[key],6)),
							else:
								print "%12s" % 'NA',
						if permp[key] != 'NA':
							f.write('%12s' % str(round(permp[key],6)))
						else:
							f.write('%12s' % 'NA')
					else:
						if PRINT:
							print "%12s" % 'NA',
						f.write('%12s' % 'NA')
				# new line
				if PRINT:
					print
				f.write('\n')
	f.close()
elif TEST == 'raw': # raw test
	keys = sorted(rs2.keys())
	for k in keys:
		if PRINT:
			print "%12s" % k,
		f.write('%12s' % k)

		fs = rs2[k].split()
		for ff in fs:
			if PRINT:
				print "%12s" % ff,
			f.write('%12s' % ff)
		### perm
		if PERM:
			if k in permp:
				if PRINT:
					if permp[k] != 'NA':
						print "%12s" % str(round(permp[k],6)),
					else:
						print "%12s" % 'NA',
				if permp[k] != 'NA':
					f.write('%12s' % str(round(permp[k],6)))
				else:
					f.write('%12s' % 'NA')
			else:
				if PRINT:
					print "%12s" % 'NA',
				f.write('%12s' % 'NA')

		if PRINT:
			print
		f.write('\n')
else: # score test U
	keys = sorted(rs2.keys())
	for k in keys:
		if PRINT:
			print "%12s" % k,
			print "%12s" % rs2[k],
		f.write('%12s' % k)
		f.write('%12s' % rs2[k])
		### perm
		if PERM:
			if k in permp:
				if PRINT:
					if permp[k] != 'NA':
						print "%12s" % str(round(permp[k],6)),
					else:
						print "%12s" % 'NA',
				if permp[k] != 'NA':
					f.write('%12s' % str(round(permp[k],6)))
				else:
					f.write('%12s' % 'NA')
			else:
				if PRINT:
					print "%12s" % 'NA',
				f.write('%12s' % 'NA')

		if PRINT:
			print
		f.write('\n')
##############################################################
import time
print "Finished at ",
print time.strftime("%H:%M:%S %d %b %Y")
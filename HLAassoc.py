#!/usr/bin/env python

import argparse
import HLAtest
import HLAregression
import pAdjust

parser = argparse.ArgumentParser(description='HLA Association Analysis', prog="HLAassoc.py")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.1')
parser.add_argument('-i', '--file', help='input file', required=True, type=str)
parser.add_argument('-d', '--digits', help='digits to test, default 4', default=4, type=int, choices=[2,4])
parser.add_argument('-m', '--model', help='genetic model, default allelic', default='allelic', type=str, choices=['allelic','dom','rec'])
parser.add_argument('-t', '--test', help='statistical test method, default chisq', default='chisq', type=str, choices=['chisq','fisher','logistic','linear'])
parser.add_argument('-c', '--covar', help='covariants file', type=str)
parser.add_argument('-n', '--covarname', help='select a particular subset of covariates', type=str)
parser.add_argument('-f', '--freq', help='minimal frequency, default 0.05', default=0.05, type=float)
parser.add_argument('-a', '--adjust', help='p value correction, default FDR', default='FDR', type=str,choices=['FDR','Bonferroni','Holm'])
parser.add_argument('-o', '--out', help='output file', default='hlaassoc.txt')
parser.add_argument('-p', '--print', help='print output to screen', type=bool, default=False, choices=[False, True])

args = vars(parser.parse_args())

INFILE = args['file']
OUTFILE = args['out']
FREQ = args['freq']
DIGIT = args['digits']
MODEL = args['model']
ADJUST = args['adjust']
TEST = args['test']
PRINT =args['print']

COVFILE = False
COVNAME = False
if 'covar' in args:
	COVFILE = args['covar']
if 'covarname' in args:
	COVNAME = args['covarname']

#####################################################################
print "@-------------------------------------------------------------@"
print "|       HLAassoc       |     v 1.1     |     12 Feb 2015      |"
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
else:
	if COVFILE:
		print "\t--covar", COVFILE
		if COVNAME:
			print "\t--covarname", COVNAME
print "\t--freq", FREQ
print "\t--adjust", ADJUST
print "\t--print", PRINT
print "\t--out", OUTFILE
print
#####################################################################
### output header
f = open(OUTFILE,"w")
if TEST == 'chisq':
	header = ("ID","A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","Chisq","DF","P_chisq","OR","L95","U95","P_adj")
elif TEST == 'fisher':
	header = ("ID","A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","P_Fisher","OR","L95","U95","P_adj")
elif TEST == 'logistic':
	header = ("ID","A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","P_logistic","OR","L95","U95","P_adj")
elif TEST == 'linear':
	header = ("ID","A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","P_linear","beta","L95","U95","P_adj")
for h in header:
	if PRINT:
		print "%12s" % h,
	f.write("%12s" % h,)
if PRINT:		
	print
f.write('\n')

############################################################
### run association analysis
rs = []
if TEST == 'chisq' or TEST == 'fisher':
	rs = HLAtest.runAssoc(INFILE, DIGIT, FREQ, MODEL, TEST)
elif TEST == 'logistic' or TEST == 'linear':
	if COVFILE:
		ans = HLAregression.regressionCov(INFILE, DIGIT, FREQ, TEST, COVFILE, COVNAME)
	else:
		ans = HLAregression.regression(INFILE, DIGIT, FREQ, TEST)
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
	pp = 9
else:
	pp = 7

for r in rs:                        ### GENE BY GENE
	keys = r.keys()
	keys = sorted(keys)
	ps = []
	### get p values for a group
	for key in keys:
		words = r[key].split()
		ps.append(float(words[pp]))
	### adjust p
	cp = pAdjust.adjustP(ps,ADJUST)
	### output
	for key in keys:
			words = r[key].split()
			for word in words:
				if PRINT:
					print "%12s" % word,
				f.write("%12s" % word,)
			tmp = cp.pop(0)
			if PRINT:
				print "%12s" % str(round(tmp,6))
			f.write('%12s\n' % str(round(tmp,6)))
f.close()

##############################################################
import time
print "Finished at ",
print time.strftime("%H:%M:%S %d %b %Y")
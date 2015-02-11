import argparse
import HLAtest
import pAdjust

parser = argparse.ArgumentParser(description='HLA Association Analysis', prog="HLAassoc.py")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
# parser.add_argument('-i','--file', help='input file', required=True, type=argparse.FileType('r'))
parser.add_argument('-i', '--file', help='input file', required=True, type=str)
parser.add_argument('-d', '--digits', help='digits to test, default 4', default=4, type=int, choices=[2,4])
parser.add_argument('-m', '--model', help='genetic model, default allelic', default='allelic', type=str, choices=['allelic','dom','rec'])
parser.add_argument('-t', '--test', help='statistical test method, default chisq', default='chisq', type=str, choices=['chisq','fisher'])
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
#####################################################################
print "@-------------------------------------------------------------@"
print "|       HLAassoc       |     v 1.0     |     23 Jan 2015      |"
print "|-------------------------------------------------------------|"
print "|  (C) 2015 Felix Yanhui Fan, GNU General Public License, v2  |"
print "|-------------------------------------------------------------|"
print "|    For documentation, citation & bug-report instructions:   |"
print "|        http://felixfan.github.io/HLAassoc                   |"
print "@-------------------------------------------------------------@"
print "\n\tOptions in effect:"
print "\t--file", INFILE
print "\t--digits", DIGIT
print "\t--model", MODEL
print "\t--test", TEST
print "\t--freq", FREQ
print "\t--adjust", ADJUST
print "\t--print", PRINT
print "\t--out", OUTFILE
print
#####################################################################
### call functions and print out results
f = open(OUTFILE,"w")
if TEST == 'chisq':
	if ADJUST:
		header = ("ID","A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","Chisq","DF","P_chisq","OR","L95","U95","P_adj")
	else:
		header = ("ID","A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","Chisq","DF","P_chisq","OR","L95","U95")
elif TEST == 'fisher':
	if ADJUST:
		header = ("ID","A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","P_Fisher","OR","L95","U95","P_adjusted")
	else:
		header = ("ID","A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","P_Fisher","OR","L95","U95")
for h in header:
	if PRINT:
		print "%10s" % h,
	f.write("%10s" % h,)
if PRINT:		
	print
f.write('\n')

### run association analysis
rs = HLAtest.runAssoc(INFILE, DIGIT, FREQ, MODEL, TEST)

### adjust p values and write output
### position of p value
if TEST == 'chisq':
	pp = 9
elif TEST == 'fisher':
	pp = 7

for r in rs:
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
					print "%10s" % word,
				f.write("%10s" % word,)
			tmp = cp.pop(0)
			if PRINT:
				print "%10s" % str(round(tmp,6))
			f.write('%10s\n' % str(round(tmp,6)))
	# print '-------'
f.close()

##############################################################
import time
print "Finished at ",
print time.strftime("%H:%M:%S %d %b %Y")
'''
This file reads .gaf files and finds out the annotation statistics
[1]--> file path
Urminder Singh
'''
from __future__ import division
import sys
import Bio.UniProt.GOA
#read .gaf file with gafiterator
prots_exp=[]
evd_exp=[]
prots_all=[]
exp_ids=(["EXP","IDA","IPI","IMP","IGI","IEP"])
with open(sys.argv[1], 'r') as handle:
	for rec in Bio.UniProt.GOA.gafiterator(handle):
		#make two different lists for experimentally verified and other type of prots
		prots_all.append(rec["DB_Object_ID"])
		if rec["Evidence"] in exp_ids:
			prots_exp.append(rec["DB_Object_ID"])
			evd_exp.append(rec["Evidence"])
#remove redundant entries
totalprots=len(set(prots_all))
num_exp_prots=len(set(prots_exp))
num_noexp_prots=totalprots-num_exp_prots
print 'Results for input file:',sys.argv[1]
print 'Total proteins:',totalprots
print 'Total proteins with experimental evidence:', num_exp_prots,'fraction:',(num_exp_prots/totalprots*100)
print 'Total proteins without experimental evidence:', num_noexp_prots,' fraction:',(num_noexp_prots/totalprots*100)


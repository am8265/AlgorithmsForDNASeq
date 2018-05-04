#!/usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd


def naive(p, t):
	occurences = []
	for i in range(len(t) - len(p) + 1):  # loop over alignments
		match = True
	  	for j in range(len(p)):  # loop over characters
	  		if t[i + j] != p[j]:  # compare chracters
	  			match = False
	  			break
	  	if match:
	  		occurences.append(i)  # all chars matched; record
	return occurences


def reverseComplement(s):
	complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
	t = ''
	for base in s:
		t = complement[base] + t
	return t

def readGenome(filename):
	genome = ''
	with open(filename, 'r') as f:
		for line in f:
			# ignore header line with genome information
			if not line[0] == '>':
				genome += line.rstrip()
	return genome


def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def naive_with_rc(p, t):
	rev_comp_p = reverseComplement(p)
 	occurence_rev_comp_p = naive(rev_comp_p, t)
 	return occurence_rev_comp_p

def strandAwareNaive(p, t):
	occurence = naive(p, t)
 	occurence.extend(naive_with_rc(p, t))
 	occurence = list(set(occurence))
 	return occurence


def naive_2mm(p, t):
	occurences = []
	for i in range(len(t) - len(p) + 1):  # loop over alignments
		match = True
		c = 0  # c to keep track of number of mismatches 
	  	for j in range(len(p)):  # loop over characters
	  		if t[i + j] != p[j]:  # compare chracters
	  			c += 1
	  			if c > 2:
	  				match = False
	  				break
	  	if match:
	  		occurences.append(i)  # all chars matched; record
	return occurences



def QtoPhred33(Q):
	return chr(Q + 33)


def phred33ToQ(qual):
	return ord(qual) - 33


def cycErrors(filename):
	sequences, qualities = readFastq(filename)
	qual_transfrm_final = []
	for qual_list in qualities:
		qual_transfrm = [phred33ToQ(qual) for qual in qual_list]
		qual_transfrm_final.append(qual_transfrm)
	
	return sequences, qual_transfrm_final, np.array(qual_transfrm_final), pd.DataFrame(qual_transfrm_final)

test = cycErrors("ERR037900_1.first1000.fastq")

df = test[3]
mean_sr = df.mean(axis = 0)

print mean_sr

print mean_sr.idxmin(axis = 0)
#mean_sr.idxmin
#print mean_sr.idxmin



#t = readGenome(sys.argv[1])
#p = 'AGGAGGTT'

#print sorted(list(set(naive_2mm(p, t))))

#print len(sorted(list(set(naive(p, t)))))


#print sorted(strandAwareNaive(p, t))
#print len(strandAwareNaive(p, t))




# Written by Ehsan Motazedi, Wageningen UR, 06-02-2017.
# Last Updated: 05-07-2017

import sys
from logprob import diffVec, Hamming
from haplotypes import Haplotypes
from libc.math cimport log, exp
from reads import Read

cdef extern from "math.h":
	float INFINITY

cpdef double loge(double x):
	""" Return log x in base e."""
	try:
		return log(x)
	except ValueError as e:
		garbage = sys.stderr.write("WARNING: {0}\n".format(e.args[0]))
		return -1*INFINITY

def veclog(vec):
	""" Logarithm function defined on vector space."""
	try:
		return map(loge, vec)
	except TypeError:
		return loge(vec)

def GetlogProb(r, Vset, eps, qual = None, counts = False, min_r_length=2):
	""" Return the log probability of a read conditional on a presumed haplotype and 
	variant error rate, i.e. P(r|Vset, eps) (Berger et al. 2014, p. 4). If counts if True,
	assign the read to the most compatible homologue (count 1 for that homologue and zero of 
	the others) to later calculate the number of reads compatible with each homologue for a set 
	of reads by the upstream functions."""
	ar = r.GetalleleSet(Vset.GetStart(), Vset.GetStop())
	if len([_allele for _allele in ar if _allele not in set(['-','.'])])<min_r_length:
		if counts:
			return [0, [0 for _h in range(0, len(Vset.GetVS()))]]
		return 0 # Corresponding to P[r|Vset, eps] = 1
	sr='-'*(r.GetBegin()-Vset.GetStart())+ar+'-'*(Vset.GetStop()-r.GetEnd()) # Make the lengths equal to calculate the distance
	#warn, veclogprobs = False, []
	warn, probs = False, []
	if qual:
		try:
			er = []
			escores = []
			for _pos in r.GetPos():
				escores.append(10**(-float(qual[_pos])/10))
			_n = 0
			for _x in ar:
				if _x=='-':
					er.append('-')
				else:
					er.append(escores[_n])
					_n+=1
			er1 = []
			er2 = []
			for _x in range(0, r.GetBegin()-Vset.GetStart()):
				er1.append('-')
			for _x in range(0, Vset.GetStop()-r.GetEnd()):
				er2.append('-')
			er = er1 + er
			er = er + er2
		except ValueError as e:
			e.args=(e.args[0]+"\nInvalid quality scores detected in the alignment file!",)+e.args[1:]
			raise
		except KeyError as e:
			warn = True
			qual = None
	#garbage=sys.stderr.write("******\nER:\n "+str(er)+"\n")
	#garbage=sys.stderr.write("\nSR:\n "+str(sr)+"\n")
	if not qual:
		for _homolo in Vset.GetVS(): # _homolo corresponds to v in Vset in Berger et al. 2014 p. 4              
			D, A = Hamming(sr, _homolo) # D corresponds to D(r,v) and A correspond to A(r,v) in Berger et al. 2014 p. 4
			#veclogprobs.append(A*veclog((1-eps)/(1-2./3*eps))+D*veclog(eps/3.)) # in case of error, any of the three other bases could have been called. So the probability for a specific erroneous base-call is eps./3. As the probability of no error is 1-eps, we have to normalize by 1-eps+eps/3 to get probs for dichotomy no-error/specific erroneous base-call. 
			probs.append(((1-eps)/(1-2./3*eps))**A*(eps/3.)**D) # in case of error, any of the three other bases could have been called. So the probability for a specific erroneous base-call is eps./3. As the probability of no error is 1-eps, we have to normalize by 1-eps+eps/3 to get probs for dichotomy no-error/specific erroneous base-call. 
	else: # Use eps(s) as a function of SNP position in the read instead of a fixed eps.
		for _homolo in Vset.GetVS():
			D, A = diffVec(sr, _homolo)
			#garbage=sys.stderr.write("\nD:\n "+str(D)+"\n")
			#garbage=sys.stderr.write("\nA:\n "+str(A)+"\n")
			#_veclogprob = 0
			_prob = 1
			for _s in A:
				#_veclogprob+=loge((1-er[_s])/(1-2./3*er[_s]))
				_prob*= ((1-er[_s])/(1-2./3*er[_s]))
			for _s in D:
				#_veclogprob+=loge(er[_s]/3.)
				_prob *= (er[_s]/3.)
			#veclogprobs.append(_veclogprob)
			probs.append(_prob)
		#garbage=sys.stderr.write("******\n")
	if warn:
		garbage = sys.stderr.write("WARNING: Quality scores specified for none or only some positions of the read! The fixed given or deafult error rate was used for all the positions!\n")
	unique_homolos = set(Vset.GetVS())
	dict_homolos = {_dh:1./len([_h for _h in Vset.GetVS() if _h == _dh]) for _dh in unique_homolos} # the counts are set equally to 1/j for j similar homologues
	h_counts=[0 for _h in Vset.GetVS()]
	for _k, _h in enumerate(Vset.GetVS()):
		#if _h==Vset.GetVS()[veclogprobs.index(min(veclogprobs))]:
		if _h==Vset.GetVS()[probs.index(min(probs))]:
			h_counts[_k]=dict_homolos[_h]
	if counts:
		#print([loge(sum(probs)/float(len(Vset.GetVS()))), h_counts])# corresponds to log(P[r|Vset, eps]) = log(1/k * sum(P(r|v in Vset))) in Berger et al. 2014 p. 4, assign the read to the most compatible homologue if counts of the reads assigned to each homologue are needed.
		return [loge(sum(probs)/float(len(Vset.GetVS()))), h_counts]  # corresponds to log(P[r|Vset, eps]) = log(1/k * sum(P(r|v in Vset))) in Berger et al. 2014 p. 4, assign the read to the most compatible homologue if counts of the reads assigned to each homologue are needed.
		#return [loge(sum(map(exp,veclogprobs)))-loge(len(Vset.GetVS())), h_counts]  # corresponds to log(P[r|Vset, eps]) = log(1/k * sum(P(r|v in Vset))) in Berger et al. 2014 p. 4, assign the read to the most compatible homologue if counts of the reads assigned to each homologue are needed.
	return loge(sum(probs)/float(len(Vset.GetVS()))) # corresponds to log(P[r|Vset, eps]) = log(1/k * sum(P(r|v in Vset))) in Berger et al. 2014 p. 4
	#return veclog(sum(map(exp,veclogprobs)))-veclog(len(Vset.GetVS())) # corresponds to veclog(P[r|Vset, eps]) = veclog(1/k * sum(P(r|v in Vset))) in Berger et al. 2014 p. 4


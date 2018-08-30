# Writen by Ehsan Motazedi, Wageningen UR, 09-11-2016.
# Last Updated: 14-10-2017

import itertools
import multiprocessing as mpthread
import numpy.random as nprnd
import random as rnd
import sys
import threading
from collections import Counter
from genotypes import Genotype, getAllelesPop 
#from hapcompare import getVER
from haplotypes import Haplotypes, getMEC
from logprob import diffVec, Hamming
from logprob2 import veclog as log, loge, GetlogProb
from math import sqrt, exp, factorial
from numpy import array, divide, isnan, vectorize, delete as npdel, argwhere as npwhere
from reads import remove_last_SNP, SubReads

#Global variableis for the module to store all parental transmissions to the child
transmit_parent_m = None # possible transmissions of ploidy_m//2 alleles from mother
transmit_parent_f = None # possible transmissions of ploidy_f//2 alleles from father
child_haploid = None # possible haploid genomes of the child 
transmit_parent = None # possible transmissions of alleles to the child from the parents
MAXCORES = mpthread.cpu_count()-1 # Max number of available cores
NCORES = 6 # desired number of cores to use
AllGenos_m = None # All possible genoptypes for mother at each position
AllGenos_f = None # All possible genoptypes for father at each position
Global_Sema = None # Global Semaphore not initialized yet 
 
def GetSemaLock(useglobal=True):
	""" Return the semaphore and lock objects needed for multi threading."""
	global NCORES
	global MAXCORES
	global Global_Sema
	if min(NCORES, MAXCORES)>1:
		if useglobal: # use a global semaphore
			if Global_Sema is None:
				Global_Sema = mpthread.BoundedSemaphore(min(NCORES, MAXCORES))
		else: # use a local semaphore at each call to BranchPop
			sema = mpthread.BoundedSemaphore(min(NCORES, MAXCORES))
		lock = mpthread.Lock()
		Parallel = True
	else:   # If parallel execution is not possible (as just one core could be used), concurrent execution is performed using threads.
		if useglobal: 
			if Global_Sema is None:
				Global_Sema = threading.BoundedSemaphore(NCORES) 
		else:
			sema = threading.BoundedSemaphore(NCORES)
		lock = threading.Lock()
		Parallel = False
	if useglobal:
		return Global_Sema, lock , Parallel 
	else:
		return sema, lock , Parallel 

def thread_func(sema, lock, q, func, *args):
	""" thread to call func with args and put the return value in q. """
	_locked = False
	try:
		b = func(*args)
		_locked = lock.acquire()
		if func==Branch:
			for _x in b:
				q.append(_x)
		else:
			for _i, _b in enumerate(b):
				for _x in _b:
					q[_i].append(_x)
	except:
		raise
	finally:
		if _locked:
			lock.release()	
		sema.release()

def SetGenos(vcffile, SampleNames, ploidy_m, ploidy_f):
	""" Set AllGenos_m, AllGenos_f, AllGenos_c from the vcffile. Supposed to be called by main.py only once."""
	global AllGenos_m
	global AllGenos_f
	AllGenos = getAllelesPop(vcffile, SampleNames)
	AllGenos_m = {_Alleles[0]: list(itertools.product(*(_Alleles[2].values() for _h in range(0, ploidy_m)))) for _Alleles in AllGenos} 
	AllGenos_f = {_Alleles[0]: list(itertools.product(*(_Alleles[2].values() for _h in range(0, ploidy_f)))) for _Alleles in AllGenos} 
	for AllGenos in (AllGenos_m, AllGenos_f): # Throw away duplicate genotypes at each position from each parent's possibility list
		for _key in AllGenos.keys():
			_value = AllGenos[_key]
			_new_value = []
			_current = 0
			while _current<len(_value):
				if Counter(_value[_current]) not in [Counter(_x) for _x in _new_value]:
					_new_value.append(_value[_current])
				_current+=1
			AllGenos[_key]=_new_value
	
def reduce(function, iterable, initializer=None):
	it = iter(iterable)
	if initializer is None:
		try:
			initializer = next(it)
		except StopIteration:
			raise TypeError('reduce() of empty sequence with no initial value')
	accum_value = initializer
	for x in it:
		accum_value = function(accum_value, x)
	return accum_value
 
def remove_last_SNP_vec(SNPvec):
	return map(remove_last_SNP, SNPvec) 

def GetLogProbH(H):
	""" Determine the number of occurences of each unique homologue in H, and hence calculate P[H|eps] to be 
	used in calculating RLs, according to Berger et al. 2014 p. 4. Return log(P[H|eps]) in base 2."""
	M=[]
	Vset = []
	for _v in H.GetVS():
		if _v not in Vset:
			Vset.append(_v)
			M.append(1)
		else:
			M[Vset.index(_v)]+=1
	log_denom = 0
	k = sum(M) # the ploidy level
	for _M in M:
		log_denom += loge(factorial(_M))
	return loge(factorial(k))-log_denom-H.GetlogCS()

class BlockException(Exception):
	def __init__(self, value):
		super(BlockException, self).__init__(value)
		self.args = (value,)
	def __str__(self):
		return "{}".format(':'.join(str(_arg) for _arg in self.args))
	def __repr__(self):
		return self.args

def Branch(H, G, SReads, rho, error, qscores = None, usecounts=True, samplename=None, geno=True):
        """ Branch current haplotype H at position G.s, using G, the semi-reads for G.s and the thresholds rho (hard) 
	and kappa (soft). Variant error rate, i.e. error, is passed to calculate the branching probabilities. (Berger et al. 2014 p. 5)."""
	if set(G.GetGenes())==set(['.']) and geno: # If the genotype is missing, extension will be skipped!
		garbage = sys.stderr.write('WARNING: {1:s}\'s genotype is missing at position {2:d}! Phasing extension will be escaped at s={0:d} for {1:s}!\n'.format(G.GetS()+1, samplename if samplename else '', G.GetPos()))
		return [H+Haplotypes(G.GetS(), G.GetS(), 0, 0, None, None, *['-' for _homologue in H.GetVS()])]  # skip extension if no extension has been possible
	if len(set(G.GetGenes()))<2 and geno: # If the genotype is homozygous, the phasing will be trivially determined
		garbage = sys.stderr.write('WARNING: {1:s}\'s genotype is homozygous at position {2:d}! Its phasing extension will be trivial at s={0:d}!\n'.format(G.GetS()+1, samplename if samplename else '', G.GetPos()))
		return  [H+Haplotypes(G.GetS(), G.GetS(), 0, 0, None, None, *G.GetGenes())]
	if all(r.isNULL() for r in SReads):
		#raise BlockException('No semi-reads exist at SNP position {0:d}{1:s}!\n'.format(G.GetS()+1, " for "+samplename if samplename else ''))
		garbage = sys.stderr.write('WARNING: No semi-reads exist at SNP position {2:d}, s={0:d}{1:s}!\n'.format(G.GetS()+1, " for "+samplename if samplename else '', G.GetPos()))
	myrho = loge(rho) # change rho to log_2 scale
        ProbH = H.GetRL()
        count_threshold = int((1./len(H.GetVS())*len([_SRead for _SRead in SReads if not _SRead.isNULL()])-2*sqrt(1./len(H.GetVS())*(1-1./len(H.GetVS()))*len([_SRead for _SRead in SReads if not _SRead.isNULL()]))))
        #garbage = sys.stdout.write('Base Haplotype, S= {1:d}:\n\t{0}\n'.format('\n\t'.join(('\t'.join(_x for _x in H.GetGenotype(_pos))) for _pos in range(H.GetStart(), H.GetStop()+1)), G.GetS()))
        extend_H_branched = []
        extend_logprobs_branched = []
        uniques, priors, logrprobs, mins = GetProbTot(H, G, SReads, error, True, qscores, usecounts)
        if not uniques:
                return [H+Haplotypes(G.GetS(), G.GetS(), 0, 0, None, None, *['-' for _homologue in H.GetVS()])]  # skip extension if no extension has been possible
        #drop_out = [_n for _n, _x in enumerate(mins) if _x < count_threshold]
        drop_out = []
        #if len(drop_out)==len(uniques):
        #        garbage = sys.stderr.write('WARNING: No extension survived the minimum read support criterion at SNP {0:d} for {1:s}! Extension will be skipped at {0:d} for {1:s}!\n'.format(G.GetS()+1, samplename))
        #        return [H+Haplotypes(G.GetS(), G.GetS(), 0, 0, None, None, *['-' for _homologue in H.GetVS()])]  # skip extension if no homologue has the minimum number of supporting reads
        #for _n in sorted(drop_out, reverse=True): # discard haplotypes that have homologues not supported by a minimum number of reads 
        #        del uniques[_n]
        #        del priors[_n]
        #        del logrprobs[_n]
        #garbage = sys.stderr.write("SNP: {}\n".format(G.GetS()))
        #garbage = sys.stderr.write("\tTotal with last SNP REMOVED = {:7.19f}\n".format(total))
        #garbage = sys.stderr.write("\tTotal WITH last SNP with EQUAL priors = {:7.19f}\n".format(sum(map(exp,logrprobs))/len(logrprobs)))
        _norm = max(logrprobs)
        logrprobs_adj = [_x + _y - _norm for _x, _y in zip(logrprobs, log(priors))] # adjust P[SR(s)|Hp, H, eps] by its prior P[Hp|H, eps]. Optionally, max(logrprobs) is subtracted to prevent underflow
        #garbage = sys.stderr.write("\tTotal WITH last SNP with ACTUAL priors = {:7.19f}\n".format(sum(map(exp, logrprobs_adj))))
        _norm = loge(sum(exp(_x) for _x in logrprobs_adj))
        logHpprobs = [_x - _norm for _x in logrprobs_adj] # obtain p[Hp|SR(s), H, eps] by P[SR(s)|Hp, H, eps].P[Hp|H, eps]/P[SR(s)|H, eps]
        for _n, Hcandid in enumerate(uniques): # remove duplicate extensions that occur due to presence of similar homologues in H
                #garbage = sys.stderr.write('\tCandidate Extension:\n\t    {0}\n'.format('\t'.join(str(_x) for _x in Hcandid.GetGenotype(Hcandid.GetStop()))))
                #garbage = sys.stderr.write("\t    prob={:7.19f}, logprob= {:7.19f}\n".format(exp(logHpprobs[_n]), logHpprobs[_n]))
                if logHpprobs[_n] >= myrho:  # cut the extensions with an adjusted reads-probability lower than the threshold
                        extend_H_branched.append(Hcandid)
                        extend_logprobs_branched.append(logHpprobs[_n])
                        #garbage = sys.stderr.write("\t    Candidate Accepted!\n")
                else:
                        #garbage = sys.stderr.write("\t    Candidate Rejected by rho or kappa!\n")
                        pass
        if not extend_H_branched:
                garbage = sys.stderr.write('WARNING: No extension survived the threshold at SNP {0:d} for {1:s}!\n'.format(G.GetS()+1, samplename))
                _maxindex = logHpprobs.index(max(logHpprobs))
                extend_H_branched.append(uniques[_maxindex])
                extend_H_branched[-1].SetRL(logHpprobs[_maxindex])
                #return [H+Haplotypes(G.GetS(), G.GetS(), 0, 0, None, None, *['-' for _homologue in H.GetVS()])] # skip extension if no extension passes the branching threshold
        for _H, _logprob in zip(extend_H_branched, extend_logprobs_branched):  # Update the stored RL value of H during branching to\ 
                _H.SetRL(_logprob) # Update the RL of Hp, as noted in Berger et al. 2014 p. 5.      
        return extend_H_branched

def BranchPop(Hpop, Gpop, SReadspop, rho, error, recombination_rate, QscoresPop=[None, None, None], GenoConstraint=True, useglobal=True):
	""" Branch current pop haplotypes (Hm, Hf, Hc) at position s = Gm.s = Gf.s = Gc[1].S = ... = Gc[n].S, using the genotypes, the semi-reads of each pop member
	active at s, the recombination rate at s and the threshold rho. Base error rate, i.e. error, or quality scores are passed to calculate the 
	branching probabilities using the reads, and recombination_rate is used to incorporate the recombination events in those probabilities."""
	global AllGenos_m
	global AllGenos_f
	Hm, Hf = Hpop[0:2] 
	Hc = Hpop[2:] # Unpack the population candidate haplotypes to maternal, paternal and child haplotypes
	Gm, Gf = Gpop[0:2] # Unpack the population candidate genotypes to maternal, paternal and child haplotypes
	Gc = Gpop[2:]
	SReadsm, SReadsf = SReadspop[0:2] # Unpack the reads to maternal, paternal and child haplotypes
	SReadsc = SReadspop[2:] 
	qscoresm, qscoresf = QscoresPop[0:2] # Assign quality scores, if present, to parental and progeny reads
	qscoresc = QscoresPop[2:]
	qscoresc.extend([None for _i in range(0, len(Gc)-len(qscoresc))])
	if all(r.isNULL() for r in SReadsm+SReadsf+reduce(lambda x, y: x+y, SReadsc)):
		garbage = sys.stderr.write('WARNING: No semi-reads exist at SNP position {1:d}, s={0:d} with this pop!\n'.format(Gm.GetS()+1, Gm.GetPos()))
	ProbH = [_H.GetRL() for _H in Hpop]
	extend_logprobs_branched = []
	extend_pop_branched = []
	myrho = log(rho) # change rho to log_2 scale
	PopHaplos, PopLogProbs = [], []
	parental_branches = [[], []] # uniques, priors, logrprobs, logHpp_probs = [], [], [], [] # unique extensions of each parent, branched using Branch function, with the read support (P(SR|Hpparent, Hparent, eps))
	_parent_id = -1
	sema, lock , Parallel = GetSemaLock(useglobal)
	#_manager_id = rnd.randint(0, 100000)
	#garbage = sys.stderr.write('New thread manager hired with ID number: {0:d}!\n'.format(_manager_id))
	if Parallel:
		manager = mpthread.Manager()
	for _parental_hap, _parental_geno, _parental_semireads, _parental_qscore in zip((Hm, Hf), (Gm, Gf), (SReadsm, SReadsf), (qscoresm, qscoresf)):
		_parent_id += 1
		#garbge = sys.stderr.write("Branching the {0}:\n".format("mother" if _parent_id==0 else "father")) 
		if GenoConstraint:
			parental_branches[_parent_id].extend(Branch(_parental_hap, _parental_geno, _parental_semireads, rho, error, _parental_qscore,  usecounts=False, samplename="Mother" if _parent_id==0 else "Father", geno=GenoConstraint))
		else:	
			AllGenos = AllGenos_m if _parent_id == 0 else AllGenos_f
			GenoProbs = [None for _G in AllGenos[_parental_geno.GetS()]]
			for _i, _G in enumerate(AllGenos[_parental_geno.GetS()]):
				GenoProbs[_i]= GetProbReads(SubReads(_parental_semireads, set([_parental_geno.GetS()]), 1), Haplotypes(_parental_geno.GetS(), _parental_geno.GetS(), 0, 0, None, None, *[str(_g) for _g in _G]), error, True, qscoresm, False, 1) # P(r|G_parent)
			_norm = loge(sum(exp(_logprob) for _logprob in GenoProbs))
			GenoProbs = [_logprob -_norm for _logprob in GenoProbs] # P(G_parent|r) = P(r|G_parent)P(G_parent)/P(r)
			threads = []
			if Parallel:
				_branches = [manager.list() for _G in AllGenos[_parental_geno.GetS()]]		
			else:
				_branches = [[] for _G in AllGenos[_parental_geno.GetS()]]
			_genoThreshold = loge(0.25*1./len(GenoProbs))
			for _i, _G in enumerate(AllGenos[_parental_geno.GetS()]):
				if GenoProbs[_i] >= _genoThreshold:
					if Parallel: 
						t = mpthread.Process(target=thread_func, args = (sema, lock, _branches[_i], Branch, _parental_hap, Genotype(_parental_geno.GetS(), _parental_geno.GetPos(), *[str(_g) for _g in _G]), _parental_semireads, rho, error, _parental_qscore,  False, "Mother" if _parent_id==0 else "Father", GenoConstraint))
					else:
						t = threading.Thread(target=thread_func, args = (sema, lock, _branches[_i], Branch, _parental_hap, Genotype(_parental_geno.GetS(), _parental_geno.GetPos(), *[str(_g) for _g in _G]), _parental_semireads, rho, error, _parental_qscore,  False, "Mother" if _parent_id==0 else "Father", GenoConstraint))
					sema.acquire()
					threads.append(t)
					#garbage=sys.stderr.write('....Next thread operated by Manager {0:d}: thread number {1:d}\n'.format(_manager_id, len(threads)))
					t.start()
			for _thread in threads:
				_thread.join()
			if Parallel:
				_branches = [list(_lprox) for _lprox in _branches] 
			for _i in range(0, len(_branches)):
				for _j in range(0, len(_branches[_i])):
					_branches[_i][_j].SetRL(_branches[_i][_j].GetRL()+GenoProbs[_i]) # P(H|r) = P(H, G|r)  = P(H|G, r) P(G|r)	
			parental_branches[_parent_id] = reduce(lambda x, y: x+y, _branches)
			_approved = []
			for _i in range(0, len(parental_branches[_parent_id])): 
				if parental_branches[_parent_id][_i].GetRL()>=myrho:
					_approved.append(_i)
			if not _approved: # if not parental branch passes the branching threshold
				parental_branches[_parent_id]  = [sorted(parental_branches[_parent_id], key= lambda _x: -1*_x.GetRL())[0]]	
			else:
				parental_branches[_parent_id] = [parental_branches[_parent_id][_i] for _i in _approved]
			del _branches, _approved
	if GenoConstraint and set(Gm.GetGenes()).union(Gf.GetGenes())==set(['.']): # if both parents have missing genotypes, just use the progeny sequence reads to extend its haplotypes
		garbage = sys.stderr.write("WARNING: As both parental genotypes are missing, the children's phasing is extended just using its own reads at SNP "+str(Gc[0].GetS()+1)+"!\n") 
		if Parallel:
			children_branches = [manager.list() for _Hc in Hc]
		else:
			children_branches = [[] for _Hc in Hc] 	# unique extensions of the progeny that survive the branching threshold
		threads = []
		#garbage=sys.stderr.write('All threads finished by Manager {0:d}!\n'.format(_manager_id))
		for _i, _Hc in enumerate(Hc):
			sema.acquire()
			if Parallel:
				t = mpthread.Process(target=thread_func, args = (sema, lock, children_branches[_i], Branch, _Hc, Gc[_i], SReadsc[_i], rho, error, qscoresc[_i], False, "offspring "+str(_i))) 
			else:
				t = threading.Thread(target=thread_func, args = (sema, lock, children_branches[_i], Branch, _Hc, Gc[_i], SReadsc[_i], rho, error, qscoresc[_i], False, "offspring "+str(_i))) 
			threads.append(t)
			#garbage=sys.stderr.write('....Next thread operated by Manager {0:d}: thread number {1:d}\n'.format(_manager_id, len(threads)))
			t.start()
		for _thread in threads:
			_thread.join()
		if Parallel:
			children_branches = [list(_child_branches) for _child_branches in children_branches] # cast multiprocess manager lists to lists
		for _i, _child_branches in enumerate(children_branches):
			for _child in _child_branches:
				_child.SetMO(Hc[_i].GetMO())
				_child.SetPO(Hc[_i].GetPO())
		for _solutions in itertools.product(*[[_x for _x in range(0, len(_child_branches))] for _child_branches in children_branches]):
			PopHaplos.append(tuple(_parental_extension[0] for _parental_extension in  parental_branches)+
					tuple(children_branches[_i][_solution] for _i, _solution in enumerate(_solutions)))
	else:
		for _parents in itertools.product(*parental_branches): # all possible ways of selecting a pair of parental haplotypes 
			block_exceptions = 0
			children_branches = [[] for _Hc in Hc] 	# unique extensions of the progeny that survive the branching threshold
			if Parallel:
				_child_uniques = [manager.list() for _Hc in Hc]
				_child_priors = [manager.list() for _Hc in Hc]
				_child_logprobs = [manager.list() for _Hc in Hc]
			else:
				_child_uniques = [[] for _Hc in Hc]
				_child_priors = [[] for _Hc in Hc]
				_child_logprobs = [[] for _Hc in Hc]
			threads = []
			#garbage=sys.stderr.write('All threads finished by Manager {0:d}!\n'.format(_manager_id))
			for _i, _Hc in enumerate(Hc):
				sema.acquire()
				if Parallel:
					t = mpthread.Process(target = thread_func, args = (sema, lock, [_child_uniques[_i], _child_priors[_i], _child_logprobs[_i]], GetProbTotChild, _parents[0], _parents[1], _Hc, Gc[_i], SReadsc[_i], error, recombination_rate, True, qscoresc[_i], GenoConstraint))
				else:
					t = threading.Thread(target = thread_func, args = (sema, lock, [_child_uniques[_i], _child_priors[_i], _child_logprobs[_i]], GetProbTotChild, _parents[0], _parents[1], _Hc, Gc[_i], SReadsc[_i], error, recombination_rate, True, qscoresc[_i], GenoConstraint))
				threads.append(t)
				#garbage=sys.stderr.write('....Next thread operated by Manager {0:d}: thread number {1:d}\n'.format(_manager_id, len(threads)))
				try:
					t.start()
				except BlockException as e: # Exception occurs if Mendelian inheritance is violated for a progeny
					garbage = sys.stderr.write("WARNING:{0} Variant {1} will be eliminated from the extension of progeny {2}!\n".format(e, Gc[0].GetS(), _i+1))
					block_exceptions+=1
					_child_uniques[_i] = [_Hc+Haplotypes(Gc[_i].GetS(), Gc[_i].GetS(), 0, 0, _Hc.GetMO(), _Hc.GetPO(), *['-' for _homologue in _Hc.GetVS()])]
					_child_priors[_i] = [1]
					_child_logprobs[_i] = [log(0.1)]
			for _thread in threads:
				_thread.join()
			if Parallel:
				_child_uniques = [list(_x) for _x in _child_uniques]
				_child_priors = [list(_x) for _x in _child_priors]
				_child_logprobs = [list(_x) for _x in _child_logprobs]
			for _i in range(0, len(Hc)):
				if not _child_uniques[_i]:
					children_branches[_i].append(Hc[_i]+Haplotypes(Gc[_i].GetS(), Gc[_i].GetS(), 0, 0, Hc[_i].GetMO(), Hc[_i].GetPO(), *['-' for _homologue in Hc[_i].GetVS()]))  # skip child extension if no extension has been possible
				else:
					#garbage = sys.stderr.write("Tansmissions to child at SNP "+str(Gc[0].GetS())+":\n"+"\t"+"  ".join(str(_x) for _x in _child_uniques)+"\n\tPrior Weights: "+"  ".join(str(_x) for _x in _child_priors)+"\n\tLog probs: "+"  ".join(str(_x) for _x in _child_logprobs)+"\n")
					_norm = max(_child_logprobs[_i])
					_child_logprobs_adj = [_x + _y - _norm for _x, _y in zip(_child_logprobs[_i], log(_child_priors[_i]))]
					_norm = loge(sum(exp(_x) for _x in _child_logprobs_adj))
					logHpc_probs = [_x - _norm for _x in _child_logprobs_adj]
					_num = -1
					for _Hcp, _logHpc_prob in zip(_child_uniques[_i], logHpc_probs): # cut the child branches with prob less than rho
						_num += 1
						if _logHpc_prob >=myrho:
							children_branches[_i].append(_Hcp)
							children_branches[_i][-1].SetRL(_logHpc_prob)
					if not children_branches[_i]:
						garbage = sys.stderr.write('WARNING: No extension of progeny {1:d} survived the threshold at SNP position {2:d}, s={0:d}!\n'.format(Gm.GetS()+1, _i+1, Gm.GetPos()))
						_maxindex = logHpc_probs.index(max(logHpc_probs))
						children_branches[_i].append(_child_uniques[_i][_maxindex])
						children_branches[_i][-1].SetRL(logHpc_probs[_maxindex])
			if block_exceptions == len(Hc): # Eliminate a SNP from all haplotypes if it does not conform to Mendelian inheritance in any child
				garbage = sys.stderr.write("WARNING: Variant {0} will be eliminated from all of the extensions as no progeny conformed to Mendelian inheritance!\n".format(Gc[0].GetS()))
				PopHaplos = [(Hm+Haplotypes(Gm.GetS(), Gm.GetS(), 0, 0, None, None, *['-' for _homologue in Hm.GetVS()]), Hf+Haplotypes(Gf.GetS(), Gf.GetS(), 0, 0, None, None, *['-' for _homologue in Hf.GetVS()]))+ tuple(_Hc+Haplotypes(Gc[_i].GetS(), Gc[_i].GetS(), 0, 0, _Hc.GetMO(), _Hc.GetPO(), *['-' for _homologue in _Hc.GetVS()]) for _i, _Hc in enumerate(Hc))]
				break
			else:
				for _solutions in itertools.product(*[[_x for _x in range(0, len(_child_branches))] for _child_branches in children_branches]):
					PopHaplos.append(_parents+tuple(children_branches[_i][_solution] for _i, _solution in enumerate(_solutions)))
	#garbage=sys.stderr.write('All threads finished by Manager {0:d}!\n'.format(_manager_id))
	PopLogProbs = [sum(_H.GetRL() for _H in _Hpop) for _Hpop in PopHaplos] # Log(P(child, Parent1, Parent2| Hc, Hp1, Hp2, error, recombination_rate)) = Log(P(child | Parent1, Parent2, HHc, Hp1, Hp2, error, recombination_rate)) + Log(P(Parent1 | Hc, Hp1, Hp2, error, recombination_rate)) + Log(P(Parent2| Hc, Hp1, Hp2, error, recombination_rate)) = Log(P(child | Parent1, Parent2, Hc, Hp1, Hp2, error, recombination_rate)) + Log(P(Parent1 | Hp1, error)) + Log(P(Parent2| Hp2, error))
	_norm = loge(sum(exp(_x) for _x in PopLogProbs)) # Normalize the pop probabilities
	PopLogProbs = [_x - _norm for _x in PopLogProbs]
	PopLogProbsToPrune = [[_H.GetRL() for _H in _Hpop] for _Hpop in PopHaplos] # Log(P(child, Parent1, Parent2| Hc, Hp1, Hp2, error, recombination_rate)) = Log(P(child | Parent1, Parent2, HHc, Hp1, Hp2, error, recombination_rate)) + Log(P(Parent1 | Hc, Hp1, Hp2, error, recombination_rate)) + Log(P(Parent2| Hc, Hp1, Hp2, error, recombination_rate)) + Nt*Log(recombinationrate) + Nl*Log(1-recombination_Rate) = Log(P(child | Parent1, Parent2, Hc, Hp1, Hp2, error, recombination_rate)) + Log(P(Parent1 | Hp1, error)) + Log(P(Parent2| Hp2, error)) + Nt*Log(recombinationrate)+ Nl*Log(1-recombination_Rate)
	for _n, _Hppop_candid in enumerate(PopHaplos):
		#garbage = sys.stderr.write("SNP number={0}, Haplotype candid {1} {2}, prob={3:5.10f}\n".format(Gm.GetS(), _n+1, _Hppop_candid, exp(PopLogProbs[_n])))
		if PopLogProbs[_n] >= myrho: # discard extensions with a probability lower than the threshold
			extend_pop_branched.append(tuple(_HH.GetCopy() for _HH in _Hppop_candid))
			extend_logprobs_branched.append(PopLogProbsToPrune[_n])
		else:
			pass
	if not extend_pop_branched:
		garbage = sys.stderr.write('WARNING: No family extension survived the threshold at SNP {0:d} for the branch!\n'.format(Gm.GetS()+1))
		Suri = sorted(PopHaplos, key = lambda x: -1*sum(_x.GetRL() for _x in x))[0]
		extend_pop_branched.append(tuple(_HH.GetCopy() for _HH in Suri))
		extend_logprobs_branched.append(PopLogProbsToPrune[PopHaplos.index(Suri)])
	for _Pop, _logprob in zip(extend_pop_branched, extend_logprobs_branched):  # Update the stored RL value of the pop members to be used for pruning 
		_member = -1
		for _H in _Pop:
			_member+=1
			_H.SetRL(ProbH[_member]+_logprob[_member]) # Update the log RL of Hpm, Hpf and Hpc
	return extend_pop_branched

def GetProbReads(Reads, Vset, eps = 0.019, pplog = False, Quals = None, getcounts=False, min_read_length=2):
	""" Probability of a set of reads, i.e. P[R|Vset, eps] = Mult(P[r|Vset, eps] for r in R), assuming independence 
	and using GetlogProb(r, Vset, eps) (Berger et al. 2014, p. 4). If getcounts if True, also calculate the number of 
	reads assigned to each homologue."""
	try:
		if getcounts:
			if Quals:
				#print(Vset.GetStart())
				#print([(_Read, _Qual) for _Read, _Qual in zip(Reads, Quals)])
				probs, counts = list(zip(*[GetlogProb(_Read, Vset, eps, _Qual, getcounts, min_read_length) for _Read, _Qual in zip(Reads, Quals)]))
				#print(probs)
				#print(counts)
			else:
				probs, counts = list(zip(*[GetlogProb(_Read, Vset, eps, None, getcounts, min_read_length) for _Read in Reads]))
			Vset_counts = [sum(_counts[_n] for _counts in counts) for _n in range(0, len(Vset.GetVS()))]
			#print(Vset_counts)
			if pplog:
				#print("HABAL:", [sum(probs), Vset_counts])
				return [sum(probs), Vset_counts]
			return [exp(sum(probs)), Vset_counts]
		if Quals:
			probs = [GetlogProb(_Read, Vset, eps, _Qual, getcounts, min_read_length) for _Read, _Qual in zip(Reads, Quals)]
		else:
			probs = [GetlogProb(_Read, Vset, eps, None, getcounts, min_read_length) for _Read in Reads]
		if pplog: # return log(P[R|Vset, eps])
			return sum(probs)
		return exp(sum(probs))
	except (IndexError, ValueError) as e: # Error that occurs at the event that the Reads set is empty
		if "index 0 is out of bounds for axis 0 with size 0" in e.args[0] or "need more than 0 values to unpack" in e.args[0]:
			if getcounts:
				if pplog:
					return [0, [0 for _h in range(0, len(Vset.GetVS()))]]
				return [1, [0 for _h in range(0, len(Vset.GetVS()))]]
			if pplog:
				return 0
			else:
				return 1
		else:
			raise

def GetProbTot(H, G, Reads, error_rate, plog = False, Qscores = None, usecounts=False):
	"""Determine the set of distinct extensions, calculate their prior weights and report the read-probabilities conditional on each extension. In case usecounts is True, also report the minimum number of reads compatible with each homologue, so that the upstream functions may set a threshold on the minimum number of reads compatible with each homologue."""
	perm = makePermutation(G) # distinct permutations of a genoptype
	probs = []   # the probability of SR(s) conditional on (Hp, H, eps)
	weights = [] # the prior pobability of Hp conditional on (H, eps)  
	Uniques = [] # distinct Hp's
	Mins = []
	for P in perm:
		Hp = H + P
		if Hp not in Uniques:
			Uniques.append(Hp)
			if usecounts:
				_prob, _counts = GetProbReads(Reads, Hp, error_rate, plog, Qscores, True)
				probs.append(_prob)
				Mins.append(min(_counts))
			else:
				probs.append(GetProbReads(Reads, Hp, error_rate, plog, Qscores))
			_npVset = []
			for _v in Hp.GetVS():
				_npv = array(_v)
				_npVset.append(npdel(_npv, npwhere(_npv=='-')).tolist())
			try:
				weights.append(exp(GetLogProbH(Haplotypes(1, 2, 1, log(len(set(itertools.permutations(tuple((
				_v[-2],_v[-1]) for _v in _npVset))))), None, None, *[(_v[-2],_v[-1]) for _v in _npVset]))))
			except IndexError:
				weights.append(1)
			#weights.append(exp(GetLogProbH(Hp)))                     
			#weights.append(1)
		else:
			pass
	wsum = float(sum(weights))
	weights = [_w/wsum for _w in weights]
	return Uniques, weights, probs, Mins

def GetProbTotChild(Hpm, Hpf, Hc, Gc, Readsc, error, recombination_rate, plog = True, Qscoresc = None, GenoConstraint=True):
	"""Determine the set of distinct extensions of a child conditional on its extended parents, calculate their prior weights and report the read AND recombination support for each child extension. It is assumed that at least one of the parental genotypes is available at s.""" 
	global transmit_parent_m 
	global transmit_parent_f 
	global child_haploid
	global transmit_parent
	probs = []   # the probability of SR(s) conditional on (Hp, H, eps)
	weights = [] # the prior pobability of Hp conditional on (H, eps, recombination_rate)  
	Uniques = [] # distinct Hp's
	if GenoConstraint and set(Gc.GetGenes())=={'.'}:
		garbage = sys.stderr.write('WARNING: Child\'s genotype is missing at position {1:d}! Its phasing extension will be skipped at s={0:d}!\n'.format(Gc.GetS()+1, Gc.GetPos()))
		return Uniques, weights, probs
	biparental = False
	current_Gm = Hpm.GetGenotype(Hpm.GetStop()) # current maternal allele (after extension)
	current_Gf = Hpf.GetGenotype(Hpf.GetStop()) # current paternal allele
	last_Gc = Hc.GetGenotype(Hc.GetStop()) # last child allele (before extension)
	ploidy_m = len(current_Gm) # It is assumed that the first haploid part of the child homologues has maternal descent and the second haploid paternal 
	ploidy_f = len(current_Gf)
	ploidy_c = (ploidy_m+ploidy_f)//2
	last_determined_SNPpos_m = Hpm.GetStop()-1 # the last SNP position with determined alleles in the maternal haplotypes, the default value
	last_determined_SNPpos_f = Hpf.GetStop()-1 # the last SNP position with determined alleles in the paternal haplotypes, the default value
	if set(Gc.GetGenes())=={'.'}: #'-' stands for complete genotype missing, '.' stands for partial missing or null allele
		Gc = Genotype(Gc.GetS(), Gc.GetPos(), *['-' for _x in range(0, ploidy_c)])
	if not transmit_parent_m:
		transmit_parent_m = list(itertools.combinations(range(0, ploidy_m), ploidy_m//2)) # Initialize transmit_parent_m 
	if not transmit_parent_f:
		transmit_parent_f = list(itertools.combinations(range(0, ploidy_f), ploidy_f//2)) # Initialize transmit_parent_f
	if not child_haploid:
		child_haploid = list(itertools.combinations(range(0, ploidy_c), ploidy_c//2)) # Initialize child_haploid
	if not transmit_parent:
		transmit_parent = list(itertools.product(transmit_parent_m, transmit_parent_f)) # Initialize transmit_parent  
	transmit_parent_G = []
	attempt_num=1
	OrigG = GenoConstraint
	while not transmit_parent_G and attempt_num<=2:
		if attempt_num>1:
                        GenoConstraint = False
                if '-' not in set(current_Gm).union(set(current_Gf)): # if no parental genotype is missing, consider all of the possible IBD potentiae
                        biparental = True
                        if GenoConstraint: # only allow transmissions that are compatible with the current genotype of the child 
                                transmit_parent_G = [_x for _x in transmit_parent if Counter([current_Gm[_y] for _y in _x[0]]+[current_Gf[_z] for _z in _x[1]])==Counter(Gc.GetGenes())]
                        else:
                                transmit_parent_G = [_x for _x in transmit_parent]
                elif attempt_num == 1: # Can only occur if GenoConstraint has been originally True 
                        if '-' in set(current_Gm): # impute maternal alleles from the child and the father
                                transmit_parent_G = [(tuple('-' for _i in range(0, ploidy_m//2)), _x) for _x in transmit_parent_f if any(Counter([current_Gf[_z] for _z in _x])==Counter([Gc.GetGenes()[_hh] for _hh in _h]) for _h in child_haploid)]
                        else:   # impute paternal alleles from the mother and the child
                                transmit_parent_G = [(_x, tuple('-' for _i in range(0, ploidy_f//2))) for _x in transmit_parent_m if any(Counter([current_Gm[_y] for _y in _x])==Counter([Gc.GetGenes()[_hh] for _hh in _h]) for _h in child_haploid)]
                else:
                        transmit_parent_G = []
		attempt_num+=1
        GenoConstraint = OrigG
        if attempt_num>2:
                error = min(0.999, error*1.1)
                recombination_rate = min(0.499, recombination_rate+0.005)
        if not transmit_parent_G:
                garbage = sys.stderr.write("WARNING: No extension of child haplotypes was compatible with the given parental and/or child genotypes assuming Mendelian inheritance! Child extension will be skipped at position {1:d}, s={0:d}!".format(Gc.GetS()+1, Gc.GetPos()))
                return Uniques, weights, probs
		#_return = GetProbTot(Hc, Gc, Readsc, error, True, Qscoresc)[:-1]
		#for _pro in range(0, len(_return[-1])):
		#	if plog:
		#		_return[-1][_pro]+=loge(0.1)
		#	else:
		#		_return[-1][_pro]*=0.1
		#for _sol in range(0, len(_return[0])):
		#	_return[0][_sol].SetMO(Hc.GetMO())
		#	_return[0][_sol].SetPO(Hc.GetPO())
		#return _return 
	for _t in transmit_parent_G:
		if biparental:
			Genotypes_to_pass = [[current_Gm[_y] for _y in _t[0]]+[current_Gf[_z] for _z in _t[1]]] # Genotype to pass will be the same as Gc.GetGenes() if GenoConstraint is True 
		elif '-' in set(current_Gm): # Try to impute the missing parent's alleles from the child alleles and the alleles of the other parent
                        fake_genotype_m = [_code for _code in Gc.GetGenes()] # If it is not possible, skip the extension of child at that position
                        try:
                                for _paternal_allele in (current_Gf[_z] for _z in _t[1]):
                                        fake_genotype_m.remove(_paternal_allele)
                                Genotypes_to_pass = [fake_genotype_m + [current_Gf[_z] for _z in _t[1]]] # Genotype to pass will be the same as Gc.GetGenes() if GenoConstraint is True. Only the order is changed so that the alleles with maternal descent come first. 
                        except ValueError as e:
                                garbage = sys.stderr.write("WARNING: "+e.args[0]+"\n")
                                if 'list.remove(x): x not in list' in e.args[0]:
                                        Genotypes_to_pass = []
                else:
                        fake_genotype_f = [_code for _code in Gc.GetGenes()]
                        try:
                                for _maternal_allele in (current_Gm[_y] for _y in _t[0]):
                                        fake_genotype_f.remove(_maternal_allele)
                                Genotypes_to_pass = [[current_Gm[_y] for _y in _t[0]] + fake_genotype_f] # Genotype to pass will be the same as Gc.GetGenes() if GenoConstraint is True 
                        except ValueError as e:
				garbage = sys.stderr.write("WARNING: "+e.args[0]+"\n")
                                if 'list.remove(x): x not in list' in e.args[0]:
                                        Genotypes_to_pass = []
			#print("$$$!\n", Genotypes_to_pass, fake_genotype_f,"$$$$$$$$$\n#############")
		for _Genotype_to_pass in Genotypes_to_pass: # Genotype_to_pass could have length >1 in case genoconstraint is False and a parental genotype is missing
			No_recom_info_m = bool(len(set(current_Gm))==1 or '-' in _t[0])
			No_recom_info_f = bool(len(set(current_Gf))==1 or '-' in _t[1])
			#print("Mother="+str(No_recom_info_m))
			#print("Father="+str(No_recom_info_f))
			P = Haplotypes(Gc.GetS(), Gc.GetS(), 0, loge(len(set(itertools.permutations(_Genotype_to_pass)))), None if No_recom_info_m else _t[0], None if No_recom_info_f else _t[1], *_Genotype_to_pass)
			#print("P="+str(P))
			#print(repr(P))
			Hpc = Hc + P
			#print("***********:\n",P, Hc, Hpc,'**************\n******************')
			_recombination_m = [] if (Hc.GetMO() is None or No_recom_info_m) else [0 if _origin == _new else 1 for _origin, _new in zip(Hc.GetMO(), _t[0])] # 0 if both the current allele and the preceding allele can be assigned to the same maternal descent, else 1. If the maternal descent if unknown, e.g. if the mother is homozygous, ignore maternal recombination. 
			_recombination_f = [] if (Hc.GetPO() is None or No_recom_info_f) else [0 if _origin == _new else 1 for _origin, _new in zip(Hc.GetPO(), _t[1])] # 0 if both the current allele and the preceding allele can be assigned to the same paternal descent, else 1. If the paternal descent if unknown, e.g. if the father is homozygous, ignore paternal recombination.
			#if plog:
			#	recombination_prob = sum([loge(recombination_rate) if _r==1 else loge(1-recombination_rate) for _r in (_recombination_m+_recombination_f)]) # Change the prior weights by penalizing recombination events with the cost -loge(recombination_rate)
			#	_prob_Hc_reads_recombination = _prob_Hc_reads + recombination_prob
			#else:
			#	recombination_prob = array([recombination_rate if _r==1 else (1-recombination_rate) for _r in (_recombination_m+_recombination_f)]).prod() # Change the prior weights by penalizing recombination events with the cost -loge(recombination_rate)
			#	_prob_Hc_reads_recombination = _prob_Hc_reads * recombination_prob
			recombination_weight = array([recombination_rate if _r==1 else (1-recombination_rate) for _r in (_recombination_m+_recombination_f)]).prod() # Change the prior weights by penalizing recombination events with the cost -loge(recombination_rate)
			#print(_recombination_m, _recombination_f, recombination_weight)
			#print("Hc:"+str(Hc.GetMO())+'\t'+str(Hc.GetPO()))
			#print("Hpc:"+str(Hpc.GetMO())+'\t'+str(Hpc.GetPO()))
			#print("T:"+str(_t[0])+'\t'+str(_t[1]))
			if Hpc not in Uniques:
				Uniques.append(Hpc)
				probs.append(GetProbReads(Readsc, Hpc, error, plog, Qscoresc))
				#probs.append(_prob_Hc_reads_recombination)
				#weights.append(1)
				_npVset = []
				for _v in Hpc.GetVS():
					_npv = array(_v)
					_npVset.append(npdel(_npv, npwhere(_npv=='-')).tolist())
				try:
					weights.append(exp(GetLogProbH(Haplotypes(1, 2, 1, log(len(set(itertools.permutations(tuple((_v[-2],_v[-1]) for _v in _npVset))))), None, None, *[(_v[-2],_v[-1]) for _v in _npVset]))))
				except IndexError:
					weights.append(1)
				#weights.append(exp(GetLogProbH(Hpc)))		
				weights[-1]*=recombination_weight
			elif _recombination_m or _recombination_f:
				#probs[Uniques.index(Hpc)] = max(probs[Uniques.index(Hpc)], _prob_Hc_reads_recombination)
				_npVset = []
				for _v in Hpc.GetVS():
					_npv = array(_v)
					_npVset.append(npdel(_npv, npwhere(_npv=='-')).tolist())
				try:
					new_weight = exp(GetLogProbH(Haplotypes(1, 2, 1, log(len(set(itertools.permutations(tuple((_v[-2],_v[-1]) for _v in _npVset))))), None, None, *[(_v[-2],_v[-1]) for _v in _npVset])))*recombination_weight
				except IndexError:
					new_weight = recombination_weight
				_Hpindx = Uniques.index(Hpc)
				if weights[_Hpindx] < new_weight:
					weights[_Hpindx] = new_weight
					if not No_recom_info_m:
						Uniques[_Hpindx].SetMO(Hpc.GetMO())
					if not No_recom_info_f:
						Uniques[_Hpindx].SetPO(Hpc.GetPO())
					#print("Updated Hpc:"+str(Hpc.GetMO())+'\t'+str(Hpc.GetPO()))
			else:
				pass
	wsum = float(sum(weights))
	weights = [_w/wsum for _w in weights]
	return Uniques, weights, probs

#def GetVERmat(Hlst):
#	""" Returns a list equal in size to Hlst, containing the sum of VER scores for each haplotype in Hlst compared to all other haplotypes in the list."""
#	VER_MATRIX = [[0 for _HH in Hlst] for _H in Hlst]
#	for _i in range(0, len(VER_MATRIX)):
#		for _j in range(_i+1, len(VER_MATRIX[0])):
#			VER_MATRIX[_i][_j]=getVER(Hlst[_i].GetVS(), Hlst[_j].GetVS())
#	for _i in range(1, len(VER_MATRIX)): # as the matrix is symmetric, fill in the lower diagonal part using the upper diagonal
#		for _j in range(0, _i):
#			VER_MATRIX[_i][_j]=VER_MATRIX[_j][_i]
#	return [sum(_VER_SCORES_FOR_H) for _VER_SCORES_FOR_H in VER_MATRIX]  

def makePermutation(G, lognum=False):
	""" Return all of the permutations possible for a genotype, e.g. (1,1,1,0), (0,1,1,1), (1,1,0,1) and 
	(1,0,1,1) for (1,1,1,0)."""
	_perm_set = set(itertools.permutations(G.GetGenes()))
	log_perm_cs = loge(len(_perm_set)) # logarithm of the number of permutations possible
	if lognum:
		return log_perm_cs
	return list(Haplotypes(G.GetS(), G.GetS(), 0, log_perm_cs, None, None, *_x) for _x in _perm_set)

def PrunePop(Poplst, kappa, error_rate, alpha = None, subreads=[[], [], []], QS=[(),(),()], Het=False):
	""" Prune a set of trio haplotypes using their log Relative Likelihoods and the pruning rate kappa."""
	try:
		#print(Poplst)
		RLs = [[_H.GetRL() for _H in _Pop] for _Pop in Poplst]
		#for _RLs in RLs:
		#	for _mem in range(0, len(_RLs)):
		#		if isnan(_RLs[_mem]):
		#			_RLs[_mem]=-1e2
		#RLs = array(list(_Pop[0].GetRL() for _Pop in Poplst))
		#RLs = array(list(sum(GetProbReads(_subreads, _H, error_rate, True, _qs) for _subreads, _H, _qs in zip(subreads, _T, QS)) for _T in Poplst)) # Calculate RL[H|R(SubReads), error_Rate] directly using formula (3) in Berger et al. 2014 p. 5.
		#RLs = list(sum(GetProbReads(_subreads, _H, error_rate, True, _qs) for _subreads, _H, _qs in zip(subreads, _T, QS)) for _T in Poplst) # Calculate RL[H|R(SubReads), error_Rate] directly using formula (3) in Berger et al. 2014 p. 5.
		#print("BEFORE NORM:", RLs)
		#print(sum(exp(rl) for rl in RLs))
		#_norm = max(RLs)
		#RLs = [_rl - _norm for _rl in RLs] # In case RLs are all very small and their sum is rounded to zero, must multiply them by a number large enough to be able to work with them
		#_norm = loge(sum(exp(_rl) for _rl in RLs)) 
		#RLs = [_rl - _norm for _rl in RLs]
		#maxprob = max(RLs)
		#maxprob_m = max(_RLs[0] for _RLs in RLs)
		#maxprob_f = max(_RLs[1] for _RLs in RLs)
		#maxprob_c = max(_RLs[2] for _RLs in RLs)
		maxprob = max(sum(_RLs) for _RLs in RLs)
		#print("AFTER NORM:", RLs)
		#print(sum(exp(rl) for rl in RLs))
	except ValueError as e:
		raise BlockException(''.join(e.args)+'\n'+"Pruning could not be done! "+'\n')
	pruned = [] # the pruned list of trios
	k = log(kappa)
	for _num, _P in enumerate(Poplst):
		if sum(RLs[_num]) >= k + maxprob: 
		#if RLs[_num][0]>= k + maxprob_m and RLs[_num][1]>=k+maxprob_f and  RLs[_num][2]>=k+maxprob_c:
		#if RLs[_num] >= k + maxprob:
			pruned.append(_P)
			#for _H in pruned[-1]:
			#	_H.SetRL(RLs[_num])
		else:
			pass
	non_empty_reads = [[_r for _r in _subreads if not _r.isNULL()] for _subreads in subreads]
	non_empty_reads += [[] for _extra in range(0, len(Poplst[0])-len(non_empty_reads))]
	QS += [() for _extra in range(0, len(Poplst[0])-len(QS))] 
	Separate_MEC_Scores = [] # Separate MEC scores for each population member
	norms = [0 for _member in Poplst[0]]
        for _num in range(0, len(pruned)):
		_Separate_MEC_Scores = [] # Separate MEC scores for each population member
                for _mem in range(0, len(pruned[_num])):
			_Separate_MEC_Scores.append(getMEC(non_empty_reads[_mem], pruned[_num][_mem]))
                        norms[_mem]=norms[_mem]+exp(pruned[_num][_mem].GetRL())
		Separate_MEC_Scores.append(_Separate_MEC_Scores)
        norms = log(norms) # Normalize the RLs to get proper probs
	for _num in range(0, len(pruned)):
		for _mem in range(0, len(pruned[_num])):
			pruned[_num][_mem].SetRL(pruned[_num][_mem].GetRL()-norms[_mem])
	HAPandMEC = sorted(zip(pruned, Separate_MEC_Scores), key = lambda x: (sum(x[1])/5.-sum(_x.GetRL() for _x in x[0]), -1*round(sum(_x.GetRL() for _x in x[0]), 9)))
	if Het: 
		HAPandMEC = [(_pruned, [getMEC(non_empty_reads[_mem], _pruned[_mem], Het=Het) for _mem in range(0, len(_pruned))]) for _pruned in [_x[0] for _x in HAPandMEC]] #MEC scores to be reported. If homozygous SNPs are filtered, they should not influence the reported MEC
	#VER_Scores = [0 for _T_pruned in pruned]
	#VER_weights = [1, 1, 1]
	#for _indx in range(0, len(pruned[0])):
	#	VER_Scores = [_x+VER_weights[_indx]*_y for _x, _y in zip(VER_Scores, GetVERmat([_T_pruned[_indx] for _T_pruned in pruned]))] # Sort using the VER sum of all the members of the pedigree
	#VER_Scores = GetVERmat([_T_pruned[2] for _T_pruned in pruned]) # Sort using the VER sum of the children
	#HAPandMEC = sorted(zip(pruned, MEC_Scores, VER_Scores), key = lambda x: (-1*round(x[0][0].GetRL(), 3), x[2], sum(x[1])), reverse = False)
	if (alpha is not None) and len(pruned)>(1-alpha)*len(Poplst):
		PrunedMEC = list(zip(*sorted(HAPandMEC[0:int((1-alpha)*len(Poplst))+2], key = lambda x: (-1*round(x[0][0].GetRL(), 4), x[1])))) #sorted(pruned, key=lambda x: x.GetRL(), reverse=True)
		return PrunedMEC[0], PrunedMEC[1]
	PrunedMEC = list(zip(*HAPandMEC))
	#print(PrunedMEC)
	return PrunedMEC[0], PrunedMEC[1]

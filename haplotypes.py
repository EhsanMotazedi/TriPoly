# Written by Ehsan Motazedi, Wageningen UR, 11-08-2016.
# Last Updated: 22-10-2017

import copy
import functools
import itertools
import math
import random as rnd
from datetime import datetime
from genotypes import Genotype
from reads import Read
from logprob import _xor, diffVec, Hamming

def getMEC(Rlst, Haplo, Het=False):
	"""Calculate the Minimum Error Correction (MEC) with regard to a haplotype for a list of reads.""" 
	try:
		if not all(isinstance(r, Read) for r in Rlst):
			raise TypeError
	except TypeError as e:
		e.args=('The first argument must be a list of Read objects!',)
		raise
	if not isinstance(Haplo, Haplotypes):
		raise TypeError('The second argument must be a Haplotypes object!')	
	return sum(getMinMatchFlip(r, Haplo, Het=Het)[0] for r in Rlst)

def getMinMatchFlip(r, Haplo, nulllist=['.','-'], Het=False):
	"""Assign a homologue to a read using the Minimum Error Correction (MEC) criterion. Return the MEC and the\n
number of the matching homologue. With Het set to true, only heterozygous sites are counted for MEC.""" 
	if not isinstance(r, Read) or r.isNULL():
		raise TypeError('The first argument must be a non-empty Read object!') 
	if not isinstance(Haplo, Haplotypes):
		raise TypeError('The second argument must be a Haplotypes object!')
	if r.GetEnd()<Haplo.GetStart() or r.GetBegin()>Haplo.GetStop():
		raise ValueError('The read coordinates do not match that of the homologues!')
	base_1h, base_nh = max(0, r.GetBegin()-Haplo.GetStart()), max(0, min(r.GetEnd()-Haplo.GetStart()+1, Haplo.GetStop()-Haplo.GetStart()+1)) 
	Min = float('inf')
	N = 0
	bases = r.GetalleleSet(Haplo.GetStart(), Haplo.GetStop())
	if len([_base for _base in bases if _base not in nulllist])<2:
		return 0, 0
	HaploVS = Haplo.GetVS()
	if Het:
		Homozygous_HaploVS = [] # homozygous sites of the haplotype
		for _l in range(0, len(HaploVS[0])):
			if len(set(_haplo_alleles[_l] for _haplo_alleles in HaploVS if _haplo_alleles[_l] not in nulllist))<2: 
				Homozygous_HaploVS.append(_l)
		HaploVS = [tuple(_allele if _l not in Homozygous_HaploVS else '-' for _l, _allele in enumerate(_haplo_alleles)) for _haplo_alleles in HaploVS]
	if len([_x for _x in HaploVS[0] if _x not in nulllist])<2:
		return 0, 0
	for _n in range(0, len(HaploVS)):
		_vec = HaploVS[_n]
		if isinstance(_vec, int) or isinstance(_vec, str):
			_vectmp = []
			_vectmp.append(_vec)
			_vec = copy.copy(_vectmp)
		#print(''.join(str(_v) for _v in _vec)[base_1:base_n])
		_homolo_score = Hamming(bases, ''.join(str(_v) for _v in _vec)[base_1h:base_nh])[0]
		if _homolo_score < Min:
			Min = _homolo_score
			N = _n
	return Min, N

class Haplotypes(object):
	""" Methods for vector set of haplotypes, i.e. records of homologue alleles (len=ploidy), Relative Likelihood RL of the haplotype, its start and stop positions."""
	def __init__(self, start, stop, RL, logcs, mo, po, *homolos): # Each homologue should be given as a vector. *(1,2,3,...) will correspond to multiple homologues, while (1,2,3,...) will be interpreted as one homologue. 
		if (not isinstance(start, int)) or (not isinstance(stop, int)):
			raise TypeError('The start and stop of the haplotypes must be integer coordinates!')
		if start>stop:
			raise ValueError('The start coordinate must be smaller than or equal to the stop!')
		self.start = start
		self.stop = stop
		self.CS = logcs
		self.mo = mo # the number of maternal homologues passed to a child haplotype at the last position in the haplotype, e.g. (1,2) for a tertploid. Set to None if not applciable or unknown.
		self.po = po # the number of paternal homologues passed to a child haplotype at the last position in the haplotype, similar to mo.
		try:
			self.RL = float(RL)
		except ValueError as e:
			e.args = ('The relative likelihood must be a real number!',)
			raise
                try:
                    if any(len(_x)!=(self.stop-self.start+1) for _x in homolos):
                            raise ValueError('Length of the homologues must be compatible with the given start and stop positions!')
                except TypeError as e:
                    if "object of type 'int' has no len"  in e.args[0]:
                        if any(self.stop-self.start+1!=1 for _x in homolos):
                            raise ValueError('Length of the homologues must be compatible with the given start and stop positions!')
		self.vecset = tuple(_x for _x in homolos)
	def __add__(self, other):
		if (not other) or other=='0' or other=='[]' or other=='()': # The null elements of addition
			return self
		if not isinstance(other, Haplotypes):
			raise TypeError('Haplotype objects could be added just to Haplotype objects!')
		_newstop, _newstart, _newRL = self.stop, self.start, self.RL+other.RL
		if (other.start<=self.stop) and other.stop>=self.start:
			raise ValueError("Can not add phasings when their variants' positions overlap!")
		if len(other.vecset)!=len(self.vecset):
			if set(other.vecset).issubset({'.','-'}): # make the addition of undetermined SNPs to the haplotypes possible
				other = copy.deepcopy(other)
				other.vecset = tuple('-' for _h in self.vecset)
			elif set(self.vecset).issubset({'.','-'}): # make the addition of undetermined SNPs to the haplotypes possible
				self = copy.deepcopy(self)
				self.vecset = tuple('-' for _h in other.vecset)
			else:
				raise ValueError("Can not add phasings when the ploidy levels are different!")
		_new_homolos=[]
		if other.start>self.stop:
			_newstop = other.stop
			for _x in self.vecset:
				_new_homolos.append(tuple(_y for _y in _x))
				for _s in range(self.stop+1, other.start):
					_new_homolos[-1]+=('-',)
			for _n, _y in enumerate(other.vecset):
				if isinstance(_y, int):
					_new_homolos[_n]+=tuple(str(_y))
				else:
					_new_homolos[_n]+=tuple(_y)
		else:
			_newstart = other.start
			for _x in other.vecset:
				_new_homolos.append(tuple(_y for _y in _x))
				for _s in range(other.stop+1, self.start):
					_new_homolos[-1]+=('-',)
			for _n, _y in enumerate(self.vecset):
				if isinstance(_y, int):
					_new_homolos[_n]+=tuple(str(_y))
				else:
					_new_homolos[_n]+=tuple(_y)
		if (self.CS is None) or (other.CS is None):
			_newlogCS = None
		else:
			_newlogCS = self.CS+other.CS
		if other.po is None: # update the paternal origin as it belongs to the last position
			newpo = self.po
		else:
			newpo = other.po
		if other.mo is None: # update the maternal origin as it belongs to the last position
			newmo = self.mo
		else:
			newmo = other.mo
		return Haplotypes(_newstart, _newstop, _newRL, _newlogCS, newmo, newpo, *_new_homolos)

	def __radd__(self, other):
		if (not other) or other=='0' or other=='[]' or other=='()': # The null elements of addition
			return self
		return self.__add__(other)
	
	def __repr__(self):
		return('['+','.join(str(_x) for _x in (self.start, self.stop, self.RL, self.CS, self.mo, self.po))+
			',*['+','.join(str(_x) for _x in self.GetVS())+']]')
	
	def __str__(self):
		return("start_pos= {0}, stop_pos= {1}, RL= {2}, CS= {3}, MO= {4}, PO= {5}, homologues={6}".format(self.start, self.stop, self.RL, self.CS, 'NA' if self.mo is None else '|'.join(str(_x) for _x in self.mo),
		'NA' if self.po is None else '|'.join(str(_x) for _x in self.po), '|'.join(''.join(str(_x) for _x in list(str(_y) if isinstance(_y, int) else _y)) for _y in self.GetVS())))
	
	def __eq__(self, other):
		return isinstance(other, self.__class__) and self.start == other.start and self.stop == other.stop and sorted(_homolo for _homolo in self.vecset) == sorted(_homolo for _homolo in other.vecset)
	
	def __ne__(self, other):
		return not self.__eq__(other)

	def ChangeVS(self, *homolos):
		_new = homolos
		self.vecset = tuple(_x for _x in _new)
	
	def GetCopy(self):
		""" make a copy of a haplotype object"""
		return Haplotypes(self.start, self.stop, self.RL, self.CS, self.mo, self.po, *self.vecset)

	def GetlogCS(self):
		""" Return the logarithm of the product of the number of permutations at all SNP sites."""
		return self.CS

	def GetGenotype(self, pos):
		""" Return the alleles at a specific position."""
		if self.start>pos or pos>self.stop:
			raise ValueError("The specified SNP number is not within the range of the haplotype!")
		return tuple(_homolo[pos-self.start] for _homolo in self.vecset)
	
	def GetLen(self):
		return self.stop-self.start+1

	def GetRL(self):
		return self.RL
	
	def GetStop(self):
		return self.stop
	
	def GetStart(self):
		return self.start
	
	def GetVS(self):
		return tuple(_homolo for _homolo in self.vecset)
	
	def SetRL(self, prob):
		""" Set the RL field of a haplotype object."""
		self.RL = prob

	def UpdateCS(self, log_cs):
		""" Update the logarithm of the product of permutations when a new SNP is added to the haplotype."""
		self.CS+=log_cs

	def DelLastSNP(self):
		""" Delete the last SNP from the haplotype."""
		self.vecset = tuple(_homolo[0:-1:1]+('-',) for _homolo in self.vecset)

	def GetMO(self):
		""" Give the maternal origin."""
		return self.mo

	def GetPO(self):
		""" Give the paternal origin."""
		return self.po

	def SetMO(self, tup):
		""" Set the maternal origin to tup."""
		if tup is None:
			self.mo = None
			return
		if self.mo is not None and len(tup)!=len(self.mo):
			raise TypeError("The number of maternally inherited chromosomes must be equal to "+str(len(self.mo))+"!")
		self.mo = tuple(int(_hnum) for _hnum in tup) 

	def SetPO(self, tup):
		""" Set the paternal origin to tup."""
		if tup is None:
			self.po = None
			return
		if self.po is not None and len(tup)!=len(self.po):
			raise TypeError("The number of maternally inherited chromosomes must be equal to "+str(len(self.po))+"!")
		self.po = tuple(int(_hnum) for _hnum in tup) 

def Gametogenesis(H, gamete_ploidy=None):
    """ Generates gametes from a haplotype H, through balanced or unbalanced meiosis (the latter in case ploidy level is odd)."""
    ploidy = len(H.GetVS())
    if gamete_ploidy is not None:
        garbage = sys.stderr.write("WARNING: The given ploidy level for the gamete will only be used if the precursor ploidy-level if odd, leading to unbalanced meiosis!\n")
        assert gamete_ploidy>=1 and isinstance(gamete_ploidy, int), "The given ploidy level of the gamete must be a positive integer number!"
        assert gamete_ploidy<=ploidy, "The given ploidy level of the gamete must be less than or equal to the ploidy level of its precursor!"
    if ploidy % 2 != 0: # unbalanced meiosis
        if (gamete_ploidy is None): # if no ploidy-level is specified for the gamete, set it to random in case the meiosis is not balanced
            a = [int(_x) for _x in repr(datetime.now()).split('(')[1:][0].split(')')[0:-1][0].split(',')]
            rnd.seed()
            rndseed = int((101*rnd.random()*sum(a[0:3])+953*rnd.random()*sum(a[3:6]))*a[-1])
            rnd.seed(rndseed)
            func = rnd.choice([math.floor, math.ceil])
            return [_gamete for _gamete in itertools.combinations(H.GetVS(), int(func(ploidy/2.)))]
        else:
            return [_gamete for _gamete in itertools.combinations(H.GetVS(), gamete_ploidy)]
    return [_gamete for _gamete in itertools.combinations(H.GetVS(), ploidy//2)]

#a = Haplotypes(0,3,0.4,9, None, None, (1,1,1,0),(1,1,1,0),(1,1,1,1),(1,1,1,1))
#b = Haplotypes(7,8,0.9,36, None, None, (1,1),(1,0),(0,0),(0,0))
#c = Haplotypes(11,12,1,3, None, None, (0,1),(0,1),(1,0),(0,1))
#sum([a,b,c])
#print(sum([a,b,c]))
#Haplo=Haplotypes(8, 12, 0, 0, (1,0,'-',1,1), (0,0,'-',0,0), (1,0,'-',1,1), (0,1,'-',0,0))
#r=Read({5:1, 7:1, 9:0, 10:0})
#getMinMatchFlip(r, Haplo)

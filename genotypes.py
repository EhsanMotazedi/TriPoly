# Written by Ehsan Motazedi, Wageningen UR, 04-11-2016.
# Last Updated: 30-08-2018

import copy
import re
import sys
from collections import Counter
from Bio import SeqIO

class Genotype(object):
	def __init__(self, numS, pos, *genotype):
		try:
			self.s = int(numS)
			self.pos = int(pos)
		except (ValueError, TypeError) as e:
			e.args+=("The SNP id and SNP position must be inetgers!",)
		self.genotype = copy.deepcopy(genotype)
	def __repr__(self):
		return '['+','.join(str(_x) for _x in (self.s, self.pos))+','+','.join(str(_x) for _x in self.genotype)+']'
	def __str__(self):
		return "SNP number = {}, SNP position = {}, Genotype = {}".format(self.s, self.pos, 
				'/'.join(str(_x) for _x in self.genotype))
        def __eq__(self, other):
            return self.pos==other.pos and self.s==other.s and Counter(str(_x) for _x in self.genotype)==Counter(str(_x) for _x in other.genotype)
        def isMISSING(self):
            return set(self.genotype).issubset(set(['.','-','NA','na','Na']))
	def GetGenes(self):
		return self.genotype
	def GetPos(self):
		return self.pos
	def GetS(self):
		return self.s

def dropHomozygous(allgenolst):
	"""Eliminate homozygous genotypes from a list of genotype objects."""
	_s=-1 # genotypes will be numbered starting from 0 (_s+=1)
	sorted_genolst = sorted(allgenolst, key = lambda x: x.GetPos()) # sort the input list upon the SNP position
	het_genolist=[]
	for _geno in sorted_genolst:
		if len(set(_geno.GetGenes()))>=2:
			_s+=1
			het_genolist.append(Genotype(_s, _geno.GetPos(), *_geno.GetGenes()))
	return het_genolist

def getGenotypes(vcffile, Pop=False, contig=False):
	"""Extract the genotypes from a VCF file. Return the genotypes as a list of Genotypes objects."""
	_s=-1 # genotypes will be numbered starting from 0 (_s+=1)
	genolist=[]
        if contig:
                contig_names = [] # contig names in the vcf file, should be just one as otherwise haplotyping will be nonsense! 
	with open(vcffile, 'rU') as vfile:
		for _geno in vfile:
			if _geno[0]!='#':
				_genolst=_geno.split()
				if Pop or len(set(re.split('/|\|', _genolst[-1].split(':')[0])))>1:  # If the vcf file belongs to a member of a population, keep all of the SNP, even the homozygous ones.\ 
					_s+=1															 # Throw away homozygous SNPs otherwise.
					genolist.append(Genotype(_s, _genolst[1], *re.split('/|\|', _genolst[-1].split(':')[0]))) # Extract genotypes
                                if contig:
                                        contig_names.append(_genolst[0])
        if contig:
               return genolist, contig_names
	return genolist

def getGenotypesPop(vcffile, SampleNames, contig=False):
	"""Extract the genotypes from a multiple-sample VCF file of a population with given sample names. Return the genotypes as lists of Genotype objects for each sample in the population."""
	if not SampleNames:
		raise ValueError('No sample names were given to getAllelesPop to extract their alleles!')
	_s=-1  # genotypes will be numbered starting from 0 (_s+=1)
        SampleNames = [_sample for _sample in set(SampleNames)]
	_header = False
	genolist, SMnames, map_dict = [[] for _sm in SampleNames], [], dict()
        if contig:
                contig_names = [] # contig names in the vcf file, should be just one as otherwise haplotyping will be nonsense! 
	with open(vcffile, 'rU') as vfile:
		for _geno in vfile:
			if _geno[0]!='#':
				_genolst=_geno.split()
				_s+=1
				_allalleles_at_s = "" # Collect all of the alleles at s to throw away variants homozygous for all population members
				sample_alleles=[[] for _sample in SMnames]
				for _sample in SMnames:
					sample_alleles[map_dict[_sample]] = re.split('/|\|', _genolst[9:][map_dict[_sample]].split(':')[0]) # Extract the alleles present for each member of the population
					_allalleles_at_s+="".join(sample_alleles[map_dict[_sample]])
				if len(set(_allalleles_at_s)-set(["."]))>1: # Throw away variants homozygous in the population
					for _n, _sample in enumerate(SMnames):
						genolist[map_dict[_sample]].append(Genotype(_s, _genolst[1], *sample_alleles[map_dict[_sample]])) # Extract the genotypes of each member in the population
                                if contig:
                                        contig_names.append(_genolst[0])
			elif _geno[1]!='#':
				if set(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']).issubset(set(_geno.lstrip('#').split())):
					_header = True
					SMnames = _geno.split()[9:] # extract sample names in the VCF header
                                        if len(set(SMnames))<len(SMnames):
                                            raise ValueError('Duplicate sample names observed in the VCF header!')
					if set(SMnames).symmetric_difference(set(SampleNames)):
                                            raise ValueError('The given sample names do not match those in the VCF header!')
					map_dict = dict(zip(SMnames, [SampleNames.index(_sample) for _sample in SMnames]))
				else:
					pass
        if contig:
               return genolist, contig_names
	return genolist

def getAlleles(vcffile, Pop=False):
	"""Extract the alleles and their numerical code at (heterozygous) sites from a VCF file. Return\n\ 
the number, position and a dictionary, for the alleles and their coding as a list of triple tuples. If\n\
the VCF file belongs to a member of a population, keep all of the SNPs. Otherwise, throw away the homozygous ones."""
	_s=-1
	alllist=[]
	with open(vcffile, 'rU') as vfile:
		for _geno in vfile:
			if _geno[0]!='#':
				_alllst=_geno.split()
				if Pop or len(set(re.split('/|\|', _alllst[-1].split(':')[0])))>1:
					_s+=1
					alleles = [_alllst[3]]+re.split("\ *,\ *|,\ *|\ *,|,", _alllst[4])
					codes = sorted(set(re.split('/|\|', _alllst[-1].split(':')[0])), key = lambda x: int(x))
					if '0' not in codes: # the reference allele should have been coded as 0
						codes = ['0'] + codes # the reference and alternative nucleotides should have been CONSISTENTLY coded in all population files to allow merging later.
					alllist.append((_s, int(_alllst[1]), {_x.upper():int(_y) for _x, _y in zip(alleles, codes)}))
	return alllist

def getAllelesPop(vcffile, SampleNames):
	"""Extract the alleles and their numerical code from a multiple-sample VCF file of a population with given sample names. Return\n\ 
the SNP number, the SNP position and the dictionary that specifies the reference and the observed alternative alleles for the SNP,
as a list of tuples."""
	_s=-1
	alllist, SMnames = [], []
	_header = False
	map_dict = dict()
	if not SampleNames:
		raise ValueError('No sample names were given to getAllelesPop to extract their alleles!')
        SampleNames = [_sample for _sample in set(SampleNames)] # get read of duplicate sample names
        with open(vcffile, 'rU') as vfile:
		for _geno in vfile:
			if _geno[0]!='#':
				if not _header:
					raise ValueError('Unexpected header detected in the VCF file!')
				_alllst = _geno.split()
				_s+=1
				_allalleles_at_s = "" # Collect all of the alleles at s to throw away variants homozygous for all population members
				sample_alleles=[[] for _sample in SMnames]
				for _sample in SMnames:
					sample_alleles[map_dict[_sample]] = re.split('/|\|', _alllst[9:][map_dict[_sample]].split(':')[0]) # Collect the alleles for each member of the population
					_allalleles_at_s+="".join(sample_alleles[map_dict[_sample]])
				if len(set(_allalleles_at_s)-set(["."]))>1: # Throw away variants homozygous in the population 
					alleles = [_alllst[3]]+re.split("\ *,\ *|,\ *|\ *,|,", _alllst[4])
					codes=[]
					for _alleles in sample_alleles:
						codes+=_alleles
					codes = sorted(set(codes)-set(['.','-']), key = lambda x: int(x)) # '.', '-' corresponds to missing genotypes
					if '0' not in codes:
						codes = ['0'] + codes
					alllist.append((_s, int(_alllst[1]), {_x.upper():int(_y) for _x, _y in zip(alleles, codes)}))
			elif _geno[1]!='#':
				if set(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']).issubset(set(_geno.lstrip('#').split())):
					_header = True
					SMnames = _geno.split()[9:] # extract sample names in the VCF header
                                        if len(set(SMnames))<len(SMnames):
                                            raise ValueError('Duplicate sample names observed in the VCF header!')
					if set(SMnames).symmetric_difference(set(SampleNames)):
                                            raise ValueError('The given sample names do not match those in the VCF header!')
					map_dict = dict(zip(SMnames, [SampleNames.index(_sample) for _sample in SMnames]))
				else:
					pass
	return alllist

def mergeVCFs(reference, output, *vcffiles):
	"""Merge several vcf files in a single multi sample one, so that all of the lists have genotypes for the same positions\n\ 
and with the same numbering, starting with zero for the allele in the reference sequence at variant positions. It is assumed that\n\
both VCF and the reference contain only one and the same contig."""
	P = [4, 2, 3] # In case a VCF file is empty, the ploidy level cannot be estimated from it. This ad hoc solution should then be used!
	with open(reference,'rU') as handle, open(reference+".fai",'rU') as index:
		for record in SeqIO.parse(handle,"fasta"):
			_seq=record.__getattribute__('seq')
			_id=record.__getattribute__('id')
			_descrip=record.__getattribute__('description')
		line_fai=index.readline()
		columns=re.split('\s*|\t',line_fai.rstrip())
		contig_id= columns[0]                   # id derived from .fai
		contig_length= int(columns[1])          # seq length defined from .fai
		offset = int(columns[2])                # offset derived from .fai
		if contig_id!=_id:
			raise ValueError("Contig name in the fasta file does not match the name in its index!")
		if contig_length!=len(_seq):
			raise ValueError("Contig length in the fasta file does not match the length specified in its index!")
	SampleNames=[]
	for _file in vcffiles:
		with open(_file, 'rU') as handle:
			for _line in handle:
				if _line[0]=='#' and _line[1]!='#' and set(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']).issubset(set(_line.lstrip('#').split())):
					SampleNames.append(_line.rstrip().split()[-1])
					break
	genolists = [getGenotypes(_file, Pop=True) for _file in vcffiles]
	allelelists = [getAlleles(_file, Pop=True) for _file in vcffiles]
	for _i in range(0, len(genolists)): # filter genotypes homozygous with the reference allele out
		_homoz_ref = []
		for _n in range(0, len(genolists[_i])):
			if set(genolists[_i][_n].GetGenes())==set(['0']):
				_homoz_ref.append(_n)
		for _n in sorted(_homoz_ref, reverse =  True):
			del genolists[_i][_n]
			del allelelists[_i][_n]
	pos=set([])
	missed_pos = []
	ploidies = []
	ngenolists=[[] for _i in range(0, len(genolists))] # genotype codes converted to actual nucleotides
	for _i, _genolst in enumerate(genolists):
		missed_pos.append(set(_G.GetPos() for _G in _genolst))
		try:
			ploidies.append(len(_genolst[0].GetGenes()))
		except IndexError:
			garbage = sys.stderr.write("WARNING: No variant was detected in {0:s}!\n".format(vcffiles[_i]))
			ploidies.append(P[_i])
		pos = pos.union(missed_pos[-1])
		ngenolists[_i]=[Genotype(_G.GetS(), _G.GetPos(), *[list(_k for _k in _alleles[-1].keys())[list(str(_v) for _v in _alleles[-1].values()).index(str(_code))] for _code in _G.GetGenes()]) for _alleles, _G in zip(allelelists[_i], _genolst)]
	for _n, _missed_pos in enumerate(missed_pos):
		missed_pos[_n] = pos.difference(_missed_pos)
	pos = sorted(pos, key = lambda x: int(x)) # sorted list of all of the SNP positions in the genotype lists
	reference_alleles = [_seq[_pos-1] for _pos in pos]
	unified_alleles=[]
	for _s, _pos in enumerate(pos):
		_unified_alleles=set([])
		for _allelelst in allelelists:
			for _alleles in _allelelst:
				if _alleles[1]<_pos:
					pass
				elif _alleles[1]==_pos:
					_unified_alleles=_unified_alleles.union(set(_alleles[-1].keys()))
					_allelelst.remove(_alleles)
					break
				else:
					break
		_unified_alleles = _unified_alleles-{reference_alleles[_s]} # drop the reference allele as it will be separately coded as 0
		_dic = {reference_alleles[_s]:0} # code the reference allele as 0
		_dic.update({_nuc:_code for _nuc, _code in zip(sorted(_unified_alleles), [_code for _code in range(1, len(_unified_alleles)+1)])}) # add the coding for the alternative alleles
		unified_alleles.append((_s, _pos, _dic))
	updated_genolists = []
	for _n, _genolst in enumerate(ngenolists): # new SNP numbers correspond to the unified position list. Missing genotypes of each list are assumed homozygous-reference.
		updated_genolists.append(sorted([Genotype(pos.index(_G.GetPos()), _G.GetPos(), *[unified_alleles[pos.index(_G.GetPos())][-1][_nuc] for _nuc in _G.GetGenes()]) for _G in _genolst]+
									   [Genotype(pos.index(_missed_pos), _missed_pos, *[0 for _x in range(0, ploidies[_n])]) for _missed_pos in missed_pos[_n]],
										key = lambda x: x.GetPos()))
	with open(output, 'w') as mergedvcf:
		mergedvcf.write('#'+'\t'.join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])+'\t'+'\t'.join(SampleNames)+'\n')
		for _s, _genos in enumerate(zip(*updated_genolists)):
			mergedvcf.write('.'+'\t'+str(pos[_s])+'\t'+'.'+'\t'+reference_alleles[_s]+'\t'+','.join(_allele for _allele in unified_alleles[_s][-1].keys() if _allele!=reference_alleles[_s])+'\t'+'\t'.join('.' for _x in ['QUAL', 'FILTER', 'INFO', 'FORMAT'])+'\t'+
				'\t'.join('/'.join(str(_code) for _code in _g.GetGenes()) for _g in _genos)+'\n')					
	return None
#a = Genotype(0,100,*[1,1,1,0])
#b = Genotype(1,100,*[1,0,1,0])
#c = Genotype(0,90,*[1,1,1,1])
#lst1=[[a],[b,c]]
#mergeGenotypes(*lst1)
#getGenotypesPop("example.trio.vcf", ['parent1','parent2','progeny'])
#getAllelesPop("example.trio.vcf", ['parent1','parent2','progeny'])

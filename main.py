#!/usr/bin/env python
# Haplotype estimation method for polyploid F1 populations (Motazedi at el. 2017).
# Written by Ehsan Motazedi, Wageningen UR, 09-11-2016.
# Last Updated: 28-03-2017

import bamprocess
import copy
import functools
import getopt
import os
import random as rnd
import sys
import textwrap
import traceback
from math import exp
from bamprocess import InputError
from branchprune import BranchPop, PrunePop, BlockException, makePermutation, SetGenos
from genotypes import Genotype, getGenotypesPop, dropHomozygous
from haplotypes import Haplotypes, getMEC
from reads import Extract_Components, getReads, SubReads, SemiRead

if __name__=='__main__':
	try:
		top = False
		_max_len = len('Read/Fragment file'+'[0,1]') # used in the formatting of the help message
		hlpmsg = (['']+textwrap.wrap("Haplotype estimation tool for polyploid populations using Next Generation Sequencing, based on a Bayesian framework incorporating inheritance pattern and recombination (Motazedi et al. 2017).", width = 95)+
		['Written by Ehsan Motazedi, Wageningen UR, 21-03-2017.']+
		['\nPositional Arguments:']+
		[' '*4+'Read/Fragment file'+' '*4+' STR '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The name of the multi-sample BAM file containing the input aligned reads of the members of the population. The first two sample names in the BAM file are considered the names of the maternal and parental samples, respectively, and the rest samples are considered as progeny.', width=64))]+
		[' '*4+'VCF file          '+' '*4+' STR '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The name of the multi-sample VCF file specifying the variants in the population, corresponding to the input BAM file.', width=64))]+
		[' '*4+'Output            '+' '*4+' STR '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The path name of the output directory to write the results, as Output/TriPolySolution_{samplename} and Output/MEC_{samplename}', width=64))]+
		['\nOptional Arguments:']+
		[' '*4+'-a, --aggressive  '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('With this option, aggressive pruning will be performed with the given rate, '+u"\u03B1"+', after each extension if the original pruning scheme is not able to rule-out at least (100'+u"\u00D7"+
u"\u03B1"+')% of the candidate haplotypes, so that at most 1 + 100'+u"\u00D7"+'(1-'+u"\u03B1"+')% of the haplotypes remain after pruning (default is the original scheme using maximum RL).', width=64))]+
		[' '*4+'-e, --error       '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The fixed base calling error rate in the read-fragments, between 0 and 1. If not specified, the software will first try to use the bam quality scores for each position within each fragment. (default = 0.015).', width=64))]+
		[' '*4+'-c, --recom       '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap("The recombination rate between 0 and 0.5, determined for example using Haldane's formula, converted to a bi-nomial recombination probability for each SNP on each homologue (default = 0.0001).", width=64))]+
		[' '*4+'-r, --rho         '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The branching threshold of the algorithm, '+u'\u03C1'+', as explained in Motazedi et al. (2017) (default = 0.3).', width=64))]+
		[' '*4+'-k, --kappa       '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The pruning rate of the algorithm, '+ u'\u03BA'+', using the relative likelihood of the haplotypes as explained in Motazedi et al. (2017). (default = 0.7).', width=64))]+
		[' '*4+'-w, --warmup      '+' '*4+' INT '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('The genomic length in base-pairs from the beginning of each haplotype block, along which no candidate extension is pruned or excluded during branching. This initial warm-up phase increases the precision and efficiency of the pruning and branching steps at the cost of memeory/speed (default = 0).', width=64))]+
		[' '*4+'--mmq             '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('Minimum read mapping quality to consider a read for phasing (default 20).', width=64))]+
		[' '*4+'--mbq             '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('Minimum phred score to consider a base for haplotype fragment (default 13).', width=64))]+
		[' '*4+'--maxIS           '+' '*4+'[0,1]'+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('Maximum insert-size for a paired-end read to be considered as a single fragment for phasing (default 3000).', width=64))]+
		[' '*4+'--nogeno          '+' '*4+'     '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('By default, parental haplotype estimates are restricted to those compatible with the detected genotypes. The progeny estimates are also compatible with the progeny genotypes unless Mendelian inheritance is violated in the current parents/progeny genotypes. This options relaxes this restriction and all dosages of the detected alleles are considred at each SNP position for each individual.', 
		width=64))]+
		[' '*4+'--filter          '+' '*4+'     '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('With this option, the phasing is only reported for the heterozygous SNPs for each individual, i.e. the homozygous SNPs of each individual are filtered out of its haplotype estimate (default is to report all SNPs).', width=64))]+
		[' '*4+'-t, --top         '+' '*4+'     '+' '*4+('\n'+' '*(3*4+_max_len)).join(textwrap.wrap('With this option, only the most likely haplotype is reported from the final set of haplotypes that have survived the prunings (default is to report all).', width=64))]+
		['\n'])
		optlist, args = getopt.gnu_getopt(sys.argv[1:], 'he:k:r:a:w:tc:', ["help", "error=", "kappa=", "rho=", "aggressive=", "warmup=", "top", "recom=", "mmq=", "mbq=", "maxIS=", "nogeno","filter"])
		valid_options = set([("-a", "--aggressive"),("-h", "--help"), ("-w", "--warmup"), ("--error", '-e'), ("--kappa", '-k'), ("--rho", '-r'), ("-t","--top"), ("--mmq", ), ("--mbq", ), ("--maxIS", ), ("-c", "-recom"), ("--nogeno", ), ("--filter", )])	
		valid_set = set(_y for _x in valid_options for _y in _x)
		if set(_v[0] for _v in optlist).issubset(valid_set):
			if set(["-h", "--help"]) & set(_v[0] for _v in optlist):
				try:
					garbage = sys.stdout.write('{}\n'.format('\n'.join(hlpmsg)))
				except UnicodeEncodeError as e:
					import codecs
					UTF8Writer = codecs.getwriter('utf8')
					sys.stdout = UTF8Writer(sys.stdout)
					for _x in hlpmsg:
						garbage = sys.stdout.write(_x+'\n')
				finally:
					sys.exit(0)
			valid_options.remove(("-h", "--help"))
		else:
			raise InputError('Unrecognized optional arguments!')
		if len(args)<1:
			raise InputError('No BAM and VCF file given!')
		if len(args)<2:
			raise InputError('No VCF file given!')
		if len(args)<3:
			raise InputError('Output file name not specified!')
		if len(args)>3:
			raise InputError('Unrecognized positional arguments found! Check the input to TriPoly! Detected positional arguments:\n'+'\t'+'\t'.join(args))
		if set(["-t", "--top"]) & set(_v[0] for _v in optlist):
			top = True
		valid_options.remove(("-t", "--top"))
		if set(["--nogeno"]) & set(_v[0] for _v in optlist):
			GenoConstraint = False
			sys.stderr.write("WARNING: all of the dosages will be estimated anew! The dosages given in the input VCF file will be ignored!\n")
		else:
			GenoConstraint = True
		valid_options.remove(("--nogeno", ))
		if set(["--filter"]) & set(_v[0] for _v in optlist):
			filter_on = True
		else:
			filter_on = False
		valid_options.remove(("--filter", ))
		valid_options = sorted(valid_options, key = lambda x: x[0].lstrip('-'))
		opts = [None, 0.0001, '1.5e-02', 0.7, 3000, 13, 20, 3e-01, 0]	# default values for aggressive pruning rate, recombination rate, base calling error, pruning rate, maxIS, minimum base quality, minimum mapping quality, branching threshold and warm-up length.
		for _n, _pair in enumerate(valid_options):
			_lst = list(_v[1] for _v in optlist if _v[0] in _pair)
			if len(_lst)>1:
				raise ValueError('Duplicate values for {0}/{1}!'.format(*_pair))
			try:
				opts[_n] = float(_lst[0]) if _n<4 or _n==7 else int(_lst[0])
			except IndexError:
				pass
			except ValueError as e:
				if "could not convert string to float" in e.args[0] or "invalid literal for int" in e.args[0] or "invalid literal for float" in e.args[0]:
					e.args=("Invalid value '"+_lst[0]+"' passed to "+','.join(_pair)+"! Use -h, --help for help!",)+e.args[1:]
					raise
				raise
		alpha, recombination_rate, error_rate, k, maxIS, mbq, mmq, rho, warmup = opts
		if alpha and (alpha>1 or alpha<0):
			raise ValueError('The aggressive pruning rate, '+u"\u03B1"+', must be between 0 and 1!')
		if isinstance(error_rate, str): # happens if no fixed error rate is specified
			error_fixed = False
			error_rate = float(error_rate)
		elif error_rate>1 or error_rate<0:
			raise ValueError('The variant calling error must be between 0 and 1!')
		else:
			error_fixed = True
		if recombination_rate>0.5 or recombination_rate<0:
			raise ValueError('The recombination rate must be between 0 and 0.5!')
		if k>1 or k<0:
			raise ValueError('The pruning rate, '+u'\u03BA'+', must be between 0 and 1!')
		if rho>1 or rho<0:
			raise ValueError('The branching threshold, '+u'\u03C1'+', must be between 0 and 1!')
		if maxIS < 0:
			raise ValueError('The maximum insert-size, maxIS, cannot be negative!')
		if warmup<0:
			raise ValueError('The warm-up length cannot be negative!')	
		if not all(isinstance(_x, str) for _x in args):
			raise InputError("The input/output file names must be valid strings!")
		try:
			Frag_lsts, Qual_lsts, SampleNames = bamprocess.get_frags(args[0], args[1], maxIS, mbq, mmq)  # Extract a separate SNP-fragment/quality score list for each member of the population
			#for _i, Frags in enumerate(Frag_lsts):
			#	garbage= sys.stdout.write("\n************** Fragments Individual {0}:\n".format(SampleNames[_i])+';;'.join(str(x.GetDict()) for x in Frags)+"\n################## Quality Score Individual {0}:\n".format(SampleNames[_i])+';;'.join(str(y) for y in Qual_lsts[_i]))
			if error_fixed: # do not use the quality scores if a fixed base calling error is given by the user
				Qual_lsts = [None for _i in range(0, len(Frag_lsts))]
		except IOError:
			raise InputError('The input BAM file was not found!')
		except:
			raise
		try:
			GenosORIGINAL, contig_names = getGenotypesPop(args[1], SampleNames, True) # extract the genotypes of each population from the VCF file
                        if len(set(contig_names))>1:
                                raise BlockException("The VCF file contains more than one contig! Please run TriPoly separately for each contig!")
                        if set(contig_names).issubset(set(['.',''])):
                                sys.stderr.write("WARNING: No contig name is specified in the VCF file! All of the variants are assumed located on the same contig!\n")
			ploidies_m_f = [] # ploidy-levels of the parents to be extracted from the vcf file                              
			for _GENOS in GenosORIGINAL[0:2]: # determine the ploidy of the parents
				for _geno in _GENOS:
					if set(['-']).intersection(set(_geno.GetGenes())) or len(_geno.GetGenes())<2:
						continue
					ploidies_m_f.append(len(_geno.GetGenes()))
					break
			if len(ploidies_m_f)!=2:
				raise ValueError("Could not detect ploidy levels of the parents from the VCF file or ploidy is less than 2 for some! Check the VCF file!")     
			if not GenoConstraint:
				garbage = SetGenos(args[1], SampleNames, *ploidies_m_f)
		except IOError:
			raise InputError('The VCF file was not found!')
		except:
			raise
		if [] in GenosORIGINAL:
			sys.exit('Program terminates as there is no genotypes to phase!')
		GenosORIGINAL_backup = copy.deepcopy(GenosORIGINAL) # make a backup of the list of the genotypes to be used in reporting of the output, as the latter gets modified later on
		BLOCKS, MEC_Scores, first_block = [], [], True
		Frags =[] # the set of all fragments in the population
		for _frags in Frag_lsts: # aggregate all of the SNP-fragments to build the SNP-fragment matrix from. As pedigree info is used, a SNP is considered connected\
			Frags+=_frags    # if it contained in at least one informative fragment, which can originate from either parent or the progeny.
		MECwritten = False
		for comp in Extract_Components(Frags):	 # phase each connected-component separately
			sys.stderr.write("{0}aplotype block started at SNP {1:d}...\n".format("H" if first_block else "Next h", min(comp)+1))
			first_block = False
			#garbage=print("COMP:",comp)
			Genos = []
			for _GenosORIGINAL in GenosORIGINAL:
				Genos.append(sorted((_geno for _geno in _GenosORIGINAL if _geno.GetS() in comp), key = lambda x: x.GetS())) # extract the genotypes that belong to each component and\
			#print(Genos)
			GenosORIGINAL = [list(set(_GenosORIGINAL)-set(_Genos)) for _GenosORIGINAL, _Genos in zip(GenosORIGINAL, Genos)] # eliminate those genotypes from the total list
			H_pruned_pop = [tuple(Haplotypes(_Genos[0].GetS(), _Genos[0].GetS(), 0, makePermutation(_Genos[0], lognum=True), None, None, *_Genos[0].GetGenes()) for _Genos in Genos)] # Start the haplotypes with the first SNP of the block
			for _id, _member_haplotype in enumerate(H_pruned_pop[0]):
				if '.' in _member_haplotype.GetVS():
					_member_haplotype.ChangeVS(*['-' for _x in range(0, max(len(_Geno.GetGenes()) for _Geno in Genos[_id]))])
			MEC = [0]
			_s, block_start = (0, 0)
			while _s < (len(Genos[0])-1):
				_s+=1
				SemiFrags = []
				
				SemiFrags = tuple(map(lambda x: SemiRead(x, Genos[0][_s].GetS()), _Frags) for _Frags in Frag_lsts) # Extract the semi-reads for position _s from the total list of the fragments
				try:
					H_pop_ext = [_H_pop_ext for _H_base_pop in H_pruned_pop for _H_pop_ext in BranchPop(_H_base_pop, tuple(_Genos[_s] for _Genos in Genos), SemiFrags, rho if Genos[0][_s].GetPos()-Genos[0][block_start].GetPos()>=warmup else 0, error_rate, recombination_rate, Qual_lsts, GenoConstraint, True)]
					subfrags_current = list(SubReads(_Frags, set(Genos[0][_x].GetS() for _x in range(block_start, _s+1))) for _Frags in Frag_lsts)
					H_pruned_pop, MEC = PrunePop(list(_H_pop_ext for _H_pop_ext in H_pop_ext), k if Genos[0][_s].GetPos()-Genos[0][block_start].GetPos()>=warmup else 0, error_rate, alpha if Genos[0][_s].GetPos()-Genos[0][block_start].GetPos()>=warmup else None, subfrags_current, Qual_lsts, filter_on)
				except BlockException as e:
					garbage = sys.stderr.write("WARNING: New Block started at SNP number {:d} due to extension error:\n{}\n".format(Genos[0][_s].GetS()+1, e))
					BLOCKS.append(H_pruned_pop)
					MEC_Scores.append(MEC)
					#MEC_Scores.append([tuple(getMEC(_memreads,_H_pruned_pop) for _memreads, _H_pruned_pop in zip(subfrags_current, H_pruned_pop[_solnum])) for _solnum in range(0, len(H_pruned_pop))])
					block_start = _s
					H_pruned_pop, MEC = [tuple(Haplotypes(_Genos[block_start].GetS(), _Genos[block_start].GetS(), 0, makePermutation(_Genos[block_start], lognum=True), None, None, *_Genos[block_start].GetGenes()) for _Genos in Genos)], [0]
			BLOCKS.append(H_pruned_pop)
			MEC_Scores.append(MEC)
			#MEC_Scores.append([tuple(getMEC(_memreads,_H_pruned_pop) for _memreads, _H_pruned_pop in zip(subfrags_current, H_pruned_pop[_solnum])) for _solnum in range(0, len(H_pruned_pop))])
		if filter_on and GenoConstraint:
			Genos = [dropHomozygous(_Genolst) for _Genolst in GenosORIGINAL_backup] # Restore the original genomes and throw away homozygous variants for each population member
		else:
			Genos = [_Genolst for _Genolst in GenosORIGINAL_backup] # Restore the original genomes
		OriginalS = [[] for _id in range(0, len(Genos))] # Store the original SNP numbers to be written to the final output
		for _id in range(0, len(Genos)):
			for _Geno in Genos[_id]:
				for _Org in GenosORIGINAL_backup[_id]:
					if _Org.GetPos()==_Geno.GetPos():
						OriginalS[_id]+=[_Org.GetS()]
						break
		if top:
			#BLOCKS_TOPS_index = [[(_X_pruned_population, _indx) for _indx, _X_pruned_population in enumerate(_H_pruned_population) if _X_pruned_population[0].GetRL()==_H_pruned_population[0][0].GetRL() and MEC_Scores[_blck][_indx]==MEC_Scores[_blck][0]] for _blck, _H_pruned_population in enumerate(BLOCKS)] # keep the top solutions and throw-away the others
			BLOCKS_TOPS_index = [[(_X_pruned_population, _indx) for _indx, _X_pruned_population in enumerate(_H_pruned_population) if _X_pruned_population[0].GetRL()==_H_pruned_population[0][0].GetRL() and sum(MEC_Scores[_blck][_indx])==sum(MEC_Scores[_blck][0])] for _blck, _H_pruned_population in enumerate(BLOCKS)] # keep the top solutions and throw-away the others
			BLOCKS_index = [[_H_pruned_population[rnd.choice(range(0, len(_H_pruned_population)))]] for _H_pruned_population in BLOCKS_TOPS_index] # randomly choose one of the top solutions as the final one
			BLOCKS = [[_bi[0] for _bi in bi] for bi in BLOCKS_index]
			MEC_Scores = [[MEC_Scores[num][_bi[1]] for _bi in bi] for num, bi in enumerate(BLOCKS_index)]
		_dir = True
		try:
			os.mkdir(args[2])
		except OSError as e:
			from errno import EEXIST as DIR_EXISTS # IF the dir already exists, just write the files in it! Otherwise, raise the error.
			if e.errno != DIR_EXISTS:
				_dir = False
				e.args = e.args[:-1]+(e.args[-1]+" '"+args[2]+"'",)
				e.args+=('Could not write the results to or make the specified output directory!',)
				raise OSError(*e.args)
		except:
			_dir = False
			raise
		finally:
			garbage = _dir and os.chdir(args[2]) # Only try to change the directory if the directory of destination exists

		BLOCKS_NEW = [[] for _sample in SampleNames]  # Transpose BLOCKS for easier processing
		for _H_pruned_block in BLOCKS:
			for _x in BLOCKS_NEW:
				_x.append([]) # first/next block for all of the population members  
			for _H_pruned_population in _H_pruned_block:
				for _id, _H_pruned in enumerate(_H_pruned_population):
					BLOCKS_NEW[_id][-1].append(_H_pruned)
		for _id in range(0, len(SampleNames)):
			blcknum = 0
			BLOCKS = BLOCKS_NEW[_id]
			line_start = 1
			with open('TriPolySolution_'+SampleNames[_id], 'w') as outhap, open('MEC_'+SampleNames[_id], 'w') as outMEC:
				for _blck, _H_pruned in enumerate(BLOCKS):
					blcknum+=1
					solution_number = -1
					blckrng = []
					n_snps_blck = []
					all_Invalid_Genos = []
					dropped_solutions = 0
					for _HH_pruned in _H_pruned:
						solution_number += 1
						n_snps_blck.append(0)
						_rng_start = [_n for _n, _s in enumerate(OriginalS[_id]) if _s >_HH_pruned.GetStart()-1]
						_rng_end = [_n for _n, _s in enumerate(OriginalS[_id]) if _s >_HH_pruned.GetStop()]
						if len(_rng_start)>0:
							_rng_start = _rng_start[0]
						else:
							_rng_start = -1
						if len(_rng_end)>0:
							_rng_end = _rng_end[0]
						else:
							_rng_end = len(OriginalS[_id])
						if _rng_start >= 0:
							rng = [_num for _num in range(_rng_start, _rng_end)]
						else:
							rng = []
						blckrng.append(rng)
						if filter_on: # filter out homozygous SNPs
							Invalid_Genos = [_s for _s in range(_HH_pruned.GetStart(), _HH_pruned.GetStop()+1) if '-' in _HH_pruned.GetGenotype(_s) or len(set(_HH_pruned.GetGenotype(_s)))<2] # throw away unphased as well as homozygous SNPs
							Invalid_Genos = [_s for _s in Invalid_Genos if _s in [OriginalS[_id][_n] for _n in rng]]
						else:
							Invalid_Genos = []
						all_Invalid_Genos.append(Invalid_Genos)
						n_snps_blck[-1] = len(rng)-len(Invalid_Genos)
						if len(rng)-len(Invalid_Genos)<2:
							dropped_solutions +=1
							continue
						if (solution_number == 0):
							garbage = outhap.write('-'*19+' BLOCK '+str(blcknum)+' '+'-'*19+'\n')
						garbage = outhap.write("{3}\t{0}\t{1}\t{2}\n".format(Genos[_id][rng[0]].GetPos(),
                                                                                                line_start,
                                                                                                line_start+n_snps_blck[-1]-1, 
											"Top Solution" if top else "Solution "+str(solution_number-dropped_solutions+1)))
						for _n in rng:
							if OriginalS[_id][_n] not in Invalid_Genos:
								_new_line = '{2}\t{0:d}\t{1}\n'.format(Genos[_id][_n].GetPos(), '\t'.join(str(_x[OriginalS[_id][_n]-_HH_pruned.GetStart()]) for _x in _HH_pruned.GetVS()), contig_names[0])
								garbage = outhap.write(_new_line)
					if all((len(_x)-len(_Invalid_Genos))<2 for _x, _Invalid_Genos in zip(blckrng, all_Invalid_Genos)):
						line_start += max(_n_snps_blck for _n_snps_blck in n_snps_blck)
						continue
					start_snp_pos = min(Genos[_id][blckrng[_x][0]].GetPos() for _x in range(0, len(blckrng)))
					stop_snp_pos = max(Genos[_id][blckrng[_x][-1]].GetPos() for _x in range(0, len(blckrng)))
					garbage = outMEC.write('-'*19+' BLOCK '+str(blcknum)+' '+'-'*19+'\n')
					garbage = outMEC.write('{0:d} SNPs over a genomic length of {1:d} nucleotides.\n'.format(max(_n_snps_blck for _n_snps_blck in n_snps_blck), stop_snp_pos-start_snp_pos+1))
					garbage = outMEC.write('Start SNP number, coordinate: {0:d}, {1:d}(bp)\tStop SNP number, coordinate:{2:d}, {3:d}(bp)\n'.format(line_start, start_snp_pos, line_start+max(_n_snps_blck for _n_snps_blck in n_snps_blck)-1, stop_snp_pos))
					_solution_num = 0
					line_start += max(_n_snps_blck for _n_snps_blck in n_snps_blck)
					for _k, _rng in enumerate(blckrng):
						if len(_rng)-len(all_Invalid_Genos[_k])<2:
							continue
						_solution_num+=1
						#garbage = outMEC.write("MEC for solution {0:d} (including all population reads and haplotypes):\t{1:3d}\tRL for solution {0:d} (accounting for all population reads and haplotypes):\t{2:5.3f}\n".format(_solution_num, MEC_Scores[_blck][_k], exp(_H_pruned[_k].GetRL())))		
						garbage = outMEC.write("MEC for {0}:\t{1:3d}\tRL of {0} (accounting for all population reads and haplotypes):\t{2:5.3f}\n".format("the top solution" if top else "solution "+str(_solution_num), MEC_Scores[_blck][_k][_id], exp(_H_pruned[_k].GetRL())))
	except getopt.GetoptError as e:
		garbage = sys.stderr.write("GetoptError: {0}\n".format(str(e)))
		sys.exit(2)
	except OSError as e:
		garbage = sys.stderr.write("OSError: {0}\n".format(str(e)))
		sys.exit(4)
	except SystemExit as e:
		sys.exit(e.args[0])
	except (TypeError, ValueError) as e:
		try:
			garbage = sys.stderr.write(repr(e).split('(')[0]+': '+e.args[0]+'\n')
		except UnicodeEncodeError:
			import unicodedata
			e.args = tuple(unicodedata.normalize('NFKD', _arg).encode('ascii','ignore') for _arg in e.args) # drop the unicode characters
			garbage = sys.stderr.write(repr(e).split('(')[0]+': '+' '.join(_x.strip() for _x in e.args[0].split(',') if _x!=' ')+'\n')	# eliminate the empty place of the dropped characters
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_tb(exc_traceback, limit=1, file=sys.stderr)
		traceback.print_exception(exc_type, exc_value, exc_traceback, limit = 20, file=sys.stderr)
		sys.exit(2)
	except InputError as e:
		garbage = sys.stderr.write(str(e))
		#exc_type, exc_value, exc_traceback = sys.exc_info()
		#traceback.print_tb(exc_traceback, limit=1, file=sys.stderr)
		#traceback.print_exception(exc_type, exc_value, exc_traceback, limit = 20, file=sys.stderr)
		sys.exit(2)
	except BlockException as e:
		garbage = sys.stderr.write("BlockException: {}\n".format(e))
	except:
		garbage = sys.stderr.write('Unexpected error:\n')
		exc_type, exc_value, exc_traceback = sys.exc_info()
		traceback.print_tb(exc_traceback, limit=1, file=sys.stderr)
		traceback.print_exception(exc_type, exc_value, exc_traceback, limit = 4, file=sys.stderr)
		sys.exit(3)

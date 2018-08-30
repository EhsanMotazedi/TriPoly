# Make the SNP-fragment list from a multi-sample bam file and its corresponding VCF file  with input options: mmq, mbq, maxIS and qoffset.
# Written by Ehsan Motazedi, Wageningen UR, 04-11-2016.
# Last mofified: 22-02-2018.

import copy
import pysam
import re
import subprocess
from genotypes import getAllelesPop, getGenotypesPop
from math import log
from reads import Read

def adjust_seq(read):
	""" Adjust the read sequence according to the mapped positions, i.e. get read of the insertions and clipped bases."""
	cig = list(_x for _x in re.split('[0-9]{1,}', read.cigarstring) if _x)
	cign = list(_x for _x in re.split('[^0-9]', read.cigarstring) if _x)
	cig = list(cig[_n]*int(cign[_n]) for _n in range(0, len(cig)))
	cig = ''.join(_x for _x in cig if 'D' not in _x) # deleted nucleotides from the reference are not present in the read sequence
        adj_seq = []
        adj_qual = []
        for _n, _x in enumerate(read.seq):
            if cig[_n]=='M': # Allow only match/mismatch of the read nucleotides, i.e. no clipping, no insertion.
                adj_qual.append(read.qual[_n])
                adj_seq.append(_x)
        return ''.join(adj_seq), adj_qual

def frag_gen(varpos, allelelst, genolst, coordinates, nucleotides, qscores):
	""" Generate SNP-fragments from (paired-end) reads.""" 
	if (not coordinates) or (not nucleotides):
		return Read(), {}
	var_codes = []
	var_num = []
	var_q = []
	for _n, _x in enumerate(varpos):
		if _x < coordinates[0]:
			continue
		if _x > coordinates[-1]:
			break
		if _x in coordinates:
			if set(genolst[_n].GetGenes()).intersection(set(['.','-'])): # Throw away missing genotypes or genotypes with one or more missing allele(s)
				continue
			if len(set(genolst[_n].GetGenes()))<2: # Do not include in the SNP-fragments belonging to a population member its homozygous alleles
				continue
			try:
				var_codes.append(allelelst[_n][2][nucleotides[coordinates.index(_x)]])
				var_num.append(allelelst[_n][0])
				var_q.append(qscores[coordinates.index(_x)])
			except KeyError: # if the called nucleotide is wrong, i.e. does not exist in VCF alleles
				pass
	try: # return the reads {SNP number: allele} and the associated quality scores {SNP number: Qscore}
		return Read({_x:str(_y) for _x, _y in zip([var_num[0], var_num[1]]+var_num[2:], var_codes)}), {_x:str(_y) for _x, _y in zip(var_num, var_q)}
	except IndexError: # throw away reads with less than 2 SNPs
		return Read(), {}

class InputError(Exception):
	""" Handle invalid input specifications.""" 
	def __init__(self, msg, *args):
		super(InputError, self).__init__(args)
		self.msg = msg
	def __str__(self):
		return "InputError: {}\n".format(self.msg)
	def __repr__(self):
		return (self.msg,)+self.args+('\n',)

def get_frags(bamfile, vcffile, maxIS=3000, mbq=13, mmq=20, qoffset=33):
	""" 
mmq : minimum read mapping quality to consider a read for phasing, default 20\n
qvoffset <33/64> : quality value offset, 33/64 depending on how quality values were encoded, default is 33\n
mbq  : minimum base quality to consider a base for haplotype fragment, default 13\n
maxIS : maximum insert-size for a paired-end read to be considered as a single fragment for phasing, default 3000.\n 
	"""
	try:
		all_reads = pysam.Samfile(bamfile, 'rb')
	except IOError:
		raise InputError('The input BAM file was not found!')
	ReadHeader = subprocess.Popen(["samtools","view","-H", bamfile], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
	header, err_header = ReadHeader.communicate()
	if ReadHeader.returncode!=0:
		raise InputError('Failed to read the header from the bam file! Original error message:\n'+err_header)
	if isinstance(header, bytes):
		header=bytes.decode(header)
	else:
		pass
	RGIDs, SMnames = [], []
	for _headerline in header.splitlines(): # pasre the header of the bam file to extract the ID and SM fields of each Read Group
		if _headerline[0:3]=='@RG':
			RGID_added, SM_added = False, False
			for _n, _field in enumerate(_headerline.split()):
				if 'ID' in _field:
					if not RGID_added:
						RGIDs.append(''.join(_headerline.split()[_n].split(':')[1:])) # add read group ID
						RGID_added = True
					else:
						raise InputError('Double ID fields detected in @RG header line!')
				elif 'SM' in _field:
					if not SM_added:
						SMnames.append(''.join(_headerline.split()[_n].split(':')[1:])) # add the sample name associated with the read group ID
						SM_added = True
					else:
						raise InputError('Double SM fields detected in @RG header line!')
			if SM_added and RGID_added:
				pass
			elif SM_added:
				raise InputError('ID field missing in @RG header line!')
			elif RGID_added:
				raise InputError('SM field missing in @RG header line!')
			else:
				raise InputError('ID and SM fields missing in @RG header line!')
	if len(RGIDs)!=len(set(RGIDs)):
		raise InputError('Duplicate read group IDs detected in the bam header!')
	GroupedReadsWithID = [[] for _id in RGIDs] # separate reads belonging to each Read Group
	for _read in all_reads:
		GroupedReadsWithID[RGIDs.index(dict(_read.get_tags())['RG'])].append(_read)
	GroupedReads, GroupedSM = [], [] # combine the reads with different RGID but the same SM as they are assumed to belong to the same sample
	for _SMname, _ReadGroup in zip(SMnames, GroupedReadsWithID):
		if _SMname not in GroupedSM:
			GroupedReads.append(_ReadGroup)
			GroupedSM.append(_SMname)
		else:
			GroupedReads[GroupedSM.index(_SMname)]+=_ReadGroup
	del GroupedReadsWithID
	try:
		genolst = getGenotypesPop(vcffile, GroupedSM)
		allelelst = getAllelesPop(vcffile, GroupedSM)
	except IOError:
		raise InputError('The VCF file was not found!')
	except:
		raise
	Frag_lsts, Q_lsts = [], []
	varpos = list(_x[1] for _x in allelelst)
	for _group, reads in enumerate(GroupedReads):
		frag_lst, q_lst = [], []
		reads = sorted(reads, key= lambda _x: (_x.qname, _x.flag & 0x900)) # sort the alignments using their names, with the primary alignments being placed first.
		_rnum = 0
		rNUM = len(reads)
		while _rnum < rNUM: # scan through the alignments to find the pairs/singles 
			break_mate = False
			is_proper_pair = False
			read = copy.deepcopy(reads[_rnum])
			if read.is_unmapped or read.is_duplicate or (read.flag & 0x900): # throw away unmapped reads, duplicates and secondary/supplementary alignments
				_rnum+=1
				continue
			try:
				if read.qname == reads[_rnum+1].qname: # means the read is paired to a mate or has multiple/supplemenatry alignments
					if reads[_rnum+1].is_unmapped or reads[_rnum+1].is_duplicate or (reads[_rnum+1].flag & 0x900): # if the next read is unmapped, a duplicate or not primarym: skip it
						pass
					else:
						is_proper_pair = True	# means the read is paired to a proper mate
						mate = copy.deepcopy(reads[_rnum+1])
					_rnum+=2
				else: # means the read is single
					_rnum+=1
			except IndexError: # could occur for the last read in the alignments' list
				_rnum+=1
			if is_proper_pair:
				if (max(mate.positions+read.positions)-min(mate.positions+read.positions)+1)>maxIS: # Check the maximum insert-size to consider the mates as a single fragment
					break_mate = True
				if read.mapping_quality >= mmq:
					try:
                                                adj_seq, adj_qual =  adjust_seq(read)
                                                coordinates, nucleotides, quals = list(zip(*[(int(_x), _y, _z) for _x, _y, _z in zip(read.positions, adj_seq.upper(), list(ord(_x)-qoffset for _x in adj_qual)) if _z>=mbq]))
					except ValueError as e:
						if e.args[0][0:len("need more than 0 values to unpack")]=="need more than 0 values to unpack" or e.args[0][0:len("not enough values to unpack")]=="not enough values to unpack":
							coordinates, nucleotides, quals= [(), (), ()] 
						else:
							raise
				else:
					coordinates, nucleotides, quals = [(), (), ()] 
				if mate.mapping_quality >= mmq:
					try: 
                                                adj_seq, adj_qual =  adjust_seq(mate)
						coordinates_mate, nucleotides_mate, quals_mate = list(zip(*[(int(_x), _y, _z) for _x, _y, _z in zip(mate.positions, adj_seq.upper(), list(ord(_x)-qoffset for _x in adj_qual)) if _z>=mbq])) 
					except ValueError as e:
						if e.args[0][0:len("need more than 0 values to unpack")]=="need more than 0 values to unpack" or e.args[0][0:len("not enough values to unpack")]=="not enough values to unpack":
							coordinates_mate, nucleotides_mate, quals_mate = [(), (), ()] 
						else:
							raise
				else:
					coordinates_mate, nucleotides_mate, quals_mate = [(), (), ()] 
				if break_mate:
					pass
				else: # merge the sub-reads if the insert-size is less than maxIS
					try: 
						coordinates, nucleotides, quals = list(zip(*sorted(zip(coordinates+coordinates_mate, nucleotides + nucleotides_mate, quals+quals_mate), key = lambda x: x[0])))
					except ValueError as e:
						if e.args[0][0:len("need more than 0 values to unpack")]=="need more than 0 values to unpack" or e.args[0][0:len("not enough values to unpack")]=="not enough values to unpack":
							coordinates, nucleotides, quals = [(), (), ()] 
						else:
							raise
			else: 
				break_mate = True
				if read.mapping_quality >= mmq:
					try:
                                                adj_seq, adj_qual =  adjust_seq(read)
                                                coordinates, nucleotides, quals = list(zip(*[(int(_x), _y, _z) for _x, _y, _z in zip(read.positions, adj_seq.upper(), list(ord(_x)-qoffset for _x in adj_qual)) if _z>=mbq]))
					except ValueError as e:
						if e.args[0][0:len("need more than 0 values to unpack")]=="need more than 0 values to unpack" or e.args[0][0:len("not enough values to unpack")]=="not enough values to unpack":
							coordinates, nucleotides, quals = [(), (), ()] 
						else:
							raise
				else:
					coordinates, nucleotides, quals = [(), (), ()] 
				coordinates_mate, nucleotides_mate, quals_mate = [(), (), ()] 
			if break_mate:
				pass
			else:
				unique_q = []
				unique_c = []
				unique_n = []
				for _n, _c in enumerate(coordinates): # remove the duplicates from overlapping positions 
					try:
						if unique_c[-1]!=_c:
							#print("\nN:\n",nucleotides, "\nQ:\n", quals)
							unique_c.append(_c)
							unique_n.append(nucleotides[_n])
							unique_q.append(quals[_n])
						elif unique_n[-1]==nucleotides[_n]:
							unique_q[-1] = min(126-qoffset, unique_q[-1]+quals[_n])
						else: # if the called nucleotides differ at overlapping sites, use the one with the highest Phred score and adjust the Phred score.
							if quals[_n]>unique_q[-1]:
								_new_q_score = round(-10*log(1-10**(-unique_q[-1]/10)*(1-10**(-quals[_n]/10)), 10), 5) # Q=-10log(p,10)
								if _new_q_score >= mbq:
									unique_n[-1] = nucleotides[_n]
									unique_q[-1] = _new_q_score
								else:
									del(unique_c[-1], unique_n[-1], unique_q[-1])
							else:
								_new_q_score = round(-10*log(1-(1-10**(-unique_q[-1]/10))*10**(-quals[_n]/10), 10), 5)
								if _new_q_score >= mbq:
									unique_q[-1] = _new_q_score
								else:
									del(unique_c[-1], unique_n[-1], unique_q[-1])
					except IndexError:
						unique_c.append(_c)
						unique_n.append(nucleotides[_n])
						unique_q.append(quals[_n])
				coordinates, nucleotides, quals = [unique_c, unique_n, unique_q]
			coordinates = list(_x+1 for _x in coordinates) # Convert the zero-based BAM coordinates to 1-based, as the coordinates are 1-based in the VCF (like the SAM format).
			new_frag, new_q = frag_gen(varpos, allelelst, genolst[_group], coordinates, nucleotides, quals)
			frag_lst.append(new_frag)
			q_lst.append(new_q)
			if break_mate:
				coordinates_mate = list(_x+1 for _x in coordinates_mate)
				new_frag_mate, new_q_mate = frag_gen(varpos, allelelst, genolst[_group], coordinates_mate, nucleotides_mate, quals_mate)
				frag_lst.append(new_frag_mate)
				q_lst.append(new_q_mate)
		try:
			frag_lst, q_lst = [_lst for _lst in zip(*[(_x, _y) for _x, _y in zip(frag_lst, q_lst) if not _x.isNULL()])]
		except ValueError as e:
			if e.args[0][0:len("need more than 0 values to unpack")]=="need more than 0 values to unpack" or e.args[0][0:len("not enough values to unpack")]=="not enough values to unpack":
				frag_lst, q_lst = [], []
		Frag_lsts.append(frag_lst)
		Q_lsts.append(q_lst)
	return Frag_lsts, Q_lsts, GroupedSM

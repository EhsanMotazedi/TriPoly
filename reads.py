# Written by Ehsan Motazedi, Wageningen UR, 27-07-2016.
# Last Updated: 22-12-2017

import ast
import networkx as nx

def remove_last_SNP(r):
	""" Remove the last SNP from a read."""
	d = r.GetDict() #ast.literal_eval(repr(r))
	try:
		del d[max(d.keys())]
	except ValueError:
		return Read(d)
	return Read(d)

def Extract_Components(fragmentlist):
	""" Make the connected Read_Graphs from the list of Read objects, so that each connected Read_Graph is phased separately."""
	Read_Graph = nx.Graph()    # The SNP-fragment graph to be grown from the input fragment file
	_fragdict, _edges, _ls = (dict(), [], []) # Dictionary and lists to add each fragment to the SNP-fragment graph
	for _frag in fragmentlist: # Add the fragments to the graph nodes, and add an edge between every two variants that lie on the same fragment
		_fragdict = _frag.GetDict() #dict(ast.literal_eval(repr(_frag))) # Convert each fragment Read object to a conventional dictionary
		Read_Graph.add_nodes_from(_fragdict.keys()) # Adding the fragments as nodes to the graph
		_ls = list((x,y) for x in _fragdict.keys() for y in _fragdict.keys()) # Considering an edge between each two variants within the same fragment
		_edges = list(edge for edge in _ls if edge[0]!=edge[1]) # Adding the edges to the graph, avoiding self-loops
		Read_Graph.add_edges_from(_edges)
	_edges, _ls, _fragdict = (None, None, None)
	return sorted(nx.connected_components(Read_Graph), key=lambda x: min(x))

def getReads(fragfile, type=0):
	"""Extract the reads from a HapTree fragment file containing a dictionary for each (paired-end) fragment (type=0),\n 
as produced by 'fragmentpoly.py' or a HapCut (Bansal & Bafina 2008) format file containing a tab-delimited row for each\n
fragment (type=1). Return the reads as Read objects in a list."""
	Rout=[]
	with open(fragfile, 'rU') as ffile:
		if type==0:
			for _strdict in ffile:
				_dict = ast.literal_eval(_strdict)
				if len(_dict)>1: 
					Rout.append(Read(_dict)) # Add the read is it has more than two variants
		else:
			for _rec in ffile:
				_reclst = _rec.split()
				del _reclst[1], _reclst[-1] # get rid of the fragment name and quality scores
				if len(''.join(_reclst[_x] for _x in range(2, len(_reclst), 2)))>1:
					Rout.append(Read([int(_reclst[_x])-1 for _x in range(1, len(_reclst)-1, 2)],# convert 1, 2, ... numbering in HapCut to 0, 1, ... 
						[_reclst[_x] for _x in range(2, len(_reclst), 2)]))
	return Rout

class Read(object):
	def __init__(self, *args):
		""" Two ways to construct a read object:\n 1) Give a fragment dictionary as input.\n 2) Give a vector input:\n 
[start positions (single int or list/tuple: first part, second part(maybe None)),\n
alleles (single string or list/tuple: first part, second part (maybe None))].\nNo input corresponds to an empty dictionary.""" 
		if len(args)==0:
			args = [dict()]
		if len(args)==1:
			if not isinstance(args[0], dict):
				raise TypeError("The single input must be a dictionary!")
			rdict = args[0]
			self.pos = sorted((int(_x) for _x in rdict.keys()), key=lambda p: int(p))
			self.alleles = []	
		else:
			self.start, self.alleles = args
			if (not isinstance(self.start, int)) and (not all(isinstance(_start, int) for _start in self.start)):
				raise TypeError("The given start positions of the (sub)reads must be integers!")
			self.pos = []
			if isinstance(self.start, int):
				self.start = [self.start,]
			if isinstance(self.alleles, str):
				self.alleles = [self.alleles,]
			for _n, _start in enumerate(self.start):
				self.pos+=range(_start, _start+len(self.alleles[_n]))
			rdict = dict(zip(self.pos, list(''.join(self.alleles))))
			self.alleles=[]
		try:
			self.start = self.pos[0]
		except IndexError: # This exception occurs if the position list, self.pos, is empty
			self.start = 0
		try:
			self.end = self.pos[-1]
		except IndexError:
			self.end = -1 
		for _pos in range(self.start, self.end+1):
			if _pos==self.pos[0]:
				self.alleles.append(rdict[_pos])
				del(self.pos[0])
			else:
				self.alleles.append('-')
		self.pos = list(range(self.start, self.end+1))
		self.alleles = self.alleles
	def __eq__(self, other):
		return isinstance(other, Read) and self.alleles == other.allele and self.pos == other.pos
	def __repr__(self):
		rdict = dict()
		for _pos, _allele in zip(self.pos, self.alleles):
			if _allele!='-':
				rdict[_pos] = _allele
		return str(rdict)
	def __str__(self):
		if self.pos==[]:
			return 'None'
		return "begin={0}, end={1}, alleles={2}".format(self.start, self.end, ''.join(str(_x) for _x in self.alleles))
	def __len__(self):
		return len([_allele for _allele in self.alleles if _allele not in set(['.','-'])]) 
	def isNULL(self):
		return self.GetDict()=={} 
		#return ast.literal_eval(repr(self))=={}
	def GetalleleSet(self, start, stop):
		if self.pos==[]:
			return ''
		return ''.join(str(_x) for _x in self.alleles[min(self.end-self.start+1, max(0,start-self.start)):max(0, min(stop+1-self.start, self.end-self.start+1))])
	def GetBegin(self):
		if self.pos==[]:
			return 'None'
		return self.start
	def GetEnd(self):
		if self.pos==[]:
			return 'None'
		return self.end
	def GetDict(self):
		rdict = dict()
		for _pos, _allele in zip(self.pos, self.alleles):
			if _allele!='-':
				rdict[_pos] = _allele
		return rdict
	def GetPos(self):
		return sorted(_pos for _pos in (self.GetDict()).keys())
		#return sorted(_pos for _pos in dict(ast.literal_eval(repr(self))).keys())

def SemiRead(R, s):
	""" Make a semi-read relative to the SNP position 's' from a read 'R' (Berger et al. 2014 p. 4)."""
	if not isinstance(R, Read):
		raise TypeError("The argument to SemiRead must be a Read object!")
	d = R.GetDict()
	if s not in d.keys():
		return Read(dict())
	rdict = dict()
	for _s, _allele in d.iteritems():
		if _s<s+1:
			rdict[_s]=str(_allele)
	if len(rdict)<2:
		return Read(dict())
	else:
		return Read(rdict)

def SemiRead_old(R, s):
	""" Make a semi-read relative to the SNP position 's' from a read 'R' (Berger et al. 2014 p. 4)."""
	if not isinstance(R, Read):
		raise TypeError("The argument to SemiRead must be a Read object!")
	if s<R.GetBegin()+1 or R.GetPos()==[]:
		return Read(dict())
	rdict = dict()
	for _n, _allele in enumerate(R.GetalleleSet(R.start, s)):
		if _allele!='-':
			rdict[R.start+_n]=_allele
	if s not in rdict.keys() or len(rdict)<2:
		return Read(dict())
	else:
		return Read(rdict)

def SubReads(Rset, Sset, min_subread_length=2):
	"""Obtain the sub-reads of a read set, i.e. those reads used 
	for phasing of Sset and only those read positions relevant to phasing Sset
	(Berger et al. 2014 p. 4)."""
	Rout=[]
	if not all(isinstance(r, Read) for r in Rset):
		raise TypeError("The argument to SubReads must be a list of Read objects!")
	for r in Rset:
		_Sset=set(_s for _s in Sset)
		nrdict = dict()
		rdict = r.GetDict() #dict(ast.literal_eval(repr(r)))
		for s in rdict.keys():
			if s in _Sset:
				_Sset.remove(s)
				nrdict[s]=rdict[s]
		if len(nrdict) >= min_subread_length:
			Rout.append(Read(nrdict))
		else:
			Rout.append(Read(dict()))	
	return Rout

#diclst=[{1:1, 2:1, 3:1, 4:1}, {3:1, 4:1, 5:0, 6:0}, {4:0, 5:1, 6:1}, 
#{4:0, 5:1,6:1, 7:0}, {5:0, 6:0, 7:1}, {5:1, 6:1, 7:0}]
#a=Read({14:1, 15:1, 16:0, 21:0, 22:1})
#b=Read({14:1, 15:1, 16:0})
#Read(ast.literal_eval(repr(a)))
#str(a)
#str(b)
#a.GetalleleSet(16,22)
#c=Read(dict())
#nlst=map(lambda rdict: Read(rdict), diclst)
#for n in range(1,8):
#	for r in nlst:
#		print("{0}->{1}".format(n, repr(SemiRead(r, n))))

#SubReads(nlst, {1,2,3,4,5})
#a='1 HWI-ST163_0386:6:27:3566:6857#P7PEM12p 14 01010 ;:;78'.split()
#b='2 HWI-ST163_0386:6:6:10650:188444#P7PEM12p_MP 5 0100 17 000 ;;7<779'.split()
#del a[1], b[1], a[-1], b[-1]
#Read(a[0], int(a[1]), a[2])
#Read(b[0], (int(b[1]), int(b[3])), (b[2], b[4]))
#Read(b[0], (b[1], b[3]), (b[2], b[4]))

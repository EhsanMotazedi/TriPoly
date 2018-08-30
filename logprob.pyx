# Written by Ehsan Motazedi, Wageningen UR, 05-02-2017.
# Last Updated: 01-09-2017

import math
import sys

cpdef int _xor(x, y):
	return 0 if str(x)==str(y) else 1

cpdef Hamming(a,  b, null_lst = [".", "-"]):
	""" Hamming distance and similarity between two strings or lists or tuples. 
	(Hamming, Richard W., Error detecting and error correcting codes, 
	Bell System technical journal 29.2 (1950), pp. 147-160)."""
	if isinstance(a, int):
		a = [a]
	if isinstance(b, int):
		b = [b]
	if len(a) != len(b):
		raise ValueError('The length of the two inputs must be equal for comparison!')
	if isinstance(a, str):
		a = list(a)
	if isinstance(b, str):
		b = list(b)
	hamvec = [_xor(x, y) for x, y in zip(a,b) if (str(x) not in null_lst) and (str(y) not in null_lst)]
	hamgelijk = hamvec.count(0)
	hamafstand = hamvec.count(1)
	return hamafstand, hamgelijk

cpdef diffVec(a, b, null_lst = [".","-"]):
	""" Obtain differences between two lists. Return two vectors: one specifying positions of similarity and one positions of difference."""
	if isinstance(a, int):
		a = [a]
	if isinstance(b, int):
		b = [b]
	if len(a) != len(b):
		raise ValueError('The length of the two inputs must be equal for comparison!')
	if isinstance(a, str):
		a = list(a)
	if isinstance(b, str):
		b = list(b)
	ongelijk, gelijk = [], []
	for n, xy in enumerate(zip(a,b)):
		if (str(xy[0]) in null_lst) or (str(xy[1]) in null_lst):
			continue
		if _xor(*xy)==1:
			ongelijk.append(n)
		else:
			gelijk.append(n)
	return ongelijk, gelijk

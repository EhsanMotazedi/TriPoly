# Break complex variants and throw out indels from the input vcf file. i
# EX: ./break_vcf.sh my.vcf out.vcf
# Written by Ehsan Motazedi, Wageningen UR, 17-10-2017
# Last Modified: 27-07-2018
file=$1
outfile=$2
tmpfile=$(mktemp /tmp/dummy.XXXXXX)
exec 3> "$tmpfile" #open tmpfile for writing on file descpritor 3
python << EOS
import re
with open("$file", 'rU') as vcf, open("$tmpfile", "w") as recoded:
        detected_SNPs = dict()
        for line in vcf:
                if line.lstrip()[0]=='#':
                        recoded.write(line)
                else:
                        fields = line.split()
			contig = fields[0]
			if contig not in detected_SNPs.keys():
				detected_SNPs[contig] = []
                        if '.' in fields[3] or '.' in fields[4]: #indels
                                pass
                        elif ',' in fields[4]: # multi-allelic variants
                                _temp_alts = fields[4].split(',')
                                if not all(len(_alt)==len(fields[3]) for _alt in _temp_alts): # do not parse complex variants that contain indels
                                        pass
                                else: # break complex variants into bi-allelic SNPs if possible
                                        difflocation = []
                                        diff_ref_alts = []
                                        diff_ref_alts_dic = []
                                        for _base in range(0, len(fields[3])):
                                                ref_alt = [fields[3][_base]] + [_x for _x in set(_alt[_base] for _alt in _temp_alts)-set([fields[3][_base]])] # check at each position within the complex variant if ref, alt corresponds to a SNP
						if len(ref_alt)>1:
							ref_alt_dic = {ref_alt[_n]:str(_n) for _n in range(0, len(ref_alt))}
							ref_alt_dic['.']='.' # missing genotype added
							difflocation.append(_base)
							diff_ref_alts.append(ref_alt)
							diff_ref_alts_dic.append(ref_alt_dic)
                                        for _n in range(0, len(difflocation)):
                                                _pos = str(int(fields[1])+difflocation[_n])
                                                if (_pos not in detected_SNPs[contig]):
                                                        recoded.write(fields[0]+'\t'+_pos+'\t'+fields[2]+'\t'+diff_ref_alts[_n][0]+'\t'+','.join(diff_ref_alts[_n][1:])+'\t'+'\t'.join(fields[5:9]))
                                                        for _genotype in fields[9:]:
                                                                _dosage = re.split('/|\|', _genotype.split(':')[0])
                                                                _old_dic = {'0':diff_ref_alts[_n][0], '.':'.'} #'.' stands for missing genotypes
                                                                _old_dic.update({str(_m+1):_temp_alts[_m][difflocation[_n]] for _m in range(0, len(_temp_alts))})
                                                                _new_dosage = '/'.join(diff_ref_alts_dic[_n][_old_dic[_allele]] for _allele in _dosage) 
                                                                recoded.write('\t'+_new_dosage+':'+':'.join(_genotype.split(':')[1:]))
                                                        recoded.write('\n')
                                                        detected_SNPs[contig].append(_pos)
                                                else:
                                                        pass
			elif len(fields[3])!=len(fields[4]): #indels
				pass
                        elif len(fields[3])==1: #biallelic SNPs
                                if (fields[1] not in detected_SNPs[contig]):
                                        recoded.write(line)
                                        detected_SNPs[contig].append(fields[1])
                                else:
                                        pass
                        else: # MNPs
                                difflocation = []
                                for _n in range(0, len(fields[3])):
                                        if fields[3][_n]!=fields[4][_n]:
                                                difflocation.append(_n)
                                for _n in range(0, len(difflocation)):
                                        _pos = str(int(fields[1])+difflocation[_n])
                                        if (_pos not in detected_SNPs[contig]):
                                                recoded.write(fields[0]+'\t'+_pos+'\t'+fields[2]+'\t'+fields[3][difflocation[_n]]+'\t'+fields[4][difflocation[_n]]+'\t'+'\t'.join(fields[5:])+'\n')
                                                detected_SNPs[contig].append(_pos)
                                        else:
                                                pass
EOS
cp $tmpfile $outfile
status=$?
if [ $status -eq 0 ]; then
	rm $tmpfile
fi

import gzip
import argparse
from textwrap import dedent
from subprocess import PIPE, Popen

class Faidx_processor:
    # samtools faidx /root/ceph/staging/volume/pub/hg19.resource/1/Annotator/referenceGenome/ucsc.hg19.fasta
    def __init__(self, path_samtools: str, path_ref: str):
        self.path_samtools = path_samtools
        self.path_ref = path_ref
        self.base_cmd = [self.path_samtools, 'faidx', self.path_ref]
    
    def getNucleotide(self, chromosome, start, end) -> str:
        if chromosome.startswith('chr') is False:
            chromosome = 'chr' + chromosome
        if chromosome == "chrMT":
            chromosome = "chrM"
        cmd = self.base_cmd + ["%s:%s-%s" % (chromosome, start, end)]
        proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
        try:
            outs, errs = proc.communicate()
            res = outs.decode('utf-8').split()
        except:   # no result from referenceGenome
            return None
        if len(res) == 2:
            return res[1].upper()
        else:
            print("getNucleotide error(%s) with cmd: %s\n" % (errs, ' '.join(cmd)))
        return None

# process = Faidx_processor('samtools', '/root/ceph/staging/volume/pub/hg19.resource/1/Annotator/referenceGenome/ucsc.hg19.fasta')
# process.getNucleotide('chrX', '129299784', '129299784')
def variantIndelConverter(proc_faidx, chromosome, start, end, ref, alt, to_dash):
    if alt == '*':
        alt = '-'
    if ref == '*':
        ref = '-'
    if to_dash is False:   # input with '-'
        if len(ref) != len(alt) and ref != '-' and alt != '-':
            return chromosome, start, end, ref, alt
        if ref != '-' and alt != '-' and len(ref) == 1 and len(alt) == 1:   # SNV
            return chromosome, start, end, ref, alt
        elif ref == '-':   # INS
            return dash_to_noDash_INS(proc_faidx, chromosome, start, end, ref, alt)
        elif alt == '-':   # DEL
            return dash_to_noDash_DEL(proc_faidx, chromosome, start, end, ref, alt)
    elif to_dash is True:   # input without '-'
        if len(ref) == 1 and len(alt) == 1:   # SNV
            return chromosome, start, end, ref, alt
        elif len(ref) < len(alt):   # INS
            return noDash_to_dash_INS(chromosome, start, end, ref, alt)
        elif len(ref) > len(alt):   # DEL
            return noDash_to_dash_DEL(chromosome, start, end, ref, alt)
    return None

def dash_to_noDash_INS(proc_faidx, chromosome, start, end, ref, alt):
    """ chr1 13417 13417 - GAGA   to   chr1 13417 13417 C CGAGA (require fetch from reference)"""
    nt = proc_faidx.getNucleotide(chromosome, start, start)
    if nt is not None:
        return chromosome, start, end, nt, nt + alt
    return chromosome, start, end, ref, alt

def dash_to_noDash_DEL(proc_faidx, chromosome, start, end, ref, alt):
    """ chr1 10145 10145 A -   to   chr1 10144 10145 TA T (require fetch from reference)"""
    nt = proc_faidx.getNucleotide(chromosome, start - 1, start - 1)
    if nt is not None:
        return chromosome, start - 1, end, nt + ref, nt
    return chromosome, start, end, ref, alt

def noDash_to_dash_INS(chromosome, start, end, ref, alt):
    """ chr1 13417 13417 C CGAGA   to   chr1 13417 13417 - GAGA """
    return chromosome, start, end, '-', alt[1:]

def noDash_to_dash_DEL(chromosome, start, end, ref, alt):
    """ chr1 10144 10145 TA T   to   chr1 10145 10145 A - """
    return chromosome, start + 1, end, ref[1:], '-'

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
Testing environment: Python 3
Require inputs:
1. in_file
2. out_file
3. in_reference
4. type (txt, vcf)
5. to_dash (ex. chr1 99 99 T TA  to  chr1 100 100 - A)

Usage:
- python3 indelConverter.py --in_file data/dash.txt --in_reference ucsc.hg19.fasta --out_file data/out_dash.txt 
"""))
    required = parser.add_argument_group('required arguments')
    required.add_argument('--in_file', type=str, help='input file', required=True)
    required.add_argument('--out_file', type=str, help='output file', required=True)
    required.add_argument('--in_reference', type=str, help='input corresponding reference fasta', required=True)
    required.add_argument('--type', type=str, default='txt', choices=['txt', 'vcf'], help='input file type')
    required.add_argument('--to_dash', action='store_true', help="if true, convert indel in input file to format with '-', else on the contrary")
    required.add_argument('--cmd_samtools', type=str, help="samtools execute path", default='samtools')
    args = parser.parse_args()

    proc_faidx = Faidx_processor(args.cmd_samtools, args.in_reference)

    out = open(args.out_file, 'w')
    f = gzip.open(args.in_file, 'rb') if args.in_file.endswith('.gz') else open(args.in_file)
    n = 0    
    for line in f:
        if args.in_file.endswith('.gz'):
            line = line.decode('utf-8')
        if line.startswith('#'):
            out.write(line)
            continue
        sep = line.strip('\n').split('\t')
        try:
            chromosome = sep[0]
            ref = sep[3]
            alt = sep[4]
            others = sep[5:]
            if args.type == 'txt':
                start = int(sep[1])
                end = int(sep[2])
            elif args.type == 'vcf':
                start = int(sep[1])
                end = start + len(ref) - 1
            # print('+', '\t'.join([chromosome, str(start), str(end), ref, alt]))
            chromosome, start, end, ref, alt = variantIndelConverter(proc_faidx, chromosome, start, end, ref, alt, args.to_dash)
            # print('-', '\t'.join([chromosome, str(start), str(end), ref, alt]), '\n')
            out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(chromosome, str(start), str(end), ref, alt, '\t'.join(others)))
        except Exception as e:
            print(e)
            out.write(line)

        n += 1
    print('processed %s variants' % n)

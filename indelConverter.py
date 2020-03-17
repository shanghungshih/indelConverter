import sys
import docker
import argparse
from textwrap import dedent

def variantIndelConverter(client, volumes, container, reference, chromosome, start, end, ref, alt, to_dash):
    if to_dash is False:   # input with '-'
        if ref != '-' and alt != '-' and len(ref) == 1 and len(alt) == 1:   # SNV
            return chromosome, start, end, ref, alt
        elif ref == '-':   # INS
            return dash_to_noDash_INS(client, volumes, container, reference, chromosome, start, end, ref, alt)
        elif alt == '-':   # DEL
            return dash_to_noDash_DEL(client, volumes, container, reference, chromosome, start, end, ref, alt)
    elif to_dash is True:   # input without '-'
        if len(ref) == 1 and len(alt) == 1:   # SNV
            return chromosome, start, end, ref, alt
        elif len(ref) < len(alt):   # INS
            return noDash_to_dash_INS(chromosome, start, end, ref, alt)
        elif len(ref) > len(alt):   # DEL
            return noDash_to_dash_DEL(chromosome, start, end, ref, alt)
    return None

def getNucleotide(client, volumes, container, reference, chromosome, start, end):
    if chromosome.startswith('chr') is False:
        chromosome = 'chr' + chromosome
    logs = container.exec_run("samtools faidx /ref/%s %s:%s-%s" %(reference.split('/')[-1], chromosome, start, end), stdout=True, stderr=True, stream=True)
    for i in logs[1]:
        res = i.decode('utf-8').split()
    if len(res) == 2:
        return res[1].upper()
    return None

def dash_to_noDash_INS(client, volumes, container, reference, chromosome, start, end, ref, alt):
    """ chr1 13417 13417 - GAGA   to   chr1 13417 13417 C CGAGA (require fetch from reference)"""
    nt = getNucleotide(client, volumes, container, reference, chromosome, start, start)
    if nt is not None:
        return chromosome, start, end, nt, nt + alt
    return chromosome, start, end, ref, alt

def dash_to_noDash_DEL(client, volumes, container, reference, chromosome, start, end, ref, alt):
    """ chr1 10145 10145 A -   to   chr1 10144 10145 TA T (require fetch from reference)"""
    nt = getNucleotide(client, volumes, container, reference, chromosome, start - 1, start - 1)
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
- python3 indelConverter.py --in_file data/dash.txt --in_reference ucsc.hg19.fasta --out_file data/out_dash.txt --type txt --to_dash false
"""))
    required = parser.add_argument_group('required arguments')
    required.add_argument('--in_file', type=str, help='input file')
    required.add_argument('--out_file', type=str, help='output file')
    required.add_argument('--in_reference', type=str, help='input corresponding reference fasta')
    required.add_argument('--type', type=str, help='input file type')
    required.add_argument('--to_dash', type=str2bool, help="if true, convert indel in input file to format with '-', else on the contrary")
    args = parser.parse_args()

    in_reference_path = '/'.join(args.in_reference.split('/')[:-1]) if '/' in args.in_reference else '$PWD'
    client = docker.from_env()
    volumes = {
        in_reference_path: {'bind': '/ref', 'mode': 'rw'}
    }
    container = client.containers.run(image='biocontainers/samtools:1.3.1',
                                    command="/bin/bash",
                                    volumes=volumes,
                                    tty=True,
                                    stdin_open=True,
                                    detach=True)
    container.start()

    out = open(args.out_file, 'w')
    with open(args.in_file) as f:
        for line in f:
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
                chromosome, start, end, ref, alt = variantIndelConverter(client, volumes, container, args.in_reference, chromosome, start, end, ref, alt, args.to_dash)
                # print('-', '\t'.join([chromosome, str(start), str(end), ref, alt]), '\n')
                out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(chromosome, str(start), str(end), ref, alt, '\t'.join(others)))
            except Exception as e:
                print(e)
                out.write(line)

    container.stop()
    container.remove()

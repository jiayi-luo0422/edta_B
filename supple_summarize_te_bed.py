#!/usr/bin/env python3
# 用法: python supple_summarize_te_bed.py SCZ_TE_all.bed SCZ.chrom.fasta.fai
import sys
from collections import defaultdict

bed, fai = sys.argv[1], sys.argv[2]
G = sum(int(x.split('\t')[1]) for x in open(fai))   # 基因组大小

def sf(cls):
    c = cls.lower()
    if 'gypsy' in c:        return 'LTR/Gypsy'
    if 'copia' in c:        return 'LTR/Copia'
    if c.startswith('ltr'): return 'LTR/unknown'
    if 'line' in c:         return 'LINE'
    if 'sine' in c:         return 'SINE'
    if 'helitron' in c or c.startswith('rc'): return 'DNA/Helitron'
    if 'mite' in c:         return 'MITE'
    if c.startswith('dna') or c.startswith('tir'): return 'DNA/TIR'
    return 'Unclassified'

bp, n = defaultdict(int), defaultdict(int)
for L in open(bed):
    f = L.rstrip('\n').split('\t')
    name = sf(f[3].split('|')[0])
    bp[name] += int(f[2]) - int(f[1])
    n[name]  += 1

order   = ['LTR/Copia','LTR/Gypsy','LTR/unknown','LINE','SINE',
           'DNA/TIR','DNA/Helitron','MITE','Unclassified']
classI  = {'LTR/Copia','LTR/Gypsy','LTR/unknown','LINE','SINE'}
classII = {'DNA/TIR','DNA/Helitron','MITE'}

def row(lab, c, l):
    print(f"{lab:<16}{c:>12,}{l:>16,}{l/1e6:>12.2f}{100*l/G:>10.2f}")

print(f"# Genome size: {G:,} bp ({G/1e6:.1f} Mb)\n")
print(f"{'Class':<16}{'Count':>12}{'Length(bp)':>16}{'Length(Mb)':>12}{'%Genome':>10}")
print('-'*66)
tb=tn=I=II=0
for k in order:
    if not bp[k]: continue
    row(k, n[k], bp[k]); tb+=bp[k]; tn+=n[k]
    if k in classI:  I  += bp[k]
    if k in classII: II += bp[k]
print('-'*66)
row('Class I (Retro)',  sum(n[k] for k in classI),  I)
row('Class II (DNA)',   sum(n[k] for k in classII), II)
row('Total TE', tn, tb)
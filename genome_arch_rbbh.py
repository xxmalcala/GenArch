#!/usr/bin/env python3

"""
Categorizes the outputs from a self BLAST ALL-vs-ALL into types of duplications.
Note that this assumes a particular protein-naming convention:
    taxon_name_XX_protein_name

Additionally, protein_name is expected to include the scaffold name as well (i.e.
the output from AUGUSTUS and BRAKER1/2's gene prediction pipelines).

Dependencies: DIAMOND, BioPython
"""

import os, sys
from collections import defaultdict

from Bio import SeqIO

def prep_fasta(prot_fasta, taxon_code):
    seq_recs = []
    for i in SeqIO.parse(prot_fasta,'fasta'):
        i.id = f'{taxon_code}_XX_{i.id}'
        i.description = ''
        seq_recs.append(i)

    SeqIO.write(seq_recs, f'{taxon_code}.CDS.AA.fas','fasta')

    return f'{taxon_code}.CDS.AA.fas', seq_recs


def self_v_self(prot_fasta, stringent = False):
    dmnd_cmd = 'diamond blastp ' \
        f'-q {prot_fasta} ' \
        f'-d {prot_fasta} ' \
        f'-o {prot_fasta.replace(".fas",".SelfHits.tsv")} ' \
        '-f 6 --ultra-sensitive'

    if stringent:
        dmnd_cmd += ' --id 30 --query-cover 40 --subject-cover 40'

    os.system(dmnd_cmd)

    return f'{prot_fasta.replace(".fas",".SelfHits.tsv")}'


def get_rbbh(tsv, ncbi = True):
    hits = defaultdict(list)
    rbbh = []
    seen = []
    data = []
    for i in open(tsv).readlines():
        if not ncbi:
            if list(set([j.rsplit(".",1)[-1] for j in i.split('\t')[:2]])) != ['t1']:
                continue
        data.append(i)

    data.sort(key = lambda x: -float(x.split('\t')[-1]))

    for i in data:
        h = i.split('\t')[:2]
        if len(set(h)) == 2:
            hits[h[0]].append((h[1], float(i.split('\t')[-1])))

    for k, v in hits.items():
        if k not in seen:
            seen.append(k)
            q = v[0][0]
            if q in hits.keys():
                s = hits[q][0][0]
                if s == k:
                    if q not in seen:
                        seen.append(q)
                        rbbh.append((k, q))
    return rbbh


def parse_gff3(gff3_file, seqs_to_eval, ncbi = True):
    chrom_pos = defaultdict(list)
    gene_scf = {}
    gene_prot_id = {}

    for line in open(gff3_file).readlines():
        if '\tmrna\t' in line.lower():

            mrna_id = line.split(';')[0].split('=')[-1].strip('rna-')
            gene_prot_id[mrna_id] = [line.split('\t')[0], int(line.split('\t')[3]), int(line.split('\t')[4])]

            if not ncbi and mrna_id in seqs_to_eval:
                gene_prot_id[mrna_id].insert(1, mrna_id)

        if ncbi and '\tcds\t' in line.lower():

                prot_id = line.split(';')[0].split('cds-')[1]

                if prot_id in seqs_to_eval:
                    mrna_id = line.split(';')[1].split('rna-')[1]

                    if prot_id not in gene_prot_id[mrna_id]:
                        gene_prot_id[mrna_id].insert(1, prot_id)


    for k, v in gene_prot_id.items():
        if len(v) == 4:
            chrom_pos[v[0]].append(tuple(v[1:]))
            gene_scf[v[1]] = [v[0], tuple(v[1:])]

    for k, v in chrom_pos.items():
        v.sort(key = lambda x: min(x[-2:]))
        chrom_pos[k] = v

    return chrom_pos, gene_scf


def eval_rbbh(rbbh, chrom_coords, gene_scf):
    smry = []
    for i in rbbh:
        dup_type = 'Dispersed'
        scf_type = 'Different'
        p1, p2 = i[0], i[1]
        g1, g2 = p1.split('_XX_')[1], p2.split('_XX_')[-1]
        scf1, scf2 = gene_scf[g1], gene_scf[g2]
        g1_pos = chrom_coords[scf1[0]].index(scf1[1])
        g2_pos = chrom_coords[scf2[0]].index(scf2[1])

        if scf1[0] == scf2[0]:
            scf_type = 'Same'
            if abs(g1_pos-g2_pos) <= 5:
                 dup_type = 'Tandem'
            elif 5 < abs(g1_pos-g2_pos) <= 10:
                 dup_type = 'Proximal'
            else:
                 dup_type = 'Dispersed'
        smry.append(f'{p1}\t{p2}\t{scf1[0]}\t{scf2[0]}\t{scf_type}\t{g1_pos+1}\t{g2_pos+1}\t{dup_type}\n')
        smry.append(f'{p2}\t{p1}\t{scf2[0]}\t{scf1[0]}\t{scf_type}\t{g2_pos+1}\t{g1_pos+1}\t{dup_type}\n')
    return smry

def save_summary(output_name, dup_summary, stringent):
    same_prox = 0
    same_disp = 0
    same_tandem = 0
    diff_disp = 0

    fout_verb = f'{output_name}.RBBH_Duplications.tsv'
    fout_smry = f'{output_name}.RBBH_Duplication.Summary.tsv'

    if stringent:
        fout_verb = fout_verb.replace("tsv","Stringent.tsv")
        fout_smry = fout_smry.replace("tsv","Stringent.tsv")

    with open(fout_verb,'w+') as w:
        w.write(f'Protein-1\tProtein-2\tScaffold-1\tScaffold-2\tScaffold-Type\tProtein-1-Position\tProtein-2-Position\tDuplication-Type\n')
        for i in dup_summary:
            w.write(i)
            if 'Same' in i:
                if 'Proximal' in i:
                    same_prox += 0.5
                elif 'Dispersed' in i:
                    same_disp += 0.5
                else:
                    same_tandem += 0.5
            else:
                diff_disp += 0.5

    with open(fout_smry,'w+') as w:
        w.write('Duplication-Type\tNumber-Pairs\n')
        w.write(f'Same-Scaffold-Tandem\t{int(same_tandem)}\n')
        w.write(f'Same-Scaffold-Proximal\t{int(same_prox)}\n')
        w.write(f'Same-Scaffold-Dispersed\t{int(same_disp)}\n')
        w.write(f'Different-Scaffold-Dispersed\t{int(diff_disp)}\n')


if __name__ == '__main__':

    if 3 <= len(sys.argv[1:]) <= 4:
        fasta_file = sys.argv[1]
        gff3_file = sys.argv[2]
        output_name = sys.argv[3]
        try:
            if sys.argv[4]: stringent = True
        except:
            stringent = False
    else:

        print('\nUsage:\n\n    python3 genome_arch_rbbh.tsv [PROTEOME] [GFF3] [Taxon-Name]\n')
        sys.exit(1)

    ncbi = ' NCBI ' in open(gff3_file).read()


    prot_fas, good_seqs = prep_fasta(fasta_file, output_name)
    seqs_to_eval = [i.name for i in good_seqs]

    tsv = self_v_self(prot_fas, stringent)
    rbbh_info = get_rbbh(tsv, ncbi)

    chrom_coords, gene_scf = parse_gff3(gff3_file, seqs_to_eval, ncbi)

    dup_summary = eval_rbbh(rbbh_info, chrom_coords, gene_scf)
    save_summary(output_name, dup_summary, stringent)

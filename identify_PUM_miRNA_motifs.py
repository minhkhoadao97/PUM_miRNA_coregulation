import pandas as pd
import os, sys
import requests
from tssv import align
import pyBigWig
from Bio.Seq import reverse_complement

##################################################################################################################################################################

# split the attribute column of annotation df to seperate columns
def splitbracket(txt, pos=2):
    spl = txt.split(' ')
    out = spl[pos].replace('"', '')
    return out

def split_info(info):
    spl = info.split(';')
    
    geneid = ''
    transcriptid = ''
    genename = ''
    transcriptname = ''
    biotype = ''
    
    for i in spl:
        if 'gene_id' in i:
            geneid = splitbracket(i, 1)
        elif 'transcript_id' in i:
            transcriptid = splitbracket(i)
        elif 'gene_name' in i:
            genename = splitbracket(i)
        elif 'gene_biotype' in i:
            biotype = splitbracket(i)
        elif 'transcript_name' in i:
            transcriptname = splitbracket(i)
        else:
            continue
            
    return geneid, transcriptid, genename, transcriptname, biotype

def apply_split_info(x): return split_info(x['info'])

# find the idx in the annotation df for the transcript to be used
# incase there are multiple transcript, use the longest one

def find_longest_isoform(annotation):
    temp = annotation.reset_index()
    temp['len'] = temp['stop'] - temp['start']
    max_len = temp['len'].max()
    anno_idx = temp.loc[temp['len']==max_len, 'index']
    return anno_idx.values.astype(int)[0]
    

def filter_by_transcript(transcriptid):
 
    annotation = anno[anno['transcriptid']==transcriptid]
    if len(annotation) == 0:
        return 'no_transcript_found'
    elif len(annotation) > 1:
        return find_longest_isoform(annotation)
    else:
        return annotation.index.values.astype(int)[0]
    
def apply_filter_by_transcript(x): return filter_by_transcript(x['Transcript'])

# retrieve sequence and genomic location of miRNA sites
def retrieve_mir_sequence(transcript_start, transcript_stop, anno_idx, extend):
    
    chrom = anno.loc[anno_idx, 'chr']
    strand = anno.loc[anno_idx, 'strand']
    
    if strand == '+':
        strandno = 1
        start_3utr = anno.loc[anno_idx, 'start']
        genomic_start = start_3utr + transcript_start 
        extended_genomic_start = genomic_start - extend
        genomic_stop = start_3utr + transcript_stop -1 
        extended_genomic_stop = genomic_stop + extend

    elif strand == '-':
        strandno = -1
        start_3utr = anno.loc[anno_idx, 'stop']
        genomic_start = start_3utr - transcript_stop + 1 
        extended_genomic_start = genomic_start - extend
        genomic_stop = start_3utr - transcript_start
        extended_genomic_stop = genomic_stop + extend
    
    return chrom, strand, genomic_start, genomic_stop, extended_genomic_start, extended_genomic_stop

def apply_retrieve_mir_sequence(x): return retrieve_mir_sequence(x['mir_start'], x['mir_stop'], x['anno_idx'], 300)

# function to retrieve sequence using ENSEMBL REST API given a set of genomic coordinate
def retrieve_sequence(chrom, start, stop, strand):
    server = "http://rest.ensembl.org"
    ext = f"/sequence/region/mouse/{chrom}:{start}..{stop}:{strand}?"

    try:
        r = requests.get(server+ext, headers={"Content-Type" : "text/plain"})
    except:
        return 'cannot_find_sequence'
        
    if not r.ok:
        return 'cannot_find_sequence'
    else:
        return r.text

# function to find PUM motif in a sequence
strand_var = {'+': ['A', 'C', 'T'], '-': ['A', 'T', 'G']}

def find_PUM(chrom, eg_start, eg_stop, strand):
    
        if strand == '+':
            PUM='TGTAXATA'
            strandno = 1
        elif strand == '-':
            PUM = 'TATXTACA'
            strandno = -1
            
        seq = retrieve_sequence(chrom, eg_start, eg_stop, strandno)

        if seq == 'cannot_find_sequence':
            return seq

        PUMpos = []
        extra = 0
        start = 0
        stop = 0
        while len(seq) >= len(PUM):
                alignment = align(seq, PUM, 1000)
                if alignment['distance'] == 1:
                        endpos = alignment['position']
                        site = seq[endpos-8:endpos] # need to check for actual position
                        if strand == '+':
                                variable_pos = site[4]
                        elif strand == '-':
                                variable_pos = site[3]

                        start = endpos - 8 + extra
                        stop = endpos - 1 + extra

                        if variable_pos in strand_var[strand]:
                                PUMpos.append(f'{start}-{stop}')

                        seq = seq[endpos:]
                        extra = stop + 1
                else:
                        break
        return ' '.join(PUMpos)
    
def apply_find_PUM(x): return find_PUM(x['chr'], x['eg_start'], x['eg_stop'], x['strand'])

# function to retrieve genomic location of found PUM motifs
def convert_PUM(PUMpos, strand, extended_genomic_start, extended_genomic_stop):
    if len(PUMpos) == 0:
        return 'no_PUM', 0
    elif PUMpos == 'cannot_find_sequence':
        return 'error', 0
    else:
        spl = PUMpos.split(' ')
        out = []
        
        for coors in spl:
            coor = coors.split('-')
            if strand == '+':
                PUM_genomic_start = extended_genomic_start + int(coor[0])
                PUM_genomic_stop = extended_genomic_start + int(coor[1])
            elif strand == '-':
                PUM_genomic_start = extended_genomic_stop - int(coor[1])
                PUM_genomic_stop = extended_genomic_stop - int(coor[0])
                       
            out.append(f'{PUM_genomic_start}-{PUM_genomic_stop}')
                
        return ' '.join(out), len(spl)

def apply_convert_PUM(x): return convert_PUM(x['PUM_transcript_coor'], x['strand'],
                                            x['eg_start'], x['eg_stop'])

# if multiple PUM motifs is found, split each motif into a separate row
def split_rows(df):
    b = pd.DataFrame(df['PUM_genomic_coor'].str.split(' ').tolist(), index=df['index'])
    b = b.reset_index()
    b = pd.melt(b, id_vars = 'index', value_vars= b.columns.to_list()[1:], var_name = 'to_drop', 
                value_name = 'PUM_genomic_coor')
    b.drop(columns=['to_drop'], inplace=True)
    b.dropna(inplace=True)
    
    c = pd.DataFrame(df['PUM_transcript_coor'].str.split(' ').tolist(), index=df['index'])
    c = c.reset_index()
    c = pd.melt(c, id_vars = 'index', value_vars= c.columns.to_list()[1:], var_name = 'to_drop', 
                value_name = 'PUM_transcript_coor', ignore_index=True)
    c.drop(columns=['to_drop'], inplace=True)
    c.dropna(inplace=True)
    
    b['PUM_transcript_coor'] = c['PUM_transcript_coor']
    
    df.drop(columns=['PUM_transcript_coor', 'PUM_genomic_coor'], inplace=True)
    out = df.merge(b, on='index', how='left')
    
    return out

# retrieve phyloP score for the PUM and miRNA site
def retrieve_phyloP(chrom, mir_start, mir_stop, PUM_genomic_coor):
    
    if PUM_genomic_coor == 'no_PUM':
        return 'no_PUM', 'no_PUM', 'no_PUM', 'no_PUM', 'no_PUM', 'no_PUM', 'no_PUM'
    elif PUM_genomic_coor == 'error':
        return 'error', 'error', 'error', 'error', 'error', 'error', 'error'

    spl = PUM_genomic_coor.split('-')

    genomic_loc = [mir_start, mir_stop, int(spl[0]), int(spl[1])]
    start = min(genomic_loc) - 5
    stop = max(genomic_loc) + 5
    chrom = 'chr' + chrom
    
    avg = bw.stats(chrom, start, stop)
    max_val = bw.stats(chrom, start, stop, type='max')
    min_val = bw.stats(chrom, start, stop, type='min')
    coverage = bw.stats(chrom, start, stop, type='coverage')
    std = bw.stats(chrom, start, stop, type='std')
    return avg[0], max_val[0], min_val[0], coverage[0], std[0], start, stop

def apply_retrieve_phyloP(x): return retrieve_phyloP(x['chr'], x['g_start'], x['g_stop'], x['PUM_genomic_coor'])

############################################################################################################################################################
def pipeline(df):
    
    df['anno_idx'] = df.apply(apply_filter_by_transcript, axis=1)
    df = df[df['anno_idx']!='no_transcript_found'].reset_index()
    
    df['chr'],\
    df['strand'],\
    df['g_start'],\
    df['g_stop'],\
    df['eg_start'],\
    df['eg_stop'] = zip(*df.apply(apply_retrieve_mir_sequence, axis = 1))
    

    df['PUM_transcript_coor'] = df.apply(apply_find_PUM, axis = 1)
    df['PUM_genomic_coor'],df['PUM_count'] = zip(*df.apply(apply_convert_PUM, axis = 1))
    df = split_rows(df)
    
    df['phyloP_avg'],\
    df['phyloP_max'],\
    df['phyloP_min'],\
    df['phyloP_covrg'],\
    df['phyloP_std'],\
    df['phyloP_start'],\
    df['phyloP_stop'] = zip(*df.apply(apply_retrieve_phyloP, axis=1))
    
    return df

###########################################################################################################################################################
# import bigwig file for phyloP score
bw = pyBigWig.open('ENSEMBL/mm39.phyloP35way.bw')

# import ENSEMBL annotation file for 3'UTR
anno = pd.read_csv('ENSEMBL/mouse_3UTR_annotation.gtf.gz',\
                    usecols = [0, 3, 4, 6, 8],\
                    sep = '\t',\
                    compression = 'gzip',\
                    names = ['chr', 'start', 'stop', 'strand', 'info'])

# split attribute column 
anno['geneid'], anno['transcriptid'], anno['genename'], anno['transcript_name'], anno['gene_biotype'] = zip(*anno.apply(apply_split_info, axis=1))

# import mirgate prediction results and run pipeline

#for fn in os.listdir('mir_prediction/'):
for fn in os.listdir('test/'):
    if fn.endswith('result.csv.gz'):
        spl = fn.split('-')
        mir = spl[1] + spl[2] 
        df = pd.read_csv(f'mir_prediction/{fn}', compression = 'gzip', skiprows=13)
        
	# keep 7mer and 8mer site only
        sites = ['7mer', '8mer']
        df = df[df['Computational Predictions (method|target-site|start|stop|score)'].str.contains('|'.join(sites))]
        df.reset_index(inplace=True)
        
        df = df[[\
                'Input Gene',\
                'Transcript',\
                'Start',\
                'Stop',\
                'Computational Predictions (method|target-site|start|stop|score)'\
               ]].copy()
        
        df.rename(columns={'Input Gene': 'gene', 'Start': 'mir_start', 'Stop': 'mir_stop'}, inplace=True)
        
        outdf = pipeline(df) # run pipeline
        outname = f'OUTPUT/{mir}_out.csv'
        outdf.to_csv(outname, index=False, sep = '\t') # save df

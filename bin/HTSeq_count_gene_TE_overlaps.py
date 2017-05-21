import HTSeq
import collections
import argparse

def get_args():
    
    parser = argparse.ArgumentParser(description='A script to count stranded paired end reads that overlap between TEs and genes.')
    parser.add_argument("-i",
                        type=str,
                        default='',
                        dest='bam_file',
                        action='store',
                        help='Path to input bam file.'
                        )
    parser.add_argument("--gene_gtf",
                        type=str,
                        default='',
                        dest='gene_gtf',
                        action='store',
                        help='Path to input bam file.'
                        )
    parser.add_argument("--TE_gtf",
                        type=str,
                        default='',
                        dest='TE_gtf',
                        action='store',
                        help='Path to input bam file.'
                        )
    parser.add_argument("-o",
                        type=str,
                        default='',
                        dest='output_prefix',
                        action='store',
                        help='Output file prefix.'
                        )
    
    args = parser.parse_args()
    return args

def extract_GTF_features(file_path,feature_type='exon',attribute_label='gene_id'):
    gtf_file = HTSeq.GFF_Reader(file_path)
    output = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    
    for feature in gtf_file:
        if feature.type == feature_type:
            output[ feature.iv ] += feature.attr[attribute_label]
    return output

def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2

def output_counts(outfile,count_dict):
    with open(outfile,'w') as f:
        for key in count_dict:
            f.write('{key}\t{count}\n'.format(key=key,count=count_dict[key]))


if __name__=='__main__':
    args = get_args()
    
    #gene_gtf="/home/daniel/local_data/hg19/annotation/Homo_sapiens.GRCh37.87.gtf"
    #TE_gtf="/home/daniel/local_data/hg19/annotation/hg19_rmsk_TE_ensemblChrmapping.gtf"
    
    
    #gene_gtf="/home/daniel/local_data/hg19/annotation/hg19_refGene.gtf"
    attribute_label='gene_id'
    feature_type='exon'
    gene_features = extract_GTF_features(args.gene_gtf,feature_type,attribute_label)
    
    #TE_gtf="/home/daniel/local_data/hg19/annotation/hg19_rmsk_TE.gtf"
    #TE_gtf="/home/daniel/local_data/hg19/annotation/hg19_rmsk_TE_exon_filtered.gtf"
    #TE_gtf="/home/daniel/local_data/hg19/annotation/hg19_rmsk_TE_fort2014_CAGEseq.gtf"
    attribute_label='transcript_id'
    feature_type='exon'
    TE_features = extract_GTF_features(args.TE_gtf,feature_type,attribute_label)
    
    TE_first_counts = collections.Counter()
    TE_second_counts = collections.Counter()
    TE_only_counts = collections.Counter()
    
    #bam_file="/home/daniel/local_data/hipsci/star/test_bam_chr19_sorted.bam"
    almnt_file = HTSeq.BAM_Reader(args.bam_file)

    nUnmapped = 0
    nMultipleAlignments = 0
    nAttributeErrors = 0
    for bundle in HTSeq.pair_SAM_alignments(almnt_file, bundle=True):
        
        if len(bundle) != 1:
            nMultipleAlignments+=1
            continue  # Skip multiple alignments
    
        first_almnt, second_almnt = bundle[0]  # extract pair
        if not (first_almnt and second_almnt):
            nUnmapped+=1
            continue
        if not first_almnt.aligned and second_almnt.aligned:
            nUnmapped+=1
            continue
        
        #reverse stranded library
        try:
            first_almnt.iv =invert_strand(first_almnt.iv)
        except AttributeError:
            nAttributeErrors+=1
            continue
        
    
        #for read pairs with TE first in transcript
        # Note that the SECOND alignment is from the FIRST part of the transcript (ISR library)
        TE_ids_first = set()
        TE_ids_second = set()
        gene_ids_first = set()
        gene_ids_second = set()
        
        for iv, val in gene_features[first_almnt.iv].steps():
            gene_ids_second |= val
        for iv, val in TE_features[first_almnt.iv].steps():
            TE_ids_second |= val
        for iv, val in gene_features[second_almnt.iv].steps():
            gene_ids_first |= val
        for iv, val in TE_features[second_almnt.iv].steps():
            TE_ids_first |= val
    
        if (len(TE_ids_first)==1) and (len(gene_ids_second) == 1) and (len(TE_ids_second)==0):
            #TE is in the first part of the fragment, not in the second
            identifier = '_'.join(list(TE_ids_first)+list(gene_ids_second))
            TE_first_counts[ identifier ] += 1
        elif (len(TE_ids_first)==0) and (len(gene_ids_first) == 1) and (len(TE_ids_second)==1):
            #TE is in the second part of the fragment, not in the first
            identifier = '_'.join(list(TE_ids_second)+list(gene_ids_first))
            TE_second_counts[ identifier ] += 1
        elif (len(TE_ids_first&TE_ids_second)==1) and (len(gene_ids_first&gene_ids_second)==0):
            identifier = list(TE_ids_first)[0]
            TE_only_counts[identifier] += 1
    
    output_counts(args.output_prefix+'_TE_only.txt',TE_only_counts)
    output_counts(args.output_prefix+'_TE_first.txt',TE_first_counts)
    output_counts(args.output_prefix+'_TE_second.txt',TE_second_counts)
    print('{} attribute errors (iv is None)'.format(nAttributeErrors))
    print('{} multi-mapped reads discarded'.format(nMultipleAlignments))
    print('{} unmapped reads'.format(nUnmapped))


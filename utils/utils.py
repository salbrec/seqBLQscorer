import pandas as pd
import os
import subprocess
from copy import deepcopy

fastqc_values = {'FAIL':0, 'WARN':1, 'PASS':2}
bl_type_map = {'LM':1, 'HSR':2}

def getFileName(file_path):
	file_name = file_path[ -file_path[::-1].find('/') : ]
	if file_name.endswith('.gz'):
		file_name = file_name[:-3]
	if file_name.endswith('.fastq'):
		file_name = file_name[:-6]
	elif file_name.endswith('.fq'):
		file_name = file_name[:-3]
	return file_name

def get_file_length(fp):
    """
    Use the linux function to get the number of lines in a file.

    Args:
        fp (str): File path

    Returns:
        int: number of lines
    """
    wc = None
    try:
        call = 'wc -l %s'%(fp)
        res = subprocess.check_output(call, shell=True, text=True)
        wc = int(res.split()[0])
    except:
        print('Could not get the wc -l here:',fp)
    return wc


def read_Bowtie_stats(fp):
    """
    Reads a Bowtie2 report and returns the feature values in a dict.

    Args:
        fp (str): File path to FastQC summary.txt.

    Returns:
        dict: a map of MAP feature names and feature values.
        str: unparsed file content or "not exist" message.
    """
    if not os.path.exists(fp):
        return None, 'File does not exist.'
    lines = open(fp,'r').read().split('\n')
    if len(lines) != 7:
        return None, lines
    stats = {'total': int(lines[0].split()[0]) }
    stats['unpaired'] = int(lines[1].split()[0])
    stats['0times'] = int(lines[2].split()[0])
    stats['1time'] = int(lines[3].split()[0])
    stats['multi'] = int(lines[4].split()[0])
    stats['overall'] = stats['1time'] + stats['multi']
    percentages = {}
    for key in stats:
        if key != 'total':
            percentages['perc_'+key] = stats[key] / stats['total'] * 100.0
    stats.update(percentages)
    return stats, lines

def read_blocklist(bl_file):
    """
    Reads in a blocklist file and uses a certain format that allows for 
    using the regions efficiently in read counting procedure

    Args:
        bl_file (str): File path to the blocklist regions file.

    Returns:
        dict: dictionary containing lists of regions for each chromosome used as key
    """
    df_blocklist = pd.read_csv(bl_file, sep='\t', names=['chr','start','end','ID'])
    blocklist = {}
    for index, row in df_blocklist.iterrows():
        if not row['chr'] in blocklist:
            blocklist[row['chr']] = []
        bl_type = row['ID'].split('_')[0][2:]
        blocklist[row['chr']].append( (index+1, row['chr'], row['start'], row['end'],
            bl_type_map[bl_type], row['ID']) )
    return blocklist

def count_reads_in_regions(summits, regions, chrom_size_map, bl_mapping):
    """
    For each region, counting the summits within the region. 
    An overlap is only given if the summit is in the region, hence, 
    if more than half of the read overlaps with the blocklist region
    Args:
        summits (dict): dictionary with chromosome as key. The values
                        are lists of summits for the chromosome. The summits 
                        describe the center of the reads in this application.
        regions (dict): dictionary with chromosome as key. The values 
                        are lists of blocklist regions. 
        chrom_size_map (dict): dictionary with chromosome as key. The values
                                describe the chromosome length.
    Returns:
        pd.DataFrame: 
    """
    bincov = {'binID':[], 'chr':[], 'start':[], 'end':[], 'count':[],
        'blID':[], 'blType':[]}

    if bl_mapping:
        for chrom in chrom_size_map:
            if not chrom in regions:
                continue
            for binID, reg_chrom, reg_start, reg_end, reg_type, reg_ID in regions[chrom]:
                if reg_ID in summits:
                    bincov['binID'].append( binID ); bincov['chr'].append( reg_chrom );
                    bincov['start'].append( reg_start ); bincov['end'].append( reg_end );
                    bincov['count'].append( len(summits[reg_ID]) ); 
                    bincov['blType'].append( reg_type ); bincov['blID'].append( reg_ID );
        return pd.DataFrame(bincov)
    else:
        for chrom in chrom_size_map:
            if not chrom in summits or not chrom in regions:
                continue

            for binID, reg_chrom, reg_start, reg_end, reg_type, reg_ID in regions[chrom]:
                if chrom != reg_chrom:
                    print('!!! Something is WRONG here: ', reg_chrom, chrom)
                    return None
            
                count = len(list(filter( lambda x: x > reg_start and x < reg_end, summits[chrom] )))
            
                if count != 0:
                    bincov['binID'].append( binID ); bincov['chr'].append( reg_chrom );
                    bincov['start'].append( reg_start ); bincov['end'].append( reg_end );
                    bincov['count'].append( count ); bincov['blType'].append( reg_type );
                    bincov['blID'].append( reg_ID )
        return pd.DataFrame(bincov)




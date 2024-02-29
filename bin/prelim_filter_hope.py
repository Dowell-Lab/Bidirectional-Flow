import os
import sys
import math
from string import Template
import argparse


def bedtools_cov(prelim, bedgraph_pos, bedgraph_neg, sample_name, outdir):
    """Wrapper for running bedtools coverage
    Parameters
    ----------
    prelim : str (path)
        path to prelim bedgraph

    bedgraph_pos : str (path)
        path to positive bedgraph file corresponding to the input prelim
        
    bedgraph_neg : str (path)
        path to negative bedgraph file corresponding to the input prelim

    outbed : str (path)
        path to output bed file

    Returns
    -------
    Output from 'bedtools coverage'

    """
    # set a base name for the output files

    base_name = sample_name

    os.system(
        "bedtools coverage -a {} -b {} > {}/{}_prelim_bedtools_coverage_+.bed".format(
            prelim, bedgraph_pos, outdir, base_name
        )
    )
    
    os.system(
        "bedtools coverage -a {} -b {} > {}/{}_prelim_bedtools_coverage_-.bed".format(
            prelim, bedgraph_neg, outdir, base_name
        )
    )


def read_bedfile(bed, bed_type="tfit_prelim"):
    """Read input BED files and return a list of coordinates

    Parameters
    ----------
    bed : str (path)
        path to bed file

    bed_type : str
        determine source of bed file, either...
        ["tfit_prelim", "bedtools_coverage", "bed3", "bed4"]

    Returns
    -------
    bed_regions : list of list
        coordinates from BED file

    """

    bed_types = ["tfit_prelim", "bedtools_coverage", "bed3", "bed4"]
    bed_regions = []

    if bed_type not in bed_types:
        raise ValueError("Invalid BED file type. Expected one of: {}".format(bed_types))

    with open(bed) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            # get individual regions from bed file
            regions = line.split("\t")

            # bed file from tfit prelim
            if bed_type == "tfit_prelim" and len(regions) != 5:
                raise IndexError(
                    "Tfit prelim BED file should have 5 columns. Double check input."
                )

            # bedtools coverage bed file
            elif bed_type == "bedtools_coverage" and len(regions) < 7:
                raise IndexError(
                    "Bedtools coverage BED file should have 7 columns (or greater for BED4> files). Double check input."
                )

            # bed3 file format
            elif bed_type == "bed3" and len(regions) != 3:
                raise IndexError("BED3 file should have 3 columns. Double check input.")

            # bed4 file format
            elif bed_type == "bed4" and len(regions) != 4:
                raise IndexError("BED4 file should have 4 columns. Double check input.")
            else:
                bed_regions.append(regions)

    return bed_regions


# def bedtools_cov_filter(bedtools_cov_bed, filter_at=9):
#     """Filter BED file from 'bedtools coverage run' based on the 4th from last column

#     Parameter
#     ---------
#     bedtools_cov_bed : list of lists
#         bed regions from bedtools coverage in list format

#     filter_at : int (default = 9)
#         minimum coverage to use as cut-off

#     Returns
#     -------
#     regions_filtered : list of lists
#         filtered bed regions with coverage greater than 'filter_at'

#     """
#     regions_filtered = []

#     # get coverage from the bedtools coverage output file
#     for region in bedtools_cov_bed:

#         # coverage is in the 4th to last column
#         cov = int(region[-4])
#         if cov > filter_at:
#             regions_filtered.append([region[0], region[1], region[2], region[3]])
#         else:
#             pass

#     return regions_filtered


def break_large_regions(chrom, width, start, stop, region_id, break_by=10000):
    """break up regions based on region length
    Parameters
    ----------
    chrom : str
        chromosome id

    width : int
        width of the region

    start : int
        start coordinate of region

    stop : int
        stop coordinate of region

    region_id : str
        region name (ME_ for prelim regions)

    break_by : int
        region chucks to be generated (default = 10,000)

    Returns
    -------
    updated_regions : list of lists
        list with all coordinates


    """

    # determine how many regions per prelim
    num_regions = math.ceil(int(width) / break_by)

    # get the widths for the smaller regions
    break_width = round(width / num_regions)

    # initalize list for adding new regions
    updated_regions = []

    for i in range(num_regions):

        # create updated region id
        _region_id = "{}:{}".format(str(region_id), str(i))

        # determine start coordinate
        _start = int(start + (break_width * i))  # (break_by *i)

        # determine stop coordinates
        _stop = int(start + (break_width * (i + 1)))  # - 1 #(break_by *(i +1))) - 1

        # update the last stop coordinate to the regions stop coordinate
        if _stop < stop:
            _stop
        else:
            _stop = stop

        # add the update coordinates to a list
        updated_regions.append([chrom, _start, _stop, _region_id])

    return updated_regions


def break_large_regions_staggered(chrom, width, start, stop, region_id, break_by=10000):
    """break up regions based on region length in a overlapping manner
    Parameters
    ----------
    chrom : str
        chromosome id

    width : int
        width of the region

    start : int
        start coordinate of region

    stop : int
        stop coordinate of region

    region_id : str
        region name (ME_ for prelim regions)

    break_by : int
        region chucks to be generatated (default = 10,000)

    Returns
    -------
    updated_regions : list of lists
        list with all coordinates

    """
    # determine number of region breaks without staggering
    num_regions = math.ceil(width / break_by)

    # determine the breaks widths
    breaks_width = math.ceil(width / num_regions)

    # determine number of region breaks WITH staggering
    num_regions_exp = int(num_regions + (num_regions - 1))

    # how much to stagger the start coordinates
    # this will vary depending on the region widths
    step = math.ceil(width / num_regions_exp)

    # initialize start coordinates
    _start = [(step * i) + start for i in range(num_regions_exp)]

    # assign new region ids
    _region_id = [
        "{}:{}".format(str(region_id), str(i)) for i in range(num_regions_exp)
    ]

    # initalize list for adding new regions
    updated_regions = []

    for i, j in zip(_start, _region_id):

        _stop = int(i + breaks_width)

        if _stop < stop:
            _stop
        else:
            _stop = stop

        # add the update coordinates to a list
        updated_regions.append([chrom, int(i), _stop, j])

    return updated_regions


def dice_prelim(prelim_bed_list, staggered=True, break_by=10000):
    """load prelim regions and break larger region into smaller parts
    determined by the 'break_by' parameter

    Parameters
    ----------
    prelim_bed_list : list of list
        list with list of bed file coordinates

    break_by : int (10000 default)
        break regions greater than this value

    Returns
    -------
    all_regions : list of lists
        containing coordinates smaller or equalt to specified length (BED4 format)

    """

    all_regions = []

    for region in prelim_bed_list:

        # split them by columns
        chrm = region[0]
        start = int(region[1])
        stop = int(region[2])
        region_id = region[3]
        # region_summary = region[4]

        # get region widths
        region_widths = stop - start

        # if regions greater than 'break_by', break them up
        # to smaller regions
        if int(region_widths) < int(break_by):
            all_regions.append([chrm, start, stop, region_id])
        else:
            if staggered == True:
                new_regions = break_large_regions_staggered(
                    chrm, region_widths, start, stop, region_id, break_by=int(break_by)
                )
                for smaller_regions in new_regions:
                    all_regions.append(smaller_regions)
            else:
                new_regions = break_large_regions(
                    chrm, region_widths, start, stop, region_id, break_by=int(break_by)
                )
                for smaller_regions in new_regions:
                    all_regions.append(smaller_regions)

    return all_regions

def parse_gtf_tss(tss_gtf, chrom_sizes, outdir, tss_width=2000, add_slop=True, slop=750):
    """load in gene gtf file and extract TSS ends into regions

    Parameters
    ----------
    tss_gtf : str (path)
        path to gene annotation file (in gtf format)

    chrom_sizes : str (path)
        path to chromosome sizes file

    outdir : str (path)
        output directory for TSS bedfiles

    tss_width : int (2000 default)
        width of region around TSS coordinate

    add_slop : boolean (True default)
        whether to add bases to TSS regions for merging

    slop : int (750 default)
        amount to add to TSS regions

    Returns
    -------
    Bed file(s) containing TSS regions and slopped regions, if specified

    Path to TSS bedfile
    """

    # Load in gtf file and get 5' coordinate of gene
    tss_coord = []

    with open(tss_gtf) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            # get individual regions from bed file
            regions = line.split("\t")

            # check if gtf file
            if len(regions) != 9:
                raise IndexError(
                    "GTF file should have 9 columns. Double check input."
                )
            # Take first coordinate if on + strand
            elif regions[6] == "+":
                meta = regions[8].split("\"")
                region_to_add = [regions[0], (int(regions[3]) - 1), regions[3], meta[1], "0,0,0"]
                if region_to_add not in tss_coord:
                    tss_coord.append([regions[0], (int(regions[3]) - 1), regions[3], meta[1], "0,0,0"])
            # Take second coordinate if on - strand
            elif regions[6] == "-":
                meta = regions[8].split("\"")
                region_to_add = [regions[0], (int(regions[4]) - 1), regions[4], meta[1], "0,0,0"]
                if region_to_add not in tss_coord:
                    tss_coord.append([regions[0], (int(regions[4]) - 1), regions[4], meta[1], "0,0,0"])
            else:
                raise IndexError(
                    "GTF file should have strand info in 7th column. Double check input."
                )

    # Write out base 5' coordinate bedfile in bed5 format
    write_bedfile(str(outdir) + "/tss_single_base.bed", tss_coord, 5)

    # Generate filepaths
    tss_bed = str(outdir) + "/tss_" + str(tss_width) + "bp.bed"
    if add_slop == True:
        tss_slop_bed = tss_bed[0:-4] + "_" + str(slop) + "overlap.bed"
    else:
        tss_slop_bed = ""

    # Run bedtools slop to get base TSS regions
    os.system(
        "bedtools slop -i {} -g {} -b {} | sort -k1,1 -k2,2n | bedtools merge -c 4,5 -o collapse,distinct > {}".format(
            (outdir + "/tss_single_base.bed"),
            chrom_sizes,
            math.ceil(tss_width/2),
            tss_bed
        )
    )

    # If adding slop for merging, run again with extra slop
    if add_slop == True:
        os.system(
            "bedtools slop -i {} -g {} -b {} | sort -k1,1 -k2,2n | bedtools merge -c 4,5 -o collapse,distinct > {}".format(
                (outdir + "/tss_single_base.bed"),
                chrom_sizes,
                (math.ceil(tss_width/2) + int(slop)),
                tss_slop_bed
            )
        )

    return tss_bed, tss_slop_bed

def slop_mu(mu_bed, chrom_sizes, outdir, mu_width=2000, add_slop=False, slop=750):
    """load in gene gtf file and extract TSS ends into regions

    Parameters
    ----------
    mu_bed : str (path)
        path to bed file with mu values of Bids want to include in prelim

    chrom_sizes : str (path)
        path to chromosome sizes file

    outdir : str (path)
        output directory for TSS bedfiles

    mu_width : int (600 default)
        width of region around mu coordinate

    add_slop : boolean (True default)
        whether to add bases to mu regions for merging

    slop : int (750 default)
        amount to add to mu regions

    Returns
    -------
    Bed file(s) containing mu regions and slopped regions, if specified

    Path to mu bedfile
    """
    # get names of files ready

    # Generate filepaths
    mu_width_bed = str(outdir) + "/mu_" + str(mu_width) + "bp.bed"
    if add_slop == True:
        mu_slop_bed = mu_width_bed[0:-4] + "_" + str(slop) + "overlap.bed"
    else:
        mu_slop_bed = ""

    # Run bedtools slop to get base mu regions
    os.system(
        "bedtools slop -i {} -g {} -b {} | sort -k1,1 -k2,2n | bedtools merge -c 4,5 -o collapse,distinct > {}".format(
            mu_bed,
            chrom_sizes,
            math.ceil(mu_width/2),
            mu_width_bed
        )
    )

    # If adding slop for merging, run again with extra slop
    if add_slop == True:
        os.system(
            "bedtools slop -i {} -g {} -b {} | sort -k1,1 -k2,2n | bedtools merge -c 4,5 -o collapse,distinct > {}".format(
                mu_bed,
                chrom_sizes,
                (math.ceil(mu_width/2) + int(slop)),
                mu_slop_bed
            )
        )

    return mu_width_bed, mu_slop_bed


def designate_tss(prelim, tss_bed, sample_name, outdir, tss_slop_bed, add_slop=True):

    """Wrapper for running bedtools to subtract TSS regions from 
    prelim and merge prelim and TSS regions

    Parameters
    ----------
    prelim : str (path)
        path to prelim bedgraph

    tss_bed : str (path)
        path to bed file containing TSS regions

    sample_name : str
        basename for output files

    outdir : str (path)
        path to output directory

    tss_slop_bed : str (path)
        path to bed file with added overlap for merging

    add_slop : boolean (True default)
        whether to add bases to TSS regions for merging

    Returns
    -------
    Sorted prelim file with TSS regions added or clearly defined 
    within previous prelim regions

    Path to prelim file

    """
    # set a base name for the output files and set paths
    base_name = str(sample_name)
    prelim_filepath = str(outdir) + "/" + base_name + "_prelim_tss.bed"
    prelim_subtract = str(outdir) + "/" + base_name + "_prelim_tss_subtract.bed"

    os.system(
        "bedtools subtract -a {} -b {} > {}".format(
            prelim, tss_bed, prelim_subtract
        )
    )

    # rename regions in subtracted file to remain distinct
    region_names = {}
    regions_to_write = []

    with open(prelim_subtract) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            # get individual regions from bed file
            regions = line.split("\t")

            if regions[3] in region_names.keys():
                region_names[regions[3]] = int(region_names[regions[3]]) + 1
                regions[3] = str(regions[3]) + ":" + str(region_names[regions[3]])
                regions_to_write.append(regions)
            else:
                region_names[regions[3]] = 0
                regions_to_write.append(regions)

    write_bedfile(prelim_subtract, regions_to_write, 5)
                
    # Concatenate and sort prelim and tss files
    if add_slop == False:
        os.system(
            "cat {} {} | bedtools sort > {}".format(
                prelim_subtract, tss_bed, prelim_filepath
            )
        )
    else:
        os.system(
            "cat {} {} | bedtools sort > {}".format(
                prelim_subtract, tss_slop_bed, prelim_filepath
            )
        )

    return prelim_filepath

def designate_mu(prelim, mu_width_bed, sample_name, outdir, mu_slop_bed, add_slop=True):

    """Wrapper for running bedtools to subtract mu regions from the original
    prelim file and merge prelim and mu regions

    Parameters
    ----------
    prelim : str (path)
        path to prelim bedgraph

    mu_width_bed : str (path)
        path to bed file containing mu regions (NOT single base resolution)

    sample_name : str
        basename for output files

    outdir : str (path)
        path to output directory

    mu_slop_bed : str (path)
        path to bed file with added overlap for merging

    add_slop : boolean (True default)
        whether to add bases to mu regions for merging

    Returns
    -------
    Sorted prelim file with mu regions added or clearly defined 
    within previous prelim regions

    Path to prelim file

    """
    # set a base name for the output files and set paths
    base_name = str(sample_name)
    prelim_filepath = str(outdir) + "/" + base_name + "_prelim_mu.bed"
    prelim_subtract = str(outdir) + "/" + base_name + "_prelim_mu_subtract.bed"

    os.system(
        "bedtools subtract -a {} -b {} > {}".format(
            prelim, mu_width_bed, prelim_subtract
        )
    )

    # rename regions in subtracted file to remain distinct
    region_names = {}
    regions_to_write = []

    with open(prelim_subtract) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            # get individual regions from bed file
            regions = line.split("\t")

            if regions[3] in region_names.keys():
                region_names[regions[3]] = int(region_names[regions[3]]) + 1
                regions[3] = str(regions[3]) + ":" + str(region_names[regions[3]])
                regions_to_write.append(regions)
            else:
                region_names[regions[3]] = 0
                regions_to_write.append(regions)

    write_bedfile(prelim_subtract, regions_to_write, 5)
                
    # Concatenate and sort prelim and tss files
    if add_slop == False:
        os.system(
            "cat {} {} | bedtools sort > {}".format(
                prelim_subtract, mu_width_bed, prelim_filepath
            )
        )
    else:
        os.system(
            "cat {} {} | bedtools sort > {}".format(
                prelim_subtract, mu_slop_bed, prelim_filepath
            )
        )

    return prelim_filepath


def write_bedfile(region_prelim, region_list, columns=4):
    """Write output from the input list as a bedfile
    Parameters
    ----------
    region_prelim : str (file name)
        name of file where bed regions will be written

    region_list : list of lists
        containing regions in BED4 or BED5 format

    columns: int (default = 4)

    Returns
    -------
    writes 'region_prelim' in BED4 or BED5 format
    """

    with open(region_prelim, "w") as output:
        for prelims in region_list:
            if int(columns) == 4:
                output.write(
                    str(
                        str(prelims[0])
                        + "\t"
                        + str(prelims[1])
                        + "\t"
                        + str(prelims[2])
                        + "\t"
                        + str(prelims[3])
                        + "\n"
                    )
                )
            elif int(columns) == 5:
                output.write(
                    str(
                        str(prelims[0])
                        + "\t"
                        + str(prelims[1])
                        + "\t"
                        + str(prelims[2])
                        + "\t"
                        + str(prelims[3])
                        + "\t"
                        + str(prelims[4])
                        + "\n"
                    )
                )
            else:
                raise IndexError(
                    "Number of columns not supported for writing. Check input bedfile."
                )

        output.close()


def main(
    prelim_bed,
    prelim_bedgraph,
    output_directory,
    sample_name,
    mu_bed,
    tss_gtf="",
    chrom_sizes="",
    bed_type="tfit_prelim",
    filter_at=9,
    break_by=10000,
    staggered=True,
    parse_tss=True,
    add_mu=True,
    tss_width=2000,
    mu_width=2000,
    add_slop=True,
    tss_slop=750
):
    """Process input prelim files"""
    print("MU_BED", mu_bed)
    # 1 : set a base name for the output files
    if sample_name:
        base_name = sample_name
    else:
        base = os.path.basename(prelim_bed)
        base_name = base[: base.index(".")]
        
    # 2A : Add mu regions and point to new prelim bed file
    if add_mu == True:
        print("adding mu")
        # add mu_width/2 to both sides of mu (default 1kb --> 2kb width)
        mu_width_bed, mu_slop_bed = slop_mu(
            mu_bed, chrom_sizes, output_directory, mu_width, False, tss_slop)
        # replace orig Tfit prelim calls with the mu calls if they overlap
        prelim_filepath = designate_mu(
            prelim_bed, mu_width_bed, base_name, output_directory, mu_slop_bed, False)
        prelim_bed = prelim_filepath

    # 2B : parse TSS regions and point to new prelim bed file 
    if parse_tss == True:
        print("parsing TSS with prelim file", prelim_bed)
        tss_bed, tss_slop_bed = parse_gtf_tss(
            tss_gtf, chrom_sizes, output_directory, tss_width, add_slop, tss_slop)
        # replace orig Tfit prelim OR mu calls with TSS if they overlap
        prelim_filepath = designate_tss(
            prelim_bed, tss_bed, base_name, output_directory, tss_slop_bed, add_slop)
        prelim_bed = prelim_filepath
        

    # 3 : load the prelim bed file
    prelim_bed_list = read_bedfile(prelim_bed, bed_type="tfit_prelim")

    # 4: break regions
    prelim_diced_list = dice_prelim(
        prelim_bed_list, staggered=staggered, break_by=break_by)

    # 5 : save the intermediate diced bed file
    write_bedfile(
        "{}/{}_prelim_diced_intermediate.bed".format(output_directory, base_name),
        prelim_diced_list,
        4)
    
    # 6. Get the positive and negative bedgraphs
    os.system(
        "grep '-' {} | sed -r 's/-//' > {}/{}_neg.bedGraph".format(prelim_bedgraph, output_directory, base_name
        ))
    os.system(
        "grep -v '-' {} > {}/{}_pos.bedGraph".format(prelim_bedgraph, output_directory, base_name
        ))
    
    # 7. Get the coverage for both strands
    os.system(
        "bedtools coverage -a {}/{}_prelim_diced_intermediate.bed -b {}/{}_neg.bedGraph > {}/{}_cov_neg.bed".format(output_directory, base_name, output_directory, base_name, output_directory, base_name)
    )
    os.system(
        "bedtools coverage -a {}/{}_prelim_diced_intermediate.bed -b {}/{}_pos.bedGraph > {}/{}_cov_pos.bed".format(output_directory, base_name, output_directory, base_name, output_directory, base_name)
    )
    
    # 8 : Filter regions with low coverage on BOTH strands
    os.system(Template("awk -v FS='\t' -v OFS='\t' '{ if ($e > 5) print $a,$b,$c,$d }' $f/${g}_cov_neg.bed > $f/${g}_cov_filt_neg.bed").safe_substitute(a="$1", b="$2", c="$3", d="$4", e="$5", f=output_directory, g=base_name))
    os.system(Template("awk -v FS='\t' -v OFS='\t' '{ if ($e > 5) print $a,$b,$c,$d }' $f/${g}_cov_pos.bed > $f/${g}_cov_filt_pos.bed").safe_substitute(a="$1", b="$2", c="$3", d="$4", e="$5", f=output_directory, g=base_name)
             )
# #     os.system("awk -v FS='\t' -v OFS='\t' '{ if ($5 > 9) print $1,$2,$3,$4 }' {}/{}_cov_pos.bed > {}/{}_cov_filt_pos.bed".format(output_directory, base_name, output_directory, base_name)
# #              )
     
    # 9 : Only keep regions with signfiicant coverage on BOTH strands
    os.system("bedtools intersect -a {}/{}_cov_filt_neg.bed -b {}/{}_cov_filt_pos.bed -u -f 0.97 -r > {}/{}_prelim_coverage_filtered_diced.bed".format(output_directory, base_name, output_directory, base_name, output_directory, base_name)
             )
              
#### Original ways of doing things        
#          # 6 : get coverage over prelim regions
#     bedtools_cov(
#         "{}/{}_prelim_diced_intermediate.bed".format(output_directory, base_name),
#         bedgraph,
#         base_name,
#         output_directory,
#     )

     # 7 : load diced region
#     prelim_coverage_list = read_bedfile(
#         "{}/{}_prelim_bedtools_coverage.bed".format(output_directory, base_name),
#         bed_type="bedtools_coverage",
#     )

    # 8 : filter regions with low coverage
#     prelim_filtered_list = bedtools_cov_filter(
#         prelim_coverage_list, filter_at=filter_at
#     )

    # 9 : write the new filtered file in bed format
#     write_bedfile(
#         "{}/{}_prelim_coverage_filtered_diced.bed".format(output_directory, base_name),
#         prelim_filtered_list,
#         4,
#     )


parser = argparse.ArgumentParser(description="Filter prelim BED files from Tfit")

parser.add_argument(
    "-p", "--prelim_bed", dest="prelm", help="Input prelim BED file", metavar="FILE"
)
parser.add_argument(
    "-b",
    "--prelim_bedgraph",
    dest="bedgraph",
    help="BedGraph file for sample",
    metavar="FILE",
)
parser.add_argument(
    "-o", "--outdirectory", dest="outdir", help="Directory for output", metavar="DIR"
)
parser.add_argument(
    "-t",
    "--bedtype",
    dest="btype",
    type=str,
    help="Type of bed file",
    default="tfit_prelim",
    metavar="STR",
)
parser.add_argument(
    "-f",
    "--cov_filt_by",
    dest="filt_cov",
    type=int,
    default=9,
    help="Minimum read coverage for regions",
    metavar="INT",
)
parser.add_argument(
    "-k",
    "--region",
    dest="brk",
    type=int,
    default=10000,
    help="Break regions into -k smaller segments",
    metavar="INT",
)

parser.add_argument(
    "-s",
    "--samp_id",
    dest="samp",
    type=str,
    help="Base name for the input sample",
    metavar="STR",
)
parser.add_argument(
    "-g",
    "--gtf_file",
    dest="gtf",
    type=str,
    help="Gene annotion GTF file",
    metavar="FILE",
)
parser.add_argument(
    "-m",
    "--mu_file",
    dest="mu_file",
    type=str,
    help="Bed file of single base resolution of mus to include in prelims",
    metavar="FILE",
)
parser.add_argument(
    "-c",
    "-chrom_file",
    dest="chrom",
    type=str,
    help="Chromosome sizes file",
    metavar="FILE",
)
parser.add_argument(
    "-w",
    "--tss_width",
    dest="tsswidth",
    type=int,
    default=2000,
    help="Width of region around TSS",
    metavar="INT",
)
parser.add_argument(
    "-mw",
    "--mu_width",
    dest="muwidth",
    type=int,
    default=2000,
    help="Width of region around TSS",
    metavar="INT",
)
parser.add_argument(
    "-l",
    "--slop_val",
    dest="slopval",
    type=int,
    default=750,
    help="Amount of extra overlap of TSS regions",
    metavar="INT",
)
parser.add_argument(
    "--stagger",
    dest="stagger",
    action="store_true",
    help="Stagger the diced regions into -k smaller segments",
)
parser.add_argument(
    "--no-stagger",
    dest="stagger",
    action="store_false",
    help="No stagger, regions broken in to -k smaller segments",
)
parser.add_argument(
    "--tss",
    dest="tss",
    action="store_true",
    help="Parse prelim with TSS file to separate out TSS regions",
)
parser.add_argument(
    "--no-tss",
    dest="tss",
    action="store_false",
    help="No parsing prelim file along with TSS file",
)

parser.add_argument(
    "--mu",
    dest="mu",
    action="store_true",
    help="Parse prelim file along with mu file",
)

parser.add_argument(
    "--no-mu",
    dest="mu",
    action="store_false",
    help="No parsing prelim file along with mu file",
)
parser.add_argument(
    "--tss-slop",
    dest="slop",
    action="store_true",
    help="Add slop around TSS before region concatenation",
)
parser.add_argument(
    "--no-tss-slop",
    dest="slop",
    action="store_false",
    help="No slop around TSS regions",
)

parser.set_defaults(stagger=True, tss=True, slop=True, mu=True)

args = parser.parse_args()

# main.run(args.freq, args.exp, args.outdir, args.len, args.num, args.sd)
main(
    args.prelm,
    args.bedgraph,
    args.outdir,
    args.samp,
    args.mu_file,
    args.gtf,
    args.chrom,
    args.btype,
    args.filt_cov,
    args.brk,
    args.stagger,
    args.tss,
    args.mu,
    args.tsswidth,
    args.muwidth,
    args.slop,
    args.slopval
)

import os
import sys
import math

import argparse


def bedtools_cov(prelim, bedgraph, sample_name, outdir):
    """Wrapper for running bedtools coverage
    Parameters
    ----------
    prelim : str (path)
        path to prelim bedgraph

    bedgraph : str (path)
        path to bedgraph file corresponding to the input prelim

    outbed : str (path)
        path to output bed file

    Returns
    -------
    Output from 'bedtools coverage'

    """
    # set a base name for the output files

    base_name = sample_name

    os.system(
        "bedtools coverage -a {} -b {} > {}/{}_prelim_bedtools_coverage.bed".format(
            prelim, bedgraph, outdir, base_name
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


def bedtools_cov_filter(bedtools_cov_bed, filter_at=9):
    """Filter BED file from 'bedtools coverage run' based on the 4th from last column

    Parameter
    ---------
    bedtools_cov_bed : list of lists
        bed regions from bedtools coverage in list format

    filter_at : int (default = 9)
        minimum coverage to use as cut-off

    Returns
    -------
    regions_filtered : list of lists
        filtered bed regions with coverage greater than 'filter_at'

    """
    regions_filtered = []

    # get coverage from the bedtools coverage output file
    for region in bedtools_cov_bed:

        # coverage is in the 4th to last column
        cov = int(region[-4])
        if cov > filter_at:
            regions_filtered.append([region[0], region[1], region[2], region[3]])
        else:
            pass

    return regions_filtered


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
        region chucks to be generatated (default = 10,000)

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
        _start = start + (break_width * i)  # (break_by *i)

        # determine stop coordinates
        _stop = start + (break_width * (i + 1))  # - 1 #(break_by *(i +1))) - 1

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
    num_regions_exp = num_regions + (num_regions - 1)

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

        _stop = i + breaks_width

        if _stop < stop:
            _stop
        else:
            _stop = stop

        # add the update coordinates to a list
        updated_regions.append([chrom, i, _stop, j])

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


def write_bedfile(region_prelim, region_list):
    """Write output from the input list as a bedfile
    Parameters
    ----------
    region_prelim : str (file name)
        name of file where bed regions will be written

    region_list : list of lists
        containing regions in BED4 format

    Returns
    -------
    writes 'region_prelim' in BED4 format
    """

    with open(region_prelim, "w") as output:
        for prelims in region_list:

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

        output.close()


def main(
    prelim_bed,
    prelim_bedgraph,
    output_directory,
    sample_name,
    bed_type="tfit_prelim",
    filter_at=9,
    break_by=10000,
    staggered=True,
):
    """Process input prelim files"""

    # 1 : set a base name for the output files
    if sample_name:
        base_name = sample_name
    else:
        base = os.path.basename(prelim_bed)
        base_name = base[: base.index(".")]

    # 1a : add TSS regions here?

    # 1b : load the new prelim + TSS regions

    # 2 : load the prelim bed file
    prelim_bed_list = read_bedfile(prelim_bed, bed_type="tfit_prelim")

    # 4: break regions
    prelim_diced_list = dice_prelim(
        prelim_bed_list, staggered=staggered, break_by=break_by
    )

    # print(prelim_diced_list)

    # 5 : save the intermediate diced bed file
    write_bedfile(
        "{}/{}_prelim_diced_intermediate.bed".format(output_directory, base_name),
        prelim_diced_list,
    )
    

    # 6 : get coverage over prelim regions
    bedtools_cov(
        "{}/{}_prelim_diced_intermediate.bed".format(output_directory, base_name),
        prelim_bedgraph,
        base_name,
        output_directory,
    )

    # 7 : load diced region
    prelim_coverage_list = read_bedfile(
        "{}/{}_prelim_bedtools_coverage.bed".format(output_directory, base_name),
        bed_type="bedtools_coverage",
    )

    # 8 : filter regions with low coverage
    prelim_filtered_list = bedtools_cov_filter(
        prelim_coverage_list, filter_at=filter_at
    )

    # 9 : write the new filtered file in bed format
    write_bedfile(
        "{}/{}_prelim_coverage_filtered_diced.bed".format(output_directory, base_name),
        prelim_filtered_list,
    )


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
parser.set_defaults(stagger=True)

args = parser.parse_args()


# main.run(args.freq, args.exp, args.outdir, args.len, args.num, args.sd)
main(
    args.prelm,
    args.bedgraph,
    args.outdir,
    args.samp,
    args.btype,
    args.filt_cov,
    args.brk,
    args.stagger,
)

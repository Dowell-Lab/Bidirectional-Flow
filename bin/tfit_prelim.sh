#!/bin/bash

##help message for for the script
usage() { echo "Usage: $0 [-t <tfit paths>] [-c <tfit config>] [-b <bedgraph file>] [-p <prefix>] [-n <threads>]" 1>&2; exit 1; }


while getopts ":t:c:b:p:n:" arg; do
    case "${arg}" in
        t)
            t=${OPTARG}
            ;;
        c)
            c=${OPTARG}
            ;;
	b)
	    b=${OPTARG}
	    ;;
	p)
	    p=${OPTARG}
	    ;;
	n)
	    n=${OPTARG}
	    ;;

        h | *) ##display help message
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${t}" ] || [ -z "${c}" ] ; then
    usage
fi

##load modules to run tfit

module load openmpi/1.6.4
module load gcc/7.1.0
module load bedtools/2.28.0

##print and export parameters
export "OMP_NUM_THREADS=1"
echo "Number of threads: ${n}"
echo "Run on           : ${SLURM_JOB_NODELIST}"

echo "Tfit path        : ${t}"
echo "Tfit config      : ${c}"
echo "Bedgraph path    : ${b}"
echo "Output prefix    : ${p}"
echo "Output directory : $(pwd)"

##run tfit!
mpirun ${t} bidir -config ${c} -ij $(pwd)/${b} -N ${p} -o $(pwd)

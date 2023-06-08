#!/bin/bash

I_SRC=""
O_SRC="src"

I_OUTPUT=""
O_OUTPUT="out"

X_KMER_SIZE=11
I_KMER_SIZE=${X_KMER_SIZE}
O_KMER_SIZE="kmer_size"

X_THREADS=1
I_THREADS=${X_THREADS}
O_THREADS="num_threads"

X_MODEL_PATH=""
I_MODEL_PATH=""
O_MODEL_PATH="model_path"

X_METH_PROB=0.5
I_METH_PROB=${X_METH_PROB}
O_METH_PROB="meth_prob"

X_UNMETH_PROB=0.5
I_UNMETH_PROB=${X_UNMETH_PROB}
O_UNMETH_PROB="unmeth_prob"

X_MOTIF="CG"
I_MOTIF=${X_MOTIF}
O_MOTIF="motif"

X_MAPQ=10
I_MAPQ=${X_MAPQ}
O_MAPQ="mapping_quality"

X_REFERENCE=""
I_REFERENCE=""
O_REFERENCE="reference"

X_PHASE=0
I_PHASE=${X_PHASE}
O_PHASE="phase"

X_BAM=""
I_BAM=""
O_BAM="bam"

function parse_arguments()
{
	ARGS=`getopt --long ${O_SRC}:,${O_OUTPUT}:,${O_KMER_SIZE}:,${O_THREADS}:,${O_MODEL_PATH}:,${O_METH_PROB}:,${O_UNMETH_PROB}:,${O_MAPQ}:,${O_REFERENCE}:,${O_MOTIF}:,${O_PHASE},${O_BAM}: -- "$@"`
	if [ $? -ne 0 ]; then
		echo [$(date)] parsing argument fail
		exit 1
	fi

	local verbose=1

	if [ ${verbose} -eq 1 ]; then
		echo ARGS=[$ARGS]
	fi
	eval set -- "${ARGS}"
	if [ ${verbose} -eq 1 ]; then
		echo formatted params=[$@]
	fi

	while true
	do
		case "$1" in
			--${O_SRC})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_SRC}], $2
				fi
				I_SRC=$2
				shift 2
				;;
			--${O_OUTPUT})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_OUTPUT}], $2
				fi
				I_OUTPUT=$2
				shift 2
				;;
			--${O_KMER_SIZE})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_KMER_SIZE}], $2
				fi
				I_KMER_SIZE=$2
				shift 2
				;;
			--${O_THREADS})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_THREADS}], $2
				fi
				I_THREADS=$2
				shift 2
				;;
			--${O_MODEL_PATH})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_MODEL_PATH}], $2
				fi
				I_MODEL_PATH=$2
				shift 2
				;;
			--${O_METH_PROB})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_METH_PROB}], $2
				fi
				I_METH_PROB=$2
				shift 2
				;;
			--${O_UNMETH_PROB})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_UNMETH_PROB}], $2
				fi
				I_UNMETH_PROB=$2
				shift 2
				;;
			--${O_MOTIF})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_MOTIF}], $2
				fi
				I_MOTIF=$2
				shift 2
				;;
			--${O_MAPQ})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_MAPQ}], $2
				fi
				I_MAPQ=$2
				shift 2
				;;
			--${O_REFERENCE})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_REFERENCE}], $2
				fi
				I_REFERENCE=$2
				shift 2
				;;
			--${O_PHASE})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_PHASE}], $2
				fi
				I_PHASE=1
				shift 1
				;;
			--${O_BAM})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_BAM}], $2
				fi
				I_BAM=$2
				shift 2
				;;
			--)
				shift
				break
				;;
			*)
				echo [$(date)] fail at parsing arguments
				exit 1
				;;
		esac
	done
}

function dump_usage()
{
	pn=$1
	echo ""
	echo "USAGE:"
	echo "${pn} [--${O_KMER_SIZE}] [--${O_THREADS}] [--${O_METH_PROB}] [--${O_UNMETH_PROB}]" 
	echo "  <--${O_SRC}> <--${O_OUTPUT}> <--${O_MODEL_PATH}>"

	echo ""
	echo "OPTION DESCRIPTIONS:"
	echo "--${O_KMER_SIZE} <integer, >0>"
	echo "    Kmer size around motif loci"
	echo "    Default = '${X_KMER_SIZE}'"
	echo "--${O_THREADS} <integer, >=0>"
	echo "    Number of CPU threads"
	echo "    Default = '${X_THREADS}'"
	echo "--${O_MODEL_PATH}"
	echo "    Path to trainned model"
	echo "--${O_SRC} <path>"
	echo "    Path to 5mc-gnn"
	echo "--${O_OUTPUT} <path>"
	echo "    Output directory"
    echo "--${O_METH_PROB} <real>\n"
    echo "    Probability cutoff for methylation"
    echo "    Default = '${X_METH_PROB}'"
    echo "--${O_UNMETH_PROB} <real>\n"
    echo "    Probability cutoff for unmethylation"
    echo "    Default = '${X_UNMETH_PROB}'"
	echo "--${O_MOTIF} <string>"
	echo "    The motif sequence"
	echo "    Default = '${X_MOTIF}'"
	echo "--${O_MAPQ} <integer>"
	echo "    Minimum mapping quality"
	echo "    Default = '${X_MAPQ}'"
	echo "--${O_REFERENCE} <fa>"
	echo "    Reference path"
	echo "--${O_PHASE}"
	echo "    Also dump paternal and maternal frequencies"
	echo "--${O_BAM} <path>"
	echo "    Path to BAM file"
}

function validate_arguments()
{
	if [ "${I_SRC}"x = x ]; then
		echo "5mc directory is not provided"
		return 1
	fi

	if [ "${I_OUTPUT}"x = x ]; then
		echo "output directory is not provided"
		return 1
	fi

	if [ "${I_REFERENCE}"x = x ]; then
		echo "reference is not provided"
		return 1
	fi

	if [ "${I_BAM}"x = x ]; then 
		echo "BAM file is not provided"
		return 1
	fi

	return 0
}

function dump_params()
{
	echo ""
	echo "========================== Params"
	echo "5mc-gnn path: ${I_SRC}"
	echo "output directory: ${I_OUTPUT}"
	echo "kmer size: ${I_KMER_SIZE}"
	echo "CPU threads: ${I_THREADS}"
    echo "model: ${I_MODEL_PATH}"
    echo "meth prob: ${I_METH_PROB}"
    echo "unmeth prob: ${I_UNMETH_PROB}"
	echo "motif: ${I_MOTIF}"
	echo "mapping quality: ${I_MAPQ}"
	echo "reference: ${I_REFERENCE}"
	echo "phase: ${I_PHASE}"
	echo "BAM: ${I_BAM}"
	echo ""
}

function run_commands()
{
    cmd=$*
    echo [$(date)] Running command
    echo "  ${cmd}"
    ${cmd}
    if [ $? -ne 0 ]; then
        echo [$(date)] Fail at running
        echo "  ${cmd}"
        exit 1
    fi
}

parse_arguments $0 $*

validate_arguments
if [ $? -ne 0 ]; then
	dump_usage $0
	exit 1
fi

dump_params

mkdir -p ${I_OUTPUT}

###############

function extract_features()
{
    local job_is_done=${I_OUTPUT}/extract-features.done
    if [ -f ${job_is_done} ]; then
        return 
    fi 

	PHASE_OPT=""
	if [ ${I_PHASE} -eq 1 ]; then
		PHASE_OPT="-p"
	fi
	cmd="samtools view -@8 ${I_BAM} | ${I_SRC}/src/cpp_tools/extract_freq_features_from_sam -q ${I_MAPQ} -M ${I_MOTIF} -s ${I_KMER_SIZE} -m ${I_METH_PROB} -u ${I_UNMETH_PROB} ${PHASE_OPT} -t ${I_THREADS} - ${I_REFERENCE} ${I_OUTPUT}"
	echo "Running command"
	echo "  ${cmd}"
	samtools view -@8 ${I_BAM} | ${I_SRC}/src/cpp_tools/extract_freq_features_from_sam -q ${I_MAPQ} -M ${I_MOTIF} -s ${I_KMER_SIZE} -m ${I_METH_PROB} -u ${I_UNMETH_PROB} ${PHASE_OPT} -t ${I_THREADS} - ${I_REFERENCE} ${I_OUTPUT}
    if [ $? -ne 0 ]; then
        echo "fail"
        exit 1
    fi 

    touch ${job_is_done}
}

function call_freq()
{
    local job_is_done=${I_OUTPUT}/freq-call-model.done
    if [ -f ${job_is_done} ]; then
        return 
    fi

	features=${I_OUTPUT}/freq-call.model.features
    freq_call=${I_OUTPUT}/freq-call.model.txt
    cmd="python ${I_SRC}/src/lstm-freq/5mc_freq.py --model_path ${I_MODEL_PATH} --kmer_size ${I_KMER_SIZE} --input ${features} --out ${freq_call} --dl_workers ${I_THREADS}"
    echo "Running command"
    echo "  ${cmd}"
    python ${I_SRC}/src/lstm-freq/5mc_freq.py --model_path ${I_MODEL_PATH} --kmer_size ${I_KMER_SIZE} --input ${features} --out ${freq_call} --dl_workers ${I_THREADS}
    if [ $? -ne 0 ]; then
        echo "fail"
        exit 1
    fi 

	if [ ${I_PHASE} -eq 1 ]; then
		features=${I_OUTPUT}/freq-call.hp1.model.features
    	freq_call=${I_OUTPUT}/freq-call.hp1.model.txt
    	cmd="python ${I_SRC}/src/lstm-freq/5mc_freq.py --model_path ${I_MODEL_PATH} --kmer_size ${I_KMER_SIZE} --input ${features} --out ${freq_call} --dl_workers ${I_THREADS}"
    	echo "Running command"
    	echo "  ${cmd}"
    	python ${I_SRC}/src/lstm-freq/5mc_freq.py --model_path ${I_MODEL_PATH} --kmer_size ${I_KMER_SIZE} --input ${features} --out ${freq_call} --dl_workers ${I_THREADS}
    	if [ $? -ne 0 ]; then
        	echo "fail"
        	exit 1
    	fi 
	fi 

	if [ ${I_PHASE} -eq 1 ]; then
		features=${I_OUTPUT}/freq-call.hp2.model.features
    	freq_call=${I_OUTPUT}/freq-call.hp2.model.txt
    	cmd="python ${I_SRC}/src/lstm-freq/5mc_freq.py --model_path ${I_MODEL_PATH} --kmer_size ${I_KMER_SIZE} --input ${features} --out ${freq_call} --dl_workers ${I_THREADS}"
    	echo "Running command"
    	echo "  ${cmd}"
    	python ${I_SRC}/src/lstm-freq/5mc_freq.py --model_path ${I_MODEL_PATH} --kmer_size ${I_KMER_SIZE} --input ${features} --out ${freq_call} --dl_workers ${I_THREADS}
    	if [ $? -ne 0 ]; then
        	echo "fail"
        	exit 1
    	fi 
	fi 

    touch ${job_is_done}
}

extract_features

if [ "${I_MODEL_PATH}"x = x ]; then
	exit 0
fi

call_freq

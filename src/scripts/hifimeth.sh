#!/bin/bash

I_SRC=""
O_SRC="src"

I_OUTPUT=""
O_OUTPUT="out"

X_MOTIF="CG"
I_MOTIF=${X_MOTIF}
O_MOTIF="motif"

X_KMER_SIZE=400
I_KMER_SIZE=${X_KMER_SIZE}
O_KMER_SIZE="kmer_size"

X_THREADS=1
I_THREADS=${X_THREADS}
O_THREADS="num_threads"

X_MODEL_PATH=""
I_MODEL_PATH=""
O_MODEL_PATH="model_path"

X_INPUT=""
I_INPUT=""
O_INPUT="input"
I_INPUT_IS_BAM=1

X_MAPQ=10
I_MAPQ=${X_MAPQ}
O_MAPQ="mapping_quality"

X_REFERENCE=""
I_REFERENCE=""
O_REFERENCE="reference"

function parse_arguments()
{
	ARGS=`getopt --long ${O_SRC}:,${O_OUTPUT}:,${O_MOTIF}:,${O_KMER_SIZE}:,${O_THREADS}:,${O_MODEL_PATH}:,${O_INPUT}:,${O_MAPQ}:,${O_REFERENCE}: -- "$@"`
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
			--${O_INPUT})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_INPUT}], $2
				fi
				I_INPUT=$2
				shift 2
				;;
			--${O_MOTIF})
				if [ ${verbose} -eq 1 ]; then
					echo match [${O_MOTIF}], $2
				fi
				I_MOTIF=$2
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
	echo "${pn} [--${O_MOTIF}] [--${O_KMER_SIZE}] [--${O_THREADS}]" 
	echo "  <--${O_SRC}> <--${O_OUTPUT}> <--${O_INPUT}> <--${O_MODEL_PATH}>"

	echo ""
	echo "OPTION DESCRIPTIONS:"
	echo "--${O_MOTIF} <string>"
	echo "    The motif sequence"
	echo "    Default = '${X_MOTIF}'"
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
	echo "--${O_INPUT} <file>"
	echo "    Path to sample reads (SAM, BAM)"
	echo "--${O_MAPQ} <integer>"
	echo "    Minimum mapping quality"
	echo "    Default = '${X_MAPQ}'"
	echo "--${O_REFERENCE} <fa>"
	echo "    Reference path"
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

	if [ "${I_MODEL_PATH}"x = x ]; then
		echo "path to model is not provided"
		return 1
	fi

	if [ "${I_INPUT}"x = x ]; then
		echo "sample reads are not provided"
		return 1
	fi

	suffix="${I_INPUT##*.}"
	if [ "${suffix}"x = x ]; then
		echo unrecognised file format: ${I_INPUT}. Should be SAM or BAM
		return 1
	fi
	suffix=$(echo ${suffix} | awk '{print toupper($0)}')
    	if [ "${suffix}"x = "SAM"x ]; then
        	I_INPUT_IS_BAM=0
    	elif [ "${suffix}"x = "BAM"x ]; then
        	I_INPUT_IS_BAM=1
    	else
        	echo "unrecognised file format: ${I_INPUT}"
        	echo "  suffix should be sam or bam"
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
	echo "motif: ${I_MOTIF}"
	echo "kmer size: ${I_KMER_SIZE}"
	echo "CPU threads: ${I_THREADS}"
	echo "minimum mapping quality: ${I_MAPQ}"
	echo "reference: ${I_REFERENCE}"

	if [ ${I_INPUT_IS_BAM} -eq 1 ]; then
		echo "sample reads (BAM): ${I_INPUT}"
	else
		echo "sample reads (SAM): ${I_INPUT}"
	fi
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
    local job_is_done=${I_OUTPUT}/extract_features.done
    if [ -f ${job_is_done} ]; then
        return
    fi 

    if [ ${I_INPUT_IS_BAM} -eq 1 ]; then
	    cmd="samtools view -@8 ${I_INPUT} | ${I_SRC}/src/cpp_tools/extract_read_features -k ${I_KMER_SIZE} -t ${I_THREADS} -m ${I_MOTIF} -q ${I_MAPQ} - ${I_OUTPUT} ${I_REFERENCE}"
	    echo "Running command"
	    echo "  ${cmd}"
	    samtools view -@8 ${I_INPUT} | ${I_SRC}/src/cpp_tools/extract_read_features -k ${I_KMER_SIZE} -t ${I_THREADS} -m ${I_MOTIF} -q ${I_MAPQ} - ${I_OUTPUT} ${I_REFERENCE}
	    if [ $? -ne 0 ]; then
		    echo fail
		    exit 1
	    fi
    else
        cmd="${I_SRC}/src/cpp_tools/extract_read_features -k ${I_KMER_SIZE} -t ${I_THREADS} -m ${I_MOTIF} -q ${I_MAPQ} ${I_INPUT} ${I_OUTPUT} ${I_REFERENCE}"
	    run_commands ${cmd}
    fi

    touch ${job_is_done}
}

function call_5mc()
{
    local job_is_done=${I_OUTPUT}/call_5mc.done
    if [ -f ${job_is_done} ]; then
        return
    fi

    cmd="python ${I_SRC}/src/py-mol/x-gcn5mc.py --model ${I_MODEL_PATH} --kmer_size ${I_KMER_SIZE} --dl_workers ${I_THREADS} --features ${I_OUTPUT}"
    echo "running commnad"
    echo "  ${cmd}"
    python ${I_SRC}/src/py-mol/x-gcn5mc.py --model ${I_MODEL_PATH} --kmer_size ${I_KMER_SIZE} --dl_workers ${I_THREADS} --features ${I_OUTPUT}
    if [ $? -ne 0 ]; then 
        exit 1
    fi

	MOD_PATH=${I_OUTPUT}/5mc-call.txt
	rm -f ${MOD_PATH}
	NP=$(cat ${I_OUTPUT}/num_chunks)
	for ((i=0; i<${NP}; i++))
	do 
		echo "merge chunk ${i}"
		paste ${I_OUTPUT}/chunk_${i}.info ${I_OUTPUT}/chunk_${i}.prob >> ${MOD_PATH}
	done 

    touch ${job_is_done}
}

function add_5mc_tag_to_sam()
{
    local job_is_done=${I_OUTPUT}/add_5mc_tag.done
    if [ -f ${job_is_done} ]; then
        return 
    fi 

    MOD_PATH=${I_OUTPUT}/5mc-call.txt
    BAM_OUTPUT=${I_OUTPUT}/5mc-call.bam

    if [ ${I_INPUT_IS_BAM} -eq 1 ]; then
	    cmd="samtools view -h -@8 ${I_INPUT} | ${I_SRC}/src/cpp_tools/add_5mc_tags ${I_MOTIF} ${MOD_PATH} - - | samtools view -@8 -bS -o ${BAM_OUTPUT}"
	    echo "running command"
	    echo "  ${cmd}"
	    samtools view -h -@8 ${I_INPUT} | ${I_SRC}/src/cpp_tools/add_5mc_tags ${I_MOTIF} ${MOD_PATH} - - | samtools view -@8 -bS -o ${BAM_OUTPUT}
	    if [ $? -ne 0 ]; then
		    echo fail
		    exit 1
	    fi
    else
	    cmd="${I_SRC}/src/cpp_tools/add_5mc_tags ${I_MOTIF} ${MOD_PATH} ${I_INPUT} - | samtools view -bS -o ${BAM_OUTPUT}"
	    echo "running command"
	    echo "  ${cmd}"
	    ${I_SRC}/src/cpp_tools/add_5mc_tags ${I_MOTIF} ${MOD_PATH} ${I_INPUT} - | samtools view -bS -o ${BAM_OUTPUT}
	    if [ $? -ne 0 ]; then
		    echo fail
		    exit 1
	    fi
    fi

    touch ${job_is_done}
}

extract_features

call_5mc

add_5mc_tag_to_sam

#rm -f ${I_OUTPUT}/chunk_*
#rm -f ${I_OUTPUT}/*.done 
#rm -f ${I_OUTPUT}/num_chunks
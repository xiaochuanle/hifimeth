#!/bin/bash

src=$(pwd)/src/cpp_src
bin=$(pwd)/src/cpp_tools

cauxsrc="${src}/hbn_aux.c ${src}/line_reader.c"
cppauxsrc="${src}/5mc_aux.cpp ${src}/kmer_signal.cpp"

BUILDDEBUG=0
if [ ${BUILDDEBUG} -eq 1 ]; then
	CXXFLAGS="-D_GLIBCXX_PARALLEL -pthread -O0 -Wall -g -std=c++11 -mpopcnt"
	CFLAGS="-D_GLIBCXX_PARALLEL -pthread -O0 -Wall -g -std=gnu99 -mpopcnt"
else
	CXXFLAGS="-D_GLIBCXX_PARALLEL -pthread -O3 -Wall -std=c++11 -DNDEBUG -mpopcnt"
	CFLAGS="-D_GLIBCXX_PARALLEL -pthread -O3 -Wall -std=gnu99 -DNDEBUG -mpopcnt"
fi

LDFLAGS="-pthread -lm -lz -lstdc++"

OSTYPE=$(echo `uname`)
echo ${OSTYPE}
if [ "${OSTYPE}" == "Linux" ]; then
	CXXFLAGS="${CXXFLAGS} -fopenmp"
	CFLAGS="${CFLAGS} -fopenmp"
	LDFLAGS="${LDFLAGS} -fopenmp"
fi

echo ${CXXFLAGS}
echo ${CFLAGS}
echo ${LDFLAGS}

mkdir -p ${bin}

echo [$(date)] build hbn_aux.o
gcc -c ${src}/hbn_aux.c -o ${bin}/hbn_aux.o ${CFLAGS}

echo [$(date)] build line_reader.o
gcc -c ${src}/line_reader.c -o ${bin}/line_reader.o ${CFLAGS}

echo [$(date)] build 5mc_aux.o
g++ -c ${src}/5mc_aux.cpp -o ${bin}/5mc_aux.o ${CXXFLAGS}

echo [$(date)] build 5mc_tags.o
g++ -c ${src}/5mc_tags.cpp -o ${bin}/5mc_tags.o ${CXXFLAGS}

echo [$(date)] build kmer_signal.o
g++ -c ${src}/kmer_signal.cpp -o ${bin}/kmer_signal.o ${CXXFLAGS}

aux="${bin}/hbn_aux.o ${bin}/line_reader.o ${bin}/5mc_aux.o ${bin}/5mc_tags.o ${bin}/kmer_signal.o"

name="add_5mc_tags"
echo [$(date)] build ${name}
g++ -c ${src}/${name}.cpp -o ${bin}/${name}.o ${CXXFLAGS}
g++ ${bin}/${name}.o ${aux} -o ${bin}/${name} ${LDFLAGS}

name="extract_freq_features_from_sam"
echo [$(date)] build ${name}
g++ -c ${src}/${name}.cpp -o ${bin}/${name}.o ${CXXFLAGS}
g++ ${bin}/${name}.o ${aux} -o ${bin}/${name} ${LDFLAGS}

name="extract_freq_features"
echo [$(date)] build ${name}
g++ -c ${src}/${name}.cpp -o ${bin}/${name}.o ${CXXFLAGS}
g++ ${bin}/${name}.o ${aux} -o ${bin}/${name} ${LDFLAGS}

name="extract_train_features"
echo [$(date)] build ${name}
g++ -c ${src}/${name}.cpp -o ${bin}/${name}.o ${CXXFLAGS}
g++ ${bin}/${name}.o ${aux} -o ${bin}/${name} ${LDFLAGS}

name="extract_read_features"
echo [$(date)] build ${name}
g++ -c ${src}/${name}.cpp -o ${bin}/${name}.o ${CXXFLAGS}
g++ ${bin}/${name}.o ${aux} -o ${bin}/${name} ${LDFLAGS}

name="merge_features"
echo [$(date)] build ${name}
g++ -c ${src}/${name}.cpp -o ${bin}/${name}.o ${CXXFLAGS}
g++ ${bin}/${name}.o ${aux} -o ${bin}/${name} ${LDFLAGS}

rm -f ${bin}/*.o

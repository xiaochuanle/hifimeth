# How to build HiFiMeth

`HiFiMeth` is compiled on a server running CentOS Linux 7 with GCC/G++ version 8.5.0. If your system compiler is outdated, please install HiFiMeth via [Conda](https://www.anaconda.com/docs/getting-started/miniconda/main) in a new environment:
```shell
$ conda create -n GCC8 -c conda-forge gcc_linux-64=8.5.0 gxx_linux-64=8.5.0 -y
$ conda activate GCC8
$ ln -sf $CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gcc $CONDA_PREFIX/bin/gcc
$ ln -sf $CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++ $CONDA_PREFIX/bin/g++
$ ln -sf $CONDA_PREFIX/bin/gcc $CONDA_PREFIX/bin/cc
$ ln -sf $CONDA_PREFIX/bin/g++ $CONDA_PREFIX/bin/c++
```
Installing GCC and G++ in an isolated Conda environment ensures that your existing system build environment remains untouched. You can verify the successful installation by running the following commands:
```shell
$ gcc --version
  gcc (conda-forge gcc 8.5.0-17) 8.5.0
$ g++ --version
  g++ (conda-forge gcc 8.5.0-17) 8.5.0
```

***Remember to activate your environment (`conda activate GCC8`) before the installation.***

1. Clone HiFiMeth repository and goto the HiFiMeth source code directory:
```shell
$ git clone https://github.com/xiaochuanle/hifimeth.git
$ cd hifimeth
$ export HiFiMethHome=$(pwd)
```
2. Download third-party dependencies:
```shell
$ mkdir -p 3rdparty
$ cd 3rdparty
$ wget --no-check-certificate https://zlib.net/fossils/zlib-1.3.1.tar.gz
$ wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2
```
3. Install [OpenVINO](https://github.com/openvinotoolkit/openvino/tree/2025.4.0). Please be aware that compiling OpenVINO from source can be an arduous process with a high risk of failure.
```shell
$ git clone -b 2025.4.0 https://github.com/openvinotoolkit/openvino.git
$ git submodule update --init --recursive
$ cd openvino/
$ mkdir build && cd build
$ cmake -DENABLE_INTEL_GPU=OFF \
      -DENABLE_INTEL_NPU=OFF \
      -DENABLE_TEMPLATE=OFF \
      -DENABLE_HETERO=OFF \
      -DENABLE_MULTI=OFF \
      -DENABLE_AUTO=OFF \
      -DENABLE_AUTO_BATCH=OFF \
      -DENABLE_OV_ONNX_FRONTEND=ON \
      -DENABLE_OV_PADDLE_FRONTEND=OFF \
      -DENABLE_OV_TF_FRONTEND=OFF \
      -DENABLE_OV_TF_LITE_FRONTEND=OFF \
      -DENABLE_OV_JAX_FRONTEND=OFF \
      -DENABLE_OV_PYTORCH_FRONTEND=OFF \
      -DENABLE_OV_JAX_FRONTEND=OFF \
      -DENABLE_INTEL_CPU=ON \
      -DENABLE_OV_IR_FRONTEND=ON \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DENABLE_PYTHON=OFF \
      -DENABLE_WHEEL=OFF \
      -DENABLE_CLANG_FORMAT=OFF \
      -DENABLE_SAMPLES=OFF \
      -DENABLE_TEMPLATE=OFF \
      -DENABLE_SSE42=ON \
      -DENABLE_AVX2=ON \
      ..
$ cmake --build . --parallel 12
$ cmake -DCMAKE_INSTALL_PREFIX=. -P cmake_install.cmake

$ cd $HiFiMethHome
$ OSTYPE=$(uname -s)
$ MACHINETYPE=$(arch=$(uname -m); echo ${arch/x86_64/amd64})
$ TARGET_DIR=src/app/${OSTYPE}-${MACHINETYPE}/bin
$ mkdir -p ${TARGET_DIR}/openvino
$ mkdir -p ${TARGET_DIR}/tbb
$ cp -P 3rdparty/openvino/build/runtime/3rdparty/tbb/lib/lib* ${TARGET_DIR}/tbb
$ cp -P 3rdparty/openvino/build/runtime/lib/intel64/libopenvino* ${TARGET_DIR}/openvino
$ cp -r models ${TARGET_DIR}
```

4. Compile `HiFiMeth`.
```shell
$ cd src/app
$ make
$ cd $HiFiMethHome
```

The binary executable `hifimeth` is presented at `${HiFiMethHome}/${TARGET_DIR}/hifimeth`.
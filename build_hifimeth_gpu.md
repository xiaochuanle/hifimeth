# Howto build HiFiMeth-GPU

`HiFiMeth-GPU` is compiled and tested on a computer running Ubuntu 22.04.4 LTS with GCC and G++ version 11.4.0.

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
$ wget https://download.pytorch.org/libtorch/cu118/libtorch-cxx11-abi-shared-with-deps-2.7.1%2Bcu118.zip
$ unzip libtorch-cxx11-abi-shared-with-deps-2.7.1%2Bcu118.zip
$ cp libtorch/lib/libnvrtc-672ee683.so.11.2 libtorch/lib/libnvrtc.so.11.2
$ cp libtorch/lib/libnvrtc-builtins-2dc4bf68.so.11.8 libtorch/lib/libnvrtc-builtins.so.11.8
```

3. Build `HiFiMeth-GPU`
```shell
$ cd $HiFiMethHome
$ OSTYPE=$(uname -s)
$ MACHINETYPE=$(arch=$(uname -m); echo ${arch/x86_64/amd64})
$ TARGET_DIR=src/app-gpu/${OSTYPE}-${MACHINETYPE}/bin
$ mkdir -p $TARGET_DIR
$ ln -s $(pwd)/3rdparty/libtorch/lib/ $(pwd)/${TARGET_DIR}/libtorch
$ cp -r models ${TARGET_DIR}
$ cd src/app-gpu
$ make
$ cd $HiFiMethHome
```
The binary executable `hifimeth` is presented at `${HiFiMethHome}/${TARGET_DIR}/hifimeth-gpu`.
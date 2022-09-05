This repository accompanies the paper **Parallel Small Polynomial Multiplication for Dilithium: A Faster Design and Implementation** accepted by ACSAC22.


Authors: 
 - Jieyu  Zheng`<jieyuzheng21@m.fudan.edu.cn>`
 - Feng He `<fhe20@fudan.edu.cn>`
 - Shiyu Shen`<shenshiyu21@m.fudan.edu.cn>`
 - Chenxi Xue`<cxxue21@m.fudan.edu.cn>` 
 - Yunlei Zhao `<ylzhao@fudan.edu.cn>`

It contains our source code for Dilithium with parallel small polynomial multiplication algorithm optimized for both C reference implementation and Arm Cortex-A72 neon implementation.

# Installation
## Cloning the code
Clone the code 

```
git clone https://github.com/zhengjieyu/Dilithium-smallpoly.git
```

## C reference implementation code

We have tested this code on an Intel(R) Core(TM) i7-10510U CPU at 2.3GHz(16 GB memory) with TurBo Boost and Hyperthreading disabled. The operating system is Ubuntu 20.04 LTS with Linux Kernel 4.4.0 and the gcc version is 9.4.0.

Before running the code, you should install openssl. For linux, you can install openssl by this instruction


```
sudo apt-get install openssl
sudo apt-get install libssl-dev
```

To compile C reference code, go to the C_reference_implementation directory,  for example, for Dilithium-3, go to the subdirectory dilithium3-smallpoly, run

```
make
```
This produces the executables 

```
test/test_dilithium3
test/test_speed3
```
The first executable is used for verifying correctness of code, the second is for benchmark.
## Arm neon implementation code

We have tested this code with a Raspberry Pi 4B (RPi  4 Model B) with ARMv8-A instruction sets, Cortex-A72(1.8 GHz) CPU and 4GB RAM, it also supports ARMv8-A Neon SIMD instructions, running Ubuntu 10.3.0.
It is essential that you are running a 64-bit OS to be able to execute aarch64 code.
The 64-bit Raspbian OS can be used, but we have not tested all code with it.

The Makefiles included in this repo assume that you are natively compiling your
code using gcc. We have tested it with gcc 10.3.0.

### Enable access to the cycle counter

For accurate benchmarking on Cortex-A72, we make use of the performance counters.
By default the access from user mode is disabled. You will need to enable it using a kernel module.

We have included on here that you can install
```
cd enable_ccr
make install
```

Alternatively, you can get one from: https://github.com/rdolbeau/enable_arm_pmu.

You may have to install the kernel headers manually:
https://www.raspberrypi.com/documentation/computers/linux_kernel.html

```
sudo apt install raspberrypi-kernel-headers
```

## Benchmarking on Cortex-A72
To compile the test programs, go to the Arm_Neon_implementation directory, run
```
make speed 
make all
```
This produces the executables 

```
test/test_speed2_test2
test/test_speed3_test2
test/test_speed5_test2
test/test_dilithium2
test/test_dilithium3
test/test_dilithium5
```
The first three executables are used for benchmark, the following three executables are used for verifying correctness of code.
# License

This repository includes code from other sources that has the following license/license waivers
- `feat.S` modified from https://github.com/bwesterb/armed-keccak: MIT
- Dilithium reference code https://github.com/pq-crystals/dilithium/blob/master/LICENSE: CC0
- `fips202.{c,h}` http://bench.cr.yp.to/supercop.html: public domain
- `fips202x2.{c,h}` https://github.com/cothan/kyber/blob/master/neon/fips202x2.c: CC0
- `m1cycles.{c, h}`: https://github.com/cothan/SABER/blob/master/Cortex-A_Implementation_KEM/m1cycles.c: public domain
- `gen_table/common` from https://github.com/multi-moduli-ntt-saber/multi-moduli-ntt-saber: CC0

All remaining code is covered by CC0.

# Build Instructions

The source runs on Linux and requires an x86 CPU supporting the AVX2 instruction set. To compile, just run `make`. This builds the binaries `test_addtion` and `test_mult` which test the integer addition and multiplication proof systems by proving relations between random integers. Also the binaries report the median and average counts of the time step counter (TSC) for the proving and verifying tasks.

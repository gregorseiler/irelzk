CC ?= gcc
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -mavx2 -mpopcnt -maes \
  -march=native -mtune=native -O3
#CFLAGS += -DUSE_RDPMC

default: test_ntt test_addition test_mult test_mult2 test_main

%_addition.o: %_addition.c params.h Makefile
	$(CC) $(CFLAGS) -DADDITION_PROOF -c $< -o $@

%_addition.o: %.c params.h Makefile
	$(CC) $(CFLAGS) -DADDITION_PROOF -c $< -o $@

%_mult2.o: %_mult2.c params.h Makefile
	$(CC) $(CFLAGS) -DMULTIPLICATION_PROOF_2 -c $< -o $@

%_mult2.o: %.c params.h Makefile
	$(CC) $(CFLAGS) -DMULTIPLICATION_PROOF_2 -c $< -o $@

%_mult.o: %_mult.c params.h Makefile
	$(CC) $(CFLAGS) -DMULTIPLICATION_PROOF -c $< -o $@

%_mult.o: %.c params.h Makefile
	$(CC) $(CFLAGS) -DMULTIPLICATION_PROOF -c $< -o $@

%_main.o: %_main.c params.h Makefile
	$(CC) $(CFLAGS) -DMAIN_PROOF -c $< -o $@

%_main.o: %.c params.h Makefile
	$(CC) $(CFLAGS) -DMAIN_PROOF -c $< -o $@

%.o: %.c params.h Makefile
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.S params.h fq.inc shuffle.inc Makefile
	$(CC) $(CFLAGS) -c $< -o $@

polyveck.o: polyvec.c polyvec.h params.h Makefile
	$(CC) $(CFLAGS) -DADDITION_PROOF -DPOLYVEC_TYPE=polyveck -DPOLYVEC_LENGTH=K -c $< -o $@

polyvecl.o: polyvec.c polyvec.h params.h Makefile
	$(CC) $(CFLAGS) -DADDITION_PROOF -DPOLYVEC_TYPE=polyvecl -DPOLYVEC_LENGTH=L -c $< -o $@

polyvecm_addition.o: polyvec.c polyvec.h params.h Makefile
	$(CC) $(CFLAGS) -DADDITION_PROOF -DPOLYVEC_TYPE=polyvecm -DPOLYVEC_LENGTH=M -c $< -o $@

polyvecm_mult.o: polyvec.c polyvec.h params.h Makefile
	$(CC) $(CFLAGS) -DMULTIPLICATION_PROOF -DPOLYVEC_TYPE=polyvecm -DPOLYVEC_LENGTH=M -c $< -o $@

polyvecm_mult2.o: polyvec.c polyvec.h params.h Makefile
	$(CC) $(CFLAGS) -DMULTIPLICATION_PROOF_2 -DPOLYVEC_TYPE=polyvecm -DPOLYVEC_LENGTH=M -c $< -o $@

polyvecm_main.o: polyvec.c polyvec.h params.h Makefile
	$(CC) $(CFLAGS) -DMAIN_PROOF -DPOLYVEC_TYPE=polyvecm -DPOLYVEC_LENGTH=M -c $< -o $@

test_ntt: test_ntt.o ntt.o invntt.o consts.o rounding.o poly.o aes256ctr.o randombytes.o cpucycles.o
	$(CC) $(CFLAGS) $^ -o $@

test_poly: test_poly.o ntt.o invntt.o consts.o rounding.o poly.o aes256ctr.o randombytes.o
	$(CC) $(CFLAGS) $^ -o $@

test_addition: test_addition.o ntt.o invntt.o consts.o rounding.o poly.o polyveck.o polyvecl.o \
  polyvecm_addition.o comm_addition.o opening_addition.o \
  product_addition.o linear_addition.o addition_addition.o \
  aes256ctr.o fips202.o randombytes.o cpucycles.o speed_print.o
	$(CC) $(CFLAGS) -DADDITION_PROOF $^ -o $@

test_mult: test_mult.o ntt.o invntt.o consts.o rounding.o poly.o polyveck.o polyvecl.o \
  polyvecm_mult.o comm_mult.o opening_mult.o \
  product_mult.o linear_mult.o apprshort_mult.o mult_mult.o \
  aes256ctr.o fips202.o randombytes.o cpucycles.o speed_print.o
	$(CC) $(CFLAGS) -DMULTIPLICATION_PROOF $^ -o $@

test_mult2: test_mult2.o ntt.o invntt.o consts.o rounding.o poly.o polyveck.o polyvecl.o \
  polyvecm_mult2.o comm_mult2.o opening_mult2.o \
  product_mult2.o linear_mult2.o apprshort_mult2.o mult2_mult2.o \
  aes256ctr.o fips202.o randombytes.o cpucycles.o speed_print.o
	$(CC) $(CFLAGS) -DMULTIPLICATION_PROOF_2 $^ -o $@

test_main: test_main.o ntt.o invntt.o consts.o rounding.o poly.o polyveck.o polyvecl.o \
  polyvecm_main.o comm_main.o opening_main.o product_main.o linear_main.o main_main.o \
  aes256ctr.o fips202.o randombytes.o cpucycles.o speed_print.o
	$(CC) $(CFLAGS) -DMAIN_PROOF $^ -o $@

clean:
	rm -f *.o
	rm -f test_ntt
	rm -f test_poly
	rm -f test_addition
	rm -f test_mult
	rm -f test_mult2
	rm -f test_main

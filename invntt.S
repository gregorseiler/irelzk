#include "consts.h"
.include "shuffle.inc"

.macro butterfly l0,l1,h0,h1,zl0=15,zh0=2,zl1=15,zh1=2,flag=0
vpsubd		%ymm\l0,%ymm\h0,%ymm8
vpsubd		%ymm\l1,%ymm\h1,%ymm9
vpaddd		%ymm\h0,%ymm\l0,%ymm\l0

vpmuldq		%ymm\zl0,%ymm8,%ymm12
vpmuldq		%ymm\zl1,%ymm9,%ymm13
vmovshdup	%ymm8,%ymm10

vmovshdup	%ymm9,%ymm11
.if \flag
vpsrlq		$32,%ymm\zl0,%ymm\zl0
.if \zl0 != \zl1
vpsrlq		$32,%ymm\zl1,%ymm\zl1
.endif
.endif

vpmuldq		%ymm\zl0,%ymm10,%ymm14
vpmuldq		%ymm\zl1,%ymm11,%ymm15
.if \flag && \zh0 != \zh1
vmovshdup	%ymm\zh0,%ymm\h0
.endif

vpmuldq		%ymm\zh0,%ymm8,%ymm8
vpmuldq		%ymm\zh1,%ymm9,%ymm9
.if \flag
vmovshdup	%ymm\zh1,%ymm\zh1
.endif

.if \flag && \zh0 != \zh1
vpmuldq		%ymm\h0,%ymm10,%ymm10
.else
vpmuldq		%ymm\zh0,%ymm10,%ymm10
.endif
vpmuldq		%ymm\zh1,%ymm11,%ymm11
vpaddd		%ymm\h1,%ymm\l1,%ymm\l1

vpmuldq		%ymm0,%ymm12,%ymm12
vpmuldq		%ymm0,%ymm13,%ymm13
vpand		%ymm1,%ymm\l0,%ymm\h0

vpsrad		$30,%ymm\l0,%ymm\l0
vpsrad		$30,%ymm\l1,%ymm\h1
vpand		%ymm1,%ymm\l1,%ymm\l1

vpmuldq		%ymm0,%ymm14,%ymm14
vpmuldq		%ymm0,%ymm15,%ymm15
vpsubd		%ymm\l0,%ymm\h0,%ymm\h0

vpsubd		%ymm\h1,%ymm\l1,%ymm\l1
vpsrlq		$32,%ymm8,%ymm8
vmovshdup	%ymm9,%ymm9

vpblendd	$0xAA,%ymm10,%ymm8,%ymm8
vpblendd	$0xAA,%ymm11,%ymm9,%ymm9
vpslld		$18,%ymm\l0,%ymm\l0

vpsrlq		$32,%ymm12,%ymm12
vmovshdup	%ymm13,%ymm13
vpaddd		%ymm\l0,%ymm\h0,%ymm\l0

vpslld		$18,%ymm\h1,%ymm\h1
vpblendd	$0xAA,%ymm14,%ymm12,%ymm12
vpblendd	$0xAA,%ymm15,%ymm13,%ymm13

vpaddd		%ymm\h1,%ymm\l1,%ymm\l1
vpsubd		%ymm12,%ymm8,%ymm\h0
vpsubd		%ymm13,%ymm9,%ymm\h1
.endm

.macro levels0t4 off
/* level0 */
vmovdqa		_PMASK*4(%rsi),%ymm3
vpermd		92(%rdx),%ymm3,%ymm14
vpermd		60(%rdx),%ymm3,%ymm15
vpermd		(_ZETAS-_ZETAS_QINV+23)*4(%rdx),%ymm3,%ymm2
vpermd		(_ZETAS-_ZETAS_QINV+15)*4(%rdx),%ymm3,%ymm3

/* load */
vmovdqa		(\off+ 0)*4(%rdi),%ymm4
vmovdqa		(\off+ 8)*4(%rdi),%ymm5
vmovdqa		(\off+16)*4(%rdi),%ymm6
vmovdqa		(\off+24)*4(%rdi),%ymm7

butterfly	4,6,5,7,14,2,15,3,1

/* level1 */
vmovdqa		_PMASK*4(%rsi),%ymm3
vpermd		28(%rdx),%ymm3,%ymm15
vpermd		(_ZETAS-_ZETAS_QINV+7)*4(%rdx),%ymm3,%ymm2

butterfly	4,5,6,7,15,2,15,2,1

/* level2 */
vpmovzxdq	12(%rdx),%ymm15
vpmovzxdq	(_ZETAS-_ZETAS_QINV+3)*4(%rdx),%ymm2
vpermq		$0x1B,%ymm15,%ymm15
vpermq		$0x1B,%ymm2,%ymm2

shuffle2	4,5,3,5
shuffle2	6,7,4,7

butterfly	3,4,5,7

/* level3 */
vpbroadcastd	8(%rdx),%ymm14
vpbroadcastd	4(%rdx),%ymm15
vpbroadcastd	(_ZETAS-_ZETAS_QINV+2)*4(%rdx),%ymm2
vpbroadcastd	(_ZETAS-_ZETAS_QINV+1)*4(%rdx),%ymm6
vpblendd	$0xF0,%ymm15,%ymm14,%ymm15
vpblendd	$0xF0,%ymm6,%ymm2,%ymm2

shuffle4	3,4,6,4
shuffle4	5,7,3,7

butterfly	6,3,4,7

/* level4 */
vpbroadcastd	(%rdx),%ymm15
vpbroadcastd	(_ZETAS-_ZETAS_QINV)*4(%rdx),%ymm2

shuffle8	6,3,5,3
shuffle8	4,7,6,7

butterfly	5,6,3,7

/* store */
vmovdqa		%ymm5,(\off+ 0)*4(%rdi)
vmovdqa		%ymm6,(\off+ 8)*4(%rdi)
vmovdqa		%ymm3,(\off+16)*4(%rdi)
vmovdqa		%ymm7,(\off+24)*4(%rdi)
.endm

.macro levels5t6 off
/* level5 */
vpbroadcastd	(_ZETAS_QINV+2)*4(%rsi),%ymm14
vpbroadcastd	(_ZETAS_QINV+1)*4(%rsi),%ymm15
vpbroadcastd	(_ZETAS+2)*4(%rsi),%ymm2
vpbroadcastd	(_ZETAS+1)*4(%rsi),%ymm3

/* load */
vmovdqa		(\off+ 0)*4(%rdi),%ymm4
vmovdqa		(\off+32)*4(%rdi),%ymm5
vmovdqa		(\off+64)*4(%rdi),%ymm6
vmovdqa		(\off+96)*4(%rdi),%ymm7

butterfly	4,6,5,7,14,2,15,3

/* level6 */
vpbroadcastd	(_ZETAS_QINV)*4(%rsi),%ymm15
vpbroadcastd	(_ZETAS)*4(%rsi),%ymm2

butterfly	4,5,6,7

/* store */
vmovdqa         %ymm4,(\off+ 0)*4(%rdi)
vmovdqa         %ymm5,(\off+32)*4(%rdi)
vmovdqa         %ymm6,(\off+64)*4(%rdi)
vmovdqa         %ymm7,(\off+96)*4(%rdi)
.endm

.text
.global cdecl(invntt_tomont_avx)
cdecl(invntt_tomont_avx):
#consts
vmovdqa		_8XQ*4(%rsi),%ymm0
vpcmpeqd	%ymm1,%ymm1,%ymm1
vpsrld		$2,%ymm1,%ymm1

/* levels0t4 */
lea		(_ZETAS_QINV+3+3*31)*4(%rsi),%rdx
levels0t4	0
sub		$124,%rdx
levels0t4	32
sub		$124,%rdx
levels0t4	64
sub		$124,%rdx
levels0t4	96

/* levels2t6 */
levels5t6	0
levels5t6	8
levels5t6	16
levels5t6	24

ret

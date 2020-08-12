#include "consts.h"
.include "shuffle.inc"

.macro butterfly l0,l1,h0,h1,zl0=2,zh0=3,zl1=2,zh1=3,flag=0
vpmuldq		%ymm\zl0,%ymm\h0,%ymm10
vpmuldq		%ymm\zl1,%ymm\h1,%ymm11
vmovshdup	%ymm\h0,%ymm8

vmovshdup	%ymm\h1,%ymm9
.if \flag
vpsrlq		$32,%ymm\zl0,%ymm\zl0
.if \zl0 != \zl1
vpsrlq		$32,%ymm\zl1,%ymm\zl1
.endif
.endif

vpmuldq		%ymm\zl0,%ymm8,%ymm12
vpmuldq		%ymm\zl1,%ymm9,%ymm13
.if \flag
vmovshdup	%ymm\zh0,%ymm\zl0
.endif

vpmuldq		%ymm\zh0,%ymm\h0,%ymm\h0
vpmuldq		%ymm\zh1,%ymm\h1,%ymm\h1
.if \flag && \zh0 != \zh1
vmovshdup	%ymm\zh1,%ymm\zl1
.endif

.if \flag
vpmuldq		%ymm\zl0,%ymm8,%ymm8
vpmuldq		%ymm\zl1,%ymm9,%ymm9
.else
vpmuldq		%ymm\zh0,%ymm8,%ymm8
vpmuldq		%ymm\zh1,%ymm9,%ymm9
.endif
vpand		%ymm1,%ymm\l0,%ymm\zl0

vpmuldq		%ymm0,%ymm10,%ymm10
vpmuldq		%ymm0,%ymm11,%ymm11
vpand		%ymm1,%ymm\l1,%ymm\zh0

vpsrad		$30,%ymm\l0,%ymm\l0
vpsrad		$30,%ymm\l1,%ymm\l1

vpmuldq		%ymm0,%ymm12,%ymm12
vpmuldq		%ymm0,%ymm13,%ymm13
vpsubd		%ymm\l0,%ymm\zl0,%ymm\zl0

vpsrlq		$32,%ymm\h0,%ymm\h0
vmovshdup	%ymm\h1,%ymm\h1
vpslld		$18,%ymm\l0,%ymm\l0

vpsubd		%ymm\l1,%ymm\zh0,%ymm\zh0
vpblendd	$0xAA,%ymm8,%ymm\h0,%ymm\h0
vpblendd	$0xAA,%ymm9,%ymm\h1,%ymm\h1

vpsrlq		$32,%ymm10,%ymm10
vmovshdup	%ymm11,%ymm11
vpslld		$18,%ymm\l1,%ymm\l1

vpaddd		%ymm\l0,%ymm\zl0,%ymm\l0
vpblendd	$0xAA,%ymm12,%ymm10,%ymm10
vpblendd	$0xAA,%ymm13,%ymm11,%ymm11

vpsubd		%ymm10,%ymm\h0,%ymm8
vpsubd		%ymm11,%ymm\h1,%ymm9
vpaddd		%ymm\l1,%ymm\zh0,%ymm\l1

vpsubd		%ymm8,%ymm\l0,%ymm\h0
vpaddd		%ymm8,%ymm\l0,%ymm\l0
vpsubd		%ymm9,%ymm\l1,%ymm\h1

vpaddd		%ymm9,%ymm\l1,%ymm\l1
.endm

.macro levels0t1 off
/* level0 */
vpbroadcastd	_ZETAS_QINV*4(%rsi),%ymm2
vpbroadcastd	_ZETAS*4(%rsi),%ymm3

/* load */
vmovdqa		(8*\off+64)*4(%rdi),%ymm6
vmovdqa		(8*\off+96)*4(%rdi),%ymm7
vmovdqa		(8*\off+ 0)*4(%rdi),%ymm4
vmovdqa		(8*\off+32)*4(%rdi),%ymm5

butterfly	4,5,6,7

/* level1 */
vpbroadcastd	(_ZETAS_QINV+1)*4(%rsi),%ymm2
vpbroadcastd	(_ZETAS_QINV+2)*4(%rsi),%ymm14
vpbroadcastd	(_ZETAS+1)*4(%rsi),%ymm3
vpbroadcastd	(_ZETAS+2)*4(%rsi),%ymm15

butterfly	4,6,5,7,2,3,14,15

/* store */
vmovdqa		%ymm5,(8*\off+32)*4(%rdi)
vmovdqa		%ymm4,(8*\off+ 0)*4(%rdi)
vmovdqa		%ymm7,(8*\off+96)*4(%rdi)
vmovdqa		%ymm6,(8*\off+64)*4(%rdi)
.endm

.macro levels2t6 off
/* level2 */
vpbroadcastd	(_ZETAS_QINV+3+31*\off)*4(%rsi),%ymm2
vpbroadcastd	(_ZETAS+3+31*\off)*4(%rsi),%ymm3

/* load */
vmovdqa		(32*\off+16)*4(%rdi),%ymm6
vmovdqa		(32*\off+24)*4(%rdi),%ymm7
vmovdqa		(32*\off+ 0)*4(%rdi),%ymm4
vmovdqa		(32*\off+ 8)*4(%rdi),%ymm5

butterfly	4,5,6,7

/* level3 */
shuffle8	4,6,3,6
shuffle8	5,7,4,7

vpbroadcastd	(_ZETAS_QINV+3+31*\off+1)*4(%rsi),%ymm2
vpbroadcastd	(_ZETAS_QINV+3+31*\off+2)*4(%rsi),%ymm5
vpbroadcastd	(_ZETAS+3+31*\off+1)*4(%rsi),%ymm14
vpbroadcastd	(_ZETAS+3+31*\off+2)*4(%rsi),%ymm15
vpblendd	$0xF0,%ymm5,%ymm2,%ymm2
vpblendd	$0xF0,%ymm15,%ymm14,%ymm5

butterfly	3,6,4,7,2,5,2,5

/* level4 */
shuffle4	3,4,5,4
shuffle4	6,7,3,7

vpmovzxdq	(_ZETAS_QINV+3+31*\off+3)*4(%rsi),%ymm2
vpmovzxdq	(_ZETAS+3+31*\off+3)*4(%rsi),%ymm6

butterfly	5,4,3,7,2,6,2,6

/* level5 */
shuffle2	5,3,6,3
shuffle2	4,7,5,7

vmovdqu		(_ZETAS_QINV+3+31*\off+7)*4(%rsi),%ymm2
vmovdqu		(_ZETAS+3+31*\off+7)*4(%rsi),%ymm4

butterfly	6,3,5,7,2,4,2,4,1

/* level6 */
vmovdqu		(_ZETAS_QINV+3+31*\off+15)*4(%rsi),%ymm2
vmovdqu		(_ZETAS_QINV+3+31*\off+23)*4(%rsi),%ymm14
vmovdqu		(_ZETAS+3+31*\off+15)*4(%rsi),%ymm4
vmovdqu		(_ZETAS+3+31*\off+23)*4(%rsi),%ymm15

butterfly	6,5,3,7,2,4,14,15,1

/* store */
vmovdqa		%ymm3,(32*\off+ 8)*4(%rdi)
vmovdqa		%ymm6,(32*\off+ 0)*4(%rdi)
vmovdqa		%ymm7,(32*\off+24)*4(%rdi)
vmovdqa		%ymm5,(32*\off+16)*4(%rdi)
.endm

.global cdecl(ntt_avx)
cdecl(ntt_avx):
vmovdqa		_8XQ*4(%rsi),%ymm0
vpcmpeqd	%ymm1,%ymm1,%ymm1
vpsrld		$2,%ymm1,%ymm1

levels0t1	0
levels0t1	1
levels0t1	2
levels0t1	3

levels2t6	0
levels2t6	1
levels2t6	2
levels2t6	3

ret

.macro red32 r,x=12
#vpcmpeqd	%ymm2,%ymm2,%ymm2
#vpslld		$30,%ymm2,%ymm2
vpsrad		$30,%ymm\r,%ymm\x
vpand		%ymm2,%ymm\r,%ymm\r
vpsubd		%ymm\x,%ymm\r,%ymm\r
vpslld		$18,%ymm\x,%ymm\x
vpaddd		%ymm\x,%ymm\r,%ymm\r
.endm

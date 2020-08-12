#ifndef CONSTS_H
#define CONSTS_H

#include "params.h"

#define QINV -1073479679 // q^-1 mod 2^32
#define MONT 1048572 // 2^32 mod q
#define MONTSQ 260045840 // 2^64 mod q
#define DIV -132153352 // mont^2/128 mod q

#define _8XQ          0
#define _8XQINV       8
#define _8XDIV       16
#define _PMASK       24
#define _ZETAS       32
#define _ZETAS_QINV 160

/* The C ABI on MacOS exports all symbols with a leading
 * underscore. This means that any symbols we refer to from
 * C files (functions) can't be found, and all symbols we
 * refer to from ASM also can't be found.
 *
 * This define helps us get around this
 */
#if defined(__WIN32__) || defined(__APPLE__)
#define decorate(s) _##s
#define _cdecl(s) decorate(s)
#define cdecl(s) _cdecl(NAMESPACE(_##s))
#else
#define cdecl(s) NAMESPACE(_##s)
#endif

#ifndef __ASSEMBLER__
#include <stdint.h>
#include "poly.h"
#define qdata NAMESPACE(_qdata)
extern const int32_t qdata[];
#define rejidx NAMESPACE(_rejidx)
extern const uint8_t rejidx[256][8];
#define nttx NAMESPACE(_nttx)
extern const poly nttx;
#define nttx2 NAMESPACE(_nttx2)
extern const poly nttx2;
#define nttx3 NAMESPACE(_nttx3)
extern const poly nttx3;
#define nttx64 NAMESPACE(_nttx64)
extern const poly nttx64;
#define twist NAMESPACE(_twist)
extern const poly twist;
#define invtwist NAMESPACE(_invtwist)
extern const poly invtwist;
#endif

#endif

/* Autogenerated */
/* curve description: p251 */
/* requested operations: (all) */
/* n = 9 (from "9") */
/* s-c = 2^251 - [(1, 9)] (from "2^251 - 9") */
/* machine_wordsize = 32 (from "32") */

/* Computed values: */
/* carry_chain = [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1] */

#include <stdint.h>
typedef unsigned char fiat_p251_uint1;
typedef signed char fiat_p251_int1;
typedef signed __int128 fiat_p251_int128;
typedef unsigned __int128 fiat_p251_uint128;

#if (-1 & 3) != 3
#error "This code only works on a two's complement system"
#endif


/*
 * The function fiat_p251_addcarryx_u28 is an addition with carry.
 * Postconditions:
 *   out1 = (arg1 + arg2 + arg3) mod 2^28
 *   out2 = ⌊(arg1 + arg2 + arg3) / 2^28⌋
 *
 * Input Bounds:
 *   arg1: [0x0 ~> 0x1]
 *   arg2: [0x0 ~> 0xfffffff]
 *   arg3: [0x0 ~> 0xfffffff]
 * Output Bounds:
 *   out1: [0x0 ~> 0xfffffff]
 *   out2: [0x0 ~> 0x1]
 */
static void fiat_p251_addcarryx_u28(uint32_t* out1, fiat_p251_uint1* out2, fiat_p251_uint1 arg1, uint32_t arg2, uint32_t arg3) {
  uint32_t x1 = ((arg1 + arg2) + arg3);
  uint32_t x2 = (x1 & UINT32_C(0xfffffff));
  fiat_p251_uint1 x3 = (fiat_p251_uint1)(x1 >> 28);
  *out1 = x2;
  *out2 = x3;
}

/*
 * The function fiat_p251_subborrowx_u28 is a subtraction with borrow.
 * Postconditions:
 *   out1 = (-arg1 + arg2 + -arg3) mod 2^28
 *   out2 = -⌊(-arg1 + arg2 + -arg3) / 2^28⌋
 *
 * Input Bounds:
 *   arg1: [0x0 ~> 0x1]
 *   arg2: [0x0 ~> 0xfffffff]
 *   arg3: [0x0 ~> 0xfffffff]
 * Output Bounds:
 *   out1: [0x0 ~> 0xfffffff]
 *   out2: [0x0 ~> 0x1]
 */
static void fiat_p251_subborrowx_u28(uint32_t* out1, fiat_p251_uint1* out2, fiat_p251_uint1 arg1, uint32_t arg2, uint32_t arg3) {
  int32_t x1 = ((int32_t)(arg2 - arg1) - (int32_t)arg3);
  fiat_p251_int1 x2 = (fiat_p251_int1)(x1 >> 28);
  uint32_t x3 = (x1 & UINT32_C(0xfffffff));
  *out1 = x3;
  *out2 = (fiat_p251_uint1)(0x0 - x2);
}

/*
 * The function fiat_p251_addcarryx_u27 is an addition with carry.
 * Postconditions:
 *   out1 = (arg1 + arg2 + arg3) mod 2^27
 *   out2 = ⌊(arg1 + arg2 + arg3) / 2^27⌋
 *
 * Input Bounds:
 *   arg1: [0x0 ~> 0x1]
 *   arg2: [0x0 ~> 0x7ffffff]
 *   arg3: [0x0 ~> 0x7ffffff]
 * Output Bounds:
 *   out1: [0x0 ~> 0x7ffffff]
 *   out2: [0x0 ~> 0x1]
 */
static void fiat_p251_addcarryx_u27(uint32_t* out1, fiat_p251_uint1* out2, fiat_p251_uint1 arg1, uint32_t arg2, uint32_t arg3) {
  uint32_t x1 = ((arg1 + arg2) + arg3);
  uint32_t x2 = (x1 & UINT32_C(0x7ffffff));
  fiat_p251_uint1 x3 = (fiat_p251_uint1)(x1 >> 27);
  *out1 = x2;
  *out2 = x3;
}

/*
 * The function fiat_p251_subborrowx_u27 is a subtraction with borrow.
 * Postconditions:
 *   out1 = (-arg1 + arg2 + -arg3) mod 2^27
 *   out2 = -⌊(-arg1 + arg2 + -arg3) / 2^27⌋
 *
 * Input Bounds:
 *   arg1: [0x0 ~> 0x1]
 *   arg2: [0x0 ~> 0x7ffffff]
 *   arg3: [0x0 ~> 0x7ffffff]
 * Output Bounds:
 *   out1: [0x0 ~> 0x7ffffff]
 *   out2: [0x0 ~> 0x1]
 */
static void fiat_p251_subborrowx_u27(uint32_t* out1, fiat_p251_uint1* out2, fiat_p251_uint1 arg1, uint32_t arg2, uint32_t arg3) {
  int32_t x1 = ((int32_t)(arg2 - arg1) - (int32_t)arg3);
  fiat_p251_int1 x2 = (fiat_p251_int1)(x1 >> 27);
  uint32_t x3 = (x1 & UINT32_C(0x7ffffff));
  *out1 = x3;
  *out2 = (fiat_p251_uint1)(0x0 - x2);
}

/*
 * The function fiat_p251_cmovznz_u32 is a single-word conditional move.
 * Postconditions:
 *   out1 = (if arg1 = 0 then arg2 else arg3)
 *
 * Input Bounds:
 *   arg1: [0x0 ~> 0x1]
 *   arg2: [0x0 ~> 0xffffffff]
 *   arg3: [0x0 ~> 0xffffffff]
 * Output Bounds:
 *   out1: [0x0 ~> 0xffffffff]
 */
static void fiat_p251_cmovznz_u32(uint32_t* out1, fiat_p251_uint1 arg1, uint32_t arg2, uint32_t arg3) {
  fiat_p251_uint1 x1 = (!(!arg1));
  uint32_t x2 = ((fiat_p251_int1)(0x0 - x1) & UINT32_C(0xffffffff));
  uint32_t x3 = ((x2 & arg3) | ((~x2) & arg2));
  *out1 = x3;
}

/*
 * The function fiat_p251_carry_mul multiplies two field elements and reduces the result.
 * Postconditions:
 *   eval out1 mod m = (eval arg1 * eval arg2) mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x1a666664]]
 *   arg2: [[0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x1a666664]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x8cccccc]]
 */
static void fiat_p251_carry_mul(uint32_t out1[9], const uint32_t arg1[9], const uint32_t arg2[9]) {
  uint64_t x1 = ((arg1[8]) * ((uint64_t)(arg2[8]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x2 = ((arg1[8]) * ((uint64_t)(arg2[7]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x3 = ((arg1[8]) * ((uint64_t)(arg2[6]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x4 = ((arg1[8]) * ((uint64_t)(arg2[5]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x5 = ((arg1[8]) * ((uint64_t)(arg2[4]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x6 = ((arg1[8]) * ((uint64_t)(arg2[3]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x7 = ((arg1[8]) * ((uint64_t)(arg2[2]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x8 = ((arg1[8]) * ((uint64_t)(arg2[1]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x9 = ((arg1[7]) * ((uint64_t)(arg2[8]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x10 = ((arg1[7]) * ((uint64_t)(arg2[7]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x11 = ((arg1[7]) * ((uint64_t)(arg2[6]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x12 = ((arg1[7]) * ((uint64_t)(arg2[5]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x13 = ((arg1[7]) * ((uint64_t)(arg2[4]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x14 = ((arg1[7]) * ((uint64_t)(arg2[3]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x15 = ((arg1[7]) * ((uint64_t)(arg2[2]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x16 = ((arg1[6]) * ((uint64_t)(arg2[8]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x17 = ((arg1[6]) * ((uint64_t)(arg2[7]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x18 = ((arg1[6]) * ((uint64_t)(arg2[6]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x19 = ((arg1[6]) * ((uint64_t)(arg2[5]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x20 = ((arg1[6]) * ((uint64_t)(arg2[4]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x21 = ((arg1[6]) * ((uint64_t)(arg2[3]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x22 = ((arg1[5]) * ((uint64_t)(arg2[8]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x23 = ((arg1[5]) * ((uint64_t)(arg2[7]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x24 = ((arg1[5]) * ((uint64_t)(arg2[6]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x25 = ((arg1[5]) * ((uint64_t)(arg2[5]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x26 = ((arg1[5]) * ((uint64_t)(arg2[4]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x27 = ((arg1[4]) * ((uint64_t)(arg2[8]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x28 = ((arg1[4]) * ((uint64_t)(arg2[7]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x29 = ((arg1[4]) * ((uint64_t)(arg2[6]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x30 = ((arg1[4]) * ((uint64_t)(arg2[5]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x31 = ((arg1[3]) * ((uint64_t)(arg2[8]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x32 = ((arg1[3]) * ((uint64_t)(arg2[7]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x33 = ((arg1[3]) * ((uint64_t)(arg2[6]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x34 = ((arg1[2]) * ((uint64_t)(arg2[8]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x35 = ((arg1[2]) * ((uint64_t)(arg2[7]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x36 = ((arg1[1]) * ((uint64_t)(arg2[8]) * ((uint32_t)0x2 * 0x9)));
  uint64_t x37 = ((uint64_t)(arg1[8]) * (arg2[0]));
  uint64_t x38 = ((uint64_t)(arg1[7]) * (arg2[1]));
  uint64_t x39 = ((uint64_t)(arg1[7]) * (arg2[0]));
  uint64_t x40 = ((uint64_t)(arg1[6]) * (arg2[2]));
  uint64_t x41 = ((uint64_t)(arg1[6]) * (arg2[1]));
  uint64_t x42 = ((uint64_t)(arg1[6]) * (arg2[0]));
  uint64_t x43 = ((uint64_t)(arg1[5]) * (arg2[3]));
  uint64_t x44 = ((uint64_t)(arg1[5]) * (arg2[2]));
  uint64_t x45 = ((uint64_t)(arg1[5]) * (arg2[1]));
  uint64_t x46 = ((uint64_t)(arg1[5]) * (arg2[0]));
  uint64_t x47 = ((uint64_t)(arg1[4]) * (arg2[4]));
  uint64_t x48 = ((uint64_t)(arg1[4]) * (arg2[3]));
  uint64_t x49 = ((uint64_t)(arg1[4]) * (arg2[2]));
  uint64_t x50 = ((uint64_t)(arg1[4]) * (arg2[1]));
  uint64_t x51 = ((uint64_t)(arg1[4]) * (arg2[0]));
  uint64_t x52 = ((uint64_t)(arg1[3]) * (arg2[5]));
  uint64_t x53 = ((uint64_t)(arg1[3]) * (arg2[4]));
  uint64_t x54 = ((uint64_t)(arg1[3]) * (arg2[3]));
  uint64_t x55 = ((uint64_t)(arg1[3]) * (arg2[2]));
  uint64_t x56 = ((uint64_t)(arg1[3]) * (arg2[1]));
  uint64_t x57 = ((uint64_t)(arg1[3]) * (arg2[0]));
  uint64_t x58 = ((uint64_t)(arg1[2]) * (arg2[6]));
  uint64_t x59 = ((uint64_t)(arg1[2]) * (arg2[5]));
  uint64_t x60 = ((uint64_t)(arg1[2]) * (arg2[4]));
  uint64_t x61 = ((uint64_t)(arg1[2]) * (arg2[3]));
  uint64_t x62 = ((uint64_t)(arg1[2]) * (arg2[2]));
  uint64_t x63 = ((uint64_t)(arg1[2]) * (arg2[1]));
  uint64_t x64 = ((uint64_t)(arg1[2]) * (arg2[0]));
  uint64_t x65 = ((uint64_t)(arg1[1]) * (arg2[7]));
  uint64_t x66 = ((uint64_t)(arg1[1]) * (arg2[6]));
  uint64_t x67 = ((uint64_t)(arg1[1]) * (arg2[5]));
  uint64_t x68 = ((uint64_t)(arg1[1]) * (arg2[4]));
  uint64_t x69 = ((uint64_t)(arg1[1]) * (arg2[3]));
  uint64_t x70 = ((uint64_t)(arg1[1]) * (arg2[2]));
  uint64_t x71 = ((uint64_t)(arg1[1]) * (arg2[1]));
  uint64_t x72 = ((uint64_t)(arg1[1]) * (arg2[0]));
  uint64_t x73 = ((uint64_t)(arg1[0]) * (arg2[8]));
  uint64_t x74 = ((uint64_t)(arg1[0]) * (arg2[7]));
  uint64_t x75 = ((uint64_t)(arg1[0]) * (arg2[6]));
  uint64_t x76 = ((uint64_t)(arg1[0]) * (arg2[5]));
  uint64_t x77 = ((uint64_t)(arg1[0]) * (arg2[4]));
  uint64_t x78 = ((uint64_t)(arg1[0]) * (arg2[3]));
  uint64_t x79 = ((uint64_t)(arg1[0]) * (arg2[2]));
  uint64_t x80 = ((uint64_t)(arg1[0]) * (arg2[1]));
  uint64_t x81 = ((uint64_t)(arg1[0]) * (arg2[0]));
  fiat_p251_uint128 x82 = (x81 + (x36 + (x35 + (x33 + (x30 + (x26 + (x21 + ((fiat_p251_uint128)x15 + x8))))))));
  uint64_t x83 = (uint64_t)(x82 >> 28);
  uint32_t x84 = (uint32_t)(x82 & UINT32_C(0xfffffff));
  uint64_t x85 = (x73 + (x65 + (x58 + (x52 + (x47 + (x43 + (x40 + (x38 + x37))))))));
  uint64_t x86 = (x74 + (x66 + (x59 + (x53 + (x48 + (x44 + (x41 + (x39 + x1))))))));
  fiat_p251_uint128 x87 = (x75 + ((fiat_p251_uint128)x67 + (x60 + (x54 + (x49 + (x45 + (x42 + (x9 + x2))))))));
  fiat_p251_uint128 x88 = (x76 + (x68 + (x61 + (x55 + (x50 + (x46 + (x16 + ((fiat_p251_uint128)x10 + x3))))))));
  fiat_p251_uint128 x89 = (x77 + (x69 + (x62 + (x56 + (x51 + (x22 + (x17 + ((fiat_p251_uint128)x11 + x4))))))));
  fiat_p251_uint128 x90 = (x78 + (x70 + (x63 + (x57 + (x27 + (x23 + (x18 + ((fiat_p251_uint128)x12 + x5))))))));
  fiat_p251_uint128 x91 = (x79 + (x71 + (x64 + (x31 + (x28 + (x24 + (x19 + ((fiat_p251_uint128)x13 + x6))))))));
  fiat_p251_uint128 x92 = (x80 + (x72 + (x34 + (x32 + (x29 + (x25 + (x20 + ((fiat_p251_uint128)x14 + x7))))))));
  fiat_p251_uint128 x93 = (x83 + x92);
  uint64_t x94 = (uint64_t)(x93 >> 28);
  uint32_t x95 = (uint32_t)(x93 & UINT32_C(0xfffffff));
  fiat_p251_uint128 x96 = (x94 + x91);
  uint64_t x97 = (uint64_t)(x96 >> 28);
  uint32_t x98 = (uint32_t)(x96 & UINT32_C(0xfffffff));
  fiat_p251_uint128 x99 = (x97 + x90);
  uint64_t x100 = (uint64_t)(x99 >> 28);
  uint32_t x101 = (uint32_t)(x99 & UINT32_C(0xfffffff));
  fiat_p251_uint128 x102 = (x100 + x89);
  uint64_t x103 = (uint64_t)(x102 >> 28);
  uint32_t x104 = (uint32_t)(x102 & UINT32_C(0xfffffff));
  fiat_p251_uint128 x105 = (x103 + x88);
  uint64_t x106 = (uint64_t)(x105 >> 28);
  uint32_t x107 = (uint32_t)(x105 & UINT32_C(0xfffffff));
  fiat_p251_uint128 x108 = (x106 + x87);
  uint64_t x109 = (uint64_t)(x108 >> 28);
  uint32_t x110 = (uint32_t)(x108 & UINT32_C(0xfffffff));
  uint64_t x111 = (x109 + x86);
  uint64_t x112 = (x111 >> 28);
  uint32_t x113 = (uint32_t)(x111 & UINT32_C(0xfffffff));
  uint64_t x114 = (x112 + x85);
  uint64_t x115 = (x114 >> 27);
  uint32_t x116 = (uint32_t)(x114 & UINT32_C(0x7ffffff));
  uint64_t x117 = (x115 * (uint64_t)0x9);
  uint64_t x118 = (x84 + x117);
  uint32_t x119 = (uint32_t)(x118 >> 28);
  uint32_t x120 = (uint32_t)(x118 & UINT32_C(0xfffffff));
  uint32_t x121 = (x119 + x95);
  uint32_t x122 = (x121 >> 28);
  uint32_t x123 = (x121 & UINT32_C(0xfffffff));
  uint32_t x124 = (x122 + x98);
  out1[0] = x120;
  out1[1] = x123;
  out1[2] = x124;
  out1[3] = x101;
  out1[4] = x104;
  out1[5] = x107;
  out1[6] = x110;
  out1[7] = x113;
  out1[8] = x116;
}

/*
 * The function fiat_p251_carry_square squares a field element and reduces the result.
 * Postconditions:
 *   eval out1 mod m = (eval arg1 * eval arg1) mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x1a666664]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x8cccccc]]
 */
static void fiat_p251_carry_square(uint32_t out1[9], const uint32_t arg1[9]) {
  uint32_t x1 = ((arg1[8]) * (uint32_t)0x9);
  uint64_t x2 = (x1 * (uint64_t)0x2);
  uint32_t x3 = ((arg1[8]) * (uint32_t)0x2);
  uint64_t x4 = ((arg1[7]) * (uint64_t)0x9);
  uint64_t x5 = (x4 * (uint64_t)0x2);
  uint32_t x6 = ((arg1[7]) * (uint32_t)0x2);
  uint64_t x7 = ((arg1[6]) * (uint64_t)0x9);
  uint64_t x8 = (x7 * (uint64_t)0x2);
  uint32_t x9 = ((arg1[6]) * (uint32_t)0x2);
  uint64_t x10 = ((arg1[5]) * (uint64_t)0x9);
  uint64_t x11 = (x10 * (uint64_t)0x2);
  uint32_t x12 = ((arg1[5]) * (uint32_t)0x2);
  uint32_t x13 = ((arg1[4]) * (uint32_t)0x2);
  uint32_t x14 = ((arg1[3]) * (uint32_t)0x2);
  uint32_t x15 = ((arg1[2]) * (uint32_t)0x2);
  uint32_t x16 = ((arg1[1]) * (uint32_t)0x2);
  uint64_t x17 = ((arg1[8]) * (x1 * (uint64_t)0x2));
  uint64_t x18 = ((arg1[7]) * (x2 * (uint64_t)0x2));
  uint64_t x19 = ((arg1[7]) * (x4 * (uint64_t)0x2));
  uint64_t x20 = ((arg1[6]) * (x2 * (uint64_t)0x2));
  fiat_p251_uint128 x21 = ((arg1[6]) * (fiat_p251_uint128)(x5 * (uint64_t)0x2));
  uint64_t x22 = ((arg1[6]) * (x7 * (uint64_t)0x2));
  uint64_t x23 = ((arg1[5]) * (x2 * (uint64_t)0x2));
  fiat_p251_uint128 x24 = ((arg1[5]) * (fiat_p251_uint128)(x5 * (uint64_t)0x2));
  fiat_p251_uint128 x25 = ((arg1[5]) * (fiat_p251_uint128)(x8 * (uint64_t)0x2));
  uint64_t x26 = ((arg1[5]) * (x10 * (uint64_t)0x2));
  uint64_t x27 = ((arg1[4]) * (x2 * (uint64_t)0x2));
  fiat_p251_uint128 x28 = ((arg1[4]) * (fiat_p251_uint128)(x5 * (uint64_t)0x2));
  fiat_p251_uint128 x29 = ((arg1[4]) * (fiat_p251_uint128)(x8 * (uint64_t)0x2));
  fiat_p251_uint128 x30 = ((arg1[4]) * (fiat_p251_uint128)(x11 * (uint64_t)0x2));
  uint64_t x31 = ((uint64_t)(arg1[4]) * (arg1[4]));
  uint64_t x32 = ((arg1[3]) * (x2 * (uint64_t)0x2));
  fiat_p251_uint128 x33 = ((arg1[3]) * (fiat_p251_uint128)(x5 * (uint64_t)0x2));
  fiat_p251_uint128 x34 = ((arg1[3]) * (fiat_p251_uint128)(x8 * (uint64_t)0x2));
  uint64_t x35 = ((uint64_t)(arg1[3]) * x12);
  uint64_t x36 = ((uint64_t)(arg1[3]) * x13);
  uint64_t x37 = ((uint64_t)(arg1[3]) * (arg1[3]));
  uint64_t x38 = ((arg1[2]) * (x2 * (uint64_t)0x2));
  fiat_p251_uint128 x39 = ((arg1[2]) * (fiat_p251_uint128)(x5 * (uint64_t)0x2));
  uint64_t x40 = ((uint64_t)(arg1[2]) * x9);
  uint64_t x41 = ((uint64_t)(arg1[2]) * x12);
  uint64_t x42 = ((uint64_t)(arg1[2]) * x13);
  uint64_t x43 = ((uint64_t)(arg1[2]) * x14);
  uint64_t x44 = ((uint64_t)(arg1[2]) * (arg1[2]));
  uint64_t x45 = ((arg1[1]) * (x2 * (uint64_t)0x2));
  uint64_t x46 = ((uint64_t)(arg1[1]) * x6);
  uint64_t x47 = ((uint64_t)(arg1[1]) * x9);
  uint64_t x48 = ((uint64_t)(arg1[1]) * x12);
  uint64_t x49 = ((uint64_t)(arg1[1]) * x13);
  uint64_t x50 = ((uint64_t)(arg1[1]) * x14);
  uint64_t x51 = ((uint64_t)(arg1[1]) * x15);
  uint64_t x52 = ((uint64_t)(arg1[1]) * (arg1[1]));
  uint64_t x53 = ((uint64_t)(arg1[0]) * x3);
  uint64_t x54 = ((uint64_t)(arg1[0]) * x6);
  uint64_t x55 = ((uint64_t)(arg1[0]) * x9);
  uint64_t x56 = ((uint64_t)(arg1[0]) * x12);
  uint64_t x57 = ((uint64_t)(arg1[0]) * x13);
  uint64_t x58 = ((uint64_t)(arg1[0]) * x14);
  uint64_t x59 = ((uint64_t)(arg1[0]) * x15);
  uint64_t x60 = ((uint64_t)(arg1[0]) * x16);
  uint64_t x61 = ((uint64_t)(arg1[0]) * (arg1[0]));
  fiat_p251_uint128 x62 = (x61 + (x45 + (x39 + (x34 + x30))));
  uint64_t x63 = (uint64_t)(x62 >> 28);
  uint32_t x64 = (uint32_t)(x62 & UINT32_C(0xfffffff));
  uint64_t x65 = (x53 + (x46 + (x40 + (x35 + x31))));
  uint64_t x66 = (x54 + (x47 + (x41 + (x36 + x17))));
  fiat_p251_uint128 x67 = ((fiat_p251_uint128)x55 + (x48 + (x42 + (x37 + x18))));
  fiat_p251_uint128 x68 = (x56 + (x49 + (x43 + ((fiat_p251_uint128)x20 + x19))));
  fiat_p251_uint128 x69 = (x57 + (x50 + (x44 + (x23 + x21))));
  fiat_p251_uint128 x70 = (x58 + (x51 + (x27 + (x24 + x22))));
  fiat_p251_uint128 x71 = (x59 + (x52 + (x32 + (x28 + x25))));
  fiat_p251_uint128 x72 = (x60 + (x38 + (x33 + (x29 + x26))));
  fiat_p251_uint128 x73 = (x63 + x72);
  uint64_t x74 = (uint64_t)(x73 >> 28);
  uint32_t x75 = (uint32_t)(x73 & UINT32_C(0xfffffff));
  fiat_p251_uint128 x76 = (x74 + x71);
  uint64_t x77 = (uint64_t)(x76 >> 28);
  uint32_t x78 = (uint32_t)(x76 & UINT32_C(0xfffffff));
  fiat_p251_uint128 x79 = (x77 + x70);
  uint64_t x80 = (uint64_t)(x79 >> 28);
  uint32_t x81 = (uint32_t)(x79 & UINT32_C(0xfffffff));
  fiat_p251_uint128 x82 = (x80 + x69);
  uint64_t x83 = (uint64_t)(x82 >> 28);
  uint32_t x84 = (uint32_t)(x82 & UINT32_C(0xfffffff));
  fiat_p251_uint128 x85 = (x83 + x68);
  uint64_t x86 = (uint64_t)(x85 >> 28);
  uint32_t x87 = (uint32_t)(x85 & UINT32_C(0xfffffff));
  fiat_p251_uint128 x88 = (x86 + x67);
  uint64_t x89 = (uint64_t)(x88 >> 28);
  uint32_t x90 = (uint32_t)(x88 & UINT32_C(0xfffffff));
  uint64_t x91 = (x89 + x66);
  uint64_t x92 = (x91 >> 28);
  uint32_t x93 = (uint32_t)(x91 & UINT32_C(0xfffffff));
  uint64_t x94 = (x92 + x65);
  uint64_t x95 = (x94 >> 27);
  uint32_t x96 = (uint32_t)(x94 & UINT32_C(0x7ffffff));
  uint64_t x97 = (x95 * (uint64_t)0x9);
  uint64_t x98 = (x64 + x97);
  uint32_t x99 = (uint32_t)(x98 >> 28);
  uint32_t x100 = (uint32_t)(x98 & UINT32_C(0xfffffff));
  uint32_t x101 = (x99 + x75);
  uint32_t x102 = (x101 >> 28);
  uint32_t x103 = (x101 & UINT32_C(0xfffffff));
  uint32_t x104 = (x102 + x78);
  out1[0] = x100;
  out1[1] = x103;
  out1[2] = x104;
  out1[3] = x81;
  out1[4] = x84;
  out1[5] = x87;
  out1[6] = x90;
  out1[7] = x93;
  out1[8] = x96;
}

/*
 * The function fiat_p251_carry reduces a field element.
 * Postconditions:
 *   eval out1 mod m = eval arg1 mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x1a666664]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x8cccccc]]
 */
static void fiat_p251_carry(uint32_t out1[9], const uint32_t arg1[9]) {
  uint32_t x1 = (arg1[0]);
  uint32_t x2 = ((x1 >> 28) + (arg1[1]));
  uint32_t x3 = ((x2 >> 28) + (arg1[2]));
  uint32_t x4 = ((x3 >> 28) + (arg1[3]));
  uint32_t x5 = ((x4 >> 28) + (arg1[4]));
  uint32_t x6 = ((x5 >> 28) + (arg1[5]));
  uint32_t x7 = ((x6 >> 28) + (arg1[6]));
  uint32_t x8 = ((x7 >> 28) + (arg1[7]));
  uint32_t x9 = ((x8 >> 28) + (arg1[8]));
  uint32_t x10 = ((x1 & UINT32_C(0xfffffff)) + ((x9 >> 27) * (uint32_t)0x9));
  uint32_t x11 = ((x10 >> 28) + (x2 & UINT32_C(0xfffffff)));
  uint32_t x12 = (x10 & UINT32_C(0xfffffff));
  uint32_t x13 = (x11 & UINT32_C(0xfffffff));
  uint32_t x14 = ((x11 >> 28) + (x3 & UINT32_C(0xfffffff)));
  uint32_t x15 = (x4 & UINT32_C(0xfffffff));
  uint32_t x16 = (x5 & UINT32_C(0xfffffff));
  uint32_t x17 = (x6 & UINT32_C(0xfffffff));
  uint32_t x18 = (x7 & UINT32_C(0xfffffff));
  uint32_t x19 = (x8 & UINT32_C(0xfffffff));
  uint32_t x20 = (x9 & UINT32_C(0x7ffffff));
  out1[0] = x12;
  out1[1] = x13;
  out1[2] = x14;
  out1[3] = x15;
  out1[4] = x16;
  out1[5] = x17;
  out1[6] = x18;
  out1[7] = x19;
  out1[8] = x20;
}

/*
 * The function fiat_p251_add adds two field elements.
 * Postconditions:
 *   eval out1 mod m = (eval arg1 + eval arg2) mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x8cccccc]]
 *   arg2: [[0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x8cccccc]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x1a666664]]
 */
static void fiat_p251_add(uint32_t out1[9], const uint32_t arg1[9], const uint32_t arg2[9]) {
  uint32_t x1 = ((arg1[0]) + (arg2[0]));
  uint32_t x2 = ((arg1[1]) + (arg2[1]));
  uint32_t x3 = ((arg1[2]) + (arg2[2]));
  uint32_t x4 = ((arg1[3]) + (arg2[3]));
  uint32_t x5 = ((arg1[4]) + (arg2[4]));
  uint32_t x6 = ((arg1[5]) + (arg2[5]));
  uint32_t x7 = ((arg1[6]) + (arg2[6]));
  uint32_t x8 = ((arg1[7]) + (arg2[7]));
  uint32_t x9 = ((arg1[8]) + (arg2[8]));
  out1[0] = x1;
  out1[1] = x2;
  out1[2] = x3;
  out1[3] = x4;
  out1[4] = x5;
  out1[5] = x6;
  out1[6] = x7;
  out1[7] = x8;
  out1[8] = x9;
}

/*
 * The function fiat_p251_sub subtracts two field elements.
 * Postconditions:
 *   eval out1 mod m = (eval arg1 - eval arg2) mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x8cccccc]]
 *   arg2: [[0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x8cccccc]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x1a666664]]
 */
static void fiat_p251_sub(uint32_t out1[9], const uint32_t arg1[9], const uint32_t arg2[9]) {
  uint32_t x1 = ((UINT32_C(0x1fffffee) + (arg1[0])) - (arg2[0]));
  uint32_t x2 = ((UINT32_C(0x1ffffffe) + (arg1[1])) - (arg2[1]));
  uint32_t x3 = ((UINT32_C(0x1ffffffe) + (arg1[2])) - (arg2[2]));
  uint32_t x4 = ((UINT32_C(0x1ffffffe) + (arg1[3])) - (arg2[3]));
  uint32_t x5 = ((UINT32_C(0x1ffffffe) + (arg1[4])) - (arg2[4]));
  uint32_t x6 = ((UINT32_C(0x1ffffffe) + (arg1[5])) - (arg2[5]));
  uint32_t x7 = ((UINT32_C(0x1ffffffe) + (arg1[6])) - (arg2[6]));
  uint32_t x8 = ((UINT32_C(0x1ffffffe) + (arg1[7])) - (arg2[7]));
  uint32_t x9 = ((UINT32_C(0xffffffe) + (arg1[8])) - (arg2[8]));
  out1[0] = x1;
  out1[1] = x2;
  out1[2] = x3;
  out1[3] = x4;
  out1[4] = x5;
  out1[5] = x6;
  out1[6] = x7;
  out1[7] = x8;
  out1[8] = x9;
}

/*
 * The function fiat_p251_opp negates a field element.
 * Postconditions:
 *   eval out1 mod m = -eval arg1 mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x8cccccc]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x34cccccb], [0x0 ~> 0x1a666664]]
 */
static void fiat_p251_opp(uint32_t out1[9], const uint32_t arg1[9]) {
  uint32_t x1 = (UINT32_C(0x1fffffee) - (arg1[0]));
  uint32_t x2 = (UINT32_C(0x1ffffffe) - (arg1[1]));
  uint32_t x3 = (UINT32_C(0x1ffffffe) - (arg1[2]));
  uint32_t x4 = (UINT32_C(0x1ffffffe) - (arg1[3]));
  uint32_t x5 = (UINT32_C(0x1ffffffe) - (arg1[4]));
  uint32_t x6 = (UINT32_C(0x1ffffffe) - (arg1[5]));
  uint32_t x7 = (UINT32_C(0x1ffffffe) - (arg1[6]));
  uint32_t x8 = (UINT32_C(0x1ffffffe) - (arg1[7]));
  uint32_t x9 = (UINT32_C(0xffffffe) - (arg1[8]));
  out1[0] = x1;
  out1[1] = x2;
  out1[2] = x3;
  out1[3] = x4;
  out1[4] = x5;
  out1[5] = x6;
  out1[6] = x7;
  out1[7] = x8;
  out1[8] = x9;
}

/*
 * The function fiat_p251_selectznz is a multi-limb conditional select.
 * Postconditions:
 *   eval out1 = (if arg1 = 0 then eval arg2 else eval arg3)
 *
 * Input Bounds:
 *   arg1: [0x0 ~> 0x1]
 *   arg2: [[0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff]]
 *   arg3: [[0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff], [0x0 ~> 0xffffffff]]
 */
static void fiat_p251_selectznz(uint32_t out1[9], fiat_p251_uint1 arg1, const uint32_t arg2[9], const uint32_t arg3[9]) {
  uint32_t x1;
  fiat_p251_cmovznz_u32(&x1, arg1, (arg2[0]), (arg3[0]));
  uint32_t x2;
  fiat_p251_cmovznz_u32(&x2, arg1, (arg2[1]), (arg3[1]));
  uint32_t x3;
  fiat_p251_cmovznz_u32(&x3, arg1, (arg2[2]), (arg3[2]));
  uint32_t x4;
  fiat_p251_cmovznz_u32(&x4, arg1, (arg2[3]), (arg3[3]));
  uint32_t x5;
  fiat_p251_cmovznz_u32(&x5, arg1, (arg2[4]), (arg3[4]));
  uint32_t x6;
  fiat_p251_cmovznz_u32(&x6, arg1, (arg2[5]), (arg3[5]));
  uint32_t x7;
  fiat_p251_cmovznz_u32(&x7, arg1, (arg2[6]), (arg3[6]));
  uint32_t x8;
  fiat_p251_cmovznz_u32(&x8, arg1, (arg2[7]), (arg3[7]));
  uint32_t x9;
  fiat_p251_cmovznz_u32(&x9, arg1, (arg2[8]), (arg3[8]));
  out1[0] = x1;
  out1[1] = x2;
  out1[2] = x3;
  out1[3] = x4;
  out1[4] = x5;
  out1[5] = x6;
  out1[6] = x7;
  out1[7] = x8;
  out1[8] = x9;
}

/*
 * The function fiat_p251_to_bytes serializes a field element to bytes in little-endian order.
 * Postconditions:
 *   out1 = map (λ x, ⌊((eval arg1 mod m) mod 2^(8 * (x + 1))) / 2^(8 * x)⌋) [0..31]
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x8cccccc]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0x7]]
 */
static void fiat_p251_to_bytes(uint8_t out1[32], const uint32_t arg1[9]) {
  uint32_t x1;
  fiat_p251_uint1 x2;
  fiat_p251_subborrowx_u28(&x1, &x2, 0x0, (arg1[0]), UINT32_C(0xffffff7));
  uint32_t x3;
  fiat_p251_uint1 x4;
  fiat_p251_subborrowx_u28(&x3, &x4, x2, (arg1[1]), UINT32_C(0xfffffff));
  uint32_t x5;
  fiat_p251_uint1 x6;
  fiat_p251_subborrowx_u28(&x5, &x6, x4, (arg1[2]), UINT32_C(0xfffffff));
  uint32_t x7;
  fiat_p251_uint1 x8;
  fiat_p251_subborrowx_u28(&x7, &x8, x6, (arg1[3]), UINT32_C(0xfffffff));
  uint32_t x9;
  fiat_p251_uint1 x10;
  fiat_p251_subborrowx_u28(&x9, &x10, x8, (arg1[4]), UINT32_C(0xfffffff));
  uint32_t x11;
  fiat_p251_uint1 x12;
  fiat_p251_subborrowx_u28(&x11, &x12, x10, (arg1[5]), UINT32_C(0xfffffff));
  uint32_t x13;
  fiat_p251_uint1 x14;
  fiat_p251_subborrowx_u28(&x13, &x14, x12, (arg1[6]), UINT32_C(0xfffffff));
  uint32_t x15;
  fiat_p251_uint1 x16;
  fiat_p251_subborrowx_u28(&x15, &x16, x14, (arg1[7]), UINT32_C(0xfffffff));
  uint32_t x17;
  fiat_p251_uint1 x18;
  fiat_p251_subborrowx_u27(&x17, &x18, x16, (arg1[8]), UINT32_C(0x7ffffff));
  uint32_t x19;
  fiat_p251_cmovznz_u32(&x19, x18, 0x0, UINT32_C(0xffffffff));
  uint32_t x20;
  fiat_p251_uint1 x21;
  fiat_p251_addcarryx_u28(&x20, &x21, 0x0, x1, (x19 & UINT32_C(0xffffff7)));
  uint32_t x22;
  fiat_p251_uint1 x23;
  fiat_p251_addcarryx_u28(&x22, &x23, x21, x3, (x19 & UINT32_C(0xfffffff)));
  uint32_t x24;
  fiat_p251_uint1 x25;
  fiat_p251_addcarryx_u28(&x24, &x25, x23, x5, (x19 & UINT32_C(0xfffffff)));
  uint32_t x26;
  fiat_p251_uint1 x27;
  fiat_p251_addcarryx_u28(&x26, &x27, x25, x7, (x19 & UINT32_C(0xfffffff)));
  uint32_t x28;
  fiat_p251_uint1 x29;
  fiat_p251_addcarryx_u28(&x28, &x29, x27, x9, (x19 & UINT32_C(0xfffffff)));
  uint32_t x30;
  fiat_p251_uint1 x31;
  fiat_p251_addcarryx_u28(&x30, &x31, x29, x11, (x19 & UINT32_C(0xfffffff)));
  uint32_t x32;
  fiat_p251_uint1 x33;
  fiat_p251_addcarryx_u28(&x32, &x33, x31, x13, (x19 & UINT32_C(0xfffffff)));
  uint32_t x34;
  fiat_p251_uint1 x35;
  fiat_p251_addcarryx_u28(&x34, &x35, x33, x15, (x19 & UINT32_C(0xfffffff)));
  uint32_t x36;
  fiat_p251_uint1 x37;
  fiat_p251_addcarryx_u27(&x36, &x37, x35, x17, (x19 & UINT32_C(0x7ffffff)));
  uint32_t x38 = (x34 << 4);
  uint32_t x39 = (x30 << 4);
  uint32_t x40 = (x26 << 4);
  uint32_t x41 = (x22 << 4);
  uint32_t x42 = (x20 >> 8);
  uint8_t x43 = (uint8_t)(x20 & UINT8_C(0xff));
  uint32_t x44 = (x42 >> 8);
  uint8_t x45 = (uint8_t)(x42 & UINT8_C(0xff));
  uint8_t x46 = (uint8_t)(x44 >> 8);
  uint8_t x47 = (uint8_t)(x44 & UINT8_C(0xff));
  uint32_t x48 = (x46 + x41);
  uint32_t x49 = (x48 >> 8);
  uint8_t x50 = (uint8_t)(x48 & UINT8_C(0xff));
  uint32_t x51 = (x49 >> 8);
  uint8_t x52 = (uint8_t)(x49 & UINT8_C(0xff));
  uint8_t x53 = (uint8_t)(x51 >> 8);
  uint8_t x54 = (uint8_t)(x51 & UINT8_C(0xff));
  uint8_t x55 = (uint8_t)(x53 & UINT8_C(0xff));
  uint32_t x56 = (x24 >> 8);
  uint8_t x57 = (uint8_t)(x24 & UINT8_C(0xff));
  uint32_t x58 = (x56 >> 8);
  uint8_t x59 = (uint8_t)(x56 & UINT8_C(0xff));
  uint8_t x60 = (uint8_t)(x58 >> 8);
  uint8_t x61 = (uint8_t)(x58 & UINT8_C(0xff));
  uint32_t x62 = (x60 + x40);
  uint32_t x63 = (x62 >> 8);
  uint8_t x64 = (uint8_t)(x62 & UINT8_C(0xff));
  uint32_t x65 = (x63 >> 8);
  uint8_t x66 = (uint8_t)(x63 & UINT8_C(0xff));
  uint8_t x67 = (uint8_t)(x65 >> 8);
  uint8_t x68 = (uint8_t)(x65 & UINT8_C(0xff));
  uint8_t x69 = (uint8_t)(x67 & UINT8_C(0xff));
  uint32_t x70 = (x28 >> 8);
  uint8_t x71 = (uint8_t)(x28 & UINT8_C(0xff));
  uint32_t x72 = (x70 >> 8);
  uint8_t x73 = (uint8_t)(x70 & UINT8_C(0xff));
  uint8_t x74 = (uint8_t)(x72 >> 8);
  uint8_t x75 = (uint8_t)(x72 & UINT8_C(0xff));
  uint32_t x76 = (x74 + x39);
  uint32_t x77 = (x76 >> 8);
  uint8_t x78 = (uint8_t)(x76 & UINT8_C(0xff));
  uint32_t x79 = (x77 >> 8);
  uint8_t x80 = (uint8_t)(x77 & UINT8_C(0xff));
  uint8_t x81 = (uint8_t)(x79 >> 8);
  uint8_t x82 = (uint8_t)(x79 & UINT8_C(0xff));
  uint8_t x83 = (uint8_t)(x81 & UINT8_C(0xff));
  uint32_t x84 = (x32 >> 8);
  uint8_t x85 = (uint8_t)(x32 & UINT8_C(0xff));
  uint32_t x86 = (x84 >> 8);
  uint8_t x87 = (uint8_t)(x84 & UINT8_C(0xff));
  uint8_t x88 = (uint8_t)(x86 >> 8);
  uint8_t x89 = (uint8_t)(x86 & UINT8_C(0xff));
  uint32_t x90 = (x88 + x38);
  uint32_t x91 = (x90 >> 8);
  uint8_t x92 = (uint8_t)(x90 & UINT8_C(0xff));
  uint32_t x93 = (x91 >> 8);
  uint8_t x94 = (uint8_t)(x91 & UINT8_C(0xff));
  uint8_t x95 = (uint8_t)(x93 >> 8);
  uint8_t x96 = (uint8_t)(x93 & UINT8_C(0xff));
  uint8_t x97 = (uint8_t)(x95 & UINT8_C(0xff));
  uint32_t x98 = (x36 >> 8);
  uint8_t x99 = (uint8_t)(x36 & UINT8_C(0xff));
  uint32_t x100 = (x98 >> 8);
  uint8_t x101 = (uint8_t)(x98 & UINT8_C(0xff));
  uint8_t x102 = (uint8_t)(x100 >> 8);
  uint8_t x103 = (uint8_t)(x100 & UINT8_C(0xff));
  out1[0] = x43;
  out1[1] = x45;
  out1[2] = x47;
  out1[3] = x50;
  out1[4] = x52;
  out1[5] = x54;
  out1[6] = x55;
  out1[7] = x57;
  out1[8] = x59;
  out1[9] = x61;
  out1[10] = x64;
  out1[11] = x66;
  out1[12] = x68;
  out1[13] = x69;
  out1[14] = x71;
  out1[15] = x73;
  out1[16] = x75;
  out1[17] = x78;
  out1[18] = x80;
  out1[19] = x82;
  out1[20] = x83;
  out1[21] = x85;
  out1[22] = x87;
  out1[23] = x89;
  out1[24] = x92;
  out1[25] = x94;
  out1[26] = x96;
  out1[27] = x97;
  out1[28] = x99;
  out1[29] = x101;
  out1[30] = x103;
  out1[31] = x102;
}

/*
 * The function fiat_p251_from_bytes deserializes a field element from bytes in little-endian order.
 * Postconditions:
 *   eval out1 mod m = bytes_eval arg1 mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0x7]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x11999999], [0x0 ~> 0x8cccccc]]
 */
static void fiat_p251_from_bytes(uint32_t out1[9], const uint8_t arg1[32]) {
  uint32_t x1 = ((uint32_t)(arg1[31]) << 24);
  uint32_t x2 = ((uint32_t)(arg1[30]) << 16);
  uint32_t x3 = ((uint32_t)(arg1[29]) << 8);
  uint8_t x4 = (arg1[28]);
  uint32_t x5 = ((uint32_t)(arg1[27]) << 20);
  uint32_t x6 = ((uint32_t)(arg1[26]) << 12);
  uint32_t x7 = ((uint32_t)(arg1[25]) << 4);
  uint32_t x8 = ((uint32_t)(arg1[24]) << 24);
  uint32_t x9 = ((uint32_t)(arg1[23]) << 16);
  uint32_t x10 = ((uint32_t)(arg1[22]) << 8);
  uint8_t x11 = (arg1[21]);
  uint32_t x12 = ((uint32_t)(arg1[20]) << 20);
  uint32_t x13 = ((uint32_t)(arg1[19]) << 12);
  uint32_t x14 = ((uint32_t)(arg1[18]) << 4);
  uint32_t x15 = ((uint32_t)(arg1[17]) << 24);
  uint32_t x16 = ((uint32_t)(arg1[16]) << 16);
  uint32_t x17 = ((uint32_t)(arg1[15]) << 8);
  uint8_t x18 = (arg1[14]);
  uint32_t x19 = ((uint32_t)(arg1[13]) << 20);
  uint32_t x20 = ((uint32_t)(arg1[12]) << 12);
  uint32_t x21 = ((uint32_t)(arg1[11]) << 4);
  uint32_t x22 = ((uint32_t)(arg1[10]) << 24);
  uint32_t x23 = ((uint32_t)(arg1[9]) << 16);
  uint32_t x24 = ((uint32_t)(arg1[8]) << 8);
  uint8_t x25 = (arg1[7]);
  uint32_t x26 = ((uint32_t)(arg1[6]) << 20);
  uint32_t x27 = ((uint32_t)(arg1[5]) << 12);
  uint32_t x28 = ((uint32_t)(arg1[4]) << 4);
  uint32_t x29 = ((uint32_t)(arg1[3]) << 24);
  uint32_t x30 = ((uint32_t)(arg1[2]) << 16);
  uint32_t x31 = ((uint32_t)(arg1[1]) << 8);
  uint8_t x32 = (arg1[0]);
  uint32_t x33 = (x32 + (x31 + (x30 + x29)));
  uint8_t x34 = (uint8_t)(x33 >> 28);
  uint32_t x35 = (x33 & UINT32_C(0xfffffff));
  uint32_t x36 = (x4 + (x3 + (x2 + x1)));
  uint32_t x37 = (x7 + (x6 + x5));
  uint32_t x38 = (x11 + (x10 + (x9 + x8)));
  uint32_t x39 = (x14 + (x13 + x12));
  uint32_t x40 = (x18 + (x17 + (x16 + x15)));
  uint32_t x41 = (x21 + (x20 + x19));
  uint32_t x42 = (x25 + (x24 + (x23 + x22)));
  uint32_t x43 = (x28 + (x27 + x26));
  uint32_t x44 = (x34 + x43);
  uint32_t x45 = (x44 & UINT32_C(0xfffffff));
  uint8_t x46 = (uint8_t)(x42 >> 28);
  uint32_t x47 = (x42 & UINT32_C(0xfffffff));
  uint32_t x48 = (x46 + x41);
  uint32_t x49 = (x48 & UINT32_C(0xfffffff));
  uint8_t x50 = (uint8_t)(x40 >> 28);
  uint32_t x51 = (x40 & UINT32_C(0xfffffff));
  uint32_t x52 = (x50 + x39);
  uint32_t x53 = (x52 & UINT32_C(0xfffffff));
  uint8_t x54 = (uint8_t)(x38 >> 28);
  uint32_t x55 = (x38 & UINT32_C(0xfffffff));
  uint32_t x56 = (x54 + x37);
  uint32_t x57 = (x56 & UINT32_C(0xfffffff));
  out1[0] = x35;
  out1[1] = x45;
  out1[2] = x47;
  out1[3] = x49;
  out1[4] = x51;
  out1[5] = x53;
  out1[6] = x55;
  out1[7] = x57;
  out1[8] = x36;
}


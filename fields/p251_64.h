/* Autogenerated */
/* curve description: p251 */
/* requested operations: (all) */
/* n = 4 (from "4") */
/* s-c = 2^251 - [(1, 9)] (from "2^251 - 9") */
/* machine_wordsize = 64 (from "64") */

/* Computed values: */
/* carry_chain = [0, 1, 2, 3, 0, 1] */

#include <stdint.h>
typedef unsigned char fiat_p251_uint1;
typedef signed char fiat_p251_int1;
typedef signed __int128 fiat_p251_int128;
typedef unsigned __int128 fiat_p251_uint128;

#if (-1 & 3) != 3
#error "This code only works on a two's complement system"
#endif


/*
 * The function fiat_p251_addcarryx_u62 is an addition with carry.
 * Postconditions:
 *   out1 = (arg1 + arg2 + arg3) mod 2^62
 *   out2 = ⌊(arg1 + arg2 + arg3) / 2^62⌋
 *
 * Input Bounds:
 *   arg1: [0x0 ~> 0x1]
 *   arg2: [0x0 ~> 0x3fffffffffffffff]
 *   arg3: [0x0 ~> 0x3fffffffffffffff]
 * Output Bounds:
 *   out1: [0x0 ~> 0x3fffffffffffffff]
 *   out2: [0x0 ~> 0x1]
 */
static void fiat_p251_addcarryx_u62(uint64_t* out1, fiat_p251_uint1* out2, fiat_p251_uint1 arg1, uint64_t arg2, uint64_t arg3) {
  uint64_t x1 = ((arg1 + arg2) + arg3);
  uint64_t x2 = (x1 & UINT64_C(0x3fffffffffffffff));
  fiat_p251_uint1 x3 = (fiat_p251_uint1)(x1 >> 62);
  *out1 = x2;
  *out2 = x3;
}

/*
 * The function fiat_p251_subborrowx_u62 is a subtraction with borrow.
 * Postconditions:
 *   out1 = (-arg1 + arg2 + -arg3) mod 2^62
 *   out2 = -⌊(-arg1 + arg2 + -arg3) / 2^62⌋
 *
 * Input Bounds:
 *   arg1: [0x0 ~> 0x1]
 *   arg2: [0x0 ~> 0x3fffffffffffffff]
 *   arg3: [0x0 ~> 0x3fffffffffffffff]
 * Output Bounds:
 *   out1: [0x0 ~> 0x3fffffffffffffff]
 *   out2: [0x0 ~> 0x1]
 */
static void fiat_p251_subborrowx_u62(uint64_t* out1, fiat_p251_uint1* out2, fiat_p251_uint1 arg1, uint64_t arg2, uint64_t arg3) {
  int64_t x1 = ((int64_t)(arg2 - (int64_t)arg1) - (int64_t)arg3);
  fiat_p251_int1 x2 = (fiat_p251_int1)(x1 >> 62);
  uint64_t x3 = (x1 & UINT64_C(0x3fffffffffffffff));
  *out1 = x3;
  *out2 = (fiat_p251_uint1)(0x0 - x2);
}

/*
 * The function fiat_p251_addcarryx_u63 is an addition with carry.
 * Postconditions:
 *   out1 = (arg1 + arg2 + arg3) mod 2^63
 *   out2 = ⌊(arg1 + arg2 + arg3) / 2^63⌋
 *
 * Input Bounds:
 *   arg1: [0x0 ~> 0x1]
 *   arg2: [0x0 ~> 0x7fffffffffffffff]
 *   arg3: [0x0 ~> 0x7fffffffffffffff]
 * Output Bounds:
 *   out1: [0x0 ~> 0x7fffffffffffffff]
 *   out2: [0x0 ~> 0x1]
 */
static void fiat_p251_addcarryx_u63(uint64_t* out1, fiat_p251_uint1* out2, fiat_p251_uint1 arg1, uint64_t arg2, uint64_t arg3) {
  uint64_t x1 = ((arg1 + arg2) + arg3);
  uint64_t x2 = (x1 & UINT64_C(0x7fffffffffffffff));
  fiat_p251_uint1 x3 = (fiat_p251_uint1)(x1 >> 63);
  *out1 = x2;
  *out2 = x3;
}

/*
 * The function fiat_p251_subborrowx_u63 is a subtraction with borrow.
 * Postconditions:
 *   out1 = (-arg1 + arg2 + -arg3) mod 2^63
 *   out2 = -⌊(-arg1 + arg2 + -arg3) / 2^63⌋
 *
 * Input Bounds:
 *   arg1: [0x0 ~> 0x1]
 *   arg2: [0x0 ~> 0x7fffffffffffffff]
 *   arg3: [0x0 ~> 0x7fffffffffffffff]
 * Output Bounds:
 *   out1: [0x0 ~> 0x7fffffffffffffff]
 *   out2: [0x0 ~> 0x1]
 */
static void fiat_p251_subborrowx_u63(uint64_t* out1, fiat_p251_uint1* out2, fiat_p251_uint1 arg1, uint64_t arg2, uint64_t arg3) {
  int64_t x1 = ((int64_t)(arg2 - (int64_t)arg1) - (int64_t)arg3);
  fiat_p251_int1 x2 = (fiat_p251_int1)((fiat_p251_int128)x1 >> 63);
  uint64_t x3 = (x1 & UINT64_C(0x7fffffffffffffff));
  *out1 = x3;
  *out2 = (fiat_p251_uint1)(0x0 - x2);
}

/*
 * The function fiat_p251_cmovznz_u64 is a single-word conditional move.
 * Postconditions:
 *   out1 = (if arg1 = 0 then arg2 else arg3)
 *
 * Input Bounds:
 *   arg1: [0x0 ~> 0x1]
 *   arg2: [0x0 ~> 0xffffffffffffffff]
 *   arg3: [0x0 ~> 0xffffffffffffffff]
 * Output Bounds:
 *   out1: [0x0 ~> 0xffffffffffffffff]
 */
static void fiat_p251_cmovznz_u64(uint64_t* out1, fiat_p251_uint1 arg1, uint64_t arg2, uint64_t arg3) {
  fiat_p251_uint1 x1 = (!(!arg1));
  uint64_t x2 = ((fiat_p251_int1)(0x0 - x1) & UINT64_C(0xffffffffffffffff));
  uint64_t x3 = ((x2 & arg3) | ((~x2) & arg2));
  *out1 = x3;
}

/*
 * The function fiat_p251_carry_mul multiplies two field elements and reduces the result.
 * Postconditions:
 *   eval out1 mod m = (eval arg1 * eval arg2) mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0xd333333333333332]]
 *   arg2: [[0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0xd333333333333332]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x4666666666666666]]
 */
static void fiat_p251_carry_mul(uint64_t out1[4], const fiat_p251_uint128 arg1[4], const fiat_p251_uint128 arg2[4]) {
  fiat_p251_uint256 x1 = ((uint64_t)(arg1[3]) * (fiat_p251_uint256)((fiat_p251_uint128)(uint64_t)(arg2[3]) * ((uint64_t)0x2 * 0x9)));
  fiat_p251_uint256 x2 = ((uint64_t)(arg1[3]) * (fiat_p251_uint256)((arg2[2]) * ((uint64_t)0x2 * 0x9)));
  fiat_p251_uint256 x3 = ((uint64_t)(arg1[3]) * (fiat_p251_uint256)((arg2[1]) * ((uint64_t)0x2 * 0x9)));
  fiat_p251_uint256 x4 = ((fiat_p251_uint256)(arg1[2]) * ((fiat_p251_uint128)(uint64_t)(arg2[3]) * ((uint64_t)0x2 * 0x9)));
  fiat_p251_uint256 x5 = ((fiat_p251_uint256)(arg1[2]) * ((arg2[2]) * ((uint64_t)0x2 * 0x9)));
  fiat_p251_uint256 x6 = ((fiat_p251_uint256)(arg1[1]) * ((fiat_p251_uint128)(uint64_t)(arg2[3]) * ((uint64_t)0x2 * 0x9)));
  fiat_p251_uint256 x7 = ((uint64_t)(arg1[3]) * (fiat_p251_uint256)(arg2[0]));
  fiat_p251_uint256 x8 = ((fiat_p251_uint256)(arg1[2]) * (arg2[1]));
  fiat_p251_uint256 x9 = ((fiat_p251_uint256)(arg1[2]) * (arg2[0]));
  fiat_p251_uint256 x10 = ((fiat_p251_uint256)(arg1[1]) * (arg2[2]));
  fiat_p251_uint256 x11 = ((fiat_p251_uint256)(arg1[1]) * (arg2[1]));
  fiat_p251_uint256 x12 = ((fiat_p251_uint256)(arg1[1]) * (arg2[0]));
  fiat_p251_uint256 x13 = ((fiat_p251_uint256)(arg1[0]) * (uint64_t)(arg2[3]));
  fiat_p251_uint256 x14 = ((fiat_p251_uint256)(arg1[0]) * (arg2[2]));
  fiat_p251_uint256 x15 = ((fiat_p251_uint256)(arg1[0]) * (arg2[1]));
  fiat_p251_uint256 x16 = ((fiat_p251_uint256)(arg1[0]) * (arg2[0]));
  fiat_p251_uint256 x17 = (x16 + (x6 + (x5 + x3)));
  fiat_p251_uint128 x18 = (fiat_p251_uint128)(x17 >> 63);
  uint64_t x19 = (uint64_t)(x17 & UINT64_C(0x7fffffffffffffff));
  fiat_p251_uint256 x20 = (x13 + (x10 + (x8 + x7)));
  fiat_p251_uint256 x21 = (x14 + (x11 + (x9 + x1)));
  fiat_p251_uint256 x22 = (x15 + (x12 + (x4 + x2)));
  fiat_p251_uint256 x23 = (x18 + x22);
  fiat_p251_uint128 x24 = (fiat_p251_uint128)(x23 >> 63);
  uint64_t x25 = (uint64_t)(x23 & UINT64_C(0x7fffffffffffffff));
  fiat_p251_uint256 x26 = (x24 + x21);
  fiat_p251_uint128 x27 = (fiat_p251_uint128)(x26 >> 63);
  uint64_t x28 = (uint64_t)(x26 & UINT64_C(0x7fffffffffffffff));
  fiat_p251_uint256 x29 = (x27 + x20);
  fiat_p251_uint128 x30 = (fiat_p251_uint128)(x29 >> 62);
  uint64_t x31 = (uint64_t)(x29 & UINT64_C(0x3fffffffffffffff));
  fiat_p251_uint128 x32 = (x30 * (fiat_p251_uint128)0x9);
  fiat_p251_uint128 x33 = (x19 + x32);
  uint64_t x34 = (uint64_t)(x33 >> 63);
  uint64_t x35 = (uint64_t)(x33 & UINT64_C(0x7fffffffffffffff));
  uint64_t x36 = (x34 + x25);
  uint64_t x37 = (x36 >> 63);
  uint64_t x38 = (x36 & UINT64_C(0x7fffffffffffffff));
  uint64_t x39 = (x37 + x28);
  out1[0] = x35;
  out1[1] = x38;
  out1[2] = x39;
  out1[3] = x31;
}

/*
 * The function fiat_p251_carry_square squares a field element and reduces the result.
 * Postconditions:
 *   eval out1 mod m = (eval arg1 * eval arg1) mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0xd333333333333332]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x4666666666666666]]
 */
static void fiat_p251_carry_square(uint64_t out1[4], const fiat_p251_uint128 arg1[4]) {
  fiat_p251_uint128 x1 = ((uint64_t)(arg1[3]) * (fiat_p251_uint128)0x9);
  fiat_p251_uint128 x2 = (x1 * (fiat_p251_uint128)0x2);
  fiat_p251_uint128 x3 = ((uint64_t)(arg1[3]) * (fiat_p251_uint128)0x2);
  fiat_p251_uint128 x4 = ((arg1[2]) * (fiat_p251_uint128)0x9);
  fiat_p251_uint128 x5 = ((arg1[2]) * (fiat_p251_uint128)0x2);
  fiat_p251_uint128 x6 = ((arg1[1]) * (fiat_p251_uint128)0x2);
  fiat_p251_uint256 x7 = ((uint64_t)(arg1[3]) * (fiat_p251_uint256)(x1 * (fiat_p251_uint128)0x2));
  fiat_p251_uint256 x8 = ((fiat_p251_uint256)(arg1[2]) * (x2 * (fiat_p251_uint128)0x2));
  fiat_p251_uint256 x9 = ((fiat_p251_uint256)(arg1[2]) * (x4 * (fiat_p251_uint128)0x2));
  fiat_p251_uint256 x10 = ((fiat_p251_uint256)(arg1[1]) * (x2 * (fiat_p251_uint128)0x2));
  fiat_p251_uint256 x11 = ((fiat_p251_uint256)(arg1[1]) * x5);
  fiat_p251_uint256 x12 = ((fiat_p251_uint256)(arg1[1]) * (arg1[1]));
  fiat_p251_uint256 x13 = ((fiat_p251_uint256)(arg1[0]) * x3);
  fiat_p251_uint256 x14 = ((fiat_p251_uint256)(arg1[0]) * x5);
  fiat_p251_uint256 x15 = ((fiat_p251_uint256)(arg1[0]) * x6);
  fiat_p251_uint256 x16 = ((fiat_p251_uint256)(arg1[0]) * (arg1[0]));
  fiat_p251_uint256 x17 = (x16 + (x10 + x9));
  fiat_p251_uint128 x18 = (fiat_p251_uint128)(x17 >> 63);
  uint64_t x19 = (uint64_t)(x17 & UINT64_C(0x7fffffffffffffff));
  fiat_p251_uint256 x20 = (x13 + x11);
  fiat_p251_uint256 x21 = (x14 + (x12 + x7));
  fiat_p251_uint256 x22 = (x15 + x8);
  fiat_p251_uint256 x23 = (x18 + x22);
  fiat_p251_uint128 x24 = (fiat_p251_uint128)(x23 >> 63);
  uint64_t x25 = (uint64_t)(x23 & UINT64_C(0x7fffffffffffffff));
  fiat_p251_uint256 x26 = (x24 + x21);
  fiat_p251_uint128 x27 = (fiat_p251_uint128)(x26 >> 63);
  uint64_t x28 = (uint64_t)(x26 & UINT64_C(0x7fffffffffffffff));
  fiat_p251_uint256 x29 = (x27 + x20);
  fiat_p251_uint128 x30 = (fiat_p251_uint128)(x29 >> 62);
  uint64_t x31 = (uint64_t)(x29 & UINT64_C(0x3fffffffffffffff));
  fiat_p251_uint128 x32 = (x30 * (fiat_p251_uint128)0x9);
  fiat_p251_uint128 x33 = (x19 + x32);
  uint64_t x34 = (uint64_t)(x33 >> 63);
  uint64_t x35 = (uint64_t)(x33 & UINT64_C(0x7fffffffffffffff));
  uint64_t x36 = (x34 + x25);
  uint64_t x37 = (x36 >> 63);
  uint64_t x38 = (x36 & UINT64_C(0x7fffffffffffffff));
  uint64_t x39 = (x37 + x28);
  out1[0] = x35;
  out1[1] = x38;
  out1[2] = x39;
  out1[3] = x31;
}

/*
 * The function fiat_p251_carry reduces a field element.
 * Postconditions:
 *   eval out1 mod m = eval arg1 mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0xd333333333333332]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x4666666666666666]]
 */
static void fiat_p251_carry(uint64_t out1[4], const fiat_p251_uint128 arg1[4]) {
  fiat_p251_uint128 x1 = (arg1[0]);
  fiat_p251_uint128 x2 = ((uint64_t)(x1 >> 63) + (arg1[1]));
  fiat_p251_uint128 x3 = ((uint64_t)(x2 >> 63) + (arg1[2]));
  uint64_t x4 = ((uint64_t)(x3 >> 63) + (uint64_t)(arg1[3]));
  uint64_t x5 = ((uint64_t)(x1 & UINT64_C(0x7fffffffffffffff)) + ((x4 >> 62) * (uint64_t)0x9));
  uint64_t x6 = ((x5 >> 63) + (uint64_t)(x2 & UINT64_C(0x7fffffffffffffff)));
  uint64_t x7 = (x5 & UINT64_C(0x7fffffffffffffff));
  uint64_t x8 = (x6 & UINT64_C(0x7fffffffffffffff));
  uint64_t x9 = ((x6 >> 63) + (uint64_t)(x3 & UINT64_C(0x7fffffffffffffff)));
  uint64_t x10 = (x4 & UINT64_C(0x3fffffffffffffff));
  out1[0] = x7;
  out1[1] = x8;
  out1[2] = x9;
  out1[3] = x10;
}

/*
 * The function fiat_p251_add adds two field elements.
 * Postconditions:
 *   eval out1 mod m = (eval arg1 + eval arg2) mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x4666666666666666]]
 *   arg2: [[0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x4666666666666666]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0xd333333333333332]]
 */
static void fiat_p251_add(fiat_p251_uint128 out1[4], const uint64_t arg1[4], const uint64_t arg2[4]) {
  fiat_p251_uint128 x1 = ((fiat_p251_uint128)(arg1[0]) + (arg2[0]));
  fiat_p251_uint128 x2 = ((fiat_p251_uint128)(arg1[1]) + (arg2[1]));
  fiat_p251_uint128 x3 = ((fiat_p251_uint128)(arg1[2]) + (arg2[2]));
  uint64_t x4 = ((arg1[3]) + (arg2[3]));
  out1[0] = x1;
  out1[1] = x2;
  out1[2] = x3;
  out1[3] = x4;
}

/*
 * The function fiat_p251_sub subtracts two field elements.
 * Postconditions:
 *   eval out1 mod m = (eval arg1 - eval arg2) mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x4666666666666666]]
 *   arg2: [[0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x4666666666666666]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0xd333333333333332]]
 */
static void fiat_p251_sub(fiat_p251_uint128 out1[4], const uint64_t arg1[4], const uint64_t arg2[4]) {
  fiat_p251_uint128 x1 = (((fiat_p251_uint128)UINT64_C(0xffffffffffffffee) + (arg1[0])) - (arg2[0]));
  fiat_p251_uint128 x2 = (((fiat_p251_uint128)UINT64_C(0xfffffffffffffffe) + (arg1[1])) - (arg2[1]));
  fiat_p251_uint128 x3 = (((fiat_p251_uint128)UINT64_C(0xfffffffffffffffe) + (arg1[2])) - (arg2[2]));
  uint64_t x4 = ((UINT64_C(0x7ffffffffffffffe) + (arg1[3])) - (arg2[3]));
  out1[0] = x1;
  out1[1] = x2;
  out1[2] = x3;
  out1[3] = x4;
}

/*
 * The function fiat_p251_opp negates a field element.
 * Postconditions:
 *   eval out1 mod m = -eval arg1 mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x4666666666666666]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0x1a666666666666664], [0x0 ~> 0xd333333333333332]]
 */
static void fiat_p251_opp(uint64_t out1[4], const uint64_t arg1[4]) {
  uint64_t x1 = (UINT64_C(0xffffffffffffffee) - (arg1[0]));
  uint64_t x2 = (UINT64_C(0xfffffffffffffffe) - (arg1[1]));
  uint64_t x3 = (UINT64_C(0xfffffffffffffffe) - (arg1[2]));
  uint64_t x4 = (UINT64_C(0x7ffffffffffffffe) - (arg1[3]));
  out1[0] = x1;
  out1[1] = x2;
  out1[2] = x3;
  out1[3] = x4;
}

/*
 * The function fiat_p251_selectznz is a multi-limb conditional select.
 * Postconditions:
 *   eval out1 = (if arg1 = 0 then eval arg2 else eval arg3)
 *
 * Input Bounds:
 *   arg1: [0x0 ~> 0x1]
 *   arg2: [[0x0 ~> 0xffffffffffffffff], [0x0 ~> 0xffffffffffffffff], [0x0 ~> 0xffffffffffffffff], [0x0 ~> 0xffffffffffffffff]]
 *   arg3: [[0x0 ~> 0xffffffffffffffff], [0x0 ~> 0xffffffffffffffff], [0x0 ~> 0xffffffffffffffff], [0x0 ~> 0xffffffffffffffff]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0xffffffffffffffff], [0x0 ~> 0xffffffffffffffff], [0x0 ~> 0xffffffffffffffff], [0x0 ~> 0xffffffffffffffff]]
 */
static void fiat_p251_selectznz(uint64_t out1[4], fiat_p251_uint1 arg1, const uint64_t arg2[4], const uint64_t arg3[4]) {
  uint64_t x1;
  fiat_p251_cmovznz_u64(&x1, arg1, (arg2[0]), (arg3[0]));
  uint64_t x2;
  fiat_p251_cmovznz_u64(&x2, arg1, (arg2[1]), (arg3[1]));
  uint64_t x3;
  fiat_p251_cmovznz_u64(&x3, arg1, (arg2[2]), (arg3[2]));
  uint64_t x4;
  fiat_p251_cmovznz_u64(&x4, arg1, (arg2[3]), (arg3[3]));
  out1[0] = x1;
  out1[1] = x2;
  out1[2] = x3;
  out1[3] = x4;
}

/*
 * The function fiat_p251_to_bytes serializes a field element to bytes in little-endian order.
 * Postconditions:
 *   out1 = map (λ x, ⌊((eval arg1 mod m) mod 2^(8 * (x + 1))) / 2^(8 * x)⌋) [0..31]
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x4666666666666666]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0x7]]
 */
static void fiat_p251_to_bytes(uint8_t out1[32], const uint64_t arg1[4]) {
  uint64_t x1;
  fiat_p251_uint1 x2;
  fiat_p251_subborrowx_u63(&x1, &x2, 0x0, (arg1[0]), UINT64_C(0x7ffffffffffffff7));
  uint64_t x3;
  fiat_p251_uint1 x4;
  fiat_p251_subborrowx_u63(&x3, &x4, x2, (arg1[1]), UINT64_C(0x7fffffffffffffff));
  uint64_t x5;
  fiat_p251_uint1 x6;
  fiat_p251_subborrowx_u63(&x5, &x6, x4, (arg1[2]), UINT64_C(0x7fffffffffffffff));
  uint64_t x7;
  fiat_p251_uint1 x8;
  fiat_p251_subborrowx_u62(&x7, &x8, x6, (arg1[3]), UINT64_C(0x3fffffffffffffff));
  uint64_t x9;
  fiat_p251_cmovznz_u64(&x9, x8, 0x0, UINT64_C(0xffffffffffffffff));
  uint64_t x10;
  fiat_p251_uint1 x11;
  fiat_p251_addcarryx_u63(&x10, &x11, 0x0, x1, (x9 & UINT64_C(0x7ffffffffffffff7)));
  uint64_t x12;
  fiat_p251_uint1 x13;
  fiat_p251_addcarryx_u63(&x12, &x13, x11, x3, (x9 & UINT64_C(0x7fffffffffffffff)));
  uint64_t x14;
  fiat_p251_uint1 x15;
  fiat_p251_addcarryx_u63(&x14, &x15, x13, x5, (x9 & UINT64_C(0x7fffffffffffffff)));
  uint64_t x16;
  fiat_p251_uint1 x17;
  fiat_p251_addcarryx_u62(&x16, &x17, x15, x7, (x9 & UINT64_C(0x3fffffffffffffff)));
  fiat_p251_uint128 x18 = ((fiat_p251_uint128)x16 << 5);
  fiat_p251_uint128 x19 = ((fiat_p251_uint128)x14 << 6);
  fiat_p251_uint128 x20 = ((fiat_p251_uint128)x12 << 7);
  uint64_t x21 = (x10 >> 8);
  uint8_t x22 = (uint8_t)(x10 & UINT8_C(0xff));
  uint64_t x23 = (x21 >> 8);
  uint8_t x24 = (uint8_t)(x21 & UINT8_C(0xff));
  uint64_t x25 = (x23 >> 8);
  uint8_t x26 = (uint8_t)(x23 & UINT8_C(0xff));
  uint64_t x27 = (x25 >> 8);
  uint8_t x28 = (uint8_t)(x25 & UINT8_C(0xff));
  uint64_t x29 = (x27 >> 8);
  uint8_t x30 = (uint8_t)(x27 & UINT8_C(0xff));
  uint64_t x31 = (x29 >> 8);
  uint8_t x32 = (uint8_t)(x29 & UINT8_C(0xff));
  uint8_t x33 = (uint8_t)(x31 >> 8);
  uint8_t x34 = (uint8_t)(x31 & UINT8_C(0xff));
  fiat_p251_uint128 x35 = (x33 + x20);
  uint64_t x36 = (uint64_t)(x35 >> 8);
  uint8_t x37 = (uint8_t)(x35 & UINT8_C(0xff));
  uint64_t x38 = (x36 >> 8);
  uint8_t x39 = (uint8_t)(x36 & UINT8_C(0xff));
  uint64_t x40 = (x38 >> 8);
  uint8_t x41 = (uint8_t)(x38 & UINT8_C(0xff));
  uint64_t x42 = (x40 >> 8);
  uint8_t x43 = (uint8_t)(x40 & UINT8_C(0xff));
  uint64_t x44 = (x42 >> 8);
  uint8_t x45 = (uint8_t)(x42 & UINT8_C(0xff));
  uint64_t x46 = (x44 >> 8);
  uint8_t x47 = (uint8_t)(x44 & UINT8_C(0xff));
  uint64_t x48 = (x46 >> 8);
  uint8_t x49 = (uint8_t)(x46 & UINT8_C(0xff));
  uint8_t x50 = (uint8_t)(x48 >> 8);
  uint8_t x51 = (uint8_t)(x48 & UINT8_C(0xff));
  fiat_p251_uint128 x52 = (x50 + x19);
  uint64_t x53 = (uint64_t)(x52 >> 8);
  uint8_t x54 = (uint8_t)(x52 & UINT8_C(0xff));
  uint64_t x55 = (x53 >> 8);
  uint8_t x56 = (uint8_t)(x53 & UINT8_C(0xff));
  uint64_t x57 = (x55 >> 8);
  uint8_t x58 = (uint8_t)(x55 & UINT8_C(0xff));
  uint64_t x59 = (x57 >> 8);
  uint8_t x60 = (uint8_t)(x57 & UINT8_C(0xff));
  uint64_t x61 = (x59 >> 8);
  uint8_t x62 = (uint8_t)(x59 & UINT8_C(0xff));
  uint64_t x63 = (x61 >> 8);
  uint8_t x64 = (uint8_t)(x61 & UINT8_C(0xff));
  uint64_t x65 = (x63 >> 8);
  uint8_t x66 = (uint8_t)(x63 & UINT8_C(0xff));
  uint8_t x67 = (uint8_t)(x65 >> 8);
  uint8_t x68 = (uint8_t)(x65 & UINT8_C(0xff));
  fiat_p251_uint128 x69 = (x67 + x18);
  uint64_t x70 = (uint64_t)(x69 >> 8);
  uint8_t x71 = (uint8_t)(x69 & UINT8_C(0xff));
  uint64_t x72 = (x70 >> 8);
  uint8_t x73 = (uint8_t)(x70 & UINT8_C(0xff));
  uint64_t x74 = (x72 >> 8);
  uint8_t x75 = (uint8_t)(x72 & UINT8_C(0xff));
  uint64_t x76 = (x74 >> 8);
  uint8_t x77 = (uint8_t)(x74 & UINT8_C(0xff));
  uint64_t x78 = (x76 >> 8);
  uint8_t x79 = (uint8_t)(x76 & UINT8_C(0xff));
  uint64_t x80 = (x78 >> 8);
  uint8_t x81 = (uint8_t)(x78 & UINT8_C(0xff));
  uint64_t x82 = (x80 >> 8);
  uint8_t x83 = (uint8_t)(x80 & UINT8_C(0xff));
  uint8_t x84 = (uint8_t)(x82 >> 8);
  uint8_t x85 = (uint8_t)(x82 & UINT8_C(0xff));
  out1[0] = x22;
  out1[1] = x24;
  out1[2] = x26;
  out1[3] = x28;
  out1[4] = x30;
  out1[5] = x32;
  out1[6] = x34;
  out1[7] = x37;
  out1[8] = x39;
  out1[9] = x41;
  out1[10] = x43;
  out1[11] = x45;
  out1[12] = x47;
  out1[13] = x49;
  out1[14] = x51;
  out1[15] = x54;
  out1[16] = x56;
  out1[17] = x58;
  out1[18] = x60;
  out1[19] = x62;
  out1[20] = x64;
  out1[21] = x66;
  out1[22] = x68;
  out1[23] = x71;
  out1[24] = x73;
  out1[25] = x75;
  out1[26] = x77;
  out1[27] = x79;
  out1[28] = x81;
  out1[29] = x83;
  out1[30] = x85;
  out1[31] = x84;
}

/*
 * The function fiat_p251_from_bytes deserializes a field element from bytes in little-endian order.
 * Postconditions:
 *   eval out1 mod m = bytes_eval arg1 mod m
 *
 * Input Bounds:
 *   arg1: [[0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0xff], [0x0 ~> 0x7]]
 * Output Bounds:
 *   out1: [[0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x8ccccccccccccccc], [0x0 ~> 0x4666666666666666]]
 */
static void fiat_p251_from_bytes(uint64_t out1[4], const uint8_t arg1[32]) {
  uint64_t x1 = ((uint64_t)(arg1[31]) << 59);
  uint64_t x2 = ((uint64_t)(arg1[30]) << 51);
  uint64_t x3 = ((uint64_t)(arg1[29]) << 43);
  uint64_t x4 = ((uint64_t)(arg1[28]) << 35);
  uint64_t x5 = ((uint64_t)(arg1[27]) << 27);
  uint64_t x6 = ((uint64_t)(arg1[26]) << 19);
  uint64_t x7 = ((uint64_t)(arg1[25]) << 11);
  uint64_t x8 = ((uint64_t)(arg1[24]) << 3);
  fiat_p251_uint128 x9 = ((fiat_p251_uint128)(arg1[23]) << 58);
  uint64_t x10 = ((uint64_t)(arg1[22]) << 50);
  uint64_t x11 = ((uint64_t)(arg1[21]) << 42);
  uint64_t x12 = ((uint64_t)(arg1[20]) << 34);
  uint64_t x13 = ((uint64_t)(arg1[19]) << 26);
  uint64_t x14 = ((uint64_t)(arg1[18]) << 18);
  uint64_t x15 = ((uint64_t)(arg1[17]) << 10);
  uint64_t x16 = ((uint64_t)(arg1[16]) << 2);
  fiat_p251_uint128 x17 = ((fiat_p251_uint128)(arg1[15]) << 57);
  uint64_t x18 = ((uint64_t)(arg1[14]) << 49);
  uint64_t x19 = ((uint64_t)(arg1[13]) << 41);
  uint64_t x20 = ((uint64_t)(arg1[12]) << 33);
  uint64_t x21 = ((uint64_t)(arg1[11]) << 25);
  uint64_t x22 = ((uint64_t)(arg1[10]) << 17);
  uint64_t x23 = ((uint64_t)(arg1[9]) << 9);
  uint64_t x24 = ((uint64_t)(arg1[8]) * 0x2);
  uint64_t x25 = ((uint64_t)(arg1[7]) << 56);
  uint64_t x26 = ((uint64_t)(arg1[6]) << 48);
  uint64_t x27 = ((uint64_t)(arg1[5]) << 40);
  uint64_t x28 = ((uint64_t)(arg1[4]) << 32);
  uint64_t x29 = ((uint64_t)(arg1[3]) << 24);
  uint64_t x30 = ((uint64_t)(arg1[2]) << 16);
  uint64_t x31 = ((uint64_t)(arg1[1]) << 8);
  uint8_t x32 = (arg1[0]);
  uint64_t x33 = (x32 + (x31 + (x30 + (x29 + (x28 + (x27 + (x26 + x25)))))));
  fiat_p251_uint1 x34 = (fiat_p251_uint1)(x33 >> 63);
  uint64_t x35 = (x33 & UINT64_C(0x7fffffffffffffff));
  uint64_t x36 = (x8 + (x7 + (x6 + (x5 + (x4 + (x3 + (x2 + x1)))))));
  fiat_p251_uint128 x37 = (x16 + (x15 + (x14 + (x13 + (x12 + (x11 + (x10 + x9)))))));
  fiat_p251_uint128 x38 = (x24 + (x23 + (x22 + (x21 + (x20 + (x19 + (x18 + x17)))))));
  fiat_p251_uint128 x39 = (x34 + x38);
  uint8_t x40 = (uint8_t)(x39 >> 63);
  uint64_t x41 = (uint64_t)(x39 & UINT64_C(0x7fffffffffffffff));
  fiat_p251_uint128 x42 = (x40 + x37);
  uint8_t x43 = (uint8_t)(x42 >> 63);
  uint64_t x44 = (uint64_t)(x42 & UINT64_C(0x7fffffffffffffff));
  uint64_t x45 = (x43 + x36);
  out1[0] = x35;
  out1[1] = x41;
  out1[2] = x44;
  out1[3] = x45;
}

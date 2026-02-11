// Spooky Hash
// A 128-bit noncryptographic hash, for checksums and table lookup
// By Bob Jenkins.  Public domain.
//   Oct 31 2010: published framework, disclaimer ShortHash isn't right
//   Nov 7 2010: disabled ShortHash
//   Oct 31 2011: replace End, ShortMix, ShortEnd, enable ShortHash again

#include <memory.h>
#include "spookyhash.h"

#define ALLOW_UNALIGNED_READS 1

// @brief spookyhash_short is used for messages under 192 bytes in length
//
// spookyhash_short has a low startup cost, the normal mode is good for long
// keys, the cost crossover is at about 192 bytes.  The two modes were
// held to the same quality bar.
//
// @param message (array of bytes, not necessarily aligned)
// @param length length of message in bytes
// @param hash1 in/out: in seed 1, out hash value 1
// @param hash2 in/out: in seed 2, out hash value 2
static void spookyhash_short(const void* message, size_t length,
                             uint64_t* hash1, uint64_t* hash2) {

  // short hash ... it could be used on any message,
  // but it's used by spooky just for short messages.

  uint64_t buf[24]; // (sc_num_vars * 2)

  union {

    const uint8_t* p8;
    uint32_t* p32;
    uint64_t* p64;
    size_t i;
  } u;

  u.p8 = (const uint8_t*)message;

  if (!ALLOW_UNALIGNED_READS && (u.i & 0x7)) {

    memcpy(buf, message, length);
    u.p64 = buf;
  }

  size_t remainder = length % 32;
  uint64_t a = *hash1;
  uint64_t b = *hash2;
  uint64_t c = sc_const;
  uint64_t d = sc_const;

  if (length > 15) {

    const uint64_t* end = u.p64 + (length / 32) * 4;

    // handle all complete sets of 32 bytes
    for (; u.p64 < end; u.p64 += 4) {

      c += u.p64[0];
      d += u.p64[1];
      spookyhash_smix(&a, &b, &c, &d);
      a += u.p64[2];
      b += u.p64[3];
    }

    // Handle the case of 16+ remaining bytes.
    if (remainder >= 16) {

      c += u.p64[0];
      d += u.p64[1];
      spookyhash_smix(&a, &b, &c, &d);
      u.p64 += 2;
      remainder -= 16;
    }
  }

  // handle the last 0..15 bytes, and its length
  d += ((uint64_t)length) << 56;

  switch (remainder) {

    case 15:
      d += ((uint64_t)u.p8[14]) << 48;
      // fall-through
    case 14:
      d += ((uint64_t)u.p8[13]) << 40;
      // fall-through
    case 13:
      d += ((uint64_t)u.p8[12]) << 32;
      // fall-through
    case 12:
      d += u.p32[2];
      c += u.p64[0];
      break;
    case 11:
      d += ((uint64_t)u.p8[10]) << 16;
      // fall-through
    case 10:
      d += ((uint64_t)u.p8[9]) << 8;
      // fall-through
    case 9:
      d += (uint64_t)u.p8[8];
      // fall-through
    case 8:
      c += u.p64[0];
      break;
    case 7:
      c += ((uint64_t)u.p8[6]) << 48;
      // fall-through
    case 6:
      c += ((uint64_t)u.p8[5]) << 40;
      // fall-through
    case 5:
      c += ((uint64_t)u.p8[4]) << 32;
      // fall-through
    case 4:
      c += u.p32[0];
      break;
    case 3:
      c += ((uint64_t)u.p8[2]) << 16;
      // fall-through
    case 2:
      c += ((uint64_t)u.p8[1]) << 8;
      // fall-through
    case 1:
      c += (uint64_t)u.p8[0];
      break;
    case 0:
      c += sc_const;
      d += sc_const;
  }

  spookyhash_short_end(&a, &b, &c, &d);

  *hash1 = a;
  *hash2 = b;
}

// do the whole hash in one call
void spookyhash128(const void* message, size_t length, uint64_t* hash1,
                   uint64_t* hash2) {

  if (length < sc_buf_size) {

    spookyhash_short(message, length, hash1, hash2);
    return;
  }

  uint64_t h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11;
  uint64_t buf[12]; // sc_num_vars
  size_t remainder;
  uint64_t* end;

  union {

    const uint8_t* p8;
    uint64_t* p64;
    size_t i;
  } u;

  h0 = h3 = h6 = h9 = *hash1;
  h1 = h4 = h7 = h10 = *hash2;
  h2 = h5 = h8 = h11 = sc_const;

  u.p8 = (const uint8_t*)message;
  end = u.p64 + (length / sc_block_size) * sc_num_vars;

  // handle all whole sc_block_size blocks of bytes
  if (ALLOW_UNALIGNED_READS || ((u.i & 0x7) == 0)) {

    while (u.p64 < end) {

      spookyhash_mix(u.p64, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9,
                     &h10, &h11);
      u.p64 += sc_num_vars;
    }
  } else {

    while (u.p64 < end) {

      memcpy(buf, u.p64, sc_block_size);
      spookyhash_mix(buf, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9,
                     &h10, &h11);

      u.p64 += sc_num_vars;
    }
  }

  // handle the last partial block of sc_block_size bytes
  remainder = (length - ((const uint8_t*)end - (const uint8_t*)message));

  memcpy(buf, end, remainder);
  memset(((uint8_t*)buf) + remainder, 0, sc_block_size - remainder);

  ((uint8_t*)buf)[sc_block_size - 1] = remainder;

  // do some final mixing
  spookyhash_end(buf, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);
  *hash1 = h0;
  *hash2 = h1;
}

// init spooky state
void spookyhash_init(uint64_t seed1, uint64_t seed2, spooky_state* state) {

  state->length = 0;
  state->remainder = 0;
  state->vars[0] = seed1;
  state->vars[1] = seed2;
}

// add a message fragment to the state
void spookyhash_update(const void* message, size_t length,
                       spooky_state* state) {

  uint64_t h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11;
  size_t new_length = length + state->remainder;
  uint8_t remainder;
  const uint64_t* end;

  union {

    const uint8_t* p8;
    uint64_t* p64;
    size_t i;
  } u;

  // is this message fragment too short, if it is, stuff it away
  if (new_length < sc_buf_size) {

    memcpy(&((uint8_t*)state->data)[state->remainder], message, length);

    state->length = length + state->length;
    state->remainder = (uint8_t)new_length;
    return;
  }

  // init the variables
  if (state->length < sc_buf_size) {

    h0 = h3 = h6 = h9 = state->vars[0];
    h1 = h4 = h7 = h10 = state->vars[1];
    h2 = h5 = h8 = h11 = sc_const;
  } else {

    h0 = state->vars[0];
    h1 = state->vars[1];
    h2 = state->vars[2];
    h3 = state->vars[3];
    h4 = state->vars[4];
    h5 = state->vars[5];
    h6 = state->vars[6];
    h7 = state->vars[7];
    h8 = state->vars[8];
    h9 = state->vars[9];
    h10 = state->vars[10];
    h11 = state->vars[11];
  }

  state->length = length + state->length;

  // if we've got anything stuffed away, use it now
  if (state->remainder) {

    uint8_t prefix = sc_buf_size - state->remainder;

    memcpy(&(((uint8_t*)state->data)[state->remainder]), message, prefix);

    u.p64 = state->data;

    spookyhash_mix(u.p64, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9,
                   &h10, &h11);
    spookyhash_mix(&u.p64[sc_num_vars], &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7,
                   &h8, &h9, &h10, &h11);

    u.p8 = ((const uint8_t*)message) + prefix;

    length -= prefix;
  } else {

    u.p8 = (const uint8_t*)message;
  }

  // handle all whole blocks of sc_block_size bytes
  end = u.p64 + (length / sc_block_size) * sc_num_vars;
  remainder = (uint8_t)(length - ((const uint8_t*)end - u.p8));

  if (ALLOW_UNALIGNED_READS || (u.i & 0x7) == 0) {

    while (u.p64 < end) {

      spookyhash_mix(u.p64, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9,
                     &h10, &h11);
      u.p64 += sc_num_vars;
    }
  } else {

    while (u.p64 < end) {

      memcpy(state->data, u.p8, sc_block_size);
      spookyhash_mix(state->data, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8,
                     &h9, &h10, &h11);
      u.p64 += sc_num_vars;
    }
  }

  // stuff away the last few bytes
  state->remainder = remainder;
  memcpy(state->data, end, remainder);

  // stuff away the variables
  state->vars[0] = h0;
  state->vars[1] = h1;
  state->vars[2] = h2;
  state->vars[3] = h3;
  state->vars[4] = h4;
  state->vars[5] = h5;
  state->vars[6] = h6;
  state->vars[7] = h7;
  state->vars[8] = h8;
  state->vars[9] = h9;
  state->vars[10] = h10;
  state->vars[11] = h11;
}

// report the hash for the concatenation of all message fragments so far
void spookyhash_final(uint64_t* hash1, uint64_t* hash2, spooky_state* state) {

  // init the variables
  if (state->length < sc_buf_size) {

    *hash1 = state->vars[0];
    *hash2 = state->vars[1];

    spookyhash_short(state->data, state->length, hash1, hash2);
    return;
  }

  const uint64_t* data = (const uint64_t*)state->data;
  uint8_t remainder = state->remainder;

  uint64_t h0 = state->vars[0];
  uint64_t h1 = state->vars[1];
  uint64_t h2 = state->vars[2];
  uint64_t h3 = state->vars[3];
  uint64_t h4 = state->vars[4];
  uint64_t h5 = state->vars[5];
  uint64_t h6 = state->vars[6];
  uint64_t h7 = state->vars[7];
  uint64_t h8 = state->vars[8];
  uint64_t h9 = state->vars[9];
  uint64_t h10 = state->vars[10];
  uint64_t h11 = state->vars[11];

  if (remainder >= sc_block_size) {

    // state->data can contain two blocks; handle any whole first block
    spookyhash_mix(data, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10,
                   &h11);
    data += sc_num_vars;
    remainder -= sc_block_size;
  }

  // mix in the last partial block, and the length mod sc_block_size
  memset(&((uint8_t*)data)[remainder], 0, sc_block_size - remainder);

  ((uint8_t*)data)[sc_block_size - 1] = remainder;

  // do some final mixing
  spookyhash_end(data, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);

  *hash1 = h0;
  *hash2 = h1;
}


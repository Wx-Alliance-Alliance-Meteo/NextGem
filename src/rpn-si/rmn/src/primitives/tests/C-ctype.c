/* Copyright (C) 1995-1999, 2000, 2001 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Contributed by Ulrich Drepper <drepper@cygnus.com>, 1995.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */


/* This table's entries are taken from POSIX.2 Table 2-6
   ``LC_CTYPE Category Definition in the POSIX Locale''.

   The `_nl_C_LC_CTYPE_width' array is a GNU extension.

   In the `_nl_C_LC_CTYPE_class' array the value for EOF (== -1)
   is set to always return 0 and the conversion arrays return EOF.  */

static const char _nl_C_LC_CTYPE_class[768] =
  /* 0x80 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x86 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x8c */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x92 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x98 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x9e */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xa4 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xaa */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xb0 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xb6 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xbc */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xc2 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xc8 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xce */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xd4 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xda */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xe0 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xe6 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xec */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xf2 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xf8 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xfe */ "\000\000" "\000\000" "\002\000" "\002\000" "\002\000" "\002\000"
  /* 0x04 */ "\002\000" "\002\000" "\002\000" "\002\000" "\002\000" "\003\040"
  /* 0x0a */ "\002\040" "\002\040" "\002\040" "\002\040" "\002\000" "\002\000"
  /* 0x10 */ "\002\000" "\002\000" "\002\000" "\002\000" "\002\000" "\002\000"
  /* 0x16 */ "\002\000" "\002\000" "\002\000" "\002\000" "\002\000" "\002\000"
  /* 0x1c */ "\002\000" "\002\000" "\002\000" "\002\000" "\001\140" "\004\300"
  /* 0x22 */ "\004\300" "\004\300" "\004\300" "\004\300" "\004\300" "\004\300"
  /* 0x28 */ "\004\300" "\004\300" "\004\300" "\004\300" "\004\300" "\004\300"
  /* 0x2e */ "\004\300" "\004\300" "\010\330" "\010\330" "\010\330" "\010\330"
  /* 0x34 */ "\010\330" "\010\330" "\010\330" "\010\330" "\010\330" "\010\330"
  /* 0x3a */ "\004\300" "\004\300" "\004\300" "\004\300" "\004\300" "\004\300"
  /* 0x40 */ "\004\300" "\010\325" "\010\325" "\010\325" "\010\325" "\010\325"
  /* 0x46 */ "\010\325" "\010\305" "\010\305" "\010\305" "\010\305" "\010\305"
  /* 0x4c */ "\010\305" "\010\305" "\010\305" "\010\305" "\010\305" "\010\305"
  /* 0x52 */ "\010\305" "\010\305" "\010\305" "\010\305" "\010\305" "\010\305"
  /* 0x58 */ "\010\305" "\010\305" "\010\305" "\004\300" "\004\300" "\004\300"
  /* 0x5e */ "\004\300" "\004\300" "\004\300" "\010\326" "\010\326" "\010\326"
  /* 0x64 */ "\010\326" "\010\326" "\010\326" "\010\306" "\010\306" "\010\306"
  /* 0x6a */ "\010\306" "\010\306" "\010\306" "\010\306" "\010\306" "\010\306"
  /* 0x70 */ "\010\306" "\010\306" "\010\306" "\010\306" "\010\306" "\010\306"
  /* 0x76 */ "\010\306" "\010\306" "\010\306" "\010\306" "\010\306" "\004\300"
  /* 0x7c */ "\004\300" "\004\300" "\004\300" "\002\000" "\000\000" "\000\000"
  /* 0x82 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x88 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x8e */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x94 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0x9a */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xa0 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xa6 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xac */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xb2 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xb8 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xbe */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xc4 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xca */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xd0 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xd6 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xdc */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xe2 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xe8 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xee */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xf4 */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
  /* 0xfa */ "\000\000" "\000\000" "\000\000" "\000\000" "\000\000" "\000\000"
;
static const int _nl_C_LC_CTYPE_toupper[384] = {
  /* 0x80 */ 0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
  /* 0x88 */ 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
  /* 0x90 */ 0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97,
  /* 0x98 */ 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
  /* 0xa0 */ 0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
  /* 0xa8 */ 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
  /* 0xb0 */ 0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7,
  /* 0xb8 */ 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf,
  /* 0xc0 */ 0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7,
  /* 0xc8 */ 0xc8, 0xc9, 0xca, 0xcb, 0xcc, 0xcd, 0xce, 0xcf,
  /* 0xd0 */ 0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7,
  /* 0xd8 */ 0xd8, 0xd9, 0xda, 0xdb, 0xdc, 0xdd, 0xde, 0xdf,
  /* 0xe0 */ 0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7,
  /* 0xe8 */ 0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef,
  /* 0xf0 */ 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7,
  /* 0xf8 */ 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xffffffff,
  /* 0x00 */ 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
  /* 0x08 */ 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
  /* 0x10 */ 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
  /* 0x18 */ 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
  /* 0x20 */ 0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27,
  /* 0x28 */ 0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
  /* 0x30 */ 0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
  /* 0x38 */ 0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f,
  /* 0x40 */ 0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47,
  /* 0x48 */ 0x48, 0x49, 0x4a, 0x4b, 0x4c, 0x4d, 0x4e, 0x4f,
  /* 0x50 */ 0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57,
  /* 0x58 */ 0x58, 0x59, 0x5a, 0x5b, 0x5c, 0x5d, 0x5e, 0x5f,
  /* 0x60 */ 0x60, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47,
  /* 0x68 */ 0x48, 0x49, 0x4a, 0x4b, 0x4c, 0x4d, 0x4e, 0x4f,
  /* 0x70 */ 0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57,
  /* 0x78 */ 0x58, 0x59, 0x5a, 0x7b, 0x7c, 0x7d, 0x7e, 0x7f,
  /* 0x80 */ 0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
  /* 0x88 */ 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
  /* 0x90 */ 0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97,
  /* 0x98 */ 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
  /* 0xa0 */ 0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
  /* 0xa8 */ 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
  /* 0xb0 */ 0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7,
  /* 0xb8 */ 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf,
  /* 0xc0 */ 0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7,
  /* 0xc8 */ 0xc8, 0xc9, 0xca, 0xcb, 0xcc, 0xcd, 0xce, 0xcf,
  /* 0xd0 */ 0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7,
  /* 0xd8 */ 0xd8, 0xd9, 0xda, 0xdb, 0xdc, 0xdd, 0xde, 0xdf,
  /* 0xe0 */ 0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7,
  /* 0xe8 */ 0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef,
  /* 0xf0 */ 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7,
  /* 0xf8 */ 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff
};
static const int _nl_C_LC_CTYPE_tolower[384] ={
  /* 0x80 */ 0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
  /* 0x88 */ 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
  /* 0x90 */ 0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97,
  /* 0x98 */ 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
  /* 0xa0 */ 0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
  /* 0xa8 */ 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
  /* 0xb0 */ 0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7,
  /* 0xb8 */ 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf,
  /* 0xc0 */ 0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7,
  /* 0xc8 */ 0xc8, 0xc9, 0xca, 0xcb, 0xcc, 0xcd, 0xce, 0xcf,
  /* 0xd0 */ 0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7,
  /* 0xd8 */ 0xd8, 0xd9, 0xda, 0xdb, 0xdc, 0xdd, 0xde, 0xdf,
  /* 0xe0 */ 0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7,
  /* 0xe8 */ 0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef,
  /* 0xf0 */ 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7,
  /* 0xf8 */ 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xffffffff,
  /* 0x00 */ 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
  /* 0x08 */ 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
  /* 0x10 */ 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
  /* 0x18 */ 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
  /* 0x20 */ 0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27,
  /* 0x28 */ 0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
  /* 0x30 */ 0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
  /* 0x38 */ 0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f,
  /* 0x40 */ 0x40, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67,
  /* 0x48 */ 0x68, 0x69, 0x6a, 0x6b, 0x6c, 0x6d, 0x6e, 0x6f,
  /* 0x50 */ 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77,
  /* 0x58 */ 0x78, 0x79, 0x7a, 0x5b, 0x5c, 0x5d, 0x5e, 0x5f,
  /* 0x60 */ 0x60, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67,
  /* 0x68 */ 0x68, 0x69, 0x6a, 0x6b, 0x6c, 0x6d, 0x6e, 0x6f,
  /* 0x70 */ 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77,
  /* 0x78 */ 0x78, 0x79, 0x7a, 0x7b, 0x7c, 0x7d, 0x7e, 0x7f,
  /* 0x80 */ 0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
  /* 0x88 */ 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
  /* 0x90 */ 0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97,
  /* 0x98 */ 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
  /* 0xa0 */ 0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
  /* 0xa8 */ 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
  /* 0xb0 */ 0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7,
  /* 0xb8 */ 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf,
  /* 0xc0 */ 0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7,
  /* 0xc8 */ 0xc8, 0xc9, 0xca, 0xcb, 0xcc, 0xcd, 0xce, 0xcf,
  /* 0xd0 */ 0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7,
  /* 0xd8 */ 0xd8, 0xd9, 0xda, 0xdb, 0xdc, 0xdd, 0xde, 0xdf,
  /* 0xe0 */ 0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7,
  /* 0xe8 */ 0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef,
  /* 0xf0 */ 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7,
  /* 0xf8 */ 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff
};


#define b(t,x,o) (((const t *) _nl_C_LC_CTYPE_##x) + o)

#pragma weak __ctype_b
const unsigned short int *__ctype_b = b (unsigned short int, class, 128);

#pragma weak __ctype_tolower
const int *__ctype_tolower = b (int, tolower, 128);

#pragma weak __ctype_toupper
const int *__ctype_toupper = b (int, toupper, 128);

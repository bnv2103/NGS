#pragma once

#define _FILE_OFFSET_BITS 64

#ifdef WIN32
#define ftell64(a)     _ftelli64(a)
#define fseek64(a,b,c) _fseeki64(a,b,c)
typedef __int64 off_type;
#else
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
#endif

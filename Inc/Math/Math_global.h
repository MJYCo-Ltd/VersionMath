#ifndef MATH_GLOBAL_H
#define MATH_GLOBAL_H

#include "../NoQt.h"

#if defined(MATH_LIBRARY)
#  define MATH_EXPORT Q_DECL_EXPORT
#else
#  define MATH_EXPORT Q_DECL_IMPORT
#endif

#endif // MATH_GLOBAL_H

#ifndef SAT_GLOBAL_H
#define SAT_GLOBAL_H

#include "../NoQt.h"

#if defined(ALGORITHM_LIBRARY)
#  define ALGORITHM_EXPORT Q_DECL_EXPORT
#else
#  define ALGORITHM_EXPORT Q_DECL_IMPORT
#endif

#endif // SAT_GLOBAL_H

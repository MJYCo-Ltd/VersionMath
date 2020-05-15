#ifndef SOFA_GLOBAL_H
#define SOFA_GLOBAL_H

#include "../NoQt.h"

#if defined(SOFA_LIBRARY)
#  define SOFASHARED_EXPORT Q_DECL_EXPORT
#else
#  define SOFASHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // SOFA_GLOBAL_H

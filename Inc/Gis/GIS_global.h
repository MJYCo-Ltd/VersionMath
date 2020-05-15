#ifndef GIS_GLOBAL_H
#define GIS_GLOBAL_H

#include "../NoQt.h"

#if defined(GIS_LIBRARY)
#  define GISSHARED_EXPORT Q_DECL_EXPORT
#else
#  define GISSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // GIS_GLOBAL_H

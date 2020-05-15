#ifndef MY_GISMATH_GLOBAL_H
#define MY_GISMATH_GLOBAL_H


#if defined(GIS_LIBRARY)
#  define GISMATHSHARED_EXPORT __declspec(dllexport)
#else
#  define GISMATHSHARED_EXPORT __declspec(dllimport)
#endif

#endif // MY_GISMATH_GLOBAL_H

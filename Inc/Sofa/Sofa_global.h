#ifndef SOFA_GLOBAL_H
#define SOFA_GLOBAL_H

#include <QtCore/QtGlobal>

#if defined(SOFA_LIBRARY)
#  define SOFASHARED_EXPORT Q_DECL_EXPORT
#else
#  define SOFASHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // SOFA_GLOBAL_H

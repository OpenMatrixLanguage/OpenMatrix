/////////////////////////////////////////////////////////////////////////
 // File : OmlPythonToolboxDefs.h
 // Copyright (c) 2017 solidThinking Inc.  All Rights Reserved
 // Contains trade secrets of solidThinking, Inc.  Copyright notice
 // does not imply publication.  Decompilation or disassembly of this
 // software is strictly prohibited.
 ////////////////////////////////////////////////////////////////////////
#ifndef OML_PYTHON_TOOLBOX_DEFS_H__
#define OML_PYTHON_TOOLBOX_DEFS_H__

#ifdef OS_WIN
#   ifdef OMLPYTHONBRIDGE_EXPORTS
#      undef  OMLPYTHONBRIDGE_DECLS
#      define OMLPYTHONBRIDGE_DECLS __declspec(dllexport)
#   else
#      undef  OMLPYTHONBRIDGE_DECLS
#      define OMLPYTHONBRIDGE_DECLS __declspec(dllimport)
#   endif  // OMLPYTHONBRIDGE_EXPORTS
#else
#   undef  OMLPYTHONBRIDGE_DECLS
#   define OMLPYTHONBRIDGE_DECLS
#endif // OS_WIN

#ifndef NULL
#   define NULL 0
#endif

#endif

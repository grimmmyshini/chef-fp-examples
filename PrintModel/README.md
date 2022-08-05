# Print Model

### Building the Model

```Nix
CladBuildDir="PATH_TO_CLAD_BUILD_DIR"
CladIncludeDir="PATH_TO_CLAD_INCLUDE_DIR"
LLVMIncludeDir="PATH_TO_LLVM_INCLUDE_DIR"
ClangIncludeDir="PATH_TO_CLANG_INCLUDE_DIR"
ToolsClangIncludeDir="PATH_TO_TOOLS_CLANG_INCLUDE_DIR"
ClangBuildDir="PATH_TO_CLANG_BUILD_DIR"
CladSrcDir="PATH_TO_CLAD_SRC_DIR"

/usr/local/bin/c++ -DCLAD_INSTDIR_INCL=\"$CladBuildDir/include\" -DCLAD_SRCDIR_INCL=\"$CladIncludeDir\" -D_GNU_SOURCE -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS -DcladCustomModelPlugin_EXPORTS -I$CladBuildDir/demos/ErrorEstimation/CustomModel -I$CladSrcDir/demos/ErrorEstimation/CustomModel -I$CladBuildDir/include -ICladBuildDir="$CladBuildDir" -IClangBuildDir -I$LLVMIncludeDir -I$ClangBuildDir/include -I$ToolsClangIncludeDir -fPIC -fvisibility-inlines-hidden -Werror=date-time -std=c++11 -w -fno-common -Woverloaded-virtual -Wcast-qual -fno-strict-aliasing -pedantic -Wno-long-long -Wall -W -Wno-unused-parameter -Wwrite-strings -g -fPIC   -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS -fno-rtti -MD -MT PrintModel.cpp.o -MF PrintModel.cpp.o.d -o PrintModel.cpp.o -c PrintModel.cpp

/usr/local/bin/c++ -fPIC  -fPIC -fvisibility-inlines-hidden -Werror=date-time -std=c++11 -w -fno-common -Woverloaded-virtual -Wcast-qual -fno-strict-aliasing -pedantic -Wno-long-long -Wall -W -Wno-unused-parameter -Wwrite-strings -g  -Wl,-z,defs -Wl,-z,nodelete -shared -Wl,-soname,libPrintModel.so -o libPrintModel.so PrintModel.cpp.o  -Wl,--unresolved-symbols=ignore-in-object-files
```
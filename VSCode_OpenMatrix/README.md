# Compose OML - OpenMatrix Language Extension for Visual Studio Code
This VS Code Extension provides support for [Open Sourced OpenMatrix Language](https://www.openmatrix.org/) and [Compose OML](https://altair.com/compose).

## Table of Contents

## Features include
Extension provide Basic Compose OML Features, without installing Altair Compose or OpenSourced OpenMatrix Language. 

1. Syntax highlighting
   
	![Syntax highlighting](https://github.com/OpenMatrixLanguage/OpenMatrix/tree/master/VSCode_OpenMatrix/images/syntax_highlighting.png))

2. Declarative language features - 
   1. Comments
   2. Brackets
   3. Indentation rules
   4. collapsible/folding
   

3. Advance Features
   1. Auto-completion
   2. Execution
   3. Plotting & Visualization

![language_features](https://github.com/OpenMatrixLanguage/OpenMatrix/tree/master/VSCode_OpenMatrix/images/code_features.gif)
      
Advanced Features are supported, if either of below is installed or build on your system 
   1. [Altair Compose](https://altairone.com/Marketplace?queryText=compose&tab=Info&app=Compose) 2023 or later
   2. [OpenMatrix 1.0.13](https://www.openmatrix.org/) or later
   



   
## Installing Compose OML - OpenMatrix Language Extension
---
you can install the Extension from within VS Code Extension (ctrl+shift+X) by searching for 'OML' or from [visual studio code Marketplace]().

### Configuring Compose OML - OpenMatrix Language with Altair Compose

To configure the extension with Altair Compose or Open Source OpenMatrix Language, go to extension settings and change Path to OML executable <span style="color:orange">
OML_EXE
</span> to compose.exe.

```
C:\Program Files\Altair\2023\Compose2023\hwx\bin\win64\Compose.exe
```

![OML_EXE](https://github.com/OpenMatrixLanguage/OpenMatrix/tree/master/VSCode_OpenMatrix/images/Configuration_oml_exe.gif)


### Configuring Compose OML - OpenMatrix Language with Open Source OpenMatrix Language
1. Follow instructions [windows](https://github.com/OpenMatrixLanguage/OpenMatrix/blob/master/INSTALL_WIN) or [Linux](https://github.com/OpenMatrixLanguage/OpenMatrix/blob/master/INSTALL_LINUX) to build OpenMatrix 
2. Set the path to the OML executable <span style="color:orange">
OML_EXE
</span> to omlcompose.exe.
for example
```
C:\Program Files\OpenMatrix_1.0.13_win64\OpenMatrix\src\bin\win64\omlconsole.exe
```

3. set the third-party library paths by going through the extensions settings.

      | Settings  | Description | Relative Path
      | ------------- | ------------- | ------------- |
      | OML_INTEL_COMP  |Path to intel compiler |/intel/compilers_and_libraries_2019.5.281/windows/redist/intel64_win/compiler
      | OML_INTEL_MKL |Path to mkl compiler |/intel/compilers_and_libraries_2019.5.281/windows/redist/intel64_win/mkl
      | OML_FFTW   |Path to fftw executable | /fftw/fftw-3.2.2/fftw-3.2.2-libs/x64/Release
      | OML_MATIO |Path to matio executable |/matio/matio-1.5.19/win64/bin
      | OML_HDF |Path to hdf5 executable |/hdf/hdf5-1.12.0/win64/bin
      | OML_QHULL |Path to qhull executable |/qhull/qhull-2015.2/bin



### Questions and Support
We encourage All feedback. If you face any issues please reach out to [Altair Compose Community Forum](https://community.altair.com/community?sys_id=d4e4e9d61b9d0c50a028542d1e4bcb47&view=sp&id=community_topic&table=sn_communities_topic)





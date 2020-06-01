# diploma
creating project structure for compiling and parallelizing fortran programs.

The solution was developed in Microsoft Visual Studio 2017.


All files of fortran program should be listed in one file (for example projectFiles.txt).

Fortran program files should be placed in the same directory.

Path to projectFiles.txt is passing in project as comand line argument.

If there is no arguments, then projectFile.txt should be in current directory.

Line size in free form files is not limited by default. In fixed form files line containes 72 characters by default.

In comand line you can pass options: -free132 and -fix132. They both will set line size to 132 characters.


You can run project in Microsoft Visual Studio.

Also you can run it in comand line.

For example from "diploma" directory (windows PowerShell):

.\x64\Debug\DependensyResolution.exe .\tests\t1\projectFiles.txt

In this directory will be created files: fileLevel.txt, moduleInit.txt, usedModules.txt, projectFilesLst.txt.


1-14 tests: in ProjectFileDescriptor.cpp uncomment line 51

17 test: for more information in ProjectFileDescriptor.cpp uncomment lines 97 and 116 and in ProjectStruct.cpp uncomment line 78


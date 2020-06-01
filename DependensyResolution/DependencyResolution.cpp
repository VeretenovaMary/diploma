#include "ProjectStruct.h"

#include <iostream>
#include <fstream>

int main(int argc, char **argv) {

	int fixPos = 72;
	int freePos = 0; // неограниченное количество символов
	std::string filePath = "tmptest/ProjectFiles.txt";
	for (int i = 1; i < argc; ++i) {
		const char *curr_arg = argv[i];
		switch (curr_arg[0]) {
		case '-':
			if (std::string(curr_arg) == "-fix132") {
				fixPos = 132;
			}
			else if (std::string(curr_arg) == "-free132") {
				freePos = 132;
			}
			break;
		default:
			filePath = std::string(curr_arg);
			break;
		}
	}

	ProjectStruct DR(filePath, fixPos, freePos);
	int result = DR.parser();
	if (result == 1) {
		return 1;
	}

	std::streambuf *coutbuf = std::cout.rdbuf(); // запоминаем старый буфер
	// projectFilesLst
	std::ofstream out("projectFilesLst.txt");	// открываем файл для вывода 
	std::cout.rdbuf(out.rdbuf());				// перенаправляем поток вывода
	DR.printProjectFilesLst();
	out.close();								// закрываем файл
	// fileLevel
	out = std::ofstream("fileLevel.txt");
	std::cout.rdbuf(out.rdbuf());
	DR.printFilesLevels();
	out.close();
	// moduleInit
	out = std::ofstream("moduleInit.txt");
	std::cout.rdbuf(out.rdbuf());
	DR.printModuleInit();
	out.close();
	// usedModules
	out = std::ofstream("usedModules.txt");
	std::cout.rdbuf(out.rdbuf());
	DR.printUsedModules();
	out.close();

	std::cout.rdbuf(coutbuf);

	return 0;
}

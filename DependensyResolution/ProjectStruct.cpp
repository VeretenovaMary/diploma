#include "ProjectStruct.h"
#include "FortranFile.h"
#include "ProjectFileDescriptor.h"

#include <list>
#include <stack>
#include <fstream>
#include <iostream>


ProjectStruct::ProjectStruct(std::string filePath, int fixLineLen, int freeLineLen) {
	_projectFiles = extractProjectFiles(filePath, fixLineLen, freeLineLen);
	_workDir = extractWorkDir(filePath);
	_filesLevels = {};
	_usedModules = {};
	_moduleInit = {};
}

ProjectStruct::~ProjectStruct() {}

void ProjectStruct::setLineLen(std::string filePath, int strLen) {
	std::map<std::string, std::pair<formType, int>>::iterator it = _projectFiles.find(filePath);
	if (it != _projectFiles.end())
	{
		it->second = std::pair<formType, int>(it->second.first, strLen);
	}
}

int ProjectStruct::parser() {
	if (_projectFiles.empty()) {
		return 0;
	}
	if (_projectFiles.begin()->first == "") {
		std::cerr << "error: can't open file with project files names" << std::endl;
		return 1;
	}
	std::stack <ProjectFileDescriptor*> stackSH;
	std::map<std::string, std::pair<formType, int>>::iterator it = _projectFiles.begin();
	for (it; it != _projectFiles.end(); it++) {
		std::string fileName = it->first;
		std::pair<formType, int> pairFP = it->second;
		ProjectFileDescriptor PFD(_workDir + fileName, pairFP.first, pairFP.second);
		if (PFD.getFilePath() == "") { return 1; }
		std::set<std::string> modLinks = {};

		std::set<std::string> filesInStack = {};

		while (true) {
			std::pair<stmtType, std::string> stmt = PFD.nextKeyWordStmt();

			stmtType firstLex = stmt.first;
			std::string nameLex = stmt.second;

			if (firstLex == module) {
				if (modLinks.count(nameLex)) {
					std::cerr << "error: in the same file module " << nameLex.c_str() << " was used before its description" << std::endl;
					return 1;
				}
				_moduleInit.insert(std::pair<std::string, std::string>(nameLex, fileName));
			}
			if (firstLex == use) {
				modLinks.insert(nameLex);
			}
			if (firstLex == include) {
				if (fileName == nameLex or filesInStack.count(_workDir + nameLex)) {
					std::cerr << "error: file " << nameLex.c_str() << " includes itself" << std::endl;
					return 1;
				}
				filesInStack.insert(PFD.getFilePath());
				if (extractFileForm(nameLex) != pairFP.first) {
					std::cerr << "error: file " << nameLex.c_str() << " includes file in different form" << std::endl;
					return 1;
				}
				ProjectFileDescriptor tmpPFD(_workDir + nameLex, pairFP.first, pairFP.second);
				if (tmpPFD.getFilePath() == "") { return 1; }
				stackSH.push(new ProjectFileDescriptor(PFD));
				PFD = tmpPFD;
				// std::cout << "include level = " << stackSH.size() << std::endl;
			}
			if (firstLex == eof) {
				if (stackSH.empty()) {
					break;
				}
				else {
					PFD = *stackSH.top();
					stackSH.pop();
					filesInStack.erase(PFD.getFilePath());
				}
			}
		}
		_usedModules.insert(std::pair<std::string, std::set<std::string>>(fileName, modLinks));
	}
	for (std::map<std::string, std::pair<formType, int>>::iterator it = _projectFiles.begin(); it != _projectFiles.end(); it++) {
		std::string fileName = it->first;
		_filesLevels.insert(std::pair<int, std::string>(defineFileLevel(fileName, {}), fileName));
	}
	return 0;
}


int ProjectStruct::defineFileLevel(std::string fileName, std::set<std::string> dependentModules) {
	int maxLvl = dependentModules.size();
	std::set<std::string> modLinks = _usedModules.at(fileName);
	for (std::set<std::string>::iterator it = modLinks.begin(); it != modLinks.end(); it++) {
		std::string modName = *it;
		if (!_moduleInit.count(modName)) {
			std::cerr << "error: module " << (modName).c_str() << " not declared" << std::endl;
			return -1;
		}
		std::string dependFile = _moduleInit.at(modName);
		if (dependFile == fileName) {
			continue; // если ссылка на модуль в том же файле
		}
		if (dependentModules.count(dependFile)) {
			std::cerr << "error: module " << modName.c_str() << " includes itself" << std::endl;
			return -1;
		}
		dependentModules.insert(dependFile);
		int resLvl = defineFileLevel(dependFile, dependentModules);
		if (resLvl < 0) {
			return resLvl;
		}
		if (resLvl > maxLvl) {
			maxLvl = resLvl;
		}
		dependentModules.erase(dependFile);
	}
	return maxLvl;
}

std::map <std::string, std::pair<formType, int>> ProjectStruct::extractProjectFiles(std::string filePath, int fixStrLen, int freeStrLen) {
	std::map <std::string, std::pair<formType, int>> result = {};
	std::ifstream file(filePath);
	std::string line;
	if (!file) {
		std::pair<formType, int> nullValue = std::pair<formType, int>(fixForm, 0);
		result.insert(std::pair<std::string, std::pair<formType, int>>("", nullValue));
		return result;
	}
	while (std::getline(file, line)) {
		if (!line.length()) { // пустая строка
			continue;
		}
		if (line.front() == '#') { // комментарий
			continue;
		}
		formType form;
		form = extractFileForm(line);
		int strLen = fixStrLen;
		if (form == freeForm) {
			strLen = freeStrLen;
		}
		std::pair<formType, int> formTypeStrLen = std::pair<formType, int>(form, strLen);
		// <имя файла, <форма записи файла, максимальное количество символов в строке>>
		result.insert(std::pair<std::string, std::pair<formType, int>>(line, formTypeStrLen));
	}
	return result;
}

formType ProjectStruct::extractFileForm(std::string filePath) {
	for (int dotPos = filePath.length() - 1; dotPos >= 0; dotPos--) {
		if (filePath[dotPos] == '.') {
			std::string extension = filePath.substr(dotPos + 1, filePath.length() - dotPos - 1);
			if (extension == "f90") {
				return freeForm;
			}
			return fixForm;
		}
	}
	return fixForm; // default
}

std::string ProjectStruct::extractWorkDir(std::string filePath) {
	for (int splitPos = filePath.length() - 1; splitPos >= 0; splitPos--) {
		if (filePath[splitPos] == '/' or filePath[splitPos] == '\\') {
			return filePath.substr(0, splitPos + 1);
		}
	}
	return "";
}


void ProjectStruct::printUsedModules() {
	if (!_usedModules.empty()) {
		std::cout << std::endl << "usedModules :" << std::endl << std::endl;
	}
	int elemNum = 1;
	std::map<std::string, std::set<std::string>>::iterator it = _usedModules.begin();
	for (it; it != _usedModules.end(); it++) {
		std::set<std::string> dependFiles = it->second;
		if (dependFiles.empty())
		{
			continue;
		}
		std::cout << elemNum << " " << it->first.c_str() << std::endl;
		elemNum++;
		std::set<std::string>::iterator str = dependFiles.begin();
		for (it; str != dependFiles.end(); str++) {
			std::cout << "    " << (*str).c_str() << std::endl;
		}
		std::cout << std::endl;
	}
}

void ProjectStruct::printModuleInit() {
	if (!_moduleInit.empty()) {
		std::cout << std::endl << "moduleInit :" << std::endl << std::endl;
	}
	int elemNum = 1;
	std::map<std::string, std::string>::iterator it = _moduleInit.begin();
	for (it; it != _moduleInit.end(); it++) {
		std::cout << elemNum << " " << it->first.c_str() << std::endl;
		elemNum++;
		std::cout << "    " << it->second.c_str() << std::endl;
		std::cout << std::endl;
	}
}

void ProjectStruct::printFilesLevels() {
	if (!_filesLevels.empty()) {
		std::cout << std::endl << "file level, file name :" << std::endl << std::endl;
	}
	std::map<int, std::string>::iterator it = _filesLevels.begin();
	for (it; it != _filesLevels.end(); it++) {
		std::cout << it->first << ", " << it->second.c_str() << std::endl;
	}
}

void ProjectStruct::printProjectFilesLst() {
	std::map<int, std::string>::iterator it = _filesLevels.begin();
	for (it; it != _filesLevels.end(); it++) {
		std::cout << it->second.c_str() << std::endl;
	}
}


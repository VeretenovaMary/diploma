#include "FortranFile.h"
#include "ProjectFileDescriptor.h"

#include <fstream>
#include <string>
#include <map>
#include <set>

#pragma once
class ProjectStruct {
	std::string _workDir;
	std::map<std::string, std::pair<formType, int>> _projectFiles;

	std::map <std::string, std::string> _moduleInit;
	std::map <std::string, std::set<std::string>> _usedModules;
	std::multimap <int, std::string> _filesLevels;
public:
	ProjectStruct(std::string s = "projectFiles.txt", int fixLineLen = 72, int freeLineLen = 0);
	~ProjectStruct();

	int parser();
	void setLineLen(std::string, int);
	int defineFileLevel(std::string, std::set<std::string>);

	void printProjectFilesLst();
	void printFilesLevels();
	void printModuleInit();
	void printUsedModules();

	std::map <std::string, std::pair<formType, int>> extractProjectFiles(std::string, int, int);
	formType extractFileForm(std::string);
	std::string extractWorkDir(std::string);
};


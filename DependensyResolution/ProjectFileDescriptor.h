#include "FortranFile.h"

#pragma once
class ProjectFileDescriptor {
private:
	char _c;
	int _pos;
	FortranFile _FF;

	int _ibLevel; // interface block level
	bool _isNewLine; // current statement on new line
	bool _isNextNewLine; // next statement on new line
public:
	ProjectFileDescriptor(std::string, formType, int);
	ProjectFileDescriptor(ProjectFileDescriptor&);
	~ProjectFileDescriptor();
	ProjectFileDescriptor & operator= (ProjectFileDescriptor&);

	std::string getFilePath() { return _FF.getFilePath(); }

	std::string nextStmt();
	std::pair<stmtType, std::string> nextKeyWordStmt();
	void nextChar() { _c = _FF.nextChar(); _pos++; }

	std::string constructName(std::string, int, bool isFileName);

};


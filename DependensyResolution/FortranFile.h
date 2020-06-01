#include <map>

enum formType { freeForm, fixForm };
enum stmtType { module, use, include, eof };

#pragma once
class FortranFile {
	std::string _filePath;
	formType _form;
	int _lineLen;
	FILE *_fd;
public:
	FortranFile(std::string, formType, int);
	FortranFile(FortranFile &);
	~FortranFile();
	FortranFile& operator= (FortranFile&);

	std::string getFilePath() { return _filePath; }
	formType getForm() { return _form; }
	int getLineLen() { return _lineLen; }


	char nextChar();
};


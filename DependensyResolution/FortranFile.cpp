#include "FortranFile.h"

#include <iostream>
#include <fstream>
#include <string>

FortranFile::FortranFile(std::string filePath, formType formType, int strLen) :
	_filePath(filePath), _form(formType), _lineLen(strLen)
{
	if (fopen_s(&_fd, _filePath.c_str(), "r")) {
		std::cerr << "error: cannot open file " << _filePath << std::endl;
		_filePath = "";
		_fd = NULL;
	}
}

FortranFile::~FortranFile() {
	if (_fd != NULL) {
		fclose(_fd);
	}
}

FortranFile::FortranFile(FortranFile &FH) {
	this->_filePath = FH._filePath;
	this->_form = FH._form;
	this->_fd = FH._fd;
	this->_lineLen = FH._lineLen;
	FH._fd = NULL;
}

FortranFile& FortranFile::operator= (FortranFile& FH) {
	this->_filePath = FH._filePath;
	this->_form = FH._form;
	this->_fd = FH._fd;
	this->_lineLen = FH._lineLen;
	FH._fd = NULL;
	return *this;
}

char FortranFile::nextChar() {
	return fgetc(_fd);
}


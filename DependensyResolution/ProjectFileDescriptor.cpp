#include <iostream>
#include "ProjectFileDescriptor.h"


ProjectFileDescriptor::ProjectFileDescriptor(std::string filePath, formType formType, int strLen) :
	_c(' '), _pos(0), _FF(filePath, formType, strLen),
	_ibLevel(0),
	_isNewLine(false), _isNextNewLine(true)
{}

ProjectFileDescriptor::ProjectFileDescriptor(ProjectFileDescriptor & PFD) :
	_c(PFD._c), _pos(PFD._pos), _FF(PFD._FF),
	_ibLevel(PFD._ibLevel),
	_isNewLine(PFD._isNewLine), _isNextNewLine(PFD._isNextNewLine)
{}

ProjectFileDescriptor& ProjectFileDescriptor::operator= (ProjectFileDescriptor &PFD)
{
	this->_c = PFD._c;
	this->_pos = PFD._pos;
	this->_FF = PFD._FF;
	this->_ibLevel = PFD._ibLevel;
	this->_isNewLine = PFD._isNewLine;
	this->_isNextNewLine = PFD._isNextNewLine;
	return *this;
}

ProjectFileDescriptor::~ProjectFileDescriptor() {}

std::string ProjectFileDescriptor::constructName(std::string str, int namePos, bool isFileName = false) {
	std::string name;
	int i = namePos;
	char c = str[i];
	if (isalpha(c)) {
		while (i < str.length() and (isalnum(c) or c == '_') or (isFileName and c == '.')) {
			if (c == '.') {
				isFileName = false;
			}
			name += tolower(c);
			i += 1;
			c = str[i];
		}
	}
	return name;
}

std::pair<stmtType, std::string> ProjectFileDescriptor::nextKeyWordStmt() {
	std::string stmt = nextStmt();
	std::string nameLex;
	while (stmt != "") {
		// std::cout << stmt.c_str() << std::endl;
		nameLex = "";
		if (stmt.compare(0, 3, "use") == 0) {
			nameLex = constructName(stmt, 3);
			if (nameLex.length() != 0 and
				(stmt.length() == 3 + nameLex.length() or stmt[3 + nameLex.length()] == ',')) {
				return std::pair<stmtType, std::string>(use, nameLex);
			}
		}
		else if (_ibLevel == 0 and stmt.compare(0, 6, "module") == 0) {
			nameLex = constructName(stmt, 6);
			if (nameLex.length() != 0 and stmt.length() == 6 + nameLex.length()) {
				return std::pair<stmtType, std::string>(module, nameLex);
			}
		}
		else if (_isNewLine and _isNextNewLine and stmt.compare(0, 7, "include") == 0) {
			int tmpi = 0;
			if (stmt.length() > 9) {
				char c = stmt[7];
				if (c == '\'' or c == '\"') {
					nameLex = constructName(stmt, 8, true);
					tmpi = 7 + 1 + nameLex.length();
					if (nameLex.length() != 0 and stmt.length() == tmpi + 1 and stmt[tmpi] == c) {
						return std::pair<stmtType, std::string>(include, nameLex);
					}
				}
			}
		}
		else if (stmt.compare(0, 9, "interface") == 0) {
			bool openBracket = false;
			if (stmt.length() == 9) {
				_ibLevel++;
			}
			else if (isalpha(stmt[9])) {
				int i = 10;
				for (i; i < stmt.length() and (isalnum(stmt[i]) or stmt[i] == '_'); i++);
				if (i == stmt.length()) {
					_ibLevel++;
				}
				else if (stmt[i] == '(') {
					for (i; i < stmt.length() and stmt[i] != ')'; i++);
					if (i == stmt.length() - 1) {
						_ibLevel++;
					}
				}
			}
			// std::cout << "ibLevel = " << _ibLevel << std::endl;
		}
		else if (stmt.compare(0, 12, "endinterface") == 0) {
			if (stmt.length() == 12) {
				_ibLevel--;
			}
			else if (isalpha(stmt[12])) {
				int i = 13;
				for (i; i < stmt.length() and (isalnum(stmt[i]) or stmt[i] == '_'); i++);
				if (i == stmt.length()) {
					_ibLevel--;
				}
				else if (stmt[i] == '(') {
					for (i; i < stmt.length() and stmt[i] != ')'; i++);
					if (i == stmt.length() - 1) {
						_ibLevel--;
					}
				}
			}
			// std::cout << "ibLevel = " << _ibLevel << std::endl;
		}
		stmt = nextStmt();
	}
	return std::pair<stmtType, std::string>(eof, "");
}

std::string ProjectFileDescriptor::nextStmt() {
	std::string stmt = "";
	_isNewLine = _isNextNewLine;
	_isNextNewLine = false;

	bool isNewStmt = false;
	bool isComment = false;
	bool isCharContext = false;
	bool isStmtContinue = false;
	char quoteType = ' ';
	while (_c != EOF) {
		if (_FF.getForm() == fixForm and _pos == 1 and (tolower(_c) == 'c' or _c == '*')) {
			isComment = true;
			nextChar();
			continue;
		}
		int mp = _FF.getLineLen();
		if (mp != 0 and _pos > mp and !(isComment or isCharContext)) {
			isComment = true;
			continue;
		}
		if (_c == ' ' and !(isComment or isCharContext)) {
			nextChar();
			continue;
		}
		if (!isComment and (_c == ' ' or _c == '0') and _pos <= 6 and _FF.getForm() == fixForm) {
			nextChar();
			continue;
		}
		if (_c == '\r') {
			nextChar();
			_pos--;
			continue;
		}
		if (_c == '\t') {
			if (_FF.getForm() == fixForm and _pos <= 6 and !(isComment or isCharContext)) {
				if (_pos == 6) {
					isStmtContinue = true;
				}
				nextChar();
				_pos = 7;
				if (isdigit(_c)) {
					_pos = 6;
					nextChar();
				}
				continue;
			}
			nextChar(); // _pos >> 6 or freeForm
			continue;
		}
		if (_c == '\n') {
			nextChar();
			_pos = 1;
			isNewStmt = true;
			_isNextNewLine = true;
			isComment = false;
			continue;
		}
		if (!isComment
			and (_FF.getForm() == freeForm and _c == '&'
				or _FF.getForm() == fixForm and _pos == 6 and _c != ' ' and _c != '0' and _c != '\t')) {
			if (isStmtContinue and !_isNextNewLine and _FF.getForm() == freeForm) {
				stmt += '&';
			}
			isStmtContinue = true;
			_isNextNewLine = false;
			nextChar();
			continue;
		}
		if (_c == ';' and !(isComment or isCharContext)) {
			nextChar();
			isNewStmt = true;
			continue;
		}
		if ((_c == '\'' or _c == '\"') and !isComment) {
			if (!isCharContext) {
				isCharContext = true;
				quoteType = _c;
				stmt += _c;
				nextChar();
				continue;
			}
			// isCharContext
			if (_c == quoteType) {
				stmt += _c;
				nextChar();
				isCharContext = false;
				quoteType = ' ';
				continue;
			}
			stmt += _c;
			nextChar();
			continue;
		}

		if (_c == '!' and !isComment) { // обработка комментария
			if (isStmtContinue and isCharContext and !isNewStmt) {
				stmt += '&';
			}
			if (isCharContext) {
				if (isStmtContinue and _isNextNewLine) {
					_isNextNewLine = false;
					isComment = true;
				}
				else {
					stmt += '!';
					isStmtContinue = false;
				}
			}
			else {
				isComment = true;
			}
			nextChar();
			isNewStmt = false;
			continue;
		}

		if (isComment) {
			nextChar();
			continue;
		}

		if (isdigit(_c) and _isNewLine and stmt.length() == 0) {
			// label
			nextChar();
			continue;
		}

		if (isStmtContinue) {
			if (!isNewStmt and _FF.getForm() == freeForm) {
				stmt += '&';
			}
			isNewStmt = false;
			_isNextNewLine = false;
			isStmtContinue = false;
		}
		if (isNewStmt and stmt.length() != 0) {
			return stmt;
		}
		isNewStmt = false;
		stmt += tolower(_c);
		nextChar();
	}
	_isNextNewLine = true;
	return stmt;
}

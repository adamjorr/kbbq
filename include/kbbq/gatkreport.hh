#ifndef KBBQ_GATKREPORT_HH
#define KBBQ_GATKREPORT_HH

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>

namespace gatkreport{

enum column_type{STRING = 0, FLOAT, INT};

class TableValue{
private:
	column_type type;
	union{
		long double f;
		unsigned long long i;
	}
	std::string s = "";
public:
	inline column_type current_type(){return type;}
	template<typename T> T get();
	template<>
	unsigned long long get<unsigned long long>(){
		return type == INT ? i : 0;
	}
	template<>
	long double get<long double>(){
		return type == FLOAT ? f : 0;
	}
	template<>
	std::string get<std::string>(){
		return type == STRING ? s : "";
	}
	inline TableValue& set(unsigned long long i){this->i = i; this->type = INT; return this;}
	inline TableValue& set(long double f){this->f = f; this->type = FLOAT; return this;}
	inline TableValue& set(std::string s){this->s = s; this->type = STRING; return this;}
};

class TableRow{
private:
	std::vector<TableValue> columns;
	std::vector<column_type> types;
public:

}

class GATKTable{
private:
	std::vector<TableRow> rows;
};

class GATKReport{
private:
	std::vector<GATKTable> tables;
	std::string version;
public:
	GATKReport(){}
	GATKReport(const std::vector<GATKTable>& tables): tables(tables) {}
	GATKReport(const std::string& filename): tables(), version(){
		std::ifstream inf(filename);
		std::string tablestr;
		std::string headerline;
		std::string versionprefix = "#:GATKReport.v";
		std::getline(inf, headerline);
		if(headerline.substr(0,versionprefix.length()) != versionprefix){
			throw std::invalid_argument("Error: Unable to parse first line of input file " +
				filename + " ; Ensure it is a valid GATKReport.");
		}
		headerline.erase(0,versionprefix.length());
		size_t colon_pos = headerline.find(':');
		version = headerline.substr(0, colon_pos);
		size_t ntables = std::stoull(headerline.substr(colon_pos+1, std::string::npos));
		for(std::string line; std::getline(inf, line);){
			if(line.empty()){ //blank line delimits tables
				tablestr.pop_back(); //get rid of ending newline
				tables.emplace_back(tablestr);
				tablestr.clear();
			} else {
				tablestr += line + '\n';
			}
		}
		if(ntables != tables.size()){
			throw std::invalid_argument("Error: Unable to parse first line of input file " +
				filename + " ; Ensure it is a valid GATKReport.");
		}
	}
	//return the header string without a trailing newline
	inline std::string headerstring() const{
		return "#:GATKReport.v" + version + ":" + std::to_string(tables.size());
	}
	inline explicit operator std::string() const {
		std::string ret = this->headerstring() + "\n";
		for(const std::string& s : tables){ret += s + "\n";}
		return ret;
	}
};

}

#endif
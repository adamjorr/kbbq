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
	std::string s;
public:
	TableValue(): type(), f(), s() {}
	TableValue(unsigned long long i): TableValue() {this->set(i);}
	TableValue(long double f): TableValue() {this->set(f);}
	TableValue(std::string s): TableValue() {this->set(s);}

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
	std::vector<std::string> headers;
public:
	TableRow(): columns(), types() {}
	TableRow(const std::string& line, const std::vector<column_type>& t, const std::vector<std::string>& h):
		types(t), columns(), headers(h)
	{
		std::istringstream is(line);
		std::string token;
		for(column_type t : types){
			is >> token;
			switch(t){
				case STRING: columns.emplace_back(token);
					break;
				case INT: columns.emplace_back(std::stoull(token));
					break;
				case FLOAT: columns.emplace_back(std::stold(token));
					break;
			}
		}
	}
	//space-delimited values in the row.
	inline explicit operator std::string() const {
		std::string ret;
		for(size_t i = 0; i < columns.size(); ++i){
			std::string column_val = "";
			switch(types[i]){
				case FLOAT: long double val = columns[i].get<long double>();
					std::string format = headers[i] == "Errors" ? "%.2f" : "%.4f"
					int outsize = std::snprintf(nullptr, 0, format.c_str(), val);
					column_val.resize(outsize + 1);
					std::snprintf(&column_val[0], column_val.size(), format.c_str(), val);
					break;
				case INT: column_val = std::to_string(columns[i].get<unsigned long long>());
					break
				case STRING:
				default:
					column_val = columns[i].get<std::string>;
			}
			ret += column_val + " ";
		}
		ret.pop_back(); //remove trailing space
		return ret;
	}
}

class GATKTable{
private:
	std::vector<TableRow> rows;
	std::string title;
	std::string description;
	std::vector<std::string> headers;
	std::vector<column_type> types;
	std::vector<size_t> widths;
	/*
	const std::map<std::string, std::string> precision = {
		{"EmpiricalQuality", ".4"},
		{"EstimatedQReported",".4"},
		{"Errors",".2"}
	};
	*/
public:
	GATKTable(): rows(), title(), description(), headers(), types() {}
	GATKTable(std::string tablestr):
		rows(), title(), description(), headers(), types()
	{
		std::istringstream in(tablestr);
		std::string headerline;
		std::string headerprefix = "#:GATKTable:"
		std::getline(in, headerline); //process first line
		if(headerline.substr(0,headerprefix.length()) != headerprefix){
			throw std::invalid_argument("Error: Unable to parse first line of table. " +
				"Ensure input is a valid GATKTable.");
		}
		headerline.erase(0,headerprefix.length());
		std::istringstream hdrin(headerline);
		std::string token;
		std::getline(hdrin, token, ':');
		size_t ncols = std::stoull(token);
		std::getline(hdrin, token, ':');
		size_t nrows = std::stoull(token);

		while(std::getline(hdrin, token, ':')){
			switch(token.pop_back()){
				case 'd': types.push_back(INT);
					break;
				case 'f': types.push_back(FLOAT);
					break;
				case ';': //do nothing
					break;
				case 's':
				default: types.push_back(STRING);
			}
		}
		if(types.size() != ncols){
			throw std::invalid_argument("Error: Number of types doesn't match header!" +
				"Ensure the type string includes all types.\n");
		}

		std::getline(in, headerline); //process second line; it's the title and description
		if(headerline.substr(0,headerprefix.length()) != headerprefix){
			throw std::invalid_argument("Error: Unable to parse second line of table. " +
				"Ensure input is a valid GATKTable.");
		}
		headerline.erase(0,headerprefix.length());
		hdrin.str(headerline);
		std::getline(hdrin, title, ':');
		std::getline(hdrin, description, ':');
		std::getline(in, headerline); //process the 3rd line; it's the column headers
		hdrin.str(headerline);
		std::copy(std::istream_iterator<std::string>(hdrin),
				std::istream_iterator<std::string>(),
				std::back_inserter(headers));
		if(types.size() != headers.size()){
			throw std::invalid_argument(
				"Error: Number of headers doesn't match the number of columns.\n");
		}
		while(std::getline(in, token)){
			rows.emplace_back(token, types, headers);
		}
		if(rows.size() != nrows){
			throw std::invalid_argument(
				"Error: Number of stated rows doesn't match number of actual rows.\n");
		}
	}
	//TODO: ctor, default and std::string
	//return the header string without a trailing newline
	inline std::string headerstring() const{
		std::string headerstr = "#:GATKTable:" + std::to_string(types.size()) +
			":" + std::to_string(rows.size()) + ":";
		for(size_t i = 0; i < types.size(); ++i){
			column_type t = types[i];
			std::string h = headers[i];
			std::string colcode = "";
			switch(t){
				case STRING: colcode = "%s";
					break;
				case INT: colcode = "%d";
					break;
				case FLOAT: colcode = h == "Errors" ? "%.2f" : "%.4f";
					break;
			}
			headerstr += colcode + ":";
		}
		return headerstr + ";";
	}

	//the title string without a trailing newline
	inline std::string titlestring() const{
		return "#:GATKTable:" + title + ":" + description;
	}

	inline std::vector<size_t> get_col_widths() const{
		std::vector<size_t> widths{headers.size(),0};
		std::transform(headers.begin(), headers.end(), widths.begin(),
			[](std::string str) -> size_t {return str.length();});
		for(size_t i = 0; i < rows.size(); ++i){
			std::istringstream rowin = rows[i];
			for(size_t j = 0; j < widths.size(); ++j){
				std::string colstr;
				rowin >> colstr;
				if(colstr.length() > widths[j]){
					widths[j] = colstr.length();
				}
			}
		}
		return widths;
	}

	inline explicit operator std::string() const {
		std::string ret = this->headerstring() + "\n";
		ret += this->titlestring() + "\n";
		std::ostringstream os();
		std::vector<size_t> widths = this->get_col_widths();
		for(size_t i = 0; i < widths.size(); ++i){
			os << std::setw(widths[i]) << headers[i] << "  "; 
		}
		ret += os.str();
		ret.erase(ret.end()-2);
		os.str(""); //reset output stream
		std::istringstream is("");
		std::string token;
		for(const std::string& s : rows){
			is.str(s);
			is >> token;
			os << std::setw(widths[0]) << token;
			for(size_t i = 1; i < widths.size(); ++i){
				is >> token;
				os << "  " << std::setw(widths[i]) << token;
			}
			os << "\n";
		}
		ret += os.str();
		return ret;
	}
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
			throw std::invalid_argument("Found " + std::to_string(ntables) + " tables in " +
				filename + ", but " + std::to_string(tables.size()) " were declared." +
				"Ensure the file is not truncated and adjust the declared number of tables.");
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
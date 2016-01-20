// this file is distributed under 
// MIT license
#include <gnuplot_wrap.h>
namespace GnuplotWrap{
	using namespace std;
	string ReplaceAll(const string&source, const string& from, const string& to) {
		string str=source;
		size_t start_pos = 0;
		while((start_pos = str.find(from, start_pos)) != string::npos) {
			str.replace(start_pos, from.length(), to);
			start_pos += to.length();
		}
		return str;
	}
	string ReplaceAll(string&& str, const string& from, const string& to){
		return ReplaceAll(str,from,to);
	}
	
	Plotter::Plotter(){
		counter=0;
		outpath="*";
	}
	Plotter::~Plotter(){
		ofstream script;
		string name=m_prefix+"_script.gnuplot";
		script.open((outpath+"/"+name).c_str());
		if(script.is_open()){
			for(string line:lines)
				script<<line<<"\n";
			script.close();
			name=string("gnuplot ")+name;
			string old=getcwd(NULL,0);
			chdir(outpath.c_str());
			system(name.c_str());
			chdir(old.c_str());
		}
		
	}
	Plotter &Plotter::Instance(){
		static Plotter m_instance;
		return m_instance;
	}
	void Plotter::SetOutput(string&&out,string&&prefix){
		if(outpath!="*")
			throw;
		outpath=out;
		m_prefix=prefix;
	}
	string&Plotter::OutPath()const{
		if(outpath=="*")
			throw;
		return const_cast<string&>(outpath);
	}
	string&Plotter::Prefix()const{
		if(outpath=="*")
			throw;
		return const_cast<string&>(m_prefix);
	}
	string Plotter::GetTerminal(){
		counter++;
		return string("set terminal pngcairo size 1024,868 font 'Verdana,18'\nset output '")+m_prefix+to_string(counter)+".png'";
	}
	Plotter &Plotter::operator<<(const string&line){
		lines.push_back(line);
		return *this;
	}
	Plotter& Plotter::operator<<(string&&line){return operator<<(line);}
};
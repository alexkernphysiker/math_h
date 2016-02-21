// this file is distributed under 
// MIT license
#include <gnuplot_wrap.h>
#include <math_h/error.h>
namespace GnuplotWrap{
	using namespace std;
	using namespace MathTemplates;
	Plotter::Plotter(){
		terminal_counter=0;
		filename_counter=0;
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
			throw Exception<Plotter>("Attempt to reset plot output settings");
		outpath=out;
		m_prefix=prefix;
	}
	const string&Plotter::OutPath()const{
		if(outpath=="*")
			throw Exception<Plotter>("Attempt to use Plotter without initializing");
		return outpath;
	}
	const string&Plotter::Prefix()const{
		if(outpath=="*")
			throw Exception<Plotter>("Attempt to use Plotter without initializing");
		return m_prefix;
	}
	string Plotter::GetTerminal(){
		terminal_counter++;
		return string("set terminal pngcairo size 1024,868 font 'Verdana,18'\nset output '")
		+m_prefix+to_string(terminal_counter)+".png'";
	}
	string Plotter::GetFileName(){
		filename_counter++;
		return m_prefix+"-data-"+to_string(filename_counter)+".txt";
	}
	Plotter &Plotter::operator<<(const string&line){
		lines.push_back(line);
		return *this;
	}
	Plotter& Plotter::operator<<(string&&line){return operator<<(line);}
};
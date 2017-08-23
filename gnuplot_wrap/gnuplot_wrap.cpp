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
		string name="!"+m_prefix+".gnuplot-script";
		ofstream str(string(outpath+"/"+name).c_str(),ofstream::out);
		for(string line:lines)
		    str<<line<<"\n";
		str.close();
		name="gnuplot "+name;
		string old=getcwd(NULL,0);
		chdir(outpath.c_str());
		system(name.c_str());
		chdir(old.c_str());
	}
	Plotter &Plotter::Instance(){
		static Plotter m_instance;
		return m_instance;
	}
	void Plotter::SetOutput(const string&out,const string&prefix){
		if(outpath!="*")
			throw Exception<Plotter>("Attempt to reset plot output settings");
		outpath=out;
		m_prefix=prefix;
	}
	const string Plotter::GetTerminal(const string&name){
	    if(outpath=="*")
		throw Exception<Plotter>("Attempt to use Plotter without initializing");
	    const string firstline="set terminal pngcairo size 1024,868 font 'Verdana,18'\n";
	    if(name==""){
		terminal_counter++;
		string cnt=to_string(terminal_counter);
		while(cnt.length()<5)cnt="0"+cnt;
		return firstline+"set output '!"+m_prefix+"-plot-"+cnt+".png'";
	    }else{
		return firstline+"set output '"+name+".png'";
	    }
	}
	const string Plotter::GetFileName( const std::string& name ){
	    if(outpath=="*")
		throw Exception<Plotter>("Attempt to use Plotter without initializing");
	    if(name==""){
		filename_counter++;
		string cnt=to_string(filename_counter);
		while(cnt.length()<6)cnt="0"+cnt;
		return string(".")+m_prefix+"-numeric-"+cnt+".txt";
	    }else{
		return name+".txt";
	    }
	}
	pair<const string,ofstream> Plotter::File(const string&name){
	    const string n=GetFileName(name);
	    return make_pair(n,ofstream((outpath+"/"+n).c_str(),ofstream::out));
	}
	ifstream Plotter::GetInput(const string&name){
	    if(name=="")throw Exception<Plotter>("Cannot get input without name");
	    return ifstream((outpath+"/"+GetFileName(name)).c_str());
	}
	Plotter &Plotter::operator<<(const string&line){
		lines.push_back(line);
		return *this;
	}
};

// this file is distributed under 
// GPL v 3.0 license
#include "gnuplot.h"
using namespace std;
string ReplaceAll(string str, const string& from, const string& to) {
	size_t start_pos = 0;
	while((start_pos = str.find(from, start_pos)) != string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
	return str;
}
Plotter::Plotter(){
	counter=0;
	outpath="*";
}
Plotter::~Plotter(){
	ofstream script;
	string name="script.gnuplot";
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
void Plotter::SetOutput(string out){
	if(outpath!="*")
		throw;
	outpath=out;
}
string Plotter::OutPath(){
	if(outpath=="*")
		throw;
	return outpath;
}
string Plotter::GetTerminal(){
	counter++;
	return string("set terminal pngcairo size 1024,868 enhanced monochrome font 'Verdana,18'\nset output '")+to_string(counter)+".png'";
}
Plotter &Plotter::operator<<(string line){
	lines.push_back(line);
	return *this;
}
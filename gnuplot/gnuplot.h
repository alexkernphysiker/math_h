// this file is distributed under 
// GPL v 3.0 license
#ifndef VIJVUSSC
#define VIJVUSSC
#include <vector>
#include <functional>
#include <memory>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <memory>
#include <unistd.h>
#include "../interpolate.h"
std::string ReplaceAll(std::string str, const std::string& from, const std::string& to);
class Plotter{
public:
	Plotter();
	~Plotter();
	static Plotter &Instance();
	void SetOutput(std::string out);
	std::string OutPath();
	std::string GetTerminal();
	Plotter &operator<<(std::string line);
private:
	std::vector<std::string> lines;
	unsigned int counter;
	std::string outpath;
};

template<class numt=double>class Plot{
public:
	typedef std::function<numt(numt)> FUNC;
	Plot();
	~Plot();
	Plot &operator<<(std::string line);
	Plot &Points(std::string name,std::vector<std::pair<numt,numt>>&&points);
	Plot &Points(std::string name,LinearInterpolation_fixedsize<numt>&&points);
	Plot &Points(std::string name,LinearInterpolation<numt>&&points);
	Plot &Points(std::string name,LinearInterpolation<numt>&&points,FUNC error);
	Plot &Function(std::string name,FUNC func,numt from,numt to,numt step);
	Plot &Function(std::string name,FUNC,FUNC error,numt from,numt to,numt step);
private:
	std::vector<std::string> lines;
	std::vector<std::string> plots;
};
template<class numt>
Plot<numt>::Plot(){
	operator<<(Plotter::Instance().GetTerminal());
}
template<class numt>
Plot<numt>& Plot<numt>::operator<<(std::string line){
	lines.push_back(line);
}
template<class numt>
Plot<numt>::~Plot(){
	for(std::string line:lines)
		Plotter::Instance()<<line;
	for(int i=0,n=plots.size();i<n;i++){
		std::string line=plots[i];
		if(i==0)
			line="plot "+line;
		if(i<(n-1))
			line+=",\\";
		Plotter::Instance()<<line;
	}
}
#define WRITE_TO(data) \
std::ofstream data;\
std::string filename=ReplaceAll(ReplaceAll(name," ","_"),"=","_")+".txt";\
data.open((Plotter::Instance().OutPath()+"/"+filename).c_str());\
if(data.is_open()){
#define END(data)\
	data.close();}
template<class numt>
Plot<numt>& Plot<numt>::Points(std::string name,std::vector<std::pair<numt,numt>>&&points){
	WRITE_TO(data)
	for(auto p:points)
		data<<p.first<<" "<<p.second<<"\n";
	plots.push_back("\""+filename+"\" using 1:2 title \""+name+"\"");
	END(data)
}
template<class numt>
Plot<numt>& Plot<numt>::Points(std::string name, LinearInterpolation_fixedsize<numt>&& points){
	WRITE_TO(data)
	for(int i=0,n=points.size();i<n;i++)
			data<<points.getX(i)<<" "<<points.getY(i)<<"\n";
	plots.push_back("\""+filename+"\" using 1:2 title \""+name+"\"");
	END(data)
}
template<class numt>
Plot<numt>& Plot<numt>::Points(std::string name, LinearInterpolation<numt>&&points){
	WRITE_TO(data)
	for(auto p:points)
		data<<p.first<<" "<<p.second<<"\n";
	plots.push_back("\""+filename+"\" using 1:2 title \""+name+"\"");
	END(data)
}
template<class numt>
Plot<numt>& Plot<numt>::Points(std::string name, LinearInterpolation<numt>&&points,FUNC error){
	WRITE_TO(data)
	for(auto p:points)
		data<<p.first<<" "<<p.second<<" "<<error(p.first)<<"\n";
	plots.push_back("\""+filename+"\" using 1:2:($2-$3):($2+$3) with yerrorbars title \""+name+"\"");
	END(data)
}
template<class numt>
Plot<numt>& Plot<numt>::Function(std::string name,FUNC func,numt from,numt to,numt step){
	WRITE_TO(out)
	for(double x=from;x<=to;x+=step)
		out<<x<<" "<<func(x)<<"\n";
	plots.push_back("\""+filename+"\" w l title \""+name+"\"");
	END(out)
}
template<class numt>
Plot<numt>& Plot<numt>::Function(std::string name, FUNC func, FUNC error, numt from, numt to, numt step){
	WRITE_TO(out)
	for(double x=from;x<=to;x+=step)
		out<<x<<" "<<func(x)<<" "<<error(x)<<"\n";
	plots.push_back("\""+filename+"\" using 1:2:($2-$3):($2+$3) with yerrorbars title \""+name+"\"");
	END(out)
}
#undef WRITE_TO
#undef END
#endif

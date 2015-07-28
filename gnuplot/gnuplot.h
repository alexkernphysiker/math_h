// this file is distributed under 
// GPL v 3.0 license
#ifndef VIJVUSSC
#define VIJVUSSC
#include <vector>
#include <functional>
#include <memory>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <memory>
#include <unistd.h>
#include <math.h>
#include <exception>
#include "../interpolate.h"
std::string ReplaceAll(std::string str, const std::string& from, const std::string& to);
class Plotter{
public:
	Plotter();
	~Plotter();
	static Plotter &Instance();
	void SetOutput(std::string out,std::string prefix="");
	std::string OutPath();
	std::string Prefix();
	std::string GetTerminal();
	Plotter &operator<<(std::string line);
private:
	std::vector<std::string> lines;
	unsigned int counter;
	std::string outpath;
	std::string m_prefix;
};
template<class numt>class Plot{
private:
	std::vector<std::string> lines;
	std::vector<std::string> plots;
public:
	typedef std::function<numt(numt)> FUNC;
	typedef std::function<void(std::ofstream&)> PLOTOUTPUT;
	Plot &operator<<(std::string line){
		lines.push_back(line);
		return *this;
	}
	Plot(){
		operator<<(Plotter::Instance().GetTerminal());
	}
	virtual ~Plot(){
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
public:
	Plot& OutputPlot(std::string name,PLOTOUTPUT delegate,std::string plotstr){
		std::ofstream data;
		std::string filename=ReplaceAll(ReplaceAll(name," ","_"),"=","_")+".txt";
		data.open((Plotter::Instance().OutPath()+"/"+filename).c_str());
		if(data.is_open()){
			delegate(data);
			data.close();
			std::string line="\""+filename+"\" ";
			line+=plotstr;
			line+=" title \"";
			line+=name;
			line+="\"";
			plots.push_back(line);
		}
		return *this;
	}
	Plot &Hist(std::string name,LinearInterpolation_fixedsize<numt>&&points){
		OutputPlot(name,[this,&points](std::ofstream&data){
			for(int i=0,n=points.size();i<n;i++)
				data<<points.getX(i)<<" "<<points.getY(i)<<"\n";
		},"using 1:2");
		return *this;
	}
	Plot &HistWithStdError(std::string name,LinearInterpolation_fixedsize<numt>&&points){
		OutputPlot(name,[&points](std::ofstream&data){
			for(int i=0,n=points.size();i<n;i++){
				double y=points.getY(i);
				if(y<0)throw std::exception();
				   double dy=sqrt(y);if(y<1)y=1;
				data<<points.getX(i)<<" "<<y<<" "<<dy<<"\n";
			}
		},"using 1:2:($2-$3):($2+$3) with yerrorbars");
		return *this;
	}
	Plot &Line(std::string name,FUNC func,numt from,numt to,numt step){
		OutputPlot(name,[func,from,to,step](std::ofstream&data){
			for(double x=from;x<=to;x+=step)
				data<<x<<" "<<func(x)<<"\n";
		},"w l");
		return *this;
	}
};
template<class numt,class Indexer>class PlotPoints:public Plot<numt>{
public:
	typedef std::function<numt(numt)> FUNC;
	typedef std::pair<numt,numt> PAIR;
	PlotPoints():Plot<numt>(){}
	virtual ~PlotPoints(){}
	PlotPoints &LineOnly(std::string name,Indexer&&points){
		Plot<numt>::OutputPlot(name,[&points](std::ofstream&data){
			for(PAIR p:points)
				data<<p.first<<" "<<p.second<<"\n";
		},"w l");
		return *this;
	}
	PlotPoints &WithoutErrors(std::string name,Indexer&&points){
		Plot<numt>::OutputPlot(name,[&points](std::ofstream&data){
			for(PAIR p:points)
				data<<p.first<<" "<<p.second<<"\n";
		},"using 1:2");
		return *this;
	}
	PlotPoints &WithErrorOnX(std::string name,Indexer&&points,FUNC error){
		Plot<numt>::OutputPlot(name,[&points,error](std::ofstream&data){
			for(auto p:points)
				data<<p.first<<" "<<p.second<<" "<<error(p.first)<<"\n";
		},"using 1:2:($2-$3):($2+$3) with yerrorbars");
		return *this;
	}
	PlotPoints &WithErrorOnY(std::string name,Indexer&&points,FUNC error){
		Plot<numt>::OutputPlot(name,[&points,error](std::ofstream&data){
			for(auto p:points)
				data<<p.first<<" "<<p.second<<" "<<error(p.second)<<"\n";
		},"using 1:2:($2-$3):($2+$3) with yerrorbars");
		return *this;
	}
};
#endif

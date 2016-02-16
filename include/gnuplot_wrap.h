// this file is distributed under 
// MIT license
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
#include <math_h/interpolate.h>
#include <math_h/sigma.h>
#include <math_h/hist.h>
namespace GnuplotWrap{
	using namespace std;
	using namespace MathTemplates;
	class Plotter{
	public:
		Plotter();
		~Plotter();
		static Plotter &Instance();
		void SetOutput(string&&out,string&&prefix="");
		string&OutPath()const;
		string&Prefix()const;
		string GetTerminal();
		string GetFileName();
		Plotter &operator<<(const string&line);
		Plotter &operator<<(string&&line);
	private:
		vector<string> lines;
		unsigned int terminal_counter,filename_counter;
		string outpath;
		string m_prefix;
	};
	template<class numt>class Plot{
	private:
		vector<string> lines;
		vector<string> plots;
	public:
		typedef function<numt(numt)> FUNC;
		typedef function<void(ofstream&)> PLOTOUTPUT;
		Plot&operator<<(const string&line){
			lines.push_back(line);
			return *this;
		}
		Plot&operator<<(string&&line){return operator<<(line);}
		Plot(){
			operator<<(Plotter::Instance().GetTerminal());
		}
		virtual ~Plot(){
			for(string&line:lines)
				Plotter::Instance()<<line;
			for(int i=0,n=plots.size();i<n;i++){
				string line=plots[i];
				if(i==0)
					line="plot "+line;
				if(i<(n-1))
					line+=",\\";
				Plotter::Instance()<<line;
			}
		}
	public:
		Plot& Object(string&&plot){
			plots.push_back(plot);
			return *this;
		}
		Plot& File(const string&name,const string&options,const string&title){
			string line="\""+name+"\" ";
			line+=options;
			line+=" title \"";
			line+=title;
			line+="\"";
			return Object(static_cast<std::string&&>(line));
		}
		Plot& File(string&&name,string&&options,string&&title=""){
			return File(name,options,title);
		}
		Plot&OutputPlot(PLOTOUTPUT delegate,string&&options,const string&title){
			ofstream data;
			string filename=Plotter::Instance().GetFileName();
			data.open((Plotter::Instance().OutPath()+"/"+filename).c_str());
			if(data.is_open()){
				delegate(data);
				File(filename,options,title);
				data.close();
			}
			return *this;
		}
		Plot&OutputPlot(PLOTOUTPUT delegate,string&&options,string&&title=""){
			return OutputPlot(delegate,static_cast<string&&>(options),title);
		}
		Plot &Func(FUNC func,numt from,numt to,numt step,string&&title=""){
			OutputPlot([func,from,to,step](ofstream&data){
				for(numt x=from;x<=to;x+=step)
					data<<x<<" "<<func(x)<<endl;
			},"w l",title);
			return *this;
		}
	};
	template<class numt,class Indexer=LinearInterpolation<numt>>
	class PlotPointsErrorless:public Plot<numt>{
	public:
		typedef pair<numt,numt> POINT;
		PlotPointsErrorless():Plot<numt>(){}
		virtual ~PlotPointsErrorless(){}
		PlotPointsErrorless &Line(const Indexer&points,string&&title=""){
			Plot<numt>::OutputPlot([&points](ofstream&data){
				for(POINT p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"w l",title);
			return *this;
		}
		PlotPointsErrorless &Points(const Indexer&points,string&&title=""){
			Plot<numt>::OutputPlot([&points](ofstream&data){
				for(POINT p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"using 1:2",title);
			return *this;
		}
	};
	template<class numt,class Indexer=hist<numt>>
	class PlotPoints:public Plot<numt>{
	public:
		PlotPoints(){}
		virtual ~PlotPoints(){}
		PlotPoints&Hist(const Indexer&data,const string&title){
			Plot<double>::OutputPlot([&data](ofstream&str){
				for(const point<numt> p:data)
					str<<p.X().val()<<" "<<p.Y().val()<<" "<<p.X().delta()<<" "<<p.Y().delta()<<endl;
			},"using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars",title);
			return *this;
		}
		PlotPoints&Hist(const Indexer&data,string&&title=""){return Hist(data,title);}
		PlotPoints&Hist(Indexer&&data,string&&title=""){return Hist(data,title);}
	};
};
#endif

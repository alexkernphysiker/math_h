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
		Plot& File(string&&name,string&&options,string&&title=""){
			string line="\""+name+"\" ";
			line+=options;
			line+=" title \"";
			line+=title;
			line+="\"";
			return Object(static_cast<std::string&&>(line));
		}
		Plot&OutputPlot(PLOTOUTPUT delegate,string&&options,string&&title=""){
			ofstream data;
			string filename=Plotter::Instance().GetFileName();
			data.open((Plotter::Instance().OutPath()+"/"+filename).c_str());
			if(data.is_open()){
				delegate(data);
				File(
					static_cast<std::string&&>(filename),
					static_cast<std::string&&>(options),
					static_cast<std::string&&>(title)
				);
				data.close();
			}
			return *this;
		}
		Plot &Hist(const LinearInterpolation_fixedsize<numt>&points,string&&title=""){
			OutputPlot([this,&points](ofstream&data){
				for(int i=0,n=points.size();i<n;i++)
					data<<points.getX(i)<<" "<<points.getY(i)<<endl;
			},"using 1:2",static_cast<string&&>(title));
			return *this;
		}
		Plot &HistWithStdError(const LinearInterpolation_fixedsize<numt>&points,string&&title=""){
			OutputPlot([&points](ofstream&data){
				for(int i=0,n=points.size();i<n;i++){
					double y=points.getY(i);
					if(y<0)throw std::exception();
				   double dy=sqrt(y);if(y<1)y=1;
				   data<<points.getX(i)<<" "<<y<<" "<<dy<<endl;
				}
			},"using 1:2:($2-$3):($2+$3) with yerrorbars",static_cast<string&&>(title));
			return *this;
		}
		Plot &Func(FUNC func,numt from,numt to,numt step,string&&title=""){
			OutputPlot([func,from,to,step](ofstream&data){
				for(double x=from;x<=to;x+=step)
					data<<x<<" "<<func(x)<<endl;
			},"w l",static_cast<string&&>(title));
			return *this;
		}
	};
	template<class numt,class Indexer>class PlotPoints:public Plot<numt>{
	public:
		typedef std::pair<numt,numt> POINT;
		PlotPoints():Plot<numt>(){}
		virtual ~PlotPoints(){}
		PlotPoints &Line(const Indexer&points,string&&title=""){
			Plot<numt>::OutputPlot([&points](ofstream&data){
				for(POINT p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"w l",static_cast<string&&>(title));
			return *this;
		}
		PlotPoints &Points(const Indexer&points,string&&title=""){
			Plot<numt>::OutputPlot([&points](ofstream&data){
				for(POINT p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"using 1:2",static_cast<string&&>(title));
			return *this;
		}
	};
	template<class numt,class Indexer>class PlotValues:public Plot<numt>{
	public:
		typedef std::pair<value<double>,value<double>> POINT;
		PlotValues():Plot<numt>(){}
		virtual ~PlotValues(){}
		PlotValues &Points(const Indexer&points,string&&title=""){
			Plot<numt>::OutputPlot([&points](ofstream&data){
				for(POINT p:points)
					data<<p.first.val()<<" "<<p.first.delta()<<" "<<p.second.val()<<" "<<p.second.delta()<<endl;
			},"using 1:3:($1-$2):($1+$2):($3-$4):($3+$4) with xyerrorbars",static_cast<string&&>(title));
			return *this;
		}
	};
};
#endif

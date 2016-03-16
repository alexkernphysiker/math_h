// this file is distributed under 
// MIT license
#ifndef VIJVUSSC
#define VIJVUSSC
#include <vector>
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
		static Plotter&Instance();
		void SetOutput(const string&&out,const string&&prefix="");
		const string&OutPath()const;
		const string&Prefix()const;
		const string GetTerminal();
		const string GetFileName();
		Plotter &operator<<(const string&line);
		Plotter &operator<<(const string&&line);
	private:
		vector<string> lines;
		unsigned int terminal_counter,filename_counter;
		string outpath;
		string m_prefix;
	};
	template<class numtX,class numtY=numtX>class Plot{
	private:
		vector<string> lines;
		vector<string> plots;
	public:
		typedef function<numtY(numtX)> FUNC;
		typedef function<void(ofstream&)> PLOTOUTPUT;
		Plot&operator<<(const string&line){
			lines.push_back(line);
			return *this;
		}
		Plot&operator<<(const string&&line){return operator<<(line);}
		Plot(){
			operator<<(Plotter::Instance().GetTerminal());
			operator<<("unset pm3d");
		}
		virtual ~Plot(){
			for(const string&line:lines)
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
		Plot& Object(const string&plot){
			plots.push_back(plot);
			return *this;
		}
		Plot& Object(const string&&plot){
			return Object(plot);
		}
		Plot& File(const string&name,const string&options,const string&title){
			string line="\""+name+"\" ";
			line+=options;
			line+=" title \"";
			line+=title;
			line+="\"";
			return Object(line);
		}
		Plot& File(const string&&name,const string&&options,const string&&title=""){
			return File(name,options,title);
		}
		Plot&OutputPlot(const PLOTOUTPUT delegate,const string&options,const string&title){
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
		Plot&OutputPlot(const PLOTOUTPUT delegate,const string&&options,const string&title){
			return OutputPlot(delegate,options,title);
		}
		Plot&OutputPlot(const PLOTOUTPUT delegate,const string&&options,const string&&title=""){
			return OutputPlot(delegate,options,title);
		}
		Plot&Line(const vector<pair<numtX,numtY>>&points,const string&title){
			Plot<numtX,numtY>::OutputPlot([&points](ofstream&data){
				for(const auto&p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"w l",title);
			return *this;
		}
		Plot&Line(const vector<pair<numtX,numtY>>&points,const string&&title=""){
			return Line(points,title);
		}
		Plot&Line(const vector<pair<numtX,numtY>>&&points,const string&&title=""){
			return Line(points,title);
		}
		Plot&Line(const initializer_list<pair<numtX,numtY>>&&points,const string&&title=""){
			Plot<numtX,numtY>::OutputPlot([&points](ofstream&data){
				for(const auto&p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"w l",title);
			return *this;
		}
		Plot&Line(const SortedPoints<numtX,numtY>&points,const string&title){
			Plot<numtX,numtY>::OutputPlot([&points](ofstream&data){
				for(const auto&p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"w l",title);
			return *this;
		}
		Plot&Line(const SortedPoints<numtX,numtY>&points,const string&&title=""){
			return Line(points,title);
		}
		Plot&Line(const SortedPoints<numtX,numtY>&&points,const string&&title=""){
			return Line(points,title);
		}
		Plot&Points(const vector<pair<numtX,numtY>>&points,const string&title){
			Plot<numtX,numtY>::OutputPlot([&points](ofstream&data){
				for(const auto&p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"using 1:2",title);
			return *this;
		}
		Plot&Points(const vector<pair<numtX,numtY>>&points,const string&&title=""){
			return Points(points,title);
		}
		Plot&Points(const vector<pair<numtX,numtY>>&&points,const string&&title=""){
			return Points(points,title);
		}
		Plot&Points(const initializer_list<pair<numtX,numtY>>&&points,const string&&title=""){
			Plot<numtX,numtY>::OutputPlot([&points](ofstream&data){
				for(const auto&p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"using 1:2",title);
			return *this;
		}
		Plot&Points(const SortedPoints<numtX,numtY>&points,const string&title){
			Plot<numtX,numtY>::OutputPlot([&points](ofstream&data){
				for(const auto&p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"using 1:2",title);
			return *this;
		}
		Plot&Points(const SortedPoints<numtX,numtY>&points,const string&&title=""){
			return Points(points,title);
		}
		Plot&Points(const SortedPoints<numtX,numtY>&&points,const string&&title=""){
			return Points(points,title);
		}
		Plot&Hist(const hist<numtX,numtY>&data,const string&title){
			Plot<numtX,numtY>::OutputPlot([&data](ofstream&str){
				for(const point<numtX,numtY> p:data)
					str<<p.X().val()<<" "<<p.Y().val()<<" "<<p.X().delta()<<" "<<p.Y().delta()<<endl;
			},"using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars",title);
			return *this;
		}
		Plot&Hist(const hist<numtX,numtY>&data,const string&&title=""){return Hist(data,title);}
		Plot&Hist(const hist<numtX,numtY>&&data,const string&&title=""){return Hist(data,title);}
	};
	
	enum TypeOf3D{normal,sp2};
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class PlotHist2d{
	private:
		vector<string> lines;
		vector<string> plots;
		string Distr2File(const hist2d<numtX,numtY,numtZ>&D)const{
			string filename=Plotter::Instance().GetFileName();
			ofstream data;
			data.open((Plotter::Instance().OutPath()+"/"+filename).c_str());
			if(data.is_open()){
				data<<D.X().size()<<" ";
				for(const auto&x:D.X())
					data<<x.val()<<" ";
				for(size_t i=0,I=D.Y().size();i<I;i++){
					data<<endl<<D.Y()[i].val();
					for(size_t j=0,J=D.X().size();j<J;j++)
						data<<" "<<D[j][i].val();
				}
				data.close();
			}
			return filename;
		}
		string Points2File(const vector<point3d<numtX,numtY,numtZ>>&points){
			string filename=Plotter::Instance().GetFileName();
			ofstream data;
			data.open((Plotter::Instance().OutPath()+"/"+filename).c_str());
			if(data.is_open()){
				for(const point3d<numtX,numtY,numtZ>&p:points)
					data<<p.X().val()<<" "<<p.Y().val()<<" "<<p.Z().val()<<endl;
				data.close();
			}
			return filename;
		}
	public:
		PlotHist2d&operator<<(const string&line){
			lines.push_back(line);
			return *this;
		}
		PlotHist2d&operator<<(const string&&line){return operator<<(line);}
		PlotHist2d(const TypeOf3D type){
			operator<<(Plotter::Instance().GetTerminal());
			if(sp2==type){
				operator<<("unset key");
				operator<<("unset surface");
				operator<<("set view map");
				operator<<("set pm3d at b");
			}else{
				operator<<("unset key");
				operator<<("unset surface");
				operator<<("unset view");
				operator<<("set pm3d");
			}
		}
		virtual ~PlotHist2d(){
			for(const string&line:lines)Plotter::Instance()<<line;
			for(int i=0,n=plots.size();i<n;i++){
				string line=plots[i];
				if(i==0)
					line="splot "+line;
				if(i<(n-1))
					line+=",\\";
				Plotter::Instance()<<line;
			}
		}
		PlotHist2d& Object(const string&&plot){
			plots.push_back(plot);
			return *this;
		}
		PlotHist2d&Distr(const hist2d<numtX,numtY,numtZ>&D,const string&title){
			return Object(string("'")+Distr2File(D)+"' matrix nonuniform title'"+title+"'");
		}
		PlotHist2d&Distr(const hist2d<numtX,numtY,numtZ>&D,const string&&title=""){return Distr(D,title);}
		PlotHist2d&Distr(const hist2d<numtX,numtY,numtZ>&&D,const string&&title=""){return Distr(D,title);}
		
		PlotHist2d&Points(const vector<point3d<numtX,numtY,numtZ>>&points,const string&title=""){
			return Object(string("'")+Points2File(points)+"' u 1:2:3 w points title'"+title+"'");
		}
		PlotHist2d&Points(const vector<point3d<numtX,numtY,numtZ>>&points,const string&&title=""){return Points(points,title);}
		PlotHist2d&Points(const vector<point3d<numtX,numtY,numtZ>>&&points,const string&&title=""){return Points(points,title);}
		PlotHist2d&Points(const initializer_list<point3d<numtX,numtY,numtZ>>&points,const string&&title=""){return Points(vector<point3d<numtX,numtY,numtZ>>(points),title);}
		
		PlotHist2d&Line(const vector<point3d<numtX,numtY,numtZ>>&points,const string&title){
			return Object(string("'")+Points2File(points)+"' u 1:2:3 w line title'"+title+"'");
		}
		PlotHist2d&Line(const vector<point3d<numtX,numtY,numtZ>>&points,const string&&title=""){return Line(points,title);}
		PlotHist2d&Line(const vector<point3d<numtX,numtY,numtZ>>&&points,const string&&title=""){return Line(points,title);}
		PlotHist2d&Line(const initializer_list<point3d<numtX,numtY,numtZ>>&&points,const string&&title=""){return Line(vector<point3d<numtX,numtY,numtZ>>(points),title);}
	};
};
#endif

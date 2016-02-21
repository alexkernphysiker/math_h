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
		Plot&operator<<(string&&line){return operator<<(line);}
		Plot(){
			operator<<(Plotter::Instance().GetTerminal());
			operator<<("unset pm3d");
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
		Plot&Line(const vector<pair<numtX,numtY>>&points,const string&title){
			Plot<numtX,numtY>::OutputPlot([&points](ofstream&data){
				for(const auto&p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"w l",title);
			return *this;
		}
		Plot&Line(const vector<pair<numtX,numtY>>&points,string&&title=""){
			return Line(points,title);
		}
		Plot&Line(vector<pair<numtX,numtY>>&&points,string&&title=""){
			return Line(points,title);
		}
		Plot&Line(initializer_list<pair<numtX,numtY>>&&points,string&&title=""){
			Plot<numtX,numtY>::OutputPlot([&points](ofstream&data){
				for(const auto&p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"w l",title);
			return *this;
		}
		Plot&Line(const LinearInterpolation<numtX,numtY>&points,string&&title=""){
			Plot<numtX,numtY>::OutputPlot([&points](ofstream&data){
				for(const auto&p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"w l",title);
			return *this;
		}
		Plot&Line(LinearInterpolation<numtX,numtY>&&points,string&&title=""){
			return Line(points,static_cast<string&&>(title));
		}
		Plot&Points(const vector<pair<numtX,numtY>>&points,const string&title){
			Plot<numtX,numtY>::OutputPlot([&points](ofstream&data){
				for(const auto&p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"using 1:2",title);
			return *this;
		}
		Plot&Points(const vector<pair<numtX,numtY>>&points,string&&title=""){
			return Points(points,title);
		}
		Plot&Points(vector<pair<numtX,numtY>>&&points,string&&title=""){
			return Points(points,title);
		}
		Plot&Points(initializer_list<pair<numtX,numtY>>&&points,string&&title=""){
			Plot<numtX,numtY>::OutputPlot([&points](ofstream&data){
				for(const auto&p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"using 1:2",title);
			return *this;
		}
		Plot&Points(const LinearInterpolation<numtX,numtY>&points,const string&title){
			Plot<numtX,numtY>::OutputPlot([&points](ofstream&data){
				for(const auto&p:points)
					data<<p.first<<" "<<p.second<<endl;
			},"using 1:2",title);
			return *this;
		}
		Plot&Points(const LinearInterpolation<numtX,numtY>&points,string&&title=""){
			return Points(points,title);
		}
		Plot&Points(LinearInterpolation<numtX,numtY>&&points,string&&title=""){
			return Points(points,title);
		}
		Plot&Hist(const hist<numtX,numtY>&data,const string&title){
			Plot<numtX,numtY>::OutputPlot([&data](ofstream&str){
				for(const point<numtX,numtY> p:data)
					str<<p.X().val()<<" "<<p.Y().val()<<" "<<p.X().delta()<<" "<<p.Y().delta()<<endl;
			},"using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars",title);
			return *this;
		}
		Plot&Hist(const hist<numtX,numtY>&data,string&&title=""){return Hist(data,title);}
		Plot&Hist(hist<numtX,numtY>&&data,string&&title=""){return Hist(data,title);}
	};
	
	enum TypeOf3D{normal,sp2};
	template<class numtX,class numtY=numtX,class numtZ=numtY>
	class PlotHist2d{
	public:
		class point{
		private:
			numtX _x;
			numtY _y;
			numtZ _z;
		public:
			point(numtX x,numtY y, numtZ z){
				_x=x;_y=y;_z=z;
			}
			point(const point&source):point(source._x,source._y,source._z){}
			~point(){}
			numtX x()const{return _x;}
			numtY y()const{return _y;}
			numtZ z()const{return _z;}
		};
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
		string Points2File(const vector<point>&points){
			string filename=Plotter::Instance().GetFileName();
			ofstream data;
			data.open((Plotter::Instance().OutPath()+"/"+filename).c_str());
			if(data.is_open()){
				for(const point&p:points)
					data<<p.x()<<" "<<p.y()<<" "<<p.z()<<endl;
				data.close();
			}
			return filename;
		}
	public:
		PlotHist2d&operator<<(const string&line){
			lines.push_back(line);
			return *this;
		}
		PlotHist2d&operator<<(string&&line){return operator<<(line);}
		PlotHist2d(TypeOf3D type){
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
			for(string&line:lines)Plotter::Instance()<<line;
			for(int i=0,n=plots.size();i<n;i++){
				string line=plots[i];
				if(i==0)
					line="splot "+line;
				if(i<(n-1))
					line+=",\\";
				Plotter::Instance()<<line;
			}
		}
		PlotHist2d& Object(string&&plot){
			plots.push_back(plot);
			return *this;
		}
		PlotHist2d&Distr(const hist2d<numtX,numtY,numtZ>&D,string&&title=""){
			return Object(string("'")+Distr2File(D)+"' matrix nonuniform title'"+title+"'");
		}
		PlotHist2d&Distr(hist2d<numtX,numtY,numtZ>&&D,string&&title=""){
			return Object(string("'")+Distr2File(D)+"' matrix nonuniform title'"+title+"'");
		}
		PlotHist2d&Points(const vector<point>&points,string&&title=""){
			return Object(string("'")+Points2File(points)+"' u 1:2:3 w points title'"+title+"'");
		}
		PlotHist2d&Points(const vector<point>&&points,string&&title=""){
			return Points(points,static_cast<string&&>(title));
		}
		PlotHist2d&Points(initializer_list<point>&&points,string&&title=""){
			return Points(vector<point>(points),static_cast<string&&>(title));
		}
		PlotHist2d&Line(const vector<point>&points,string&&title=""){
			return Object(string("'")+Points2File(points)+"' u 1:2:3 w line title'"+title+"'");
		}
		PlotHist2d&Line(const vector<point>&&points,string&&title=""){
			return Line(points,static_cast<string&&>(title));
		}
		PlotHist2d&Line(initializer_list<point>&&points,string&&title=""){
			return Line(vector<point>(points),static_cast<string&&>(title));
		}
	};
};
#endif

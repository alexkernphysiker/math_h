// this file is distributed under 
// MIT license
#ifndef VIJVUSSC
#define VIJVUSSC
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <utility>
#include <sstream>
#include <memory>
#include <unistd.h>
#include <math.h>
#include <exception>
#include <math_h/hists.h>
namespace GnuplotWrap{
	class Plotter{
	public:
		Plotter();
		~Plotter();
		static Plotter&Instance();
		void SetOutput(const std::string&out,const std::string&prefix="");
		const std::string GetTerminal(const std::string&name="");
		std::pair<const std::string,std::ofstream> File(const std::string&name="");
		std::ifstream GetInput(const std::string&name);
		template<class numt=double>
		const MathTemplates::SortedPoints<numt> GetPoints2(const std::string&name=""){
		    MathTemplates::SortedPoints<numt> res;
		    auto str=GetInput(name);
		    numt x,y;
		    while(str){
			str>>x>>y;
			res<<MathTemplates::point<numt>(x,y);
		    }
		    return res;
		}
		template<class numt=double>
		const MathTemplates::SortedPoints<MathTemplates::value<numt>> GetPoints4(const std::string&name=""){
		    MathTemplates::SortedPoints<MathTemplates::value<numt>> res;
		    auto str=GetInput(name);
		    numt x,y,dx,dy;
		    while(str){
			str>>x>>y>>dx>>dy;
			res<<MathTemplates::point<MathTemplates::value<numt>>({x,dx},{y,dy});
		    }
		    return res;
		}
		Plotter &operator<<(const std::string&line);
	private:
		const std::string GetFileName(const std::string&name);
		std::vector<std::string> lines;
		unsigned int terminal_counter,filename_counter;
		std::string outpath;
		std::string m_prefix;
	};
	template<class numtX=double,class numtY=numtX>class Plot{
	private:
		std::vector<std::string> lines;
		std::vector<std::string> plots;
	public:
		typedef std::function<numtY(numtX)> FUNC;
		typedef std::function<void(std::ofstream&)> PLOTOUTPUT;
		Plot&operator<<(const std::string&line){
			lines.push_back(line);
			return *this;
		}
		Plot(const std::string&name=""){
			operator<<(Plotter::Instance().GetTerminal(name));
			operator<<("unset pm3d");
			operator<<("unset title");
			operator<<("unset key");
			operator<<("unset surface");
			operator<<("unset view");
			operator<<("unset xrange");
			operator<<("unset yrange");
			operator<<("unset xlabel");
			operator<<("unset ylabel");
		}
		virtual ~Plot(){
			for(const std::string&line:lines)
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
		Plot& Object(const std::string&plot){
			plots.push_back(plot);
			return *this;
		}
		Plot& File(const std::string&name,const std::string&options,const std::string&title){
			std::string line="\""+name+"\" ";
			line+=options;
			line+=" title \"";
			line+=title;
			line+="\"";
			return Object(line);
		}
		Plot&OutputPlot(const PLOTOUTPUT delegate,const std::string&options,const std::string&title="",const std::string&name=""){
			auto data=Plotter::Instance().File(name);
			if(data.second){
				delegate(data.second);
				File(data.first,options,title);
				data.second.close();
			}
			return *this;
		}
		Plot&Line(const std::vector<MathTemplates::point<numtX,numtY>>&points,const std::string&title="",const std::string&name=""){
			Plot<numtX,numtY>::OutputPlot([&points](std::ofstream&data){
				for(const auto&p:points)
					data<<p.X()<<" "<<p.Y()<<std::endl;
			},"w l",title,name);
			return *this;
		}
		Plot&Line(const MathTemplates::SortedPoints<numtX,numtY>&points,const std::string&title="",const std::string&name=""){
			Plot<numtX,numtY>::OutputPlot([&points](std::ofstream&data){
				for(const auto&p:points)
					data<<p.X()<<" "<<p.Y()<<std::endl;
			},"w l",title,name);
			return *this;
		}
		Plot&Points(const std::vector<MathTemplates::point<numtX,numtY>>&points,const std::string&title="",const std::string&name=""){
			Plot<numtX,numtY>::OutputPlot([&points](std::ofstream&data){
				for(const auto&p:points)
					data<<p.X()<<" "<<p.Y()<<std::endl;
			},"using 1:2",title,name);
			return *this;
		}
		Plot&Points(const MathTemplates::SortedPoints<numtX,numtY>&points,const std::string&title="",const std::string&name=""){
			Plot<numtX,numtY>::OutputPlot([&points](std::ofstream&data){
				for(const auto&p:points)
					data<<p.X()<<" "<<p.Y()<<std::endl;
			},"using 1:2",title,name);
			return *this;
		}
		Plot&Hist(const MathTemplates::SortedPoints<MathTemplates::value<numtX>,MathTemplates::value<numtY>>&data,
			  const std::string&title="",const std::string&name=""){
			Plot<numtX,numtY>::OutputPlot([&data](std::ofstream&str){
				for(const auto p:data)
					str<<p.X().val()<<" "<<p.Y().val()<<" "<<p.X().uncertainty()<<" "<<p.Y().uncertainty()<<std::endl;
			},"using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars",title,name);
			return *this;
		}
	};
	
	enum TypeOf3D{normal,sp2};
	template<class numtX=double,class numtY=numtX,class numtZ=numtY>
	class PlotHist2d{
	private:
		std::vector<std::string> lines;
		std::vector<std::string> plots;
		std::string Surf2File(const MathTemplates::BiSortedPoints<numtX,numtY,numtZ>&D,const std::string&name="")const{
			auto data=Plotter::Instance().File(name);
			if(data.second){
				data.second<<D.X().size()<<" ";
				for(const auto&x:D.X())
					data.second<<x<<" ";
				for(size_t i=0,I=D.Y().size();i<I;i++){
					data.second<<std::endl<<D.Y()[i];
					for(size_t j=0,J=D.X().size();j<J;j++)
						data.second<<" "<<D[j][i];
				}
				data.second.close();
			}
			return data.first;
		}
		std::string Points2File(const std::vector<MathTemplates::point3d<numtX,numtY,numtZ>>&points,const std::string&name=""){
			auto data=Plotter::Instance().File(name);
			if(data.second){
				for(const auto&p:points)
					data.second<<p.X()<<" "<<p.Y()<<" "<<p.Z()<<std::endl;
				data.second.close();
			}
			return data.first;
		}
		std::string Distr2File(const MathTemplates::hist2d<numtX,numtY,numtZ>&D,const std::string&name="")const{
			auto data=Plotter::Instance().File(name);
			if(data.second){
				data.second<<D.X().size()<<" ";
				for(const auto&x:D.X())
					data.second<<x.val()<<" ";
				for(size_t i=0,I=D.Y().size();i<I;i++){
					data.second<<std::endl<<D.Y()[i].val();
					for(size_t j=0,J=D.X().size();j<J;j++)
						data.second<<" "<<D[j][i].val();
				}
				data.second.close();
			}
			return data.first;
		}
		std::string Points2File(const std::vector<
		MathTemplates::point3d<MathTemplates::value<numtX>,
		MathTemplates::value<numtY>,MathTemplates::value<numtZ>>>&points,
		const std::string&name=""){
			auto data=Plotter::Instance().File(name);
			if(data.second){
				for(const auto&p:points)
					data.second<<p.X().val()<<" "<<p.Y().val()<<" "<<p.Z().val()<<std::endl;
				data.second.close();
			}
			return data.first;
		}
	public:
		PlotHist2d&operator<<(const std::string&line){
			lines.push_back(line);
			return *this;
		}
		PlotHist2d(const TypeOf3D type,const std::string&name=""){
			operator<<(Plotter::Instance().GetTerminal(name));
			operator<<("unset title");
			operator<<("unset key");
			operator<<("unset surface");
			if(sp2==type){
				operator<<("set view map");
				operator<<("set pm3d at b");
			}else{
				operator<<("unset view");
				operator<<("set pm3d");
			}
			operator<<("unset xrange");
			operator<<("unset yrange");
			operator<<("unset zrange");
			operator<<("unset xlabel");
			operator<<("unset ylabel");
			operator<<("unset zlabel");
		}
		virtual ~PlotHist2d(){
			for(const std::string&line:lines)Plotter::Instance()<<line;
			for(int i=0,n=plots.size();i<n;i++){
				std::string line=plots[i];
				if(i==0)
					line="splot "+line;
				if(i<(n-1))
					line+=",\\";
				Plotter::Instance()<<line;
			}
		}
		PlotHist2d& Object(const std::string&&plot){
			plots.push_back(plot);
			return *this;
		}
		PlotHist2d&Surface(const MathTemplates::BiSortedPoints<numtX,numtY,numtZ>&D,const std::string&title=""){
			return Object(std::string("'")+Surf2File(D)+"' matrix nonuniform title'"+title+"'");
		}
		
		PlotHist2d&Distr(const MathTemplates::hist2d<numtX,numtY,numtZ>&D,const std::string&title=""){
			return Object(std::string("'")+Distr2File(D)+"' matrix nonuniform title'"+title+"'");
		}
		PlotHist2d&Points(const std::vector<MathTemplates::point3d<
		MathTemplates::value<numtX>,MathTemplates::value<numtY>,
		MathTemplates::value<numtZ>>>&points,const std::string&title=""){
			return Object(std::string("'")+Points2File(points)+"' u 1:2:3 w points title'"+title+"'");
		}
		PlotHist2d&Points(const std::initializer_list<
		MathTemplates::point3d<MathTemplates::value<numtX>,
		MathTemplates::value<numtY>,MathTemplates::value<numtZ>>>&points,const std::string&title=""){
			return Points(std::vector<MathTemplates::point3d<
			MathTemplates::value<numtX>,MathTemplates::value<numtY>,
			MathTemplates::value<numtZ>>>(points),title);
		}
		
		PlotHist2d&Line(const std::vector<MathTemplates::point3d<
		MathTemplates::value<numtX>,MathTemplates::value<numtY>,
		MathTemplates::value<numtZ>>>&points,const std::string&title=""){
			return Object(std::string("'")+Points2File(points)+"' u 1:2:3 w line title'"+title+"'");
		}

		PlotHist2d&Points(const std::vector<MathTemplates::point3d<numtX,numtY,numtZ>>&points,const std::string&title=""){
			return Object(std::string("'")+Points2File(points)+"' u 1:2:3 w points title'"+title+"'");
		}
		PlotHist2d&Points(const std::initializer_list<MathTemplates::point3d<numtX,numtY,numtZ>>&points,const std::string&title=""){
			return Points(std::vector<MathTemplates::point3d<numtX,numtY,numtZ>>(points),title);
		}
		
		PlotHist2d&Line(const std::vector<MathTemplates::point3d<numtX,numtY,numtZ>>&points,const std::string&title=""){
			return Object(std::string("'")+Points2File(points)+"' u 1:2:3 w line title'"+title+"'");
		}
		
	};
	template<class numt=double>
	class PlotDistr1D:public MathTemplates::Distribution1D<numt>{
	private:
	    std::string m_title,m_axis;
	public:
	    PlotDistr1D(
		    const std::string&title,const std::string&axis,
		    const MathTemplates::SortedChain<MathTemplates::value<numt>>&data
	    ):MathTemplates::Distribution1D<numt>(data),m_title(title),m_axis(axis){}
	    virtual ~PlotDistr1D(){
		Plot<numt>().Hist(*this)<<"set title '"+m_title+"'"<<"set yrange [0:]"
		<<"set xlabel '"+m_axis+"'"<<"set ylabel 'counts'";
	    }
	};
	template<class numt=double>
	class PlotDistr2D:public MathTemplates::Distribution2D<numt>{
	private:
	    std::string m_title;
	public:
	    PlotDistr2D(
		const std::string&title,
		const MathTemplates::SortedChain<MathTemplates::value<numt>>&X,
		const MathTemplates::SortedChain<MathTemplates::value<numt>>&Y
	    ):MathTemplates::Distribution2D<numt>(X,Y),m_title(title){}
	    virtual ~PlotDistr2D(){
		PlotHist2d<numt>(sp2).Distr(*this)<<"set title '"+m_title+"'";
	    }
	};
};
#endif

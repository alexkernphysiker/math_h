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
#include <math_h/error.h>
namespace GnuplotWrap
{
#define TDecl_init size_t width = 1024,size_t height=868
#define TDecl size_t width,size_t height
#define TUse width,height
template<TDecl_init> class PlotHist2d;
template<TDecl_init> class Plot;
template<TDecl_init> class Plotter
{
    friend class Plot<TUse>;
    friend class PlotHist2d<TUse>;
private:
    std::vector<std::string> lines;
    unsigned int terminal_counter, filename_counter;
    std::string outpath;
    std::string m_prefix;
    const std::string GetFileName(const std::string &name)
    {
        if (name == "") {
            filename_counter++;
            auto cnt = std::to_string(filename_counter);
            while (cnt.length() < 6)cnt = "0" + cnt;
            return "." + m_prefix + ".numeric-data."
	+std::to_string(width)+"."+std::to_string(height)+"." + cnt + ".txt";
        } else {
            return name + ".txt";
        }
    }
protected:
    Plotter &operator<<(const std::string &line)
    {
        lines.push_back(line);
        return *this;
    }
    const std::string File(std::ofstream&str,const std::string &name = "")
    {
        const auto n = GetFileName(name);
	str.open((outpath + "/" + n).c_str());
        return n;
    }
    std::ifstream GetInput(const std::string &name)
    {
        if (name == "")throw MathTemplates::Exception<Plotter>("Cannot get input without name");
        return std::ifstream((outpath + "/" + GetFileName(name)).c_str());
    }
    const std::string GetTerminal(const std::string &name = "")
    {
        const std::string firstline = "set terminal pngcairo size "
	+std::to_string(width)+","+std::to_string(height)+" font 'Verdana,18'\n";
        if (name == "") {
            terminal_counter++;
            auto cnt = std::to_string(terminal_counter);
            while (cnt.length() < 5)cnt = "0" + cnt;
            return firstline + "set output '"
	+std::to_string(width)+"-"+std::to_string(height)+"."+
	m_prefix + "-plot-" + cnt + ".png'";
        } else {
            return firstline + "set output '" + name + ".png'";
        }
    }
public:
    Plotter()
    {
        terminal_counter = 0;
        filename_counter = 0;
        outpath = ".";
    }
    ~Plotter()
    {
        std::string name = m_prefix + ".gnuplot-script";
        std::ofstream str(std::string(outpath + "/" + name).c_str(), std::ofstream::out);
        for (const auto &line : lines)
            str << line << "\n";
        str.close();
        name = "gnuplot " + name;
        const std::string old = getcwd(NULL, 0);
        chdir(outpath.c_str());
        system(name.c_str());
        chdir(old.c_str());
    }
    static Plotter &Instance()
    {
        static Plotter m_instance;
        return m_instance;
    }
    void SetOutput(const std::string &out, const std::string &prefix = "")
    {
        if (lines.size() > 0)
            throw MathTemplates::Exception<Plotter>("Attempt to reset plot output settings");
        outpath = out;
        m_prefix = prefix;
    }
    template<class numtX = double,class numtY = numtX>
    const MathTemplates::Points<numtX,numtY> GetPoints(const std::string &name)
    {
        MathTemplates::Points<numtX,numtY> res;
        auto str = GetInput(name);
        numtX x;numtY y;
        while (str >> x >> y) {
            res.push_back(make_point(x, y));
        }
        return res;
    }
    template<class numtX = double,class numtY = numtX>
    const std::string SavePoints(const MathTemplates::Points<numtX,numtY>&data,const std::string &name = "")
    {
	std::ofstream str;
        auto n = Plotter::Instance().File(str,name);
        if (str) {
            for(const auto&P:data)
		str<<P.X()<<"\t"<<P.Y()<<std::endl;
            str.close();
        }
        return n;
    }
};
template<TDecl> class Plot
{
private:
    std::vector<std::string> lines;
    std::vector<std::string> plots;
    Plot &Object(const std::string &plot)
    {
        plots.push_back(plot);
        return *this;
    }
public:
    Plot &operator<<(const std::string &line)
    {
        lines.push_back(line);
        return *this;
    }
    Plot(const std::string &name = "")
    {
        operator<<(Plotter<TUse>::Instance().GetTerminal(name));
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
    virtual ~Plot()
    {
        for (const std::string &line : lines)
            Plotter<TUse>::Instance() << line;
        for (int i = 0, n = plots.size(); i < n; i++) {
            std::string line = plots[i];
            if (i == 0)
                line = "plot " + line;
            if (i < (n - 1))
                line += ",\\";
            Plotter<TUse>::Instance() << line;
        }
    }
    Plot &File(const std::string &name, const std::string &options, const std::string &title)
    {
        std::string line = "\"" + name + "\" ";
        line += options;
        line += " title \"";
        line += title;
        line += "\"";
        return Object(line);
    }
    template<class numtX = double,class numtY = numtX>
    inline Plot &Output(const MathTemplates::Points<numtX,numtY>&data,const std::string &options, const std::string &title = "", const std::string &name = "")
    {
        File(Plotter<TUse>::Instance().SavePoints(data,name), options, title);
        return *this;
    }
    template<class numtX = double,class numtY = numtX>
    inline Plot &Line(const MathTemplates::Points<numtX,numtY> &points, const std::string &title = "", const std::string &name = "")
    {
        return Output(points,"w l",title,name);
    }
    template<class numtX = double,class numtY = numtX>
    inline Plot &Line(const MathTemplates::SortedPoints<numtX,numtY> &points, const std::string &title = "", const std::string &name = "")
    {
        return Output(points(),"w l",title,name);
    }
    template<class numtX = double,class numtY = numtX>
    inline Plot &Points(const MathTemplates::Points<numtX,numtY> &points, const std::string &title = "", const std::string &name = "")
    {
        return Output<numtX,numtY>(points,"using 1:2",title,name);
    }
    template<class numtX = double,class numtY = numtX>
    inline Plot &Points(const MathTemplates::SortedPoints<numtX,numtY> &points, const std::string &title = "", const std::string &name = "")
    {
        return Output<numtX,numtY>(points(),"using 1:2",title,name);
    }
    template<class numtX = double,class numtY = numtX>
    inline Plot &Hist(const MathTemplates::SortedPoints<MathTemplates::value<numtX>, MathTemplates::value<numtY>> &data,
               const std::string &title = "", const std::string &name = "")
    {
        return Output(data(),"using 1:3:($1-$2):($1+$2):($3-$4):($3+$4) with xyerrorbars", title, name);
    }
    template<class numtX = double,class numtY = numtX>
    inline Plot &XYUncertainties(const MathTemplates::Points<MathTemplates::value<numtX>, MathTemplates::value<numtY>> &data,
               const std::string &title = "", const std::string &name = "")
    {
        return Output(data,"using 1:3:($1-$2):($1+$2):($3-$4):($3+$4) with xyerrorbars", title, name);
    }
    template<class numtX = double,class numtY = numtX>
    inline Plot &YUncertainties(const MathTemplates::Points<numtX, MathTemplates::value<numtY>> &data,
               const std::string &title = "", const std::string &name = "")
    {
        return Output(data,"using 1:2:($2-$3):($2+$3) with yerrorbars", title, name);
    }
    template<class numtX = double,class numtY = numtX>
    inline Plot &XUncertainties(const MathTemplates::Points<MathTemplates::value<numtX>,numtY> &data,
               const std::string &title = "", const std::string &name = "")
    {
        return Output(data,"using 1:3:($1-$2):($1+$2) with xerrorbars", title, name);
    }

};

enum TypeOf3D {normal, sp2};
template<TDecl>class PlotHist2d
{
private:
    std::vector<std::string> lines;
    std::vector<std::string> plots;
    template<class numtX = double,class numtY = numtX,class numtZ = numtY>
    std::string Surf2File(const MathTemplates::BiSortedPoints<numtX,numtY,numtZ> &D, const std::string &name = "")const
    {
        std::ofstream str;
        auto n = Plotter<TUse>::Instance().File(str,name);
        if (str) {
            str << D.X().size() << " ";
            for (const auto &x : D.X())
                str << x << " ";
            for (size_t i = 0, I = D.Y().size(); i < I; i++) {
                str << std::endl << D.Y()[i];
                for (size_t j = 0, J = D.X().size(); j < J; j++)
                    str << " " << D[j][i];
            }
            str.close();
        }
        return n;
    }
    template<class numtX = double,class numtY = numtX,class numtZ = numtY>
    std::string Points2File(const std::vector<MathTemplates::point3d<numtX,numtY,numtZ>> &points, const std::string &name = "")
    {
        std::ofstream str;
        auto n = Plotter<TUse>::Instance().File(str,name);
        if (str) {
            for (const auto &p : points)
                str << p.X() << " " << p.Y() << " " << p.Z() << std::endl;
            str.close();
        }
        return n;
    }
    template<class numtX = double,class numtY = numtX,class numtZ = numtY>
    std::string Distr2File(const MathTemplates::hist2d<numtX,numtY,numtZ> &D, const std::string &name = "")const
    {
        std::ofstream str;
        auto n = Plotter<TUse>::Instance().File(str,name);
        if (str) {
            str << D.X().size() << " ";
            for (const auto &x : D.X())
                str << x.val() << " ";
            for (size_t i = 0, I = D.Y().size(); i < I; i++) {
                str << std::endl << D.Y()[i].val();
                for (size_t j = 0, J = D.X().size(); j < J; j++)
                    str << " " << D[j][i].val();
            }
            str.close();
        }
        return n;
    }
    template<class numtX = double,class numtY = numtX,class numtZ = numtY>
    std::string Points2File(const std::vector <
                            MathTemplates::point3d<MathTemplates::value<numtX>,
                            MathTemplates::value<numtY>, MathTemplates::value<numtZ> >> &points,
                            const std::string &name = "")
    {
        std::ofstream str;
        auto n = Plotter<TUse>::Instance().File(str,name);
        if (str) {
            for (const auto &p : points)
                str << p.X().val() << " " << p.Y().val() << " " << p.Z().val() << std::endl;
            str.close();
        }
        return n;
    }
public:
    PlotHist2d &operator<<(const std::string &line)
    {
        lines.push_back(line);
        return *this;
    }
    PlotHist2d(const TypeOf3D type, const std::string &name = "")
    {
        operator<<(Plotter<TUse>::Instance().GetTerminal(name));
        operator<<("unset title");
        operator<<("unset key");
        operator<<("unset surface");
        if (sp2 == type) {
            operator<<("set view map");
            operator<<("set pm3d at b");
        } else {
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
    virtual ~PlotHist2d()
    {
        for (const std::string &line : lines)Plotter<TUse>::Instance() << line;
        for (int i = 0, n = plots.size(); i < n; i++) {
            std::string line = plots[i];
            if (i == 0)
                line = "splot " + line;
            if (i < (n - 1))
                line += ",\\";
            Plotter<TUse>::Instance() << line;
        }
    }
    PlotHist2d &Object(const std::string&plot)
    {
        plots.push_back(plot);
        return *this;
    }
    template<class numtX = double,class numtY = numtX,class numtZ = numtY>
    inline PlotHist2d &Surface(const MathTemplates::BiSortedPoints<numtX,numtY,numtZ> &D, const std::string &title = "")
    {
        return Object(std::string("'") + Surf2File(D) + "' matrix nonuniform title'" + title + "'");
    }
    template<class numtX = double,class numtY = numtX,class numtZ = numtY>
    inline PlotHist2d &Distr(const MathTemplates::hist2d<numtX,numtY,numtZ> &D, const std::string &title = "")
    {
        return Object(std::string("'") + Distr2File(D) + "' matrix nonuniform title'" + title + "'");
    }
    template<class numtX = double,class numtY = numtX,class numtZ = numtY>
    inline PlotHist2d &Points(const std::vector<MathTemplates::point3d<
                       MathTemplates::value<numtX>, MathTemplates::value<numtY>,
                       MathTemplates::value<numtZ>>> &points, const std::string &title = "")
    {
        return Object(std::string("'") + Points2File(points) + "' u 1:2:3 w points title'" + title + "'");
    }
    template<class numtX = double,class numtY = numtX,class numtZ = numtY>
    inline PlotHist2d &Line(const std::vector<MathTemplates::point3d<
                     MathTemplates::value<numtX>, MathTemplates::value<numtY>,
                     MathTemplates::value<numtZ>>> &points, const std::string &title = "")
    {
        return Object(std::string("'") + Points2File(points) + "' u 1:2:3 w line title'" + title + "'");
    }
    template<class numtX = double,class numtY = numtX,class numtZ = numtY>
    inline PlotHist2d &Points(const std::vector<MathTemplates::point3d<numtX,numtY,numtZ>> &points, const std::string &title = "")
    {
        return Object(std::string("'") + Points2File(points) + "' u 1:2:3 w points title'" + title + "'");
    }
    template<class numtX = double,class numtY = numtX,class numtZ = numtY>
    inline PlotHist2d &Line(const std::vector<MathTemplates::point3d<numtX,numtY,numtZ>> &points, const std::string &title = "")
    {
        return Object(std::string("'") + Points2File(points) + "' u 1:2:3 w line title'" + title + "'");
    }

};
template<class numtX = double,class numtY = numtX,TDecl_init>
class PlotDistr1D: public MathTemplates::Distribution1D<numtX,numtY>
{
private:
    std::string m_title, m_axis,m_imgname;
public:
    PlotDistr1D(
	const std::string &title, const std::string &axis,
	const MathTemplates::SortedChain<MathTemplates::value<numtX>> &data,
	const std::string&imgname=""
    ): MathTemplates::Distribution1D<numtX,numtY>(data), m_title(title), m_axis(axis),m_imgname(imgname) {}
    virtual ~PlotDistr1D()
    {
        Plot<TUse>(m_imgname).template Hist<numtX,numtY>(*this) << "set title '" + m_title + "'" << "set yrange [0:]"
                                 << "set xlabel '" + m_axis + "'" << "set ylabel 'counts'";
    }
};
template<class numtX = double,class numtY = numtX,class numtZ = numtY,TDecl_init>
class PlotDistr2D: public MathTemplates::Distribution2D<numtX,numtY,numtZ>
{
private:
    std::string m_title,m_imgname;
public:
    PlotDistr2D(
	const std::string &title,
	const MathTemplates::SortedChain<MathTemplates::value<numtX>> &X,
	const MathTemplates::SortedChain<MathTemplates::value<numtY>> &Y,
	const std::string&imgname=""
    ): MathTemplates::Distribution2D<numtX,numtY,numtZ>(X, Y), m_title(title),m_imgname(imgname){}
    virtual ~PlotDistr2D()
    {
        PlotHist2d<TUse>(sp2,m_imgname).template Distr<numtX,numtY,numtZ>(*this) << "set title '" + m_title + "'";
    }
};
#undef TDecl_init
#undef TDecl
#undef TUse
};
#endif

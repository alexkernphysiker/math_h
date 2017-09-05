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
template<class numt = double>
class Plotter
{
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
            return "." + m_prefix + ".numeric-data." + cnt + ".txt";
        } else {
            return name + ".txt";
        }
    }
public:
    Plotter()
    {
        terminal_counter = 0;
        filename_counter = 0;
        outpath = ".";
    }
    Plotter &operator<<(const std::string &line)
    {
        lines.push_back(line);
        return *this;
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
    const std::string GetTerminal(const std::string &name = "")
    {
        const std::string firstline = "set terminal pngcairo size 1024,868 font 'Verdana,18'\n";
        if (name == "") {
            terminal_counter++;
            auto cnt = std::to_string(terminal_counter);
            while (cnt.length() < 5)cnt = "0" + cnt;
            return firstline + "set output 'z_" + m_prefix + "-plot-" + cnt + ".png'";
        } else {
            return firstline + "set output '" + name + ".png'";
        }
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
    const MathTemplates::SortedPoints<numt> GetPoints2(const std::string &name = "")
    {
        MathTemplates::SortedPoints<numt> res;
        auto str = GetInput(name);
        numt x, y;
        while (str >> x >> y) {
            res << MathTemplates::point<numt>(x, y);
        }
        return res;
    }
    const MathTemplates::SortedPoints<MathTemplates::value<numt>> GetPoints4(const std::string &name = "")
    {
        MathTemplates::SortedPoints<MathTemplates::value<numt>> res;
        auto str = GetInput(name);
        numt x, y, dx, dy;
        while (str >> x >> y >> dx >> dy) {
            res << MathTemplates::point<MathTemplates::value<numt>>({x, dx}, {y, dy});
        }
        return res;
    }
};
template<class numt = double>class Plot
{
private:
    std::vector<std::string> lines;
    std::vector<std::string> plots;
public:
    typedef std::function<numt(numt)> FUNC;
    typedef std::function<void(std::ofstream &)> PLOTOUTPUT;
    Plot &operator<<(const std::string &line)
    {
        lines.push_back(line);
        return *this;
    }
    Plot(const std::string &name = "")
    {
        operator<<(Plotter<numt>::Instance().GetTerminal(name));
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
            Plotter<numt>::Instance() << line;
        for (int i = 0, n = plots.size(); i < n; i++) {
            std::string line = plots[i];
            if (i == 0)
                line = "plot " + line;
            if (i < (n - 1))
                line += ",\\";
            Plotter<numt>::Instance() << line;
        }
    }
    Plot &Object(const std::string &plot)
    {
        plots.push_back(plot);
        return *this;
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
    Plot &OutputPlot(const PLOTOUTPUT delegate, const std::string &options, const std::string &title = "", const std::string &name = "")
    {
	std::ofstream str;
        auto n = Plotter<numt>::Instance().File(str,name);
        if (str) {
            delegate(str);
            File(n, options, title);
            str.close();
        }
        return *this;
    }
    Plot &Line(const std::vector<MathTemplates::point<numt>> &points, const std::string &title = "", const std::string &name = "")
    {
        Plot<numt>::OutputPlot([&points](std::ofstream & data) {
            for (const auto &p : points)
                data << p.X() << " " << p.Y() << std::endl;
        }, "w l", title, name);
        return *this;
    }
    Plot &Line(const MathTemplates::SortedPoints<numt> &points, const std::string &title = "", const std::string &name = "")
    {
        Plot<numt>::OutputPlot([&points](std::ofstream & data) {
            for (const auto &p : points)
                data << p.X() << " " << p.Y() << std::endl;
        }, "w l", title, name);
        return *this;
    }
    Plot &Points(const std::vector<MathTemplates::point<numt>> &points, const std::string &title = "", const std::string &name = "")
    {
        Plot<numt>::OutputPlot([&points](std::ofstream & data) {
            for (const auto &p : points)
                data << p.X() << " " << p.Y() << std::endl;
        }, "using 1:2", title, name);
        return *this;
    }
    Plot &Points(const MathTemplates::SortedPoints<numt> &points, const std::string &title = "", const std::string &name = "")
    {
        Plot<numt>::OutputPlot([&points](std::ofstream & data) {
            for (const auto &p : points)
                data << p.X() << " " << p.Y() << std::endl;
        }, "using 1:2", title, name);
        return *this;
    }
    Plot &Hist(const MathTemplates::SortedPoints<MathTemplates::value<numt>, MathTemplates::value<numt>> &data,
               const std::string &title = "", const std::string &name = "")
    {
        Plot<numt>::OutputPlot([&data](std::ofstream & str) {
            for (const auto p : data)
                str << p.X().val() << " " << p.Y().val() << " " << p.X().uncertainty() << " " << p.Y().uncertainty() << std::endl;
        }, "using 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars", title, name);
        return *this;
    }
};

enum TypeOf3D {normal, sp2};
template<class numt=double>
class PlotHist2d
{
private:
    std::vector<std::string> lines;
    std::vector<std::string> plots;
    std::string Surf2File(const MathTemplates::BiSortedPoints<numt> &D, const std::string &name = "")const
    {
        std::ofstream str;
        auto n = Plotter<numt>::Instance().File(str,name);
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
    std::string Points2File(const std::vector<MathTemplates::point3d<numt>> &points, const std::string &name = "")
    {
        std::ofstream str;
        auto n = Plotter<numt>::Instance().File(str,name);
        if (str) {
            for (const auto &p : points)
                str << p.X() << " " << p.Y() << " " << p.Z() << std::endl;
            str.close();
        }
        return n;
    }
    std::string Distr2File(const MathTemplates::hist2d<numt> &D, const std::string &name = "")const
    {
        std::ofstream str;
        auto n = Plotter<numt>::Instance().File(str,name);
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
    std::string Points2File(const std::vector <
                            MathTemplates::point3d<MathTemplates::value<numt>,
                            MathTemplates::value<numt>, MathTemplates::value<numt> >> &points,
                            const std::string &name = "")
    {
        std::ofstream str;
        auto n = Plotter<numt>::Instance().File(str,name);
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
        operator<<(Plotter<numt>::Instance().GetTerminal(name));
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
        for (const std::string &line : lines)Plotter<numt>::Instance() << line;
        for (int i = 0, n = plots.size(); i < n; i++) {
            std::string line = plots[i];
            if (i == 0)
                line = "splot " + line;
            if (i < (n - 1))
                line += ",\\";
            Plotter<numt>::Instance() << line;
        }
    }
    PlotHist2d &Object(const std::string&&plot)
    {
        plots.push_back(plot);
        return *this;
    }
    PlotHist2d &Surface(const MathTemplates::BiSortedPoints<numt> &D, const std::string &title = "")
    {
        return Object(std::string("'") + Surf2File(D) + "' matrix nonuniform title'" + title + "'");
    }

    PlotHist2d &Distr(const MathTemplates::hist2d<numt> &D, const std::string &title = "")
    {
        return Object(std::string("'") + Distr2File(D) + "' matrix nonuniform title'" + title + "'");
    }
    PlotHist2d &Points(const std::vector<MathTemplates::point3d<
                       MathTemplates::value<numt>, MathTemplates::value<numt>,
                       MathTemplates::value<numt>>> &points, const std::string &title = "")
    {
        return Object(std::string("'") + Points2File(points) + "' u 1:2:3 w points title'" + title + "'");
    }
    PlotHist2d &Points(const std::initializer_list <
                       MathTemplates::point3d<MathTemplates::value<numt>,
                       MathTemplates::value<numt>, MathTemplates::value<numt> >> &points, const std::string &title = "")
    {
        return Points(std::vector<MathTemplates::point3d<
                      MathTemplates::value<numt>, MathTemplates::value<numt>,
                      MathTemplates::value<numt>>>(points), title);
    }

    PlotHist2d &Line(const std::vector<MathTemplates::point3d<
                     MathTemplates::value<numt>, MathTemplates::value<numt>,
                     MathTemplates::value<numt>>> &points, const std::string &title = "")
    {
        return Object(std::string("'") + Points2File(points) + "' u 1:2:3 w line title'" + title + "'");
    }

    PlotHist2d &Points(const std::vector<MathTemplates::point3d<numt>> &points, const std::string &title = "")
    {
        return Object(std::string("'") + Points2File(points) + "' u 1:2:3 w points title'" + title + "'");
    }
    PlotHist2d &Points(const std::initializer_list<MathTemplates::point3d<numt>> &points, const std::string &title = "")
    {
        return Points(std::vector<MathTemplates::point3d<numt>>(points), title);
    }

    PlotHist2d &Line(const std::vector<MathTemplates::point3d<numt>> &points, const std::string &title = "")
    {
        return Object(std::string("'") + Points2File(points) + "' u 1:2:3 w line title'" + title + "'");
    }

};
template<class numt = double>
class PlotDistr1D: public MathTemplates::Distribution1D<numt>
{
private:
    std::string m_title, m_axis,m_imgname;
public:
    PlotDistr1D(
	const std::string &title, const std::string &axis,
	const MathTemplates::SortedChain<MathTemplates::value<numt>> &data,
	const std::string&imgname=""
    ): MathTemplates::Distribution1D<numt>(data), m_title(title), m_axis(axis),m_imgname(imgname) {}
    virtual ~PlotDistr1D()
    {
        Plot<numt>(m_imgname).Hist(*this) << "set title '" + m_title + "'" << "set yrange [0:]"
                                 << "set xlabel '" + m_axis + "'" << "set ylabel 'counts'";
    }
};
template<class numt = double>
class PlotDistr2D: public MathTemplates::Distribution2D<numt>
{
private:
    std::string m_title,m_imgname;
public:
    PlotDistr2D(
	const std::string &title,
	const MathTemplates::SortedChain<MathTemplates::value<numt>> &X,
	const MathTemplates::SortedChain<MathTemplates::value<numt>> &Y,
	const std::string&imgname=""
    ): MathTemplates::Distribution2D<numt>(X, Y), m_title(title),m_imgname(imgname){}
    virtual ~PlotDistr2D()
    {
        PlotHist2d<numt>(sp2,m_imgname).Distr(*this) << "set title '" + m_title + "'";
    }
};
};
#endif

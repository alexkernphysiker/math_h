// this file is distributed under
// LGPLv3 license
#ifndef VIJVUSSC
#define VIJVUSSC
#if __cplusplus<201100L
#error c++>=11 is needed for using math_h headers
#endif
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <unistd.h>
#include <math.h>
#include "math_h/error.h"
#include "math_h/tabledata.h"
#include "math_h/sigma.h"
#include "math_h/sigma2.h"
#include "math_h/sigma3.h"
namespace GnuplotWrap
{
    class PlotHist2d;
    class Plot;
    class Plotter
    {
        friend class Plot;
        friend class PlotHist2d;
    private:
        MathTemplates::Chain<std::string> lines;
        unsigned int terminal_counter, filename_counter;
        std::string outpath;
        std::string m_prefix;
        inline std::string GetFileName(const std::string& name)
        {
            if (name == "") {
                filename_counter++;
                auto cnt = std::to_string(filename_counter);
                while (cnt.length() < 6)cnt = "0" + cnt;
                return "." + m_prefix + ".numeric-data." + cnt + ".txt";
            }
            else {
                return name + ".txt";
            }
        }
    protected:
        inline Plotter& operator<<(const std::string& line)
        {
            lines.push_back(line);
            return *this;
        }
        inline std::string File(std::ofstream& str, const std::string& name = "")
        {
            const auto n = GetFileName(name);
            str.open((outpath + "/" + n).c_str());
            return n;
        }
        inline std::string GetTerminal(const std::string& name, int fontsize, int imgsize)
        {
            static const std::vector<std::string> IMGSIZE{ "640,480","800,600","1000,750","1600,1200","2000,1500" };

            const std::string firstline = "set terminal pngcairo size " + IMGSIZE[imgsize] + " font 'Verdana," + std::to_string(int(14) + 4 * imgsize + 2 * fontsize) + "'\n"
                + "set termoption enhanced\nset encoding utf8\n";
            if (name == "") {
                terminal_counter++;
                auto cnt = std::to_string(terminal_counter);
                while (cnt.length() < 5)cnt = "0" + cnt;
                return firstline + "set output '" +
                    m_prefix + "-plot-" + cnt + ".png'";
            }
            else {
                return firstline + "set output '" + name + ".png'";
            }
        }
    public:
        inline Plotter()
        {
            terminal_counter = 0;
            filename_counter = 0;
            outpath = ".";
        }
        inline ~Plotter()
        {
            std::string name = m_prefix + ".gnuplot-script";
            std::ofstream str(std::string(outpath + "/" + name).c_str(), std::ofstream::out);
            for (const auto& line : lines)
                str << line << "\n";
            str.close();
            name = "gnuplot " + name;
            const std::string old = getcwd(NULL, 0);
            chdir(outpath.c_str());
            system(name.c_str());
            chdir(old.c_str());
        }
        inline static Plotter& Instance()
        {
            static Plotter m_instance;
            return m_instance;
        }
        inline void SetOutput(const std::string& out, const std::string& prefix = "")
        {
            if (lines.size() > 0)
                throw MathTemplates::Exception<Plotter>("Attempt to reset plot output settings");
            outpath = out;
            m_prefix = prefix;
        }
        template<class numtX = double, class numtY = numtX>
        MathTemplates::Points<numtX, numtY> GetPoints(const std::string& name)
        {
            MathTemplates::Points<numtX, numtY> res;
            std::ifstream str((outpath + "/" + GetFileName(name)).c_str());
            numtX x;numtY y;
            while (str >> x >> y) {
                res.push_back(MathTemplates::make_point(x, y));
            }
            return res;
        }
        template<class numtX = double, class numtY = numtX>
        std::string SavePoints(const MathTemplates::Points<numtX, numtY>& data, const std::string& name = "")
        {
            std::ofstream str;
            auto n = Plotter::Instance().File(str, name);
            if (str) {
                for (const auto& P : data)
                    str << P.X() << "\t" << P.Y() << std::endl;
                str.close();
            }
            return n;
        }
    };
    class Plot
    {
    private:
        MathTemplates::Chain<std::string> lines;
        MathTemplates::Chain<std::string> plots;
        MathTemplates::Chain<std::string> finalization;
        inline Plot& Object(const std::string& plot)
        {
            plots.push_back(plot);
            return *this;
        }
    public:
        inline Plot& operator<<(const std::string& line)
        {
            lines.push_back(line);
            return *this;
        }
        inline Plot& operator>>(const std::string& line)
        {
            finalization.push_back(line);
            return *this;
        }
        inline Plot(const std::string& name = "", int fontsize = 0, int imgsize = 2)
        {
            operator<<(Plotter::Instance().GetTerminal(name, fontsize, imgsize));
            operator<<("unset pm3d");
            operator<<("unset title");
            operator<<("unset key");
            operator<<("unset surface");
            operator<<("unset view");
            operator<<("unset xrange");
            operator<<("unset yrange");
            operator<<("unset xlabel");
            operator<<("unset ylabel");
            operator<<("set xtics");
            operator<<("set ytics");
        }
        inline ~Plot()
        {
            for (const std::string& line : lines)
                Plotter::Instance() << line;
            for (int i = 0, n = plots.size(); i < n; i++) {
                std::string line = plots[i];
                if (i == 0)
                    line = "plot " + line;
                if (i < (n - 1))
                    line += ",\\";
                Plotter::Instance() << line;
            }
            for (const std::string& line : finalization)
                Plotter::Instance() << line;
        }
        inline Plot& File(const std::string& name, const std::string& options, const std::string& title)
        {
            std::string line = "\"" + name + "\" ";
            line += options;
            line += " title \"";
            line += title;
            line += "\"";
            return Object(line);
        }
        template<class numtX = double, class numtY = numtX>
        inline Plot& Output(const MathTemplates::Points<numtX, numtY>& data, const std::string& options, const std::string& title = "", const std::string& name = "")
        {
            File(Plotter::Instance().SavePoints(data, name), options, title);
            return *this;
        }
        template<class numtX = double, class numtY = numtX>
        inline Plot& Line(const MathTemplates::Points<numtX, numtY>& points, const std::string& title = "", const std::string& name = "", const std::string& mod = "")
        {
            return Output(points, "w l " + mod, title, name);
        }
        template<class numtX = double, class numtY = numtX>
        inline Plot& Line(const MathTemplates::SortedPoints<numtX, numtY>& points, const std::string& title = "", const std::string& name = "", const std::string& mod = "")
        {
            return Output(points(), "w l " + mod, title, name);
        }
        template<class numtX = double, class numtY = numtX>
        inline Plot& Points(const MathTemplates::Points<numtX, numtY>& points, const std::string& title = "", const std::string& name = "", const std::string& mod = "")
        {
            return Output<numtX, numtY>(points, "using 1:2 " + mod, title, name);
        }
        template<class numtX = double, class numtY = numtX>
        inline Plot& Points(const MathTemplates::SortedPoints<numtX, numtY>& points, const std::string& title = "", const std::string& name = "", const std::string& mod = "")
        {
            return Output<numtX, numtY>(points(), "using 1:2 " + mod, title, name);
        }
        template<class numtX = double, class numtY = numtX>
        inline Plot& Hist(const MathTemplates::SortedPoints<MathTemplates::value<numtX>, MathTemplates::value<numtY>>& data,
            const std::string& title = "", const std::string& name = "", const std::string& mod = "")
        {
            return Output(data(), "using 1:3:($1-$2):($1+$2):($3-$4):($3+$4) with xyerrorbars " + mod, title, name);
        }
        template<size_t index, size_t sz, class numtX = double, class numtY = numtX>
        inline Plot& Hist(const MathTemplates::SortedPoints<MathTemplates::value<numtX>, MathTemplates::Uncertainties<sz, numtY>>& data,
            const std::string& title = "", const std::string& name = "", const std::string& mod = "")
        {
            return Output(data(), "using 1:3:($1-$2):($1+$2):($3-$" + std::to_string(3 + index) + "):($3+$" + std::to_string(3 + index) + ") with xyerrorbars " + mod, title, name);
        }
        template<size_t index, size_t index2, size_t sz, class numtX = double, class numtY = numtX>
        inline Plot& Hist_2bars(const MathTemplates::SortedPoints<MathTemplates::value<numtX>, MathTemplates::Uncertainties<sz, numtY>>& data,
            const std::string& title1 = "", const std::string& title2 = "", const std::string& name = "", const std::string& mod = "")
        {
            return Output(data(), "using 1:3:($1-$2):($1+$2):($3-$" + std::to_string(3 + index) + "):($3+$" + std::to_string(3 + index) + ") with xyerrorbars title'" + title1 + "'," +
                "'' using 1:3:($3-$" + std::to_string(3 + index2) + "):($3+$" + std::to_string(3 + index2) + ") with yerrorbars " + mod, title2, name);
        }
        template<class numtX = double, class numtY = numtX>
        inline Plot& XYUncertainties(const MathTemplates::Points<MathTemplates::value<numtX>, MathTemplates::value<numtY>>& data,
            const std::string& title = "", const std::string& name = "", const std::string& mod = "")
        {
            return Output(data, "using 1:3:($1-$2):($1+$2):($3-$4):($3+$4) with xyerrorbars " + mod, title, name);
        }
        template<class numtX = double, class numtY = numtX>
        inline Plot& YUncertainties(const MathTemplates::Points<numtX, MathTemplates::value<numtY>>& data,
            const std::string& title = "", const std::string& name = "", const std::string& mod = "")
        {
            return Output(data, "using 1:2:($2-$3):($2+$3) with yerrorbars " + mod, title, name);
        }
        template<class numtX = double, class numtY = numtX>
        inline Plot& XUncertainties(const MathTemplates::Points<MathTemplates::value<numtX>, numtY>& data,
            const std::string& title = "", const std::string& name = "", const std::string& mod = "")
        {
            return Output(data, "using 1:3:($1-$2):($1+$2) with xerrorbars " + mod, title, name);
        }

    };

    enum TypeOf3D { normal, sp2 };
    class PlotHist2d
    {
    private:
        MathTemplates::Chain<std::string> lines;
        MathTemplates::Chain<std::string> plots;
        MathTemplates::Chain<std::string> finalization;
        template<class numtX = double, class numtY = numtX, class numtZ = numtY>
        std::string Surf2File(const MathTemplates::BiSortedPoints<numtX, numtY, numtZ>& D, const std::string& name = "")const
        {
            std::ofstream str;
            auto n = Plotter::Instance().File(str, name);
            if (str) {
                str << D.X().size() << " ";
                for (const auto& x : D.X())
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
        template<class numtX = double, class numtY = numtX, class numtZ = numtY>
        std::string Points2File(const MathTemplates::Chain<MathTemplates::point3d<numtX, numtY, numtZ>>& points, const std::string& name = "")
        {
            std::ofstream str;
            auto n = Plotter::Instance().File(str, name);
            if (str) {
                for (const auto& p : points)
                    str << p.X() << " " << p.Y() << " " << p.Z() << std::endl;
                str.close();
            }
            return n;
        }
        template<class numtX = double, class numtY = numtX, class numtZ = numtY>
        std::string Distr2File(const MathTemplates::hist2d<numtX, numtY, numtZ>& D, const std::string& name = "")const
        {
            std::ofstream str;
            auto n = Plotter::Instance().File(str, name);
            if (str) {
                str << D.X().size() << " ";
                for (const auto& x : D.X())
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
        template<class numtX = double, class numtY = numtX, class numtZ = numtY>
        std::string Points2File(const MathTemplates::Chain <
            MathTemplates::point3d<MathTemplates::value<numtX>,
            MathTemplates::value<numtY>, MathTemplates::value<numtZ> >> &points,
            const std::string& name = "")
        {
            std::ofstream str;
            auto n = Plotter::Instance().File(str, name);
            if (str) {
                for (const auto& p : points)
                    str << p.X().val() << " " << p.Y().val() << " " << p.Z().val() << std::endl;
                str.close();
            }
            return n;
        }
    public:
        inline PlotHist2d& operator<<(const std::string& line)
        {
            lines.push_back(line);
            return *this;
        }
        inline PlotHist2d& operator>>(const std::string& line)
        {
            finalization.push_back(line);
            return *this;
        }
        inline PlotHist2d(const TypeOf3D type, const std::string& name = "", int fontsize = 0, int imgsize = 2)
        {
            operator<<(Plotter::Instance().GetTerminal(name, fontsize, imgsize));
            operator<<("unset title");
            operator<<("unset key");
            operator<<("unset surface");
            if (sp2 == type) {
                operator<<("set view map");
                operator<<("set pm3d at b");
                operator<<("unset colorbox");
            }
            else {
                operator<<("unset view");
                operator<<("set pm3d");
            }
            operator<<("unset xrange");
            operator<<("unset yrange");
            operator<<("unset zrange");
            operator<<("unset xlabel");
            operator<<("unset ylabel");
            operator<<("unset zlabel");
            operator<<("set xtics");
            operator<<("set ytics");
            operator<<("set ztics");
        }
        inline ~PlotHist2d()
        {
            for (const std::string& line : lines)
                Plotter::Instance() << line;
            for (int i = 0, n = plots.size(); i < n; i++) {
                std::string line = plots[i];
                if (i == 0)
                    line = "splot " + line;
                if (i < (n - 1))
                    line += ",\\";
                Plotter::Instance() << line;
            }
            for (const std::string& line : finalization)
                Plotter::Instance() << line;
        }
        inline PlotHist2d& Object(const std::string& plot)
        {
            plots.push_back(plot);
            return *this;
        }
        template<class numtX = double, class numtY = numtX, class numtZ = numtY>
        inline PlotHist2d& Surface(const MathTemplates::BiSortedPoints<numtX, numtY, numtZ>& D, const std::string& title = "")
        {
            return Object(std::string("'") + Surf2File(D) + "' matrix nonuniform title'" + title + "'");
        }
        template<class numtX = double, class numtY = numtX, class numtZ = numtY>
        inline PlotHist2d& Distr(const MathTemplates::hist2d<numtX, numtY, numtZ>& D, const std::string& title = "")
        {
            return Object(std::string("'") + Distr2File(D) + "' matrix nonuniform title'" + title + "'");
        }
        template<class numtX = double, class numtY = numtX, class numtZ = numtY>
        inline PlotHist2d& Points(const MathTemplates::Chain<MathTemplates::point3d<
            MathTemplates::value<numtX>, MathTemplates::value<numtY>,
            MathTemplates::value<numtZ>>>& points, const std::string& title = "")
        {
            return Object(std::string("'") + Points2File(points) + "' u 1:2:3 w points title'" + title + "'");
        }
        template<class numtX = double, class numtY = numtX, class numtZ = numtY>
        inline PlotHist2d& Line(const MathTemplates::Chain<MathTemplates::point3d<
            MathTemplates::value<numtX>, MathTemplates::value<numtY>,
            MathTemplates::value<numtZ>>>& points, const std::string& title = "")
        {
            return Object(std::string("'") + Points2File(points) + "' u 1:2:3 w line title'" + title + "'");
        }
        template<class numtX = double, class numtY = numtX, class numtZ = numtY>
        inline PlotHist2d& Points(const MathTemplates::Chain<MathTemplates::point3d<numtX, numtY, numtZ>>& points, const std::string& title = "")
        {
            return Object(std::string("'") + Points2File(points) + "' u 1:2:3 w points title'" + title + "'");
        }
        template<class numtX = double, class numtY = numtX, class numtZ = numtY>
        inline PlotHist2d& Line(const MathTemplates::Chain<MathTemplates::point3d<numtX, numtY, numtZ>>& points, const std::string& title = "")
        {
            return Object(std::string("'") + Points2File(points) + "' u 1:2:3 w line title'" + title + "'");
        }

    };
    template<class numtX = double, class numtY = numtX>
    class PlotDistr1D : public MathTemplates::Distribution1D<numtX, numtY>
    {
    private:
        std::string m_title, m_axis, m_imgname;
        int m_fontsize, m_imgsize;
    public:
        PlotDistr1D(
            const std::string& title, const std::string& axis,
            MathTemplates::SortedChain<MathTemplates::value<numtX>>&& data,
            const std::string& imgname = "", int fontsize = 0, int imgsize = 2
        ) : MathTemplates::Distribution1D<numtX, numtY>(std::move(data)), m_title(title), m_axis(axis), m_imgname(imgname), m_fontsize(fontsize), m_imgsize(imgsize) {}
        virtual ~PlotDistr1D()
        {
            Plot(m_imgname, m_fontsize, m_imgsize).template Hist<numtX, numtY>(*this)
                << "set title '" + m_title + "'" << "set yrange [0:]"
                << "set xlabel '" + m_axis + "'" << "set ylabel 'counts'";
        }
    };
    template<class numtX = double, class numtY = numtX, class numtZ = numtY>
    class PlotDistr2D : public MathTemplates::Distribution2D<numtX, numtY, numtZ>
    {
    private:
        std::string m_title, m_imgname;
        int m_fontsize, m_imgsize;
    public:
        PlotDistr2D(
            const std::string& title,
            MathTemplates::SortedChain<MathTemplates::value<numtX>>&& X,
            MathTemplates::SortedChain<MathTemplates::value<numtY>>&& Y,
            const std::string& imgname = "", int fontsize = 0, int imgsize = 2
        ) : MathTemplates::Distribution2D<numtX, numtY, numtZ>(std::move(X), std::move(Y)), m_title(title), m_imgname(imgname), m_fontsize(fontsize), m_imgsize(imgsize) {}
        virtual ~PlotDistr2D()
        {
            PlotHist2d(sp2, m_imgname, m_fontsize, m_imgsize).template Distr<numtX, numtY, numtZ>(*this) << "set title '" + m_title + "'";
        }
    };
};
#endif

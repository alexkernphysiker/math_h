INCLUDEPATH +=../..
QMAKE_CXXFLAGS+= -std=c++11
QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
TARGET = mandelbrot.example
TEMPLATE = app
SOURCES += mandelbrotwindow.cpp mandelbrot_main.cpp
HEADERS  += mandelbrotwindow.h
FORMS    += mandelbrotwindow.ui

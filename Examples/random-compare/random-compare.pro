INCLUDEPATH +=../..
QMAKE_CXXFLAGS+= -std=c++11
QT       -= gui
TARGET = random-compare.example
TEMPLATE = app
HEADERS += randcmp.h
SOURCES += device.cpp pseudo.cpp randcmp.cpp

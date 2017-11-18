TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
INCLUDEPATH +=  /usr/local/Cellar/armadillo/8.200.0/include/
LIBS += -L /usr/local/Cellar/armadillo/8.200.0/lib -larmadillo

SOURCES += main_analytical_results.cpp

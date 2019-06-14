#-------------------------------------------------
#
# Project created by QtCreator 2015-03-18T16:29:16
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = EngCalcProgramm
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    qcustomplot.cpp

HEADERS  += mainwindow.h \
    qcustomplot.h

FORMS    += mainwindow.ui

QMAKE_LFLAGS += -static-libgcc

RESOURCES += \
    degs/degs.qrc

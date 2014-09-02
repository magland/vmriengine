TEMPLATE = app

OBJECTS_DIR = build
DESTDIR = bin
TARGET = vmriengine

INCLUDEPATH += .
DEPENDPATH += .
SOURCES += vmrimain.cpp

INCLUDEPATH += core
DEPENDPATH += core
HEADERS += vmriengine.h abstractsimsequence.h abstractphantom.h vmriinterface.h
SOURCES += vmriengine.cpp abstractsimsequence.cpp abstractphantom.cpp vmriinterface.cpp
HEADERS += textfile.h qjson.h parser.h qjson_export.h serializer.h qwait.h
SOURCES += textfile.cpp qjson.cpp qwait.cpp
LIBS += -lqjson

INCLUDEPATH += phantoms
DEPENDPATH += phantoms
HEADERS += staticphantom.h phantom1.h
SOURCES += staticphantom.cpp phantom1.cpp

INCLUDEPATH += sequences
DEPENDPATH += sequences
HEADERS += sequence1.h
SOURCES += sequence1.cpp

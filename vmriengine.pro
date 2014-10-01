TEMPLATE = app

OBJECTS_DIR = build
DESTDIR = bin
TARGET = vmriengine

INCLUDEPATH += .
DEPENDPATH += .
SOURCES += vmrimain.cpp

INCLUDEPATH += ./core
DEPENDPATH += ./core
HEADERS += core/vmriengine.h core/abstractsimsequence.h core/abstractphantom.h core/vmriinterface.h
SOURCES += core/vmriengine.cpp core/abstractsimsequence.cpp core/abstractphantom.cpp core/vmriinterface.cpp
HEADERS += core/textfile.h core/qjson.h core/parser.h core/qjson_export.h core/serializer.h core/qwait.h
SOURCES += core/textfile.cpp core/qjson.cpp core/qwait.cpp
LIBS += -lqjson

INCLUDEPATH += ./phantoms
DEPENDPATH += ./phantoms
HEADERS += phantoms/staticphantom.h phantoms/phantom1.h
SOURCES += phantoms/staticphantom.cpp phantoms/phantom1.cpp

INCLUDEPATH += ./sequences
DEPENDPATH += ./sequences
HEADERS += sequences/sequence1.h
SOURCES += sequences/sequence1.cpp

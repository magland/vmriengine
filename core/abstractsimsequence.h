#ifndef abstractsimsequence_h
#define abstractsimsequence_h

#include <QString>
#include <QList>
#include "vmriinterface.h" //for Matrix44
#include <QMap>
#include <QVariant>

#define COMMAND_TYPE_addRFOperator 1
#define COMMAND_TYPE_initialize 2
#define COMMAND_TYPE_evolve 3
#define COMMAND_TYPE_excite 4
#define COMMAND_TYPE_readout 5


struct SimCommand {
	QString command;
	QMap<QString,QVariant> data;
};

class AbstractSimSequence {
public:
	AbstractSimSequence() {}
	virtual ~AbstractSimSequence();
	
	virtual bool atEnd()=0;
	virtual SimCommand nextSimCommand()=0;
};

#endif

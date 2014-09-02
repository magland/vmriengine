#include "sequence1.h"
#include <QList>
#include <stdio.h>
#include "../core/qjson.h"
#include "../core/textfile.h"
#include <QDebug>

class Sequence1Private {
public:
	
	QList<SimCommand> m_commands;
	int m_command_index;
};

Sequence1::Sequence1() {
	d=new Sequence1Private;
	d->m_command_index=0;
}

void Sequence1::loadFromFile(const QString &fname) {
	d->m_commands.clear();
	
	QMap<QString,QVariant> simblocks=parseJSON(read_text_file(fname)).toMap();
	
	QList<QVariant> rf_waveforms=simblocks["rf_waveforms"].toList();
	QList<QVariant> blocks=simblocks["blocks"].toList();
	
	{
		SimCommand CC;
		CC.command="initialize";
		d->m_commands << CC;
	}
	for (int i=0; i<rf_waveforms.count(); i++) {
		SimCommand CC;
		CC.command="add_RF_waveform";
		CC.data=rf_waveforms[i].toMap();
		d->m_commands << CC;
	}
	for (int i=0; i<blocks.count(); i++) {
		SimCommand CC;
		CC.command="run_block";
		CC.data=blocks[i].toMap();
		d->m_commands << CC;
	}
	{
		SimCommand CC;
		CC.command="finalize";
		d->m_commands << CC;
	}
}

Sequence1::~Sequence1() {
	delete d;
}
	
bool Sequence1::atEnd() {
	return (d->m_command_index>=d->m_commands.count());
}
SimCommand Sequence1::nextSimCommand() {
	return d->m_commands[d->m_command_index++];
}

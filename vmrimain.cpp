#include "core/vmriengine.h"
#include "phantoms/phantom1.h"
#include "sequences/sequence1.h"
#include <QMap>
#include <QString>
#include <QStringList>
#include <stdio.h>
#include "qwait.h"
#include <QTime>
#include <QDebug>
#include <QDir>
#include <QDateTime>

typedef QList<VmriReadout> VmriReadoutList;

void usage() {
	printf("vmriengine -simblocks [fname.json] -out [fname.dat] -iso [num] -threads [num]\n");
	printf("vmriengine -combine [fname_prefix] -out [fname.dat]\n");
}

QStringList get_combine_paths(const QString &prefix) {
	
	
	QString str1=".";
	QString str2=prefix;
	int ind=prefix.lastIndexOf("/");
	if (ind>=0) {
		str1=prefix.mid(0,ind);
		str2=prefix.mid(ind+1);
	}
	QDir dir(str1);
	QStringList list=dir.entryList(QStringList(str2+"*.dat"),QDir::Files,QDir::Name);
	
	QStringList ret;
	for (int i=0; i<list.count(); i++) {
		if (str1==".") ret << list[i];
		else ret << str1+"/"+list[i];
	}
	
	return ret;
}

QList<VmriReadout> combine_readout_lists(const QList<VmriReadoutList> &X) {
	QList<VmriReadout> readouts;
	if (X.count()==0) return readouts;
	readouts=X[0];
	for (int tt=1; tt<X.count(); tt++) {
		QList<VmriReadout> readouts2=X[tt];
		for (int j=0; j<readouts.count(); j++) {
			VmriReadout *RO=&readouts[j];
			VmriReadout *RO2=&readouts2[j];
			for (int k=0; k<RO->N; k++) {
				RO->data_real[k]+=RO2->data_real.value(k);
				RO->data_imag[k]+=RO2->data_imag.value(k);
			}
		}
	}
	return readouts;
}

void write_readouts_to_file(const QString &fname,const QList<VmriReadout> &readouts) {
	FILE *outf=fopen(fname.toAscii().data(),"wb");
	for (int i=0; i<readouts.count(); i++) {
		VmriReadout RO=readouts[i];
		fwrite(&RO.N,sizeof(quint32),1,outf);
		fwrite(&RO.readout_index,sizeof(quint32),1,outf);
		fwrite(&RO.dwell_time,sizeof(float),1,outf);
		for (unsigned int i=0; i<29*4; i++) {
			unsigned char c=0;
			fwrite(&c,1,1,outf);
		}
		for (unsigned int i=0; i<RO.N; i++) {
			float re,im;
			re=RO.data_real.value(i);
			im=RO.data_imag.value(i);
		
			fwrite(&re,sizeof(float),1,outf);
			fwrite(&im,sizeof(float),1,outf);
		}
	}
	fclose(outf);
}

QList<VmriReadout> load_readouts_from_file(const QString &fname) {
	QList<VmriReadout> readouts;
	FILE *inf=fopen(fname.toAscii().data(),"rb");
	
	bool done=false;
	while (!done) {
		quint32 N;
		quint32 readout_index;
		float dwell_time;
		if ((!feof(inf))
				&&(fread(&N,sizeof(quint32),1,inf))
				&&(fread(&readout_index,sizeof(quint32),1,inf))
				&&(fread(&dwell_time,sizeof(float),1,inf))
			 ) {
			for (unsigned int i=0; i<29*4; i++) {
				unsigned char c;
				fread(&c,1,1,inf);
			}
			if ((N>0)&&(N<100000)) {
				VmriReadout RO;
				RO.N=N; RO.readout_index=readout_index; RO.dwell_time=dwell_time;
				for (unsigned int i=0; i<N; i++) {
					float re,im;
					fread(&re,sizeof(float),1,inf);
					fread(&im,sizeof(float),1,inf);
					RO.data_real << re;
					RO.data_imag << im;
				}
				readouts << RO;
			}
			else {
				qWarning() << "Unexpected value of N:" << N;
			}
		}
		else done=true;
	}
	
	fclose(inf);
	return readouts;
}



int main(int argc,char *argv[]) {
	
	qsrand(QDateTime::currentDateTime().toTime_t());
	
	QStringList args;
	for (int i=1; i<argc; i++) args << QString(argv[i]);
	
	QMap<QString,QString> params;
	for (int i=0; i<args.count(); i++) {
		QString key=args[i];
		if ((key.indexOf("-")==0)&&(i+1<args.count())) {
			params[key.mid(1)]=args[i+1];
		}
	}
	
	QString output_file=params.value("out","out.dat");
	
	QString combine_prefix=params.value("combine","");
	if (!combine_prefix.isEmpty()) {
		QStringList paths=get_combine_paths(combine_prefix);
		qDebug() << paths;
		QList<VmriReadoutList> readout_lists;
		for (int i=0; i<paths.count(); i++) {
			qDebug()  << QString("Loading %1").arg(paths[i]);
			readout_lists << load_readouts_from_file(paths[i]);
		}
		qDebug()  << "Combining...";
		QList<VmriReadout> readouts=combine_readout_lists(readout_lists);
		qDebug()  << "Writing...";
		write_readouts_to_file(output_file,readouts);
		qDebug()  << "Done.";
		return 0;
	}

	long num_iso=params.value("iso","10000").toInt();
	long num_threads=params.value("threads","1").toInt();
	
	QList<VmriEngine *> engines;
	
	int total_iso=0;
	for (int tt=0; tt<num_threads; tt++) {
		Sequence1 *SEQ=new Sequence1;
		QString simblocks_fname=params.value("simblocks","");
		if (!simblocks_fname.isEmpty()) {
			qDebug()  << QString("Loading sequence %1...").arg(simblocks_fname);
			SEQ->loadFromFile(simblocks_fname);
		}
		else {
			usage();
			return -1;
		}
		
		Phantom1 *PP=new Phantom1;
		
		int tmp=num_iso/num_threads;
		if (tt+1==num_threads) {
			tmp=num_iso-total_iso;
		}
		PP->setNumIsochromats(tmp);
		PP->setDensityFactor(1.0/num_threads);
		PP->initialize();
		total_iso+=tmp;
		
		VmriEngine *EE=new VmriEngine;
		EE->setSimSequence(SEQ);
		EE->setPhantom(PP);
		engines << EE;
	}
	
	QTime timer; timer.start();
	for (int tt=0; tt<num_threads; tt++) {
		engines[tt]->start();
	}

	bool done=false;
	while (!done) {
		done=true;
		for (int tt=0; tt<num_threads; tt++) {
			if (!engines[tt]->isFinished()) {
				done=false;
			}
		}
		qWait(100);
	}
	
	qDebug()  << QString("ELAPSED TIME: %1 ms").arg(timer.elapsed());
	
	//add the readouts together
	QList<VmriReadoutList> readout_lists;
	for (int tt=0; tt<num_threads; tt++) {
		readout_lists << engines[tt]->getReadouts();
	}
	QList<VmriReadout> readouts=combine_readout_lists(readout_lists);
	
	write_readouts_to_file(output_file,readouts);
	
	return 0;
}

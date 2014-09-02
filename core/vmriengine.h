#ifndef vmriengine_h
#define vmriengine_cpp
#include <QThread>

#include <QString>
#include "abstractsimsequence.h"
#include "abstractphantom.h"

struct VmriReadout {
	quint32 N;
	quint32 readout_index;
	float dwell_time;
	QList<double> data_real;
	QList<double> data_imag;
};

class VmriEnginePrivate;
class VmriEngine : public QThread {
public:
	friend class VmriEnginePrivate;
	
	VmriEngine();
	virtual ~VmriEngine();
	
	void setSimSequence(AbstractSimSequence *seq);
	void setPhantom(AbstractPhantom *phantom);
	QList<VmriReadout> getReadouts();
	
	void run();
	
private:
	VmriEnginePrivate *d;
};

#endif
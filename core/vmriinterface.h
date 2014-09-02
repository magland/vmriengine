#ifndef vmriinterface_H
#define vmriinterface_H

#include <QList>

struct Matrix44 {
	double data[4][4];
};

class VmriInterfacePrivate;
class VmriInterface {
public:
	friend class VmriInterfacePrivate;
	
	VmriInterface();
	virtual ~VmriInterface();
		
	void setT1(double T1);
	void setT2(double T2);
	void setIsochromats(int N,double *x,double *y,double *z,double *f,double *d);
	void addRFOperator(double freq_min,double freq_step,int num_freqs,const QList<Matrix44> &operators); //operators is of length num_freqs x 16
	void initialize();
	
	void evolve(const QList<double> &M,double t);
	void excite(int RF_index,const QList<double> &A,double phase,double frequency);
	void readout(double *ret_real,double *ret_imag,int num_readout_points,const QList<double> &A,double dt,double phase,double frequency);
	
private:
	VmriInterfacePrivate *d;
};

#endif

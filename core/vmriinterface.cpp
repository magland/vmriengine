#include "vmriinterface.h"
#include <QList>
#include <stdio.h>
#include <QDebug>

#define PI 3.141592653589793
#include <math.h>

struct Isochromat2 {
	double x,y,z,f;
	double Mx,My,Mz;
	double density;
};
struct RF_operator {
	double freq_min,freq_step;
	QList<Matrix44> matrices;
};

class VmriInterfacePrivate {
public:
	QList<Isochromat2> m_isochromats;
	QList<RF_operator> m_RF_operators;
	double m_T1,m_T2;
	double m_gamma;
	
	void rotateZ(double *M,double theta);
	void decayMag(double *M,double E1,double E2);
	double dot_product(const QList<double> &V1,double *V2);
	void apply_matrix(double *Mag,const Matrix44 &MM);
};

VmriInterface::VmriInterface() {
	d=new VmriInterfacePrivate;
	d->m_T1=300;
	d->m_T2=40;
	d->m_gamma=42.57;
}

VmriInterface::~VmriInterface() {
	delete d;
}

void VmriInterface::setT1(double T1) {
	d->m_T1=T1;
}
void VmriInterface::setT2(double T2) {
	d->m_T2=T2;
}

void VmriInterface::setIsochromats(int N,double *x,double *y,double *z,double *f,double *density) {
	d->m_isochromats.clear();
	for (int i=0; i<N; i++) {
		Isochromat2 II;
		II.x=x[i];
		II.y=y[i];
		II.z=z[i];
		II.f=f[i];
		II.density=density[i];
		II.Mx=0;
		II.My=0;
		II.Mz=1;
		d->m_isochromats << II;
	}
}
void VmriInterface::addRFOperator(double freq_min,double freq_step,int num_freqs,const QList<Matrix44> &operators) {
	RF_operator XX;
	XX.freq_min=freq_min;
	XX.freq_step=freq_step;
	for (int i=0; i<num_freqs; i++) {
		Matrix44 MM;
		int ct=0;
		for (int r=0; r<4; r++)
		for (int c=0; c<4; c++) {
			MM.data[r][c]=operators[i].data[r][c];
			ct++;
		}
		XX.matrices << MM;
	}
	d->m_RF_operators << XX;
}
void VmriInterface::initialize() {
}

double randu() {
	return ((qrand()%10000)+0.5)*1.0/10000;
}

void VmriInterface::evolve(const QList<double> &M,double t) {
	for (int i=0; i<d->m_isochromats.count(); i++) {
		Isochromat2 *II=&d->m_isochromats[i];
		double pos0[3];
		pos0[0]=II->x;
		pos0[1]=II->y;
		pos0[2]=II->z;
		double theta=d->dot_product(M,pos0)*d->m_gamma*2*PI/(1000*1000); //radians (uT/mm-us*mm*Hz/uT)
		theta+=(II->f)*t*2*PI/(1000*1000);
		double Mag[3];
		Mag[0]=II->Mx;
		Mag[1]=II->My;
		Mag[2]=II->Mz;
		d->rotateZ(Mag,theta);
		double E1=exp(-t/(d->m_T1*1000));
		double E2=exp(-t/(d->m_T2*1000));
		d->decayMag(Mag,E1,E2);
		II->Mx=Mag[0];
		II->My=Mag[1];
		II->Mz=Mag[2];
		
		/*
		//test some crazy diffusion
		for (int a=0; a<t/1000; a++) {
			II->x+=(randu()*2-1)*0.01;
			II->y+=(randu()*2-1)*0.01;
			II->z+=(randu()*2-1)*0.01;
		}
		*/
		
	}
}
void VmriInterface::excite(int RF_index,const QList<double> &A,double phase,double frequency) {
	for (int i=0; i<d->m_isochromats.count(); i++) {
		Isochromat2 *II=&d->m_isochromats[i];
		double pos0[3];
		pos0[0]=II->x;
		pos0[1]=II->y;
		pos0[2]=II->z;
		double freq=II->f + d->dot_product(A,pos0)*d->m_gamma+frequency; //Hz
		
		RF_operator *RF=&d->m_RF_operators[RF_index];
		int ind=round((freq-RF->freq_min)/RF->freq_step);
		if (ind<0) return;
		if (ind>=RF->matrices.count()) return;
		Matrix44 MM=RF->matrices[ind];
		
		double Mag[3];
		Mag[0]=II->Mx;
		Mag[1]=II->My;
		Mag[2]=II->Mz;
		double theta=phase*PI/180;
		d->rotateZ(Mag,theta);
		d->apply_matrix(Mag,MM);
		d->rotateZ(Mag,-theta);
		II->Mx=Mag[0];
		II->My=Mag[1];
		II->Mz=Mag[2];
	}
}
void VmriInterface::readout(double *output_real,double *output_imag,int N,const QList<double> &A,double dt,double phase,double frequency) {
	
	//we are going to accumulate the spins in the Fourier domain
	//this has a couple advantages... first it may be faster, second we can filter out the high freqs
	
	QList<double> frequencies; //Hz
	QList<double> signals_real;
	QList<double> signals_imag;
	
	for (int ii=0; ii<d->m_isochromats.count(); ii++) {
		Isochromat2 *II=&d->m_isochromats[ii];
		double pos0[3];
		pos0[0]=II->x;
		pos0[1]=II->y;
		pos0[2]=II->z;
		frequencies << II->f+d->dot_product(A,pos0)*d->m_gamma+frequency; //this will be in Hz
		signals_real << II->Mx*II->density;
		signals_imag << II->My*II->density;
	}
	
	int oversamp=4; //oversampling factor
	QList<double> X_real;
	QList<double> X_imag;
	for (int i=0; i<N*oversamp; i++) {
		X_real << 0;
		X_imag << 0;
	}
	//The sampling time is dt
	//so the sampling bandwidth is 1/dt
	double sampling_bandwidth=1.0/dt*(1000*1000); //Hz
	double freq_step=sampling_bandwidth/(N*oversamp);
	double min_freq=-freq_step*(N*oversamp*1.0/2);
	double max_freq=freq_step*(N*oversamp*1.0/2-1);
	double min0=0;
	double max0=0;
	for (int i=0; i<frequencies.count(); i++) {
		double f0=frequencies[i];
		min0=qMin(f0,min0); max0=qMax(f0,max0);
		if ((min_freq<=f0)&&(f0<=max_freq)) {
			int ind0=round((f0-min_freq)/freq_step);
			if ((0<=ind0)&&(ind0<N*oversamp)) {
				//in the future, consider spreading this using sinc interpolation
				X_real[ind0]+=signals_real[i];
				X_imag[ind0]+=signals_imag[i];
			}
		}
	}
	
	for (int t=0; t<N; t++) {
		output_real[t]=0;
		output_imag[t]=0;
	}
	
	for (int i=0; i<N*oversamp; i++) {
		double re0=X_real[i];
		double im0=X_imag[i];
		if ((re0)||(im0)) {
			double freq=min_freq+freq_step*i;
			double dtheta=freq*dt/(1000*1000)*2*PI;
			double cosdtheta=cos(dtheta);
			double sindtheta=sin(dtheta);
			double cosdtheta2=cos(dtheta/2);
			double sindtheta2=sin(dtheta/2);
		
			//evolve by half of a step
			double re1=cosdtheta2*re0-sindtheta2*im0;
			double im1=sindtheta2*re0+cosdtheta2*im0;
			re0=re1;
			im0=im1;

			for (int t=0; t<N; t++) {
				output_real[t]+=re0;
				output_imag[t]+=im0;
				
				//evolve by full step
				re1=cosdtheta*re0-sindtheta*im0;
				im1=sindtheta*re0+cosdtheta*im0;
				re0=re1;
				im0=im1;
				
			}
		}
	}
	
	//finally, apply global phase and decay
	double cosphase=cos(phase*PI/180);
	double sinphase=sin(phase*PI/180);
	for (int t=0; t<N; t++) {
		double re1=output_real[t]*cosphase-output_imag[t]*sinphase;
		double im1=output_real[t]*sinphase+output_imag[t]*cosphase;
		re1*=exp(-(t+0.5)*dt/(1000)/d->m_T2);
		im1*=exp(-(t+0.5)*dt/(1000)/d->m_T2);
		output_real[t]=re1;
		output_imag[t]=im1;
	}
	
	double dur0=dt*N;
	QList<double> moment0;
	moment0 << A[0]*dur0;
	moment0 << A[1]*dur0;
	moment0 << A[2]*dur0;
	evolve(moment0,dur0);
}

void VmriInterfacePrivate::rotateZ(double *M,double theta) {
	double costheta=cos(theta);
	double sintheta=sin(theta);
	double ret0=costheta*M[0]-sintheta*M[1];
	double ret1=sintheta*M[0]+costheta*M[1];
	M[0]=ret0;
	M[1]=ret1;
}
void VmriInterfacePrivate::decayMag(double *M,double E1,double E2) {
	M[0]*=E2;
	M[1]*=E2;
	M[2]=1-(1-M[2])*E1;
}

double VmriInterfacePrivate::dot_product(const QList<double> &V1,double *V2) {
	return V1[0]*V2[0]+V1[1]*V2[1]+V1[2]*V2[2];
}
void VmriInterfacePrivate::apply_matrix(double *Mag,const Matrix44 &MM) {
	double ret0=Mag[0]*MM.data[0][0]+Mag[1]*MM.data[0][1]+Mag[2]*MM.data[0][2]+MM.data[0][3];
	double ret1=Mag[0]*MM.data[1][0]+Mag[1]*MM.data[1][1]+Mag[2]*MM.data[1][2]+MM.data[1][3];
	double ret2=Mag[0]*MM.data[2][0]+Mag[1]*MM.data[2][1]+Mag[2]*MM.data[2][2]+MM.data[2][3];
	Mag[0]=ret0;
	Mag[1]=ret1;
	Mag[2]=ret2;
}

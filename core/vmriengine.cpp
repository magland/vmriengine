#include "vmriengine.h"
#include <QList>
#include "vmriinterface.h"
#include <stdio.h>
#include <math.h>
#include <QDebug>
#include "qjson.h"
#include <QDateTime>

#define COMMAND_TYPE_addRFOperator 1
#define COMMAND_TYPE_initialize 2
#define COMMAND_TYPE_evolve 3
#define COMMAND_TYPE_excite 4
#define COMMAND_TYPE_readout 5

#ifndef PI
#define PI 3.141592
#endif

class VmriEnginePrivate {
public:
	AbstractPhantom *m_phantom;
	AbstractSimSequence *m_sim_sequence;
	double m_T1,m_T2;
	double m_gamma;
	
	VmriInterface m_interface;
	
	QList<VmriReadout> m_readouts;
	int m_thread_id;
	int m_readout_number;
	
	void load_isochromats();
	QList<SimCommand> load_sim_commands();
	
	Matrix44 compute_RF_operator(int N,double dt,const QList<double> &waveform_real,const QList<double> &waveform_imag,double dtheta);
	void rotateZ(double *M,double theta);
	void rotateY(double *M,double theta);
	void decayMag(double *M,double E1,double E2);
};

VmriEngine::VmriEngine() {
	d=new VmriEnginePrivate;
	d->m_phantom=0;
	d->m_sim_sequence=0;
	d->m_T1=300; //just applies to RF pulses
	d->m_T2=40; //just applies to RF pulses
	d->m_gamma=42.57;
	static int thread_id=0;
	d->m_thread_id=thread_id++;
	d->m_readout_number=0;
}

VmriEngine::~VmriEngine() {
	delete d;
}

void VmriEngine::setPhantom(AbstractPhantom *phantom) {
	d->m_phantom=phantom;
	d->m_T1=phantom->T1();
	d->m_T2=phantom->T2();
	d->m_interface.setT1(d->m_T1);
	d->m_interface.setT2(d->m_T2);
}
void VmriEngine::setSimSequence(AbstractSimSequence *sim_sequence) {
	d->m_sim_sequence=sim_sequence;
}

void VmriEnginePrivate::load_isochromats() {
	if (!m_phantom) return;
	
	qint32 N=m_phantom->isochromatCount();
	double *x,*y,*z,*f,*d;
	x=(double *)malloc(sizeof(double)*N);
	y=(double *)malloc(sizeof(double)*N);
	z=(double *)malloc(sizeof(double)*N);
	f=(double *)malloc(sizeof(double)*N);
	d=(double *)malloc(sizeof(double)*N);
	
	for (long i=0; i<N; i++) {
		Isochromat II=m_phantom->getIsochromat(i);
		
		x[i]=II.x;
		y[i]=II.y;
		z[i]=II.z;
		f[i]=II.f;
		d[i]=II.d;
	}
	
	m_interface.setIsochromats(N,x,y,z,f,d);
	
	free(x);
	free(y);
	free(z);
	free(f);
	free(d);
}

QList<VmriReadout> VmriEngine::getReadouts() {
	return d->m_readouts;
}

QList<SimCommand> VmriEnginePrivate::load_sim_commands() {
	QList<SimCommand> ret;
	
	if (!m_sim_sequence) return ret;
	
	while (!m_sim_sequence->atEnd()) {
		ret << m_sim_sequence->nextSimCommand();
	}
	
	return ret;
}

QList<double> to_double_list(const QList<QVariant> &L) {
	QList<double> ret;
	for (int i=0; i<L.count(); i++) ret << L[i].toDouble();
	return ret;
}

Matrix44 VmriEnginePrivate::compute_RF_operator(int N,double dt,const QList<double> &waveform_real,const QList<double> &waveform_imag,double dtheta) {
	
	double E1=exp(-(dt/1000)/m_T1);
	double E2=exp(-(dt/1000)/m_T2);
	double E1b=exp(-(dt/2/1000)/m_T1);
	double E2b=exp(-(dt/2/1000)/m_T2);
	
	double M0[3]; M0[0]=0; M0[1]=0; M0[2]=0;
	double M1[3]; M1[0]=1; M1[1]=0; M1[2]=0;
	double M2[3]; M2[0]=0; M2[1]=1; M2[2]=0;
	double M3[3]; M3[0]=0; M3[1]=0; M3[2]=1;
	
	rotateZ(M0,dtheta/2);
	rotateZ(M1,dtheta/2);
	rotateZ(M2,dtheta/2);
	rotateZ(M3,dtheta/2);
	
	decayMag(M0,E1b,E2b);
	decayMag(M1,E1b,E2b);
	decayMag(M2,E1b,E2b);
	decayMag(M3,E1b,E2b);
	
	for (int i=0; i<N; i++) {
		if (i>0) {
			rotateZ(M0,dtheta);
			rotateZ(M1,dtheta);
			rotateZ(M2,dtheta);
			rotateZ(M3,dtheta);
			
			decayMag(M0,E1,E2);
			decayMag(M1,E1,E2);
			decayMag(M2,E1,E2);
			decayMag(M3,E1,E2);
		}
		double re0=waveform_real[i]; //uT
		double im0=waveform_imag[i]; //uT
		double mag0=sqrt(re0*re0+im0*im0); //uT
		double phi=atan2(im0,re0); //radians
		double theta=mag0*m_gamma*dt/(1000*1000)*2*PI; //radians
		
		rotateZ(M0,phi);
		rotateZ(M1,phi);
		rotateZ(M2,phi);
		rotateZ(M3,phi);
		
		rotateY(M0,theta);
		rotateY(M1,theta);
		rotateY(M2,theta);
		rotateY(M3,theta);
		
		rotateZ(M0,-phi);
		rotateZ(M1,-phi);
		rotateZ(M2,-phi);
		rotateZ(M3,-phi);
	}
	
	rotateZ(M0,dtheta/2);
	rotateZ(M1,dtheta/2);
	rotateZ(M2,dtheta/2);
	rotateZ(M3,dtheta/2);
	
	decayMag(M0,E1b,E2b);
	decayMag(M1,E1b,E2b);
	decayMag(M2,E1b,E2b);
	decayMag(M3,E1b,E2b);
	
	Matrix44 ret;
	ret.data[0][0]=M1[0]-M0[0]; ret.data[0][1]=M2[0]-M0[0]; ret.data[0][2]=M3[0]-M0[0]; ret.data[0][3]=M0[0];
	ret.data[1][0]=M1[1]-M0[1]; ret.data[1][1]=M2[1]-M0[1]; ret.data[1][2]=M3[1]-M0[1]; ret.data[1][3]=M0[1];
	ret.data[2][0]=M1[2]-M0[2]; ret.data[2][1]=M2[2]-M0[2]; ret.data[2][2]=M3[2]-M0[2]; ret.data[2][3]=M0[2];
	ret.data[3][0]=0; ret.data[3][1]=0; ret.data[3][2]=0; ret.data[3][3]=1;
		
	return ret;
}
void VmriEnginePrivate::rotateZ(double *M,double theta) {
	double costheta=cos(theta);
	double sintheta=sin(theta);
	double ret0=costheta*M[0]-sintheta*M[1];
	double ret1=sintheta*M[0]+costheta*M[1];
	M[0]=ret0;
	M[1]=ret1;
}
void VmriEnginePrivate::rotateY(double *M,double theta) {
	double costheta=cos(theta);
	double sintheta=sin(theta);
	double ret0=costheta*M[0]+sintheta*M[2];
	double ret2=-sintheta*M[0]+costheta*M[2];
	M[0]=ret0;
	M[2]=ret2;
}
void VmriEnginePrivate::decayMag(double *M,double E1,double E2) {
	M[0]*=E2;
	M[1]*=E2;
	M[2]=1-(1-M[2])*E1;
}

void display_Matrix44(const Matrix44 &M) {
	for (int i=0; i<4; i++) {
		for (int j=0; j<4; j++) {
			printf("%.4f ",M.data[i][j]);
		}
		printf("\n");
	}
}


void VmriEngine::run() {
	
	qsrand(QDateTime::currentDateTime().toMSecsSinceEpoch());
	
	d->load_isochromats();
	QList<SimCommand> CClist=d->load_sim_commands();
	
	for (int i=0; i<CClist.count(); i++) {
		SimCommand CC=CClist[i];
		if (CC.command=="initialize") {
			d->m_interface.initialize();
		}
		else if (CC.command=="finalize") {
		}
		else if (CC.command=="add_RF_waveform") {
			double dt=CC.data.value("dt",10).toDouble();
			QList<double> data_real=to_double_list(CC.data.value("data_real").toList());
			QList<double> data_imag=to_double_list(CC.data.value("data_imag").toList());
			
			int N=data_real.count();
			
			int oversamp=4;
			double sampling_bandwidth=1.0/dt*(1000*1000); //Hz
			double freq_step=sampling_bandwidth/(N*oversamp);
			double min_freq=-freq_step*(N*oversamp/2);
			double num_freqs=N*oversamp;
			
			QList<Matrix44> operators;
			for (int ii=0; ii<N*oversamp; ii++) {
				double freq0=min_freq+ii*freq_step;
				double dtheta=freq0*dt/(1000*1000)*2*PI;
				Matrix44 MM=d->compute_RF_operator(N,dt,data_real,data_imag,dtheta);
				operators << MM;
			}
			
			d->m_interface.addRFOperator(min_freq,freq_step,num_freqs,operators);
		}
		else if (CC.command=="run_block") {
			QString block_type=CC.data.value("block_type").toString();
			if (block_type=="evolve") {
				QList<double> gradient_moment=to_double_list(CC.data.value("gradient_moment").toList());
				double duration=CC.data.value("duration",0).toDouble();
				d->m_interface.evolve(gradient_moment,duration);
			}
			else if (block_type=="rf_pulse") {
				int RF_index=CC.data.value("rf_waveform_index",0).toInt();
				QList<double> A=to_double_list(CC.data.value("gradient_amplitude").toList());
				double phase=CC.data.value("phase",0).toDouble();
				double frequency=CC.data.value("frequency",0).toDouble();
				d->m_interface.excite(RF_index,A,phase,frequency);
			}
			else if (block_type=="readout") {
				int num_readout_points=CC.data.value("N",1).toInt();
				//qDebug()  << "READOUT" << QString("thread=%1").arg(d->m_thread_id) << QString("readout=%1").arg(d->m_readout_number++) << QString("N=%1").arg(num_readout_points);
				QList<double> A=to_double_list(CC.data.value("gradient_amplitude").toList());
				double dwell_time=CC.data.value("dt",10).toDouble();
				double phase=CC.data.value("phase",0).toDouble();
				double frequency=CC.data.value("frequency",0).toDouble();
				int readout_index=CC.data.value("readout_index",0).toInt();
				
				double *ret_real=(double *)malloc(sizeof(double)*num_readout_points);
				double *ret_imag=(double *)malloc(sizeof(double)*num_readout_points);
							
				d->m_interface.readout(ret_real,ret_imag,num_readout_points,A,dwell_time,phase,frequency);
				
				quint32 N=num_readout_points;
				quint32 readout_index0=readout_index;
				float dwell_time0=(float)dwell_time;
				
				VmriReadout RO;
				RO.N=N;
				RO.readout_index=readout_index0;
				RO.dwell_time=dwell_time0;
				for (int i=0; i<N; i++) {
					RO.data_real << ret_real[i];
					RO.data_imag << ret_imag[i];
				}
				d->m_readouts << RO;
				
				free(ret_real);
				free(ret_imag);
			}
			else {
				qWarning() << "Unrecognized block_type: "+block_type;
			}
		}
		else {
			qWarning() << "Unrecognized command: "+CC.command;
		}
	}
}

#include "staticphantom.h"
#include <math.h>
#include <stdlib.h>
#include <QDebug>

class StaticPhantomPrivate {
public:
	StaticPhantom *q;
	long m_num_isochromats;
	double m_density_factor; //for isochromat density calculation
	double m_volume;
	
	void compute_volume();
};

StaticPhantom::StaticPhantom() {
	d=new StaticPhantomPrivate;
	d->q=this;
	d->m_num_isochromats=1000;
	d->m_volume=0;
	d->m_density_factor=1;
}

StaticPhantom::~StaticPhantom() {
	delete d;
}

long StaticPhantom::isochromatCount() {
	return d->m_num_isochromats;
}

void StaticPhantom::setNumIsochromats(long num) {
	d->m_num_isochromats=num;
}

void StaticPhantom::setDensityFactor(double factor) {
	d->m_density_factor=factor;
}

void StaticPhantom::initialize() {
	d->compute_volume();
}

void make_rand(double &x,double xmin,double xmax) {
	double ret=(qrand()%100000)*1.0/100000;
	x=ret*(xmax-xmin)+xmin;
}

void make_rand(double &x,double &y,double &z,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax) {
	make_rand(x,xmin,xmax);
	make_rand(y,ymin,ymax);
	make_rand(z,zmin,zmax);
}

Isochromat StaticPhantom::getIsochromat(long i) {
	
	double xmin,xmax,ymin,ymax,zmin,zmax;
	getBoundingBox(xmin,xmax,ymin,ymax,zmin,zmax);
	
	Isochromat II;
	
	bool done=false;
	while (!done) {
		double x,y,z;
		make_rand(x,y,z,xmin,xmax,ymin,ymax,zmin,zmax);
		SpinProperties SP=spinPropertiesAt(x,y,z);
		if (SP.density) {
			done=true;
			II.x=x;
			II.y=y;
			II.z=z;
			if (d->m_volume) II.d=SP.density*d->m_volume/d->m_num_isochromats*d->m_density_factor;
			else II.d=0;
			II.f=SP.chemical_shift*42.57*1.5; //ppm * Hz/uT * T
			//figure out what to do with T2star
		}
	}
	
	return II;
}
void StaticPhantomPrivate::compute_volume() {
	double xmin,xmax,ymin,ymax,zmin,zmax;
	q->getBoundingBox(xmin,xmax,ymin,ymax,zmin,zmax);
	
	long num_hits=0;
	long num_trials=0;
	long num_trials_per_region=100;
	for (int iz=0; iz<10; iz++)
	for (int iy=0; iy<10; iy++)
	for (int ix=0; ix<10; ix++) {
		double xmin0=xmin+ix*(xmax-xmin)/10;
		double xmax0=xmin+(ix+1)*(xmax-xmin)/10;
		double ymin0=ymin+iy*(ymax-ymin)/10;
		double ymax0=ymin+(iy+1)*(ymax-ymin)/10;
		double zmin0=zmin+iz*(zmax-zmin)/10;
		double zmax0=zmin+(iz+1)*(zmax-zmin)/10;
		for (long i=0; i<num_trials_per_region; i++) {
			double x0,y0,z0;
			make_rand(x0,y0,z0,xmin0,xmax0,ymin0,ymax0,zmin0,zmax0);
			SpinProperties PP=q->spinPropertiesAt(x0,y0,z0);
			if (PP.density) {
				num_hits++;
			}
			num_trials++;
		}
	}
	m_volume=num_hits*1.0/num_trials*(xmax-xmin)*(ymax-ymin)*(zmax-zmin);
	qDebug()  << "VOLUME OF PHANTOM" << m_volume;
}

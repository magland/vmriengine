#ifndef StaticPhantom_h
#define StaticPhantom_h

#include "../core/abstractphantom.h"

struct SpinProperties {
	double T2star; //ms
	double chemical_shift; //ppm
	double density;
};

class StaticPhantomPrivate;
class StaticPhantom : public AbstractPhantom {
public:
	friend class StaticPhantomPrivate;
	StaticPhantom();
	virtual ~StaticPhantom();
	
	void setNumIsochromats(long num);
	void setDensityFactor(double factor);
	void initialize();
	
	virtual long isochromatCount();
	virtual Isochromat getIsochromat(long i);
	
	virtual SpinProperties spinPropertiesAt(double x,double y,double z) const=0;
	virtual void getBoundingBox(double &xmin,double &xmax,double &ymin,double &ymax,double &zmin,double &zmax) const=0;
	
public:
	StaticPhantomPrivate *d;
};


#endif
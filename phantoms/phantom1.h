#ifndef Phantom1_h
#define Phantom1_h

#include "staticphantom.h"

class Phantom1Private;
class Phantom1 : public StaticPhantom {
public:
	friend class Phantom1Private;
	Phantom1();
	virtual ~Phantom1();
	
	virtual double T1();
	virtual double T2();
	
	virtual SpinProperties spinPropertiesAt(double x,double y,double z) const;
	virtual void getBoundingBox(double &xmin,double &xmax,double &ymin,double &ymax,double &zmin,double &zmax) const;
	
public:
	Phantom1Private *d;
};


#endif
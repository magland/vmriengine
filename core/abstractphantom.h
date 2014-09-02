#ifndef abstractphantom_h
#define abstractphantom_h

struct Isochromat {
	double x,y,z;
	double f,d;
	double T1,T2;
};

class AbstractPhantom {
public:
	AbstractPhantom() {};
	virtual ~AbstractPhantom();
	
	virtual long isochromatCount()=0;
	virtual Isochromat getIsochromat(long i)=0;
	virtual double T1()=0;
	virtual double T2()=0;
};

#endif

#ifndef sequence1_h
#define sequence1_h

#include <QString>
#include "../core/abstractsimsequence.h"

class Sequence1Private;
class Sequence1 : public AbstractSimSequence {
public:
	friend class Sequence1Private;
	Sequence1();
	virtual ~Sequence1();
	
	void loadFromFile(const QString &fname);
	
	bool atEnd();
	SimCommand nextSimCommand();
	
	
private:
	Sequence1Private *d;
};

#endif

// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: forceField.sip,v 1.5 2004/12/27 17:24:57 amoll Exp $
//

class ForceField
{
%TypeHeaderCode
	#include <BALL/MOLMEC/COMMON/forceField.h>
	#include <BALL/COMMON/exception.h>
	using namespace BALL;
	using namespace BALL::Exception;
%End
	public:
	//typedef vector<pair<Atom*, Atom*> > PairVector;

	ForceField();
	ForceField(System&);
	ForceField(System&, Options&);
	ForceField(const ForceField&);

	// ForceField& operator = (const ForceField& force_field);
	bool isValid();
	bool setup(System& /system/) throw(FileNotFound);
	bool setup(System& /system/, Options& /options/) throw(FileNotFound);
	virtual bool specificSetup() throw(TooManyErrors);
	void setName(const String&);
	String getName() const;
	int getNumberOfAtoms() const;
	int getNumberOfMovableAtoms() const;
	const AtomVector& getAtoms() const;
	System*  getSystem();
	bool getUseSelection();
	void disableSelection() throw();
	void enableSelection() throw();
	// ForceFieldParameters& getParameters();
	int countComponents() const;
	void insertComponent(ForceFieldComponent*);
	void removeComponent(const ForceFieldComponent*);
	void removeComponent(const String&);
	ForceFieldComponent* getComponent(int) const;
	ForceFieldComponent* getComponent(String&) const;
	double getEnergy() const;
	double updateEnergy();
	void  updateForces();
	double getRMSGradient() const;
	int getUpdateFrequency() const;
	void update() throw (TooManyErrors);
	Options options;
	PeriodicBoundary periodic_boundary;
};

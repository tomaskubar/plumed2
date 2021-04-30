#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>

typedef double dvec[3];

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR QMCHARGE
/*
Consider the charge of an atom that is being considered as a QM atom in a QM/MM simulation.

\par Examples
In the following, the QM zone in the simulation consists of atoms 1 through 14 as well as 29,
and the charge of atom #2 is considered for the collective variable:
the atoms from 1-10 and print it every 5 steps
\verbatim
q: QMCHARGE ATOM=2 QMATOMS=1-14,29
\endverbatim

\attention
Attention!

*/
//+ENDPLUMEDOC

class QMcharge : public Colvar {
  unsigned iQMatom;
  int QMnr;
  vector<int> QMix;
  vector<double> QMq;
  vector<double> QMdqdx;

public:
  explicit QMcharge(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
  void setQMq(void *val);
  void setQMdqdx(void *val);
};

PLUMED_REGISTER_ACTION(QMcharge,"QMCHARGE")

void QMcharge::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOM","the atom for which the QM charges is being considered");
  keys.add("atoms","QMATOMS","all of the QM atoms -- cumbersome but important!");
}

QMcharge::QMcharge(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOM", atoms);

  if(atoms.size()==1) {
    log.printf("  of QM atom no. %d\n",atoms[0].serial());
  } else error("Number of specified atoms for QMCHARGE should be 1");

  // put the atom[0] into memory!
  iQMatom = atoms[0].index();

  addValueWithDerivatives(); setNotPeriodic();

  // We need to use *all* of the QM atoms.
  vector<AtomNumber> QMatoms;
  parseAtomList("QMATOMS", QMatoms);
  QMnr = QMatoms.size();
  QMix.resize(QMnr);
  for (int i=0; i<QMnr; i++) {
    QMix[i] = QMatoms[i].index();
  }
  QMq.resize(QMnr);
  QMdqdx.resize(3*QMnr);

  requestAtoms(QMatoms);

  checkRead();
}

// calculator
void QMcharge::calculate() {

  for (int i=0; i<QMnr; i++) {
    Vector temp_dvec;
    temp_dvec[0] = QMdqdx[3*i  ];
    temp_dvec[1] = QMdqdx[3*i+1];
    temp_dvec[2] = QMdqdx[3*i+2];
    setAtomsDerivatives(i, temp_dvec);
  }
  setValue(QMq[iQMatom]);
}

// Grab the charges and charge derivatives
//   already on the level of a command in mdrun.cpp.
//   We will iterate over all of the actions, and whenever
//   the action is QMcharges, we will obtain the data.
//   The data structures of type QMcharge (inherited from
//   Action) are accessible from there.

void QMcharge::setQMq(void *val) {
  for (int i=0; i<QMnr; i++) {
    QMq[i] = ((double *) val)[i];
  }
}

void QMcharge::setQMdqdx(void *val) {
  // fill it with content
  dvec *source = (dvec *) val;
  for (int i=0; i<QMnr; i++) {
    QMdqdx[3*i  ] = source[iQMatom * QMnr + i][0];
    QMdqdx[3*i+1] = source[iQMatom * QMnr + i][1];
    QMdqdx[3*i+2] = source[iQMatom * QMnr + i][2];
  }
}

}
}

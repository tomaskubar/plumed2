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
atoms nos. 35 through 75 are considered as MM atoms affecting the QM charges,
and the charge of atom #2 constitutes the collective variable:
\verbatim
q: QMCHARGE ATOM=2 QMATOMS=1-14,29 MMATOMS=35-75
\endverbatim

\attention
Attention!

*/
//+ENDPLUMEDOC

class QMcharge : public Colvar {
  unsigned iQMatom;
  int QMnr;
  int MMnr;
  vector<int> QMix;
  vector<double> QMq;
  vector<double> QMdqdx;
  vector<int> MMix;
  vector<double> QMdqdxMM;

public:
  explicit QMcharge(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
  void setQMq(void *val);
  void setQMdqdx(void *val);
  void setQMdqdxMM(void *val);
};

PLUMED_REGISTER_ACTION(QMcharge,"QMCHARGE")

void QMcharge::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOM","the atom for which the QM charge is being considered");
  keys.add("atoms","QMATOMS","all of the QM atoms -- cumbersome but important!");
  keys.add("atoms","MMATOMS","the MM atoms considered to be affecting the QM charges -- cumbersome but important!");
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

  // The MM atoms that are considered to affect the QM charges.
  vector<AtomNumber> MMatoms;
  parseAtomList("MMATOMS", MMatoms);
  MMnr = MMatoms.size();
  MMix.resize(MMnr);
  for (int i=0; i<MMnr; i++) {
    MMix[i] = MMatoms[i].index();
  }
  QMdqdxMM.resize(3*MMnr);

  // it seems that all of the atoms (QM+MM) need to be requested at the same time,
  //   so we first append MMatoms to QMatoms, and then request the whole vector
  QMatoms.insert(QMatoms.end(), MMatoms.begin(), MMatoms.end());
  requestAtoms(QMatoms);
  // this will then be unnecessary
//requestAtoms(MMatoms);

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
  for (int i=0; i<MMnr; i++) {
    Vector temp_dvec;
    temp_dvec[0] = QMdqdxMM[3*i  ];
    temp_dvec[1] = QMdqdxMM[3*i+1];
    temp_dvec[2] = QMdqdxMM[3*i+2];
    setAtomsDerivatives(QMnr + i, temp_dvec);
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

void QMcharge::setQMdqdxMM(void *val) {
  // fill it with content
  dvec *source = (dvec *) val;
  for (int i=0; i<MMnr; i++) {
    QMdqdxMM[3*i  ] = source[iQMatom * MMnr + i][0];
    QMdqdxMM[3*i+1] = source[iQMatom * MMnr + i][1];
    QMdqdxMM[3*i+2] = source[iQMatom * MMnr + i][2];
  }
}

}
}

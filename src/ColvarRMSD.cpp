#include "ColvarRMSD.h"
#include "PlumedMain.h"
#include "ActionRegister.h"
#include "PDB.h"
#include "RMSD.h"
#include "Atoms.h"


using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR RMSD
/*
Calculate the RMSD with respect to a reference structure.  To perform
an ??optimal?? (what does this mean algorithmical speed wise?) alignment
using the Kearsley algorithm then use TYPE=OPTIMAL.  Otherwise
use TYPE=SIMPLE, which will not perform optimal alignment and will only 
remove the translation of the center of mass.

\attention
The documentation here needs some work as it is not very clear to me 
sorry GAT. 

\par Examples

The following tells plumed to calculate the RMSD distance between
the positions of the atoms in the reference file and their instantaneous
position.  The Kearseley algorithm is used so this is done optimally.

\verbatim
RMSD REFERENCE=file.pdb TYPE=OPTIMAL
\endverbatim

...

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(ColvarRMSD,"RMSD")

void ColvarRMSD::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV. " + PDB::documentation() );
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.addFlag("SQUARED",false," This should be setted if you want MSD instead of RMSD ");
}

ColvarRMSD::ColvarRMSD(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),rmsd(log),squared(false)
{
  string reference;
  parse("REFERENCE",reference);
  string type;	
  type.assign("SIMPLE");
  parse("TYPE",type);
  parseFlag("SQUARED",squared);

  checkRead();


  addValueWithDerivatives(); setNotPeriodic();
  PDB pdb;

  // read everything in ang and transform to nm if we are not in natural units
  pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().length);

  rmsd.set(pdb,type);

  requestAtoms(pdb.getAtomNumbers());

  derivs.resize(getNumberOfAtoms());
  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %d atoms\n",getNumberOfAtoms());
  log.printf("  method for alignment : %s \n",rmsd.getMethod().c_str() );
  if(squared)log.printf("  chosen to use SQARED option for MSD instead of RMSD\n");

}


// calculator
void ColvarRMSD::calculate(){
  double r=rmsd.calculate(getPositions(),derivs,squared);
  setValue(r);
  for(unsigned i=0;i<derivs.size();i++) setAtomsDerivatives(i,derivs[i]);
  Tensor virial;
  for(unsigned i=0;i<derivs.size();i++) virial=virial+(-1.0*Tensor(getPosition(i),derivs[i]));
  setBoxDerivatives(virial);
}

}




/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace adjmat {

class Sprint : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Sprint(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Sprint,"SPRINT")

void Sprint::registerKeywords(Keywords& keys){
  ActionShortcut::registerKeywords( keys );
  keys.add("numbered","GROUP","specifies the list of atoms that should be assumed indistinguishable");
  keys.add("numbered","SWITCH","specify the switching function to use between two sets of indistinguishable atoms");
}

Sprint::Sprint(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  std::vector<std::string> grp_str; std::string grp_inpt;
  for(unsigned i=1;; ++i) {
      if( !parseNumbered("GROUP",i,grp_inpt) ) break; 
      grp_str.push_back( grp_inpt );
  }
  if( grp_str.size()>9 ) error("cannot handle more than 9 groups");
  std::vector<unsigned> nin_group; unsigned ntot_atoms=0;
  for(unsigned i=0;i<grp_str.size(); ++i) {
    std::string sw_str, num; Tools::convert( i+1, num ); parseNumbered("SWITCH", (i+1)*10 + 1 + i,  sw_str );
    if( sw_str.length()==0 ) error("missing SWITCH" + num + num + " keyword");
    readInputLine( getShortcutLabel() + "_mat" + num +  num + ": CONTACT_MATRIX GROUP=" + grp_str[i] + " SWITCH={" + sw_str + "}" );
    // Get number of atoms in each group
    std::vector<std::string> words=Tools::getWords(grp_str[i],"\t\n ,"); Tools::interpretRanges(words); 
    nin_group.push_back( words.size() ); ntot_atoms += words.size();
    for(unsigned j=0; j<i; ++j) {
        std::string sw_str2, jnum; Tools::convert( j+1, jnum ); parseNumbered("SWITCH", (j+1)*10 + 1 + i, sw_str2);
        if( sw_str2.length()==0 ) error("missing SWITCH" + jnum + num + " keyword");
        readInputLine( getShortcutLabel() + "_mat" + jnum + num + ": CONTACT_MATRIX GROUPA=" + grp_str[j] + " GROUPB=" + grp_str[i] + " SWITCH={" + sw_str2 +"}");
        readInputLine( getShortcutLabel() + "_mat" + num +  jnum + ": TRANSPOSE ARG=" + getShortcutLabel() + "_mat" + jnum + num + ".w");
    }
  }
  std::string join_matrices = getShortcutLabel() + "_jmat: COMBINE_MATRICES"; 
  for(unsigned i=0; i<grp_str.size(); ++i) {
    std::string inum; Tools::convert(i+1,inum);
    for(unsigned j=0; j<grp_str.size(); ++j) {
      std::string jnum; Tools::convert(j+1,jnum);
      if( i>j ) join_matrices += " MATRIX" + inum + jnum + "=" + getShortcutLabel() + "_mat" + inum +  jnum;
      else join_matrices += " MATRIX" + inum + jnum + "=" + getShortcutLabel() + "_mat" + inum +  jnum + ".w";
    }
  }
  readInputLine( join_matrices );
  // Diagonalization
  readInputLine( getShortcutLabel() + "_diag: DIAGONALIZE ARG=" + getShortcutLabel() + "_jmat VECTORS=1");
  // Compute sprint coordinates as product of eigenvalue and eigenvector times square root of number of atoms in all groups
  std::string str_natoms; Tools::convert( ntot_atoms, str_natoms );
  readInputLine( getShortcutLabel() + "_sp: MATHEVAL ARG1=" + getShortcutLabel() + "_diag.vals-1 ARG2=" + getShortcutLabel() + 
                 "_diag.vecs-1 FUNC=sqrt(" + str_natoms + ")*x*y PERIODIC=NO"); 
  // Sort sprint coordinates for each group of atoms
  unsigned k=0;
  for(unsigned j=0; j<nin_group.size(); ++j) {
    std::string jnum, knum; Tools::convert( j+1, jnum ); Tools::convert(k+1, knum); k++;
    std::string sort_act = getShortcutLabel() + jnum + ": SORT ARG=" + getShortcutLabel() + "_sp." + knum;
    for(unsigned n=1; n<nin_group[j]; ++n) {
      Tools::convert( k+1, knum ); sort_act += "," + getShortcutLabel() + "_sp." + knum; k++;
    }
    readInputLine( sort_act );
  }
}

}
}

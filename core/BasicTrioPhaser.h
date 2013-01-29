
/* This file is made available under the terms of the GNU GPL 2.0 */

#ifndef __BASIC_TRIO_PHASER_H__
#define __BASIC_TRIO_PHASER_H__

/* This function tries to phase a set of alleles observed in a parent-child
   trio. It returns true if successful and rearranges alleles so that the
   father and mother's transmitted alleles are assigned to the slots
   father1 and mother1. The function will generate an assertion failure if
   the genotypes are not Mendelian consistent */

bool PhaseTrioTransmissions(int & father1, int & father2,
                           int & mother1, int & mother2,
                           int & child1, int & child2);

#endif


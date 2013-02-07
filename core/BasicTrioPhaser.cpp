
/* This file is made available under the terms of the GNU GPL 2.0 */

#include "BasicTrioPhaser.h"

#include <stdio.h>

bool PhaseTrioTransmissions(int & father1, int & father2,
                           int & mother1, int & mother2,
                           int & child1, int & child2)
   {
   /* List up to four potential trio configurations, conditional on parental genotypes */
   int phasedFather[4][2];
   int phasedMother[4][2];
   int phasedChild[4][2];

   for (int i = 0; i < 4; i++)
      {
      phasedFather[i][0] = (i & 1) ? father1 : father2;
      phasedFather[i][1] = (i & 1) ? father2 : father1;
      phasedChild[i][0] = phasedFather[i][0];
      }

   int cycle = father1 == father2 ? 1 : 2;

   for (int i = 0; i < 4; i++)
      {
      phasedMother[i][0] = (i & cycle) ? mother1 : mother2;
      phasedMother[i][1] = (i & cycle) ? mother2 : mother1;
      phasedChild[i][1] = phasedMother[i][0];
      }

   /* Count the number of distinct phases, ignoring homozygous parents */

   int possiblePhases = 1;

   if (father1 != father2)
      possiblePhases *= 2;

   if (mother1 != mother2)
      possiblePhases *= 2;

   /* Search for phases that match the observed child genotypes */
   int matchingPhase = -1;

   for (int i = 0; i < possiblePhases; i++)
      {
      if (child1 == phasedChild[i][0] && child2 == phasedChild[i][1] ||
          child1 == phasedChild[i][1] && child2 == phasedChild[i][0])
          if (matchingPhase == -1)
            matchingPhase = i;
          else
            /* Multiple different phases match observed genotype, can't be phased */
            return false;
      }

   /* If there is no matching phase, we have a Mendelian inconsistency */
   if (matchingPhase < 0)
      {
      // These calls often occur when one makes a marginal call, which picks
      // the best genotype for each individual ...

      // printf("Odd situation: %d/%d x %d/%d = %d/%d\n", father1, father2,
      //                                                  mother1, mother2,
      //                                                  child1, child2);
      return false;
      }

   /* A single phase matches observed genotypes, trio is phased */
   father1 = phasedFather[matchingPhase][0];
   father2 = phasedFather[matchingPhase][1];
   mother1 = phasedMother[matchingPhase][0];
   mother2 = phasedMother[matchingPhase][1];
   child1 = phasedChild[matchingPhase][0];
   child2 = phasedChild[matchingPhase][1];

   return true;
   }

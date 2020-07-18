#pragma once

#include "../utils/FromBias.h"

real NormalizeAngle(real);
void ScaleTo2Pi(interval&);
void iSin(interval&, interval&);
void iCos(interval&, interval&);
void iTan(interval&, interval&);
void iCot(interval&, interval&);
void iExp(interval&, interval&);
void iLog (interval&, interval&);
interval iSqr (const interval&);
interval iSqrt (const interval&);

/*BiasPi         = ;
  BiasTwoPi      = ;
  BiasPiHalf     = ;
  BiasPiQuarter  = ;
  BiasE          = ;
  BiasSqrt2      = ;
  BiasInvSqrt2   = ;
  BiasLn10       = ;*/

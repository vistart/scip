/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_sol.h,v 1.2 2004/02/04 17:27:46 bzfpfend Exp $"

/**@file   struct_sol.h
 * @brief  datastructures for storing primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __STRUCT_SOL_H__
#define __STRUCT_SOL_H__


#include "def.h"
#include "type_misc.h"
#include "type_sol.h"
#include "type_heur.h"



/** primal CIP solution
 *  For reasons of efficiency, a working solution only stores values that have been accessed at least once,
 *  or that have been changed from the value in the solution's source.
 *  The user has to call SCIPsolUnlink() in order to retrieve all non-cached elements from the solution's source
 *  and to store the values in the solution's own array. This changes the solution's origin to SCIP_SOLORIGIN_ZERO.
 *  A linked solution with origin SCIP_SOLORIGIN_LPSOL or SCIP_SOLORIGIN_PSEUDOSOL becomes invalid after the
 *  next node is activated (i.e. the LP and pseudo solutions changed) and cannot be accessed anymore.
 */
struct Sol
{
   REALARRAY*       vals;               /**< solution values for variables */
   BOOLARRAY*       valid;              /**< for solutions originating from LPSOL or PSEUDOSOL: TRUE iff variable's val
                                         *   is valid; otherwise the value has to be retrieved from the origin */
   HEUR*            heur;               /**< heuristic that found the solution (or NULL if it's an LP solution) */
   Real             obj;                /**< objective value of solution */
   Real             time;               /**< clock time, when the solution was discovered */
   Longint          nodenum;            /**< last node number, where this solution was modified */
   SOLORIGIN        solorigin;          /**< origin of solution: where to retrieve uncached elements */
   int              depth;              /**< depth at which the solution was found */
};


#endif

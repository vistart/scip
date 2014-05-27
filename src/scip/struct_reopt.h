/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_reopt.h
 * @brief  datastructures for collecting reoptimization information
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_REOPT_H__
#define __SCIP_STRUCT_REOPT_H__


#include "scip/def.h"
#include "scip/type_reopt.h"

#ifdef __cplusplus
extern "C" {
#endif

/** primal data and solution storage */
struct SCIP_Reopt
{
   SCIP_SOL***           sols;               /**< solutions of the reoptimization runs */

   int                   run;                /**< current position in the sols array*/
   int                   runsize;            /**< allocated memory for runs */
   int*                  solssize;           /**< size of sols[x] arrays */
   int*                  nsols;              /**< number of solutions stored in sols[x] array */

   SCIP_Bool*            solsused;           /**< True or False if the solutions in run x were used at least once */

   SCIP_HASHMAP*         varnamehash;        /**< hashmap which hashed varnames to indices */
   SCIP_Real**           objs;               /**< list of objective coefficients */
};

#ifdef __cplusplus
}
#endif

#endif

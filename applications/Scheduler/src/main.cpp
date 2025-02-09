/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   Scheduler/src/main.cpp
 * @brief  Main file for C++ compilation
 * @author Stefan Heinz
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/scipshell.h"

#include "cons_optcumulative.h"
#include "heur_optcumulative.h"
#include "heur_listscheduling.h"
#include "reader_cmin.h"
#include "reader_sch.h"
#include "reader_sm.h"
#include "reader_rcp.h"

/** runs the shell */
static
SCIP_RETCODE runShell(
   int                        argc,               /**< number of shell parameters */
   char**                     argv,               /**< array with shell parameters */
   const char*                defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* we explicitly enable the use of a debug solution for this main SCIP instance */
   SCIPenableDebugSol(scip);

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include problem reader */
   SCIP_CALL( SCIPincludeReaderCmin(scip) );
   SCIP_CALL( SCIPincludeReaderSch(scip) );
   SCIP_CALL( SCIPincludeReaderSm(scip) );
   SCIP_CALL( SCIPincludeReaderRcp(scip) );

   /* include problem specific heuristic */
   SCIP_CALL( SCIPincludeHeurListScheduling(scip) );
   SCIP_CALL( SCIPincludeHeurOptcumulative(scip) );

   /* include cumulative constraint handler with optional activities */
   SCIP_CALL( SCIPincludeConshdlrOptcumulative(scip) );

#ifdef WITH_CPOPTIMIZER
   SCIP_CALL( SCIPsetSolveCumulative(scip, cpoptimizer) );
#endif

   /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   /* check block memory */
   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main method */
int main(
   int                   argc,          /**< number of arguments */
   char**                argv           /**< string array with arguments */
   )
{
  SCIP_RETCODE retcode;

  retcode = runShell(argc, argv, "scip.set");

  if( retcode != SCIP_OKAY )
  {
     SCIPprintError(retcode);
     return -1;
  }

  return 0;
}

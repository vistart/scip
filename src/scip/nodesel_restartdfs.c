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
#pragma ident "@(#) $Id: nodesel_restartdfs.c,v 1.12 2004/04/19 17:08:37 bzfpfend Exp $"

/**@file   nodesel_restartdfs.c
 * @brief  node selector for depth first search with periodical selection of the best node
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "nodesel_restartdfs.h"


#define NODESEL_NAME             "restartdfs"
#define NODESEL_DESC             "depth first search with periodical selection of the best node"
#define NODESEL_STDPRIORITY       50000
#define NODESEL_MEMSAVEPRIORITY   50000
#define NODESEL_LOWESTFIRST       FALSE   /**< are the nodes sorted such that the lowest bound node comes first? */




/*
 * Default parameter settings
 */

#define SELECTBESTFREQ             1000 /**< frequency for selecting the best node instead of the deepest one */



/** node selector data for best first search node selection */
struct NodeselData
{
   int              selectbestfreq;     /**< frequency for selecting the best node instead of the deepest one */
};




/*
 * Callback methods
 */

static
DECL_NODESELFREE(nodeselFreeRestartdfs)
{  /*lint --e{715}*/
   NODESELDATA* nodeseldata;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   /* free user data of node selector */
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);
   SCIPfreeMemory(scip, &nodeseldata);
   SCIPnodeselSetData(nodesel, nodeseldata);

   return SCIP_OKAY;
}

static
DECL_NODESELSELECT(nodeselSelectRestartdfs)
{  /*lint --e{715}*/
   NODESELDATA* nodeseldata;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);
   assert(selnode != NULL);

   /* get node selector user data */
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert(nodeseldata != NULL);

   /* decide if we want to select the node with lowest bound or the deepest node */
   if( nodeseldata->selectbestfreq >= 1 && SCIPgetNodenum(scip) % nodeseldata->selectbestfreq == 0 )
      *selnode = SCIPgetBestboundNode(scip);
   else
      *selnode = SCIPgetBestNode(scip);

   return SCIP_OKAY;
}

static
DECL_NODESELCOMP(nodeselCompRestartdfs)
{  /*lint --e{715}*/
   int depth1;
   int depth2;

   assert(nodesel != NULL);
   assert(strcmp(SCIPnodeselGetName(nodesel), NODESEL_NAME) == 0);
   assert(scip != NULL);

   depth1 = SCIPnodeGetDepth(node1);
   depth2 = SCIPnodeGetDepth(node2);
   if( depth1 > depth2 )
      return -1;
   else if( depth1 < depth2 )
      return +1;
   else
   {
      Real lowerbound1;
      Real lowerbound2;

      /**@todo create a better node selection rule for nodes at equal depth (to find primal solutions quickly) */
      lowerbound1 = SCIPnodeGetLowerbound(node1);
      lowerbound2 = SCIPnodeGetLowerbound(node2);
      if( lowerbound1 < lowerbound2 )
         return -1;
      else if( lowerbound1 > lowerbound2 )
         return +1;
      else
         return 0;
   }
}





/*
 * restartdfs specific interface methods
 */

/** creates the node selector for restarting depth first search and includes it in SCIP */
RETCODE SCIPincludeNodeselRestartdfs(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   NODESELDATA* nodeseldata;

   /* allocate and initialize node selector data; this has to be freed in the destructor */
   CHECK_OKAY( SCIPallocMemory(scip, &nodeseldata) );
   nodeseldata->selectbestfreq = SELECTBESTFREQ;

   /* include node selector */
   CHECK_OKAY( SCIPincludeNodesel(scip, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
                  NODESEL_LOWESTFIRST,
                  nodeselFreeRestartdfs, NULL, NULL, nodeselSelectRestartdfs, nodeselCompRestartdfs,
                  nodeseldata) );

   /* add node selector parameters */
   CHECK_OKAY( SCIPaddIntParam(scip,
                  "nodeselection/restartdfs/selectbestfreq",
                  "frequency for selecting the best node instead of the deepest one (0: never)",
                  &nodeseldata->selectbestfreq, SELECTBESTFREQ, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}


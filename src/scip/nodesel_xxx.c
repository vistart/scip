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
#pragma ident "@(#) $Id: nodesel_xxx.c,v 1.4 2004/02/04 17:27:29 bzfpfend Exp $"

/**@file   nodesel_xxx.c
 * @brief  xxx node selector
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "nodesel_xxx.h"


#define NODESEL_NAME            "xxx"
#define NODESEL_DESC            "node selector template"
#define NODESEL_STDPRIORITY     0
#define NODESEL_MEMSAVEPRIORITY 0
#define NODESEL_LOWESTFIRST     FALSE   /**< are the nodes sorted such that the lowest bound node comes first? */




/*
 * Data structures
 */

/* TODO: fill in the necessary node selector data */

/** node selector data */
struct NodeselData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of node selector
 */

/* TODO: Implement all necessary node selector methods. The methods with an #if 0 ... #else #define ... are optional */

/** destructor of node selector to free user data (called when SCIP is exiting) */
#if 0
static
DECL_NODESELFREE(nodeselFreeXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx node selector not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselFreeXxx NULL
#endif


/** initialization method of node selector (called when problem solving starts) */
#if 0
static
DECL_NODESELINIT(nodeselInitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx node selector not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselInitXxx NULL
#endif


/** deinitialization method of node selector (called when problem solving exits) */
#if 0
static
DECL_NODESELEXIT(nodeselExitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx node selector not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselExitXxx NULL
#endif


/** node selection method of node selector */
static
DECL_NODESELSELECT(nodeselSelectXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx node selector not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** node comparison method of node selector */
static
DECL_NODESELCOMP(nodeselCompXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx node selector not implemented yet\n");
   abort(); /*lint --e{527}*/

   return 0;
}




/*
 * node selector specific interface methods
 */

/** creates the xxx node selector and includes it in SCIP */
RETCODE SCIPincludeNodeselXxx(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   NODESELDATA* nodeseldata;

   /* create xxx node selector data */
   nodeseldata = NULL;
   /* TODO: (optional) create node selector specific data here */

   /* include node selector */
   CHECK_OKAY( SCIPincludeNodesel(scip, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
                  NODESEL_LOWESTFIRST,
                  nodeselFreeXxx, nodeselInitXxx, nodeselExitXxx, nodeselSelectXxx, nodeselCompXxx,
                  nodeseldata) );

   /* add xxx node selector parameters */
   /* TODO: (optional) add node selector specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

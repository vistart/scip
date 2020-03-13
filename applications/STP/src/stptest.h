/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   stptest.h
 * @brief  includes various testing methods for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_STPTEST_H_
#define APPLICATIONS_STP_SRC_STPTEST_H_

#include <stdlib.h>
#include "graph.h"
#include "reduce.h"
#include "scip/scip.h"

#define STPTEST_ASSERT_MSG_ARGS(cond, msg, ...)  \
   do                                            \
   {                                             \
      if (!(cond))                               \
      {                                          \
         printf("[%s:%d] unit test failed: ", __FILE__, __LINE__); \
         printf((msg), __VA_ARGS__);             \
         abort();                                \
      }                                          \
   } while(0)


#define STPTEST_ASSERT_MSG(cond, msg)  \
   do                                  \
   {                                   \
      if (!(cond))                     \
      {                                \
         printf("[%s:%d] unit test failed: ", __FILE__, __LINE__); \
         printf((msg));                \
         abort();                      \
      }                                \
   } while(0)


#define STPTEST_ASSERT(cond)           \
   do                                  \
   {                                   \
      if (!(cond))                     \
      {                                \
         printf("[%s:%d] unit test failed: ", __FILE__, __LINE__); \
         printf("\n");                 \
         abort();                      \
      }                                \
   } while(0)


/* stptest_base.c
 */
extern SCIP_RETCODE    stptest_testAll(SCIP*);

/* stptest_graph.c
 */
extern SCIP_RETCODE    stptest_csrdepo(SCIP*);
extern void            stptest_graphTearDown(SCIP*, GRAPH*);
extern SCIP_RETCODE    stptest_graphSetUp(SCIP*, GRAPH*);
extern SCIP_RETCODE    stptest_graphSetUpPcOrg(SCIP*, GRAPH*, int*, int*);
extern SCIP_RETCODE    stptest_graphSetUpPcExtended(SCIP*, GRAPH*, int*, int*);

/* stptest_extreduce.c
 */
extern SCIP_RETCODE    stptest_extreduce(SCIP*);
extern void            stptest_extreduceTearDown(SCIP*, GRAPH*, REDCOST*);

/* stptest_extutils.c
 */
extern SCIP_RETCODE    stptest_extmldists(SCIP*);

/* stptest_reduce.c
 */
extern SCIP_RETCODE    stptest_dcmst(SCIP*);
extern SCIP_RETCODE    stptest_reduce_sdpcmw(SCIP*);

/* stptest_misc.c
 */
extern SCIP_RETCODE    stptest_pseudoAncestors(SCIP*);
extern SCIP_RETCODE    stptest_pseudoDel(SCIP*);
extern SCIP_RETCODE    stptest_dheap(SCIP*);
extern SCIP_RETCODE    stptest_completegraph(SCIP*);

/* stptest_heurlocal.c
 */
extern SCIP_RETCODE    stptest_testHeurLocal(SCIP*);


/* stptest_heurtm.c
 */
extern SCIP_RETCODE    stptest_testSolPrune(SCIP*);
extern SCIP_RETCODE    stptest_testHeurTm(SCIP*);


#endif /* APPLICATIONS_STP_SRC_STPTEST_H_ */

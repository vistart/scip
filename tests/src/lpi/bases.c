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

/**@file   bases.c
 * @brief  unit test for checking the settings of slack variables in a basis of the lpi
 * @author Marc Pfetsch
 * @author Franziska Schloesser
 * @author Felipe Serrano
 *
 * The behavior of different LP solvers w.r.t. the slack variables should not differ, if interfaced by LPI.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <scip/scip.h>
#include <lpi/lpi.h>
#include "include/scip_test.h"
#include "stdio.h"

#define EPS 1e-6

/* global variable for LPI */
static SCIP_LPI* lpi;


/*** TEST SUITE SIMPLE ***/
static
void setup_simple(void)
{
   int nrows;
   int ncols;
   int beg = 0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 3.0;
   SCIP_Real lhs = 1.0;
   SCIP_Real rhs = 2.0;
   SCIP_Real obj = 1.0;
   SCIP_Real val = 1.0;
   int ind = 0;

   lpi = NULL;

   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

   /* use the following LP as base:
    *   max x
    *       1 <= x <= 2  (linear constraint)
    *       0 <= x <= 3  (bounds)
    */
   /* add one column */
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

   /* add one row */
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 1, &beg, &ind, &val) );

   /* check size */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   cr_assert( nrows == 1 );
   cr_assert( ncols == 1 );

#ifdef SCIP_DEBUG
   /* turn on output */
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPINFO, 1) );
#endif
}

static
void teardown(void)
{
   SCIP_CALL( SCIPlpiFree(&lpi) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}
TestSuite(simple, .init = setup_simple, .fini = teardown);

/*** TESTS ***/
Test(simple, test1)
{
   int cstat;
   int rstat;

   /* solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   /* get basis */
   SCIP_CALL( SCIPlpiGetBase(lpi, &cstat, &rstat) );

   /* the variable should be basic and the slack variable at the upper bound */
   assert( cstat == SCIP_BASESTAT_BASIC );
   assert( rstat == SCIP_BASESTAT_UPPER );
}

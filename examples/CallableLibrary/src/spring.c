/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   spring.c
 * @brief  Coil Compression Spring Design model
 * @author Stefan Vigerske
 *
 * This example shows how to setup quadratic and nonlinear constraints in SCIP when using SCIP as callable library.
 * The example implements a model for the design of a coil compression spring as it can be found in the GAMS model library:
 * https://www.gams.com/latest/gamslib_ml/libhtml/gamslib_spring.html
 *
 * The task is to find a minimum volume of a wire for the production of a coil compression spring.
 *
 * Original model source:
 * @par
 *    E. Sangren@n
 *    Nonlinear Integer and Discrete Programming in Mechanical Design Optimization@n
 *    Journal of Mechanical Design, Trans. ASME 112 (1990), 223-229
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <math.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_var.h"

#ifndef M_PI
#define M_PI           3.141592653589793238462643
#endif


/* Model parameters */

/** number of possible wire types */
#define nwires 11

/** diameters of available diameters (in) */
static const SCIP_Real diameters[] = { 0.207, 0.225, 0.244, 0.263, 0.283, 0.307, 0.331, 0.362, 0.394, 0.4375, 0.500 };

/** preload (lb) */
static const SCIP_Real preload = 300;

/** maximal working load (lb) */
static const SCIP_Real maxworkload = 1000;

/** maximal deflection (in) */
static const SCIP_Real maxdeflect = 6;

/** deflection from preload (in) */
static const SCIP_Real deflectpreload = 1.25;

/** maximal free length of spring (in) */
static const SCIP_Real maxfreelen = 14.0;

/** maximal coil diameter (in) */
static const SCIP_Real maxcoildiam = 3.0;

/** maximal shear stress */
static const SCIP_Real maxshearstress = 189000.0;

/** shear modulus of material */
static const SCIP_Real shearmod = 11500000.0;


/** sets up problem */
static
SCIP_RETCODE setupProblem(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_VAR* coil;        /* coil diameter */
   SCIP_VAR* wire;        /* wire diameter */
   SCIP_VAR* defl;        /* deflection */
   SCIP_VAR* ncoils;      /* number of coils (integer) */
   SCIP_VAR* const1;      /* a constant */
   SCIP_VAR* const2;      /* another constant */
   SCIP_VAR* volume;      /* total volume */
   SCIP_VAR* y[nwires];   /* wire choice (binary) */

   SCIP_CONSHDLR* consexprhdlr;
   SCIP_CONSEXPR_EXPR* coilexpr;
   SCIP_CONSEXPR_EXPR* wireexpr;
   SCIP_CONSEXPR_EXPR* deflexpr;
   SCIP_CONSEXPR_EXPR* ncoilsexpr;
   SCIP_CONSEXPR_EXPR* const1expr;
   SCIP_CONSEXPR_EXPR* const2expr;
   SCIP_CONSEXPR_EXPR* volumeexpr;

   SCIP_CONS* voldef;
   SCIP_CONS* defconst1;
   SCIP_CONS* defconst2;
   SCIP_CONS* shear;
   SCIP_CONS* defdefl;
   SCIP_CONS* freel;
   SCIP_CONS* coilwidth;
   SCIP_CONS* defwire;
   SCIP_CONS* selectwire;

   char name[SCIP_MAXSTRLEN];
   int i;

   /* find expression constraint handler */
   consexprhdlr = SCIPfindConshdlr(scip, "expr");
   assert(consexprhdlr != NULL);

   /* create empty problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "spring") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &coil, "coildiam", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &wire, "wirediam", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &defl, "deflection", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &ncoils, "ncoils", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &const1, "const1", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &const2, "const2", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &volume, "volume", 0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   for( i = 0; i < nwires; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "wire%d", i+1);
      SCIP_CALL( SCIPcreateVarBasic(scip, &y[i], name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   }

   /* set nonstandard variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, defl, deflectpreload / (maxworkload - preload)) );
   SCIP_CALL( SCIPchgVarUb(scip, defl, maxdeflect / preload) );

   /* add variables to problem */
   SCIP_CALL( SCIPaddVar(scip, coil) );
   SCIP_CALL( SCIPaddVar(scip, wire) );
   SCIP_CALL( SCIPaddVar(scip, defl) );
   SCIP_CALL( SCIPaddVar(scip, ncoils) );
   SCIP_CALL( SCIPaddVar(scip, const1) );
   SCIP_CALL( SCIPaddVar(scip, const2) );
   SCIP_CALL( SCIPaddVar(scip, volume) );
   for( i = 0; i < nwires; ++i )
   {
      SCIP_CALL( SCIPaddVar(scip, y[i]) );
   }

   /* create variable expressions */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &coilexpr, coil) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &wireexpr, wire) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &deflexpr, defl) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &ncoilsexpr, ncoils) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &const1expr, const1) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &const2expr, const2) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &volumeexpr, volume) );

   /* nonlinear constraint voldef: PI/2 * (ncoils+2)*coil*wire^2 - volume == 0 */
   {
      SCIP_CONSEXPR_EXPR* exprs[3];
      SCIP_CONSEXPR_EXPR* powexpr;
      SCIP_CONSEXPR_EXPR* prodexpr;
      SCIP_CONSEXPR_EXPR* sumexpr;
      SCIP_CONSEXPR_EXPR* expr;
      SCIP_Real coefs[2];

      /* create wire^2 */
      SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, &powexpr, wireexpr, 2.0) );

      /* create (ncoils+2) */
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, &sumexpr, 1, &ncoilsexpr, NULL, 2.0) );

      /* create (ncoils+2)*coil*wire^2 */
      exprs[0] = sumexpr;
      exprs[1] = coilexpr;
      exprs[2] = powexpr;
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prodexpr, 3, exprs, 1.0) );

      /* create PI/2 * (ncoils+2)*coil*wire^2 - volume */
      exprs[0] = prodexpr;
      coefs[0] = M_PI / 2.0;
      exprs[1] = volumeexpr;
      coefs[1] = -1.0;
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, &expr, 2, exprs, coefs, 0.0) );

      /* create nonlinear constraint */
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &voldef, "voldef", expr, 0.0, 0.0) );

      /* release expressions */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexpr) );
   }

   /* nonlinear constraint defconst1: coil / wire - const1 == 0.0 */
   {
      SCIP_CONSEXPR_EXPR* exprs[2];
      SCIP_CONSEXPR_EXPR* powexpr;
      SCIP_CONSEXPR_EXPR* prodexpr;
      SCIP_CONSEXPR_EXPR* sumexpr;
      SCIP_Real coefs[2];

      /* create 1 / wire */
      SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, &powexpr, wireexpr, -1.0) );

      /* create coil / wire */
      exprs[0] = coilexpr;
      exprs[1] = powexpr;
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prodexpr, 2, exprs, 1.0) );

      /* create coil / wire - const1 */
      exprs[0] = prodexpr;
      coefs[0] = 1.0;
      exprs[1] = const1expr;
      coefs[1] = -1.0;
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, &sumexpr, 2, exprs, coefs, 0.0) );

      /* create nonlinear constraint */
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &defconst1, "defconst1", sumexpr, 0.0, 0.0) );

      /* release expressions */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexpr) );
   }

   /* nonlinear constraint defconst2: (4.0 * const1 - 1.0) / (4.0 * const1 - 4.0) + 0.615 / const1 - const2 == 0.0 */
   {
      SCIP_CONSEXPR_EXPR* exprs[3];
      SCIP_CONSEXPR_EXPR* sumexpr1;
      SCIP_CONSEXPR_EXPR* sumexpr2;
      SCIP_CONSEXPR_EXPR* powexpr1;
      SCIP_CONSEXPR_EXPR* powexpr2;
      SCIP_CONSEXPR_EXPR* prodexpr;
      SCIP_CONSEXPR_EXPR* expr;
      SCIP_Real coefs[3];

      /* create (4.0 * const1 - 1.0) */
      coefs[0] = 4.0;
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, &sumexpr1, 1, &const1expr, coefs, -1.0) );

      /* create (4.0 * const1 - 4.0) */
      coefs[0] = 4.0;
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, &sumexpr2, 1, &const1expr, coefs, -4.0) );

      /* create 1 / (4.0 * const1 - 4.0) */
      SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, &powexpr1, sumexpr2, -1.0) );

      /* create (4.0 * const1 - 1.0) / (4.0 * const1 - 4.0) */
      exprs[0] = sumexpr1;
      exprs[1] = powexpr1;
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prodexpr, 2, exprs, 1.0) );

      /* create 1 / const1 */
      SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, &powexpr2, const1expr, -1.0) );

      /* create (4.0 * const1 - 1.0) / (4.0 * const1 - 4.0) + 0.615 / const1 - const2 */
      exprs[0] = prodexpr;
      coefs[0] = 1.0;
      exprs[1] = powexpr2;
      coefs[1] = 0.615;
      exprs[2] = const2expr;
      coefs[2] = -1.0;
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, &expr, 3, exprs, coefs, 0.0) );

      /* create nonlinear constraint */
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &defconst2, "defconst2", expr, 0.0, 0.0) );

      /* release expressions */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexpr2) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexpr1) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr2) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr1) );
   }

   /* quadratic constraint shear: 8.0*maxworkload/PI * const1 * const2 - maxshearstress * wire^2 <= 0.0 */
   {
      SCIP_VAR* quadvars1[2] = {const1, wire};
      SCIP_VAR* quadvars2[2] = {const2, wire};
      SCIP_Real quadcoefs[2] = {8.0 * maxworkload / M_PI, -maxshearstress};

      /* create empty quadratic constraint with right-hand-side 0.0 */
      SCIP_CALL( SCIPcreateConsExprQuadratic(scip, &shear, "shear", 0, NULL, NULL, 2, quadvars1, quadvars2, quadcoefs,
         -SCIPinfinity(scip), 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   }

   /* nonlinear constraint defdefl: 8.0/shearmod * ncoils * const1^3 / wire - defl == 0.0 */
   {
      SCIP_CONSEXPR_EXPR* exprs[3];
      SCIP_CONSEXPR_EXPR* prodexpr;
      SCIP_CONSEXPR_EXPR* powexpr1;
      SCIP_CONSEXPR_EXPR* powexpr2;
      SCIP_CONSEXPR_EXPR* expr;
      SCIP_Real coefs[3];

      /* create const1^3 */
      SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, &powexpr1, const1expr, 3.0) );

      /* create 1 / wire */
      SCIP_CALL( SCIPcreateConsExprExprPow(scip, consexprhdlr, &powexpr2, wireexpr, -1.0) );

      /* create ncoils * const1^3 / wire */
      exprs[0] = ncoilsexpr;
      exprs[1] = powexpr1;
      exprs[2] = powexpr2;
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prodexpr, 3, exprs, 1.0) );

      /* create 8.0/shearmod * ncoils * const1^3 / wire - defl */
      exprs[0] = prodexpr;
      coefs[0] = 8.0 / shearmod;
      exprs[1] = deflexpr;
      coefs[1] = -1.0;
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, &expr, 2, exprs, coefs, 0.0) );

      /* create nonlinear constraint */
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &defdefl, "defdefl", expr, 0.0, 0.0) );

      /* release expressions */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexpr2) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &powexpr1) );
   }

   /* quadratic constraint freel: maxworkload * defl + 1.05 * ncoils * wire + 2.1 * wire <= maxfreelen */
   {
      SCIP_VAR* linvars[2] = {defl, wire};
      SCIP_Real lincoefs[2] = {maxworkload, 2.1};
      SCIP_Real one05 = 1.05;

      /* create quadratic constraint maxworkload * defl + 1.05 * ncoils * wire <= maxfreelen */
      SCIP_CALL( SCIPcreateConsExprQuadratic(scip, &freel, "freel", 2, linvars, lincoefs, 1, &ncoils, &wire, &one05,
         -SCIPinfinity(scip), maxfreelen, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
         FALSE, FALSE) );
   }

   /* linear constraint coilwidth: coil + wire <= maxcoildiam */
   {
      /* create empty linear constraint with right-hand-side maxcoildiam */
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &coilwidth, "coilwidth", 0, NULL, NULL, -SCIPinfinity(scip), maxcoildiam) );

      /* add linear term coil + wire */
      SCIP_CALL( SCIPaddCoefLinear(scip, coilwidth, coil, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, coilwidth, wire, 1.0) );
   }

   /* linear constraint defwire: sum_i b[i]*y[i] - wire == 0.0 */
   {
      /* create linear constraint sum_i b[i]*y[i] == 0.0 */
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &defwire, "defwire", nwires, y, (SCIP_Real*)diameters, 0.0, 0.0) );

      /* add term -wire */
      SCIP_CALL( SCIPaddCoefLinear(scip, defwire, wire, -1.0) );
   }

   /* specialized linear constraint selectwire: sum_i y[i] == 1.0 */
   {
      SCIP_CALL( SCIPcreateConsBasicSetpart(scip, &selectwire, "selectwire", nwires, y) );
   }

   /* add constraints to problem */
   SCIP_CALL( SCIPaddCons(scip, voldef) );
   SCIP_CALL( SCIPaddCons(scip, defconst1) );
   SCIP_CALL( SCIPaddCons(scip, defconst2) );
   SCIP_CALL( SCIPaddCons(scip, shear) );
   SCIP_CALL( SCIPaddCons(scip, defdefl) );
   SCIP_CALL( SCIPaddCons(scip, freel) );
   SCIP_CALL( SCIPaddCons(scip, coilwidth) );
   SCIP_CALL( SCIPaddCons(scip, defwire) );
   SCIP_CALL( SCIPaddCons(scip, selectwire) );

   /* release variable expressions */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &volumeexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &const2expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &const1expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &ncoilsexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &deflexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &wireexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &coilexpr) );

   /* release variables and constraints
    * the problem has them captured, and we do not require them anymore
    */
   SCIP_CALL( SCIPreleaseVar(scip, &coil) );
   SCIP_CALL( SCIPreleaseVar(scip, &wire) );
   SCIP_CALL( SCIPreleaseVar(scip, &defl) );
   SCIP_CALL( SCIPreleaseVar(scip, &ncoils) );
   SCIP_CALL( SCIPreleaseVar(scip, &const1) );
   SCIP_CALL( SCIPreleaseVar(scip, &const2) );
   SCIP_CALL( SCIPreleaseVar(scip, &volume) );
   for( i = 0; i < nwires; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &y[i]) );
   }

   SCIP_CALL( SCIPreleaseCons(scip, &voldef) );
   SCIP_CALL( SCIPreleaseCons(scip, &defconst1) );
   SCIP_CALL( SCIPreleaseCons(scip, &defconst2) );
   SCIP_CALL( SCIPreleaseCons(scip, &shear) );
   SCIP_CALL( SCIPreleaseCons(scip, &defdefl) );
   SCIP_CALL( SCIPreleaseCons(scip, &freel) );
   SCIP_CALL( SCIPreleaseCons(scip, &coilwidth) );
   SCIP_CALL( SCIPreleaseCons(scip, &defwire) );
   SCIP_CALL( SCIPreleaseCons(scip, &selectwire) );

   return SCIP_OKAY;
}

/** runs spring example */
static
SCIP_RETCODE runSpring(void)
{
   SCIP* scip;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "************************************************\n");
   SCIPinfoMessage(scip, NULL, "* Running Coil Compression Spring Design Model *\n");
   SCIPinfoMessage(scip, NULL, "************************************************\n");
   SCIPinfoMessage(scip, NULL, "\n");

   SCIP_CALL( setupProblem(scip) );

   SCIPinfoMessage(scip, NULL, "Original problem:\n");
   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );

   SCIPinfoMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPpresolve(scip) );

   /* SCIPinfoMessage(scip, NULL, "Reformulated problem:\n");
   SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
   */

   SCIPinfoMessage(scip, NULL, "\nSolving...\n");
   SCIP_CALL( SCIPsolve(scip) );

   if( SCIPgetNSols(scip) > 0 )
   {
      SCIPinfoMessage(scip, NULL, "\nSolution:\n");
      SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );
   }

   SCIP_CALL( SCIPfree(&scip) );

   return SCIP_OKAY;
}


/** main method starting SCIP */
int main(
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< array of shell arguments */
   )
{  /*lint --e{715}*/
   SCIP_RETCODE retcode;

   retcode = runSpring();

   /* evaluate return code of the SCIP process */
   if( retcode != SCIP_OKAY )
   {
      /* write error back trace */
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}

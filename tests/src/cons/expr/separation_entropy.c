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

/**@file   separation_entropy.c
 * @brief  tests separation of entropy()
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_entropy.c"
#include "separation.h"

Test(separation, entropy, .init = setup, .fini = teardown,
   .description = "test separation for an exponential expression"
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_ROW* cut;
   int i;

   SCIP_CALL( SCIPcreateConsExprExprEntropy(scip, conshdlr, &expr, zexpr) );

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   /* compute a cut for which we need an overestimation (linearization) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 0.0) );

   cut = NULL;
   SCIP_CALL( separatePointEntropy(scip,  conshdlr, expr, sol, TRUE, &cut) );

   cr_assert(cut != NULL);
   cr_assert_eq(SCIProwGetNNonz(cut), 2);
   cr_assert(SCIPisEQ(scip, SCIProwGetLhs(cut), -2.0));
   cr_assert_eq(SCIProwGetRhs(cut), SCIPinfinity(scip));

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];

      if( var == SCIPvarGetTransVar(z) )
         cr_assert(SCIPisEQ(scip, coef, -log(2.0) - 1.0));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_assert_eq(coef, -1.0);
      else
         cr_assert(FALSE, "found an unknown variable");
   }

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   /* compute a cut for which we need an underestimation (linearization) */
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, -10.0) );

   cut = NULL;
   SCIP_CALL( separatePointEntropy(scip,  conshdlr, expr, sol, FALSE, &cut) );

   cr_assert(cut != NULL);
   cr_assert_eq(SCIProwGetNNonz(cut), 2);
   cr_assert_eq(SCIProwGetLhs(cut), -SCIPinfinity(scip));
   cr_assert(SCIPisEQ(scip, SCIProwGetRhs(cut), -1.5 * log(3.0) + 1.5 * log(1.0)));

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];

      if( var == SCIPvarGetTransVar(z) )
         cr_assert_eq(coef, 0.5 * (-3.0 * log(3.0) + log(1.0)));
      else if( var == SCIPvarGetTransVar(auxvar) )
         cr_assert_eq(coef, -1.0);
      else
         cr_assert(FALSE, "found an unknown variable");
   }

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

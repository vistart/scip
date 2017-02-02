/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    nlpi_worhp.c
 * @ingroup NLPIS
 * @brief   Worhp NLP interface
 * @author  Benjamin Mueller
 * @author  Renke Kuhlmann
 *
 * @todo So far, Worhp can not handle the case that variables have been fixed before warm-starting. Remove the code in
 * nlpiChgVarBoundsWorhp when this has changed.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "nlpi/nlpi_worhp.h"
#include "nlpi/nlpi.h"
#include "nlpi/nlpioracle.h"
#include "nlpi/exprinterpret.h"
#include "scip/interrupt.h"
#include "scip/pub_misc.h"

#include <stdio.h>
#include <stdlib.h>

#include "worhp/worhp.h"

#define NLPI_NAME              "worhp"                      /**< short concise name of solver */
#define NLPI_DESC              "Worhp interface"            /**< description of solver */
#define NLPI_PRIORITY          1                            /**< priority of NLP solver */

#define DEFAULT_ALGORITHM      1                            /**< default algorithm to solve NLP (1: SQP 2: interior point) */
#define DEFAULT_VERBLEVEL      0                            /**< default verbosity level (0: normal 1: full 2: debug >2: more debug) */
#define DEFAULT_SCALEDKKT      TRUE                         /**< default whether KKT conditions are allowed to be scaled in the solver */

/*
 * Data structures
 */

struct SCIP_NlpiData
{
   BMS_BLKMEM*                 blkmem;       /**< block memory */
   SCIP_MESSAGEHDLR*           messagehdlr;  /**< message handler */
   SCIP_Real                   infinity;     /**< initial value for infinity */
};

struct SCIP_NlpiProblem
{
   SCIP_NLPIORACLE*            oracle;       /**< Oracle-helper to store and evaluate NLP */
   BMS_BLKMEM*                 blkmem;       /**< block memory */

   SCIP_NLPTERMSTAT            lasttermstat; /**< termination status from last run */
   SCIP_NLPSOLSTAT             lastsolstat;  /**< solution status from last run */
   SCIP_Real                   lasttime;     /**< time spend in last run */
   int                         lastniter;    /**< number of iterations in last run */

   SCIP_Real*                  lastprimal;   /**< primal solution from last run, if available */
   SCIP_Real*                  lastdualcons; /**< dual solution from last run, if available */
   SCIP_Real*                  lastduallb;   /**< dual solution for lower bounds from last run, if available */
   SCIP_Real*                  lastdualub;   /**< dual solution for upper bounds from last run, if available */
   int                         lastprimalsize; /**< size of lastprimal array */
   int                         lastdualconssize; /**< size of lastdualcons array */
   int                         lastduallbsize; /**< size of lastduallb array */
   int                         lastdualubsize; /**< size of lastdualub array */

   SCIP_Bool                   firstrun;     /**< whether the next NLP solve will be the first one (with the current problem structure) */
   SCIP_Real*                  initguess;    /**< initial values for primal variables, or NULL if not known */

   /* Worhp data structures */
   OptVar*                     opt;          /**< Worhp variables  */
   Workspace*                  wsp;          /**< Worhp working space */
   Params*                     par;          /**< Worhp parameters */
   Control*                    cnt;          /**< Worhp control */

   /* parameters */
   SCIP_Real                   feastol;      /**< feasibility tolerance for primal variables and slacks */
   SCIP_Real                   relobjtol;    /**< relative objective tolerance */
   SCIP_Real                   lobjlim;      /**< lower objective limit (cutoff) */
   SCIP_Real                   timelim;      /**< NLP time limit */
   int                         fromscratch;  /**< solver should start from scratch at next call?: 0 no, 1 yes */
   int                         verblevel;    /**< verbosity level of output of NLP solver to the screen: 0 off, 1 normal, 2 debug, > 2 more debug */
   int                         itlim;        /**< NLP iteration limit */
   int                         fastfail;     /**< should the NLP solver stop early if convergence is slow?: 0 no, 1 yes */
};

/*
 * Local methods
 */

/** clears the last solution information */
static
void invalidateSolution(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(problem != NULL);
   assert(problem->blkmem != NULL);

   BMSfreeBlockMemoryArrayNull(problem->blkmem, &(problem->lastprimal), problem->lastprimalsize);
   BMSfreeBlockMemoryArrayNull(problem->blkmem, &(problem->lastdualcons), problem->lastdualconssize);
   BMSfreeBlockMemoryArrayNull(problem->blkmem, &(problem->lastduallb), problem->lastduallbsize);
   BMSfreeBlockMemoryArrayNull(problem->blkmem, &(problem->lastdualub), problem->lastdualubsize);

   problem->lastprimal = NULL;
   problem->lastdualcons = NULL;
   problem->lastduallb = NULL;
   problem->lastdualub = NULL;
   problem->lastprimalsize = 0;
   problem->lastdualconssize = 0;
   problem->lastduallbsize = 0;
   problem->lastdualubsize = 0;
   problem->lastsolstat = SCIP_NLPSOLSTAT_UNKNOWN;
   problem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
}

/** evaluate last Worhp run */
static
SCIP_RETCODE evaluateWorhpRun(
   SCIP_NLPIPROBLEM*     problem,            /**< pointer to problem data structure */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler (might be NULL) */
   )
{
   int i;

   assert(problem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->par != NULL);
   assert(problem->cnt != NULL);

  /* initialization error */
  if( problem->cnt->status == InitError )
  {
     SCIPdebugMessage("Worhp failed because of initialization error!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_MEMERR;
  }
  /* data error */
  else if( problem->cnt->status == DataError )
  {
     SCIPdebugMessage("Worhp failed because of data error!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
  }
  /* license error */
  else if( problem->cnt->status == LicenseError )
  {
     SCIPerrorMessage("Worhp failed because of license error!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_LICERR;
  }

  /* evaluation errors */
  else if( problem->cnt->status == evalsNaN )
  {
     SCIPdebugMessage("Worhp failed because of a NaN value in an evaluation!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_EVALERR;
  }

  /*
   * catch unsuccessful states
   */

  /* numerical errors during solution of NLP */
  else if( (problem->cnt->status == QPerror) || (problem->cnt->status == MinimumStepsize) || (problem->cnt->status == TooBig) )
  {
     SCIPdebugMessage("Worhp failed because of a numerical error during optimization!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_NUMERR;
  }
  /* maximal number of calls or iteration */
  else if( (problem->cnt->status == MaxCalls) || (problem->cnt->status == MaxIter) )
  {
     SCIPdebugMessage("Worhp failed because maximal number of calls or iterations is reached!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_ITLIM;
  }
  /* time limit reached */
  else if( problem->cnt->status == Timeout )
  {
     SCIPdebugMessage("Worhp failed because time limit is reached!\n");
     invalidateSolution(problem);
     problem->lastsolstat = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_TILIM;
  }
  /* infeasible stationary point found */
  else if( problem->cnt->status == LocalInfeas || problem->cnt->status == LocalInfeasOptimal )
  {
     SCIPdebugMessage("Worhp failed because of convergence against infeasible stationary point!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_LOCINFEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* regularization of Hessian matrix failed */
  else if( problem->cnt->status == RegularizationFailed )
  {
     SCIPdebugMessage("Worhp failed because of regularization of Hessian matrix failed!\n");
     invalidateSolution(problem);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_NUMERR;
  }

  /*
   * catch successful states
   */

  /* everything went fine */
  else if( problem->cnt->status == OptimalSolution || problem->cnt->status == OptimalSolutionBoxEqual )
  {
     SCIPdebugMessage("Worhp terminated successfully at a local optimum!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_LOCOPT;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* feasible point, KKT conditions are not satisfied, and Worhp thinks that there is no objective function */
  else if( problem->cnt->status == OptimalSolutionConstantF )
  {
     SCIPdebugMessage("Worhp terminated successfully with a feasible point but KKT are not met!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* feasible point but KKT conditions are violated in unscaled space */
  else if( problem->cnt->status == AcceptableSolutionSKKT || problem->cnt->status == AcceptableSolutionScaled
     || problem->cnt->status == AcceptablePreviousScaled )
  {
     SCIPdebugMessage("Worhp terminated successfully with a feasible point but KKT are violated in unscaled space!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* feasible and no further progress */
  else if( problem->cnt->status == LowPassFilterOptimal )
  {
     SCIPdebugMessage("Worhp terminated at feasible solution without further progress!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* feasible and in feasibility mode, i.e., optimality not required */
  else if( problem->cnt->status == FeasibleSolution )
  {
     SCIPdebugMessage("Worhp terminated at feasible solution, optimality was not required!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* acceptable solution found, but stopped due to limit or error */
  else if( problem->cnt->status == AcceptableSolution || problem->cnt->status == AcceptableSolutionConstantF )
  {
     SCIPdebugMessage("Worhp terminated at acceptable solution due to limit or error!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* previously acceptable solution was found, but stopped due to limit or error */
  else if( problem->cnt->status == AcceptablePrevious || problem->cnt->status == AcceptablePreviousConstantF )
  {
     SCIPdebugMessage("Worhp previously found acceptable solution but terminated due to limit or error!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* acceptable solution found, and no further progress */
  else if( problem->cnt->status == LowPassFilterAcceptable )
  {
     SCIPdebugMessage("Worhp found acceptable solution but terminated due no further progress!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* acceptable solution found, but search direction is small or zero */
  else if( (problem->cnt->status == SearchDirectionZero) || (problem->cnt->status == SearchDirectionSmall) )
  {
     SCIPdebugMessage("Worhp found acceptable solution but search direction is small or zero!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }
  /* acceptable solution found, but not optimal */
  else if( (problem->cnt->status == FritzJohn) || (problem->cnt->status == NotDiffable) || (problem->cnt->status == Unbounded) )
  {
     SCIPdebugMessage("Worhp found acceptable solution but terminated perhaps due to nondifferentiability, unboundedness or at Fritz John point!\n");
     problem->lastsolstat  = SCIP_NLPSOLSTAT_FEASIBLE;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OKAY;
  }

  /*
   * catch other states
   */
  else
  {
     SCIPerrorMessage("Worhp returned with unknown solution status %d\n", problem->cnt->status);
     problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
     problem->lasttermstat = SCIP_NLPTERMSTAT_OTHER;
  }

  /* store solution */
  if( problem->lastprimal == NULL )
  {
     SCIP_ALLOC( BMSduplicateBlockMemoryArray(problem->blkmem, &problem->lastprimal, problem->opt->X, problem->opt->n) );
     problem->lastprimalsize = problem->opt->n;

     SCIP_ALLOC( BMSduplicateBlockMemoryArray(problem->blkmem, &problem->lastdualcons, problem->opt->Mu, problem->opt->m) );
     problem->lastdualconssize = problem->opt->m;

     SCIP_ALLOC( BMSallocBlockMemoryArray(problem->blkmem, &problem->lastduallb, problem->opt->n) );
     problem->lastduallbsize = problem->opt->n;

     SCIP_ALLOC( BMSallocBlockMemoryArray(problem->blkmem, &problem->lastdualub, problem->opt->n) );
     problem->lastdualubsize = problem->opt->n;
  }
  else
  {
     BMScopyMemoryArray(problem->lastprimal, problem->opt->X, problem->opt->n);
     BMScopyMemoryArray(problem->lastdualcons, problem->opt->Mu, problem->opt->m);
  }

  for( i = 0; i < problem->opt->n; ++i )
  {
     if( problem->opt->Lambda[i] <= 0.0 )
     {
        problem->lastduallb[i] = -problem->opt->Lambda[i];
        problem->lastdualub[i] = 0.0;
     }
     else
     {
        problem->lastduallb[i] = 0.0;
        problem->lastdualub[i] = problem->opt->Lambda[i];
     }
  }

  assert(problem->lastprimal != NULL);
  assert(problem->lastdualcons != NULL);
  assert(problem->lastduallb != NULL);
  assert(problem->lastdualub != NULL);

  return SCIP_OKAY;
}

/** evaluates objective function and store the result in the corresponding Worhp data fields */
static
SCIP_RETCODE userF(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   SCIP_Real objval;

   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->blkmem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n == SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m == SCIPnlpiOracleGetNConstraints(problem->oracle));

   SCIP_CALL( SCIPnlpiOracleEvalObjectiveValue(problem->oracle, problem->opt->X, &objval) );
   problem->opt->F = problem->wsp->ScaleObj * objval;

#ifdef SCIP_DEBUG_USERF
   {
      int i;

      printf("userF()\n");
      for( i = 0; i < problem->opt->n; ++i )
         printf("  x[%d] = %g\n", i, problem->opt->X[i]);
      printf("  obj = %g\n", problem->opt->F);
   }
#endif

   return SCIP_OKAY;
}

/** evaluates constraints and store the result in the corresponding Worhp data fields */
static
SCIP_RETCODE userG(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->blkmem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n == SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m == SCIPnlpiOracleGetNConstraints(problem->oracle));

   SCIP_CALL( SCIPnlpiOracleEvalConstraintValues(problem->oracle, problem->opt->X, problem->opt->G) );

#ifdef SCIP_DEBUG_USERG
   {
      int i;

      printf("userG()\n");
      for( i = 0; i < problem->opt->n; ++i )
         printf("  x[%d] = %g\n", i, problem->opt->X[i]);

      for( i = 0; i < problem->opt->m; ++i )
         printf("  cons[%d] = %g\n", i, problem->opt->G[i]);
   }
#endif

   return SCIP_OKAY;
}

/** computes objective gradient and store the result in the corresponding Worhp data fields */
static
SCIP_RETCODE userDF(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   SCIP_Real objval;

   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->blkmem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n == SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m == SCIPnlpiOracleGetNConstraints(problem->oracle));

   /* TODO this needs to be changed if we store the gradient of the objective function in a sparse format */
   SCIP_CALL( SCIPnlpiOracleEvalObjectiveGradient(problem->oracle, problem->opt->X, TRUE, &objval,
         problem->wsp->DF.val) );

   /* scale gradient if necessary */
   if( problem->wsp->ScaleObj != 1.0 )
   {
      int i;
      for( i = 0; i < problem->opt->n; ++i )
         problem->wsp->DF.val[i] *= problem->wsp->ScaleObj;
   }

#ifdef SCIP_DEBUG_USERDF
   {
      int i;

      printf("userDF()\n");
      for( i = 0; i < problem->opt->n; ++i )
         printf("  x[%d] = %g\n", i, problem->opt->X[i]);

      for( i = 0; i < problem->opt->n; ++i )
         printf("  DF[%d] = %g\n", i, problem->wsp->DF.val[i]);
   }
#endif

   return SCIP_OKAY;
}

/** computes jacobian matrix and store the result in the corresponding Worhp data fields */
static
SCIP_RETCODE userDG(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   SCIP_RETCODE retcode;
   SCIP_Real* jacvals;
   int i;

   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->blkmem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n == SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m == SCIPnlpiOracleGetNConstraints(problem->oracle));

   SCIP_ALLOC( BMSallocBlockMemoryArray(problem->blkmem, &jacvals, problem->wsp->DG.nnz) );
   retcode = SCIPnlpiOracleEvalJacobian(problem->oracle, problem->opt->X, TRUE, NULL, jacvals);

   if( retcode == SCIP_OKAY )
   {
      /* map values with DG indices */
      for( i = 0; i < problem->wsp->DG.nnz; ++i )
      {
         problem->wsp->DG.val[i] = jacvals[ problem->wsp->DG.perm[i]-1 ];
      }

#ifdef SCIP_DEBUG_USERDG
      printf("userDG()\n");
      for( i = 0; i < problem->opt->n; ++i )
         printf("  x[%d] = %g\n", i, problem->opt->X[i]);
      for( i = 0; i < problem->wsp->DG.nnz; ++i )
         printf("  DG[%d] = %g\n", i, problem->wsp->DG.val[i]);
#endif
   }

   /* free memory */
   BMSfreeBlockMemoryArray(problem->blkmem, &jacvals, problem->wsp->DG.nnz);

   return retcode;
}

/** computes hessian matrix and store the result in the corresponding Worhp data fields */
static
SCIP_RETCODE userHM(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   const int* offset;
   SCIP_Real* hessianvals;
   SCIP_RETCODE retcode;
   int nnonz;
   int i;

   assert(problem != NULL);
   assert(problem->oracle != NULL);
   assert(problem->blkmem != NULL);
   assert(problem->opt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->opt->n == SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m == SCIPnlpiOracleGetNConstraints(problem->oracle));

   /* get nonzero entries in HM of SCIP (excludes unused diagonal entries) */
   SCIP_CALL( SCIPnlpiOracleGetHessianLagSparsity(problem->oracle, &offset, NULL) );
   nnonz = offset[problem->opt->n];

   /* evaluate hessian */
   SCIP_ALLOC( BMSallocBlockMemoryArray(problem->blkmem, &hessianvals, problem->wsp->HM.nnz) );
   retcode = SCIPnlpiOracleEvalHessianLag(problem->oracle, problem->opt->X, TRUE, problem->wsp->ScaleObj,
         problem->opt->Mu, hessianvals);

   if( retcode == SCIP_OKAY )
   {
      assert(problem->wsp->HM.nnz >= nnonz);
      for( i = 0; i < problem->wsp->HM.nnz; ++i )
      {
         /* an entry i with HM.perm[i] - 1 >= nnonz corresponds to an in SCIP non-existing diagonal element */
         if( problem->wsp->HM.perm[i] - 1 >= nnonz )
            problem->wsp->HM.val[i] = 0.0;
         else
            problem->wsp->HM.val[i] = hessianvals[ problem->wsp->HM.perm[i] - 1 ];
      }

#ifdef SCIP_DEBUG_HM
      printf("userHM()\n");
      for( i = 0; i < problem->opt->n; ++i )
         printf("  x[%d] = %g\n", i, problem->opt->X[i]);
      for( i = 0; i < problem->wsp->HM.nnz; ++i )
         printf("  HM[%d] = %g\n", i, problem->wsp->HM.val[i]);
#endif
   }

   /* free memory */
   BMSfreeBlockMemoryArray(problem->blkmem, &hessianvals, problem->wsp->HM.nnz);

   return retcode;
}

/** initialize Worhp data */
static
SCIP_RETCODE initWorhp(
   SCIP_NLPI*            nlpi,               /**< pointer to NLPI datastructure */
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   Workspace* wsp = problem->wsp;
   Control* cnt = problem->cnt;
   OptVar* opt = problem->opt;
   Params* par = problem->par;

   const SCIP_Real* lbs;
   const SCIP_Real* ubs;
   const int* offset;
   const int* cols;
   SCIP_NLPIDATA* nlpidata;
   int i;
   int j;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->cnt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->par != NULL);
   assert(problem->opt != NULL);
   assert(problem->firstrun);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   /* properly zeros everything */
   WorhpPreInit(opt, wsp, par, cnt);

   /* set problem dimensions */
   opt->n = SCIPnlpiOracleGetNVars(problem->oracle);
   opt->m = SCIPnlpiOracleGetNConstraints(problem->oracle);
   SCIPdebugMessage("nvars %d nconss %d\n", opt->n, opt->m);

   /* assume that objective function is dense; TODO use sparse representation */
   wsp->DF.nnz = opt->n;

   /* get number of non-zero entries in Jacobian */
   SCIP_CALL( SCIPnlpiOracleGetJacobianSparsity(problem->oracle, &offset, NULL) );
   wsp->DG.nnz = offset[opt->m];
   SCIPdebugMessage("nnonz jacobian %d\n", wsp->DG.nnz);

   /* get number of non-zero entries in hessian
    *
    * note that Worhp wants to have the full diagonal in ANY case
    */
   SCIP_CALL( SCIPnlpiOracleGetHessianLagSparsity(problem->oracle, &offset, &cols) );
   wsp->HM.nnz = 0;

   j = offset[0];
   for( i = 0; i < opt->n; ++i )
   {
      /* diagonal element */
      ++(wsp->HM.nnz);

      /* strict lower triangle elements */
      for( ; j < offset[i+1]; ++j )
      {
         if( i != cols[j] )
         {
            assert(i > cols[j]);
            ++(wsp->HM.nnz);
         }
      }
   }
   assert(offset[opt->n] <= wsp->HM.nnz);
   SCIPdebugMessage("nnonz hessian %d\n", wsp->HM.nnz);

   /* initialize data in Worhp */
   WorhpInit(opt, wsp, par, cnt);
   if (cnt->status != FirstCall)
   {
      SCIPerrorMessage("Initialisation failed.\n");
      return SCIP_ERROR;
   }

   /* set variable bounds */
   lbs = SCIPnlpiOracleGetVarLbs(problem->oracle);
   ubs = SCIPnlpiOracleGetVarUbs(problem->oracle);
   for( i = 0; i < opt->n; ++i )
   {
      opt->XL[i] = lbs[i];
      opt->XU[i] = ubs[i];
      SCIPdebugMessage("bounds %d [%g,%g]\n", i, opt->XL[i], opt->XU[i]);
   }

   /* set constraint sides */
   for( i = 0; i < opt->m; ++i )
   {
      opt->GL[i] = SCIPnlpiOracleGetConstraintLhs(problem->oracle, i);
      opt->GU[i] = SCIPnlpiOracleGetConstraintRhs(problem->oracle, i);
      SCIPdebugMessage("sides %d [%g,%g]\n", i, opt->GL[i], opt->GU[i]);
   }

   /* set column indices of objective function; note that indices go from 1 to n */
   /* if( wsp->DF.NeedStructure ) evaluates to FALSE if DF is dense */
   {
      SCIPdebugMessage("column indices of objective function:");
      for( i = 0; i < opt->n; ++i )
      {
         wsp->DF.row[i] = i + 1;
         SCIPdebugMessage(" %d", wsp->DF.row[i]);
      }
      SCIPdebugMessage("\n");
   }

   /* set column and row indices of non-zero entries in Jacobian matrix */
   /* if( wsp->DG.NeedStructure ) evaluates to FALSE if DG is dense */
   {
      int nnonz;

      SCIP_CALL( SCIPnlpiOracleGetJacobianSparsity(problem->oracle, &offset, &cols) );
      assert(offset[opt->m] == wsp->DG.nnz);

      nnonz = 0;
      j = offset[0];
      for( i = 0; i < opt->m; ++i )
      {
         for( ; j < offset[i+1]; ++j )
         {
            wsp->DG.row[nnonz] = i + 1;
            wsp->DG.col[nnonz] = cols[j] + 1;
            ++nnonz;
         }
      }
      assert(nnonz == wsp->DG.nnz);

      /* sort arrays w.r.t the column-major order */
      SortWorhpMatrix(&wsp->DG);
   }

   /* set column and row indices of non-zero entries in hessian matrix */
   if( problem->par->UserHM || problem->par->FidifHM || problem->par->BFGSmethod > 1 )
   {
      int nnonz;
      int k;

      SCIP_CALL( SCIPnlpiOracleGetHessianLagSparsity(problem->oracle, &offset, &cols) );
      assert(offset[opt->n] <= wsp->HM.nnz);

      k = offset[opt->n];
      nnonz = 0;
      j = offset[0];
      for( i = 0; i < opt->n; ++i )
      {
         SCIP_Bool adddiag = TRUE;

         for( ; j < offset[i+1]; ++j )
         {
            problem->wsp->HM.row[nnonz] = i + 1;
            problem->wsp->HM.col[nnonz] = cols[j] + 1;
            ++nnonz;

            if( i == cols[j] )
               adddiag = FALSE;
         }

         /* Worhp wants to have each diagonal element */
         if( adddiag )
         {
            problem->wsp->HM.row[k] = i + 1;
            problem->wsp->HM.col[k] = i + 1;
            ++k;
         }
      }
      assert(nnonz == offset[opt->n]);
      assert(k == wsp->HM.nnz);

      /* sort arrays w.r.t the LT column-major order */
      SortWorhpMatrix(&wsp->HM);

#ifdef SCIP_DEBUG
      SCIPdebugMessage("column and row indices of hessian:\n");
      for( i = 0; i < wsp->HM.nnz; ++i )
      {
         SCIPdebugMessage("entry %d: (row,col) = (%d,%d)\n", i, wsp->HM.row[i], wsp->HM.col[i]);
      }
#endif
   }

   return SCIP_OKAY;
}

/** update Worhp data */
static
SCIP_RETCODE updateWorhp(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   const SCIP_Real* lbs;
   const SCIP_Real* ubs;
   int i;

   assert(problem != NULL);
   assert(problem->cnt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->par != NULL);
   assert(problem->opt != NULL);
   assert(problem->oracle != NULL);
   assert(problem->opt->n == SCIPnlpiOracleGetNVars(problem->oracle));
   assert(problem->opt->m == SCIPnlpiOracleGetNConstraints(problem->oracle));

   WorhpRestart(problem->opt, problem->wsp, problem->par, problem->cnt);

   /* update variable bounds */
   lbs = SCIPnlpiOracleGetVarLbs(problem->oracle);
   ubs = SCIPnlpiOracleGetVarUbs(problem->oracle);
   for( i = 0; i < problem->opt->n; ++i )
   {
      problem->opt->XL[i] = lbs[i];
      problem->opt->XU[i] = ubs[i];
   }

   /* update constraint sides */
   for( i = 0; i < problem->opt->m; ++i )
   {
      problem->opt->GL[i] = SCIPnlpiOracleGetConstraintLhs(problem->oracle, i);
      problem->opt->GU[i] = SCIPnlpiOracleGetConstraintRhs(problem->oracle, i);
   }

   return SCIP_OKAY;
}

/** frees Worhp data */
static
SCIP_RETCODE freeWorhp(
   SCIP_NLPIPROBLEM*     problem             /**< pointer to problem data structure */
   )
{
   assert(problem != NULL);
   assert(problem->cnt != NULL);
   assert(problem->wsp != NULL);
   assert(problem->par != NULL);
   assert(problem->opt != NULL);

   if( problem->opt->initialised )
      WorhpFree(problem->opt, problem->wsp, problem->par, problem->cnt);

   return SCIP_OKAY;
}

/*
 * Callback methods of NLP solver interface
 */

/* TODO: Implement all necessary NLP interface methods. The methods with an #if 0 ... #else #define ... are optional
 * (currently, all methods are required) */

/** copy method of NLP interface (called when SCIP copies plugins)
 *
 * input:
 *  - blkmem block memory in target SCIP
 *  - sourcenlpi the NLP interface to copy
 *  - targetnlpi buffer to store pointer to copy of NLP interface
 */
static
SCIP_DECL_NLPICOPY( nlpiCopyWorhp )
{
   SCIP_NLPIDATA* sourcedata;
   SCIP_NLPIDATA* targetdata;

   assert(sourcenlpi != NULL);
   assert(targetnlpi != NULL);

   sourcedata = SCIPnlpiGetData(sourcenlpi);
   assert(sourcedata != NULL);

   SCIP_CALL( SCIPcreateNlpSolverWorhp(blkmem, targetnlpi) );
   assert(*targetnlpi != NULL);

   SCIP_CALL( SCIPnlpiSetRealPar(*targetnlpi, NULL, SCIP_NLPPAR_INFINITY, sourcedata->infinity) );
   SCIP_CALL( SCIPnlpiSetMessageHdlr(*targetnlpi, sourcedata->messagehdlr) );

   targetdata = SCIPnlpiGetData(*targetnlpi);
   assert(targetdata != NULL);

   targetdata->blkmem = sourcedata->blkmem;
   targetdata->messagehdlr = sourcedata->messagehdlr;

   /* copy parameter */
   targetdata->infinity = sourcedata->infinity;

   return SCIP_OKAY;
}  /*lint !e715*/

/** destructor of NLP interface to free nlpi data
 *
 * input:
 *  - nlpi datastructure for solver interface
 */
static
SCIP_DECL_NLPIFREE( nlpiFreeWorhp )
{
   SCIP_NLPIDATA* data;

   assert(nlpi != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);
   assert(data->blkmem != NULL);

   BMSfreeMemory(&data);

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets pointer for NLP solver
 *
 *  to do dirty stuff
 *
 * input:
 *  - nlpi datastructure for solver interface
 *
 * return: void pointer to solver
 */
static
SCIP_DECL_NLPIGETSOLVERPOINTER(nlpiGetSolverPointerWorhp)
{
   assert(nlpi != NULL);

   return NULL;
}  /*lint !e715*/

/** creates a problem instance
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem pointer to store the problem data
 *  - name name of problem, can be NULL
 */
static
SCIP_DECL_NLPICREATEPROBLEM(nlpiCreateProblemWorhp)
{
   SCIP_NLPIDATA* data;

   assert(nlpi    != NULL);
   assert(problem != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(data->blkmem, problem) );
   if( *problem == NULL )
      return SCIP_NOMEMORY;

   /* initialize problem */
   BMSclearMemory((*problem));
   (*problem)->blkmem = data->blkmem;
   (*problem)->firstrun = TRUE;
   SCIP_CALL( SCIPnlpiOracleCreate(data->blkmem, &(*problem)->oracle) );
   SCIP_CALL( SCIPnlpiOracleSetInfinity((*problem)->oracle, data->infinity) );
   SCIP_CALL( SCIPnlpiOracleSetProblemName((*problem)->oracle, name) );

   /* allocate memory for Worhp data */
   SCIP_ALLOC( BMSallocBlockMemory(data->blkmem, &(*problem)->opt) );
   SCIP_ALLOC( BMSallocBlockMemory(data->blkmem, &(*problem)->wsp) );
   SCIP_ALLOC( BMSallocBlockMemory(data->blkmem, &(*problem)->par) );
   SCIP_ALLOC( BMSallocBlockMemory(data->blkmem, &(*problem)->cnt) );
   WorhpPreInit((*problem)->opt, (*problem)->wsp, (*problem)->par, (*problem)->cnt);

   /* set default parameters */
   (*problem)->feastol = 1e-9;
   (*problem)->relobjtol = 1e-9;
   (*problem)->lobjlim = SCIP_INVALID;
   (*problem)->timelim = SCIP_DEFAULT_INFINITY;
   (*problem)->fromscratch = 0;
   (*problem)->verblevel = DEFAULT_VERBLEVEL;
   (*problem)->itlim = INT_MAX;
   (*problem)->fastfail = 0;

   return SCIP_OKAY;
}  /*lint !e715*/

/** free a problem instance
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem pointer where problem data is stored
 */
static
SCIP_DECL_NLPIFREEPROBLEM(nlpiFreeProblemWorhp)
{
   SCIP_NLPIDATA* data;

   assert(nlpi     != NULL);
   assert(problem  != NULL);
   assert(*problem != NULL);

   data = SCIPnlpiGetData(nlpi);
   assert(data != NULL);

   if( (*problem)->opt != NULL )
   {
      /* free memory for last solution information */
      invalidateSolution(*problem);

      assert((*problem)->wsp != NULL);
      assert((*problem)->par != NULL);
      assert((*problem)->cnt != NULL);

      /* free Worhp data */
      SCIP_CALL( freeWorhp(*problem) );
      BMSfreeBlockMemory(data->blkmem, &(*problem)->cnt);
      BMSfreeBlockMemory(data->blkmem, &(*problem)->par);
      BMSfreeBlockMemory(data->blkmem, &(*problem)->wsp);
      BMSfreeBlockMemory(data->blkmem, &(*problem)->opt);
   }

   if( (*problem)->oracle != NULL )
   {
      SCIP_CALL( SCIPnlpiOracleFree(&(*problem)->oracle) );
   }

   BMSfreeMemoryArrayNull(&(*problem)->initguess);
   BMSfreeBlockMemory(data->blkmem, problem);
   *problem = NULL;

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets pointer to solver-internal problem instance
 *
 *  to do dirty stuff
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 * return: void pointer to problem instance
 */
static
SCIP_DECL_NLPIGETPROBLEMPOINTER(nlpiGetProblemPointerWorhp)
{
   assert(nlpi    != NULL);
   assert(problem != NULL);

   return NULL;
}  /*lint !e715*/

/** add variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nvars number of variables
 *  - lbs lower bounds of variables, can be NULL if -infinity
 *  - ubs upper bounds of variables, can be NULL if +infinity
 *  - varnames names of variables, can be NULL
 */
static
SCIP_DECL_NLPIADDVARS( nlpiAddVarsWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleAddVars(problem->oracle, nvars, lbs, ubs, varnames) );

   BMSfreeMemoryArrayNull(&problem->initguess);
   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/


/** add constraints
 * quadratic coefficiens: row oriented matrix for each constraint
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - ncons number of added constraints
 *  - lhss left hand sides of constraints
 *  - rhss right hand sides of constraints
 *  - nlininds number of linear coefficients for each constraint
 *    may be NULL in case of no linear part
 *  - lininds indices of variables for linear coefficients for each constraint
 *    may be NULL in case of no linear part
 *  - linvals values of linear coefficient for each constraint
 *    may be NULL in case of no linear part
 *  - nquadrows number of columns in matrix of quadratic part for each constraint
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadrowidxs indices of variables for which a quadratic part is specified
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadoffsets start index of each rows quadratic coefficients in quadinds[.] and quadvals[.]
 *    indices are given w.r.t. quadrowidxs., i.e., quadoffsets[.][i] gives the start index of row quadrowidxs[.][i] in quadvals[.]
 *    quadoffsets[.][nquadrows[.]] gives length of quadinds[.] and quadvals[.]
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadinds column indices w.r.t. quadrowidxs, i.e., quadrowidxs[quadinds[.][i]] gives the index of the variable corresponding
 *    to entry i, entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - quadvals coefficient values
 *    entry of array may be NULL in case of no quadratic part
 *    may be NULL in case of no quadratic part in any constraint
 *  - exprvaridxs indices of variables in expression tree, maps variable indices in expression tree to indices in nlp
 *    entry of array may be NULL in case of no expression tree
 *    may be NULL in case of no expression tree in any constraint
 *  - exprtrees expression tree for nonquadratic part of constraints
 *    entry of array may be NULL in case of no nonquadratic part
 *    may be NULL in case of no nonquadratic part in any constraint
 *  - names of constraints, may be NULL or entries may be NULL
 */
static
SCIP_DECL_NLPIADDCONSTRAINTS( nlpiAddConstraintsWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleAddConstraints(problem->oracle,
         ncons, lhss, rhss,
         nlininds, lininds, linvals,
         nquadelems, quadelems,
         exprvaridxs, exprtrees, names) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets or overwrites objective, a minimization problem is expected
 *  May change sparsity pattern.
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nlins number of linear variables
 *  - lininds variable indices
 *    may be NULL in case of no linear part
 *  - linvals coefficient values
 *    may be NULL in case of no linear part
 *  - nquadcols number of columns in matrix of quadratic part
 *  - quadcols indices of variables for which a quadratic part is specified
 *    may be NULL in case of no quadratic part
 *  - quadoffsets start index of each rows quadratic coefficients in quadinds and quadvals
 *    quadoffsets[.][nquadcols] gives length of quadinds and quadvals
 *    may be NULL in case of no quadratic part
 *  - quadinds column indices
 *    may be NULL in case of no quadratic part
 *  - quadvals coefficient values
 *    may be NULL in case of no quadratic part
 *  - exprvaridxs indices of variables in expression tree, maps variable indices in expression tree to indices in nlp
 *    may be NULL in case of no expression tree
 *  - exprtree expression tree for nonquadratic part of objective function
 *    may be NULL in case of no nonquadratic part
 *  - constant objective value offset
 */
static
SCIP_DECL_NLPISETOBJECTIVE( nlpiSetObjectiveWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleSetObjective(problem->oracle,
         constant, nlins, lininds, linvals,
         nquadelems, quadelems,
         exprvaridxs, exprtree) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change variable bounds
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nvars number of variables to change bounds
 *  - indices indices of variables to change bounds
 *  - lbs new lower bounds
 *  - ubs new upper bounds
 */
static
SCIP_DECL_NLPICHGVARBOUNDS( nlpiChgVarBoundsWorhp )
{
   const SCIP_Real* oldlbs = SCIPnlpiOracleGetVarLbs(problem->oracle);
   const SCIP_Real* oldubs = SCIPnlpiOracleGetVarUbs(problem->oracle);
   int index;
   int i;

   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   /* TODO: So far, Worhp can not handle fixed variables (and fixed variables that have been unfixed) when applying a
    * restart. The following code needs to be removed when this has changed.
    */
   for( i = 0; i < nvars; ++i )
   {
      index = indices[i];
      SCIPdebugMessage("change bounds of %d from [%g,%g] -> [%g,%g]\n", index, oldlbs[index], oldubs[index],
         lbs[i], ubs[i]);

      if( REALABS(lbs[i] - ubs[i]) <= problem->feastol )
         problem->firstrun = TRUE;
      else
         if( REALABS(oldlbs[index] - oldubs[index]) <= problem->feastol )
            problem->firstrun = TRUE;
   }

   SCIP_CALL( SCIPnlpiOracleChgVarBounds(problem->oracle, nvars, indices, lbs, ubs) );

   invalidateSolution(problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change constraint bounds
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - nconss number of constraints to change sides
 *  - indices indices of constraints to change sides
 *  - lhss new left hand sides
 *  - rhss new right hand sides
 */
static
SCIP_DECL_NLPICHGCONSSIDES( nlpiChgConsSidesWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

#ifdef SCIP_DEBUG
   {
      SCIP_Real oldlhs;
      SCIP_Real oldrhs;
      int index;
      int i;

      for( i = 0; i < nconss; ++i )
      {
         index = indices[i];
         oldlhs = SCIPnlpiOracleGetConstraintLhs(problem->oracle, index);
         oldrhs = SCIPnlpiOracleGetConstraintRhs(problem->oracle, index);
         SCIPdebugMessage("change constraint side of %d from [%g,%g] -> [%g,%g]\n", index, oldlhs, oldrhs, lhss[i], rhss[i]);
      }
   }
#endif

   SCIP_CALL( SCIPnlpiOracleChgConsSides(problem->oracle, nconss, indices, lhss, rhss) );

   invalidateSolution(problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** delete a set of variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - dstats deletion status of vars; 1 if var should be deleted, 0 if not
 *
 * output:
 *  - dstats new position of var, -1 if var was deleted
 */
static
SCIP_DECL_NLPIDELVARSET( nlpiDelVarSetWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelVarSet(problem->oracle, dstats) );

   BMSfreeMemoryArrayNull(&problem->initguess); // @TODO keep initguess for remaining variables

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** delete a set of constraints
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - dstats deletion status of rows; 1 if row should be deleted, 0 if not
 *
 * output:
 *  - dstats new position of row, -1 if row was deleted
 */
static
SCIP_DECL_NLPIDELCONSSET( nlpiDelConstraintSetWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleDelConsSet(problem->oracle, dstats) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** changes (or adds) linear coefficients in a constraint or objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idx index of constraint or -1 for objective
 *  - nvals number of values in linear constraint to change
 *  - varidxs indices of variables which coefficient to change
 *  - vals new values for coefficients
 */
static
SCIP_DECL_NLPICHGLINEARCOEFS( nlpiChgLinearCoefsWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgLinearCoefs(problem->oracle, idx, nvals, varidxs, vals) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** changes (or adds) coefficients in the quadratic part of a constraint or objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idx index of constraint or -1 for objective
 *  - nentries number of entries in quadratic matrix to change
 *  - rows row indices of entries in quadratic matrix where values should be changed
 *  - cols column indices of entries in quadratic matrix where values should be changed
 *  - values new values for entries in quadratic matrix
 */
static
SCIP_DECL_NLPICHGQUADCOEFS( nlpiChgQuadraticCoefsWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgQuadCoefs(problem->oracle, idx, nquadelems, quadelems) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** replaces the expression tree of a constraint or objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idxcons index of constraint or -1 for objective
 *  - exprvaridxs indices of variables in expression tree, maps variable indices in expression tree to indices in nlp, or NULL
 *  - exprtree new expression tree for constraint or objective, or NULL to only remove previous tree
 */
static
SCIP_DECL_NLPICHGEXPRTREE( nlpiChgExprtreeWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgExprtree(problem->oracle, idxcons, exprvaridxs, exprtree) );

   invalidateSolution(problem);
   problem->firstrun = TRUE;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change one coefficient in the nonlinear part
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - idxcons index of constraint or -1 for objective
 *  - idxparam index of parameter
 *  - value new value for nonlinear parameter
 *
 * return: Error if parameter does not exist
 */
static
SCIP_DECL_NLPICHGNONLINCOEF( nlpiChgNonlinCoefWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgExprParam(problem->oracle, idxcons, idxparam, value) );

   invalidateSolution(problem);

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** change the constant offset in the objective
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - objconstant new value for objective constant
 */
static
SCIP_DECL_NLPICHGOBJCONSTANT( nlpiChgObjConstantWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   SCIP_CALL( SCIPnlpiOracleChgObjConstant(problem->oracle, objconstant) );

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets initial guess for primal variables
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues initial primal values for variables, or NULL to clear previous values
 *  - consdualvalues initial dual values for constraints, or NULL to clear previous values
 *  - varlbdualvalues  initial dual values for variable lower bounds, or NULL to clear previous values
 *  - varubdualvalues  initial dual values for variable upper bounds, or NULL to clear previous values
 */
static
SCIP_DECL_NLPISETINITIALGUESS( nlpiSetInitialGuessWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);
   assert(problem->oracle != NULL);

   if( primalvalues != NULL )
   {
      if( !problem->initguess )
      {
         if( BMSduplicateMemoryArray(&problem->initguess, primalvalues, SCIPnlpiOracleGetNVars(problem->oracle)) == NULL )
            return SCIP_NOMEMORY;
      }
      else
      {
         BMScopyMemoryArray(problem->initguess, primalvalues, SCIPnlpiOracleGetNVars(problem->oracle));
      }
   }
   else
   {
      BMSfreeMemoryArrayNull(&problem->initguess);
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** tries to solve NLP
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 */
static
SCIP_DECL_NLPISOLVE( nlpiSolveWorhp )
{
   SCIP_NLPIDATA* nlpidata = SCIPnlpiGetData(nlpi);
   Workspace* wsp = problem->wsp;
   Control* cnt = problem->cnt;
   OptVar* opt = problem->opt;
   Params* par = problem->par;
   int status;
   int i;

   problem->lastniter = -1;
   problem->lasttime  = -1.0;

   /* initialize Worhp data if necessary */
   if( problem->firstrun )
   {
      SCIP_CALL( freeWorhp(problem) );
      SCIP_CALL( initWorhp(nlpi, problem) );
      problem->firstrun = FALSE;
   }
   else
   {
      SCIP_CALL( updateWorhp(problem) );
   }

   /* set parameters */
   InitParams(&status, par);

   if( status != OK )
      return SCIP_INVALIDCALL;

   par->algorithm = DEFAULT_ALGORITHM;
   par->ScaledKKT = DEFAULT_SCALEDKKT;
   par->Infty = nlpidata->infinity;
   par->TolFeas = problem->feastol;
   par->TolOpti = problem->relobjtol;
   par->TolComp = problem->relobjtol;
   par->Timeout = problem->timelim;
   par->MaxIter = problem->itlim;
   par->NLPprint = problem->verblevel - 1; /* Worhp verbosity levels: -1 = off, 0 = normal, 1 = debug, >1 = more debug */

#ifdef CHECKFUNVALUES
   /* activate gradient and hessian check */
   par->CheckValuesDF = TRUE;
   par->CheckValuesDG = TRUE;
   par->CheckValuesHM = TRUE;
#endif

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPnlpiOraclePrintProblem(problem->oracle, nlpidata->messagehdlr, NULL) );
#endif

   /* set initial guess (if available) */
   if( problem->initguess != NULL )
   {
      for( i = 0; i < problem->opt->n; ++i )
      {
         problem->opt->X[i] = problem->initguess[i];
      }
   }

   /*
    * Worhp Reverse Communication loop.
    * In every iteration poll GetUserAction for the requested action, i.e. one
    * of {callWorhp, iterOutput, evalF, evalG, evalDF, evalDG, evalHM, fidif}.
    *
    * Make sure to reset the requested user action afterwards by calling
    * DoneUserAction, except for 'callWorhp' and 'fidif'.
    */
   while( cnt->status < TerminateSuccess && cnt->status > TerminateError )
   {
      /*
       * Worhp's main routine.
       * Do not manually reset callWorhp, this is only done by the FD routines.
       */
      if( GetUserAction(cnt, callWorhp) )
      {
         Worhp(opt, wsp, par, cnt);
         /* No DoneUserAction! */
      }

      /*
       * Show iteration output.
       * The call to IterationOutput() may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, iterOutput) )
      {
         IterationOutput(opt, wsp, par, cnt);
         DoneUserAction(cnt, iterOutput);
      }

      /*
       * Evaluate the objective function.
       * The call to UserF may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalF) )
      {
         if( userF(problem) != SCIP_OKAY )
            break;
         DoneUserAction(cnt, evalF);
      }

      /*
       * Evaluate the constraints.
       * The call to UserG may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalG) )
      {
         if( userG(problem) != SCIP_OKAY )
            break;
         DoneUserAction(cnt, evalG);
      }

      /*
       * Evaluate the gradient of the objective function.
       * The call to UserDF may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalDF) )
      {
         if( userDF(problem) != SCIP_OKAY )
            break;
         DoneUserAction(cnt, evalDF);
      }

      /*
       * Evaluate the Jacobian of the constraints.
       * The call to UserDG may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalDG) )
      {
         if( userDG(problem) != SCIP_OKAY )
            break;
         DoneUserAction(cnt, evalDG);
      }

      /*
       * Evaluate the Hessian matrix of the Lagrange function (L = f + mu*g)
       * The call to UserHM may be replaced by user-defined code.
       */
      if( GetUserAction(cnt, evalHM) )
      {
         if( userHM(problem) != SCIP_OKAY)
            break;
         DoneUserAction(cnt, evalHM);
      }

      /*
       * Use finite differences with RC to determine derivatives
       * Do not reset fidif, this is done by the FD routine.
       */
      if( GetUserAction(cnt, fidif) )
      {
         WorhpFidif(opt, wsp, par, cnt);
         /* No DoneUserAction! */
      }
   }

   /* interpret Worhp result */
   if( cnt->status < TerminateSuccess && cnt->status > TerminateError )
   {
      SCIPmessagePrintWarning(nlpidata->messagehdlr, "Worhp failed because of an invalid function evaluation!\n");
      problem->lastsolstat  = SCIP_NLPSOLSTAT_UNKNOWN;
      problem->lasttermstat = SCIP_NLPTERMSTAT_NUMERR;
   }
   else
   {
      SCIP_CALL( evaluateWorhpRun(problem, nlpidata->messagehdlr) );
   }

   /* free memory */
   StatusMsg(opt, wsp, par, cnt);

   /* store statistics */
   problem->lastniter = wsp->MajorIter;
   problem->lasttime = GetTimerCont(&cnt->Timer);

   return SCIP_OKAY;
}  /*lint !e715*/

/** gives solution status
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 * return: Solution Status
 */
static
SCIP_DECL_NLPIGETSOLSTAT( nlpiGetSolstatWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   return problem->lastsolstat;
}  /*lint !e715*/

/** gives termination reason
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *
 * return: Termination Status
 */
static
SCIP_DECL_NLPIGETTERMSTAT( nlpiGetTermstatWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   return problem->lasttermstat;
}  /*lint !e715*/

/** gives primal and dual solution values
 *
 * solver can return NULL in dual values if not available
 * but if solver provides dual values for one side of variable bounds, then it must also provide those for the other side
 *
 * for a ranged constraint, the dual variable is positive if the right hand side is active and negative if the left hand side is active
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - primalvalues buffer to store pointer to array to primal values, or NULL if not needed
 *  - consdualvalues buffer to store pointer to array to dual values of constraints, or NULL if not needed
 *  - varlbdualvalues buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed
 *  - varubdualvalues buffer to store pointer to array to dual values of variable lower bounds, or NULL if not needed
 */
static
SCIP_DECL_NLPIGETSOLUTION( nlpiGetSolutionWorhp )
{
   assert(problem != NULL);

   if( primalvalues != NULL )
      *primalvalues = problem->lastprimal;

   if( consdualvalues != NULL )
      *consdualvalues = problem->lastdualcons;

   if( varlbdualvalues != NULL )
      *varlbdualvalues = problem->lastduallb;

   if( varubdualvalues != NULL )
      *varubdualvalues = problem->lastdualub;

   return SCIP_OKAY;
}  /*lint !e715*/

/** gives solve statistics
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - statistics pointer to store statistics
 *
 * output:
 *  - statistics solve statistics
 */
static
SCIP_DECL_NLPIGETSTATISTICS( nlpiGetStatisticsWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   SCIPnlpStatisticsSetNIterations(statistics, problem->lastniter);
   SCIPnlpStatisticsSetTotalTime(statistics, problem->lasttime);

   return SCIP_OKAY;
}  /*lint !e715*/

/** gives required size of a buffer to store a warmstart object
 *
 *  input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - size pointer to store required size for warmstart buffer
 *
 * output:
 *  - size required size for warmstart buffer
 */
static
SCIP_DECL_NLPIGETWARMSTARTSIZE( nlpiGetWarmstartSizeWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** stores warmstart information in buffer
 *
 * required size of buffer should have been obtained by SCIPnlpiGetWarmstartSize before
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - buffer memory to store warmstart information
 *
 * output:
 *  - buffer warmstart information in solver specific data structure
 */
static
SCIP_DECL_NLPIGETWARMSTARTMEMO( nlpiGetWarmstartMemoWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets warmstart information in solver
 *
 * write warmstart to buffer
 *
 * input:
 *  - nlpi datastructure for solver interface
 *  - problem datastructure for problem instance
 *  - buffer warmstart information
 */
static
SCIP_DECL_NLPISETWARMSTARTMEMO( nlpiSetWarmstartMemoWorhp )
{
   /* TODO */

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gets integer parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival pointer to store the parameter value
 *
 * output:
 *  - ival parameter value
 */
static
SCIP_DECL_NLPIGETINTPAR( nlpiGetIntParWorhp )
{
   assert(nlpi != NULL);
   assert(ival != NULL);
   assert(problem != NULL);

   switch( type )
   {
   case SCIP_NLPPAR_FROMSCRATCH:
   {
      *ival = 1;
      break;
   }

   case SCIP_NLPPAR_VERBLEVEL:
   {
      *ival = problem->verblevel;
      break;
   }

   case SCIP_NLPPAR_FEASTOL:
   {
      SCIPerrorMessage("feasibility tolerance parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_RELOBJTOL:
   {
      SCIPerrorMessage("relative objective tolerance parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_LOBJLIM:
   {
      SCIPerrorMessage("objective limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_INFINITY:
   {
      SCIPerrorMessage("infinity parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_ITLIM:
   {
      *ival = problem->itlim;
      break;
   }

   case SCIP_NLPPAR_TILIM:
   {
      SCIPerrorMessage("time limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_OPTFILE:
   {
      SCIPerrorMessage("optfile parameter is of type string.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FASTFAIL:
   {
      *ival = problem->fastfail ? 1 : 0;
      break;
   }

   default:
   {
      SCIPerrorMessage("Parameter %d not known to Worhp interface.\n", type);
      return SCIP_PARAMETERUNKNOWN;
   }
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** sets integer parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - ival parameter value
 */
static
SCIP_DECL_NLPISETINTPAR( nlpiSetIntParWorhp )
{
   assert(nlpi != NULL);
   assert(problem != NULL);

   switch( type )
   {
   case SCIP_NLPPAR_FROMSCRATCH:
   {
      if( ival == 0 || ival == 1 )
      {
         SCIP_NLPIDATA* data;

         data = SCIPnlpiGetData(nlpi);
         assert(data != NULL);

         SCIPmessagePrintWarning(data->messagehdlr, "from scratch parameter not supported by Worhp interface yet. Ignored.\n");
      }
      else
      {
         SCIPerrorMessage("Value %d for parameter from scratch out of range {0, 1}\n", ival);
         return SCIP_PARAMETERWRONGVAL;
      }
      break;
   }

   case SCIP_NLPPAR_VERBLEVEL:
   {
      assert(ival >= 0);
      problem->verblevel = ival;
      break;
   }

   case SCIP_NLPPAR_FEASTOL:
   {
      SCIPerrorMessage("feasibility tolerance parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_RELOBJTOL:
   {
      SCIPerrorMessage("relative objective tolerance parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_LOBJLIM:
   {
      SCIPerrorMessage("objective limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_INFINITY:
   {
      SCIPerrorMessage("infinity parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_ITLIM:
   {
      if( ival >= 0 )
         problem->itlim = ival;
      else
      {
         SCIPerrorMessage("Value %d for parameter iteration limit is negative\n", ival);
         return SCIP_PARAMETERWRONGVAL;
      }
      break;
   }

   case SCIP_NLPPAR_TILIM:
   {
      SCIPerrorMessage("time limit parameter is of type real.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_OPTFILE:
   {
      SCIPerrorMessage("optfile parameter is of type string.\n");
      return SCIP_PARAMETERWRONGTYPE;
   }

   case SCIP_NLPPAR_FASTFAIL:
   {
      if( ival == 0 || ival == 1 )
      {
         problem->fastfail = ival;
      }
      else
      {
         SCIPerrorMessage("Value %d for parameter fastfail out of range {0, 1}\n", ival);
         return SCIP_PARAMETERWRONGVAL;
      }
      break;
   }

   default:
   {
      SCIPerrorMessage("Parameter %d not known to Worhp interface.\n", type);
      return SCIP_PARAMETERUNKNOWN;
   }
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets floating point parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval pointer to store the parameter value
 *
 * output:
 *  - dval parameter value
 */
static
SCIP_DECL_NLPIGETREALPAR( nlpiGetRealParWorhp )
{
   SCIP_NLPIDATA* data = SCIPnlpiGetData(nlpi);

   assert(data != NULL);
   assert(dval != NULL);

   switch( type )
   {
      case SCIP_NLPPAR_FROMSCRATCH:
      {
         SCIPerrorMessage("fromscratch parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_VERBLEVEL:
      {
         SCIPerrorMessage("verblevel parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_FEASTOL:
      {
         *dval = problem->feastol;
         break;
      }

      case SCIP_NLPPAR_RELOBJTOL:
      {
         *dval = problem->relobjtol;
         break;
      }

      case SCIP_NLPPAR_LOBJLIM:
      {
         *dval = problem->lobjlim;
         break;
      }

      case SCIP_NLPPAR_INFINITY:
      {
         *dval = data->infinity;
         break;
      }

      case SCIP_NLPPAR_ITLIM:
      {
         SCIPerrorMessage("itlim parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_TILIM:
      {
         *dval = problem->timelim;
         break;
      }

      case SCIP_NLPPAR_OPTFILE:
      {
         SCIPerrorMessage("optfile parameter is of type string.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_FASTFAIL:
      {
         SCIPerrorMessage("fastfail parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      default:
      {
         break;
      }
   }

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets floating point parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance, can be NULL only if type == SCIP_NLPPAR_INFINITY
 *  - type parameter number
 *  - dval parameter value
 */
static
SCIP_DECL_NLPISETREALPAR( nlpiSetRealParWorhp )
{
   SCIP_NLPIDATA* data = SCIPnlpiGetData(nlpi);

   assert(data != NULL);

   switch( type )
   {
      case SCIP_NLPPAR_FROMSCRATCH:
      {
         SCIPerrorMessage("fromscratch parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_VERBLEVEL:
      {
         SCIPerrorMessage("verblevel parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_FEASTOL:
      {
         problem->feastol = dval;
         break;
      }

      case SCIP_NLPPAR_RELOBJTOL:
      {
         problem->relobjtol = dval;
         break;
      }

      case SCIP_NLPPAR_LOBJLIM:
      {
         problem->lobjlim = dval;
         break;
      }

      case SCIP_NLPPAR_INFINITY:
      {
         data->infinity = dval;
         break;
      }

      case SCIP_NLPPAR_ITLIM:
      {
         SCIPerrorMessage("itlim parameter is of type real.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_TILIM:
      {
         problem->timelim = dval;
         break;
      }

      case SCIP_NLPPAR_OPTFILE:
      {
         SCIPerrorMessage("optfile parameter is of type string.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      case SCIP_NLPPAR_FASTFAIL:
      {
         SCIPerrorMessage("optfile parameter is of type int.\n");
         return SCIP_PARAMETERWRONGTYPE;
      }

      default:
      {
         break;
      }
   }

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** gets string parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - sval pointer to store the string value, the user must not modify the string
 *
 * output:
 *  - sval parameter value
 */
static
SCIP_DECL_NLPIGETSTRINGPAR( nlpiGetStringParWorhp )
{
   SCIP_NLPIDATA* nlpidata = SCIPnlpiGetData(nlpi);

   if( type == SCIP_NLPPAR_OPTFILE )
   {
      SCIPmessagePrintWarning(nlpidata->messagehdlr, "optfile parameter not supported by Worhp interface yet. Ignored.\n");
   }
   else
   {
      SCIPerrorMessage("parameter %d is not of type string.\n", type);
      return SCIP_PARAMETERWRONGTYPE;
   }

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets string parameter of NLP
 *
 * input:
 *  - nlpi NLP interface structure
 *  - problem datastructure for problem instance
 *  - type parameter number
 *  - sval parameter value
 */
static
SCIP_DECL_NLPISETSTRINGPAR( nlpiSetStringParWorhp )
{
   SCIP_NLPIDATA* nlpidata = SCIPnlpiGetData(nlpi);

   if( type == SCIP_NLPPAR_OPTFILE )
   {
      SCIPmessagePrintWarning(nlpidata->messagehdlr, "optfile parameter not supported by Worhp interface yet. Ignored.\n");
   }
   else
   {
      SCIPerrorMessage("parameter %d is not of type string.\n", type);
      return SCIP_PARAMETERWRONGTYPE;
   }

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/** sets message handler for message output
 *
 * input:
 *  - nlpi NLP interface structure
 *  - messagehdlr SCIP message handler, or NULL to suppress all output
 */
static
SCIP_DECL_NLPISETMESSAGEHDLR( nlpiSetMessageHdlrWorhp )
{
   SCIP_NLPIDATA* nlpidata;

   assert(nlpi != NULL);

   nlpidata = SCIPnlpiGetData(nlpi);
   assert(nlpidata != NULL);

   nlpidata->messagehdlr = messagehdlr;

   return SCIP_OKAY;  /*lint !e527*/
}  /*lint !e715*/

/*
 * NLP solver interface specific interface methods
 */

/** create solver interface for Worhp solver */
SCIP_RETCODE SCIPcreateNlpSolverWorhp(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi                /**< pointer to buffer for nlpi address */
   )
{
   SCIP_NLPIDATA* nlpidata;

   assert(blkmem != NULL);
   assert(nlpi   != NULL);

   /* create Worhp solver interface data */
   BMSallocMemory(&nlpidata);
   BMSclearMemory(nlpidata);

   nlpidata->blkmem = blkmem;

   /* initialize parameter */
   nlpidata->infinity = SCIP_DEFAULT_INFINITY;

   /* checks the version of the library and header files */
   CHECK_WORHP_VERSION

   /* create solver interface */
   SCIP_CALL( SCIPnlpiCreate(nlpi,
         NLPI_NAME, NLPI_DESC, NLPI_PRIORITY,
         nlpiCopyWorhp, nlpiFreeWorhp, nlpiGetSolverPointerWorhp,
         nlpiCreateProblemWorhp, nlpiFreeProblemWorhp, nlpiGetProblemPointerWorhp,
         nlpiAddVarsWorhp, nlpiAddConstraintsWorhp, nlpiSetObjectiveWorhp,
         nlpiChgVarBoundsWorhp, nlpiChgConsSidesWorhp, nlpiDelVarSetWorhp, nlpiDelConstraintSetWorhp,
         nlpiChgLinearCoefsWorhp, nlpiChgQuadraticCoefsWorhp, nlpiChgExprtreeWorhp, nlpiChgNonlinCoefWorhp,
         nlpiChgObjConstantWorhp, nlpiSetInitialGuessWorhp, nlpiSolveWorhp, nlpiGetSolstatWorhp, nlpiGetTermstatWorhp,
         nlpiGetSolutionWorhp, nlpiGetStatisticsWorhp,
         nlpiGetWarmstartSizeWorhp, nlpiGetWarmstartMemoWorhp, nlpiSetWarmstartMemoWorhp,
         nlpiGetIntParWorhp, nlpiSetIntParWorhp, nlpiGetRealParWorhp, nlpiSetRealParWorhp, nlpiGetStringParWorhp, nlpiSetStringParWorhp,
         nlpiSetMessageHdlrWorhp,
         nlpidata) );

   return SCIP_OKAY;
}

/** gets string that identifies Worhp (version number) */
const char* SCIPgetSolverNameWorhp(
   void
   )
{
   return "WORHP " WORHP_VERSION;
}

/** gets string that describes Worhp (version number) */
const char* SCIPgetSolverDescWorhp(
   void
   )
{
   return "Sequential Quadratic Programming developed at Research Institute Steinbeis (www.worhp.de)";
}

/** returns whether Worhp is available, i.e., whether it has been linked in */
SCIP_Bool SCIPisWorhpAvailableWorhp(
   void
   )
{
   return TRUE;
}

#include <iostream>
#include <scip/message_default.h>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#define EPS 1e-6

/** macro to keep control of infinity
 *
 *  Some LPIs use std::numeric_limits<SCIP_Real>::infinity() as finity value. Comparing two infinty values then yields
 *  nan. This is a workaround. */
#define cr_assert_float_eq(Actual, Expected, Epsilon, FormatString, ...) \
      assert(ABS(Actual - Expected) < Epsilon);
#define cr_assert_float_eq_inf(Actual, Expected, Epsilon, FormatString, ...) \
   if ( ABS(Actual) < 1e30 && ABS(Expected) < 1e30 )                       \
      cr_assert_float_eq(Actual, Expected, Epsilon, FormatString, __VA_ARGS__);

/* GLOBAL VARIABLES */
static SCIP_LPI* lpi = NULL;
static SCIP_MESSAGEHDLR* messagehdlr = NULL;

/** expected feasibility status for primal or dual problem */
enum SCIPfeasStatus
{
    SCIPfeas = 0,    /**< the problem is feasible */
    SCIPunbounded = 1,    /**< the problem is unbounded */
    SCIPinfeas = 2     /**< the problem is infeasible */
};
typedef enum SCIPfeasStatus SCIPFEASSTATUS;

/** solve problem */
static
SCIP_RETCODE solveTest(
    SCIP_Bool             solveprimal,        /**< use primal simplex */
    int                   ncols,              /**< number of columns */
    int                   nrows,              /**< number of rows */
    SCIPFEASSTATUS        exp_primalfeas,     /**< expected primal feasibility status */
    SCIPFEASSTATUS        exp_dualfeas,       /**< expected primal feasibility status */
    const SCIP_Real * exp_primsol,        /**< expected primal optimal solution or primal ray if primal is unbounded or NULL */
    const SCIP_Real * exp_dualsol,        /**< expected dual optimal solution or dual ray if dual is unbounded or NULL */
    const SCIP_Real * exp_activity,       /**< expected activity of optimal solution or NULL */
    const SCIP_Real * exp_redcost         /**< expected reduced cost of optimal solution or NULL */
)
{
    /* solution data */
    SCIP_Real objval;
    SCIP_Real* primsol;
    SCIP_Real* dualsol;
    SCIP_Real* activity;
    SCIP_Real* redcost;

    /* auxiliary data */
    SCIP_Bool primalfeasible;
    SCIP_Bool dualfeasible;
    int ntmprows;
    int ntmpcols;
    int i;
    int j;

    /* check size */
    SCIP_CALL(SCIPlpiGetNRows(lpi, &ntmprows));
    SCIP_CALL(SCIPlpiGetNCols(lpi, &ntmpcols));
    assert(nrows == ntmprows);
    assert(ncols == ntmpcols);

    /* solve problem */
    if (solveprimal)
    {
        SCIP_CALL(SCIPlpiSolvePrimal(lpi));
    }
    else
    {
        SCIP_CALL(SCIPlpiSolveDual(lpi));
    }

    /* check status */
    assert(SCIPlpiWasSolved(lpi));
    assert(!SCIPlpiIsObjlimExc(lpi));
    assert(!SCIPlpiIsIterlimExc(lpi));
    assert(!SCIPlpiIsTimelimExc(lpi));

    /* check feasibility status */
    SCIP_CALL(SCIPlpiGetSolFeasibility(lpi, &primalfeasible, &dualfeasible));

    /* if we are feasible, we should be optimal */
    if (exp_primalfeas == SCIPfeas && exp_dualfeas == SCIPfeas)
    {
        assert(SCIPlpiIsOptimal(lpi));
    }

    /* check more primal statuses */
    switch (exp_primalfeas)
    {
    case SCIPfeas:
        assert(primalfeasible);
        assert(!SCIPlpiExistsPrimalRay(lpi));
        assert(!SCIPlpiHasPrimalRay(lpi));
        assert(!SCIPlpiIsPrimalUnbounded(lpi));
        assert(!SCIPlpiIsPrimalInfeasible(lpi));
        assert(SCIPlpiIsPrimalFeasible(lpi));
        break;

    case SCIPunbounded:
        /* Because of SoPlex, cannot always determine feasibility status here, even if we want to apply the primal
         * simplex. In any case, the results of primalfeasible and SCIPlpiIsPrimalFeasible(lpi) should coincide. */
        assert(primalfeasible == SCIPlpiIsPrimalFeasible(lpi));

        /* It seems that we cannot guarantee that the primal is shown to be unbounded. */
        /* cr_assert( SCIPlpiIsPrimalUnbounded(lpi) ); */

        /* primal ray should exist if the primal simplex ran */
        assert(!solveprimal || SCIPlpiExistsPrimalRay(lpi));
        assert(!SCIPlpiIsPrimalInfeasible(lpi));
        break;

    case SCIPinfeas:
        assert(!primalfeasible);
        /* It seems that we cannot always prove that primal is infeasible. */
        /* cr_assert( SCIPlpiIsPrimalInfeasible(lpi) ); */

        /* It seems that we cannot always prove that primal is not unbounded. */
        /* cr_assert( ! SCIPlpiIsPrimalUnbounded(lpi) ); */
        assert(!SCIPlpiIsPrimalFeasible(lpi));
        break;

    default:
        abort();
    }

    /* check more dual statuses */
    switch (exp_dualfeas)
    {
    case SCIPfeas:
        assert(dualfeasible);
        assert(!SCIPlpiExistsDualRay(lpi));
        assert(!SCIPlpiHasDualRay(lpi));
        assert(!SCIPlpiIsDualUnbounded(lpi));
        assert(!SCIPlpiIsDualInfeasible(lpi));
        assert(SCIPlpiIsDualFeasible(lpi));
        break;

    case SCIPunbounded:
        /* Because of SoPlex, cannot always determine feasibility status here, even if we want to apply the dual
         * simplex. In any case, the results of dualfeasible and SCIPlpiIsDualFeasible(lpi) should coincide. */
        assert(dualfeasible == SCIPlpiIsDualFeasible(lpi));

        /* It seems that we cannot guarantee that the dual is shown to be unbounded. */
        /* cr_assert( SCIPlpiIsDualUnbounded(lpi) ); */

        /* dual ray should exist if the dual simplex ran */
        assert(solveprimal || SCIPlpiExistsDualRay(lpi));
        assert(!SCIPlpiIsDualInfeasible(lpi));
        break;

    case SCIPinfeas:
        assert(!dualfeasible);
        assert(!SCIPlpiIsDualUnbounded(lpi));
        /* It seems that we cannot always prove that dual is infeasible. */
        /* cr_assert( SCIPlpiIsDualInfeasible(lpi) ); */
        assert(!SCIPlpiIsDualFeasible(lpi));
        break;

    default:
        abort();
    }

    /* allocate storage for solution */
    BMSallocMemoryArray(&primsol, ncols);
    BMSallocMemoryArray(&dualsol, nrows);
    BMSallocMemoryArray(&activity, nrows);
    BMSallocMemoryArray(&redcost, ncols);

    /* check solution */
    if (exp_primalfeas == SCIPfeas)
    {
        /* get solution */
        SCIP_CALL(SCIPlpiGetSol(lpi, &objval, primsol, dualsol, activity, redcost));

        assert(exp_primsol != NULL && exp_redcost != NULL);
        for (j = 0; j < ncols; ++j)
        {
            cr_assert_float_eq(primsol[j], exp_primsol[j], EPS, "Violation of primal solution %d: %g != %g\n", j, primsol[j], exp_primsol[j]);
            cr_assert_float_eq(redcost[j], exp_redcost[j], EPS, "Violation of reduced cost of solution %d: %g != %g\n", j, redcost[j], exp_redcost[j]);
        }
    }
    else if (exp_primalfeas == SCIPunbounded)
    {
        assert(exp_primsol != NULL);

        if (SCIPlpiHasPrimalRay(lpi))
        {
            SCIP_Real scalingfactor = 1.0;

            SCIP_CALL(SCIPlpiGetPrimalRay(lpi, primsol));

            /* loop until scaling factor can be determined */
            for (j = 0; j < ncols; ++j)
            {
                if (REALABS(exp_primsol[j]) < EPS)
                    assert(ABS(primsol[j] - exp_primsol[j]) < EPS);
                    //cr_assert_float_eq(primsol[j], exp_primsol[j], EPS, "Violation of primal ray %d: %g != %g\n", j, primsol[j], exp_primsol[j]);
                else
                {
                    scalingfactor = primsol[j] / exp_primsol[j];
                    break;
                }
            }

            /* again loop over ray */
            for (j = 0; j < ncols; ++j)
            {
                cr_assert_float_eq(primsol[j], scalingfactor * exp_primsol[j], EPS, "Violation of primal ray %d: %g != %g\n", j, primsol[j], scalingfactor * exp_primsol[j]);
            }
        }
    }

    if (exp_dualfeas == SCIPfeas)
    {
        /* get solution */
        SCIP_CALL(SCIPlpiGetSol(lpi, &objval, primsol, dualsol, activity, redcost));

        assert(exp_dualsol != NULL && exp_activity != NULL);
        for (i = 0; i < nrows; ++i)
        {
            cr_assert_float_eq(dualsol[i], exp_dualsol[i], EPS, "Violation of dual solution %d: %g != %g\n", i, dualsol[i], exp_dualsol[i]);
            cr_assert_float_eq(activity[i], exp_activity[i], EPS, "Violation of activity of solution %d: %g != %g\n", i, activity[i], exp_activity[i]);
        }
    }
    else if (exp_dualfeas == SCIPunbounded)
    {
        assert(exp_dualsol != NULL);

        if (SCIPlpiHasDualRay(lpi))
        {
            SCIP_Real scalingfactor = 1.0;
            SCIP_Real* lhs;
            SCIP_Real* rhs;

            /* get lhs/rhs for check of dual ray */
            BMSallocMemoryArray(&lhs, nrows);
            BMSallocMemoryArray(&rhs, nrows);
            SCIP_CALL(SCIPlpiGetSides(lpi, 0, nrows - 1, lhs, rhs));

            /* get dual ray */
            SCIP_CALL(SCIPlpiGetDualfarkas(lpi, dualsol));

            /* loop until scaling factor can be determined */
            for (i = 0; i < nrows; ++i)
            {
                if (REALABS(exp_dualsol[i]) < EPS)
                    assert(ABS(dualsol[i] - exp_dualsol[i]) < EPS);
                    //cr_assert_float_eq(dualsol[i], exp_dualsol[i], EPS, "Violation of dual ray %d: %g != %g\n", i, dualsol[i], exp_dualsol[i]);
                else
                {
                    scalingfactor = dualsol[i] / exp_dualsol[i];
                    break;
                }
            }

            /* again loop over ray */
            for (i = 0; i < nrows; ++i)
            {
                cr_assert_float_eq(dualsol[i], scalingfactor * exp_dualsol[i], EPS, "Violation of dual ray %d: %g != %g\n", i, dualsol[i], scalingfactor * exp_dualsol[i]);
                assert(!SCIPlpiIsInfinity(lpi, -lhs[i]) || dualsol[i] <= -EPS);
                assert(!SCIPlpiIsInfinity(lpi, rhs[i]) || dualsol[i] >= EPS);
            }

            BMSfreeMemoryArray(&rhs);
            BMSfreeMemoryArray(&lhs);
        }
    }

    BMSfreeMemoryArray(&primsol);
    BMSfreeMemoryArray(&dualsol);
    BMSfreeMemoryArray(&activity);
    BMSfreeMemoryArray(&redcost);

    return SCIP_OKAY;
}

/** perform basic test for the given problem */
static
SCIP_RETCODE performTest(
    SCIP_Bool             solveprimal,        /**< use primal simplex */
    SCIP_OBJSEN           objsen,             /**< objective sense */
    int                   ncols,              /**< number of columns */
    const SCIP_Real * obj,                /**< objective function values of columns */
    const SCIP_Real * lb,                 /**< lower bounds of columns */
    const SCIP_Real * ub,                 /**< upper bounds of columns */
    int                   nrows,              /**< number of rows */
    const SCIP_Real * lhs,                /**< left hand sides of rows */
    const SCIP_Real * rhs,                /**< right hand sides of rows */
    int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
    const int* beg,                /**< start index of each column in ind- and val-array */
    const int* ind,                /**< row indices of constraint matrix entries */
    const SCIP_Real * val,                /**< values of constraint matrix entries */
    SCIPFEASSTATUS        exp_primalfeas,     /**< expected primal feasibility status */
    SCIPFEASSTATUS        exp_dualfeas,       /**< expected primal feasibility status */
    const SCIP_Real * exp_primsol,        /**< expected primal optimal solution or primal ray if primal is unbounded or NULL */
    const SCIP_Real * exp_dualsol,        /**< expected dual optimal solution or dual ray if dual is unbounded or NULL */
    const SCIP_Real * exp_activity,       /**< expected activity of optimal solution or NULL */
    const SCIP_Real * exp_redcost         /**< expected reduced cost of optimal solution or NULL */
)
{
    /* load problem */
    SCIP_CALL(SCIPlpiLoadColLP(lpi, objsen, 2, obj, lb, ub, NULL, 2, lhs, rhs, NULL, 4, beg, ind, val));
    assert(!SCIPlpiWasSolved(lpi));

    /* solve problem */
    SCIP_CALL(solveTest(solveprimal, ncols, nrows, exp_primalfeas, exp_dualfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost));

    return SCIP_OKAY;
}

/** check whether data in LP solver aggrees with original data */
static
SCIP_RETCODE checkData(
    SCIP_OBJSEN           objsen,             /**< objective sense */
    int                   ncols,              /**< number of columns */
    const SCIP_Real * obj,                /**< objective function values of columns */
    const SCIP_Real * lb,                 /**< lower bounds of columns */
    const SCIP_Real * ub,                 /**< upper bounds of columns */
    int                   nrows,              /**< number of rows */
    const SCIP_Real * lhs,                /**< left hand sides of rows */
    const SCIP_Real * rhs,                /**< right hand sides of rows */
    int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
    const int* beg,                /**< start index of each column in ind- and val-array */
    const int* ind,                /**< row indices of constraint matrix entries */
    const SCIP_Real * val                 /**< values of constraint matrix entries */
)
{
    SCIP_OBJSEN lpiobjsen;
    SCIP_Real* lpival;
    SCIP_Real* lpilb;
    SCIP_Real* lpiub;
    SCIP_Real* lpiobj;
    SCIP_Real* lpilhs;
    SCIP_Real* lpirhs;
    int* lpibeg;
    int* lpiind;
    int lpincols;
    int lpinrows;
    int lpinnonz;
    int lpinnonz2;
    int i;
    int j;

    /* check number of rows and columns */
    SCIP_CALL(SCIPlpiGetNRows(lpi, &lpinrows));
    SCIP_CALL(SCIPlpiGetNCols(lpi, &lpincols));
    assert(lpinrows == nrows);
    assert(lpincols == ncols);

    /* check objective sense */
    SCIP_CALL(SCIPlpiGetObjsen(lpi, &lpiobjsen));
    assert(objsen == lpiobjsen);

    /* get number of nonzeros in matrix */
    SCIP_CALL(SCIPlpiGetNNonz(lpi, &lpinnonz));
    assert(lpinnonz == nnonz);

    /* allocate storage for data */
    BMSallocMemoryArray(&lpilb, ncols);
    BMSallocMemoryArray(&lpiub, ncols);
    BMSallocMemoryArray(&lpibeg, ncols);
    BMSallocMemoryArray(&lpiind, lpinnonz);
    BMSallocMemoryArray(&lpival, lpinnonz);
    BMSallocMemoryArray(&lpiobj, ncols);

    /* get matrix data */
    SCIP_CALL(SCIPlpiGetCols(lpi, 0, ncols - 1, lpilb, lpiub, &lpinnonz2, lpibeg, lpiind, lpival));
    SCIP_CALL(SCIPlpiGetObj(lpi, 0, ncols - 1, lpiobj));

    /* compare data */
    for (j = 0; j < ncols; ++j)
    {
        cr_assert_float_eq_inf(lpilb[j], lb[j], EPS, "Violation of lower bound %d: %g != %g\n", j, lpilb[j], lb[j]);
        cr_assert_float_eq_inf(lpiub[j], ub[j], EPS, "Violation of upper bound %d: %g != %g\n", j, lpiub[j], ub[j]);

        cr_assert_float_eq(lpiobj[j], obj[j], EPS, "Violation of objective coefficient %d: %g != %g\n", j, lpiobj[j], obj[j]);

        assert(lpibeg[j] == beg[j]);
    }

    /* compare matrix */
    for (j = 0; j < nnonz; ++j)
    {
        assert(lpiind[j] == ind[j]);
        cr_assert_float_eq(lpival[j], val[j], EPS, "Violation of matrix entry (%d, %d): %g != %g\n", ind[j], j, lpival[j], val[j]);
    }

    BMSfreeMemoryArray(&lpiobj);
    BMSfreeMemoryArray(&lpival);
    BMSfreeMemoryArray(&lpiind);
    BMSfreeMemoryArray(&lpibeg);
    BMSfreeMemoryArray(&lpiub);
    BMSfreeMemoryArray(&lpilb);

    /* compare lhs/rhs */
    BMSallocMemoryArray(&lpilhs, nrows);
    BMSallocMemoryArray(&lpirhs, nrows);

    SCIP_CALL(SCIPlpiGetSides(lpi, 0, nrows - 1, lpilhs, lpirhs));

    for (i = 0; i < nrows; ++i)
    {
        cr_assert_float_eq_inf(lpilhs[i], lhs[i], EPS, "Violation of lhs %d: %g != %g\n", i, lpilhs[i], lhs[i]);
        cr_assert_float_eq_inf(lpirhs[i], rhs[i], EPS, "Violation of rhs %d: %g != %g\n", i, lpirhs[i], rhs[i]);
    }

    BMSfreeMemoryArray(&lpirhs);
    BMSfreeMemoryArray(&lpilhs);

    return SCIP_OKAY;
}

SCIP_RETCODE execmain_test1(int argc, const char** argv) {
    /* data with fixed values: */
    SCIP_Real obj[2] = { 3, 1 };
    SCIP_Real lb[2] = { 0, 0 };
    SCIP_Real rhs[2] = { 10, 15 };
    int beg[2] = { 0, 2 };
    int ind[4] = { 0, 1, 0, 1 };
    SCIP_Real val[4] = { 2, 1, 1, 3 };

    /* data to be filled */
    SCIP_Real ub[2];
    SCIP_Real lhs[2];

    /* expected solutions */
    SCIP_Real exp_primsol[2] = { 5, 0 };
    SCIP_Real exp_dualsol[2] = { 1.5, 0 };
    SCIP_Real exp_activity[2] = { 10, 5 };
    SCIP_Real exp_redcost[2] = { 0, -0.5 };

    /* fill variable data */
    ub[0] = SCIPlpiInfinity(lpi);
    ub[1] = SCIPlpiInfinity(lpi);
    lhs[0] = -SCIPlpiInfinity(lpi);
    lhs[1] = -SCIPlpiInfinity(lpi);

    /* solve problem with primal simplex */
    SCIP_CALL(performTest(TRUE, SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
        SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost));

    /* check that data stored in lpi is still the same */
    SCIP_CALL(checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));

    /* clear basis status */
    SCIP_CALL(SCIPlpiClearState(lpi));

    /* solve problem with dual simplex */
    SCIP_CALL(solveTest(FALSE, 2, 2, SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost));

    /* check that data stored in lpi is still the same */
    SCIP_CALL(checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));

    /* clear basis status */
    SCIP_CALL(SCIPlpiClearState(lpi));

    /* change objective */
    obj[0] = 1.0;
    SCIP_CALL(SCIPlpiChgObj(lpi, 1, ind, obj));

    /* change expected solution */
    exp_primsol[0] = 3;
    exp_primsol[1] = 4;
    exp_dualsol[0] = 0.4;
    exp_dualsol[1] = 0.2;
    exp_activity[0] = 10;
    exp_activity[1] = 15;
    exp_redcost[1] = 0;

    /* check changed problem with primal simplex */
    SCIP_CALL(performTest(TRUE, SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
        SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost));

    /* check that data stored in lpi is still the same */
    SCIP_CALL(checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));
    return SCIP_OKAY;
}

SCIP_RETCODE execmain_test2(int argc, const char** argv) {
    /* data with fixed values: */
    SCIP_Real obj[2] = { 3, 1 };
    SCIP_Real rhs[2] = { 10, 15 };
    int beg[2] = { 0, 2 };
    int ind[4] = { 0, 1, 0, 1 };
    SCIP_Real val[4] = { 2, 1, 1, 3 };

    /* data to be filled */
    SCIP_Real lb[2];
    SCIP_Real ub[2];
    SCIP_Real lhs[2];

    /* expected ray for first run */
    SCIP_Real exp_primray[2] = { 0.5, -1 };

    /* expected solutions for second run */
    SCIP_Real exp_primsol[2] = { 3, 4 };
    SCIP_Real exp_dualsol[2] = { 0.4, 0.2 };
    SCIP_Real exp_activity[2] = { 10, 15 };
    SCIP_Real exp_redcost[2] = { 0, 0 };

    /* fill variable data */
    lb[0] = -SCIPlpiInfinity(lpi);
    lb[1] = -SCIPlpiInfinity(lpi);
    ub[0] = SCIPlpiInfinity(lpi);
    ub[1] = SCIPlpiInfinity(lpi);
    lhs[0] = -SCIPlpiInfinity(lpi);
    lhs[1] = -SCIPlpiInfinity(lpi);

    /* solve problem with primal simplex */
    SCIP_CALL(performTest(TRUE, SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
        SCIPunbounded, SCIPinfeas, exp_primray, NULL, NULL, NULL));

    /* check that data stored in lpi is still the same */
    SCIP_CALL(checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));

    /* clear basis status */
    SCIP_CALL(SCIPlpiClearState(lpi));

    /* solve problem with dual simplex */
    SCIP_CALL(solveTest(FALSE, 2, 2, SCIPunbounded, SCIPinfeas, exp_primray, NULL, NULL, NULL));

    /* check that data stored in lpi is still the same */
    SCIP_CALL(checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));

    /* clear basis status */
    SCIP_CALL(SCIPlpiClearState(lpi));

    /* change objective */
    obj[0] = 1.0;
    SCIP_CALL(SCIPlpiChgObj(lpi, 1, ind, obj));

    /* solve with primal simplex */
    SCIP_CALL(solveTest(TRUE, 2, 2, SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost));

    /* check that data stored in lpi is still the same */
    SCIP_CALL(checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val));
    return SCIP_OKAY;
}

int main(int argc, const char* argv[]) {
    printf("Hello, SCIP! This problem would be solved by using SCIP integrated with SCS.\n");
    /**
    SCIP_CALL(SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE));
    if (execmain_test1(argc, argv) == SCIP_OKAY) {
        printf("The test 1 passed!\n");
    }
    SCIP_CALL(SCIPlpiFree(&lpi));
    assert(BMSgetMemoryUsed() == 0);
    */
    SCIP_CALL(SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE));
    if (execmain_test2(argc, argv) == SCIP_OKAY) {
        printf("The test 2 passed!\n");
    }
    SCIP_CALL(SCIPlpiFree(&lpi));
    assert(BMSgetMemoryUsed() == 0);
    return 0;
}
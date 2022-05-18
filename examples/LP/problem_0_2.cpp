//
// Created by vistart on 2022/1/25.
//

#include <iostream>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#define EPS 1e-6

/** macro to keep control of infinity
 *
 *  Some LPIs use std::numeric_limits<SCIP_Real>::infinity() as finity value. Comparing two infinty values then yields
 *  nan. This is a workaround. */
#define cr_assert_float_eq(Actual, Expected, Epsilon, FormatString, ...) \
      assert(ABS(Actual - Expected) < Epsilon);
#define cr_expect_float_eq(Actual, Expected, Epsilon, FormatString, ...) \
      if (ABS(Actual - Expected) < Epsilon) printf("Warning!\n");
#define cr_assert_float_eq_inf(Actual, Expected, Epsilon, FormatString, ...) \
   if ( ABS(Actual) < 1e30 && ABS(Expected) < 1e30 )                       \
      cr_assert_float_eq(Actual, Expected, Epsilon, FormatString, __VA_ARGS__);

SCIP_RETCODE execmain(int argc, const char** argv) {
    int ncols;
    int beg = 0;
    int nrows;
    int inds[3];
    SCIP_Real vals[2];
    SCIP_Real lb;
    SCIP_Real ub;
    SCIP_Real obj;
    SCIP_Real lhs;
    SCIP_Real rhs;
    SCIP_LPI* lpi;

    /* create LPI */
    SCIP_CALL(SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE));

    /* use the following LP:
     * max 1 x0 + 1 x1 + 1 x2
     *       -8 <= -x0 -          x2 <= -1
     *       -7 <= -x0 -   x1        <= -1
     *              x0 + 2 x1        <= 12
     *              x0,    x1,    x2 >= 0
     */
     /* add columns */
    lb = 0.0;
    ub = SCIPlpiInfinity(lpi);
    obj = 1.0;

    SCIP_CALL(SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL));
    SCIP_CALL(SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL));
    SCIP_CALL(SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL));

    /* add rows */
    lhs = -8.0;
    rhs = -1.0;
    inds[0] = 0;
    inds[1] = 2;
    vals[0] = -1.0;
    vals[1] = -1.0;
    SCIP_CALL(SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 2, &beg, inds, vals));

    lhs = -7.0;
    rhs = -1.0;
    inds[0] = 0;
    inds[1] = 1;
    vals[0] = -1.0;
    vals[1] = -1.0;
    SCIP_CALL(SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 2, &beg, inds, vals));

    lhs = -SCIPlpiInfinity(lpi);
    rhs = 12.0;
    inds[0] = 0;
    inds[1] = 1;
    vals[0] = 1.0;
    vals[1] = 2.0;
    SCIP_CALL(SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 2, &beg, inds, vals));

    /* check size */
    SCIP_CALL(SCIPlpiGetNRows(lpi, &nrows));
    SCIP_CALL(SCIPlpiGetNCols(lpi, &ncols));

    printf("nrows, ncols: %d, %d\n", nrows, ncols);

    SCIP_CALL(SCIPlpiSolvePrimal(lpi));

    SCIP_Real objval;
    SCIP_Real primsol[3];
    SCIP_CALL(SCIPlpiGetSol(lpi, &objval, primsol, NULL, NULL, NULL));
    printf("objval: %8.2f\n", objval);
    printf("primsol[0, 1, 2]: (%8.2f, %8.2f, %8.2f)\n", primsol[0], primsol[1], primsol[2]);
    /**
    SCIP_Real binvrow[3];
    SCIP_Real binvcol[3];
    SCIP_Real coef[3];
    SCIP_Real coeftwo[3];
    int cstats[3];
    int rstats[3];
    int basinds[3];
    int ninds;
    int idx;
    int entry;
    
    SCIP_CALL(SCIPlpiGetBase(lpi, cstats, rstats));
    printf("(%d, %d, %d), (%d, %d, %d)\n", cstats[0], cstats[1], cstats[2], rstats[0], rstats[1], rstats[2]);
    int bind[3];
    SCIP_CALL(SCIPlpiGetBasisInd(lpi, bind));
    printf("(%d, %d, %d)\n", bind[0], bind[1], bind[2]);

    for (int i = 0; i < 3; i++) {
        SCIP_Real coef[3];
        int inds[3];
        int ninds[3];
        SCIP_CALL(SCIPlpiGetBInvRow(lpi, i, coef, inds, ninds));
        printf("row[%d]: %f, %d, %d\n", i, coef[0], inds[0], ninds[0]);
    }*/

    SCIP_Real binvrow[3];
    SCIP_Real binvcol[3];
    SCIP_Real coef[3];
    SCIP_Real coeftwo[3];
    //SCIP_Real objval;
    int cstats[3];
    //int nrows;
    int rstats[3];
    int basinds[3];
    //int inds[3];
    int ninds;
    int idx;
    int entry;
    int i;

    /* expected values for the first column of BInv with corresponding variables */
    int exp_vars[] = { -2, 1, 2 };
    float exp_vals[] = { 0.0, 0.0, -1.0 };

    /* expected values for the first column of BAInv with corresponding variables */
    float exp_avals[] = { -0.5, 0.5, 1.0 };

    /* ------------------------------------- */
    /* first solve problem */
    SCIP_CALL(SCIPlpiSolvePrimal(lpi));

    SCIP_CALL(SCIPlpiGetObjval(lpi, &objval));
    cr_assert_float_eq(objval, 14.0, EPS, "objval: %f does not equal to 14.0\n", objval);

    /* the optimal basis should be: {x2, x3, slack for second row} */
    SCIP_CALL(SCIPlpiGetBase(lpi, cstats, rstats));
    printf("cstat[0, 1, 2]: (%d, %d, %d), rstat[0, 1, 2]: (%d, %d, %d)\n", cstats[0], cstats[1], cstats[2], rstats[0], rstats[1], rstats[2]);
    assert(cstats[0] == SCIP_BASESTAT_LOWER);
    assert(cstats[1] == SCIP_BASESTAT_BASIC);
    assert(cstats[2] == SCIP_BASESTAT_BASIC);

    assert(rstats[0] == SCIP_BASESTAT_LOWER);
    assert(rstats[1] == SCIP_BASESTAT_BASIC);
    assert(rstats[2] == SCIP_BASESTAT_UPPER);

    /* get basis indices */
    SCIP_CALL(SCIPlpiGetBasisInd(lpi, basinds));
    printf("basinds[0, 1, 2]: (%d, %d, %d)\n", basinds[0], basinds[1], basinds[2]);

    /* search for slack variable in basis */
    SCIP_CALL(SCIPlpiGetNRows(lpi, &nrows));
    for (i = 0; i < nrows; ++i)
    {
        if (basinds[i] < 0)
            break;
    }
    /* assert that we found the slack variable in the basis */
    assert(i < nrows);

    /* check basis inverse for the row corresponding to the basic slack variable */
    SCIP_CALL(SCIPlpiGetBInvRow(lpi, i, binvrow, NULL, NULL));
    for (int i = 0; i < 3; i++) {
        printf("binvrow[%d]: %f, %d, %d\n", i, binvrow[i]);
    }

    for (int i = 0; i < 3; i++) {
        SCIP_Real coef[3];
        int inds[3];
        int ninds[3];
        SCIP_CALL(SCIPlpiGetBInvRow(lpi, i, coef, inds, ninds));
        for (int j = 0; j < 3; j++) {
            printf("row[%d]: coef[%d]:%f, inds[%d]:%d, ninds[%d]:%d\n", i, j, coef[j], j, inds[j], j, ninds[j]);
        }
    }

    /* row of basis inverse should be (0, 1, 0.5) */
    cr_expect_float_eq(binvrow[0], 0.0, EPS, "BInvRow[%d] = %g != %g\n", 0, binvrow[0], 0.0);
    cr_expect_float_eq(binvrow[1], 1.0, EPS, "BInvRow[%d] = %g != %g\n", 1, binvrow[1], 1.0);
    cr_expect_float_eq(binvrow[2], 0.5, EPS, "BInvRow[%d] = %g != %g\n", 2, binvrow[2], 0.5);

    /* check whether sparse version is available and the same */
    SCIP_CALL(SCIPlpiGetBInvRow(lpi, i, coef, inds, &ninds));
    if (ninds >= 0)
    {
        assert(ninds == 2);
        for (entry = 0; entry < ninds; ++entry)
        {
            idx = inds[entry];
            assert(0 <= idx && idx < 3);
            cr_expect_float_eq(coef[idx], binvrow[idx], EPS, "coef[idx] = %g != %g\n", 0, coef[idx], binvrow[idx]);
        }
    }

    /* check first column of basis inverse */
    SCIP_CALL(SCIPlpiGetBInvCol(lpi, 0, binvcol, NULL, NULL));

    /* The columns will be in the same order, however, the rows might be permuted.
     * For each row/entry we check that it corresponds to the value of the corresponding variable.
     * The correspondance variable to row/entry is given by basinds. */
    for (entry = 0; entry < nrows; entry++)
    {
        /* for the given entry try each variable in exp_vars */
        for (idx = 0; idx < nrows; idx++)
        {
            /* Check that the value is the expected one if the column corresponds to the current variable given in exp_vars. */
            if (exp_vars[idx] == basinds[entry])
            {
                cr_expect_float_eq(binvcol[entry], exp_vals[idx], EPS, "binvcol[entry] = %g != %g\n", 0, binvcol[entry], exp_vals[idx]);
            }
        }
    }

    /* check whether number of nonzeros fits */
    SCIP_CALL(SCIPlpiGetBInvCol(lpi, 0, coef, inds, &ninds));
    if (ninds >= 0)
    {
        assert(ninds == 1);
    }

    /* check basis inverse times nonbasic matrix for row corresponding to the basic slack variable */
    assert(i >= 0);
    assert(i < nrows);
    SCIP_CALL(SCIPlpiGetBInvARow(lpi, i, NULL, coef, NULL, NULL));

    /* row of basis inverse times nonbasic matrix should be (-0.5, 0, 0) */
    cr_expect_float_eq(coef[0], -0.5, EPS, "BInvARow[%d] = %g != %g\n", 0, coef[0], -0.5);
    cr_expect_float_eq(coef[1], 0.0, EPS, "BInvARow[%d] = %g != %g\n", 1, coef[1], 0.0);
    cr_expect_float_eq(coef[2], 0.0, EPS, "BInvARow[%d] = %g != %g\n", 2, coef[2], 0.0);

    /* check nonzeros */
    SCIP_CALL(SCIPlpiGetBInvARow(lpi, i, NULL, coeftwo, inds, &ninds));
    if (ninds >= 0)
    {
        assert(ninds == 1);
        for (entry = 0; entry < ninds; ++entry)
        {
            idx = inds[entry];
            assert(0 <= idx && idx < 3);
            cr_expect_float_eq(coeftwo[idx], coef[idx], EPS, "coeftwo[idx] = %g != %g\n", coeftwo[idx], coef[idx]);
        }
    }

    /* check first column of basis inverse times nonbasic matrix */
    SCIP_CALL(SCIPlpiGetBInvACol(lpi, 0, coef, NULL, NULL));

    /* The columns will be in the same order, however, the rows will be permuted.
     * For each row/entry we check that it corresponds to the value of the corresponding variable.
     * The correspondance variable to row/entry is given by basinds. */
    for (entry = 0; entry < nrows; entry++)
    {
        /* for the given entry try each variable in exp_vars */
        for (idx = 0; idx < nrows; idx++)
        {
            /* Check that the value is the expected one if the column corresponds to the current variable given in exp_vars. */
            if (exp_vars[idx] == basinds[entry])
            {
                cr_expect_float_eq(coef[entry], exp_avals[idx], EPS, "coef[entry] = %g != %g\n", coef[entry], exp_avals[idx]);
            }
        }
    }

    /* check nonzeros */
    SCIP_CALL(SCIPlpiGetBInvACol(lpi, 0, coef, inds, &ninds));
    if (ninds >= 0)
    {
        assert(ninds == 3);
    }
    return SCIP_OKAY;
}

int main(int argc, const char* argv[]) {
    printf("Hello, SCIP! This problem would be solved by using SCIP integrated with SCS.\n");
    return execmain(argc, argv) != SCIP_OKAY ? 1 : 0;
}
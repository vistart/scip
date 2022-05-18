#include <iostream>
#include <scip/message_default.h>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#define EPS 1e-6

#define cr_assert_float_eq_inf(Actual, Expected, Epsilon) \
   if ( fabs(Actual) > 1e30 || fabs(Expected) > 1e30 )    \
      assert( Actual == Expected );                    \
   else                                                   \
      assert(ABS(Actual - Expected) < Epsilon);

/* GLOBAL VARIABLES */
static SCIP_LPI* lpi = NULL;
static SCIP_MESSAGEHDLR* messagehdlr = NULL;

SCIP_RETCODE setup() {
    /* this is necessary because the theories don't setup and teardown after each call, but once before and after */
    if (lpi != NULL)
    {
        SCIP_CALL(SCIPlpiFree(&lpi));
    }

    /* create message handler if necessary */
    if (messagehdlr == NULL)
    {
        /* create default message handler */
        SCIP_CALL(SCIPcreateMessagehdlrDefault(&messagehdlr, TRUE, NULL, FALSE));
    }

    /* This name is necessary because if CPLEX reads a problem from a file its problemname will be the filename. */
    SCIP_CALL(SCIPlpiCreate(&lpi, messagehdlr, "lpi_change_test_problem.lp", SCIP_OBJSEN_MAXIMIZE));
    return SCIP_OKAY;
}

SCIP_RETCODE teardown() {
    assert(!SCIPlpiWasSolved(lpi));
    SCIP_CALL(SCIPlpiFree(&lpi));

    if (messagehdlr != NULL)
    {
        SCIP_CALL(SCIPmessagehdlrRelease(&messagehdlr));
    }
    long long mu = BMSgetMemoryUsed();
    if (mu > 0) {
        printf("There is a memory leak! Actual %lld\n", mu);
        BMSdisplayMemory();
    }
    else {
        printf("There is no memory leak.\n");
    }
    return SCIP_OKAY;
}

SCIP_RETCODE execmain_test1(int argc, const char** argv) {
    /* problem data */
    SCIP_Real obj[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
    SCIP_Real  lb[5] = { -1.0, -SCIPlpiInfinity(lpi), 0.0, -SCIPlpiInfinity(lpi), 0.0 };
    SCIP_Real  ub[5] = { 10.0, SCIPlpiInfinity(lpi), SCIPlpiInfinity(lpi), 29.0, 0.0 };
    int ncolsbefore, ncolsafter;
    int nrowsbefore, nrowsafter;
    SCIP_Real lhsvals[6] = { -SCIPlpiInfinity(lpi), -1.0,   -3e-10, 0.0, 1.0,  3e10 };
    SCIP_Real rhsvals[6] = { -1.0,                  -3e-10, 0.0,    1.0, 3e10, SCIPlpiInfinity(lpi) };
    int     nnonzs[6] = { 1, 10, -1, 6, -1 };
    int    begvals[6] = { 0, 2, 3, 5, 8, 9 };
    int    indvals[10] = { 0, 1, 3, 2, 1, 1, 2, 4, 0, 3 };
    SCIP_Real vals[10] = { 1.0, 5.0, -1.0, 3e5, 2.0, 1.0, 20, 10, -1.9, 1e-2 };

    int iterations = 5;
    int k[5] = { 1, 6, -1, 4, -2 };
    int nnonzsdiff[5] = { 1, 10, -1, 6, -3 };
    int i;
    int j;

    /* create original lp */
    SCIP_CALL(SCIPlpiAddCols(lpi, 5, obj, lb, ub, NULL, 0, NULL, NULL, NULL)); SCIPlpiGetInternalStatus(lpi);
    SCIP_CALL(SCIPlpiGetNCols(lpi, &ncolsbefore));
    for (i = 0; i < iterations; i++)
    {
        int nrows;
        int nnonzsbefore;
        int nnonzsafter;

        nrows = k[i];

        SCIP_CALL(SCIPlpiGetNNonz(lpi, &nnonzsbefore));
        SCIP_CALL(SCIPlpiGetNRows(lpi, &nrowsbefore));

        if (nrows < 0)
        {
            SCIP_CALL(SCIPlpiDelRows(lpi, 0, -(1 + nrows)));
        }
        else
        {
            SCIP_Real lhs[100];
            SCIP_Real rhs[100];
            int beg[100];

            int nnonz = nnonzs[i];
            int ind[100];
            SCIP_Real val[100];

            SCIP_Real newlhs[100];
            SCIP_Real newval[100];
            SCIP_Real newrhs[100];
            int newbeg[100];
            int newind[100];
            int newnnonz;
            int indold;
            int indnew;

            assert(nrows < 100);
            for (j = 0; j < nrows; j++)
            {
                lhs[j] = lhsvals[j];
                rhs[j] = rhsvals[j];
                beg[j] = begvals[j];
            }

            assert(nnonz < 100);
            for (j = 0; j < nnonz; j++)
            {
                ind[j] = indvals[j];
                val[j] = vals[j];
            }
            SCIP_CALL(SCIPlpiAddRows(lpi, nrows, lhs, rhs, NULL, nnonz, beg, ind, val)); SCIPlpiGetInternalStatus(lpi);

            SCIP_CALL(SCIPlpiGetRows(lpi, nrowsbefore, nrowsbefore - 1 + nrows, newlhs, newrhs, &newnnonz, newbeg, newind, newval));
            assert(nnonz == newnnonz);
            SCIPdebugMessage("`beg` and `newbeg` should be equal.\n");
            for (int i = 0; i < nrows; i++) {
                SCIPdebugMessage("row[%d]: (beg, newbeg), (%d, %d)\n", i, beg[i], newbeg[i]);
                assert(beg[i] == newbeg[i]);
            }

            beg[nrows] = nnonz;
            newbeg[nrows] = newnnonz;

            for (j = 0; j < nrows; j++)
            {
                cr_assert_float_eq_inf(lhs[j], newlhs[j], 1e-16);
                cr_assert_float_eq_inf(rhs[j], newrhs[j], 1e-16);

                for (indold = beg[j]; indold < beg[j + 1]; indold++)
                {
                    int occurrences = 0;

                    for (indnew = beg[j]; indnew < beg[j + 1]; indnew++)
                    {
                        if (ind[indold] == newind[indnew])
                        {
                            occurrences = occurrences + 1;
                            assert(ABS(val[indold] - newval[indnew]) < 1e-16);
                        }
                    }
                    assert(occurrences == 1);
                }
            }
        }

        SCIP_CALL(SCIPlpiGetNRows(lpi, &nrowsafter));
        assert(nrowsbefore + nrows == nrowsafter);

        SCIP_CALL(SCIPlpiGetNNonz(lpi, &nnonzsafter));
        SCIPdebugMessage("nnonzsbefore %d, nnonzsafter %d, nnonzsdiff[i] %d, in iteration %d\n", nnonzsbefore, nnonzsafter, nnonzsdiff[i], i);
        assert(nnonzsbefore + nnonzsdiff[i] == nnonzsafter);

        SCIP_CALL(SCIPlpiGetNCols(lpi, &ncolsafter));
        assert(ncolsbefore == ncolsafter);
    }

    SCIP_CALL(SCIPlpiGetNRows(lpi, &nrowsbefore));
    assert(8 == nrowsbefore);
    for (i = 3; i > 0; i--)
    {
        int rows[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

        for (j = 0; j < i; j++)
            rows[(2 * j) + 1] = 1;

        SCIP_CALL(SCIPlpiGetNRows(lpi, &nrowsbefore));
        SCIP_CALL(SCIPlpiDelRowset(lpi, rows)); SCIPlpiGetInternalStatus(lpi);
        SCIP_CALL(SCIPlpiGetNRows(lpi, &nrowsafter));
        assert(nrowsbefore - i == nrowsafter);
    }
    return SCIP_OKAY;
}

SCIP_RETCODE execmain_test2(int argc, const char** argv) {
    /* problem data */
    SCIP_Real obj[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
    SCIP_Real lhs[5] = { -1.0, -SCIPlpiInfinity(lpi), 0.0, -SCIPlpiInfinity(lpi), 0.0 };
    SCIP_Real rhs[5] = { 10.0, SCIPlpiInfinity(lpi), SCIPlpiInfinity(lpi), 29.0, 0.0 };
    int ncolsbefore, ncolsafter;
    int nrowsbefore, nrowsafter;
    SCIP_Real lbvals[6] = { -SCIPlpiInfinity(lpi), -1.0, -3e-10, 0.0, 1.0, 3e10 };
    SCIP_Real ubvals[6] = { -1.0, -3e-10, 0.0, 1.0, 3e10, SCIPlpiInfinity(lpi) };
    SCIP_Real   vals[10] = { 1.0, 5.0, -1.0, 3e5, 2.0, 1.0, 20, 10, -1.9, 1e-2 };
    int  nnonzs[6] = { 1, 10, -1, 6, -1 };
    int begvals[6] = { 0, 2, 3, 5, 8, 9 };
    int indvals[10] = { 0, 1, 3, 2, 1, 1, 2, 4, 0, 3 };

    int iterations = 5;
    int k[5] = { 1, 6, -1, 4, -2 };
    int nnonzsdiff[5] = { 1, 10, -1, 6, -3 };
    int i;
    int j;

    /* create original lp */
    SCIP_CALL(SCIPlpiAddRows(lpi, 5, lhs, rhs, NULL, 0, NULL, NULL, NULL));
    SCIP_CALL(SCIPlpiGetNRows(lpi, &nrowsbefore));

    for (i = 0; i < iterations; i++)
    {
        BMSdisplayMemory();
        /* setup col values */
        int ncols;
        int nnonzsbefore;
        int nnonzsafter;

        ncols = k[i];

        /* get data before modification */
        SCIP_CALL(SCIPlpiGetNNonz(lpi, &nnonzsbefore));
        SCIP_CALL(SCIPlpiGetNCols(lpi, &ncolsbefore));

        if (ncols < 0)
        {
            SCIP_CALL(SCIPlpiDelCols(lpi, 0, -(1 + ncols)));
            BMSdisplayMemory();
        }
        else
        {  /* ncols >= 0 */
            SCIP_Real lb[100];
            SCIP_Real ub[100];
            int beg[100];

            int nnonz = nnonzs[i];
            int ind[100];
            SCIP_Real val[100];

            SCIP_Real newlb[100];
            SCIP_Real newval[100];
            SCIP_Real newub[100];
            int newbeg[100];
            int newind[100];
            int newnnonz;

            assert(ncols < 100);
            for (j = 0; j < ncols; j++)
            {
                lb[j] = lbvals[j];
                ub[j] = ubvals[j];
                beg[j] = begvals[j];
            }

            assert(nnonz < 100);
            for (j = 0; j < nnonz; j++)
            {
                ind[j] = indvals[j];
                val[j] = vals[j];
            }
            //printf("ncols, j: %d, %d\n", ncols, j);
            SCIP_CALL(SCIPlpiAddCols(lpi, ncols, obj, lb, ub, NULL, nnonz, beg, ind, val));
            BMSdisplayMemory();

            /* checks */
            SCIP_CALL(SCIPlpiGetCols(lpi, ncolsbefore, ncolsbefore - 1 + ncols, newlb, newub, &newnnonz, newbeg, newind, newval));
            assert(nnonz == newnnonz);
            //cr_assert_eq(nnonz, newnnonz, "expecting %d, got %d\n", nnonz, newnnonz);
            for (int i = 0; i < ncols; i++) {
                assert(lb[i] == newlb[i]);
                assert(ub[i] == newub[i]);
                assert(beg[i] == newbeg[i]);
            }
            for (int i = 0; i < nnonz; i++) {
                assert(ind[i] == newind[i]);
                assert(val[i] == newval[i]);
            }
            /**
            cr_assert_arr_eq(lb, newlb, ncols * sizeof(SCIP_Real));
            cr_assert_arr_eq(ub, newub, ncols * sizeof(SCIP_Real));
            cr_assert_arr_eq(beg, newbeg, ncols * sizeof(int));
            cr_assert_arr_eq(ind, newind, nnonz * sizeof(int));
            cr_assert_arr_eq(val, newval, nnonz * sizeof(SCIP_Real));*/
        }

        /* checks */
        SCIP_CALL(SCIPlpiGetNRows(lpi, &nrowsafter));
        assert(nrowsbefore == nrowsafter);
        //cr_assert_eq(nrowsbefore, nrowsafter);

        SCIP_CALL(SCIPlpiGetNNonz(lpi, &nnonzsafter));
        assert(nnonzsbefore + nnonzsdiff[i] == nnonzsafter);
        //cr_assert_eq(nnonzsbefore + nnonzsdiff[i], nnonzsafter, "nnonzsbefore %d, nnonzsafter %d, nnonzsdiff[i] %d, in iteration %d\n",
        //    nnonzsbefore, nnonzsafter, nnonzsdiff[i], i);

        SCIP_CALL(SCIPlpiGetNCols(lpi, &ncolsafter));
        assert(ncolsbefore + ncols == ncolsafter);
        //cr_assert_eq(ncolsbefore + ncols, ncolsafter);
    }

    /* delete rowsets */
    /* should have 8 rows now */
    SCIP_CALL(SCIPlpiGetNCols(lpi, &ncolsbefore));
    assert(ncolsbefore == 8);
    //cr_assert_eq(8, ncolsbefore);
    for (i = 3; i > 0; i--)
    {
        int cols[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

        for (j = 0; j < i; j++)
            cols[(2 * j) + 1] = 1;

        SCIP_CALL(SCIPlpiGetNCols(lpi, &ncolsbefore));
        SCIP_CALL(SCIPlpiDelColset(lpi, cols));
        SCIP_CALL(SCIPlpiGetNCols(lpi, &ncolsafter));
        assert(ncolsbefore - i == ncolsafter);
        //cr_assert_eq(ncolsbefore - i, ncolsafter);
        /* assert that the rows that are left are the ones I intended */
    }
    return SCIP_OKAY;
}

SCIP_RETCODE execmain_test3(int argc, const char** argv) {
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

    //lpi = NULL;

    /* create LPI */
    //SCIP_CALL(SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE));

    /* use the following LP as base:
     *   max x
     *       1 <= x <= 2  (linear constraint)
     *       0 <= x <= 3  (bounds)
     */
     /* add one column */
    SCIP_CALL(SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL));

    /* add one row */
    SCIP_CALL(SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 1, &beg, &ind, &val));

    /* check size */
    SCIP_CALL(SCIPlpiGetNRows(lpi, &nrows));
    SCIP_CALL(SCIPlpiGetNCols(lpi, &ncols));
    assert(nrows == 1);
    assert(ncols == 1);

#ifdef SCIP_DEBUG
    /* turn on output */
    SCIP_CALL(SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPINFO, 1));
#endif
    // unittests-lpi-bases, Test(simple, test1)
    int cstat;
    int rstat;

    /* solve problem */
    SCIP_CALL(SCIPlpiSolvePrimal(lpi));

    /* get basis */
    SCIP_CALL(SCIPlpiGetBase(lpi, &cstat, &rstat));

    /* the variable should be basic and the slack variable at the upper bound */
    assert(cstat == SCIP_BASESTAT_BASIC);
    assert(rstat == SCIP_BASESTAT_UPPER);

    SCIP_CALL(SCIPlpiFree(&lpi));

    if (messagehdlr != NULL)
    {
        SCIP_CALL(SCIPmessagehdlrRelease(&messagehdlr));
    }
    long long mu = BMSgetMemoryUsed();
    if (mu > 0) {
        printf("There is a memory leak! Actual %lld\n", mu);
        BMSdisplayMemory();
    }
    else {
        printf("There is no memory leak.\n");
    }

    return SCIP_OKAY;
}

SCIP_RETCODE execmain_test4(int argc, const char** argv) {
    int ncols;
    int beg = 0;
    int inds[2];
    int nrows;
    SCIP_Real vals[2];
    SCIP_Real lb;
    SCIP_Real ub;
    SCIP_Real obj;
    SCIP_Real lhs;
    SCIP_Real rhs;

    /* create LPI */
    //SCIP_CALL(SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE));

    /* use the following LP:
     * max 1 x1 + 1 x2 + 1 x3
     *       -8 <= -x1 -          x3 <= -1
     *       -7 <= -x1 -   x2        <= -1
     *              x1 + 2 x2        <= 12
     *              x1,    x2,    x3 >= 0
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
    assert(nrows == 3);
    assert(ncols == 3);

#ifdef SCIP_DEBUG
    /* turn on output */
    SCIP_CALL(SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPINFO, 1));
#endif
    return SCIP_OKAY;
}

int main(int argc, const char* argv[]) {
    printf("Hello, SCIP! This problem would be solved by using SCIP integrated with SCS.\n");
    if (!setup()) {
        printf("The test 1 failed to initialize.\n");
    }
    if (execmain_test1(argc, argv) == SCIP_OKAY) {
        printf("The test 1 passed!\n");
    }
    teardown();
    if (!setup()) {
        printf("The test 2 failed to initialize.\n");
    }
    if (execmain_test2(argc, argv) == SCIP_OKAY) {
        printf("The test 2 passed!\n");
    }
    teardown();
    if (!setup()) {
        printf("The test 3 failed to initialize.\n");
    }
    if (execmain_test3(argc, argv) == SCIP_OKAY) {
        printf("The test 3 passed!\n");
    }
    // test3 does not require `teardown()`.
    if (!setup()) {
        printf("The test 4 failed to initialize.\n");
    }
    if (execmain_test4(argc, argv) == SCIP_OKAY) {
        printf("The test 4 passed!\n");
    }
    teardown();
    return 0;
}
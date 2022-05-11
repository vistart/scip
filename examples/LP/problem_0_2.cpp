//
// Created by vistart on 2022/1/25.
//

#include <iostream>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

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
    }
    return SCIP_OKAY;
}

int main(int argc, const char* argv[]) {
    printf("Hello, SCIP! This problem would be solved by using SCIP integrated with SCS.\n");
    return execmain(argc, argv) != SCIP_OKAY ? 1 : 0;
}
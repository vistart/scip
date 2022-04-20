//
// Created by vistart on 2022/1/25.
//

#include <iostream>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

SCIP_RETCODE execmain(int argc, const char** argv) {
    int nrows;
    int ncols;
    int beg = 0;
    SCIP_Real* lb = (SCIP_Real*)calloc(sizeof(SCIP_Real), 2);
    lb[0] = 0;
    lb[1] = 0;
    SCIP_Real* ub = (SCIP_Real*)calloc(sizeof(SCIP_Real), 2);
    ub[0] = 2;
    ub[1] = 2;
    SCIP_Real lhs = 0.0;
    SCIP_Real rhs = 4.0;
    SCIP_Real* obj = (SCIP_Real*)calloc(sizeof(SCIP_Real), 2);
    obj[0] = 1;
    obj[1] = 1;
    SCIP_Real* val = (SCIP_Real*)calloc(sizeof(SCIP_Real), 2);
    val[0] = 1;
    val[1] = 2;
    int* ind = (int*)calloc(sizeof(int), 2);
    ind[0] = 0;
    ind[1] = 1;
    SCIP_LPI* lpi;

    SCIP_CALL(SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE));

    /* use the following LP as base:
     *   max x
     *       1 <= x <= 2  (linear constraint)
     *       0 <= x <= 3  (bounds)
     */
     /* add one column */
    SCIP_CALL(SCIPlpiAddCols(lpi, 2, obj, lb, ub, NULL, 0, NULL, NULL, NULL));

    /* add one row */
    SCIP_CALL(SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 2, &beg, ind, val));

    SCIP_CALL(SCIPlpiSolvePrimal(lpi));

    SCIP_CALL(SCIPlpiGetSol(lpi, NULL, NULL, NULL, NULL, NULL));

    return SCIP_OKAY;
}

int main(int argc, const char * argv[]) {
    printf("Hello, SCIP! This problem would be solved by using SCIP integrated with SCS.\n");
    return execmain(argc, argv) != SCIP_OKAY ? 1 : 0;
}
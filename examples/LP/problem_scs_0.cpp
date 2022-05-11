//
// Created by vistart on 2022/1/25.
//

#include <iostream>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <vector>
#include <map>

#include "problem_scs_random.h"

SCIP_RETCODE execmain(int argc, const char** argv) {
    /**
    ScsCone* k = (ScsCone*)scs_calloc(1, sizeof(ScsCone));
    ScsData* d = (ScsData*)scs_calloc(1, sizeof(ScsData));
    ScsSettings* stgs = (ScsSettings*)scs_calloc(1, sizeof(ScsSettings));
    ScsSolution* sol = (ScsSolution*)scs_calloc(1, sizeof(ScsSolution));
    ScsSolution* opt_sol = (ScsSolution*)scs_calloc(1, sizeof(ScsSolution));
    ScsInfo info = { 0 };
    const scs_float p_f = 0.1;
    int seed = 1234;
    const scs_int n = 4000;
    const scs_int m = 8000;
    const scs_int col_nnz = (scs_int)ceil(sqrt(n));
    const scs_int nnz = n * col_nnz;
    scs_int exitflag;
    scs_float perr, derr;
    scs_int success;
    const char* fail;

    k->z = m; // (scs_int)floor(m * p_f);
    k->l = m - k->z;

    d->m = m;
    d->n = n;
    gen_random_prob_data(nnz, col_nnz, d, k, opt_sol, seed);
    print_d(d, nnz);
    */
    SCIP_LPI* lpi;
    SCIP_CALL(SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MINIMIZE));
    SCIP_CALL(SCIPlpiFree(&lpi));
    long long mu = BMSgetMemoryUsed();
    printf("The actual remaining: %lld\n", mu);
    return SCIP_OKAY;
}

int main(int argc, const char* argv[]) {
    printf("Hello, SCIP! This problem would be solved by using SCIP integrated with SCS.\n");
    return execmain(argc, argv) != SCIP_OKAY ? 1 : 0;
}
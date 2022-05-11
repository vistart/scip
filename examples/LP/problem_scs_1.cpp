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

    SCIP_LPI* lpi;
    SCIP_CALL(SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MINIMIZE));
    double lb = -SCIPlpiInfinity(lpi);
    double ub = SCIPlpiInfinity(lpi);
    for (int i = 0; i < n; i++) {
        SCIP_CALL(SCIPlpiAddCols(lpi, 1, &d->c[i], &lb, &ub, NULL, 0, NULL, NULL, NULL));
    }
    
    std::map<int, std::vector<std::pair<int, double>>> cons;
    for (int i = 0; i < d->A->p[d->n]; i += col_nnz) {
        int v = i / col_nnz;
        for (int j = i; j < i + col_nnz; j++) {
            auto it = cons.find(d->A->i[j]);
            if (it == cons.end()) {
                cons.insert({ d->A->i[j], {{ v, d->A->x[j] }} });
                continue;
            }
            it->second.push_back({v, d->A->x[j]});
        }
    }
    for (auto i = cons.begin(); i != cons.end(); i++) {
        int beg = 0;
        int ind[i->second.size()];
        double val[i->second.size()];
        int p = 0;
        std::cout << i->first << ":";
        for (auto j : i->second) {
            std::cout << "(" << j.first << "," << j.second << ") ";
            ind[p] = j.first;
            val[p] = j.second;
            p++;
        }
        std::cout << std::endl;
        SCIP_CALL(SCIPlpiAddRows(lpi, 1, &lb, &d->b[i->first], NULL, i->second.size(), &beg, ind, val));
    }
    /**
    for (int i = 0; i < d->A->p[d->n]; i++) {
        auto it = cons.find(d->A->i[i]);
        if (it != cons.end()) {
            //std::cout << it->first << " size:" << it->second.size() << std::endl;
            it->second.push_back(d->A->x[i]);
        }
        else {
            cons.insert({ d->A->i[i], {d->A->x[i]} });
        }
    }
    /**
    std::cout << cons.size() << std::endl;
    for (auto it = cons.begin(); it != cons.end(); it++) {
        std::cout << it->first << ":";
        for (auto i : it->second) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
    /**
    for (int i = 0; i < n; i += col_nnz) {
        int beg = 0;
        int ind[col_nnz];
        double val[col_nnz];
        for (int j = 0; j < col_nnz; j++) {
            ind[j] = d->A->i[i + j];
            val[j] = d->A->x[i + j];
        }
        SCIP_CALL(SCIPlpiAddRows(lpi, 1, &lb, &d->b[i], NULL, col_nnz, &beg, ind, val));
    }
    */
    SCIP_CALL(SCIPlpiSolvePrimal(lpi));
    double objval[1];
    double primsol[n];
    double dualsol[m];
    SCIP_CALL(SCIPlpiGetSol(lpi, objval, primsol, dualsol, NULL, NULL));
    scs_printf("Objective: %8.4f\n", *objval);
    print_sol_prim(primsol, n);
    print_sol_dual(dualsol, m);
    /**
    scs_set_default_settings(stgs);
    stgs->eps_abs = 1e-5;
    stgs->eps_rel = 1e-5;

    exitflag = scs(d, k, stgs, sol, &info);

    perr = SCS(dot)(d->c, sol->x, d->n) - SCS(dot)(d->c, opt_sol->x, d->n);
    derr = -SCS(dot)(d->b, sol->y, d->m) + SCS(dot)(d->b, opt_sol->y, d->m);
    scs_printf("true obj %4e\n", SCS(dot)(d->c, opt_sol->x, d->n));
    scs_printf("primal obj error %4e\n", perr);
    scs_printf("dual obj error %4e\n", derr);

    success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

    //mu_assert("small_lp: SCS failed to produce outputflag SCS_SOLVED", success);
    fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
    SCS(free_data)(d);
    SCS(free_cone)(k);
    SCS(free_sol)(sol);
    SCS(free_sol)(opt_sol);
    scs_free(stgs);*/
    return SCIP_OKAY;
}

int main(int argc, const char* argv[]) {
    printf("Hello, SCIP! This problem would be solved by using SCIP integrated with SCS.\n");
    return execmain(argc, argv) != SCIP_OKAY ? 1 : 0;
}
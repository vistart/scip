//
// Created by vistart on 2022/2/18.
//

#include <iostream>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#include "scip/set.h"

SCIP_RETCODE execmain(int argc, const char** argv) {
    SCIP* scip = nullptr;
    SCIP_CALL(SCIPcreate(&scip));
    SCIP_CALL(SCIPsetBoolParam(scip, "lp/checkdualfeas", FALSE));
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));
    SCIP_CALL(SCIPcreateProbBasic(scip, "SCIP_scs_example"));
    SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE));

    SCIP_VAR* x1 = nullptr;
    SCIP_CALL(SCIPcreateVarBasic(scip, &x1, "x1", 0.0, 40.0, 1.0, SCIP_VARTYPE_CONTINUOUS));
    SCIP_CALL(SCIPaddVar(scip, x1));

    SCIP_VAR* x2 = nullptr;
    SCIP_CALL(SCIPcreateVarBasic(scip, &x2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), 2.0, SCIP_VARTYPE_CONTINUOUS));
    SCIP_CALL(SCIPaddVar(scip, x2));

    SCIP_VAR* x3 = nullptr;
    SCIP_CALL(SCIPcreateVarBasic(scip, &x3, "x3", -SCIPinfinity(scip), SCIPinfinity(scip), 3.0, SCIP_VARTYPE_CONTINUOUS));
    SCIP_CALL(SCIPaddVar(scip, x3));

    SCIP_VAR* x4 = nullptr;
    SCIP_CALL(SCIPcreateVarBasic(scip, &x4, "x4", 2.0, 3.0, 1.0, SCIP_VARTYPE_INTEGER));
    SCIP_CALL(SCIPaddVar(scip, x4));

    SCIP_CONS* cons1 = nullptr;
    SCIP_CALL(SCIPcreateConsBasicLinear(scip, &cons1, "cons1", 0, nullptr, nullptr, -SCIPinfinity(scip), 20.0));
    SCIP_CALL(SCIPaddCoefLinear(scip, cons1, x1, -1.0));
    SCIP_CALL(SCIPaddCoefLinear(scip, cons1, x2, 1.0));
    SCIP_CALL(SCIPaddCoefLinear(scip, cons1, x3, 1.0));
    SCIP_CALL(SCIPaddCoefLinear(scip, cons1, x4, 10.0));
    SCIP_CALL( SCIPaddCons(scip, cons1));

    SCIP_CONS* cons2 = nullptr;
    SCIP_CALL(SCIPcreateConsBasicLinear(scip, &cons2, "cons2", 0, nullptr, nullptr, -SCIPinfinity(scip), 30));
    SCIP_CALL(SCIPaddCoefLinear(scip, cons2, x1, 1.0));
    SCIP_CALL(SCIPaddCoefLinear(scip, cons2, x2, -3.0));
    SCIP_CALL(SCIPaddCoefLinear(scip, cons2, x3, 1.0));
    SCIP_CALL( SCIPaddCons(scip, cons2));

    SCIP_CONS* cons3 = nullptr;;
    SCIP_CALL(SCIPcreateConsBasicLinear(scip,&cons3,"cons3",0,nullptr,nullptr,0,0.0));
    SCIP_CALL( SCIPaddCoefLinear(scip, cons3, x2, 1.0));
    SCIP_CALL( SCIPaddCoefLinear(scip, cons3, x4, -3.5));
    SCIP_CALL( SCIPaddCons(scip, cons3));

    //Scip releasing all the constraints
    SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
    SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
    SCIP_CALL( SCIPreleaseCons(scip, &cons3) );

    SCIP_CALL( SCIPsolve(scip) );

    SCIP_SOL* sol;
    sol = SCIPgetBestSol(scip);
    std::cout << "The solution(s):" << std::endl;
    std::cout << "x1: " << SCIPgetSolVal(scip, sol, x1) << " "
              << "x2: " << SCIPgetSolVal(scip, sol, x2) << " "
              << "x3: " << SCIPgetSolVal(scip, sol, x3) << " "
              << "x4: " << SCIPgetSolVal(scip, sol, x4) << std::endl;
    SCIP_CALL( (SCIPwriteOrigProblem(scip, "problem_2_example.lp", nullptr, FALSE)));
    //Freeing the variables
    SCIP_CALL( SCIPreleaseVar(scip, &x1));
    SCIP_CALL( SCIPreleaseVar(scip, &x2));
    SCIP_CALL( SCIPreleaseVar(scip, &x3));
    SCIP_CALL( SCIPreleaseVar(scip, &x4));
    SCIP_CALL( SCIPfree(&scip) );
    return SCIP_OKAY;
}

int main(int argc, const char * argv[]) {
    printf("Hello, SCIP! This problem would be solved by using SCIP integrated with SCS.\n");
    return execmain(argc, argv) != SCIP_OKAY ? 1 : 0;
}
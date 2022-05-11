//
// Created by vistart on 2022/1/25.
//

#include <iostream>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

SCIP_RETCODE execmain(int argc, const char** argv) {
    SCIP* scip = nullptr;
    SCIP_CALL( SCIPcreate(&scip));
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));
    SCIP_CALL(SCIPcreateProbBasic(scip, "SCIP_scs_example"));
    SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE));

    SCIP_VAR* x1 = nullptr;
    SCIP_CALL(SCIPcreateVarBasic(scip, &x1, "x1", 0.0, SCIPinfinity(scip), 3.0, SCIP_VARTYPE_INTEGER));
    SCIP_CALL(SCIPaddVar(scip, x1));

    SCIP_VAR* x2 = nullptr;
    SCIP_CALL(SCIPcreateVarBasic(scip, &x2, "x2", 0.0, SCIPinfinity(scip), 3.0, SCIP_VARTYPE_INTEGER));
    SCIP_CALL(SCIPaddVar(scip, x2));

    SCIP_CONS* cons1 = nullptr;
    SCIP_CALL(SCIPcreateConsBasicLinear(scip, &cons1, "cons1", 0, nullptr, nullptr, -SCIPinfinity(scip), 4));
    SCIP_CALL(SCIPaddCoefLinear(scip, cons1, x1, 1.0));
    SCIP_CALL(SCIPaddCoefLinear(scip, cons1, x2, 1.0));
    SCIP_CALL( SCIPaddCons(scip, cons1));

    SCIP_CONS* cons2 = nullptr;
    SCIP_CALL(SCIPcreateConsBasicLinear(scip, &cons2, "cons2", 0, nullptr, nullptr, -SCIPinfinity(scip), 5));
    SCIP_CALL(SCIPaddCoefLinear(scip, cons2, x1, 2.0));
    SCIP_CALL(SCIPaddCoefLinear(scip, cons2, x2, 1.0));
    SCIP_CALL( SCIPaddCons(scip, cons2));

    SCIP_CONS* cons3 = nullptr;
    SCIP_CALL(SCIPcreateConsBasicLinear(scip,&cons3,"cons3",0,nullptr,nullptr,-SCIPinfinity(scip),-2.0));
    SCIP_CALL( SCIPaddCoefLinear(scip, cons3, x1, 1.0));
    SCIP_CALL( SCIPaddCoefLinear(scip, cons3, x2, -4.0));
    SCIP_CALL( SCIPaddCons(scip, cons3));

    //Scip releasing all the constraints
    SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
    SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
    SCIP_CALL( SCIPreleaseCons(scip, &cons3) );

    SCIP_CALL( SCIPsolve(scip) );

    SCIP_SOL* sol;
    sol = SCIPgetBestSol(scip);
    std::cout << "x1: " << SCIPgetSolVal(scip, sol, x1) << " " << "x2: " << SCIPgetSolVal(scip, sol, x2) << std::endl;
    SCIP_CALL( (SCIPwriteOrigProblem(scip, "problem_1_example.lp", nullptr, FALSE)));
    //Freeing the variables
    SCIP_CALL( SCIPreleaseVar(scip, &x1));
    SCIP_CALL( SCIPreleaseVar(scip, &x2));
    SCIP_CALL( SCIPfree(&scip) );
    return SCIP_OKAY;
}

int main(int argc, const char * argv[]) {
    printf("Hello, SCIP! This problem would be solved by using SCIP integrated with SCS.\n");
    return execmain(argc, argv) != SCIP_OKAY ? 1 : 0;
}
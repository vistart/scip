/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch_conffullstrong.h,v 1.2 2004/02/04 17:27:16 bzfpfend Exp $"

/**@file   branch_conffullstrong.h
 * @brief  full strong LP branching rule, that creates infeasible children to give input to conflict analysis
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __BRANCH_CONFFULLSTRONG_H__
#define __BRANCH_CONFFULLSTRONG_H__


#include "scip.h"


/** creates the conflict full strong LP braching rule and includes it in SCIP */
extern
RETCODE SCIPincludeBranchruleConffullstrong(
   SCIP*            scip                /**< SCIP data structure */
   );

#endif

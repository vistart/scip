/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_xyz.h
 * @ingroup NLHDLRS
 * @brief  xyz nonlinear handler
 * @author Jane Doe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLHDLR_XYZ_H__
#define __SCIP_NLHDLR_XYZ_H__

#include "scip/scip.h"
#include "scip/pub_nlhdlr.h"

#ifdef __cplusplus
extern "C" {
#endif


/**@addtogroup NLHDLRS
 *
 * @{
 *
 * @name Xyz nonlinear handler.
 *
 * @{
 */

/** @}
  * @}
  */

/** includes xyz nonlinear handler to nonlinear constraint handler
 *
 * @ingroup NlhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrXyz(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLHDLR_XZY_H__ */

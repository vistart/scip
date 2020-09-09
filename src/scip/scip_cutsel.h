/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_cutsel.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for cut selector plugins
 * @author Felipe Serrano
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_CUTSEL_H__
#define __SCIP_SCIP_CUTSEL_H__


#include "scip/def.h"
#include "scip/type_cutsel.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicCutSelectorMethods
 *
 * @{
 */

/** creates a cut selector and includes it in SCIP.
 *
 *  @note method has all cut selector callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludeCutselBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeCutsel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of cut selector */
   const char*           desc,               /**< description of cut selector */
   int                   priority,           /**< priority of the cut selector */
   SCIP_DECL_CUTSELCOPY ((*cutselcopy)),     /**< copy method of cut selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CUTSELFREE ((*cutselfree)),     /**< destructor of cut selector */
   SCIP_DECL_CUTSELINIT ((*cutselinit)),     /**< initialize cut selector */
   SCIP_DECL_CUTSELEXIT ((*cutselexit)),     /**< deinitialize cut selector */
   SCIP_DECL_CUTSELINITSOL((*cutselinitsol)),/**< solving process initialization method of cut selector */
   SCIP_DECL_CUTSELEXITSOL((*cutselexitsol)),/**< solving process deinitialization method of cut selector */
   SCIP_DECL_CUTSELSELECT((*cutselselect)),  /**< cut selection method */
   SCIP_CUTSELDATA*     cutseldata           /**< cut selector data */
   );

/** Creates a cut selector and includes it in SCIP with its most fundamental callbacks. All non-fundamental
 *  (or optional) callbacks as, e.g., init and exit callbacks, will be set to NULL.
 *  Optional callbacks can be set via specific setter functions, see SCIPsetCutselCopy(), SCIPsetCutselFree(),
 *  SCIPsetCutselInit(), SCIPsetCutselExit(), SCIPsetCutselInitsol(), and SCIPsetCutselExitsol()
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeCutsel() instead
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeCutselBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTSEL**         cutsel,             /**< reference to a cut selector, or NULL */
   const char*           name,               /**< name of cut selector */
   const char*           desc,               /**< description of cut selector */
   int                   priority,           /**< priority of the cut selector in standard mode */
   SCIP_DECL_CUTSELSELECT((*cutselselect)),  /**< cut selection method */
   SCIP_CUTSELDATA*      cutseldata          /**< cut selector data */
   );

/** sets copy method of cut selector */
SCIP_EXPORT
SCIP_RETCODE SCIPsetCutselCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELCOPY  ((*cutselcopy))  /**< copy method of cut selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of cut selector */
SCIP_RETCODE SCIPsetCutselFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELFREE  ((*cutselfree))    /**< destructor of cut selector */
   );

/** sets initialization method of cut selector */
SCIP_RETCODE SCIPsetCutselInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELINIT ((*cutselinit))     /**< initialize node selector */
   );

/** sets deinitialization method of cut selector */
SCIP_RETCODE SCIPsetCutselExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELEXIT  ((*cutselexit))     /**< deinitialize cut selector */
   );

/** sets solving process initialization method of cut selector */
SCIP_RETCODE SCIPsetCutselInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELINITSOL  ((*cutselinitsol))/**< solving process initialization method of cut selector */
   );

/** sets solving process deinitialization method of cut selector */
SCIP_RETCODE SCIPsetCutselExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTSEL*          cutsel,             /**< cut selector */
   SCIP_DECL_CUTSELEXITSOL  ((*cutselexitsol))/**< solving process deinitialization method of cut selector */
   );

/** returns the cut selector of the given name, or NULL if not existing */
SCIP_EXPORT
SCIP_CUTSEL* SCIPfindCutsel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of cut selector */
   );

/** returns the array of currently available cut selectors */
SCIP_CUTSEL** SCIPgetCutsels(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available cut selectors */
int SCIPgetNCutsels(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** @} */

#ifdef __cplusplus
}
#endif

#endif

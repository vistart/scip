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
#pragma ident "@(#) $Id: pub_sepa.h,v 1.3 2004/02/04 17:27:38 bzfpfend Exp $"

/**@file   pub_sepa.h
 * @brief  public methods for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_SEPA_H__
#define __PUB_SEPA_H__


#include "def.h"
#include "type_misc.h"
#include "type_sepa.h"



/** compares two separators w. r. to their priority */
extern
DECL_SORTPTRCOMP(SCIPsepaComp);

/** gets user data of separator */
extern
SEPADATA* SCIPsepaGetData(
   SEPA*            sepa                /**< separator */
   );

/** sets user data of separator; user has to free old data in advance! */
extern
void SCIPsepaSetData(
   SEPA*            sepa,               /**< separator */
   SEPADATA*        sepadata            /**< new separator user data */
   );

/** gets name of separator */
extern
const char* SCIPsepaGetName(
   SEPA*            sepa                /**< separator */
   );

/** gets description of separator */
extern
const char* SCIPsepaGetDesc(
   SEPA*            sepa                /**< separator */
   );

/** gets priority of separator */
extern
int SCIPsepaGetPriority(
   SEPA*            sepa                /**< separator */
   );

/** gets frequency of separator */
extern
int SCIPsepaGetFreq(
   SEPA*            sepa                /**< separator */
   );

/** gets time in seconds used in this separator */
extern
Real SCIPsepaGetTime(
   SEPA*            sepa                /**< separator */
   );

/** gets the total number of times, the separator was called */
extern
Longint SCIPsepaGetNCalls(
   SEPA*            sepa                /**< separator */
   );

/** gets the number of times, the separator was called at the current node */
extern
int SCIPsepaGetNCallsAtNode(
   SEPA*            sepa                /**< separator */
   );

/** gets the total number of cutting planes found by this separator */
extern
Longint SCIPsepaGetNCutsFound(
   SEPA*            sepa                /**< separator */
   );

/** gets the number of cutting planes found by this separator at the current node */
extern
Longint SCIPsepaGetNCutsFoundAtNode(
   SEPA*            sepa                /**< separator */
   );

/** is separator initialized? */
extern
Bool SCIPsepaIsInitialized(
   SEPA*            sepa                /**< separator */
   );


#endif

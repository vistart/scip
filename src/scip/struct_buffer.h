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
#pragma ident "@(#) $Id: struct_buffer.h,v 1.2 2004/02/04 17:27:43 bzfpfend Exp $"

/**@file   struct_buffer.h
 * @brief  datastructures for memory buffers for temporary objects
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_BUFFER_H__
#define __STRUCT_BUFFER_H__


#include <assert.h>

#include "def.h"



/** memory buffer storage for temporary objects */
struct Buffer
{
   void**           data;               /**< allocated memory chunks for arbitrary data */
   int*             size;               /**< sizes of buffers in bytes */
   Bool*            used;               /**< TRUE iff corresponding buffer is in use */
   int              ndata;              /**< number of memory chunks */
   int              firstfree;          /**< first unused memory chunk */
};


#endif

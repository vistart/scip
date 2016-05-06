/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dirreduce.c
 * @brief  several basic reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements basic redution techniques for several Steiner problems.
 * All tests are described in "A Generic Approach to Solving the Steiner Tree Problem and Variants" by Daniel Rehfeldt.
 *
 * A list of all interface methods can be found in grph.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "grph.h"
#include "portab.h"
#include "scip/scip.h"


/** is there no vertex of higher prize? */
static
SCIP_Bool maxprize(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int                   i                   /**< the terminal to be checked */
   )
{
   int k;
   int t = -1;
   SCIP_Real max = -1.0;

   for( k = 0; k < g->knots; k++ )
   {
      if( Is_term(g->term[k]) && g->mark[k] && g->grad[k] > 0 )
      {
	 assert(k != g->source[0]);
	 if( SCIPisGT(scip, g->prize[k], max) )
	 {
            max = g->prize[k];
            t = k;
	 }
	 else if( t == i && SCIPisGE(scip, g->prize[k], max) )
	 {
            t = k;
	 }
      }
   }

   SCIPdebugMessage("maxprize: %f (from %d) \n", g->prize[t], t );
   return (t == i);
}

/** try to eliminate a terminal of degree one */
static
SCIP_RETCODE trydg1edgepc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            offset,             /**< pointer to store the offset */
   int*                  count,              /**< pointer storing number of eliminated edges */
   int                   i,                  /**< the terminal to be checked */
   int                   iout,               /**< outgoing arc */
   SCIP_Bool*            rerun               /**< further eliminations possible? */
   )
{
   int i1;
   int degsum;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(count  != NULL);
   assert(Is_term(g->term[i]));

   if( maxprize(scip, g, i) )
      return SCIP_OKAY;

   i1 = g->head[iout];

   if( SCIPisLE(scip, g->prize[i], g->cost[iout]) && g->stp_type != STP_MAX_NODE_WEIGHT )
   {
      /* delete terminal */

      if( (i1 < i) && (Is_term(g->term[i1]) || g->grad[i1] == 2 || g->grad[i1] == 3) )
         (*rerun) = TRUE;
      SCIPdebugMessage("Delete (degree 1) terminal %d \n", i);
      (*offset) += g->prize[i];
      *count += deleteterm(scip, g, i);
   }
   else
   {
      /* contract terminal */

      (*rerun) = TRUE;
      assert(SCIPisGT(scip, g->prize[i], 0.0 ));

      if( g->stp_type == STP_MAX_NODE_WEIGHT )
      {
	 if( SCIPisLT(scip, g->prize[i], -g->prize[i1]) )
            *offset += g->prize[i];
	 else
            *offset -= g->prize[i1];
      }
      else
      {
         *offset += g->cost[iout];
      }

      if( g->source[0] == i1 )
      {
         if( g->pcancestors[i] != NULL )
         {
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[i1]), g->pcancestors[i]) );
            SCIPintListNodeFree(scip, &(g->pcancestors[i]));
         }
	 SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[i1]), g->ancestors[iout]) );
         *count += deleteterm(scip, g, i);
         return SCIP_OKAY;
      }

      degsum = g->grad[i] + g->grad[i1];

      SCIP_CALL( graph_knot_contractpc(scip, g, i, i1, i) );

      degsum -= g->grad[i];

      assert(degsum >= 1);

      if( g->stp_type == STP_MAX_NODE_WEIGHT )
      {
	 int e;
	 int t = UNKNOWN;
	 int e2 = UNKNOWN;

	 if( SCIPisLT(scip, g->prize[i], 0.0 ) )
	 {
            i1 = UNKNOWN;
            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
            {
               i1 = g->head[e];
               if( Is_pterm(g->term[i1]) && g->source[0] != i1 )
                  t = i1;
	       else if( g->source[0] == i1 )
                  e2 = e;
            }

            assert(t != UNKNOWN);
            assert(e2 != UNKNOWN);

            /* delete artifical terminal */
            graph_knot_chg(g, t, -1);
            while( g->outbeg[t] != EAT_LAST )
            {
	       e = g->outbeg[t];
	       g->cost[e] = 0.0;
	       g->cost[flipedge(e)] = 0.0;
               graph_edge_del(scip, g, e, TRUE);
               count++;
            }

            assert(g->grad[t] == 0);

	    /* i is not a terminal anymore */
	    graph_knot_chg(g, i, -1);
	    graph_edge_del(scip, g, e2, TRUE);

            for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
               if( g->mark[g->tail[e]] )
                  g->cost[e] = -g->prize[i];

            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
	    {
               i1 = g->head[e];
               if( g->mark[i1]  )
               {
                  if( !Is_term(g->term[i1]) )
                  {
                     g->cost[e] = -g->prize[i1];
                  }
                  else
                  {
                     g->cost[e] = 0.0;
                  }
               }
	    }
	 }
         else
         {
            for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
               if( g->mark[g->tail[e]] )
                  g->cost[e] = 0.0;

	    for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
	    {
               i1 = g->head[e];
               if( g->mark[i1]  )
               {
                  if( !Is_term(g->term[i1]) )
                  {
                     assert(SCIPisLE(scip, g->prize[i1], 0.0 ));
                     g->cost[e] = -g->prize[i1];
                  }
                  else
                  {
                     assert(SCIPisGE(scip, g->prize[i1], 0.0 ));
                     g->cost[e] = 0.0;
                  }
               }
               else if( Is_pterm(g->term[i1]) && g->source[0] != i1 )
	       {
		  t = i1;
	       }

	    }
	    assert(t != UNKNOWN);

	    for( e = g->inpbeg[t]; e != EAT_LAST; e = g->ieat[e] )
               if( g->tail[e] == g->source[0] )
                  break;
            assert(e != EAT_LAST);
	    g->cost[e] = g->prize[i];
	 }
      }
      *count += degsum;
   }
   return SCIP_OKAY;
}


/** traverse one side of a chain (MWCSP) */
static
SCIP_RETCODE traverseChain(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  length,             /**< pointer to store length of chain */
   int*                  final,              /**< pointer to store final vertex */
   int                   i,                  /**< start vertex */
   int                   i1,                 /**< first vertex */
   int                   i2,                 /**< last vertex */
   int                   e1                  /**< first edge */
   )
{
   IDX* ancestors = NULL;
   IDX* revancestors = NULL;
   SCIP_Real sum;
   int k;
   int e;

   assert(g != NULL);
   assert(scip != NULL);
   assert(length != NULL);

   k = i1;
   e = e1;
   sum = 0.0;

   while( g->grad[k] == 2 && !Is_term(g->term[k]) && k != i2 )
   {
      assert(g->mark[k]);

      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors), g->ancestors[e]) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors), g->ancestors[flipedge(e)]) );

      if( e != e1 )
         graph_edge_del(scip, g, e, TRUE);

      e = g->outbeg[k];
      sum += g->prize[k];
      (*length)++;

      if( e == flipedge(e1) )
         e = g->oeat[e];

      assert(e != EAT_LAST);
      assert(SCIPisLE(scip, g->prize[k], 0.0));

      k = g->head[e];
   }
   if( k != i1 )
   {
      int ne;

      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors), g->ancestors[e]) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors), g->ancestors[flipedge(e)]) );

      graph_edge_del(scip, g, e, TRUE);

      g->prize[i] += sum;
      ne = graph_edge_redirect(scip, g, e1, i, k, 1.0);

      if( ne != -1 )
      {
         e1 = ne;

         SCIPintListNodeFree(scip, &(g->ancestors[e1]));
         SCIPintListNodeFree(scip, &(g->ancestors[flipedge(e1)]));

         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[e1]), ancestors) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[flipedge(e1)]), revancestors) );

      }
      else
      {
	 for( e1 = g->outbeg[i]; e1 != EAT_LAST; e1 = g->oeat[e1] )
            if( g->head[e1] == k )
               break;
         assert(e1 != EAT_LAST);
      }

      SCIPintListNodeFree(scip, &(ancestors));
      SCIPintListNodeFree(scip, &(revancestors));

      if( SCIPisGE(scip, g->prize[k], 0.0) )
         g->cost[e1] = 0.0;
      else
         g->cost[e1] = -g->prize[k];
      assert(SCIPisLE(scip, g->prize[i], 0.0) );
   }

   *final = k;

   return SCIP_OKAY;
}


/** delete a terminal for a (rooted) prize-collecting problem */
int deleteterm(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int                   i                   /**< index of the terminal */
   )
{
   int e;
   int t;
   int i1;
   int count;

   assert(g != NULL);
   assert(scip != NULL);
   assert(Is_term(g->term[i]));

   t = UNKNOWN;
   count = g->grad[i] + 2;

   /* delete terminal */

   graph_knot_chg(g, i, -1);
   g->mark[i] = FALSE;

   while( (e = g->outbeg[i]) != EAT_LAST )
   {
      i1 = g->head[e];

      if( Is_pterm(g->term[i1]) && g->source[0] != i1 )
         t = g->head[e];
      graph_edge_del(scip, g, e, TRUE);
   }

   assert(t != UNKNOWN);

   /* delete artifical terminal */

   graph_knot_chg(g, t, -1);

   while( g->outbeg[t] != EAT_LAST )
      graph_edge_del(scip, g, g->outbeg[t], TRUE);


   return count;
}


/** basic reduction tests for the STP */
SCIP_RETCODE degree_test(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  nelims               /**< pointer to number of reductions */
   )
{
   int i;
   int i1;
   int i2;
   int e1;
   int e2;
   int nnodes;
   int rerun = TRUE;
   int done  = TRUE;
   int count = 0;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(nelims != NULL);

   nnodes = g->knots;

   SCIPdebugMessage("Degree Test: ");

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      for( i = 0; i < nnodes; i++ )
      {
         assert(g->grad[i] >= 0);

         if( g->grad[i] == 1 )
         {
            e1  = g->outbeg[i];
            i1  = g->head[e1];

            assert(e1 >= 0);
            assert(e1 == Edge_anti(g->inpbeg[i]));
            assert(g->oeat[e1] == EAT_LAST);
            assert(g->ieat[g->inpbeg[i]] == EAT_LAST);

            if( Is_term(g->term[i]) )
            {
               *fixed += g->cost[e1];
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1]) );
            }

            SCIP_CALL( graph_knot_contract(scip, g, i1, i) );
            count++;

            assert(g->grad[i] == 0);

            /* the last node in the graph? */
            if( g->grad[i1] == 0 )
            {
               rerun = FALSE;
               break;
            }
            if( (i1 < i) && (g->grad[i1] < 3) )
               rerun = TRUE;

            continue;
         }
         if( g->grad[i] == 2 )
         {
            e1 = g->outbeg[i];
            e2 = g->oeat[e1];
            i1 = g->head[e1];
            i2 = g->head[e2];

            assert(e1 >= 0);
            assert(e2 >= 0);

            do
            {
               done = TRUE;

               if( !Is_term(g->term[i]) )
               {
                  assert(EQ(g->cost[e2], g->cost[Edge_anti(e2)]));

                  g->cost[e1]            += g->cost[e2];
                  g->cost[Edge_anti(e1)] += g->cost[e2];
                  SCIP_CALL( graph_knot_contract(scip, g, i2, i) );
                  count++;

                  break;
               }
               assert(Is_term(g->term[i]));

               if( Is_term(g->term[i1]) && Is_term(g->term[i2]) )
               {

                  if( SCIPisLT(scip, g->cost[e1], g->cost[e2]) )
                  {
                     *fixed += g->cost[e1];
                     SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1]) );
                     SCIP_CALL( graph_knot_contract(scip, g, i1, i) );
                  }
                  else
                  {
                     *fixed += g->cost[e2];
                     SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e2]) );
                     SCIP_CALL( graph_knot_contract(scip, g, i2, i) );
                  }
                  count++;

                  break;
               }
               if( Is_term(g->term[i1]) && !Is_term(g->term[i2]) && SCIPisLE(scip, g->cost[e1], g->cost[e2]) )
               {
                  *fixed += g->cost[e1];
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1]) );
                  SCIP_CALL( graph_knot_contract(scip, g, i1, i) );
                  count++;
                  break;
               }
               if( Is_term(g->term[i2]) && !Is_term(g->term[i1]) && SCIPisLE(scip, g->cost[e2], g->cost[e1]) )
               {
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e2]) );
                  *fixed += g->cost[e2];
                  SCIP_CALL( graph_knot_contract(scip, g, i2, i) );
                  count++;
                  break;
               }
               done = FALSE;
            }
            while( FALSE );

            if (done
               && (((i1 < i) && (g->grad[i1] < 3))
                  || ((i2 < i) && (g->grad[i2] < 3))))
               rerun = TRUE;
         }
         if( Is_term(g->term[i]) && g->grad[i] > 2 )
         {
            SCIP_Real mincost = FARAWAY;
            int ett = UNKNOWN;
            for( e1 = g->outbeg[i]; e1 != EAT_LAST; e1 = g->oeat[e1] )
            {
               i1 = g->head[e1];

               if( SCIPisLT(scip, g->cost[e1], mincost) )
               {
                  mincost = g->cost[e1];
                  if( Is_term(g->term[i1]) )
                     ett = e1;
               }
               else if( Is_term(g->term[i1]) && SCIPisLE(scip, g->cost[e1], mincost) )
               {
                  ett = e1;
               }
            }
            if( ett != UNKNOWN && SCIPisLE(scip, g->cost[ett], mincost) )
            {
               *fixed += g->cost[ett];
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[ett]) );
               SCIP_CALL( graph_knot_contract(scip, g, i, g->head[ett]) );
               rerun = TRUE;
            }
         }
      }
   }
   SCIPdebugMessage(" %d Knots deleted\n", count);
   assert(graph_valid(g));

   *nelims += count;
   return SCIP_OKAY;
}



/** basic reduction tests for the SAP */
SCIP_RETCODE degree_test_sap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  count               /**< pointer to number of reductions */
   )
{
   SCIP_QUEUE* queue;
   int i;
   int e;
   int i1;
   int i2;
   int e1;
   int e2;
   int root;
   int nnodes;
   int* pnode;
   char rerun;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count != NULL);

   root = g->source[0];
   rerun = TRUE;
   nnodes = g->knots;

   *count = 0;
   SCIPdebugMessage("Degree Test: ");

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      for( i = 0; i < nnodes; i++ )
      {
         assert(g->grad[i] >= 0);

         if( g->grad[i] == 1 )
         {
            e1  = g->inpbeg[i];
            i1  = g->tail[e1];

            assert(e1 >= 0);
            assert(e1 == Edge_anti(g->outbeg[i]));
            assert(g->ieat[e1] == EAT_LAST);
            assert(g->oeat[g->outbeg[i]] == EAT_LAST);

            if( Is_term(g->term[i]) )
            {
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1]) );
               *fixed += g->cost[e1];
               SCIP_CALL( graph_knot_contract(scip, g, i1, i) );
            }
            else
            {
               graph_edge_del(scip, g, e1, TRUE);
            }

            assert(g->grad[i] == 0);

            if ((i1 < i) && (g->grad[i1] < 3))
               rerun = TRUE;

            (*count)++;

            continue;
         }

         if( g->grad[i] == 2 )
         {
            e1 = g->outbeg[i];
            e2 = g->oeat[e1];
            i1 = g->head[e1];
            i2 = g->head[e2];

            assert(e1 >= 0);
            assert(e2 >= 0);

            if( !Is_term(g->term[i]) )
            {
               if( (!Is_term(g->term[i2]) && !Is_term(g->term[i1])) )
               {
                  g->cost[e1] += g->cost[Edge_anti(e2)];
                  g->cost[Edge_anti(e1)] += g->cost[e2];
		  if( SCIPisGT(scip, g->cost[e1], FARAWAY) )
                     g->cost[e1] = FARAWAY;
		  if( SCIPisGT(scip, g->cost[Edge_anti(e1)], FARAWAY) )
                     g->cost[Edge_anti(e1)] = FARAWAY;
                  SCIP_CALL( graph_knot_contract(scip, g, i2, i) );
                  (*count)++;
                  if( ((i1 < i) && (g->grad[i1] < 3))
                     || ((i2 < i) && (g->grad[i2] < 3)) )
                     rerun = TRUE;
               }
            }
            /* CONSTCOND */
            /*lint -save -e717 */
            /*lint -restore */
         }
      }
   }

   /* delete all arcs in \delta^-(root) */
   for( e = g->inpbeg[root]; e != EAT_LAST; e = g->ieat[e] )
      g->cost[e] = FARAWAY;

   /* delete all arcs in not connected to a terminal other than the root by forward arcs */

   /* BFS until all terminals are reached */
   SCIP_CALL( SCIPqueueCreate(&queue, nnodes, 2.0) );

   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) && i != root )
      {
         g->mark[i] = TRUE;
	 SCIP_CALL( SCIPqueueInsert(queue, &(g->tail[g->outbeg[i]])) );
      }
      else
      {
	 g->mark[i] = FALSE;
      }
   }

   g->mark[root] = TRUE;

   while( !SCIPqueueIsEmpty(queue) )
   {
      pnode = (SCIPqueueRemove(queue));
      for( e = g->inpbeg[*pnode]; e != EAT_LAST; e = g->ieat[e] )
      {
         if( !g->mark[g->tail[e]] )
         {
            g->mark[g->tail[e]] = TRUE;
            SCIP_CALL( SCIPqueueInsert(queue, &(g->tail[e])) );
         }
      }
   }

   SCIPqueueFree(&queue);

   for( i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] )
      {
	 while( g->inpbeg[i] != EAT_LAST )
	 {
	    SCIPdebugMessage("remove edge to node %d \n", i);
            graph_edge_del(scip, g, g->inpbeg[i], TRUE);
	 }
      }
#if 0
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
	 if( SCIPisGE(scip, g->cost[e], FARAWAY) &&  SCIPisGE(scip, g->cost[flipedge(e)], FARAWAY) )
	 {
	    printf("remove high cost edge to node %d \n", i);
            graph_edge_del(scip, g, e, TRUE);
	 }
      }
#endif
   }

   SCIPdebugMessage("dirdeg %d Knots deleted\n", *count);
   assert(graph_valid(g));

   return SCIP_OKAY;
}


/** root proximity terminal test (SAP) */
SCIP_RETCODE rptReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offset value */
   int*                  count               /**< pointer to number of reductions */
   )
{
   SCIP_Real pathcost;
   SCIP_Real* dijkdist;
   int i;
   int e;
   int i1;
   int e1;
   int old;
   int root;
   int nnodes;
   int* dijkedge;

   assert(scip != NULL);
   assert(g != NULL);
   assert(fixed != NULL);
   assert(count != NULL);

   root = g->source[0];
   nnodes = g->knots;
   *count = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &dijkdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dijkedge, nnodes) );

   graph_path_execX(scip, g, root, g->cost, dijkdist, dijkedge);

   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) && i != root && g->grad[i] > 0 )
      {
	 e1 = dijkedge[i];
	 pathcost = dijkdist[i];

	 for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
	 {
	    if( e == e1 )
               continue;

	    if( SCIPisGT(scip, pathcost, g->cost[e]) )
               break;
	 }
	 if( e == EAT_LAST )
	 {
            i1 = g->tail[e1];
            old = g->grad[i] + g->grad[i1] - 1;

            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1]) );
            *fixed += g->cost[e1];
            SCIP_CALL( graph_knot_contract(scip, g, i1, i) );

            assert(old - g->grad[i1] > 0);
            *count += old - g->grad[i1];
            SCIPdebugMessage("contract %d\n", old - g->grad[i] - g->grad[i1]);
	 }

      }
   }

   SCIPfreeBufferArray(scip, &dijkedge);
   SCIPfreeBufferArray(scip, &dijkdist);

   return SCIP_OKAY;
}


/** basic reduction tests for the MWCS problem */
SCIP_RETCODE degree_test_mw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  count               /**< pointer to number of reductions */
   )
{
   int i;
   int e;
   int i1;
   int i2;
   int e1;
   int e2;
   int nnodes;
   int nedges;
   SCIP_Bool rerun = TRUE;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count != NULL);
   assert(g->stp_type == STP_MAX_NODE_WEIGHT);

   nnodes = g->knots;
   nedges = g->edges;

   SCIPdebugMessage("MW degree test: \n");

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      /* contract adjacent positive vertices */
      for( e = 0; e < nedges; e += 2 )
      {
	 i1 = g->tail[e];
	 i2 = g->head[e];
	 if( g->mark[i1] && g->mark[i2] && Is_term(g->term[i1]) && Is_term(g->term[i2]) )
	 {
            SCIPdebugMessage("contract tt %d->%d\n ", i1, i2);
            (*count)++;
            SCIP_CALL( graph_knot_contractpc(scip, g, i1, i2, i1) );
	 }
      }

      /* main loop for remaining tests */
      for( i = 0; i < nnodes; i++ )
      {
         assert(g->grad[i] >= 0);
         if( !g->mark[i] || g->grad[i] == 0 )
            continue;

         assert( !SCIPisEQ(scip, g->prize[i], 0.0) );

	 /* non-positive vertex? */
         if( !Is_term(g->term[i]) )
         {
            if( g->grad[i] == 1 )
            {
               e1 = g->inpbeg[i];
               i1 = g->tail[e1];
               assert(e1 >= 0);
               assert(e1 == Edge_anti(g->outbeg[i]));
               assert(g->ieat[e1] == EAT_LAST);
               assert(g->oeat[g->outbeg[i]] == EAT_LAST);
	       assert(SCIPisLE(scip, g->prize[i], 0.0));

               graph_edge_del(scip, g, e1, TRUE);
               SCIPdebugMessage("delete negative vertex of degree 1 (%d)\n ",  i);
               assert(g->grad[i] == 0);

               if( (i1 < i) && (g->grad[i1] < 3 || (g->grad[i1] == 3 && Is_term(g->term[i1]))) )
                  rerun = TRUE;

               (*count)++;
               continue;
            }

            /* contract non-positive chains */
            if( g->grad[i] == 2 )
            {
               int f1 = -1;
               int f2 = -1;
               int length = 0;

               e1 = g->outbeg[i];
               e2 = g->oeat[e1];
               i1 = g->head[e1];
               i2 = g->head[e2];

               assert(e1 >= 0);
               assert(e2 >= 0);
	       assert(i1 != i2);
               assert(g->mark[i1]);
               assert(g->mark[i2]);

	       SCIP_CALL( traverseChain(scip, g, &length, &f1, i, i1, i2, e1) );
	       SCIP_CALL( traverseChain(scip, g, &length, &f2, i, i2, i1, e2) );

	       if( f1 == f2 )
	       {
		  while( g->outbeg[i] != EAT_LAST )
		     graph_edge_del(scip, g, g->outbeg[i], TRUE);
	       }
	       else if( length > 0 )
	       {
                  assert(g->grad[i] <= 2);

                  for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
                     g->cost[e] = -g->prize[i];

                  e1 = g->outbeg[i];
                  e2 = g->oeat[e1];

                  (*count) += length;
	       }
	    }
            continue;
         }

         /* node i is of positive weight (terminal): */

         /* terminal of (real) degree 0? */
         if( g->grad[i] == 2 )
         {
            /* if terminal node i is not the one with the highest prize, delete */
            if( !maxprize(scip, g, i) )
            {
               SCIPdebugMessage("delete degree 0 term %d prize: %f count:%d\n ", i, g->prize[i], *count);
               (*fixed) += g->prize[i];
               (*count) += deleteterm(scip, g, i);
            }
         }
         /* terminal of (real) degree 1? */
         else if( g->grad[i] == 3 )
         {
            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
               if( g->mark[g->head[e]] )
                  break;
            assert(e != EAT_LAST);
            assert(g->head[e] != g->source[0]);
            if( !Is_term(g->term[g->head[e]]) )
            {
               SCIP_CALL( trydg1edgepc(scip, g, fixed, count, i, e, &rerun) );
               continue;
            }
         }
      }
   }

   SCIPdebugMessage("MW basic reduction package has deleted %d edges\n", *count);

   return SCIP_OKAY;
}


/** basic reduction tests for the HCDSTP */
SCIP_RETCODE degree_test_hc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  count               /**< pointer to number of reductions */
   )
{
   int i;
   int e;
   int e2;
   int root;
   int nnodes;
   SCIP_Bool rerun = TRUE;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count  != NULL);
   assert(g->stp_type == STP_HOP_CONS);

   nnodes = g->knots;
   root = g->source[0];

   SCIPdebugMessage("basic HC test: \n");

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      /* delete incoming arcs of the root */
      e = g->inpbeg[root];
      while( e != EAT_LAST )
      {
         e2 = g->ieat[e];

         if( SCIPisGE(scip, g->cost[flipedge(e)], FARAWAY) )
         {
            SCIPdebugMessage("delete incoming root arc \n");
            (*count)++;
            graph_edge_del(scip, g, e, TRUE);
         }
         else if( SCIPisLT(scip, g->cost[e], FARAWAY) )
         {
            SCIPdebugMessage("delete anti-parallel root arcs \n");
            g->cost[e] = FARAWAY;
         }

         e = e2;
      }

      /* delete outgoing arcs of the terminals (other than the root) */
      for( i = 0; i < nnodes; i++ )
      {
         if( Is_term(g->term[i]) && i != root )
         {
            e = g->outbeg[i];
            while( e != EAT_LAST )
            {
               e2 = g->oeat[e];

               if( SCIPisGE(scip, g->cost[flipedge(e)], FARAWAY) )
               {
                  SCIPdebugMessage("delete anti-parallel terminal arcs \n");
                  (*count)++;
                  graph_edge_del(scip, g, e, TRUE);
               }

               e = e2;
            }
         }
      }
   }

   SCIPdebugMessage("HC basic reduction package has deleted %d edges\n", *count);

   return SCIP_OKAY;
}


/** basic reductions for RPCSTP and PCSPG */
SCIP_RETCODE degree_test_pc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  count               /**< pointer to number of reductions */
   )
{
   int* edges2;
   int* nodes2;
   int i;
   int i1;
   int i2;
   int e;
   int e1;
   int e2;
   int nnodes;
   SCIP_Bool pc;
   SCIP_Bool rerun = TRUE;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count != NULL);
   assert(g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_ROOTED_PRIZE_COLLECTING);

   pc = (g->stp_type == STP_PRIZE_COLLECTING);

   nnodes = g->knots;
   *count = 0;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &edges2, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes2, 2) );

   SCIPdebugMessage("Degree Test: ");

   if( !pc )
      g->mark[g->source[0]] = FALSE;

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      for( i = 0; i < nnodes; i++ )
      {
         assert(g->grad[i] >= 0);
         if( !g->mark[i] || g->grad[i] == 0 )
            continue;

         if( !Is_term(g->term[i]) )
         {
            if( g->grad[i] == 1 )
            {
               e1 = g->inpbeg[i];
               i1 = g->tail[e1];
               assert(e1 >= 0);
               assert(e1 == Edge_anti(g->outbeg[i]));
               assert(g->ieat[e1] == EAT_LAST);
               assert(g->oeat[g->outbeg[i]] == EAT_LAST);

               graph_edge_del(scip, g, e1, TRUE);
               SCIPdebugMessage("delete NT %d\n ",  i);
               assert(g->grad[i] == 0);

               /* the last node? */
               if( g->grad[i1] == 0 )
               {
                  rerun = FALSE;
                  break;
               }
               if( (i1 < i) && (g->grad[i1] < 3 || Is_term(g->term[i1])) )
                  rerun = TRUE;

               (*count)++;
               continue;
            }

            /* contract non terminals of degree 2 */
            if( g->grad[i] == 2 )
            {
               e1 = g->outbeg[i];
               e2 = g->oeat[e1];
               i1 = g->head[e1];
               i2 = g->head[e2];

               assert(e1 >= 0);
               assert(e2 >= 0);
               assert(g->mark[i1] || i1 == g->source[0]);
               assert(g->mark[i2] || i2 == g->source[0]);
               assert(EQ(g->cost[e2], g->cost[Edge_anti(e2)]));

               g->cost[e1]            += g->cost[e2];
               g->cost[Edge_anti(e1)] += g->cost[e2];

	       SCIPdebugMessage("contract NT %d %d \n ", i2, i);
               SCIP_CALL( graph_knot_contract(scip, g, i2, i) );
               (*count)++;

               if( (Is_term(g->term[i2]) && (i2 < i)) || (Is_term(g->term[i1]) && (i1 < i)) )
                  rerun = TRUE;
            }
            continue;
         }

         /* node i is a terminal: */

         /* terminal of (real) degree 0? */
         if( ( (g->grad[i] == 2 && pc) || (g->grad[i] == 1 && !pc) ) )
         {
            /* if terminal node i is node the one with the highest prize, delete*/
            if( !maxprize(scip, g, i) )
            {
               SCIPdebugMessage("delete 0 term %d prize: %f count:%d\n ", i, g->prize[i], *count);
               (*fixed) += g->prize[i];
               (*count) += deleteterm(scip, g, i);
            }
         }
         /* terminal of (real) degree 1? */
         else if( ( (g->grad[i] == 3 && pc) || (g->grad[i] == 2 && !pc) ) )
         {
            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
               if( g->mark[g->head[e]] || (!pc && g->head[e] == g->source[0]) )
                  break;
            assert(e != EAT_LAST);
            assert(g->head[e] != g->source[0] || !pc);

            SCIP_CALL( trydg1edgepc(scip, g, fixed, count, i, e, &rerun) );
         }
         /* terminal of (real) degree 2? */
         else if( ( (g->grad[i] == 4 && pc) || (g->grad[i] == 3 && !pc)  )  )
         {
            if( !maxprize(scip, g, i) )
            {
               i2 = 0;
               for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
               {
                  i1 = g->head[e];
                  if( g->mark[i1] )
                  {
                     if( i2 >= 2 )
                        assert(i2 < 2);
                     edges2[i2] = e;
                     nodes2[i2++] = i1;
                  }
               }
               if( SCIPisLE(scip, g->prize[i], g->cost[edges2[0]]) && SCIPisLE(scip, g->prize[i], g->cost[edges2[1]]) )
               {
                  int n1;
                  IDX* ancestors = NULL;
                  IDX* revancestors = NULL;

                  e = edges2[0];
                  e1 = edges2[1];
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors), g->ancestors[e]) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors), g->ancestors[Edge_anti(e1)]) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors), g->ancestors[Edge_anti(e)]) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors), g->ancestors[e1]) );
                  SCIPdebugMessage("delete - term - %d\n ", i);

                  /* contract edge */
                  n1 = graph_edge_redirect(scip, g, e, nodes2[1], nodes2[0], g->cost[e] + g->cost[e1] - g->prize[i]);

                  /* new edge inserted? */
                  if( n1 >= 0)
                  {
                     /* add ancestors */
                     SCIPintListNodeFree(scip, &(g->ancestors[n1]));
                     SCIPintListNodeFree(scip, &(g->ancestors[Edge_anti(n1)]));
                     SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[n1]), ancestors) );
                     SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), revancestors) );
                  }
                  (*count) += deleteterm(scip, g, i);
                  (*fixed) += g->prize[i];
                  SCIPintListNodeFree(scip, &(ancestors));
                  SCIPintListNodeFree(scip, &(revancestors));
               }
            }
         }

         /* try to contract adjacent terminals */
         if( g->grad[i] > 0 )
         {
            SCIP_Real mincost = FARAWAY;
            int ett = UNKNOWN;

            for( e1 = g->outbeg[i]; e1 != EAT_LAST; e1 = g->oeat[e1] )
            {
               i1 = g->head[e1];
               if( !g->mark[i1] )
                  continue;
               if( SCIPisLT(scip, g->cost[e1], mincost) )
               {
                  mincost = g->cost[e1];
                  if( Is_term(g->term[i1]) )
                     ett = e1;
               }
               else if( Is_term(g->term[i1]) && SCIPisLE(scip, g->cost[e1], mincost) )
               {
                  assert(SCIPisLT(scip, g->cost[e1], FARAWAY));
                  assert(SCIPisEQ(scip, g->cost[e1], mincost));
                  ett = e1;
               }
            }
            if( ett != UNKNOWN && SCIPisLE(scip, g->cost[ett], mincost) && SCIPisLE(scip, g->cost[ett], g->prize[i])
               && SCIPisLE(scip, g->cost[ett], g->prize[g->head[ett]]) )
            {
               i1 = g->head[ett];
               SCIPdebugMessage("contract tt %d->%d\n ", i, i1);
               assert(SCIPisLT(scip, mincost, FARAWAY));
               *fixed += g->cost[ett];
               (*count)++;
               SCIP_CALL( graph_knot_contractpc(scip, g, i, i1, i) );
               rerun = TRUE;
            }
         }
      }
   }

   if( !pc )
      g->mark[g->source[0]] = TRUE;
   SCIPdebugMessage("dirdeg %d Knots deleted\n", *count);

   /* free memory */
   SCIPfreeBufferArray(scip, &nodes2);
   SCIPfreeBufferArray(scip, &edges2);

   return SCIP_OKAY;
}

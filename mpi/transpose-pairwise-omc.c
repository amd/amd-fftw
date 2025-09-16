/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 * Copyright (C) 2022, Advanced Micro Devices, Inc. All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* Distributed transposes that overlap the memcpy operation with the MPI
   communication of asynchronous send-recv calls.
   This MPI transpose method builds on top of original transpose-pairwise
   and continues to use a sequence of carefully scheduled pairwise exchanges.
   This has the advantage that it can be done in-place, or out-of-place
   while preserving the input, using buffer space proportional to the
   local size divided by the number of processes (i.e. to the total array
   size divided by the number of processes squared). */

#include <string.h>
#include "mpi-transpose.h"

typedef struct {
     solver super;
     int preserve_input; /* preserve input even if DESTROY_INPUT was passed */
} S;

typedef struct {
     plan_mpi_transpose super;

     plan *cld1, *cld2, *cld2rest, *cld3;
     INT rest_Ioff, rest_Ooff;
     
     int n_pes, my_pe, *sched;
     INT *send_block_sizes, *send_block_offsets;
     INT *recv_block_sizes, *recv_block_offsets;
     R *send_block_bufs[2];
     MPI_Comm comm;
     int preserve_input;
} P;

static void transpose_chunks(int *sched, int n_pes, int my_pe,
                 INT *sbs, INT *sbo, INT *rbs, INT *rbo,
                 MPI_Comm comm,
                 R *I, R *O, R **bufs)
{
    if (sched) {
      int i;
      MPI_Status status;

      if (I == O) {
           //Two temporary block memory buffers: When MPI communication is 
           //performed on one buffer, memcpy is performed on the second.
           R *buf[2];
#if !defined(AMD_MPI_MALLOC_ONCE)
           buf[0] = (R*) MALLOC(sizeof(R) * sbs[0], BUFFERS);
           buf[1] = (R*) MALLOC(sizeof(R) * sbs[0], BUFFERS);
#else
           buf[0] = bufs[0];
           buf[1] = bufs[1];
#endif
           int pe = sched[0], pe2, j=0;
           MPI_Request send_request, recv_request;

#ifdef AMD_MPI_TRANSPOSE_LOGS
           printf("TRANSPOSE-PAIRWISE: n_pes[%d], my_pe[%d], first_pe[%d]\n", n_pes, my_pe, pe);
#endif
           //prolog : Memcpy of first data is done initially that needs to be sent to its remote rank
           i = 0;
           if (my_pe == pe) {
#ifdef AMD_MPI_TRANSPOSE_LOGS
               printf("PROLOG:(my_pe == pe)=> my_pe[%d], i[%d]:pe[%d], src_offset[%d], dst_offset[%d], size[%d]\n", my_pe, i, pe, sbo[pe], rbo[pe], sbs[pe]);
#endif
               if (rbo[pe] != sbo[pe])
                   memmove(O + rbo[pe], O + sbo[pe], sbs[pe] * sizeof(R));
               pe = sched[1];
               i++;
           }
           memcpy(buf[j], O + sbo[pe], sbs[pe] * sizeof(R));

           //overlapped execution region of memcpy and MPI communication : Asynchronous Send-Recv MPI opearations are used
           for (; i < (n_pes-1); ++i) {
               pe2 = sched[i+1];
               if (my_pe == pe2) {
#ifdef AMD_MPI_TRANSPOSE_LOGS
               printf("OVERLAP REGION:(my_pe == pe)=> my_pe[%d], i[%d]:pe[%d], src_offset[%d], dst_offset[%d], size[%d]\n", my_pe, i, pe2, sbo[pe2], rbo[pe2], sbs[pe2]);
#endif
                   if (rbo[pe2] != sbo[pe2])
                       memmove(O + rbo[pe2], O + sbo[pe2],
                               sbs[pe2] * sizeof(R));
                   continue;
               }
               else {
#ifdef AMD_MPI_TRANSPOSE_LOGS
               printf("OVERLAP REGION:(my_pe != pe)=> my_pe[%d], i[%d]:pe[%d], src_offset[%d], dst_offset[%d], size[%d]\n", my_pe, i, pe2, sbo[pe2], rbo[pe2], sbs[pe2]);
#endif
                   MPI_Isend(buf[j&0x1], (int) (sbs[pe]), FFTW_MPI_TYPE, pe, (my_pe * n_pes + pe) & 0xffff, comm, &send_request);
                   MPI_Irecv(O + rbo[pe], (int) (rbs[pe]), FFTW_MPI_TYPE, pe, (pe * n_pes + my_pe) & 0xffff, comm, &recv_request);
                   memcpy(buf[(j+1)&0x1], O + sbo[pe2], sbs[pe2] * sizeof(R));
                   pe = pe2;
                   //barrier
                   MPI_Wait(&send_request, MPI_STATUS_IGNORE);
                   MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
               }
               j++;
           }

           //epilog : Last send-recv to the remote rank is performed as normal blocking Send-Recv operation
           if (my_pe == pe) {
#ifdef AMD_MPI_TRANSPOSE_LOGS
               printf("EPILOG:(my_pe == pe)=> my_pe[%d], i[%d]:pe[%d], src_offset[%d], dst_offset[%d], size[%d]\n", my_pe, i, pe, sbo[pe], rbo[pe], sbs[pe]);
#endif
               if (rbo[pe] != sbo[pe])
                   memmove(O + rbo[pe], O + sbo[pe], sbs[pe] * sizeof(R));
           }
           else {
#ifdef AMD_MPI_TRANSPOSE_LOGS
               printf("EPILOG:(my_pe != pe)=> my_pe[%d], i[%d]:pe[%d], src_offset[%d], dst_offset[%d], size[%d]\n", my_pe, i, pe, sbo[pe], rbo[pe], sbs[pe]);
#endif
               MPI_Sendrecv(buf[j&0x1], (int) (sbs[pe]), FFTW_MPI_TYPE,
                       pe, (my_pe * n_pes + pe) & 0xffff,
                       O + rbo[pe], (int) (rbs[pe]),
                       FFTW_MPI_TYPE,
                       pe, (pe * n_pes + my_pe) & 0xffff,
                       comm, &status);
           }
#if !defined(AMD_MPI_MALLOC_ONCE)
           X(ifree)(buf[0]);
           X(ifree)(buf[1]);
#endif
      }
      else { /* I != O */
           for (i = 0; i < n_pes; ++i) {
            int pe = sched[i];
            if (my_pe == pe)
             memcpy(O + rbo[pe], I + sbo[pe], sbs[pe] * sizeof(R));
            else
             MPI_Sendrecv(I + sbo[pe], (int) (sbs[pe]),
                      FFTW_MPI_TYPE,
                      pe, (my_pe * n_pes + pe) & 0xffff,
                      O + rbo[pe], (int) (rbs[pe]),
                      FFTW_MPI_TYPE,
                      pe, (pe * n_pes + my_pe) & 0xffff,
                      comm, &status);
           }
      }
    }
}

static void apply(const plan *ego_, R *I, R *O)
{
     const P *ego = (const P *) ego_;
     plan_rdft *cld1, *cld2, *cld2rest, *cld3;

     /* transpose locally to get contiguous chunks */
     cld1 = (plan_rdft *) ego->cld1;
     if (cld1) {
      cld1->apply(ego->cld1, I, O);
      
      if (ego->preserve_input) I = O;

      /* transpose chunks globally */
      transpose_chunks(ego->sched, ego->n_pes, ego->my_pe,
               ego->send_block_sizes, ego->send_block_offsets,
               ego->recv_block_sizes, ego->recv_block_offsets,
               ego->comm, O, I, ego->send_block_bufs);
     }
     else if (ego->preserve_input) {
      /* transpose chunks globally */
      transpose_chunks(ego->sched, ego->n_pes, ego->my_pe,
               ego->send_block_sizes, ego->send_block_offsets,
               ego->recv_block_sizes, ego->recv_block_offsets,
               ego->comm, I, O, ego->send_block_bufs);

      I = O;
     }
     else {
      /* transpose chunks globally */
      transpose_chunks(ego->sched, ego->n_pes, ego->my_pe,
               ego->send_block_sizes, ego->send_block_offsets,
               ego->recv_block_sizes, ego->recv_block_offsets,
               ego->comm, I, I, ego->send_block_bufs);
     }

     /* transpose locally, again, to get ordinary row-major;
    this may take two transposes if the block sizes are unequal
    (3 subplans, two of which operate on disjoint data) */
     cld2 = (plan_rdft *) ego->cld2;
     cld2->apply(ego->cld2, I, O);
     cld2rest = (plan_rdft *) ego->cld2rest;
     if (cld2rest) {
      cld2rest->apply(ego->cld2rest,
              I + ego->rest_Ioff, O + ego->rest_Ooff);
      cld3 = (plan_rdft *) ego->cld3;
      if (cld3)
           cld3->apply(ego->cld3, O, O);
      /* else TRANSPOSED_OUT is true and user wants O transposed */
     }
}

static int applicable(const S *ego, const problem *p_,
              const planner *plnr, int n_pes)
{
     const problem_mpi_transpose *p = (const problem_mpi_transpose *) p_;
     /* Note: this is *not* UGLY for out-of-place, destroy-input plans;
    the planner often prefers transpose-pairwise to transpose-alltoall,
    at least with LAM MPI on my machine. */
     return (1
         && (!ego->preserve_input || (!NO_DESTROY_INPUTP(plnr)
                      && p->I != p->O))
         && ONLY_TRANSPOSEDP(p->flags)
         && (n_pes > 1));
}

static void awake(plan *ego_, enum wakefulness wakefulness)
{
     P *ego = (P *) ego_;
     X(plan_awake)(ego->cld1, wakefulness);
     X(plan_awake)(ego->cld2, wakefulness);
     X(plan_awake)(ego->cld2rest, wakefulness);
     X(plan_awake)(ego->cld3, wakefulness);
}

static void destroy(plan *ego_)
{
     P *ego = (P *) ego_;
     X(ifree0)(ego->sched);
     X(ifree0)(ego->send_block_sizes);
#if defined(AMD_MPI_MALLOC_ONCE)
     X(ifree)(ego->send_block_bufs[0]);
     X(ifree)(ego->send_block_bufs[1]);
#endif
     MPI_Comm_free(&ego->comm);
     X(plan_destroy_internal)(ego->cld3);
     X(plan_destroy_internal)(ego->cld2rest);
     X(plan_destroy_internal)(ego->cld2);
     X(plan_destroy_internal)(ego->cld1);
}

static void print(const plan *ego_, printer *p)
{
     const P *ego = (const P *) ego_;
     p->print(p, "(mpi-transpose-pairwise-omc%s%(%p%)%(%p%)%(%p%)%(%p%))", 
          ego->preserve_input==2 ?"/p":"",
          ego->cld1, ego->cld2, ego->cld2rest, ego->cld3);
}

/* Given a process which_pe and a number of processes npes, fills
   the array sched[npes] with a sequence of processes to communicate
   with for a deadlock-free, optimum-overlap all-to-all communication.
   (All processes must call this routine to get their own schedules.)
   The schedule can be re-ordered arbitrarily as long as all processes
   apply the same permutation to their schedules.

   The algorithm here is based upon the one described in:
       J. A. M. Schreuder, "Constructing timetables for sport
       competitions," Mathematical Programming Study 13, pp. 58-67 (1980). 
   In a sport competition, you have N teams and want every team to
   play every other team in as short a time as possible (maximum overlap
   between games).  This timetabling problem is therefore identical
   to that of an all-to-all communications problem.  In our case, there
   is one wrinkle: as part of the schedule, the process must do
   some data transfer with itself (local data movement), analogous
   to a requirement that each team "play itself" in addition to other
   teams.  With this wrinkle, it turns out that an optimal timetable
   (N parallel games) can be constructed for any N, not just for even
   N as in the original problem described by Schreuder.
*/
static void fill1_comm_sched(int *sched, int which_pe, int npes)
{
     int pe, i, n, s = 0;
     A(which_pe >= 0 && which_pe < npes);
     if (npes % 2 == 0) {
      n = npes;
      sched[s++] = which_pe;
     }
     else
      n = npes + 1;
     for (pe = 0; pe < n - 1; ++pe) {
      if (npes % 2 == 0) {
           if (pe == which_pe) sched[s++] = npes - 1;
           else if (npes - 1 == which_pe) sched[s++] = pe;
      }
      else if (pe == which_pe) sched[s++] = pe;

      if (pe != which_pe && which_pe < n - 1) {
           i = (pe - which_pe + (n - 1)) % (n - 1);
           if (i < n/2)
            sched[s++] = (pe + i) % (n - 1);
           
           i = (which_pe - pe + (n - 1)) % (n - 1);
           if (i < n/2)
            sched[s++] = (pe - i + (n - 1)) % (n - 1);
      }
     }
     A(s == npes);
}

/* Sort the communication schedule sched for npes so that the schedule
   on process sortpe is ascending or descending (!ascending).  This is
   necessary to allow in-place transposes when the problem does not
   divide equally among the processes.  In this case there is one
   process where the incoming blocks are bigger/smaller than the
   outgoing blocks and thus have to be received in
   descending/ascending order, respectively, to avoid overwriting data
   before it is sent. */
#ifdef AMD_FMV_AUTO
__attribute__((target_clones(TARGET_STRINGS)))
#endif
static void sort1_comm_sched(int *sched, int npes, int sortpe, int ascending)
{
     int *sortsched, i;
     sortsched = (int *) MALLOC(npes * sizeof(int) * 2, OTHER);
     fill1_comm_sched(sortsched, sortpe, npes);
     if (ascending)
      for (i = 0; i < npes; ++i)
           sortsched[npes + sortsched[i]] = sched[i];
     else
      for (i = 0; i < npes; ++i)
           sortsched[2*npes - 1 - sortsched[i]] = sched[i];
     for (i = 0; i < npes; ++i)
      sched[i] = sortsched[npes + i];
     X(ifree)(sortsched);
}

static plan *mkplan(const solver *ego_, const problem *p_, planner *plnr)
{
     const S *ego = (const S *) ego_;
     const problem_mpi_transpose *p;
     P *pln;
     plan *cld1 = 0, *cld2 = 0, *cld2rest = 0, *cld3 = 0;
     INT b, bt, vn, rest_Ioff, rest_Ooff;
     INT *sbs, *sbo, *rbs, *rbo;
     int pe, my_pe, n_pes, sort_pe = -1, ascending = 1;
     R *I, *O;
     static const plan_adt padt = {
          XM(transpose_solve), awake, print, destroy
     };

     UNUSED(ego);

     p = (const problem_mpi_transpose *) p_;
     vn = p->vn;
     I = p->I; O = p->O;

     MPI_Comm_rank(p->comm, &my_pe);
     MPI_Comm_size(p->comm, &n_pes);

     if (!applicable(ego, p_, plnr, n_pes))
          return (plan *) 0;

     b = XM(block)(p->nx, p->block, my_pe);
     
     if (!(p->flags & TRANSPOSED_IN)) { /* b x ny x vn -> ny x b x vn */
      cld1 = X(mkplan_f_d)(plnr, 
                   X(mkproblem_rdft_0_d)(X(mktensor_3d)
                             (b, p->ny * vn, vn,
                              p->ny, vn, b * vn,
                              vn, 1, 1),
                             I, O),
                   0, 0, NO_SLOW);
      if (XM(any_true)(!cld1, p->comm)) goto nada;
     }
     if (ego->preserve_input || NO_DESTROY_INPUTP(plnr)) I = O;

     if (XM(any_true)(!XM(mkplans_posttranspose)(p, plnr, I, O, my_pe,
                         &cld2, &cld2rest, &cld3,
                         &rest_Ioff, &rest_Ooff),
              p->comm)) goto nada;

     pln = MKPLAN_MPI_TRANSPOSE(P, &padt, apply);

     pln->cld1 = cld1;
     pln->cld2 = cld2;
     pln->cld2rest = cld2rest;
     pln->rest_Ioff = rest_Ioff;
     pln->rest_Ooff = rest_Ooff;
     pln->cld3 = cld3;
     pln->preserve_input = ego->preserve_input ? 2 : NO_DESTROY_INPUTP(plnr);

     MPI_Comm_dup(p->comm, &pln->comm);

     n_pes = (int) X(imax)(XM(num_blocks)(p->nx, p->block),
               XM(num_blocks)(p->ny, p->tblock));

     /* Compute sizes/offsets of blocks to exchange between processors */
     sbs = (INT *) MALLOC(4 * n_pes * sizeof(INT), PLANS);
     sbo = sbs + n_pes;
     rbs = sbo + n_pes;
     rbo = rbs + n_pes;
     b = XM(block)(p->nx, p->block, my_pe);
     bt = XM(block)(p->ny, p->tblock, my_pe);
     for (pe = 0; pe < n_pes; ++pe) {
      INT db, dbt; /* destination block sizes */
      db = XM(block)(p->nx, p->block, pe);
      dbt = XM(block)(p->ny, p->tblock, pe);

      sbs[pe] = b * dbt * vn;
      sbo[pe] = pe * (b * p->tblock) * vn;
      rbs[pe] = db * bt * vn;
      rbo[pe] = pe * (p->block * bt) * vn;

      if (db * dbt > 0 && db * p->tblock != p->block * dbt) {
           A(sort_pe == -1); /* only one process should need sorting */
           sort_pe = pe;
           ascending = db * p->tblock > p->block * dbt;
      }
     }
     pln->n_pes = n_pes;
     pln->my_pe = my_pe;
     pln->send_block_sizes = sbs;
     pln->send_block_offsets = sbo;
     pln->recv_block_sizes = rbs;
     pln->recv_block_offsets = rbo;

     if (my_pe >= n_pes) {
      pln->sched = 0; /* this process is not doing anything */
     }
     else {
      pln->sched = (int *) MALLOC(n_pes * sizeof(int), PLANS);
      fill1_comm_sched(pln->sched, my_pe, n_pes);
      if (sort_pe >= 0)
           sort1_comm_sched(pln->sched, n_pes, sort_pe, ascending);
     }

#if defined(AMD_MPI_MALLOC_ONCE)
     //Malloc temporary block-sized buffers for in-place transpose operation
     pln->send_block_bufs[0] = (R*) MALLOC(sizeof(R) * sbs[0], BUFFERS);
     pln->send_block_bufs[1] = (R*) MALLOC(sizeof(R) * sbs[0], BUFFERS);
#endif

     X(ops_zero)(&pln->super.super.ops);
     if (cld1) X(ops_add2)(&cld1->ops, &pln->super.super.ops);
     if (cld2) X(ops_add2)(&cld2->ops, &pln->super.super.ops);
     if (cld2rest) X(ops_add2)(&cld2rest->ops, &pln->super.super.ops);
     if (cld3) X(ops_add2)(&cld3->ops, &pln->super.super.ops);
     /* FIXME: should MPI operations be counted in "other" somehow? */

     return &(pln->super.super);

 nada:
     X(plan_destroy_internal)(cld3);
     X(plan_destroy_internal)(cld2rest);
     X(plan_destroy_internal)(cld2);
     X(plan_destroy_internal)(cld1);
     return (plan *) 0;
}

static solver *mksolver(int preserve_input)
{
     static const solver_adt sadt = { PROBLEM_MPI_TRANSPOSE, mkplan, 0 };
     S *slv = MKSOLVER(S, &sadt);
     slv->preserve_input = preserve_input;
     return &(slv->super);
}

void XM(transpose_pairwise_omc_register)(planner *p)
{
     int preserve_input;
     for (preserve_input = 0; preserve_input <= 1; ++preserve_input)
        REGISTER_SOLVER(p, mksolver(preserve_input));
}

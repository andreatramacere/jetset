

//=========================================================================================
//                   THREADS HANDLER
//=========================================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "Blazar_SED.h"
#include <pthread.h>


void threaded_j_evaluation(struct blob * pt, void *(*eval_j)(void *data), 
    double * j_nu_array, double * nu_array, double nu_start, double nu_stop, unsigned int I_MAX, unsigned int N_THREADS){
    unsigned int THREAD, CHUNK_SIZE, NU_INT_MAX, ACTUAL_N_THREADS;
    
    unsigned int capped = N_THREADS;
    if (capped < 1) capped = 1;
    if (capped > I_MAX + 1) capped = I_MAX + 1;
    
    if (N_THREADS>1){
        CHUNK_SIZE= (I_MAX+1)/capped;
        CHUNK_SIZE = (CHUNK_SIZE + 7) & ~7;  // round up to multiple of 8

    }else{
        N_THREADS = 1;
        CHUNK_SIZE= 0;
    }

    // derive ACTUAL_N_THREADS based on chunk_size
    if (CHUNK_SIZE < 1) {
        ACTUAL_N_THREADS = 1;
    } else {
        unsigned int needed = (I_MAX + 1 + CHUNK_SIZE - 1) / CHUNK_SIZE;
        ACTUAL_N_THREADS = (capped < needed) ? capped : needed;
    }

    struct j_args *thread_args = (struct j_args*) malloc(ACTUAL_N_THREADS * sizeof(struct j_args));
    if (thread_args == NULL) {
        struct j_args args;
        args.blob_pt = pt;
        args.NU_INT_START = 0;
        args.NU_INT_STOP = I_MAX;
        args.nu_array = nu_array;
        args.j_array = j_nu_array;
        eval_j(&args);
        return;
    }

    if ((ACTUAL_N_THREADS < 2) ){
        THREAD =0 ;
        thread_args[THREAD].blob_pt = pt;
        thread_args[THREAD].NU_INT_START = 0;
        thread_args[THREAD].NU_INT_STOP =I_MAX;
        thread_args[THREAD].nu_array = nu_array;
        thread_args[THREAD].j_array = j_nu_array;
        if (pt->core.verbose>0) {
            printf("NO THREAD, nu_start_int=%d, nu_stop_int=%d\n",thread_args[0].NU_INT_START,thread_args[0].NU_INT_STOP);
        }
        eval_j(thread_args);
    }else{

        pthread_t *threads = (pthread_t*) malloc(ACTUAL_N_THREADS * sizeof(pthread_t));
        if (threads == NULL) {
            thread_args[0].blob_pt = pt;
            thread_args[0].NU_INT_START = 0;
            thread_args[0].NU_INT_STOP =I_MAX;
            thread_args[0].nu_array = nu_array;
            thread_args[0].j_array = j_nu_array;
            eval_j(thread_args);
            free(thread_args);
            return;
        }
        unsigned int created = 0;

        for (THREAD = 0; THREAD <ACTUAL_N_THREADS; THREAD++) {
            
            thread_args[THREAD].blob_pt = pt;
            thread_args[THREAD].NU_INT_START = (THREAD*CHUNK_SIZE);
            if (THREAD == ACTUAL_N_THREADS -1){
                NU_INT_MAX = I_MAX;
            }else{
                NU_INT_MAX = thread_args[THREAD].NU_INT_START + CHUNK_SIZE -1;
            }
            thread_args[THREAD].NU_INT_STOP = min(NU_INT_MAX,I_MAX);
            thread_args[THREAD].nu_array = nu_array;
            thread_args[THREAD].j_array = j_nu_array;
            if (pt->core.verbose>0) {
                printf("THREAD=%d, nu_start_int<=%d, nu_stop_int<=%d\n",THREAD,thread_args[THREAD].NU_INT_START,thread_args[THREAD].NU_INT_STOP);
            }
            int result = pthread_create(&threads[created], NULL, eval_j, &thread_args[THREAD]);
            if (result != 0) {
                printf("Error creating thread %d\n", THREAD);
                eval_j(&thread_args[THREAD]);
            } else {
                created++;
            }
           
        }

        for (THREAD = 0; THREAD <created; THREAD++) {
            int result = pthread_join(threads[THREAD], NULL);
            if (result != 0) {
                  printf("Error joining thread %d\n", THREAD);
            }
        }
        free(threads);   
    }
    free(thread_args);
}

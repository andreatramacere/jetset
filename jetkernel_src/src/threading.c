

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
    unsigned int THREAD, CHUNK_SIZE, NU_INT_MAX;
    

    if (N_THREADS>1){
        CHUNK_SIZE= (I_MAX+1)/N_THREADS;
    }else{
        N_THREADS = 1;
        CHUNK_SIZE= 0;
    }

   
    struct j_args *thread_args = (struct j_args*) malloc(N_THREADS * sizeof(struct j_args));

    if ((N_THREADS < 2) ){
        THREAD =0 ;
        thread_args[THREAD].blob_pt = pt;
        thread_args[THREAD].NU_INT_START = 0;
        thread_args[THREAD].NU_INT_STOP =I_MAX;
        thread_args[THREAD].nu_array = nu_array;
        if (pt->verbose>0) {
            printf("NO THREAD, nu_start_int=%d, nu_stop_int=%d\n",thread_args[0].NU_INT_START,thread_args[0].NU_INT_STOP);
        }
        eval_j(thread_args);
    }else{
        //struct j_args  thread_args[N_THREADS];
        
        //pthread_t threads[N_THREADS];
        pthread_t *threads = (pthread_t*) malloc(N_THREADS * sizeof(pthread_t));

        for (THREAD = 0; THREAD <N_THREADS; THREAD++) {
            
            thread_args[THREAD].blob_pt = pt;
            thread_args[THREAD].NU_INT_START = (THREAD*CHUNK_SIZE);
            if (THREAD == N_THREADS -1){
                NU_INT_MAX = I_MAX;
            }else{
                NU_INT_MAX = thread_args[THREAD].NU_INT_START + CHUNK_SIZE -1;
            }
            thread_args[THREAD].NU_INT_STOP = min(NU_INT_MAX,I_MAX);

            thread_args[THREAD].nu_array = nu_array;
            if (pt->verbose>0) {
                printf("THREAD=%d, nu_start_int<=%d, nu_stop_int<=%d\n",THREAD,thread_args[THREAD].NU_INT_START,thread_args[THREAD].NU_INT_STOP);
            }
            int result = pthread_create(&threads[THREAD], NULL, eval_j, &thread_args[THREAD]);
            if (result != 0) {
                printf("Error creating thread %d\n", THREAD);
            }
           
        }

        for (THREAD = 0; THREAD <N_THREADS; THREAD++) {
            int result = pthread_join(threads[THREAD], NULL);
            if (result != 0) {
                  printf("Error joining thread %d\n", THREAD);
            }
        }
        free(threads);
        free(thread_args);
    }
}


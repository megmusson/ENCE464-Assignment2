#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#include <unistd.h>

/**
 * poisson.c
 * Implementation of a Poisson solver with Neumann boundary conditions.
 *
 * This template handles the basic program launch, argument parsing, and memory
 * allocation required to implement the solver *at its most basic level*. You
 * will likely need to allocate more memory, add threading support, account for
 * cache locality, etc...
 *
 * BUILDING:
 * gcc -o poisson poisson.c -lpthread
 *
 * [note: linking pthread isn't strictly needed until you add your
 *        multithreading code]
 *
 * TODO:
 * 1 - Read through this example, understand what it does and what it gives you
 *     to work with.
 * 2 - Implement the basic algorithm and get a correct output.
 * 3 - Add a timer to track how long your execution takes.
 * 4 - Profile your solution and identify weaknesses.
 * 5 - Improve it!
 * 6 - Remember that this is now *your* code and *you* should modify it however
 *     needed to solve the assignment.
 *
 * See the lab notes for a guide on profiling and an introduction to
 * multithreading (see also threads.c which is reference by the lab notes).
 */

#define ONE_SIX 1.0/6.0

// Global flag
// Set to true when operating in debug mode to enable verbose logging
static bool debug = false;


/**
 * @brief Solve Poissons equation for a given cube with Neumann boundary
 * conditions on all sides.
 *
 * @param n             The edge length of the cube. n^3 number of elements.
 * @param source        Pointer to the source term cube, a.k.a. forcing function.
 * @param iterations    Number of iterations to perform.
 * @param threads       Number of threads to use for solving.
 * @param delta         Grid spacing.
 * @return double*      Solution to Poissons equation.  Caller must free.
 */
 

typedef struct
{
    int start;          // Start index of the worker thread
    int end;            // End index of the worker thread
    int n;
    double *curr;
    double *next;
    float delta;
    double *source;
    int iterations;
} WorkerArgs;



void* worker (void* pargs)
{
    WorkerArgs* args = (WorkerArgs*)pargs;
    
    double i_next = 0;
    double i_prev = 0;
    double j_next = 0;
    double j_prev = 0;
    double k_next = 0;
    double k_prev = 0;
    
    int i = 0;
    int j = 0;
    int k = 0;

    int n = args->n;
    double *curr = args->curr;
    double *next = args->next;
    int start = args->start;
    float delta = args->delta;
    int end = args->end;
    double *source = args->source;

    for (k = start; k < end; k++)
    {
        for (j = 1; j < n-1; j++)
        {
            for (i = 1; i < n-1; i++)
            {
                i_next = curr[n*n*k + n*j + i+1];
                i_prev = curr[n*n*k + n*j + i-1];
        
                j_next = curr[n*n*k + n*(j+1) + i];
                j_prev = curr[n*n*k + n*(j-1) + i];
                
                k_next = curr[n*n*(k+1) + n*j + i];
                k_prev = curr[n*n*(k-1) + n*j + i];

                next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
            }
        }
    }     

    return NULL;
}
 
 
double* poisson_neumann (int n, double *source, int iterations, int threads, float delta)
{	
    double i_next = 0;
    double i_prev = 0;
    double j_next = 0;
    double j_prev = 0;
    double k_next = 0;
    double k_prev = 0;
    
    int i = 0;
    int j = 0;
    int k = 0;
    int t = 0;
    double range = n;

    if (debug)
    {
        printf ("Starting solver with:\n"
               "n = %i\n"
               "iterations = %i\n"
               "threads = %i\n"
               "delta = %f\n",
               n, iterations, threads, delta);
    }

    // Allocate some buffers to calculate the solution in
    double *curr = (double*)calloc (n * n * n, sizeof (double));
    double *next = (double*)calloc (n * n * n, sizeof (double));

    // Ensure we haven't run out of memory
    if (curr == NULL || next == NULL)
    {
        fprintf (stderr, "Error: ran out of memory when trying to allocate %i sized cube\n", n);
        exit (EXIT_FAILURE);
    }

    // TODO: solve Poisson's equation for the given inputs
    
    // Storage for the thread handles and arguments
    // will exist for the entire lifetime of the program.
    for (t = 0; t < iterations; t++)
    {

        pthread_t threads_index[threads];
        WorkerArgs args[threads];
        
        for (int i = 0; i < threads; i++)
        {
            // Fill in the arguments to the worker
            
            args[i].n = n;
            args[i].curr = curr; // setting the pointer
            args[i].next = next;
            args[i].start = (range * i) / threads;
            args[i].end = (range * (i + 1)) / threads;
            args[i].delta = delta;
            args[i].iterations = iterations;
            args[i].source = source;

            // Create the worker thread
            if (pthread_create (&threads_index[i], NULL, &worker, &args[i]) != 0)
            {
                fprintf (stderr, "Error creating worker thread!\n");
                // return EXIT_FAILURE;
                exit(EXIT_FAILURE);
            }
        }

        // ****************** Initial i plane**************
        i = 0;
        for (k = 1; k < n-1; k++)
        {
            for (j = 1; j < n-1; j++)
            {
                i_next = curr[n*n*k + n*j + i+1];
                i_prev = i_next;
                j_next = curr[n*n*k + n*(j+1) + i];
                j_prev = curr[n*n*k + n*(j-1) + i];
                k_next = curr[n*n*(k+1) + n*j + i];
                k_prev = curr[n*n*(k-1) + n*j + i];
                next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
            }
        }
        // Edge cases on initial i plane
        k = 0;
        for (j = 1; j < n-1; j++)
        {
            k_next = curr[n*n*(k+1) + n*j + i];
            k_prev = k_next;
            j_next = curr[n*n*k + n*(j+1) + i];
            j_prev = curr[n*n*k + n*(j-1) + i];
            i_next = curr[n*n*k + n*j + i+1];
            i_prev = i_next;
            next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        }
        k = n-1;
        for (j = 1; j < n-1; j++)
        {
            k_prev = curr[n*n*(k-1) + n*j + i];
            k_next = k_prev;
            j_next = curr[n*n*k + n*(j+1) + i];
            j_prev = curr[n*n*k + n*(j-1) + i];
            i_next = curr[n*n*k + n*j + i+1];
            i_prev = i_next;
            next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        }
        j = 0;
        for (k = 1; k < n-1; k++)
        {
            j_next = curr[n*n*k + n*(j+1) + i];
            j_prev = j_next;
            k_next = curr[n*n*(k+1) + n*j + i];
            k_prev = curr[n*n*(k-1) + n*j + i];
            i_next = curr[n*n*k + n*j + i+1];
            i_prev = i_next;
            next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        }
        j = n-1;
        for (k = 1; k < n-1; k++)
        {
            j_prev = curr[n*n*k + n*(j-1) + i];
            j_next = j_prev;
            k_next = curr[n*n*(k+1) + n*j + i];
            k_prev = curr[n*n*(k-1) + n*j + i];
            i_next = curr[n*n*k + n*j + i+1];
            i_prev = i_next;
            next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        }
        // Corner cases for inital i plane
        j = 0;
        k = 0;
        i_next = curr[n*n*k + n*j + i+1];
        i_prev = i_next;
        j_next = curr[n*n*k + n*(j+1) + i];
        j_prev = j_next;
        k_next = curr[n*n*(k+1) + n*j + i];
        k_prev = k_next;
        next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);

        j = n-1;
        k = 0;
        i_next = curr[n*n*k + n*j + i+1];
        i_prev = i_next;
        j_prev = curr[n*n*k + n*(j-1) + i];
        j_next = j_prev;
        k_next = curr[n*n*(k+1) + n*j + i];
        k_prev = k_next;
        next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);

        j = 0;
        k = n-1;
        i_next = curr[n*n*k + n*j + i+1];
        i_prev = i_next;
        j_next = curr[n*n*k + n*(j+1) + i];
        j_prev = j_next;
        k_prev = curr[n*n*(k-1) + n*j + i];
        k_next = k_prev;
        next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);

        j = n-1;
        k = n-1;
        i_next = curr[n*n*k + n*j + i+1];
        i_prev = i_next;
        j_prev = curr[n*n*k + n*(j-1) + i];
        j_next = j_prev;
        k_prev = curr[n*n*(k-1) + n*j + i];
        k_next = k_prev;
        next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);


        // ****************** Final i plane **************
        i = n-1;
        for (k = 1; k < n-1; k++)
        {
            for (j = 1; j < n-1; j++)
            {
                i_prev = curr[n*n*k + n*j + i-1];
                i_next = i_prev;
                j_next = curr[n*n*k + n*(j+1) + i];
                j_prev = curr[n*n*k + n*(j-1) + i];
                k_next = curr[n*n*(k+1) + n*j + i];
                k_prev = curr[n*n*(k-1) + n*j + i];
                next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
            }
        }
        // Edge cases for final i plane
        k = 0;
        for (j = 1; j < n-1; j++)
        {
            k_next = curr[n*n*(k+1) + n*j + i];
            k_prev = k_next;
            j_next = curr[n*n*k + n*(j+1) + i];
            j_prev = curr[n*n*k + n*(j-1) + i];
            i_prev = curr[n*n*k + n*j + i-1];
            i_next = i_prev;
            next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        }
        k = n-1;
        for (j = 1; j < n-1; j++)
        {
            k_prev = curr[n*n*(k-1) + n*j + i];
            k_next = k_prev;
            j_next = curr[n*n*k + n*(j+1) + i];
            j_prev = curr[n*n*k + n*(j-1) + i];
            i_prev = curr[n*n*k + n*j + i-1];
            i_next = i_prev;
            next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        }
        j = 0;
        for (k = 1; k < n-1; k++)
        {
            j_next = curr[n*n*k + n*(j+1) + i];
            j_prev = j_next;
            k_next = curr[n*n*(k+1) + n*j + i];
            k_prev = curr[n*n*(k-1) + n*j + i];
            i_prev = curr[n*n*k + n*j + i-1];
            i_next = i_prev;
            next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        }
        j = n-1;
        for (k = 1; k < n-1; k++)
        {
            j_prev = curr[n*n*k + n*(j-1) + i];
            j_next = j_prev;
            k_next = curr[n*n*(k+1) + n*j + i];
            k_prev = curr[n*n*(k-1) + n*j + i];
            i_prev = curr[n*n*k + n*j + i-1];
            i_next = i_prev;
            next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        }
        // Corner cases for final i plane
        j = 0;
        k = 0;
        i_prev = curr[n*n*k + n*j + i-1];
        i_next = i_prev;
        j_next = curr[n*n*k + n*(j+1) + i];
        j_prev = j_next;
        k_next = curr[n*n*(k+1) + n*j + i];
        k_prev = k_next;
        next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);

        j = n-1;
        k = 0;
        i_prev = curr[n*n*k + n*j + i-1];
        i_next = i_prev;
        j_prev = curr[n*n*k + n*(j-1) + i];
        j_next = j_prev;
        k_next = curr[n*n*(k+1) + n*j + i];
        k_prev = k_next;
        next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);

        j = 0;
        k = n-1;
        i_prev = curr[n*n*k + n*j + i-1];
        i_next = i_prev;
        j_next = curr[n*n*k + n*(j+1) + i];
        j_prev = j_next;
        k_prev = curr[n*n*(k-1) + n*j + i];
        k_next = k_prev;
        next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);

        j = n-1;
        k = n-1;
        i_prev = curr[n*n*k + n*j + i-1];
        i_next = i_prev;
        j_prev = curr[n*n*k + n*(j-1) + i];
        j_next = j_prev;
        k_prev = curr[n*n*(k-1) + n*j + i];
        k_next = k_prev;
        next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);

        

        // ****************** Initial j plane**************
        j = 0; 
        for (k = 1; k < n-1; k++)
        {
            for (i = 1; i < n-1; i++)
            {
                j_next = curr[n*n*k + n*(j+1) + i];
                j_prev = j_next;
                i_next = curr[n*n*k + n*j + i+1];
                i_prev = curr[n*n*k + n*j + i-1];
                k_next = curr[n*n*(k+1) + n*j + i];
                k_prev = curr[n*n*(k-1) + n*j + i];
                next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
            }
        }
        // Edge cases for initial j plane
        k = 0;
        for (i = 1; i < n-1; i++)
        {
            k_next = curr[n*n*(k+1) + n*j + i];
            k_prev = k_next;
            i_next = curr[n*n*k + n*j + i+1];
            i_prev = curr[n*n*k + n*j + i-1];
            j_next = curr[n*n*k + n*(j+1) + i];
            j_prev = j_next;
            next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        }
        k = n-1;
        for (i = 1; i < n-1; i++)
        {
            k_prev = curr[n*n*(k-1) + n*j + i];
            k_next = k_prev;
            i_next = curr[n*n*k + n*j + i+1];
            i_prev = curr[n*n*k + n*j + i-1];
            j_next = curr[n*n*k + n*(j+1) + i];
            j_prev = j_next;
            next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        }
        // i = 0;
        // for (k = 1; k < n-1; k++)
        // {
        //     k_next = curr[n*n*(k+1) + n*j + i];
        //     k_prev = curr[n*n*(k-1) + n*j + i];
        //     i_next = curr[n*n*k + n*j + i+1];
        //     i_prev = i_next;
        //     j_next = curr[n*n*k + n*j + i+1];
        //     j_prev = j_next;
        //     next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        // }
        // i = n-1;
        // for (k = 1; k < n-1; k++)
        // {
        //     k_next = curr[n*n*(k+1) + n*j + i];
        //     k_prev = curr[n*n*(k-1) + n*j + i];
        //     i_prev = curr[n*n*k + n*j + i+1];
        //     i_next = i_prev;
        //     j_next = curr[n*n*k + n*j + i+1];
        //     j_prev = j_next;
        //     next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        // }
        // // Corner cases for inital i plane
        // i = 0;
        // k = 0;
        // i_next = curr[n*n*k + n*j + i+1];
        // i_prev = i_next;
        // j_next = curr[n*n*k + n*j + i+1];
        // j_prev = j_next;
        // k_next = curr[n*n*k + n*j + i+1];
        // k_prev = k_next;
        // next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);

        // i = n-1;
        // k = 0;
        // i_prev = curr[n*n*k + n*j + i+1];
        // i_next = i_prev;
        // j_next = curr[n*n*k + n*j + i+1];
        // j_prev = j_next;
        // k_next = curr[n*n*k + n*j + i+1];
        // k_prev = k_next;
        // next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);

        // i = 0;
        // k = n-1;
        // i_next = curr[n*n*k + n*j + i+1];
        // i_prev = i_next;
        // j_next = curr[n*n*k + n*j + i+1];
        // j_prev = j_next;
        // k_prev = curr[n*n*k + n*j + i-1];
        // k_next = k_prev;
        // next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);

        // i = n-1;
        // k = n-1;
        // i_prev = curr[n*n*k + n*j + i+1];
        // i_next = i_prev;
        // j_next = curr[n*n*k + n*j + i+1];
        // j_prev = j_next;
        // k_prev = curr[n*n*k + n*j + i+1];
        // k_next = k_prev;
        // next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);

        




        // ****************** Final j plane**************
        j = n-1; 
        for (k = 1; k < n-1; k++)
        {
            for (i = 1; i < n-1; i++)
            {
                j_prev = curr[n*n*k + n*(j-1) + i];
                j_next = j_prev;
                i_next = curr[n*n*k + n*j + i+1];
                i_prev = curr[n*n*k + n*j + i-1];
                k_next = curr[n*n*(k+1) + n*j + i];
                k_prev = curr[n*n*(k-1) + n*j + i];
                next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
            }
        }
        // Edge cases for Final j plane
        k = 0;
        for (i = 1; i < n-1; i++)
        {
            k_next = curr[n*n*(k+1) + n*j + i];
            k_prev = k_next;
            i_next = curr[n*n*k + n*j + i+1];
            i_prev = curr[n*n*k + n*j + i-1];
            j_prev = curr[n*n*k + n*(j-1) + i];
            j_next = j_prev;
            next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        }
        k = n-1;
        for (i = 1; i < n-1; i++)
        {
            k_prev = curr[n*n*(k-1) + n*j + i];
            k_next = k_prev;
            i_next = curr[n*n*k + n*j + i+1];
            i_prev = curr[n*n*k + n*j + i-1];
            j_prev = curr[n*n*k + n*(j-1) + i];
            j_next = j_prev;
            next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        }
        
        // ****************** Inital k plane**************
        k = 0;
        // for (j = 1; j < n-1; j++)
        // {
        //     for (i = 1; i < n-1; i++)
        //     {
        //         k_next = curr[n*n*(k+1) + n*j + i];
        //         k_prev = k_next;
        //         i_next = curr[n*n*k + n*j + i+1];
        //         i_prev = curr[n*n*k + n*j + i-1];
        //         j_next = curr[n*n*k + n*(j+1) + i];
        //         j_prev = curr[n*n*k + n*(j-1) + i];
        //         next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
        //     }
        // }

        for (j = 1; j < n-1; j++)
        {
            for (i = 1; i < n-1; i++)
            {
                k_next = curr[n*n*(k+1) + n*j + i];
                k_prev = k_next;
                i_next = curr[n*n*k + n*j + i+1];
                i_prev = curr[n*n*k + n*j + i-1];
                j_next = curr[n*n*k + n*(j+1) + i];
                j_prev = curr[n*n*k + n*(j-1) + i];
                next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
            }
        }
        // ****************** Final k plane**************
        k = n-1;
        for (j = 1; j < n-1; j++)
        {
            for (i = 1; i < n-1; i++)
            {
                k_prev = curr[n*n*(k-1) + n*j + i];
                k_next = k_prev;
                i_next = curr[n*n*k + n*j + i+1];
                i_prev = curr[n*n*k + n*j + i-1];
                j_next = curr[n*n*k + n*(j+1) + i];
                j_prev = curr[n*n*k + n*(j-1) + i];
                next[n*n*k + n*j + i] = ONE_SIX*(i_next + i_prev + j_next + j_prev + k_next + k_prev - delta*delta*source[n*n*k + n*j + i]);
            }
        }
        
        double *temp = next;
        next = curr;
        curr = temp;
        // Wait for all the threads to finish using join ()
        for (int i = 0; i < threads; i++)
        {
            pthread_join (threads_index[i], NULL);
        }
    }

    
	
    // Free one of the buffers and return the correct answer in the other.
    // The caller is now responsible for free'ing the returned pointer.
    free (next);

    if (debug)
    {
        printf ("Finished solving.\n");
    }

    return curr;
}



int main (int argc, char **argv)
{
    // Default settings for solver
    int iterations = 1;
    int n = 5;
    int threads = 1;
    float delta = 1;
    
    
    

    // parse the command line arguments
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp (argv[i], "-h") == 0 || strcmp (argv[i], "--help") == 0)
        {
            printf ("Usage: poisson [-n size] [-i iterations] [-t threads] [--debug]\n");
            return EXIT_SUCCESS;
        }

        if (strcmp (argv[i], "-n") == 0)
        {
            if (i == argc - 1)
            {
                fprintf (stderr, "Error: expected size after -n!\n");
                return EXIT_FAILURE;
            }

            n = atoi (argv[++i]);
        }

        if (strcmp (argv[i], "-i") == 0)
        {
            if (i == argc - 1)
            {
                fprintf (stderr, "Error: expected iterations after -i!\n");
                return EXIT_FAILURE;
            }

            iterations = atoi (argv[++i]);
        }

        if (strcmp (argv[i], "-t") == 0)
        {
            if (i == argc - 1)
            {
                fprintf (stderr, "Error: expected threads after -t!\n");
                return EXIT_FAILURE;
            }

            threads = atoi (argv[++i]);
        }

        if (strcmp (argv[i], "--debug") == 0)
        {
            debug = true;
        }
    }

    // Ensure we have an odd sized cube
    if (n % 2 == 0)
    {
        fprintf (stderr, "Error: n should be an odd number!\n");
        return EXIT_FAILURE;
    }

    // Create a source term with a single point in the centre
    double *source = (double*)calloc (n * n * n, sizeof (double));
    if (source == NULL)
    {
        fprintf (stderr, "Error: failed to allocated source term (n=%i)\n", n);
        return EXIT_FAILURE;
    }

    source[(n * n * n) / 2] = 1;

    // Calculate the resulting field with Neumann conditions
    double *result = poisson_neumann (n, source, iterations, threads, delta);

    // Print out the middle slice of the cube for validation
    // for (int x = 0; x < n; ++x)
    // {
    //     for (int y = 0; y < n; ++y)
    //     {
    //         printf ("%0.5f ", result[((n / 2) * n + y) * n + x]);
    //     }
    //     printf ("\n");
    // }

    for (int x = 0; x < n; ++x)
    {
        for (int y = 0; y < n; ++y)
        {
            printf ("%0.5f ", result[(y * n + x) * n + (n/2)]);
        }
        printf ("\n");
    }


    free (source);
    free (result);

    return EXIT_SUCCESS;
}

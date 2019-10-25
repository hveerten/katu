/* The MIT License
 * 
 * Copyright (c) 2019 Bruno Jiménez <brunojimen@gmail.com>
 *
 * Adapted from the code by John Schember <john@nachtimwald.com> used
 * in the poddown program.
 * References:
 *  - https://github.com/user-none/poddown
 *  - https://nachtimwald.com/2019/04/12/thread-pool-in-c/
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE
 */

#include <stdlib.h>
#include <pthread.h>

#include <stdio.h>

#include "tpool.h"

struct _thread_pool_work_item_t
{
    thread_fun_t      fun;
    void              *args;
    struct _thread_pool_work_item_t *next;
};

static thread_pool_work_item_t *thread_pool_create_work_item(thread_fun_t f, void *args)
{
    if(f == NULL)
        return NULL;

    thread_pool_work_item_t *item = NULL;

    item       = calloc(1, sizeof(thread_pool_work_item_t));
    item->fun  = f;
    item->args = args;
    item->next = NULL;

    return item;
}

static void thread_pool_destroy_work_item(thread_pool_work_item_t *item)
{
    if(item == NULL)
        return;

    free(item);
}

static thread_pool_work_item_t *thread_pool_get_work_item(thread_pool_t *pool)
{
    thread_pool_work_item_t *item;

    if(pool == NULL)
        return NULL;

    item = pool->head;
    if(item == NULL)
        return NULL;

    if(item->next == NULL)
    {
        pool->head = NULL;
        pool->tail = NULL;
    }
    else
        pool->head = item->next;

    return item;
}

static void *tpool_worker(void *args)
{
    thread_pool_t *pool = args;
    thread_pool_work_item_t *item = NULL;

    while(1)
    {
        pthread_mutex_lock(&(pool->work_mutex));
        /* Keep running until told to stop. */
        if (pool->stop)
            break;

        /* If there is no work in the queue wait in the conditional until
         * there is work to take. */
        if (pool->head == NULL)
            pthread_cond_wait(&(pool->work_cond), &(pool->work_mutex));

        /* Try to pull work from the queue. */
        item = thread_pool_get_work_item(pool);
        pool->num_working++;
        pthread_mutex_unlock(&(pool->work_mutex));

        /* Call the work function and let it process.
         *
         * item can legitimately be NULL. Since multiple threads from the pool
         * will wake when there is work, a thread might not get any work. 1
         * piece of work and 2 threads, both will wake but with 1 only work 1
         * will get the work and the other won't.
         *
         * num_working has been increment and item could be NULL. While it's
         * not true there is work processing the thread is considered working
         * because it's not waiting in the conditional. Pedantic but...
         */
        if(item != NULL)
        {
            item->fun(item->args);
            thread_pool_destroy_work_item(item);
        }

        pthread_mutex_lock(&(pool->work_mutex));
        pool->num_working--;
        /* Since we're in a lock no work can be added or removed form the queue.
         * Also, the working_cnt can't be changed (except the thread holding the lock).
         * At this point if there isn't any work processing and if there is no work
         * signal this is the case. */
        if (!pool->stop && pool->num_working == 0 && pool->head == NULL)
            pthread_cond_signal(&(pool->working_cond));
        pthread_mutex_unlock(&(pool->work_mutex));
    }

    pool->num_alive--;
    if(pool->num_alive == 0)
        pthread_cond_signal(&(pool->work_cond));
    pthread_mutex_unlock(&(pool->work_mutex));
    return NULL;
}

void thread_pool_create(thread_pool_t *pool, size_t num)
{
    unsigned int i;
    pthread_t thread;

    pool->stop = false;
    pool->num_alive = num;
    pool->num_working = 0;

    pthread_mutex_init(&(pool->work_mutex), NULL);
    pthread_cond_init(&(pool->work_cond), NULL);
    pthread_cond_init(&(pool->working_cond), NULL);

    pool->head = NULL;
    pool->tail = NULL;

    /* Create the requested number of thread and detach them. */
    for(i = 0; i < num; i++)
    {
        pthread_create(&thread, NULL, tpool_worker, pool);
        pthread_detach(thread);
    }

    thread_pool_wait(pool);
}

void thread_pool_destroy(thread_pool_t *pool)
{
    thread_pool_work_item_t *item;
    thread_pool_work_item_t *item2;

    if (pool == NULL)
        return;

    /* Take all work out of the queue and destroy it. */
    pthread_mutex_lock(&(pool->work_mutex));
    item = pool->head;
    while(item != NULL)
    {
        item2 = item->next;
        thread_pool_destroy_work_item(item);
        item = item2;
    }

    /* Tell the worker threads to stop. */
    pool->stop = true;
    pthread_cond_broadcast(&(pool->work_cond));
    pthread_mutex_unlock(&(pool->work_mutex));

    /* Wait for all threads to stop. */
    thread_pool_wait(pool);

    pthread_mutex_destroy(&(pool->work_mutex));
    pthread_cond_destroy(&(pool->work_cond));
    pthread_cond_destroy(&(pool->working_cond));

    free(pool);
}

bool thread_pool_add_work(thread_pool_t *pool, thread_fun_t f, void *args)
{
    thread_pool_work_item_t *item = NULL;

    if(pool == NULL)
        return false;

    item = thread_pool_create_work_item(f, args);

    if(item == NULL)
        return false;
    /*fprintf(stderr,"ping\n");*/

    pthread_mutex_lock(&(pool->work_mutex));
    if(pool->head == NULL)
    {
        pool->head = item;
        pool->tail = pool->head;
    }
    else
    {
        pool->tail->next = item;
        pool->tail       = item;
    }

    pthread_cond_broadcast(&(pool->work_cond));
    pthread_mutex_unlock(&(pool->work_mutex));

    return true;
}

void thread_pool_wait(thread_pool_t *pool)
{
    if(pool == NULL)
        return;

    pthread_mutex_lock(&(pool->work_mutex));
    while(1)
    {
        /* working_cond is dual use. It signals when we're not stopping but the
         * working_cnt is 0 indicating there isn't any work processing. If we
         * are stopping it will trigger when there aren't any threads running. */
        if((!pool->stop && (pool->num_working != 0 || pool->head != NULL)) ||
           ( pool->stop && pool->num_alive   != 0))
            pthread_cond_wait(&(pool->working_cond), &(pool->work_mutex));
        else
            break;
    }
    pthread_mutex_unlock(&(pool->work_mutex));

    /*fprintf(stderr,"%zd\n", pool->num_working);*/
    /*fprintf(stderr,"%p\n", pool->head);*/
    /*fprintf(stderr,"\n");*/
}

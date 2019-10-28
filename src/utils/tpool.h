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

#ifndef TREAD_POOL_H
#define TREAD_POOL_H
#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <pthread.h>


typedef void *(*thread_fun_t)(void *args);

struct _thread_pool_work_item_t;
typedef struct _thread_pool_work_item_t thread_pool_work_item_t;

typedef struct _thread_pool_t
{
    thread_pool_work_item_t *head;
    thread_pool_work_item_t *tail;

    pthread_mutex_t work_mutex; 
    pthread_cond_t  work_cond;
    pthread_cond_t  working_cond;

    size_t num_working; 
    size_t num_alive;  

    bool   stop;        
} thread_pool_t;

void thread_pool_init(thread_pool_t *pool, size_t num);
void thread_pool_clear(thread_pool_t *pool);

bool thread_pool_add_work(thread_pool_t *pool, thread_fun_t f, void *args);

void thread_pool_wait(thread_pool_t *pool);

#endif /* end of include guard: TREAD_POOL_H */

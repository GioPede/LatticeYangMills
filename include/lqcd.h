/******************************************************************************
* 
* MIT License
* 
* Copyright (c) 2018 Giovanni Pederiva
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
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
******************************************************************************/

/*! \file      lqcd.h 
 *  \brief     Main include file for all headers
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 * 
 *  This file provides a simple access to all of the classes 
 *  defined in the project. Includes all actions, apps, observables
 *  as well as base objects.
 */

#pragma once

#include "Actions/action.h"
#include "Actions/actionlist.h"
#include "Apps/applist.h"
#include "InputOutput/io.h"
#include "Math/latticemath.h"
#include "Observables/observable.h"
#include "Observables/observablelist.h"
#include "ParallelTools/parallel.h"
#include "Utils/latticeunits.h"

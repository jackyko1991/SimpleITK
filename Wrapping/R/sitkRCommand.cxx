/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

// The python header defines _POSIX_C_SOURCE without a preceding #undef
#undef _POSIX_C_SOURCE
#undef _XOPEN_SOURCE
#include <Rinternals.h>

#include "sitkRCommand.h"
#include "sitkExceptionObject.h"

#include <iostream>


namespace itk
{
namespace simple
{


RCommand::RCommand()
  : m_Object(NULL)
{
}

RCommand::~RCommand()
{
  this->m_Object = NULL;
}

void RCommand::SetCallbackRCallable(SEXP o)
{
  if (o != this->m_Object)
    {
    // Don't think we have to do anything
    // fancy with references - R handles that internally??
    // store the new object
    this->m_Object = o;

    }
}

SEXP RCommand::GetCallbackRCallable()
{
  return this->m_Object;
}

void RCommand::Execute()
{
  // if null do nothing
  if (!this->m_Object)
    {
    return;
    }

  // make sure that the CommandCallable is in fact callable
  if (!isFunction(this->m_Object))
    {
    // we throw a standard ITK exception: this makes it possible for
    // our standard CableSwig exception handling logic to take this
    // through to the invoking R process
    sitkExceptionMacro(<<"R Callable is not a callable R object, "
                       <<"or it has not been set.");
    }
  else
    {
    SEXP result;

    if ( !isEnvironment(FRAME(this->m_Object))) {
    sitkExceptionMacro(<<"FRAME is not the right way to get a function environment,");

    }

    result = eval(this->m_Object, FRAME(this->m_Object));

    }
}


} // namespace simple
} // namespace itk

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
#include <Rinternals.h>

#include "sitkRCommand.h"
#include "sitkExceptionObject.h"

#include <iostream>


namespace itk
{
namespace simple
{


RCommand::RCommand()
  : m_Object(R_NilValue)
{
}

RCommand::~RCommand()
{
  this->m_Object = NULL;
  this->m_Environ = NULL;
}

void RCommand::SetCallbackRCallable(SEXP obj)
{
  if (obj != this->m_Object)
    {
    // Don't think we have to do anything
    // fancy with references - R handles that internally??
    // store the new object
    this->m_Object = obj;

    }
}
void RCommand::SetCallbackREnviron(SEXP rho)
{
  if (rho != this->m_Environ)
    {
    // Don't think we have to do anything
    // fancy with references - R handles that internally??
    // store the new object
    this->m_Environ = rho;

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

  else
    {
    SEXP result;
    // retrieve the environment for passing to eval
    result = eval(this->m_Object, this->m_Environ);
    if (result == NULL)
      sitkExceptionMacro(<<"There was an error executing the "
                         <<"R Callable.");

    }
}


} // namespace simple
} // namespace itk

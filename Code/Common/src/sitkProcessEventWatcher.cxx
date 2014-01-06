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
#include "sitkProcessEventWatcher.h"
#include "itkProcessObject.h"


namespace itk
{
namespace simple
{

ProcessEventWatcher::ProcessEventWatcher( ProcessObject &po )
  : m_Process(po),
    m_Progress(0.0),
    m_IterationNumber(0)

{
}

std::string ProcessEventWatcher::GetName( ) const
{
  return m_Process.GetName();
}

float ProcessEventWatcher::GetProgress( void ) const
{
  return this->m_Progress;
}

unsigned int ProcessEventWatcher::GetIterationNumber( void ) const
{
  return this->m_IterationNumber;
}

ProcessObject &ProcessEventWatcher::GetProcessObject()
{
  return this->m_Process;
}

void ProcessEventWatcher::OnAbortEvent( )
{
}

void ProcessEventWatcher::OnEndEvent( )
{
}

void ProcessEventWatcher::OnStartEvent( )
{
  this->m_Progress = 0.0;
  this->m_IterationNumber = 0;
}

void ProcessEventWatcher::OnProgressEvent( )
{
  this->m_Progress = m_Process.GetProgress();
}


void ProcessEventWatcher::OnIterationEvent( )
{
  ++this->m_IterationNumber;
}

} // end namespace simple
} // end namespace itk

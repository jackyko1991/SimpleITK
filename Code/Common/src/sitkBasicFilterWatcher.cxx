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
#include "sitkBasicFilterWatcher.h"

#include "itkTimeProbe.h"

#include <iomanip>


namespace itk
{
namespace simple
{

BasicFilterWatcher::BasicFilterWatcher(ProcessObject& o, const std::string &comment)
  : ProcessEventWatcher(o),
    m_TimeProbe(new TimeProbe()),
    m_Comment(comment)
{
  // Initialize state
  m_Quiet = false;

  BasicFilterWatcher *self = this;
  m_MemberCommands[0].SetCallbackFunction(self, &Self::OnProgressEvent);
  o.AddCommand(sitkProgressEvent, &m_MemberCommands[0]);
  m_MemberCommands[1].SetCallbackFunction(self, &Self::OnAbortEvent);
  o.AddCommand(sitkAbortEvent, &m_MemberCommands[1]);
  m_MemberCommands[2].SetCallbackFunction(self, &Self::OnIterationEvent);
  o.AddCommand(sitkIterationEvent, &m_MemberCommands[2]);
  m_MemberCommands[3].SetCallbackFunction(self, &Self::OnStartEvent);
  o.AddCommand(sitkStartEvent, &m_MemberCommands[3]);
  m_MemberCommands[4].SetCallbackFunction(self, &Self::OnEndEvent);
  o.AddCommand(sitkEndEvent, &m_MemberCommands[4]);
}

BasicFilterWatcher::~BasicFilterWatcher(void)
{
  delete m_TimeProbe;
}

/** Callback method to show the ProgressEvent */
void BasicFilterWatcher::OnProgressEvent( )
{
  this->Superclass::OnProgressEvent( );

  this->ShowProgress();
}

/** Callback method to show the AbortEvent */
void BasicFilterWatcher::OnAbortEvent( )
{
  this->Superclass::OnAbortEvent( );
  m_TimeProbe->Stop();
    std::cout << std::endl << "-------Aborted" << std::endl << std::flush;
}

/** Callback method to show the IterationEvent */
void BasicFilterWatcher::OnIterationEvent( )
{
  this->Superclass::OnIterationEvent( );
  this->ShowProgress();
}

/** Callback method to show the StartEvent */
void BasicFilterWatcher::OnStartEvent( )
{
  this->Superclass::OnStartEvent( );
  m_TimeProbe->Start();
  this->ShowProgress();
}

/** Callback method to show the EndEvent */
void BasicFilterWatcher::OnEndEvent( )
{
  this->Superclass::OnEndEvent( );

  this->ShowProgress();

  m_TimeProbe->Stop();

  std::cout << std::endl << this->GetName() << " took " << m_TimeProbe->GetMean() << std::endl;
}

void BasicFilterWatcher::ShowProgress( void ) const
{
  std::streamsize oldWidth = std::cout.width();
  std::streamsize oldPrecision = std::cout.precision();
  std::ios_base::fmtflags oldFlags = std::cout.flags();

  const unsigned int precision = 1;


  std::cout << "\r-------- ";

  std::cout << this->GetName()
              << " Progress: ";

  if ( this->GetProgress() != 0.0 )
      {
      std::cout << "[ ";
      std::cout << std::setprecision( precision ) << std::setw( precision+4 ) << std::right << std::fixed << this->GetProgress()*100;
      std::cout << std::setiosflags( oldFlags ) << std::setw( oldWidth ) << std::setprecision( oldPrecision );
      std::cout << " ] ";
      }

  if ( this->GetIterationNumber() != 0 )
    {
    std::cout << " ( ";
    std::cout << std::setw( 6 ) << std::right << this->GetIterationNumber();
    std::cout << std::setiosflags( oldFlags ) << std::setw( oldWidth );
    std::cout << " ) ";
    }

  std::cout << "\r" << std::flush;
}

} // end namespace simple
} // end namespace itk

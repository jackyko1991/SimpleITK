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
#ifndef __sitkProcessEventWatcher_h
#define __sitkProcessEventWatcher_h

#include "sitkProcessObject.h"

namespace itk
{

class ProcessObject;

namespace simple {

/** \class ProcessEventWatcher
 * \brief Provides a mechanism for displaying progress of multiple filters.
 *
 */
class SITKCommon_EXPORT ProcessEventWatcher
{
public:

  typedef ProcessEventWatcher Self;

  ProcessEventWatcher(ProcessObject &po );

  virtual std::string GetName( ) const;
  virtual float GetProgress( ) const;
  virtual unsigned int GetIterationNumber( ) const;

protected:

  virtual ProcessObject &GetProcessObject();

  virtual void OnAbortEvent( );
  virtual void OnEndEvent( );
  virtual void OnStartEvent( );
  virtual void OnProgressEvent( );
  virtual void OnIterationEvent( );

private:

  ProcessObject& m_Process;

  float        m_Progress;
  unsigned int m_IterationNumber;

};

} // end namespace simple
} // end namespace itk

#endif // __sitkProcessEventWatcher_h

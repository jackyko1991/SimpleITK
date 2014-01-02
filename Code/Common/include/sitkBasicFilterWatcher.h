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
#ifndef __sitkBasicFilterWatcher_h
#define __sitkBasicFilterWatcher_h

#include "sitkProcessEventWatcher.h"
#include <string>
#include <vector>

#include "sitkMemberCommand.h"

namespace itk
{

class TimeProbe;

namespace simple
{

/** \class BasicFilterWatcher
 */
class  SITKCommon_EXPORT BasicFilterWatcher
  : public ProcessEventWatcher
{
public:
  typedef BasicFilterWatcher  Self;
  typedef ProcessEventWatcher Superclass;

  /** Constructor. Takes a ProcessObject to monitor and an optional
   * comment string that is prepended to each event message. */
  BasicFilterWatcher(ProcessObject& o, const std::string &comment="");

  /** Deconstructor. To clean up the queue, just to be safe */
  ~BasicFilterWatcher(void);

  /** Methods to control the verbosity of the messages. Quiet
   * reporting limits the information emitted at a ProgressEvent. */
  void QuietOn( void ) {m_Quiet = true;};
  void QuietOff( void ) {m_Quiet = false;};

  /** Set/Get the quiet mode boolean. If true, verbose progess is
    * reported. */
  void SetQuiet(bool val) {m_Quiet=val;};
  bool GetQuiet( void ) {return m_Quiet;};

  /** Get the comment for the watcher. */
  const std::string &GetComment() const {return m_Comment;};

protected:

  /** Callback method to show the ProgressEvent */
  virtual void OnProgressEvent( );

  /** Callback method to show the AbortEvent */
  virtual void OnAbortEvent( );

  /** Callback method to show the IterationEvent */
  virtual void OnIterationEvent( );

  /** Callback method to show the StartEvent */
  virtual void OnStartEvent( );

  /** Callback method to show the EndEvent */
  virtual void OnEndEvent( );

  virtual void ShowProgress( void ) const;

private:

  TimeProbe*          m_TimeProbe;
  std::string         m_Comment;
  bool                m_Quiet;

  MemberCommand m_MemberCommands[5];

  /** Copy constructor */
  BasicFilterWatcher(const BasicFilterWatcher& ); // not implemented

  /** operator=  */
  void operator=(const BasicFilterWatcher& ); // not implemented

};

} // end namespace simple
} // end namespace itk

#endif //

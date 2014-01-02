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
#ifndef __sitkCommand_h
#define __sitkCommand_h


#include "sitkCommon.h"
#include "sitkNonCopyable.h"

#include <set>

namespace itk {
namespace simple {

class ProcessObject;

/** \class Command
 * \brief Base class for user overrideables
 */
class SITKCommon_EXPORT Command:
    protected NonCopyable
{
public:

  typedef Command Self;

  Command();

  /** Destruct or. */
  virtual ~Command(void);

  virtual void Execute(void);

protected:
  friend class itk::simple::ProcessObject;

  #ifndef SWIG
  virtual size_t AddProcessObject(itk::simple::ProcessObject *o);
  virtual size_t RemoveProcessObject(const itk::simple::ProcessObject *o);
  #endif

private:

  /** Copy constructor */
  Command(const Command& );

  /** operator=  */
  Command &operator=(const Command& );

  // a set of objects who use us
  std::set<itk::simple::ProcessObject*> m_ReferencedObjects;
};

} // end namespace simple
} // end namespace itk

#endif

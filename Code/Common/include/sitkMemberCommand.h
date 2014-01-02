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
#ifndef __sitkMemberCommand_h
#define __sitkMemberCommand_h


#include "sitkCommand.h"


#if defined SITK_HAS_STLTR1_TR1_FUNCTIONAL
#include <tr1/functional>
#elif defined SITK_HAS_STLTR1_FUNCTIONAL
#include <functional>
#else
#error "No system tr1 functional available"
#endif

namespace itk {
namespace simple {

/** \class MemberCommand
 * \brief Base class for user overrideables
 */
class SITKCommon_EXPORT MemberCommand:
    public Command
{
public:

  typedef MemberCommand Self;

  MemberCommand();

  /** Destruct or. */
    virtual ~MemberCommand(void);

  virtual void Execute(void);

  template <class T>
    void SetCallbackFunction ( T *object, void(T::* pMemberFunction )() )
  {
    m_Function = std::tr1::bind(pMemberFunction, object);
  }

  void SetCallbackFunction ( void(* pFunction )() )
  {
    m_Function = pFunction;
  }

private:

  /** Copy constructor */
  MemberCommand(const Command& );

  /** operator=  */
  MemberCommand &operator=(const Command& );

  typedef std::tr1::function<void()> FunctionObjectType;
  FunctionObjectType m_Function;

};

} // end namespace simple
} // end namespace itk

#endif

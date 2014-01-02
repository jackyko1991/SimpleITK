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

#ifndef __sitkPyCommand_h
#define __sitkPyCommand_h

#include "sitkCommand.h"



#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

namespace itk
{
namespace simple
{

/** \class PyCommand
 *  \brief Command subclass that calls a Python callable object, e.g.
 *  a Python function.
 *
 * With this class, arbitrary Python callable objects (e.g. functions)
 * can be associated with an instance to be used in AddObserver calls.
 *
 * Based of the WrapITK itkPyCommand class originally contributed by
 * Charl P. Botha <cpbotha |AT| ieee.org>.
 */
class PyCommand
  : public itk::simple::Command
{
public:
  // Standard "Self" typedef.
  typedef PyCommand         Self;


  PyCommand();
  ~PyCommand();


  /**
   * Assign a Python callable object to be used.  You don't have to keep
   * a binding to the callable, PyCommand will also take out a reference
   * to make sure the Callable sticks around.
   */
  void SetCommandCallable(PyObject *obj);

  PyObject * GetCommandCallable();

  virtual void Execute(void);

  #ifndef SWIG
  inline void SetOwnedByProcessObjects(bool o) {this->m_OwnedByProcessObjects = o;}
  inline bool GetOwnedByProcessObjects() const {return this->m_OwnedByProcessObjects;}
  inline void OwnedByProcessObjectsOn() {this->SetOwnedByProcessObjects(true);}
  inline void OwnedByProcessObjectsOff() {this->SetOwnedByProcessObjects(false);}
  #endif

protected:
  void PyExecute();

  virtual size_t RemoveProcessObject(const itk::simple::ProcessObject *o);

  PyCommand(const Self&);
  PyCommand & operator=(const Self&);

private:
  PyObject *m_Object;

  bool m_OwnedByProcessObjects;
};

} // namespace simple
} // namespace itk

#endif // _sitkPyCommand_h

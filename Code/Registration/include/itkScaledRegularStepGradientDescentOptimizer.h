/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkRegularStepGradientDescentOptimizer.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkScaledRegularStepGradientDescentOptimizer_h
#define __itkScaledRegularStepGradientDescentOptimizer_h

#include "itkRegularStepGradientDescentBaseOptimizer.h"

namespace itk
{


/** \class ScaledRegularStepGradientDescentOptimizer
 * \brief Implement a gradient descent optimizer
 *
 * \ingroup Numerics  Optimizers
 *
 */
class  ScaledRegularStepGradientDescentOptimizer :
    public RegularStepGradientDescentBaseOptimizer
{
public:
  /** Standard class typedefs. */
  typedef ScaledRegularStepGradientDescentOptimizer   Self;
  typedef RegularStepGradientDescentBaseOptimizer     Superclass;
  typedef SmartPointer<Self>                          Pointer;
  typedef SmartPointer<const Self>                    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ScaledRegularStepGradientDescentOptimizer,
                RegularStepGradientDescentBaseOptimizer );

  /** Cost function typedefs. */
  typedef Superclass::CostFunctionType        CostFunctionType;
  typedef CostFunctionType::Pointer           CostFunctionPointer;


protected:
  ScaledRegularStepGradientDescentOptimizer() {};

  /** Advance one step along the corrected gradient taking into
   * account the steplength represented by factor.
   * This method is invoked by AdvanceOneStep. It is expected
   * to be overrided by optimization methods in non-vector spaces
   * \sa AdvanceOneStep */
  virtual void StepAlongGradient( double factor, const DerivativeType & transformedGradient )
  {

    itkDebugMacro(<<"factor = " << factor << "  transformedGradient= " << transformedGradient );

    const unsigned int spaceDimension =
      m_CostFunction->GetNumberOfParameters();

    ParametersType newPosition( spaceDimension );
    ParametersType currentPosition = this->GetCurrentPosition();

    ScalesType scales = this->GetScales();

    for(unsigned int j=0; j<spaceDimension; j++)
      {
      newPosition[j] = currentPosition[j] + transformedGradient[j] * factor / scales[j];
      }

    itkDebugMacro(<<"new position = " << newPosition );

    this->SetCurrentPosition( newPosition );

    }

private:
  ScaledRegularStepGradientDescentOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#endif

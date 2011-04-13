#ifndef __sitkShapeLabelMapFilter_h
#define __sitkShapeLabelMapFilter_h

#include "sitkMacro.h"
#include "sitkMemberFunctionFactory.h"
#include "sitkImage.h"


namespace itk
{
namespace simple
{

/** \class ShapeLabelMapFilter
 * \brief Compute the shape attributes for a label map
 */
class ShapeLabelMapFilter {
public:
  typedef ShapeLabelMapFilter Self;

  // function pointer type
  typedef void (Self::*MemberFunctionType)( const Image& );

  ShapeLabelMapFilter();

  std::string GetName() const { return "ShapleLabelMapFilter"; }

  // Print ourselves out
  std::string ToString() const;

  void Execute ( const Image& image );


private:

  template <class TImageType> void ExecuteInternal ( const Image& inImage );

  // friend to get access to executeInternal member
  friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

  std::auto_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;
};

void ShapeLabelMap( const Image& image );

}
}
#endif

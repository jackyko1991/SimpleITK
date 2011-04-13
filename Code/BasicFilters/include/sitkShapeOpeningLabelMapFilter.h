#ifndef __sitkShapeOpeningLabelMapFilter_h
#define __sitkShapeOpeningLabelMapFilter_h

#include "sitkMacro.h"
#include "sitkMemberFunctionFactory.h"
#include "sitkImage.h"


namespace itk
{
namespace simple
{

/** \class ShapeOpeningLabelMapFilter
 * \brief Convert an binary image into multiple labels
 */
class ShapeOpeningLabelMapFilter {
public:
  typedef ShapeOpeningLabelMapFilter Self;

  // function pointer type
  typedef Image (Self::*MemberFunctionType)( const Image& );

  ShapeOpeningLabelMapFilter();

  // Print ourselves out
  std::string ToString() const;

  Image Execute ( const Image& image );


private:

  template <class TImageType> Image ExecuteInternal ( const Image& inImage );

  // friend to get access to executeInternal member
  friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

  std::auto_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;
};

Image ShapeOpeningLabelMap( const Image& image );

}
}
#endif

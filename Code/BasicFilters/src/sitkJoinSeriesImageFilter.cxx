#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkLabelMap.h"
#include "itkLabelObject.h"
#include "itkNumericTraits.h"
#include "itkNumericTraitsVariableLengthVectorPixel.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageToVectorImageFilter.h"

#include "itkJoinSeriesImageFilter.h"
#include "sitkJoinSeriesImageFilter.h"

namespace itk {
namespace simple {

//----------------------------------------------------------------------------

  Image JoinSeries ( const std::vector<Image> &images, double spacing, double origin ) {
    JoinSeriesImageFilter filter;
    filter.SetSpacing( spacing );
    filter.SetOrigin( origin );
    return filter.Execute ( images );
  }


//
// Default constructor that initializes parameters
//
JoinSeriesImageFilter::JoinSeriesImageFilter ()
{
  this->m_Spacing = 1.0;
  this->m_Origin = 0.0;
  this->m_MemberFactory.reset( new detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();
}

//
// ToString
//
std::string JoinSeriesImageFilter::ToString() const
  {
  std::ostringstream out;
  out << "itk::simple::JoinSeriesImageFilter\n"
      << "\tSpacing: " << this->m_Spacing << std::endl
      << "\tOrigin: " << this->m_Origin << std::endl;
  return out.str();
  }

JoinSeriesImageFilter::Self& JoinSeriesImageFilter::SetSpacing ( double s )
{
  this->m_Spacing = s;
  return *this;
}

double JoinSeriesImageFilter::GetSpacing() const
{
  return this->m_Spacing;
}

JoinSeriesImageFilter::Self& JoinSeriesImageFilter::JoinSeriesImageFilter::SetOrigin ( double s )
{
  this->m_Origin = s;
  return *this;
}

double JoinSeriesImageFilter::GetOrigin () const
{
  return this->m_Origin;
}


std::string JoinSeriesImageFilter::GetName() const
{
  return std::string ( "JoinSeries");
}



//
// Execute
//
Image JoinSeriesImageFilter::Execute ( const std::vector<Image> &images )
  {

    if ( images.empty() )
      {
      sitkExceptionMacro( << this->GetName() << " requires at least one image " );
      }

    PixelIDValueType type = images.front().GetPixelIDValue();
    unsigned int dimension = images.front().GetDimension();

    std::cout << "type: " << type << " dimension: " << dimension << std::endl;

    for ( std::vector<Image>::const_iterator i = images.begin(); i != images.end(); ++i )
      {
      // todo check that all the images are the same type and same size
      }

    return this->m_MemberFactory->GetMemberFunction( type, dimension )( images );
  }

//----------------------------------------------------------------------------

//
// ExecuteInternal
//
template <class TImageType>
Image JoinSeriesImageFilter::ExecuteInternal ( const std::vector<Image> &inImages )
{
  std::vector< typename TImageType::ConstPointer > images( inImages.size() );

  // cast all the images to raw pointers
  for ( unsigned int i = 0; i < inImages.size(); ++i )
    {
    images[i] = dynamic_cast<const TImageType*> ( inImages[i].GetImageBase() );

    if ( images[i].IsNull() )
      {
      sitkExceptionMacro( << "Could not cast input image " << i << " to proper type" );
      }
   }

  typedef itk::Image< typename TImageType::PixelType, TImageType::ImageDimension+1 > OutputImageType;

  typedef itk::JoinSeriesImageFilter<TImageType, OutputImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  for (typename  std::vector<typename TImageType::ConstPointer >::const_iterator i = images.begin(); i != images.end(); ++i )
    {
    filter->PushBackInput ( i->GetPointer() );
    }

  filter->SetSpacing( m_Spacing );
  filter->SetOrigin( m_Origin );
  filter->Update();

  return Image( filter->GetOutput() );
  }

} // end namespace simple
} // end namespace itk

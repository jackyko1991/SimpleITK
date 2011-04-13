#include "sitkShapeLabelMapFilter.h"

#include "itkShapeLabelMapFilter.h"
#include "itkConvertLabelMapFilter.h"

#include "itkLabelMap.h"
#include "itkLabelObject.h"
#include "itkShapeLabelObject.h"

namespace itk {
namespace simple {

ShapeLabelMapFilter::ShapeLabelMapFilter()
{

  typedef LabelPixelIDTypeList PixelIDTypeList;
  this->m_MemberFactory.reset( new detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();
}

// Print ourselves out
std::string ShapeLabelMapFilter::ToString() const
{
  std::ostringstream out;
  out << "itk::simple::ShapeLabelMapFilter\n";
  return out.str();
}

void ShapeLabelMapFilter::Execute ( const Image& image )
{
  PixelIDValueType type = image.GetPixelIDValue();
  unsigned int dimension = image.GetDimension();

  return this->m_MemberFactory->GetMemberFunction( type, dimension )( image );
}


template <class TImageType>
void ShapeLabelMapFilter::ExecuteInternal ( const Image& inImage )
{
  typedef TImageType InputImageType;
  typedef itk::ShapeLabelObject<typename InputImageType::LabelType, InputImageType::ImageDimension > ShapeLabelObjectType;
  typedef itk::LabelMap< ShapeLabelObjectType  > ShapeLabelMapType;

  typename InputImageType::ConstPointer image =
    dynamic_cast <const InputImageType*> ( inImage.GetImageBase() );

  if ( image.IsNull() )
    {
    sitkExceptionMacro( << "Could not cast input image to proper type" );
    }

  typedef itk::ConvertLabelMapFilter<  InputImageType, ShapeLabelMapType  > ConvertFilterType;
  typename ConvertFilterType::Pointer convert = ConvertFilterType::New();
  convert->SetInput( image );

  typedef itk::ShapeLabelMapFilter<ShapeLabelMapType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(convert->GetOutput() );
  filter->InPlaceOn();
  filter->Update();

  typename ShapeLabelMapType::Pointer labelMap = filter->GetOutput();
  for ( size_t label = 1; label < labelMap->GetNumberOfLabelObjects(); ++label )
    {
    const ShapeLabelObjectType *lo = labelMap->GetLabelObject( label );
    std::cout << "label: " << label
              << "\t" << lo->GetNumberOfPixels()
              << "\t" << lo->GetPhysicalSize()
              << "\t" << lo->GetRegionElongation()
              << "\t" << lo->GetSizeRegionRatio()
      //<< "\t" << lo->GetCentroid()
      //<< "\t" << lo->GetBoundingBox()
              << "\t" << lo->GetNumberOfPixelsOnBorder()
              << "\t" << lo->GetPerimeterOnBorder()
              << "\t" << lo->GetFeretDiameter()
      //<< "\t" << lo->GetPrincipalMoments()
      //<< "\t" << lo->GetPrincipalAxes()
              << "\t" << lo->GetElongation()
              << "\t" << lo->GetPerimeter()
              << "\t" << lo->GetRoundness()
              << "\t" << lo->GetEquivalentSphericalRadius()
              << "\t" << lo->GetEquivalentSphericalPerimeter()
      //<< "\t" << lo->GetEquivalentEllipsoidDiameter()
              << "\t" << lo->GetFlatness()
              << "\t" << lo->GetPerimeterOnBorderRatio()
              << std::endl;
    }

  //return out;
}


void ShapeLabelMap( const Image& image )
{
  return ShapeLabelMapFilter().Execute( image );
}
}
}


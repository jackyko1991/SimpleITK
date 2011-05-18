#ifndef __sitkJoinSeriesImageFilter_h
#define __sitkJoinSeriesImageFilter_h

#include "sitkMultiImageFilter.h"

#include <memory>

namespace itk {
  namespace simple {

    /** \class JoinSeriesImageFilter
     * \brief Join a series of slices into a volume
     */
    class JoinSeriesImageFilter
      : public MultiImageFilter {
    public:
      typedef JoinSeriesImageFilter Self;

      // Type List Setup
      typedef BasicPixelIDTypeList PixelIDTypeList;

      /**
       * Default Constructor that takes no arguments and initializes
       * default parameters
       */
      JoinSeriesImageFilter();

      Self& SetSpacing ( double s );
      double GetSpacing() const;

      Self& SetOrigin ( double s );
      double GetOrigin () const;

      /** Name of this class */
      virtual std::string GetName() const;

      // Print ourselves out
      virtual std::string ToString() const;

      virtual Image Execute ( const std::vector<Image> &);

    private:

      /** Spacing in the joining direction */
      double m_Spacing;

      /** Dimesion along which to extract the slice */
      double m_Origin;

      // Macro that instantiate the member function dispatching
      typedef Image (Self::*MemberFunctionType)(  const std::vector<Image> & );
      template <class TImageType> Image ExecuteInternal (  const std::vector<Image>& images );
      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;
      std::auto_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;
    };

  Image JoinSeries ( const std::vector<Image>&, double spacing = 1.0, double origin = 0.0 );
  }
}
#endif

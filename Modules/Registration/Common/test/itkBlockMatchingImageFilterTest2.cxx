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
// itkBlockMatchingImageFilterTest2 tp1_brain_regularized.nii.gz tp2_brain_regularized.nii.gz tp1_brain_mask_regularized.nii.gz displ.txt

#include <ctime>
#include <iostream>
#include "itkIndex.h"
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLineIterator.h"
#include "itkMaskFeaturePointSelectionFilter.h"
#include "itkBlockMatchingImageFilter.h"

#include "preselected_points.h"

#include <iostream>
#include <fstream>


int itkBlockMatchingImageFilterTest2( int argc, char * argv[] )
{
itk::MultiThreader::SetGlobalMaximumNumberOfThreads( 1 );

  std::cerr << "This is Test2 " << std::endl;
  if( argc < 5 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << " itkitkBlockMatchingImageFilterTest2 fixedImageFile[1] movingImageFile[2] maskImageFile[3] displacementsOutput[4]";
    return EXIT_FAILURE;
    }

  typedef unsigned short                  InputPixelType;
  typedef itk::RGBPixel<InputPixelType>  OutputPixelType;

  typedef itk::Image< InputPixelType,  3 >       InputImageType;
  typedef itk::Image< OutputPixelType, 3 >       OutputImageType;

  typedef itk::ImageFileReader< InputImageType >  ReaderType;

  typedef itk::MaskFeaturePointSelectionFilter< InputImageType >  FilterType;
  typedef FilterType::FeaturePointsType                           PointSetType;

  typedef FilterType::PointType       PointType;
  typedef FilterType::InputImageType  ImageType;

  // parameters
  FilterType::SizeType blockRadius = { 3, 3, 3 }; // 7x7x7
  FilterType::SizeType searchRadius = { 10, 10, 10 }; // 21x21x21
  FilterType::SizeType blockStep = { 1, 1, 1 };

  //Set up the readers
  ReaderType::Pointer readerFixed = ReaderType::New();
  readerFixed->SetFileName( argv[1] );
  readerFixed->Update();
  InputImageType::Pointer fixedImage = readerFixed->GetOutput();

  ReaderType::Pointer readerMoving = ReaderType::New();
  readerMoving->SetFileName( argv[2] );
  readerMoving->Update();
  InputImageType::Pointer movingImage = readerMoving->GetOutput();

  ReaderType::Pointer readerMask = ReaderType::New();
  readerMask->SetFileName( argv[3] );
  readerMask->Update();
  InputImageType::Pointer maskImage = readerMask->GetOutput();

  //// set origin to 0s
  //const float origin[3] = { 0.0, 0.0, 0.0 };
  //fixedImage->SetOrigin( origin );
  //movingImage->SetOrigin( origin );
  //maskImage->SetOrigin( origin );



  // Set up filter
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( fixedImage );
  filter->SetMaskImage( maskImage );

  // Set parameters
  filter->SetNonConnectivity( FilterType::VERTEX_CONNECTIVITY ); //26
  filter->SetBlockRadius( blockRadius );
  filter->SetSelectFraction( 0.05 ); // reject 95%
  filter->ComputeStructureTensorsOff(); // no tensors

  //s.Fill( 7 );
  //filter->SetBlockHalfWindow( s );

  std::cout << "Feature selection: " << filter << std::endl;

  try
    {
      const clock_t begin = std::clock();

//      filter->Update();

      const clock_t end = std::clock();

      std::cout << "FS execution time: "<< ( end - begin ) /  CLOCKS_PER_SEC << "sec" << std::endl << std::endl;
    }
  catch ( itk::ExceptionObject &err )
    {
      ( &err )->Print( std::cerr );
      return EXIT_FAILURE;
    }

  // create PointSet from preselected feature points array
  PointSetType::Pointer preSelectedPointSet = PointSetType::New();
  PointSetType::PointsContainer::Pointer points = PointSetType::PointsContainer::New();

  for ( unsigned i = 0; i < preselected_points_count; i++ )
    {
      PointSetType::PointType point;
      // -0.5 to convert from clatz's mapping
      point[0] = preselected_points[i * 3] - 0.5;
      point[1] = preselected_points[i * 3 + 1] - 0.5;
      point[2] = preselected_points[i * 3 + 2] - 0.5;

      points->InsertElement( i, point );
    }
  preSelectedPointSet->SetPoints(points);

  // at this time we have feature points
  typedef itk::BlockMatchingImageFilter< InputImageType >  BMFilterType;

  BMFilterType::Pointer BMFilter = BMFilterType::New();

  // inputs (all required)
  BMFilter->SetFixedImage( fixedImage );
  BMFilter->SetMovingImage( movingImage );
  BMFilter->SetFeaturePoints( preSelectedPointSet ); //  BMFilter->SetFeaturePoints( filter->GetOutput() );

  // BM parameters
  BMFilter->SetBlockRadius( blockRadius );
  BMFilter->SetSearchRadius( searchRadius );
  BMFilter->SetBlockStep( blockStep );

  std::cout << "Block matching: " << BMFilter << std::endl;
  try
    {
      const clock_t begin = std::clock();

      BMFilter->Update();

      const clock_t end = std::clock();

      std::cout << "BM execution time: "<< ( end - begin ) /  CLOCKS_PER_SEC << "sec" << std::endl << std::endl;
    }
  catch ( itk::ExceptionObject &err )
    {
      ( &err )->Print( std::cerr );
      return EXIT_FAILURE;
    }

// __DEBUG ONLY__ < < <
  // write displacements to file argv[4]
  std::cout << "BM done, wrtiting displacements to file <" << argv[4] << ">" << std::endl;
  std::ofstream fdisp(argv[4]);
  fdisp << "" << std::endl;

  unsigned pcount = BMFilter->GetOutput()->GetPoints()->size();
  std::cout << pcount << " points will be written to file" << std::endl;

  fdisp << pcount << std::endl;
  for ( unsigned i = 0; i < pcount; i++ )
    {
      BMFilterType::FeaturePointsPhysicalCoordinates coor = BMFilter->GetOutput()->GetPoints()->ElementAt(i);
      BMFilterType::DisplacementsVector              disp = BMFilter->GetOutput()->GetPointData()->ElementAt(i);
      fdisp << coor[0] << " " << coor[1] << " " << coor[2] << "\t"
            << disp[0] << " " << disp[1] << " " << disp[2] << std::endl;
    }

  fdisp <<  std::endl;
  fdisp.close();
// __DEBUG ONLY__ > > >


abort();

  //Set up the writer
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  WriterType::Pointer writer = WriterType::New();

  typedef itk::ImageRegionConstIterator< InputImageType >         InputIteratorType;
  InputIteratorType inputIterator( fixedImage, fixedImage->GetBufferedRegion() );
  typedef itk::ImageRegionIterator< OutputImageType > OutputIteratorType;

  OutputImageType::Pointer outputImage = OutputImageType::New();
  outputImage->CopyInformation( fixedImage );
  outputImage->SetBufferedRegion( fixedImage->GetBufferedRegion() );
  outputImage->SetRequestedRegion( fixedImage->GetRequestedRegion() );
  outputImage->Allocate();

  OutputIteratorType outputIterator( outputImage, outputImage->GetBufferedRegion() );
  inputIterator.GoToBegin();
  outputIterator.GoToBegin();

  // Copy input image to output image
  while (!outputIterator.IsAtEnd())
    {
    OutputPixelType rgbPixel;
    rgbPixel.SetRed( inputIterator.Get() );
    rgbPixel.SetGreen( inputIterator.Get());
    rgbPixel.SetBlue( inputIterator.Get() );
    outputIterator.Set( rgbPixel );
    ++outputIterator;
    ++inputIterator;
    }

  //Highlight the feature points identified in the output image
  typedef PointSetType::PointsContainer::ConstIterator                        PointIteratorType;
  typedef BMFilterType::DisplacementsType::PointDataContainer::ConstIterator  PointDataIteratorType;

  PointIteratorType pointItr =
          filter->GetOutput()->GetPoints()->Begin();
  PointIteratorType pointEnd =
          filter->GetOutput()->GetPoints()->End();
  PointDataIteratorType displItr =
          BMFilter->GetOutput()->GetPointData()->Begin();

  // define colors
  OutputPixelType red;
  red.SetRed( 255.0 );
  red.SetGreen( 0.0 );
  red.SetBlue( 0.0 );

  OutputPixelType green;
  green.SetRed( 0.0 );
  green.SetGreen( 255.0 );
  green.SetBlue( 0.0 );

  OutputPixelType blue;
  blue.SetRed( 0.0 );
  blue.SetGreen( 0.0 );
  blue.SetBlue( 255.0 );

  OutputImageType::IndexType index;
  while ( pointItr != pointEnd )
    {
    if ( outputImage->TransformPhysicalPointToIndex(pointItr.Value(), index) )
      {
        OutputImageType::IndexType displ;
        outputImage->TransformPhysicalPointToIndex( pointItr.Value() + displItr.Value(), displ );

        // draw line between old and new location of a point in blue
        itk::LineIterator< OutputImageType > lineIter( outputImage, index, displ );
        for ( lineIter.GoToBegin(); !lineIter.IsAtEnd(); ++lineIter )
          {
            lineIter.Set( blue );
          }

        // mark old location of a point in green
        outputImage->SetPixel(index, green);

        // mark new location of a point in red
        outputImage->SetPixel(displ, red);
      }
    pointItr++;
    displItr++;
    }

  writer->SetFileName( argv[2] );
  writer->SetInput( outputImage );
  writer->Update();


  return EXIT_SUCCESS;
}

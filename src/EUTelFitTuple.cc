
// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Revised by Martin Klassen
// Date 2014.08.06

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// this processor is built only if USE_AIDA and USE_GEAR are defined
#if defined(USE_GEAR) && ( defined(USE_AIDA) || defined(MARLIN_USE_AIDA) )

// eutelescope inlcudes
#include "EUTelFitTuple.h"
#include "EUTelVirtualCluster.h"
#include "EUTelSparseClusterImpl.h"
#include "EUTelFFClusterImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelUtility.h"
#include "EUTelGeometryTelescopeGeoDescription.h"


// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/ITupleFactory.h>

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>


#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

//include ROOT<.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TMath.h"


#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <complex> 

using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;

// definition of static 		members mainly used to name histograms
std::string EUTelFitTuple::_FitTupleName  = "EUFit";
std::string EUTelFitTuple::_DistTupleName  = "DistFit";


EUTelFitTuple::EUTelFitTuple() : Processor("EUTelFitTuple")
 {

  // modify processor description
  _description = "Prepare n-tuple with track fit results" ;


  // register steering parameters:
  //       name, description, class-variable, default value

  // input collection first:

  registerInputCollection( LCIO::TRACK,
                           "InputCollectionName" ,
                           "Name of the input Track collection"  ,
                           _inputColName ,
                           std::string("telescopetracks") ) ;

  registerInputCollection( LCIO::TRACKERHIT,
                           "InputHitCollectionName" ,
                           "Name of the input hit collection"  ,
                           _inputHitColName ,
                           std::string("hit") ) ;

/*registerInputCollection( LCIO::ALIGNEDHIT,
                           "InputaligneDHitCollectionName" ,
                           "Name of the input alignedhit collection"  ,
                           _inputalignedHitColName ,
                           std::string("alignedhit") ) ;*/

  // other processor parameters:

  registerProcessorParameter ("MissingValue",
                              "Value used for missing measurements",
                              _missingValue,  static_cast < double > (-100.));


  registerProcessorParameter ("UseManualDUT",
                              "Flag for manual DUT selection",
                              _useManualDUT,  static_cast < bool > (false));

  registerProcessorParameter ("DUTid",
                              "Id of sensor layer which should be used as DUT",
                              _DUTid,  static_cast < int > (20));

  registerProcessorParameter ("DistMax",
                              "Maximum allowed distance between fit and matched DUT hit",
                              _distMax,  static_cast < double > (0.1));
  registerProcessorParameter ("FinderRadius",
							  "Needed for distinguishing the files for diffrent Radii",
							  _FinderRadius, static_cast < double >(500));


  std::vector<float > initAlign;
  initAlign.push_back(0.);
  initAlign.push_back(0.);
  initAlign.push_back(0.);

  registerProcessorParameter ("DUTalignment",
                              "Alignment corrections for DUT: shift in X, Y and rotation around Z",
                              _DUTalign, initAlign);

}


void EUTelFitTuple::init() {

  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;
  _tluTimeStamp = 0 ;


  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    message<ERROR5> ( "The GearMgr is not available, for an unknown reason." );
    exit(-1);
  }

  // Read geometry information from GEAR

  message<MESSAGE5> ( log() << "Reading telescope geometry description from GEAR ") ;

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));


// Take all layers defined in GEAR geometry
  _nTelPlanes = _siPlanesLayerLayout->getNLayers();

// Check for DUT

  if( _siPlanesParameters->getSiPlanesType()==_siPlanesParameters->TelescopeWithDUT )
    {
      _iDUT = _nTelPlanes ;
      _nTelPlanes++;
    }
  else
    _iDUT = -1 ;

// Read position in Z (for sorting)

  _planeSort = new int[_nTelPlanes];
  _planePosition   = new double[_nTelPlanes];

  for(int ipl=0; ipl <  _siPlanesLayerLayout->getNLayers(); ipl++)
    {
      _planePosition[ipl]=_siPlanesLayerLayout->getLayerPositionZ(ipl);
      _planeSort[ipl]=ipl;
    }

  if(_iDUT>0)
    {
      _planePosition[_iDUT]=_siPlanesLayerLayout->getDUTPositionZ();
      _planeSort[_iDUT]=_iDUT;
    }

  // Binary sorting

  bool sorted;
  do{
    sorted=false;
    for(int iz=0; iz<_nTelPlanes-1 ; iz++)
      if(_planePosition[iz]>_planePosition[iz+1])
        {
          double _posZ = _planePosition[iz];
          _planePosition[iz] = _planePosition[iz+1];
          _planePosition[iz+1] = _posZ;

          int _idZ = _planeSort[iz];
          _planeSort[iz] = _planeSort[iz+1];
          _planeSort[iz+1] = _idZ;

          sorted=true;
        }

  }while(sorted);

// Book local geometry arrays

  _planeID = new vector<int>(_nTelPlanes);
  _isActive = new bool[_nTelPlanes];

// Fill remaining layer parameters

	for(int iz=0; iz < _nTelPlanes ; iz++)
    {
      int ipl=_planeSort[iz];

      double resolution;

		if(ipl != _iDUT )
        {
          _planeID->at(ipl) = (_siPlanesLayerLayout->getID(ipl));
          resolution = _siPlanesLayerLayout->getSensitiveResolution(ipl);
        }
		else
        {
          _planeID->at(ipl) = (_siPlanesLayerLayout->getID(ipl));
          resolution = _siPlanesLayerLayout->getDUTSensitiveResolution();
        }

		_isActive[iz] = (resolution > 0);
	}

  // Get new DUT position (after sorting)

  for(int iz=0;iz< _nTelPlanes ; iz++)
    if(_planeSort[iz]==_iDUT)
      {
        _iDUT=iz;
        break;
      }

  // DUT position can be changed by processor parameter

  if(_useManualDUT)
    {
      bool _manualOK=false;

      for(int iz=0; iz < _nTelPlanes ; iz++)
        if(_planeID->at(iz)==_DUTid)
          {
            _iDUT=iz;
            _manualOK=true;
          }

      if(!_manualOK)
        {
          message<ERROR5> ( log() << "Manual DUT flag set, layer not found ID = "
                           << _DUTid
                           << "\n Program will terminate! Correct geometry description!");
          exit(-1);
        }
    }

  // Print out geometry information

  message<MESSAGE5> ( log() << "Telescope configuration with " << _nTelPlanes << " planes" );


  for(int ipl=0; ipl < _nTelPlanes; ipl++)
    {
      stringstream ss ;
      if(ipl == _iDUT)
        ss << "D.U.T.  plane" ;
      else
        if(_isActive[ipl])
          ss << "Active  plane" ;
        else
          ss << "Passive plane" ;

      ss << "  ID = " << _planeID->at(ipl)
         << "  at Z [mm] = " << _planePosition[ipl];

      message<MESSAGE5> ( log() << ss.str() );
    }


  // Allocate arrays for track fitting

  _isMeasured      = new bool[_nTelPlanes];
  _isFitted        = new bool[_nTelPlanes];

  _measuredX     = new double[_nTelPlanes];
  _measuredY     = new double[_nTelPlanes];
  _measuredZ     = new double[_nTelPlanes];
  _measuredQ     = new double[_nTelPlanes];
  _fittedX       = new double[_nTelPlanes];
  _fittedY       = new double[_nTelPlanes];
  _fittedZ       = new double[_nTelPlanes];
  _dutHitcoordinatesLocal = new vector< double* > [_nTelPlanes];
  _dut_tracknumberintern = new vector< int > [_nTelPlanes];
  _distXmin		 = new double[_nTelPlanes];
  _distYmin		 = new double[_nTelPlanes];
  _distRmin		 = new double[_nTelPlanes];

	//Fill distancearry
	for(int ipl=0;ipl<_nTelPlanes;ipl++)
	{ 
	 _distXmin[ipl]= 20;
	 _distYmin[ipl]= 15;
	 _distRmin[ipl]= 25;
	}

	//EUTelUtility for identifying sensor ID from hits
	

	// Book histograms
  	bookHistos();
  	bookHistosNeighbours();

	//initialise EUTelGeometryTelescopeGeoDescription (for coordinate conversion)
	std::string name("EUTelFitTupleGeo.root");
	geo::gGeometry().initializeTGeoDescription(name,false);


}

void EUTelFitTuple::processRunHeader( LCRunHeader* runHeader) {

  auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl ( runHeader ) );
  eutelHeader->addProcessor( type() );

  _nRun++ ;

  // Decode and print out Run Header information - just a check

  _runNr = runHeader->getRunNumber();

  message<MESSAGE5> ( log() << "Processing run header " << _nRun
                     << ", run nr " << _runNr );

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();
  const std::vector<std::string> * subDets = runHeader->getActiveSubdetectors();

  message<MESSAGE5> ( log() << detectorName << " : " << detectorDescription ) ;

  int nDet = subDets->size();

  if(nDet)message<MESSAGE5> ( log() << nDet << " subdetectors defined :" );
  stringstream ss;
  for(int idet=0;idet<nDet;idet++)  message<MESSAGE5> (log()  << idet+1 << " : " << subDets->at(idet) );


}

void EUTelFitTuple::processEvent( LCEvent * event )
{
	EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*> ( event );
	if ( euEvent->getEventType() == kEORE )
	{
		message<DEBUG5> ( "EORE found: nothing else to do." );
		return;
	}

	_nEvt ++ ;
	_evtNr        = event->getEventNumber();
	_tluTimeStamp = static_cast<long int> (event->getTimeStamp());


	LCCollection* col;
	try
	{
		col = event->getCollection( _inputColName ) ;
	} catch (lcio::DataNotAvailableException& e)
	{
    	streamlog_out( DEBUG5 ) << "Not able to get collection " << _inputColName << "from event " << event->getEventNumber() << " in run " << event->getRunNumber() << endl;
    	return;
  	}

	LCCollection* hitcol = NULL;
	bool _DUTok=true;

	try
	{
		hitcol = event->getCollection( _inputHitColName ) ;
	} catch (lcio::DataNotAvailableException& e)
	{
//		message<ERROR5> ( log() << "Not able to get collection "
//								<< _inputHitColName
//								<< "\nfrom event " << event->getEventNumber()_measuredQ[hitPlane]
//								<< " in run " << event->getRunNumber()  );
		_DUTok=false;
	}

	//get plane ID of DUT
	unsigned int dutPlaneID = distance(_planeID->begin(),find(_planeID->begin(), _planeID->end(), _DUTid));

	// clear distance variables
	_dut_neighbor.clear();
	_dut_Q.clear();
	_nearestrackdistRadius.clear();
	_dut_ClusterSize.clear();

	// Loop over tracks in input collections

	int nTrack = col->getNumberOfElements()  ;
	_TOTnTrack= _TOTnTrack + nTrack;
	

	message<DEBUG5> ( log() << "Total of " << nTrack << " tracks in input collection " );
	
	_nHitCol = 0;
	if(_DUTok) _nHitCol = hitcol->getNumberOfElements()  ;

	message<DEBUG5> ( log() << "Total of " << _nHitCol << " hits in input collection " );

	//clear vectors of local hits and distances
	for(int ipl=0; ipl<_nTelPlanes; ipl++)
	{
		_dutHitcoordinatesLocal[ipl].clear();
	}

	//clear Tracknumberintern
	_Tracknumberintern=0;

	for(int itrack=0; itrack< nTrack ; itrack++)
    {
	  _numberalltracks++;
		Track * fittrack = dynamic_cast<Track*>( col->getElementAt(itrack) ) ;

		// Hit list assigned to track

		std::vector<EVENT::TrackerHit*>  trackhits = fittrack->getTrackerHits();

		// Copy hits assign to the track to local table
		// Assign hits to sensor planes


		int nHit =   trackhits.size();

		message<DEBUG5> ( log() << "Track " << itrack << " with " << nHit << " hits, Chi2 = "
                                << fittrack->getChi2() << "/" << fittrack->getNdf());


		// Clear plane tables

		for(int ipl=0;ipl<_nTelPlanes;ipl++)
        {
			_isMeasured[ipl]=false;
			_isFitted[ipl]=false;

			_measuredX[ipl]=_missingValue;
			_measuredY[ipl]=_missingValue;
			_measuredZ[ipl]=_missingValue;
			_measuredQ[ipl]=_missingValue;

			_fittedX[ipl]=_missingValue;
			_fittedY[ipl]=_missingValue;
			_fittedZ[ipl]=_missingValue;

			_dut_tracknumberintern[ipl].clear();
        }
		
		// Clear DUT variables
		double dutX=_missingValue;
		double dutY=_missingValue;
		double _dutQ=_missingValue;
		int _dutClusterSizeX=-10;
		int _dutClusterSizeY=-10;
		

		// setup cellIdDecoder to decode the hit properties
		CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);

		int _numberneighbor=0;

		// Loop over hits and fill hit tables
		for(int ihit=0; ihit< nHit ; ihit++)
        {
			TrackerHit * measHit = trackhits.at(ihit);

			//get hit position
			const double * pos = measHit->getPosition();

			//find plane number of the hit
			int sensorID = Utility::getSensorIDfromHit(measHit);
			unsigned int hitPlane = std::distance(_planeID->begin(),find(_planeID->begin(), _planeID->end(), sensorID));
			if (hitPlane==_planeID->size())
			{
				message<DEBUG5> ( log() << "hit not matched to plane, sensorID=" << sensorID);
				continue;
			}

			if( (hitCellDecoder(measHit)["properties"] & kFittedHit) == 0 )
		    {
				// Measured hits
				_isMeasured[hitPlane]=true;

				_measuredX[hitPlane]=pos[0];
				_measuredY[hitPlane]=pos[1];
				_measuredZ[hitPlane]=pos[2];
				
				// get local coordinates
				double* posLocal = new double[3];
				geo::gGeometry().master2Localtwo(sensorID,pos,posLocal);
				_dutHitcoordinatesLocal[hitPlane].push_back( posLocal );

				//save extra hit information if DUT
				if (hitPlane==dutPlaneID)
				{
					dutX = posLocal[0];
					dutY = posLocal[1];
				}

				// Get cluster charge

				EVENT::LCObjectVec rawdata = measHit->getRawHits();
				if(rawdata.size()>0 && rawdata.at(0)!=NULL )
				{
				EUTelVirtualCluster* cluster = new EUTelSparseClusterImpl<EUTelGenericSparsePixel>(static_cast<TrackerDataImpl*>(rawdata.at(0)));
				_measuredQ[hitPlane]=cluster->getTotalCharge();

					//save extra cluster information if DUT
					if (hitPlane==dutPlaneID)
					{
					_dutQ = _measuredQ[hitPlane];
					cluster->getClusterSize(_dutClusterSizeX, _dutClusterSizeY);
					message<DEBUG5> ( log() << "hitplane in 1 loop " << hitPlane << " at  X = "
		                                << pos[0] << ", Y = " << pos[1] << ", Q = " << _dutQ << " clustersizecx " <<_dutClusterSizeX << " clustersizey " 			
										<<_dutClusterSizeY );

						IMPL::TrackerDataImpl* clusterContent = cluster->trackerData();
						std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> >
						apixData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(clusterContent));
						EUTelGenericSparsePixel apixPixel;
						_nPixHits=0;
						for( unsigned int iHit = 0; iHit < apixData->size(); iHit++ )
						{
							//apixData->getSparsePixelAt( iHit, &apixPixel);
							_nPixHits++;
							//_p_col->push_back( apixPixel.getXCoord() );
							//_p_row->push_back( apixPixel.getYCoord() );
							//_TOT+=static_cast< int >(apixPixel.getSignal());
							//_p_lv1->push_back( static_cast< int >(apixPixel.getTime()) );			
						}




					}
				}
				message<DEBUG5> ( log() << "Measured hit in plane " << hitPlane << " at  X = "
		                                << pos[0] << ", Y = " << pos[1] << ", Q = " << _measuredQ[hitPlane] );

				message<DEBUG5> ( log() << "hitplane " << hitPlane << " at  X = "
		                                << pos[0] << ", Y = " << pos[1] << ", Q = " << _measuredQ[hitPlane]  << " clustersizecx " <<_dutClusterSizeX << " clustersizey " 			
										<<_dutClusterSizeY );
				
			}
			else
		    {
				// Fitted hits

				_isFitted[hitPlane]=true;

				_fittedX[hitPlane]=pos[0];
				_fittedY[hitPlane]=pos[1];
				_fittedZ[hitPlane]=pos[2];

				message<DEBUG5> ( log() << "Fitted  hit  in plane " << hitPlane << " at  X = "
		                                << pos[0] << ", Y = " << pos[1] );
		    }
        

		//find neighbours from tracks from hitcollection
		
			//_measuredCharge  =  new double [_nTelPlanes]; 
			//_measuredCharge[_nTelPlanes] = 0;

			//for (int sensorID=0;sensorID<_nTelPlanes;sensorID++)
			//TrackerHit * _FitHitwithneigbours = 0;
			
			for(int ihitL=0; ihitL< _nHitCol ; ihitL++)
		    {	double distmininX=250; //in mircons
			  	double distmininY=50; //in mircons
			  
			 

				TrackerHit * measHit2 = dynamic_cast<TrackerHit*>( hitcol->getElementAt(ihitL) );
				int sensorID2 = Utility::getSensorIDfromHit(measHit2);
				if (sensorID2!=_DUTid)//select hits on plane unequal the dutplane of hitK in tracks
				{
					
					const double * pos = measHit2->getPosition();
					//calculate distance between fitted track hit and this hit
				
						//hit in global coordinates, simply subtract
						double distXglobal = abs(pos[0]-_measuredX[sensorID2])*1000;
						if(distXglobal<distmininX && distXglobal!=0 )
						{	
							double distYglobal = abs(pos[1]-_measuredY[sensorID2]*1000);
							if(distYglobal<distmininY && distYglobal!=0)
							{	
								_numberneighbor++;
								message<DEBUG5> ( log() << "Sensoridfor neigbhorcounter" << sensorID2 << " from track "<< _TOTnTrack );

							}
						}
				}
		     
			}
		
		}// End of loop over hits in track 

		if (_numberneighbor!=0)
		{ 	/*_measuredQ[hitPlane]=0.;
			EVENT::LCObjectVec rawdata =  measHit->getRawHits();

			if(rawdata.size()>0 && rawdata.at(0)!=NULL )
			{
			
				EUTelVirtualCluster * cluster = new EUTelSparseClusterImpl<EUTelGenericSparsePixel>(static_cast<TrackerDataImpl*>(rawdata.at(0)));
				_measuredQ[hitPlane]=cluster->getTotalCharge();
				_dutCharge=_measuredQ[_DUTid];					
				cluster->getClusterSize(dutClusterSizeXdist,dutClusterSizeYdist);
				//message<MESSAGE5> ( log() << "Q= " << _measuredQ[_DUTid]<< " ClustersizeX "
		                           // << dutClusterSizeX << ", clustersizeY " << dutClusterSizeY << "with n neighbours" << _numberneighbor);
		

				IMPL::TrackerDataImpl* clusterContent = cluster->trackerData();
				std::auto_ptr<EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel> >
				apixData( new EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(clusterContent));
				EUTelGenericSparsePixel apixPixel;
				_TOT=0;
				for( unsigned int iHit = 0; iHit < apixData->size(); iHit++ )
				{
					apixData->getSparsePixelAt( iHit, &apixPixel);
					_nPixHits++;
					//_p_col->push_back( apixPixel.getXCoord() );
					//_p_row->push_back( apixPixel.getYCoord() );
					_TOT+=static_cast< int >(apixPixel.getSignal());
					//_p_lv1->push_back( static_cast< int >(apixPixel.getTime()) );

			
				}
			
			}*/
			message<DEBUG5> ( log() << "Q= " << _dutQ <<", Clustersize: " 
		                            << _nPixHits << ", clustersizeY: " << _dutClusterSizeY  << ", with " << _numberneighbor << " neighbours, tracknumber: " << _TOTnTrack);
					_Tracknumberintern++;
					_dut_Q.push_back(_dutQ);
					_dut_ClusterSize.push_back(_nPixHits);
					_dut_ClusterSizeY.push_back(_dutClusterSizeY);
					_dut_neighbor.push_back(_numberneighbor);
					_dut_tracknumber.push_back(_numberalltracks);
					/*for (int ipl=0; ipl<_nTelPlanes; ipl++)
					{
						if (_isMeasured[ipl])
						{
							_dut_tracknumberintern[ipl].push_back(_dutHitcoordinatesLocal[ipl].size()-1);
						} else
						{
							_dut_tracknumberintern[ipl].push_back(-1);
						}
					}*/
		}

			// Try to get cluster charge and size of the DUT from the track with neighbors
		 

		// Fill n-tuple

		int icol=0;
		_FitTuple->fill(icol++,_nEvt);
		_FitTuple->fill(icol++,_runNr);
		_FitTuple->fill(icol++,_evtNr);
		_FitTuple->fill(icol++,_tluTimeStamp);
		_FitTuple->fill(icol++,nTrack);
		_FitTuple->fill(icol++,fittrack->getNdf());
		_FitTuple->fill(icol++,fittrack->getChi2());

		for(int ipl=0; ipl<_nTelPlanes;ipl++)
        {
			_FitTuple->fill(icol++,_measuredX[ipl]);
			_FitTuple->fill(icol++,_measuredY[ipl]);
			_FitTuple->fill(icol++,_measuredZ[ipl]);
			_FitTuple->fill(icol++,_measuredQ[ipl]);
			_FitTuple->fill(icol++,_fittedX[ipl]);
			_FitTuple->fill(icol++,_fittedY[ipl]);
			_FitTuple->fill(icol++,_fittedZ[ipl]);
        }
		_FitTuple->fill(icol++,dutX);
		_FitTuple->fill(icol++,dutY);
		_FitTuple->fill(icol++,_dutQ);
		_FitTuple->fill(icol++,_dutClusterSizeX);
		_FitTuple->fill(icol++,_dutClusterSizeY);
		//_FitTuple->fill(icol++,dutZ);_FitTuple->fill(icol++,dutR);
		_FitTuple->addRow();
    }// End of loop over tracks

	//fill histograms for distance measurement first run over planes
	
	unsigned int nHits=0;
	
		
	// fill nTuple with distances
	double distX;
	double distY;
	int ijcolX=0;
	int ijcolY=1;
	double distRsquarttot;
	double distRsquart;
	double distRminsquart;
	//_nearestrackdistRadius.resize(_Tracknumberintern);
	/*for(int m=0; m < _Tracknumberintern;m++)
	{
		_nearestrackdistRadius.at(m) = 100;
	}*/
	//message<DEBUG5> ( log() << "_nearestrackdistRadius.size(): " << _nearestrackdistRadius.size() << ", _Tracknumberintern: " << _Tracknumberintern <<", _dut_tracknumberintern[0].size(): "<< _dut_tracknumberintern[0].size());
    for(int ipl=0; ipl <_nTelPlanes; ipl++)
	{ 
		nHits = _dutHitcoordinatesLocal[ipl].size();
		//if( _dutHitcoordinatesLocal[ipl].size() == 0 ) continue;

		message<DEBUG4> ( log() << "distX: i: " << ipl << " " << _dutHitcoordinatesLocal[ipl].size() );
		distRsquarttot=_distXmin[ipl]*_distXmin[ipl]+_distYmin[ipl]*_distYmin[ipl];
	   	for(unsigned int i=0; i< nHits; i++ )
		{	
			double distRnear;
			
			distRminsquart= 100;
			for(unsigned int j=i+1; j< nHits;j++)
			{
			 

				message<DEBUG5> ( log() << "distX: " << (_dutHitcoordinatesLocal[ipl][i][0]-_dutHitcoordinatesLocal[ipl][j][0])<<"_for Plane_" <<ipl);
				distX = abs(_dutHitcoordinatesLocal[ipl][i][0]-_dutHitcoordinatesLocal[ipl][j][0]);
				distY = abs(_dutHitcoordinatesLocal[ipl][i][1]-_dutHitcoordinatesLocal[ipl][j][1]);
				_DistTuple[ipl]->fill(ijcolX,distX);
				_DistTuple[ipl]->fill(ijcolY,distY);
				
				_DistTuple[ipl]->addRow();	
				//get minimal distance!=0 between two tracks 

				distRsquart=distX*distX+distY*distY;
				if (distRsquart < distRminsquart)
				{
					distRminsquart = distRsquart;
				}
				if (distRsquart<distRsquarttot && distRsquart !=0 )
				{ if (distX>1e-3)//for excluding to low values which come from other reasons
					{			 
					_distXmin[ipl] = distX;
					message<DEBUG5> ( log() << "distXmin: " << _distXmin[ipl] << "at Plane" << ipl <<"Evtnumber"<<_nEvt);
					_distYmin[ipl] = distY;
					message<DEBUG5> ( log() << "distYmin: " << _distYmin[ipl] << "at Plane" << ipl <<"Evtnumber"<<_nEvt);
					_distRmin[ipl] = sqrt(distRsquart);
					message<DEBUG5> ( log() << "distRmin: " << _distRmin[ipl] << "at Plane" << ipl <<"Evtnumber"<<_nEvt);
					distRsquarttot = distRsquart;
					}
				}
				/*if(distX<_distXmin[ipl] && distX !=0)
				{ if (distX>1e-5){			 //for excluding to low values which come from other reasons
				  _distXmin[ipl] = distX;
				message<MESSAGE5> ( log() << "distXmin: " << _distXmin[ipl] << "at Plane" << ipl <<"Evtnumber"<<_nEvt);}
				}
				if(distY<_distYmin[ipl] && distY !=0)
				{ if (distY>1e-5){			//for excluding to low values which come from other reasons
				  _distYmin[ipl] = distY;
				message<DEBUG5> ( log() << "distYmin: " << _distYmin[ipl] << "at Plane" << ipl <<"Evtnumber"<<_nEvt);}
				}	*/					
			}
			/*unsigned int index = distance(_dut_tracknumberintern[ipl].begin(),find(_dut_tracknumberintern[ipl].begin(), _dut_tracknumberintern[ipl].end(), i));
			if (index<_dut_tracknumberintern[ipl].size())
			{
			 	distRnear=sqrt(distRminsquart);
				if(distRnear < _nearestrackdistRadius.at(index) && distRnear!=0 )
				{
					_nearestrackdistRadius.at(index) = distRnear;
					message<MESSAGE5> ( log() << "nearestTrack at Radius: " << _nearestrackdistRadius.at(index) << " at Index " << index <<" Tracknumber "<< _dut_tracknumber.at(index));
				}
			 	
			}*/
			
		}	
		ijcolX=ijcolX+2;
		ijcolY=ijcolY+2;
	
	}

	treeneighbourw();


	
	return;
}



void EUTelFitTuple::check( LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelFitTuple::end(){

  //   std::cout << "EUTelFitTuple::end()  " << name()
  //        << " processed " << _nEvt << " events in " << _nRun << " runs "
  //        << std::endl ;


  message<MESSAGE5> ( log() << "fit N-tuple with "
                     << _FitTuple->rows() << " rows created" );
  message<MESSAGE5> ( log() << "dist0 N-tuple with "
                     << _DistTuple[0]->rows() << " rows created" );
  message<MESSAGE5> ( log() << "dist1 N-tuple with "
                     << _DistTuple[1]->rows() << " rows created" );
  message<MESSAGE5> ( log() << "dist2 N-tuple with "
                     << _DistTuple[2]->rows() << " rows created" );
  message<MESSAGE5> ( log() << "dist3 N-tuple with "
                     << _DistTuple[3]->rows() << " rows created" );
  message<MESSAGE5> ( log() << "dist4 N-tuple with "
                     << _DistTuple[4]->rows() << " rows created" );
  message<MESSAGE5> ( log() << "dist5 N-tuple with "
                     << _DistTuple[5]->rows() << " rows created" );
  message<MESSAGE5> ( log() << "dist6 N-tuple with "
                     << _DistTuple[6]->rows() << " rows created" );
  
  for(int ipl=0; ipl<_nTelPlanes;ipl++)
  {
	message<MESSAGE5> ( log() << "distXminEND: " << _distXmin[ipl] << "at Plane_" << ipl <<"_FinderRadius_"<<_FinderRadius);
	message<MESSAGE5> ( log() << "distYminEND: " << _distYmin[ipl] << "at Plane_" << ipl <<"_FinderRadius_"<<_FinderRadius);
	message<MESSAGE5> ( log() << "distRminEND: " << _distRmin[ipl] << "at Plane_" << ipl <<"_FinderRadius_"<<_FinderRadius);
  /*if (_distXmin[ipl]<_distYmin[ipl])
	{
 	 message<MESSAGE5> ( log() << "distminEND: " << _distXmin[ipl] << "at Plane_" << ipl <<"_FinderRadius_"<<_FinderRadius);
	} 
	else
	{
  	message<MESSAGE5> ( log() << "distminEND: " <<_distYmin[ipl] << "at Plane_" << ipl <<"_FinderRadius_"<<_FinderRadius);
	}*/
  }

 // Write Tree with min distances in all Planes, resolved tracks in comparison to the FinderRadius

	//call tree for histos
	//
	//treedistr();
 treedistw();
	
  // Clean memory

  delete [] _planeSort ;
  delete [] _planePosition ;
  delete _planeID ;
  delete [] _isActive ;

  delete [] _isMeasured ;
  delete [] _isFitted ;
  delete [] _measuredX  ;
  delete [] _measuredY  ;
  delete [] _measuredZ  ;
  delete [] _measuredQ  ;
  delete [] _fittedX ;
  delete [] _fittedY ;
  delete [] _fittedZ ;

  delete [] _dutHitcoordinatesLocal;
}



void EUTelFitTuple::treeneighbourw()
{
	//create a Tree file tree_neigh.root
   	message<DEBUG5> ( log() << "filling the treeneighbourw()" );
   
   //create the file, the Tree and a few branches
	//fill tree 
	int icol;


	for(int i=0;i<_Tracknumberintern;i++)
	{
      // 	message<MESSAGE5> ( log() << "i: " << i << " of " << _Tracknumberintern << " with " << _dut_neighbor.at(i) << " neighbours" << " Q:" << _dut_Q.at(i) << ", clustersize: " << _dut_ClusterSize.at(i) );
		icol=0;
    	_NeighbourTuple->fill(icol++, static_cast<int>( _dut_neighbor.at(i) ));
    	_NeighbourTuple->fill(icol++, static_cast<int> (_dut_Q.at(i)) );	
		_NeighbourTuple->fill(icol++, static_cast<int>( _dut_ClusterSize.at(i)));
    //	_NeighbourTuple->fill(icol++, static_cast<double>( _nearestrackdistRadius.at(i)));
	
	_NeighbourTuple->addRow();
	

//	dutclustersizeX=_dut_ClusterSizeX.at(i);
//	dutclustersizeY=_dut_ClusterSizeY.at(i);
//	dutclusterTOT=_dut_Q.at(i);
//	nearestrackdistRadius=_nearestrackdistRadius.at(i);
    	
	}	

	


}


void EUTelFitTuple::treedistw()
{
   //create a Tree file tree_dist.root
   
   //create the file, the Tree and a few branches
   TFile* tdist= new TFile("treedist2.root","update");
   //TTree treedist("treedist","Tree with distances");
   TTree* treedist;
	treedist = (TTree*)tdist->Get("treedist");
   
  
 
	Double_t averageXmin,averageYmin,distXmin[7], distYmin[7], distRmin[7];			
	Int_t ToTnTrack, runnumber, FitterRadius;

	treedist->SetBranchAddress("averageXmin",&averageXmin);		//,"averageXmin/D");
	treedist->SetBranchAddress("distRmin",&distRmin);		//,"distRmin[7]/D");
	treedist->SetBranchAddress("averageYmin",&averageYmin);		//,"averageYmin/D");
	treedist->SetBranchAddress("distXmin",&distXmin);		//,"distXmin[7]/D");
	treedist->SetBranchAddress("distYmin",&distYmin);		//,"distYmin[7]/D");				
	treedist->SetBranchAddress("runnumber",&runnumber);		//,"runnumber/I");
   treedist->SetBranchAddress("ToTnTrack",&ToTnTrack);		//,"ToTnTrack/I");
   treedist->SetBranchAddress("FitterRadius",&FitterRadius);		//,"FitterRadius/I");
   
   //fill the tree
  
	/*for(int ipl=0; ipl<_nTelPlanes;ipl++){
	distXmin[ipl] = _distXmin[ipl];
	distYmin[ipl] = _distYmin[ipl];
    treedist->Fill();
	}*/
   distXmin[0]= _distXmin[0];distXmin[1]= _distXmin[1];distXmin[2]= _distXmin[2];distXmin[3]= _distXmin[3];distXmin[4]=_distXmin[4];distXmin[5]= _distXmin[5];distXmin[6]= _distXmin[6];
   distRmin[0]= _distRmin[0];distRmin[1]= _distRmin[1];distRmin[2]= _distRmin[2];distRmin[3]= _distRmin[3];distRmin[4]=_distRmin[4];distRmin[5]= _distRmin[5];distRmin[6]= _distRmin[6];
   distYmin[0]= _distYmin[0];distYmin[1]= _distYmin[1];distYmin[2]= _distYmin[2];distYmin[3]= _distYmin[3];distYmin[4]=_distYmin[4];distYmin[5]= _distYmin[5];distYmin[6]= _distYmin[6];
	double avXmin=0, avYmin=0;	
	for(int ipl=0; ipl<_nTelPlanes;ipl++){
	avXmin= avXmin +_distXmin[ipl];
	avYmin= avYmin +_distYmin[ipl];
	}
	averageXmin=(avXmin/_nTelPlanes); 
	averageYmin=(avYmin/_nTelPlanes);
    runnumber = _nRun ;
	FitterRadius = _FinderRadius;
	ToTnTrack = _TOTnTrack;
	treedist->Fill();
  
  
  //save the Tree header. The file will be automatically closed
  //when going out of the function scope

  
  	treedist->Write();
	tdist->Close();
  //treedist.SetFile(0);
  //tdist->Close();
  //delete tdist;
  //treedist.Print();
}


/*void EUTelFitTuple::treedistr()
{
   //read the Tree generated by treedistw and fill histograms
   
   //note that we use "new" to create the TFile and TTree objects !
   //because we want to keep these objects alive when we leave this function.
   
   TFile *tdist = new TFile("treedist.root");
   TTree *treedist = (TTree*)tdist->Get("treedist");

   
   Int_t runnumber, FitterRadius;
   Double_t averageXmin,averageYmin,distXmin[7],distYmin[7];
   Int_t ToTnTrack;

   treedist->SetBranchAddress("FitterRadius",&FitterRadius);
   treedist->SetBranchAddress("ToTnTrack",&ToTnTrack);
   treedist->SetBranchAddress("averageXmin",&averageXmin);
   treedist->SetBranchAddress("averageYmin",&averageYmin);
   treedist->SetBranchAddress("distXmin",&distXmin); 
   treedist->SetBranchAddress("distYmin",&distYmin);
   treedist->SetBranchAddress("runnumber",&runnumber);

   
   //create two histograms
   TH1F *hTrackn   = new TH1F("Trackn","Number of Tracks vs FinderRadius",100,0,10000);
   TH2F *hTracknvsRadius = new TH2F("TracknvsRadius","Track vs Radius",50,0,1000,100,0,10000);
   TH2F *hdistXavaragevsRadius = new TH2F("distXavaragevsRadius","DistanceXmin Avarage vs Radius",50,0,1000,10000,1e-4,1e-2);
   TH2F *hdistYavaragevsRadius = new TH2F("distYavaragevsRadius","DistanceYmin Avarage vs Radius",50,0,1000,10000,1e-4,1e-2);
	TH2F *hdistX0vsRadius = new TH2F("distX0vsRadius","DistanceXmin Plane 0 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistX1vsRadius = new TH2F("distX1vsRadius","DistanceXmin Plane 1 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistX2vsRadius = new TH2F("distX2vsRadius","DistanceXmin Plane 1 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistX3vsRadius = new TH2F("distX3vsRadius","DistanceXmin Plane 3 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistX4vsRadius = new TH2F("distX4vsRadius","DistanceXmin Plane 4 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistX5vsRadius = new TH2F("distX5vsRadius","DistanceXmin Plane 5 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistX6vsRadius = new TH2F("distX6vsRadius","DistanceXmin Plane 6 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistY0vsRadius = new TH2F("distY0vsRadius","DistanceYmin Plane 0 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistY1vsRadius = new TH2F("distY1vsRadius","DistanceYmin Plane 1 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistY2vsRadius = new TH2F("distY2vsRadius","DistanceXmin Plane 1 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistY3vsRadius = new TH2F("distY3vsRadius","DistanceYmin Plane 3 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistY4vsRadius = new TH2F("distY4vsRadius","DistanceYmin Plane 4 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistY5vsRadius = new TH2F("distY5vsRadius","DistanceYmin Plane 5 vs Radius",50,0,1000,10000,1e-5,1e-2);
	TH2F *hdistY6vsRadius = new TH2F("distY6vsRadius","DistanceYmin Plane 6 vs Radius",50,0,1000,10000,1e-5,1e-2);
	//TH2F *hdistXTOTALvsRadius = new TH2F("distXTotalvsRadius","DistanceXmin for all 7 Planes vs Radius",50,0,700,10000,1e-5,1e-3);
	//TH2F *hdistYTOTALvsRadius = new TH2F("distYTotalvsRadius","DistanceYmin for all 7 Planes vs Radius",50,0,700,10000,1e-5,1e-3); 
  
   //read all entries and fill the histograms
   Long64_t nentries = treedist->GetEntries();
   for (Long64_t i=0;i<nentries;i++) {
     treedist->GetEntry(i);
	hTrackn->Fill(ToTnTrack);
	 hTracknvsRadius->Fill(FitterRadius,ToTnTrack);
	for (int j=0; j<=6; j++)
	{hdistYTOTALvsRadius->Fill(FitterRadius,distYmin[j]);
	hdistYTOTALvsRadius->SetMarkerColor(j);
	}
	 hdistXavaragevsRadius->Fill(FitterRadius,averageXmin);
	 hdistYavaragevsRadius->Fill(FitterRadius,averageYmin);
 	hdistX0vsRadius->Fill(FitterRadius,distXmin[0]);
 	hdistX1vsRadius->Fill(FitterRadius,distXmin[1]);
	hdistX2vsRadius->Fill(FitterRadius,distXmin[2]);	
	hdistX3vsRadius->Fill(FitterRadius,distXmin[3]);
	hdistX4vsRadius->Fill(FitterRadius,distXmin[4]);
	hdistX5vsRadius->Fill(FitterRadius,distXmin[5]);
	hdistX6vsRadius->Fill(FitterRadius,distXmin[6]);	
	hdistX0vsRadius->Fill(FitterRadius,distYmin[0]);
	hdistY1vsRadius->Fill(FitterRadius,distYmin[1]);
	hdistY2vsRadius->Fill(FitterRadius,distYmin[2]);
	hdistY3vsRadius->Fill(FitterRadius,distYmin[3]);
	hdistY4vsRadius->Fill(FitterRadius,distYmin[4]);
	hdistY5vsRadius->Fill(FitterRadius,distYmin[5]);
	hdistY6vsRadius->Fill(FitterRadius,distYmin[6]); 
	
    

  }

  //hTracknvsRadius->SaveAs("tracknumbervsRadius.pdf");

  TCanvas c_avarageXmin;
	c_avarageXmin.SetLogy();
	c_avarageXmin.SetFillColor(42);
	//gPad->SetFillColor(12);
	hdistXavaragevsRadius->Draw("AP");
	c_avarageXmin.SaveAs("AvarageXmin.pdf");

  TCanvas c_avarageYmin;
	c_avarageYmin.SetLogy();
	hdistXavaragevsRadius->SetMarkerStyle(8);
	hdistYavaragevsRadius->SetMarkerColor(9);
	hdistYavaragevsRadius->Draw();
	c_avarageYmin.SaveAs("AvarageYmin.pdf");

TCanvas c_AllXmin;
	c_AllXmin.SetLogy();
	hdistX0vsRadius->SetMarkerStyle(8);
 	hdistX1vsRadius->SetMarkerStyle(8);
	hdistX2vsRadius->SetMarkerStyle(8);	
	hdistX3vsRadius->SetMarkerStyle(8);
	hdistX4vsRadius->SetMarkerStyle(8);
	hdistX5vsRadius->SetMarkerStyle(8);
	hdistX6vsRadius->SetMarkerStyle(8);
	hdistX0vsRadius->SetMarkerColor(1);
 	hdistX1vsRadius->SetMarkerColor(2);
	hdistX2vsRadius->SetMarkerColor(9);	
	hdistX3vsRadius->SetMarkerColor(4);
	hdistX4vsRadius->SetMarkerColor(5);
	hdistX5vsRadius->SetMarkerColor(6);
	hdistX6vsRadius->SetMarkerColor(7);
	hdistX0vsRadius->Draw();
 	hdistX1vsRadius->Draw("same");
	hdistX2vsRadius->Draw("same");	
	hdistX3vsRadius->Draw("same");
	hdistX4vsRadius->Draw("same");
	hdistX5vsRadius->Draw("same");
	hdistX6vsRadius->Draw("same");
	c_avarageYmin.SaveAs("AllX.pdf");

TCanvas c_AllYmin;
	c_AllYmin.SetLogy();
	hdistY0vsRadius->SetMarkerStyle(8);
 	hdistY1vsRadius->SetMarkerStyle(8);
	hdistY2vsRadius->SetMarkerStyle(8);
	hdistY3vsRadius->SetMarkerStyle(8);
	hdistY4vsRadius->SetMarkerStyle(8);
	hdistY5vsRadius->SetMarkerStyle(8);
	hdistY6vsRadius->SetMarkerStyle(8);
	hdistY0vsRadius->SetMarkerColor(1);
 	hdistY1vsRadius->SetMarkerColor(2);
	hdistY2vsRadius->SetMarkerColor(9);
	hdistY3vsRadius->SetMarkerColor(4);
	hdistY4vsRadius->SetMarkerColor(5);
	hdistY5vsRadius->SetMarkerColor(6);
	hdistY6vsRadius->SetMarkerColor(7);
	hdistY0vsRadius->Draw("same");
 	hdistY1vsRadius->Draw("same");
	hdistY2vsRadius->Draw("same");	
	hdistY3vsRadius->Draw("same");
	hdistY4vsRadius->Draw("same");
	hdistY5vsRadius->Draw("same");
	hdistY6vsRadius->Draw("same");
	c_AllYmin.SaveAs("AllY.pdf");



 


    TCanvas c_TrackRadius;
    c_TrackRadius.cd();
	hTracknvsRadius->Draw();
	c_TrackRadius.SaveAs("tracknumbervsRadius.pdf");

    TCanvas c_tot;
	c_tot.cd();
	hTrackn->Draw();
	c_tot.SaveAs("TOD.pdf");

	TCanvas c_X0min;
	hdistX0vsRadius->Draw();
	c_X0min.SaveAs("distX0vsRadius.pdf");

	TCanvas c_X1min;
	hdistX1vsRadius->Draw();
	c_X1min.SaveAs("distX1vsRadius.pdf");

}   
*/
void EUTelFitTuple::bookHistos()
{


  message<MESSAGE5> ( log() << "Booking fit n-tuple \n" );

  std::vector<std::string> _columnNames;
  std::vector<std::string> _columnType;

  _columnNames.push_back("Event");
  _columnType.push_back("int");

  _columnNames.push_back("RunNr");
  _columnType.push_back("int");

  _columnNames.push_back("EvtNr");
  _columnType.push_back("int");

  _columnNames.push_back("TLUtime");
  _columnType.push_back("long int");

  _columnNames.push_back("Track");
  _columnType.push_back("int");

  _columnNames.push_back("Ndf");
  _columnType.push_back("int");

  _columnNames.push_back("Chi2");
  _columnType.push_back("float");

  const char * _varName[] = { "measX", "measY" , "measZ", "measQ", "fitX", "fitY", "fitZ" };

  for(int ipl=0; ipl<_nTelPlanes;ipl++)
  {
    for(int ivar=0; ivar<7;ivar++)
      {
        stringstream ss;
        ss << _varName[ivar] << "_" <<  _planeID->at(ipl);
        _columnNames.push_back(ss.str());
        _columnType.push_back("double");
      }
  }

  // DUT variables

  _columnNames.push_back("dutX");
  _columnType.push_back("double");

  _columnNames.push_back("dutY");
  _columnType.push_back("double");

  //_columnNames.push_back("dutZ");
  //_columnType.push_back("double");

  //_columnNames.push_back("dutR");
  //_columnType.push_back("double");

  _columnNames.push_back("dutQ");
  _columnType.push_back("double");

  _columnNames.push_back("dutClusterSizeX");
  _columnType.push_back("int");

  _columnNames.push_back("dutClusterSizeY");
  _columnType.push_back("int");

  _FitTuple=AIDAProcessor::tupleFactory(this)->create(_FitTupleName, _FitTupleName, _columnNames, _columnType, "");

  //Tuple for the distances from hit to hit in x and y
  message<MESSAGE5> ( log() << "Booking dist n-tuple \n" );
  //AIDA::ITuple _DistTuple[_nTelPlanes];
  //_DistTuple.resize(_nTelPlanes);

  std::vector<std::string> _distColumnNames;
  std::vector<std::string> _distColumnType;

	for(int ipl=0; ipl<_nTelPlanes;ipl++)
	{
		stringstream ssx;
        ssx << "xDist_" << _planeID->at(ipl);
		_distColumnNames.push_back(ssx.str());
		_distColumnType.push_back("double");

		stringstream ssy;
        ssy << "yDist_" << _planeID->at(ipl);
		_distColumnNames.push_back(ssy.str());
		_distColumnType.push_back("double");

		stringstream sName;
		sName << _DistTupleName << "_" << _planeID->at(ipl);
		_DistTuple.push_back( AIDAProcessor::tupleFactory(this)->create(sName.str(), sName.str(), _distColumnNames, _distColumnType, ""));
	}
	

    message<DEBUG5> ( log() << "Booking completed \n\n");
    message<MESSAGE5> ( log() << "Total of " << _TOTnTrack << " tracks" );
    return;
}

void EUTelFitTuple::bookHistosNeighbours()
{


  message<MESSAGE5> ( log() << "Booking fit n-tuple \n" );

  std::vector<std::string> _columnNames;
  std::vector<std::string> _columnType;



  // DUT variables

  _columnNames.push_back("nneighbours");
  _columnType.push_back("int");

  _columnNames.push_back("dutclusterTOT");
  _columnType.push_back("int");

  _columnNames.push_back("dutClusterSize");
  _columnType.push_back("int");

  //_columnNames.push_back("nearestrackdistRadius");
 // _columnType.push_back("double");

  
  _NeighbourTuple=AIDAProcessor::tupleFactory(this)->create("treeneigh", "treeneigh", _columnNames, _columnType, "");

  //Tuple for the distances from hit to hit in x and y
  message<MESSAGE5> ( log() << "Booking dist n-tuple \n" );
  //AIDA::ITuple _DistTuple[_nTelPlanes];
  //_DistTuple.resize(_nTelPlanes);

	
	//book Treeneighbours
/*
	tneigh= new TFile("treeneigh.root","recreate");
	treeneigh = new TTree("treeneigh","Tree with neighbours");

	treeneigh->Branch("nneighbours",&nneighbours,"nneighbours/I");
	treeneigh->Branch("dutclusterTOT",&dutclusterTOT,"dutclusterTOT/I");
	treeneigh->Branch("dutclustersizeY",&dutclustersizeY,"dutclustersizeY/I");
	treeneigh->Branch("dutclustersizeX",&dutclustersizeX,"dutclustersizeX/I");
	treeneigh->Branch("nearestrackdistRadius",&nearestrackdistRadius,"nearestrackdistRadius/D");
	treeneigh->Branch("Tracknumber",&Tracknumber,"Tracknumber/I");
*/

    message<DEBUG5> ( log() << "Booking completed \n\n");
    message<MESSAGE5> ( log() << "Total of " << _TOTnTrack << " tracks" );
    return;
}






#endif // GEAR && AIDA


// Author: A.F.Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Revised by Simon De Ridder
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
#include "EUTelFFClusterImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelHistogramManager.h"
#include "EUTelExceptions.h"
#include "EUTelUtility.h"


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

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

using namespace std;
using namespace lcio ;
using namespace marlin ;
using namespace eutelescope;

// definition of static 		members mainly used to name histograms
std::string EUTelFitTuple::_FitTupleName  = "EUFit";


EUTelFitTuple::EUTelFitTuple() : Processor("EUTelFitTuple") {

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

  // other processor parameters:

  registerProcessorParameter ("MissingValue",
                              "Value used for missing measurements",
                              _missingValue,  static_cast < double > (-100.));


  registerProcessorParameter ("UseManualDUT",
                              "Flag for manual DUT selection",
                              _useManualDUT,  static_cast < bool > (false));

  registerProcessorParameter ("DUTid",
                              "Id of sensor layer which should be used as DUT",
                              _DUTid,  static_cast < int > (0));

  registerProcessorParameter ("DistMax",
                              "Maximum allowed distance between fit and matched DUT hit",
                              _distMax,  static_cast < double > (0.1));


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

  _planeID         = new int[_nTelPlanes];
  _isActive        = new bool[_nTelPlanes];

// Fill remaining layer parameters

	for(int iz=0; iz < _nTelPlanes ; iz++)
    {
      int ipl=_planeSort[iz];

      double resolution;

		if(ipl != _iDUT )
        {
          _planeID[iz]=_siPlanesLayerLayout->getID(ipl);
          resolution = _siPlanesLayerLayout->getSensitiveResolution(ipl);
        }
		else
        {
          _planeID[iz]=_siPlanesLayerLayout->getDUTID();
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
        if(_planeID[iz]==_DUTid)
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

      ss << "  ID = " << _planeID[ipl]
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

	//EUTelUtility for identifying sensor ID from hits
	

	// Book histograms
  	bookHistos();
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
//								<< "\nfrom event " << event->getEventNumber()
//								<< " in run " << event->getRunNumber()  );
		_DUTok=false;
	}


	// Loop over tracks in input collections

	int nTrack = col->getNumberOfElements()  ;

	message<DEBUG5> ( log() << "Total of " << nTrack << " tracks in input collection " );

	int nHitCol = 0;
	if(_DUTok) nHitCol = hitcol->getNumberOfElements()  ;

	message<DEBUG5> ( log() << "Total of " << nHitCol << " hits in input collection " );


	for(int itrack=0; itrack< nTrack ; itrack++)
    {
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
        }
		
		// Clear DUT variables
		double dutX=_missingValue;
		double dutY=_missingValue;
		double dutZ=_missingValue;
		double dutAlpha=_missingValue;
		double dutBeta=_missingValue;
		double dutGamma=_missingValue;
		double dutR=_missingValue;
		double dutQ=_missingValue;
		int dutClusterSizeX=-10;
		int dutClusterSizeY=-10;

		// setup cellIdDecoder to decode the hit properties
		CellIDDecoder<TrackerHit>  hitCellDecoder(EUTELESCOPE::HITENCODING);

		// Loop over hits and fill hit tables
		for(int ihit=0; ihit< nHit ; ihit++)
        {
			TrackerHit * measHit = trackhits.at(ihit);

			// Hit position
			const double * pos = measHit->getPosition();

			
			// find plane number of the hit
			int hitPlane = Utility::getSensorIDfromHit(measHit);

			//set DUT hits in array position behind telescope planes
			if(hitPlane==_DUTid)
			{
				hitPlane = _nTelPlanes-1;
			}
			// Ignore hits not matched to any plane
			if(hitPlane<0 || hitPlane>=_nTelPlanes)
            {
				continue;
            }
			if( (hitCellDecoder(measHit)["properties"] & kFittedHit) == 0 )
		    {
				// Measured hits
				_isMeasured[hitPlane]=true;

				_measuredX[hitPlane]=pos[0];
				_measuredY[hitPlane]=pos[1];
				_measuredZ[hitPlane]=pos[2];

				// Get cluster charge
				_measuredQ[hitPlane]=0.;

				EVENT::LCObjectVec rawdata =  measHit->getRawHits();

				if(rawdata.size()>0 && rawdata.at(0)!=NULL )
		        {
					EUTelVirtualCluster * cluster = new EUTelFFClusterImpl ( static_cast<TrackerDataImpl*> (rawdata.at(0))) ;
					_measuredQ[hitPlane]=cluster->getTotalCharge();
		        }
				message<DEBUG5> ( log() << "Measured hit in plane " << hitPlane << " at  X = "
		                                << pos[0] << ", Y = " << pos[1] << ", Q = " << _measuredQ[hitPlane] );
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
        }// End of loop over hits in track

		//Look for closest DUT hit
		if(_DUTok)
		{
			double distmin=_distMax*_distMax;
			TrackerHit * bestFitHit = 0;
			double dist = distmin;

			for(int ihit=0; ihit< nHitCol ; ihit++)
		    {
				TrackerHit * measHit = dynamic_cast<TrackerHit*>( hitcol->getElementAt(ihit) );

				if (Utility::getSensorIDfromHit(measHit)==_DUTid)//select hits on DUT
				{
					const double * pos = measHit->getPosition();
					//calculate distance between fitted track hit and this hit
					if((hitCellDecoder(measHit)["properties"] & kHitInGlobalCoord) == 0)
					{
						//hit in pseudo-local coordinates
						message<ERROR1> ( log() << "Hits in non-global coordinates not yet implemented");
					}
					else
					{
						//hit in global coordinates, simply subtract
						double resX = pos[0]-_fittedX[_nTelPlanes-1];
						double resY = pos[0]-_fittedX[_nTelPlanes-1];
						double resZ = pos[0]-_fittedX[_nTelPlanes-1];
						dist = resX*resX + resY*resY + resZ*resZ;
					}
					//change minimum if this hit is closer than last minimum
					if(dist < distmin)
					{
						distmin=dist;
						dutX=pos[0];
						dutY=pos[1];
						dutZ=pos[2];
						bestFitHit = measHit;
					}
				}
		    }

			// Try to get DUT cluster charge and size
			if (bestFitHit!=0)
			{
				EVENT::LCObjectVec rawdata =  bestFitHit->getRawHits();

				if(rawdata.size()>0 && rawdata.at(0)!=NULL )
				{
					EUTelVirtualCluster * cluster = new EUTelFFClusterImpl( static_cast<TrackerDataImpl*> (rawdata.at(0)));
					dutQ=cluster->getTotalCharge();
					cluster->getClusterSize(dutClusterSizeX,dutClusterSizeY);
				}
				dutR=sqrt(distmin);
				message<DEBUG5> ( log() << "Matched DUT hit at X = " << dutX << "   Y = " << dutY
			                            << "   Z = " << dutZ << "   Dist = " << dutR << "   Q = " << dutQ );
			}
		}// End of if(_DUTok)

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
		_FitTuple->fill(icol++,dutZ);
		_FitTuple->fill(icol++,dutAlpha);
		_FitTuple->fill(icol++,dutBeta);
		_FitTuple->fill(icol++,dutGamma);
		_FitTuple->fill(icol++,dutR);
		_FitTuple->fill(icol++,dutQ);
		_FitTuple->fill(icol++,dutClusterSizeX);
		_FitTuple->fill(icol++,dutClusterSizeY);

		_FitTuple->addRow();
    }// End of loop over tracks
	return;
}



void EUTelFitTuple::check( LCEvent * /* evt */ ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void EUTelFitTuple::end(){

  //   std::cout << "EUTelFitTuple::end()  " << name()
  //        << " processed " << _nEvt << " events in " << _nRun << " runs "
  //        << std::endl ;


  message<MESSAGE5> ( log() << "N-tuple with "
                     << _FitTuple->rows() << " rows created" );


  // Clean memory

  delete [] _planeSort ;
  delete [] _planePosition ;
  delete [] _planeID ;
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


}



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

  for(int ipl=0; ipl<_nTelPlanes-1;ipl++)
  {
    for(int ivar=0; ivar<7;ivar++)
      {
        stringstream ss;
        ss << _varName[ivar] << "_" << ipl;
        _columnNames.push_back(ss.str());
        _columnType.push_back("double");
      }
  }
  //give DUT (last in array) its proper sensor ID
  for(int ivar=0; ivar<7;ivar++)
  {
    stringstream ss;
    ss << _varName[ivar] << "_" << _DUTid;
    _columnNames.push_back(ss.str());
    _columnType.push_back("double");
  }

  // DUT variables

  _columnNames.push_back("dutX");
  _columnType.push_back("double");

  _columnNames.push_back("dutY");
  _columnType.push_back("double");

  _columnNames.push_back("dutZ");
  _columnType.push_back("double");

  _columnNames.push_back("dutAlpha");
  _columnType.push_back("double");

  _columnNames.push_back("dutBeta");
  _columnType.push_back("double");

  _columnNames.push_back("dutGamma");
  _columnType.push_back("double");

  _columnNames.push_back("dutR");
  _columnType.push_back("double");

  _columnNames.push_back("dutQ");
  _columnType.push_back("double");

  _columnNames.push_back("dutClusterSizeX");
  _columnType.push_back("int");

  _columnNames.push_back("dutClusterSizeY");
  _columnType.push_back("int");

  _FitTuple=AIDAProcessor::tupleFactory(this)->create(_FitTupleName, _FitTupleName, _columnNames, _columnType, "");


  message<DEBUG5> ( log() << "Booking completed \n\n");

  return;
}

#endif // GEAR && AIDA

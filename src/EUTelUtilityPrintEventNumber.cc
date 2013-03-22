#include "EUTelUtilityPrintEventNumber.h"

// C++
#include <iostream>
#include <iomanip>

// Aida
#ifdef MARLIN_USE_AIDA
//AIDA
#include <AIDA/AIDA.h>
#include <marlin/AIDAProcessor.h>
#endif

// LCIO
#include "UTIL/LCTOOLS.h"

using namespace lcio;
using namespace marlin;
using namespace eutelescope;

EUTelUtilityPrintEventNumber aPrintEventNumber ;
    

EUTelUtilityPrintEventNumber::EUTelUtilityPrintEventNumber() : 
  Processor("EUTelUtilityPrintEventNumber") {
  /* the constructor call first the super constructor which the name 
   * of the processor
   * this must by the same as the class name.
   */
	
  /* description of the processor which will be displayed in the
   * steering file autometicly generated by marlin
   */
  _description = "EUTelUtilityPrintEventNumber prints event number to screen"
    " depending on the verbosity level" ;	
	
  /* register steering parameters: name, description, class-variable, default value
   * the type will by definded by the type of the default value
   * string, double, float and int are possible
   */
  registerProcessorParameter( "EveryNEvents", 
			      "Print event number for every n-th event",
			      _everyNEvents, static_cast< int >(1000));
  registerOptionalParameter( "printTimestamp",
			      "print the event timestamp as read from LCIO",
			      _printTimestamp,static_cast< bool >(0));	
}
    
    
void EUTelUtilityPrintEventNumber::init() { 
  // this method is called only once even when the rewind is active

  // It isusually a good idea to
  printParameters ();	
}
    
void EUTelUtilityPrintEventNumber::processRunHeader( LCRunHeader* run) { 
	
      
  run->parameters().setValue( _processorName + "_revision", "$Rev: 699 $");

  for ( ProcParamMap::iterator i = _map.begin(); i != _map.end(); i++ ) {
    if ( ! i->second->isOptional() || i->second->valueSet() ) {
      run->parameters().setValue( _processorName + "_" + i->second->name(), 
				  i->second->value() );
    }
  }

  // Increment the total runs counter:
  totalruns++;
} 
    
void EUTelUtilityPrintEventNumber::processEvent( LCEvent * evt ) { 

  // The regular output in case verbosity MESSAGE is set:
  if ( evt->getEventNumber() <= 10 ||
       (evt->getEventNumber() <= 100 && evt->getEventNumber()%10 == 0) ||
       evt->getEventNumber()%_everyNEvents == 0) {
    streamlog_out(MESSAGE5) << "Processing event " 
			    << std::setw(7) << evt->getEventNumber() 
			    << " in run " << std::setw(6) << std::setfill('0') << evt->getRunNumber();
    if(_printTimestamp) streamlog_out(MESSAGE5) << ", timestamp " << evt->getTimeStamp();
    streamlog_out(MESSAGE5) << std::endl;
  }
  // Additional output if we have the DEBUG level: print every event.
  else {
    streamlog_out(DEBUG5) << "Processing event " 
			    << std::setw(7) << evt->getEventNumber() 
			    << " in run " << evt->getRunNumber();
    if(_printTimestamp) streamlog_out(DEBUG5) << ", timestamp " << evt->getTimeStamp();
    streamlog_out(DEBUG5) << std::endl;
  }

  // Increment the total event counter.
  totalevents++;
}
    
    
void EUTelUtilityPrintEventNumber::check( LCEvent * /* evt */ ) { 
  /* Nothing to do here... */
}

    
    
void EUTelUtilityPrintEventNumber::end(){ 
  streamlog_out(MESSAGE4) << "Finished. Processed " << totalevents 
    			  << " events in " << totalruns << " runs in total."
			  << std::endl ;
}

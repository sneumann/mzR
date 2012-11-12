/***************************************************************************
cramp.cpp

/***************************************************************************
cramp.hpp -- renamed cramp.h to avoid R checker warning 

  A C++ wrapper for the RAMP code.

  Use this library to parse an mzXML file in a non-sequential way, by
  taking advantage of the index element.  

  (C) 2004 by Brian Pratt, Insilicos LLC 

  Based on mzXML2Other, which has this copyright:
	 -------------------
	 begin					 : Wed Apr 2
	 copyright				 : (C) 2002 by Pedrioli Patrick, ISB, Proteomics
	 email					 : ppatrick@student.ethz.ch
 ***************************************************************************/

/***************************************************************************
*																								  *
*	 This program is free software; you can redistribute it and/or modify  *
*	 it under the terms of the GNU Library or "Lesser" General Public 	  *
*	 License (LGPL) as published by the Free Software Foundation;			  *
*	 either version 2 of the License, or (at your option) any later		  *
*	 version.																				  *
*																								  *
***************************************************************************/



#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "stdio.h"
#if !defined(_MSC_VER) && !defined(__MINGW32__)
#include "sys/errno.h"
#endif
#include "cramp.h"
#include<R.h>

/**
 * This function performs a non-sequential parsing operation on an indexed
 * msxml file.
 *
 * @param fileName: Name of the msxml file
 * @param startSCan: Number of the scan we want to read from
 * @param what: -HEADER will return num, msLevel and retentionTime
 *              -SCAN will return only the peaks
 *              -ALL will return everything found in scan, precursorMz and peaks
 *
 * @return pData is dynamically allocate and becomes property of the caller, who
 *         is responsible for its deallocation!!
 */



cRamp::cRamp( const char* fileName,bool declaredScansOnly ) : 
  m_filename(fileName), m_declaredScansOnly(declaredScansOnly), m_runInfo()
{
   m_handle = rampOpenFile(fileName);
   m_scanOffsets = NULL;
   m_runInfo = NULL;
   m_lastScan = 0;   
   if (!OK()) {
     // HENRY -- I would prefer this to be silent, and let the caller deals with it
     // cout << "Error: Could not open file " << fileName << ": " << strerror(errno) << endl;
     // END HENRY
   } else {

      m_runInfo = getRunInfo();

      // HENRY -- always read index to set scan count, since scan count
      // declared at the top of the mzXML file is unreliable now that
      // there are missing scans.
      // This will also set the structs m_scanOffsets, and the value m_lastScan
      
      //      if (m_runInfo->m_data.scanCount < 0) { // undeclared scan count
         // this will provoke reading of index, which sets scan count
         rampScanInfo* tmp = getScanHeaderInfo ( 1 );
         free(tmp);
	 // }
      // END HENRY
   }
}

cRamp::~cRamp() {
   rampCloseFile(m_handle);
   
   // NB: these pointers may be null on file open failure, 
   // but free/delete of NULL is OK per C++ standard 
   free(m_scanOffsets);
   delete m_runInfo; // was free() - but allocated with new 
}

//
// here are the private guts
//
rampInfo* cRamp::do_ramp( ramp_fileoffset_t arg , eWhatToRead	what )
{
   
   switch( what ) {
   case RAMP_RUNINFO:
   case RAMP_HEADER:
   case RAMP_PEAKS:
   case RAMP_INSTRUMENT:
      break; // OK
   default:
     Rf_error("unknown read type!\n");
      return NULL;
      break;
   }	
   
   rampInfo* returnPtr=NULL;
   
   if ((RAMP_RUNINFO != what) && (RAMP_INSTRUMENT != what) && !m_scanOffsets) {
      int iLastScan = 0; 
     // we need the index to get anything besides the header
      //      std::cerr << "in: getIndexOffset(m_handle);!\n";
      ramp_fileoffset_t indexOffset = getIndexOffset(m_handle);
      //      std::cerr << "out: getIndexOffset(m_handle);!\n";

      // std::cerr << "in: readIndex();!\n";
      m_scanOffsets = readIndex(m_handle, indexOffset, &iLastScan);
      //std::cerr << "out: readIndex();!\n";
      if (iLastScan >= m_runInfo->m_data.scanCount) {
		 if (!m_declaredScansOnly) {
           m_runInfo->m_data.scanCount = iLastScan;
		 } else { // get rid of all the fake entries created

		   //  std::cerr << "get rid of all the fake entries created\n";
			 for (int n=1;n<=iLastScan;n++) { // ramp is 1 based
				 if (m_scanOffsets[n]==-1) {
					// find a run of fakes
				    int m;
					for (m=n+1;(m<=iLastScan)&&(m_scanOffsets[m]==-1);m++);
					if (m<=iLastScan) {
						memmove(m_scanOffsets+n,m_scanOffsets+m,
						  sizeof(ramp_fileoffset_t)*((iLastScan-m)+1));
					}
					iLastScan-=(m-n);
				 }
			 }
		 }
      }
      // HENRY - store last scan explicitly.
      m_lastScan = iLastScan;
      // END HENRY
   }

   
   // HENRY -- arg is out of bounds. instead of creating havoc in RAMP, let's just kill it here.
   if (RAMP_RUNINFO != what && (RAMP_INSTRUMENT != what) && (arg > m_runInfo->m_data.scanCount || arg < 1)) {
     return (NULL);
   }
     
   if (m_scanOffsets || (RAMP_RUNINFO == what) || (RAMP_INSTRUMENT == what)) {
      ramp_fileoffset_t scanOffset=-1;
      if (RAMP_RUNINFO == what || RAMP_INSTRUMENT == what) {
	scanOffset = 0; // read from head of file
      } else {
	scanOffset = m_scanOffsets[arg]; // ramp is one-based
      }
      
      if (scanOffset >= 0) {
         
         // -----------------------------------------------------------------------
         // And now we can parse the info we were looking for
         // -----------------------------------------------------------------------
         
         
         // Ok now we have to copy everything in our structure
         switch( what )
         {
         case RAMP_RUNINFO:
	   //	   std::cerr << "in: rampRunInfo( m_handle )\n";
            returnPtr = new rampRunInfo( m_handle );
	    //  std::cerr << "out: rampRunInfo( m_handle )\n";
            break;
         case RAMP_HEADER:
	   //  std::cerr <<  "in: rampScanInfo( m_handle, scanOffset, (int)arg )\n" ;
            returnPtr = new rampScanInfo( m_handle, scanOffset, (int)arg );
	    //std::cerr <<  "out: rampScanInfo( m_handle, scanOffset, (int)arg )\n" ;
            if (returnPtr) {
#ifdef HAVE_PWIZ_MZML_LIB
			   if (!m_handle->mzML) // rampadapter already set this for us
#endif
              ((rampScanInfo *)returnPtr)->m_data.filePosition = scanOffset; // for future reference
            
              // HENRY -- error checking here
              if (((rampScanInfo*)returnPtr)->m_data.acquisitionNum < 0) {
                // something failed in RAMP, possibly because it's a missing scan
                delete ((rampScanInfo*)returnPtr);
                returnPtr = NULL;
              }
            }
            break;           
         case RAMP_PEAKS:
            returnPtr = new rampPeakList( m_handle, scanOffset);
            
            // HENRY -- error checking here
            if (returnPtr && ((rampPeakList*)returnPtr)->getPeakCount() <= 0) {
              // something failed in RAMP, possibly because it's a missing scan
              delete ((rampPeakList*)returnPtr);
              returnPtr = NULL;
            }
            break;
            
         // HENRY -- add the instrument info reading functionality (present in RAMP, but not provided in cRAMP before)
         case RAMP_INSTRUMENT:
	   //	   std::cerr <<  "in: rampInstrumentInfo(m_handle)\n" ;
            returnPtr = new rampInstrumentInfo(m_handle);
	    //std::cerr <<  "out: rampInstrumentInfo(m_handle)\n" ;
            if (((rampInstrumentInfo*)returnPtr)->m_instrumentStructPtr == NULL) {
              delete ((rampInstrumentInfo*)returnPtr);
              returnPtr = NULL;
            }
            break;
         }
         
      }
   }
   
   
   
   return returnPtr;
}


/**
 * This function performs a non-sequential parsing operation on an indexed
 * msxml file to obtain minimal info on the msRun contained in the file.  

 *
 * @return rapRunInfo* is dynamically allocate and becomes property of the caller, who
 *         is responsible for its deallocation!!
 */

rampRunInfo* cRamp::getRunInfo (  ) {
   rampRunInfo* result;
   if (m_runInfo) { // did we derive this already?
      result = new rampRunInfo(*m_runInfo);
   } else {
      result = (rampRunInfo*) do_ramp(0, RAMP_RUNINFO);
   }
   return result;
}

/**
 * This function performs a non-sequential parsing operation on an indexed
 * msxml file to obtain minimal header info for a numbered scan (thus minimizing parse time).  
 *
 * @param fileName: Name of the msxml file
 * @param startSCan: Number of the scan we want to read from
 * @return rapHeaderInfo* is dynamically allocate and becomes property of the caller, who
 *         is responsible for its deallocation!! returns just the minimal header info num, msLevel and retentionTime
 */

rampScanInfo* cRamp::getScanHeaderInfo ( int whichScan  ) {
   return (rampScanInfo*) do_ramp((ramp_fileoffset_t)whichScan, RAMP_HEADER);
}


/**
 * This function performs a non-sequential parsing operation on an indexed
 * msxml file to obtain peak info for a numbered scan.
 *
 * @param fileName: Name of the msxml file
 * @param startSCan: Number of the scan we want to read from
 * @return rapPeakList* is dynamically allocate and becomes property of the caller, who
 *         is responsible for its deallocation!! returns everything found in scan, precursorMz and peaks
 */

rampPeakList* cRamp::getPeakList ( int whichScan ) {
   return (rampPeakList*) do_ramp((ramp_fileoffset_t)whichScan, RAMP_PEAKS);
}

// HENRY - provides instrument info getting method
rampInstrumentInfo* cRamp::getInstrumentInfo () {
  return (rampInstrumentInfo*) do_ramp(0, RAMP_INSTRUMENT);
}
// END HENRY


// HENRY - sequential access parser that skips over missing scans. This version only reads scan header.
bool cRampIterator::nextScan(rampScanInfo** scanInfo) {
	while (++m_currentScan <= m_cramp.getLastScan() && m_cramp.getScanOffset(m_currentScan) <= 0);
  if (m_currentScan > m_cramp.getLastScan()) {
    return (false);
  }
  
  *scanInfo = (rampScanInfo*)m_cramp.do_ramp((ramp_fileoffset_t)(m_currentScan), RAMP_HEADER);
  return (true);
}
// END HENRY

// HENRY - sequential access parser that skips over missing scans. This version reads both scan header and peak list.
bool cRampIterator::nextScan(rampScanInfo** scanInfo, rampPeakList** peakList) {
  while (++m_currentScan <= m_cramp.getLastScan() && m_cramp.getScanOffset(m_currentScan) <= 0);
  if (m_currentScan > m_cramp.getLastScan()) {
    return (false);
  }
  
  *scanInfo = (rampScanInfo*)m_cramp.do_ramp((ramp_fileoffset_t)(m_currentScan), RAMP_HEADER);
  *peakList = (rampPeakList*)m_cramp.do_ramp((ramp_fileoffset_t)(m_currentScan), RAMP_PEAKS);

  return (true);

}
// END HENRY

// HENRY - resets the sequential access parser to the first scan.
void cRampIterator::reset() {
  m_currentScan = 1;
}
// END HENRY

/**
 * populate from a file handle
 **/
rampPeakList::rampPeakList(RAMPFILE *handle, ramp_fileoffset_t index) {
   init();
   m_peaksCount = readPeaksCount(handle,index);
   m_pPeaks = (rampPeakInfoStruct *)readPeaks(handle,index);
}

/**
 * populate from a file handle
 **/
rampScanInfo::rampScanInfo(RAMPFILE *handle, ramp_fileoffset_t index, int seqNum) {
   init();
   readHeader(handle,index,&m_data);
   m_data.seqNum = seqNum;
}

/**
 * populate from a file handle
 **/
rampRunInfo::rampRunInfo(RAMPFILE *handle) {
   init();
   readMSRun(handle,&m_data);
}

// HENRY - provides instrument info reading functionality
rampInstrumentInfo::rampInstrumentInfo(RAMPFILE *handle) {
  init();
  m_instrumentStructPtr = getInstrumentStruct(handle);
}
// END HENRY

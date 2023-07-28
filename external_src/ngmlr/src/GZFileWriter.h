/**
 * Contact: philipp.rescheneder@gmail.com
 */

#ifndef GZFILEWRITER_H_
#define GZFILEWRITER_H_

#include <stdio.h>
#include <stdarg.h>
#include <iostream>
#include <cstring>

#include "NGMThreads.h"
#include "ILog.h"
#include "zlib.h"

class GZFileWriter : public FileWriter {

public:

	gzFile m_Output;

	GZFileWriter(char const * const filename) {
		std::string strFilename = std::string(filename);
		if(!hasSuffix(strFilename, ".gz")) {
			strFilename += ".gz";
		}
		if (!(m_Output = gzopen(strFilename.c_str(), "wb"))) {
			Log.Error("Unable to open output file %s", filename);
		}
	}

	~GZFileWriter() {
		gzclose(m_Output);
	}

	void doFlush(int & bufferPosition, int const BUFFER_LIMIT, char * writeBuffer, bool last = false) {

		if (bufferPosition > BUFFER_LIMIT || last) {
			if(gzwrite(m_Output, writeBuffer, bufferPosition) < 0) {
				Log.Error("Writing");
			}
			bufferPosition = 0;

		}

	}

private:
	bool hasSuffix(const std::string &str, const std::string &suffix)
	{
	    return str.size() >= suffix.size() &&
	           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
	}

};

#endif /* GZFILEWRITER_H_ */

// random_access_gzFile.cpp
// presents C-style calls based on an ifstream-derived C++ class
//
// like gzio.c from the zlib distro,
// but for dealing with reading >2GB files in win32, and efficient random access.
// Copyright (C) Insilicos LLC 2008, ALl Rights Reserved.
// For conditions of distribution and use, see copyright notice in zlib.h
//
// based on:
/* gzio.c -- IO on .gz files
* Copyright (C) 1995-2005 Jean-loup Gailly.
* For conditions of distribution and use, see copyright notice in zlib.h
*/

// efficient random access stuff based on
/* zran.c -- example of zlib/gzip stream indexing and random access
* Copyright (C) 2005 Mark Adler
* For conditions of distribution and use, see copyright notice in zlib.h
Version 1.0  29 May 2005  Mark Adler */

#include "random_access_gzFile.h"
#include "zlib.h"

struct random_access_gzFile {
	pwiz::util::random_access_compressed_ifstream *istrm;
};

/* ===========================================================================
Opens a gzip (.gz) file for reading or writing. The mode parameter
is as in fopen ("rb" or "wb"). The file is given either by file descriptor
or path name (if fd == -1).
random_access_gzopen returns NULL if the file could not be opened or if there was
insufficient memory to allocate the (de)compression state; errno
can be checked to distinguish the two cases (if errno is zero, the
zlib error is Z_MEM_ERROR).
*/
random_access_gzFile *random_access_gzopen (const char *path)
{
	random_access_gzFile *ret = new random_access_gzFile;
	if (ret) {
		ret->istrm = new pwiz::util::random_access_compressed_ifstream(path);
		if (ret->istrm && !*(ret->istrm)) {
			delete ret->istrm;
			delete ret;
			ret = NULL;
		} 
	}
	return ret;
}

/* ===========================================================================
Reads the given number of uncompressed bytes from the compressed file.
gzread returns the number of bytes actually read (0 for end of file).
*/
int random_access_gzread (random_access_gzFile *file,
				 char *buf,
				 unsigned len)
{
	file->istrm->clear(); // clear any stale failbit
	file->istrm->read((char *)buf,len);
	return file->istrm->gcount();
}
/* ===========================================================================
Reads one byte from the compressed file. gzgetc returns this byte
or -1 in case of end of file or error.
*/
int random_access_gzgetc(random_access_gzFile *file)
{
	char c;
	file->istrm->clear(); // clear any stale failbit
	file->istrm->read( &c, 1);
	return (file->istrm->gcount() == 1) ? (int)c : -1;
}

/* ===========================================================================
Reads bytes from the compressed file until len-1 characters are
read, or a newline character is read and transferred to buf, or an
end-of-file condition is encountered.  The string is then terminated
with a null character.
gzgets returns buf, or Z_NULL in case of error.
*/
char * random_access_gzgets(random_access_gzFile *file, char *buf, int len) {
	file->istrm->clear(); // clear any stale failbit
	file->istrm->getline(buf,len); // this does NOT retain the \n
	int got = file->istrm->gcount();
	if (got && (got < (len-1)) && !file->istrm->fail() && !file->istrm->eof()) {
		// looks like a newline was read then discarded
		buf[(got-1)] = '\n';
		buf[got] = 0;
	}
	return (*buf?buf:NULL);
}

/* ===========================================================================
Sets the starting position for the next gzread on the given
compressed file. The offset represents a number of bytes in the
gzseek returns the resulting offset location as measured in bytes from
the beginning of the uncompressed stream, or -1 in case of error.
NOTE - bpratt borrows zran code here for better performance
*/
random_access_gzFile_off_t random_access_gzseek (random_access_gzFile *file,
							   random_access_gzFile_off_t offset,
							   int whence)
{
	if (file == NULL) {
			return -1L;
	}
	file->istrm->clear(); // clear any stale error flags
	switch (whence) {
		case SEEK_SET:
			file->istrm->seekg(boost::iostreams::offset_to_position(offset));
			break;
		case SEEK_CUR:
			file->istrm->seekg(boost::iostreams::offset_to_position(offset),std::ios_base::cur);
			break;
		case SEEK_END:
			file->istrm->seekg(boost::iostreams::offset_to_position(offset),std::ios_base::end);
			break;
	}
	return file->istrm->tellg();
}

/* ===========================================================================
Returns the starting position for the next gzread on the
given compressed file. This position represents a number of bytes in the
uncompressed data stream.
*/
random_access_gzFile_off_t random_access_gztell(random_access_gzFile *file)
{
	return file?(random_access_gzFile_off_t)(file->istrm->tellg()):-1;
}

/* ===========================================================================
Returns 1 when EOF has previously been detected reading the given
input stream, otherwise zero.
*/
int random_access_gzeof (random_access_gzFile *s)
{
	return s?s->istrm->eof():0;
}

/* ===========================================================================
Returns 1 if reading and doing so transparently, otherwise zero.
*/
int random_access_gzdirect(random_access_gzFile *file)
{
	if ((file == NULL) || (file->istrm == NULL)) {
		return 0;
	}
	return !file->istrm->getCompressionType();
}

/* ===========================================================================
Flushes all pending output if necessary, closes the compressed file
and deallocates all the (de)compression state.
*/
int random_access_gzclose (random_access_gzFile *file)
{
	if (file == NULL) {
		return Z_STREAM_ERROR;
	}
	delete file->istrm;
	delete file;
	return 0;
}

CC = gcc
CXX = g++

CXXFLAGS = -I./boost_aux/ -I./boost -I. -DHAVE_PWIZ_MZML_LIB  -DWITHOUT_MZ5 -fpermissive -DWINDOWS_NATIVE -s -O2 -m64

.c.o:
	$(CC) $(CXXFLAGS) -c -o $@ $<

.cpp.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

PWIZOBJECTS=\
./boost/system/src/error_code.o \
./boost/regex/src/posix_api.o \
./boost/regex/src/fileiter.o \
./boost/regex/src/regex_raw_buffer.o \
./boost/regex/src/cregex.o \
./boost/regex/src/regex_debug.o \
./boost/regex/src/instances.o \
./boost/regex/src/icu.o \
./boost/regex/src/usinstances.o \
./boost/regex/src/regex.o \
./boost/regex/src/wide_posix_api.o \
./boost/regex/src/regex_traits_defaults.o \
./boost/regex/src/winstances.o \
./boost/regex/src/wc_regex_traits.o \
./boost/regex/src/c_regex_traits.o \
./boost/regex/src/cpp_regex_traits.o \
./boost/regex/src/static_mutex.o \
./boost/regex/src/w32_regex_traits.o \
./boost/iostreams/src/zlib.o \
./boost/iostreams/src/file_descriptor.o \
./boost/filesystem/src/operations.o \
./boost/filesystem/src/path.o \
./boost/filesystem/src/path_traits.o \
./boost/filesystem/src/utf8_codecvt_facet.o \
./boost/filesystem/src/windows_file_codecvt.o \
./boost/filesystem/src/codecvt_error_category.o \
./boost/thread/src/win32/thread.o \
./boost_aux/boost/nowide/iostream.o \
./pwiz/data/common/MemoryIndex.o \
./pwiz/data/common/CVTranslator.o \
./pwiz/data/common/cv.o \
./pwiz/data/common/ParamTypes.o \
./pwiz/data/common/BinaryIndexStream.o \
./pwiz/data/common/diff_std.o \
./pwiz/data/common/Unimod.o \
./pwiz/data/msdata/Version.o \
./pwiz/data/msdata/SpectrumList_MGF.o \
./pwiz/data/msdata/DefaultReaderList.o \
./pwiz/data/msdata/ChromatogramList_mzML.o \
./pwiz/data/msdata/examples.o \
./pwiz/data/msdata/Serializer_mzML.o \
./pwiz/data/msdata/Serializer_MSn.o \
./pwiz/data/msdata/Reader.o \
./pwiz/data/msdata/Serializer_MGF.o \
./pwiz/data/msdata/Serializer_mzXML.o \
./pwiz/data/msdata/SpectrumList_mzML.o \
./pwiz/data/msdata/SpectrumList_MSn.o \
./pwiz/data/msdata/BinaryDataEncoder.o \
./pwiz/data/msdata/Diff.o \
./pwiz/data/msdata/MSData.o \
./pwiz/data/msdata/References.o \
./pwiz/data/msdata/SpectrumList_mzXML.o \
./pwiz/data/msdata/IO.o \
./pwiz/data/msdata/SpectrumList_BTDX.o \
./pwiz/data/msdata/SpectrumInfo.o \
./pwiz/data/msdata/RAMPAdapter.o \
./pwiz/data/msdata/LegacyAdapter.o \
./pwiz/data/msdata/SpectrumIterator.o \
./pwiz/data/msdata/MSDataFile.o \
./pwiz/data/msdata/MSNumpress.o \
./pwiz/data/msdata/SpectrumListCache.o \
./pwiz/data/msdata/Index_mzML.o \
./pwiz/data/msdata/ramp/wglob.o \
./pwiz/data/identdata/IdentDataFile.o \
./pwiz/data/identdata/IdentData.o \
./pwiz/data/identdata/DefaultReaderList.o \
./pwiz/data/identdata/Reader.o \
./pwiz/data/identdata/Serializer_protXML.o \
./pwiz/data/identdata/Serializer_pepXML.o \
./pwiz/data/identdata/Serializer_mzid.o \
./pwiz/data/identdata/IO.o \
./pwiz/data/identdata/Version.o \
./pwiz/data/identdata/References.o \
./pwiz/data/proteome/Modification.o \
./pwiz/data/proteome/Digestion.o \
./pwiz/data/proteome/Peptide.o \
./pwiz/data/proteome/AminoAcid.o \
./pwiz/utility/chemistry/Chemistry.o \
./pwiz/utility/chemistry/ChemistryData.o \
./pwiz/utility/chemistry/MZTolerance.o \
./pwiz/utility/misc/IntegerSet.o \
./pwiz/utility/misc/Base64.o \
./pwiz/utility/misc/IterationListener.o \
./pwiz/utility/misc/MSIHandler.o \
./pwiz/utility/misc/Filesystem.o \
./pwiz/utility/misc/TabReader.o \
./pwiz/utility/misc/random_access_compressed_ifstream.o \
./pwiz/utility/misc/SHA1.o \
./pwiz/utility/misc/SHA1Calculator.o \
./pwiz/utility/minimxml/XMLWriter.o \
./pwiz/utility/minimxml/SAXParser.o \
./random_access_gzFile.o \
adler32.o compress.o crc32.o deflate.o gzclose.o gzlib.o gzread.o gzwrite.o infback.o inffast.o inflate.o inftrees.o trees.o uncompr.o zutil.o


all: libpwiz.a

libpwiz.a: $(PWIZOBJECTS)
	ar -rc $@ $(PWIZOBJECTS)
	
clean:
	rm -rf $(PWIZOBJECTS)
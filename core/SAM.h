#ifndef __SAMSTRUCTURES_H__
#define __SAMSTRUCTURES_H__

#include "StringArray.h"
#include "StringHash.h"
#include "MathVector.h"
#include "InputFile.h"
#include "LongHash.h"
#include "IntArray.h"

struct samRecordHeader
   {
   public:
      int            blockSize;
      int            referenceID;
      int            position;
      unsigned int   readNameLength : 8, mapQuality : 8, bin : 16;
      int            cigarLength : 16, flag : 16;
      int            readLength;
      int            mateReferenceID;
      int            matePosition;
      int            insertSize;             // Outer fragment length
   };

struct samRecordStruct : public samRecordHeader
   {
   public:
      char  data[1];
   };

class samRecord
   {
   public:
      samRecordStruct * header;

      String         readName;
      String         cigar;
      String         sequence;
      String         quality;

      LongHash<int>  extras;
      StringArray    strings;
      IntArray       integers;
      Vector         doubles;

      samRecord();
      ~samRecord();

      void allocateStructure(int size);
      void clearExtras();

      String & getString(const char * tag);
      int &    getInteger(const char * tag);
      double & getDouble(const char * tag);

      bool     checkString(const char * tag)    { return checkTag(tag, 'Z'); }
      bool     checkInteger(const char * tag)   { return checkTag(tag, 'i'); }
      bool     checkDouble(const char * tag)    { return checkTag(tag, 'f'); }
      bool     checkTag(const char * tag, char type);

      static int MAKEKEY(char ch1, char ch2, char type)
            { return (type << 16) + (ch2 << 8) + ch1; }

   private:
      int allocatedSize;
   };

class samFile
   {
   public:
      String         header;
      StringArray    referenceContigs;
      StringIntHash  referenceHash;
      IntArray       referenceLengths;

      samRecord      data;

      bool  isOpen;
      bool  isBinary;

      bool  Open(const char * filename);
      bool  ReadRecord();
      bool  ReadBAM();
      bool  ReadSAM();

      int   GetReferenceID(const String & referenceName);
      const String & GetReferenceLabel(int id);

   private:
      IFILE       input;
      String      buffer;
      StringArray tokens;

      void ParseHeaderLine(StringIntHash & tags, StringArray & values);
   };

#endif


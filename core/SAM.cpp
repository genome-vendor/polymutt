#include "SAM.h"

#include <stdlib.h>

samRecord::samRecord()
   {
   header = (samRecordStruct *) malloc(sizeof(samRecordStruct) + sizeof(int));

   allocatedSize = header->blockSize = sizeof(samRecordStruct);
   }

samRecord::~samRecord()
   {
   free(header);
   }

void samRecord::allocateStructure(int size)
   {
   if (allocatedSize < size)
      {
      free(header);

      header = (samRecordStruct *) malloc(size + sizeof(int));

      allocatedSize = size;
      }

   header->blockSize = size;
   }

void samRecord::clearExtras()
   {
   extras.Clear();
   strings.Clear();
   integers.Clear();
   doubles.Clear();
   }

bool samFile::Open(const char * filename)
   {
   if (input != NULL)
      ifclose(input);

   referenceContigs.Clear();
   referenceHash.Clear();
   referenceLengths.Clear();

   input = ifopen(filename, "rb");

   if (input == NULL)
      return false;

   char magic[4];
   ifread(input, magic, 4);

   if (magic[0] == 'B' && magic[1] == 'A' && magic[2] == 'M' && magic[3] == 1)
      isBinary = true;
   else
      {
      isBinary = false;
      ifrewind(input);
      }

   isOpen = true;

   if (isBinary) // Read BAM Header
      {
      int headerLength;
      ifread(input, &headerLength, sizeof(int));

      if (headerLength > 0)
         {
         ifread(input, header.LockBuffer(headerLength + 1), headerLength);
         header[headerLength] = 0;
         header.UnlockBuffer();
         }

      int referenceCount;
      ifread(input, &referenceCount, sizeof(int));

      referenceContigs.Dimension(referenceCount);
      referenceLengths.Dimension(referenceCount);
      for (int i = 0; i < referenceCount; i++)
         {
         int nameLength;
         ifread(input, &nameLength, sizeof(int));

         ifread(input, referenceContigs[i].LockBuffer(nameLength), nameLength);
         ifread(input, &referenceLengths[i], sizeof(int));
         referenceContigs[i].UnlockBuffer();

         referenceHash.Add(referenceContigs[i], i);
         }
      }
   else // Read SAM Header
      {
      do {
         StringIntHash tags;
         StringArray   values;

         buffer.ReadLine(input);

         if (ifeof(input) || buffer[0] != '@') break;

         if (buffer[1] != 'S' || buffer[2] != 'Q') continue;

         ParseHeaderLine(tags, values);

         int name = tags.Integer("SN");
         int length = tags.Integer("LN");

         if (name < 0 || length < 0) continue;

         referenceHash.Add(values[name], referenceContigs.Length());

         referenceContigs.Push(values[name]);
         referenceLengths.Push(values[length].AsInteger());

      } while (1);

      ifrewind(input);
      }

   return true;
   }

void samFile::ParseHeaderLine(StringIntHash & tags, StringArray & values)
   {
   tags.Clear();
   values.Clear();

   tokens.AddColumns(buffer, '\t');

   for (int i = 1; i < tokens.Length(); i++)
      {
      tags.Add(tokens[i].Left(2), i - 1);
      values.Push(tokens[i].SubStr(3));
      }
   }

bool samFile::ReadSAM()
   {
   do {
      buffer.Clear();
      buffer.ReadLine(input);
      if (ifeof(input)) return false;
   } while (buffer[0] == '@' || buffer[0] == 0);

   tokens.ReplaceColumns(buffer, '\t');

   if (tokens.Length() < 11)
      return false;

   data.readName = tokens[0];
   data.header->flag = tokens[1].AsInteger();
   data.header->referenceID = GetReferenceID(tokens[2]);
   data.header->position = tokens[3].AsInteger() - 1;
   data.header->mapQuality = tokens[4].AsInteger();
   data.cigar = tokens[5];
   data.header->mateReferenceID = GetReferenceID(tokens[6]);
   data.header->matePosition = tokens[7].AsInteger() - 1;
   data.header->insertSize = tokens[8].AsInteger();
   data.sequence = tokens[9];
   data.quality = tokens[10];

   data.clearExtras();

   for (int i = 11; i < tokens.Length(); i++)
      {
      String & nugget = tokens[i];

      if (nugget.Length() < 6 || nugget[2] != ':' || nugget[4] != ':') continue;

      int key = data.MAKEKEY(nugget[0], nugget[1], nugget[3]);
      int value;

      switch (nugget[3])
         {
         case 'i' :
            value = data.integers.Length();
            data.integers.Push(atoi((const char *) nugget + 5));
            break;
         case 'Z' :
            value = data.strings.Length();
            data.strings.Push((const char *) nugget + 5);
            break;
         case 'f' :
            value = data.doubles.Length();
            data.doubles.Push(atof((const char *) nugget + 5));
            break;
         default :
            printf("samFile::ReadSAM() - Unknown custom field of type %.4s\n",
                   (const char *) tokens[i]);
         }

      data.extras.Add(key, value);
      }

   return true;
   }

bool samFile::ReadBAM()
   {
   int recordSize;

   ifread(input, &recordSize, sizeof(int));

   if (ifeof(input))
      return false;

   data.allocateStructure(recordSize);

   unsigned char * ptr = (unsigned char *) (void *) data.header;

   ifread(input, ptr + 4, recordSize);

   data.readName = (char *) ptr + sizeof(samRecordHeader);

   int length = data.header->readLength;

   data.sequence.SetLength(length);
   data.quality.SetLength(length);

   unsigned char * readNameEnds = ptr + sizeof(samRecordHeader) + data.header->readNameLength;
   unsigned char * packedSequence = readNameEnds + data.header->cigarLength * sizeof(int);
   unsigned char * packedQuality = packedSequence + (length + 1) / 2;

   const char * asciiBases = "0AC.G...T......N";

   for (int i = 0; i < data.header->readLength; i++) {
      data.sequence[i] = i & 1 ?
            asciiBases[packedSequence[i / 2] & 0xF] :
            asciiBases[packedSequence[i / 2] >> 4];

      data.quality[i] = packedQuality[i] + 33;
   }

   data.cigar.Clear();

   unsigned int * packedCigar = (unsigned int *) (void *) readNameEnds;
   const char * asciiCigar = "MIDNSHP";

   for (int i = 0; i < data.header->cigarLength; i++)
      data.cigar.catprintf("%d%c",
                           packedCigar[i] >> 4,
                           asciiCigar[packedCigar[i] & 0xF]);

   unsigned char * extraPtr = packedQuality + length;

   data.clearExtras();
   while (recordSize + 4 - (extraPtr - ptr) > 0)
      {
      int key;
      int value;

      void * content = extraPtr + 3;

      switch (extraPtr[2])
         {
         case 'A' :
         case 'c' :
            key = data.MAKEKEY(extraPtr[0], extraPtr[1], 'i');
            value = data.integers.Length();
            data.integers.Push(* (char *) content);
            extraPtr += 4;
            break;
         case 'C' :
            key = data.MAKEKEY(extraPtr[0], extraPtr[1], 'i');
            value = data.integers.Length();
            data.integers.Push(* (unsigned char *) content);
            extraPtr += 4;
            break;
         case 's' :
            key = data.MAKEKEY(extraPtr[0], extraPtr[1], 'i');
            value = data.integers.Length();
            data.integers.Push(* (short *) content);
            extraPtr += 5;
            break;
         case 'S' :
            key = data.MAKEKEY(extraPtr[0], extraPtr[1], 'i');
            value = data.integers.Length();
            data.integers.Push(* (unsigned short *) content);
            extraPtr += 5;
            break;
         case 'i' :
            key = data.MAKEKEY(extraPtr[0], extraPtr[1], 'i');
            value = data.integers.Length();
            data.integers.Push(* (int *) content);
            extraPtr += 7;
            break;
         case 'I' :
            key = data.MAKEKEY(extraPtr[0], extraPtr[1], 'i');
            value = data.integers.Length();
            data.integers.Push((int) * (unsigned int *) content);
            extraPtr += 7;
            break;
         case 'Z' :
            key = data.MAKEKEY(extraPtr[0], extraPtr[1], 'Z');
            value = data.strings.Length();
            data.strings.Push((const char *) content);
            extraPtr += 4 + data.strings.Last().Length();
            break;
         case 'f' :
            key = data.MAKEKEY(extraPtr[0], extraPtr[1], 'f');
            value = data.doubles.Length();
            data.doubles.Push(* (float *) content);
            extraPtr += 7;
            break;
         default :
            printf("samFile::ReadBAM() - Unknown custom field of type %c%c:%c\n",
                   extraPtr[0], extraPtr[1], extraPtr[2]);
         }

      data.extras.Add(key, value);
      }

   return true;
   }

bool samFile::ReadRecord()
   {
   if (isBinary)
      return ReadBAM();
   else
      return ReadSAM();
   }

int samFile::GetReferenceID(const String & referenceName)
   {
   if (referenceName == "*")
      return -1;

   int id = referenceHash.Find(referenceName);

   if (id >= 0)
      return referenceHash.Integer(id);

   id = referenceContigs.Length();
   referenceContigs.Push(referenceName);
   referenceLengths.Push(0);
   referenceHash.Add(referenceName, id);

   return id;
   }

const String & samFile::GetReferenceLabel(int id)
   {
   static String noname("*");

   if (id == -1) return noname;

   return referenceContigs[id];
   }

bool samRecord::checkTag(const char * tag, char type)
   {
   int key = MAKEKEY(tag[0], tag[1], type);

   return (extras.Find(key) != LH_NOTFOUND);
   }

String & samRecord::getString(const char * tag)
   {
   int key = MAKEKEY(tag[0], tag[1], 'Z');
   int offset = extras.Find(key);

   int value;
   if (offset < 0)
      {
      value = strings.Length();
      strings.Push("");
      extras.Add(key, value);
      }
   else
      value = extras[offset];

   return strings[value];
   }

int & samRecord::getInteger(const char * tag)
   {
   int key = MAKEKEY(tag[0], tag[1], 'i');
   int offset = extras.Find(key);

   int value;
   if (offset < 0)
      {
      value = integers.Length();
      integers.Push(0);
      extras.Add(key, value);
      }
   else
      value = extras[offset];

   return integers[value];
   }

double & samRecord::getDouble(const char * tag)
   {
   int key = MAKEKEY(tag[0], tag[1], 'f');
   int offset = extras.Find(key);

   int value;
   if (offset < 0)
      {
      value = doubles.Length();
      doubles.Push(0);
      extras.Add(key, value);
      }
   else
      value = extras[offset];

   return doubles[value];
   }



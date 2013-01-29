#ifndef __FASTALIKE__
#define __FASTALIKE__

#include "SequenceBasics.h"
#include "StringBasics.h"
#include "InputFile.h"

class FastaLike : public SequenceBasics
   {
   public:
      IFILE    input;
      String   label;
      String   shortLabel;
      String   sequence;
      String   quality;

      FastaLike()  { }
      ~FastaLike() { CloseArchive(); }

      bool OpenArchive(const char * filename);
      void CloseArchive();

      int Load(IFILE & source, String & shortLabel, String & label, String & sequence);
      int Load(IFILE & source)  { return Load(source, shortLabel, label, sequence); }
      int Load()                { return Load(input); }

      int Length()     { return sequence.Length(); }

   private:
      String buffer;
   };

#endif


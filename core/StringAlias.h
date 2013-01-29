#ifndef __STRINGALIAS_H__
#define __STRINGALIAS_H__

#include "StringArray.h"
#include "StringHash.h"

class StringAlias
   {
   public:
      StringAlias()               {}
      virtual ~StringAlias()      {}

      void SetAlias(String & string, String & alias);

      const String & GetAlias(const String & string) const;

      bool  ReadFromFile(const char * filename);
      bool  ReadFromFile(IFILE & input);

   private:
      StringIntHash  lookup;
      StringArray    aliases;
   };

#endif


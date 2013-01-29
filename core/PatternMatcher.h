#ifndef __PATTERNMATCHER_H__
#define __PATTERNMATCHER_H__

#include "StringBasics.h"
#include "IntArray.h"

class PatternMatcher
   {
   public:
      bool CompilePattern(const char * input);
      int  ComparePattern(const char * string);

   private:
      // Convenience function to compare one character in a pattern
      static bool CompareCharacter(const int * mask, const char * string);

      // Convenience function that allows for internal pattern matches
      static int  ComparePattern(const int * pattern, const char * input);

      // Convenience function to advance pattern pointer
      static void AdvancePointer(const int * & pointer);

      // Convenience functions for processing character sets
      static bool CompileMask(int * mask, const char * & definition);

      // Convenience function for handling escape codes
      static int  CompileCharacter(const char * & definition);

      // Convenience functions for dealing with hex character codes
      static bool isHexDigit(char ch);
      static int  convertHexDigit(char ch);

      // Convenience functions for dealing with bitmaps
      static void SetMask(int * mask, int bit);
      static bool GetMask(const int * mask, int bit);

      IntArray compiledPattern;
   };

#endif

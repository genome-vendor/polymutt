#include "PatternMatcher.h"

#include <ctype.h>
#include <limits.h>

// TODO : Write function to evaluate matches to compiled pattern
// TODO : Write function to evaluate single characters in compiled pattern

#define  BITS_PER_INT         (sizeof(int) * CHAR_BIT)
#define  MAP_SIZE             (256 / BITS_PER_INT)

#define  PAT_BEGIN_MASK      '['
#define  PAT_END_MASK        ']'
#define  PAT_REVERSE_MASK    '^'
#define  PAT_OPTIONAL        '?'
#define  PAT_ONE_OR_MORE     '+'
#define  PAT_ZERO_OR_MORE    '*'
#define  PAT_ANY             '.'
#define  PAT_START           '^'
#define  PAT_END             '$'

#define  PM_MASK             0x1000
#define  PM_OPTIONAL         0x2000
#define  PM_ONE_OR_MORE      0x3000
#define  PM_ZERO_OR_MORE     0x4000
#define  PM_ANY              0x5000
#define  PM_START            0x6000
#define  PM_END              0x7000

bool PatternMatcher::CompilePattern(const char * pattern)
   {
   compiledPattern.Clear();

   if (*pattern == 0 || *pattern == PAT_OPTIONAL ||
       *pattern == PAT_ONE_OR_MORE || *pattern == PAT_ZERO_OR_MORE )
       return false;

   int previous;
   while (*pattern)
      {
      switch (*pattern)
         {
         case PAT_ANY :
            compiledPattern.Push(previous = PM_ANY);
            pattern++;
            break;
         case PAT_BEGIN_MASK :
            {
            compiledPattern.Push(previous = PM_MASK);
            int offset = compiledPattern.Length();
            compiledPattern.Push(0);
            for (unsigned int i = 1; i < MAP_SIZE; i++)
               compiledPattern.Push(0);
            if (!CompileMask(&compiledPattern[offset], ++pattern))
               {
               compiledPattern.Clear();
               return false;
               }
            if (*(pattern++) != PAT_END_MASK)
               {
               compiledPattern.Clear();
               return false;
               }
            break;
            }
         case PAT_OPTIONAL :
         case PAT_ONE_OR_MORE :
         case PAT_ZERO_OR_MORE :
            {
            switch (previous)
               {
               case PM_OPTIONAL :
               case PM_ONE_OR_MORE :
               case PM_ZERO_OR_MORE :
                  compiledPattern.Clear();
                  return false;
               }

            int offset = previous == PM_MASK ? MAP_SIZE + 1 : 1;
            compiledPattern.InsertAt(compiledPattern.Length() - offset,
               *pattern == PAT_OPTIONAL ? PM_OPTIONAL :
               *pattern == PAT_ONE_OR_MORE ? PM_ONE_OR_MORE : PM_ZERO_OR_MORE);
            pattern++;
            break;
            }
         default :
            compiledPattern.Push(previous = CompileCharacter(pattern));
         }
      }

   compiledPattern.Push(0);
   return true;
   }

int PatternMatcher::ComparePattern(const char * string)
   {
   if (compiledPattern.Length() == 0)
      return 0;

   const int * pattern = &compiledPattern[0];

   return ComparePattern(pattern, string);
   }

int PatternMatcher::ComparePattern(const int * pattern, const char * string)
   {
   const char * start = string;

   while (*pattern)
      {
      if (*pattern == PM_OPTIONAL)
         {
         // Match is optional, so no failure for mismatch
         if (CompareCharacter(pattern + 1, string))
            string++;
         AdvancePointer(pattern);
         AdvancePointer(pattern);
         }
      else if (*pattern == PM_ONE_OR_MORE || *pattern == PM_ZERO_OR_MORE)
         {
         // Check for required match
         if (*pattern == PM_ONE_OR_MORE)
            if (CompareCharacter(pattern + 1, string))
               string++;
            else
               return 0;

         // Additional matches are optional, so no failure for mismatch
         const char * patternStart = string;
         while ( *string && CompareCharacter(pattern + 1, string))
            string++;

         AdvancePointer(pattern);
         AdvancePointer(pattern);

         if (*pattern != 0)
            {
            for (int match = 0 ; string >= patternStart; string--)
               if ((match = ComparePattern(pattern, string)) >  0)
                  return string - start + match;
            return 0;
            }
         }
      else if (CompareCharacter(pattern, string))
         {
         string++;
         AdvancePointer(pattern);
         }
      else
         return 0;
      }
   return string - start;
   }

bool PatternMatcher::CompareCharacter(const int * pattern, const char * string)
   {
   switch (*pattern)
      {
      case PM_ANY :
         return (*string != 0 && *string != '\n');
      case PM_MASK :
         return GetMask(pattern + 1, *string);
      default :
         return *pattern == (int) *string;
      }
   }

bool PatternMatcher::CompileMask(int * mask, const char * & definition)
   {
   bool reverse_mask = false;
   if (*definition == PAT_REVERSE_MASK)
      {
      reverse_mask = true;
      definition++;
      }

   int previous;
   int previousIsValid = false;
   while (*definition && *definition != PAT_END_MASK)
      {
      if (*definition != '-' || definition[1] == PAT_END_MASK || definition[1] == 0 || !previousIsValid)
         {
         int character = CompileCharacter(definition);

         SetMask(mask, character);
         previous = character;
         previousIsValid = true;
         }
      else
         {
         int end = CompileCharacter(++definition);

         if (previous < end)
            while (previous++ != end)
               SetMask(mask, previous);
         else
            while (previous-- != end)
               SetMask(mask, previous);

         previousIsValid = false;
         }
      }

   if (reverse_mask)
      for (int i = 0; i < 4; i++)
         mask[i] ^= 0xFFFF;

   return true;
   }


bool PatternMatcher::isHexDigit(char ch)
   {
   if (ch >= '0' && ch <= '9' ||
       ch >= 'a' && ch <= 'f' ||
       ch >= 'A' && ch <= 'F')
       return true;
   return false;
   }

int PatternMatcher::convertHexDigit(char ch)
   {
   if (ch >= '0' && ch <= '9') return ch - '0';
   if (ch >= 'a' && ch <= 'f') return ch - 'a' + 10;
   if (ch >= 'A' && ch <= 'F') return ch - 'A' + 10;
   return 0;
   }

void PatternMatcher::SetMask(int * mask, int bit)
   {
   mask[bit / BITS_PER_INT] |= 1 << (bit % BITS_PER_INT);
   }

bool PatternMatcher::GetMask(const int * mask, int bit)
   {
   return mask[bit / BITS_PER_INT] & (1 << (bit % BITS_PER_INT));
   }

int PatternMatcher::CompileCharacter(const char * & definition)
   {
   int character = *(definition++);

   if (character != '\\')
      return character;

   character = *(definition++);
   switch (character)
      {
      case 0 :
         definition--;
         return '\\';
      case 'a' :
         return '\a';
      case 'b' :
         return '\b';
      case 'f' :
         return '\f';
      case 'n' :
         return '\n';
      case 'r' :
         return '\r';
      case 't' :
         return '\t';
      case 'v' :
         return '\v';
      case 'x' :
         if (isHexDigit(*definition)) character = convertHexDigit(*definition++);
         if (isHexDigit(*definition)) character = character * 16 +
                                                  convertHexDigit(*definition++);
         return character;
      case '^' :
         if (isalpha(*definition))
            return toupper(*definition++) - 'A' + 1;
      }

   return character;
   }

void PatternMatcher::AdvancePointer(const int * & pattern)
   {
   pattern += (*pattern == PM_MASK) ? MAP_SIZE + 1 : 1;
   }

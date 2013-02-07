#include "StringAlias.h"
#include "InputFile.h"

void StringAlias::SetAlias(String & string, String & alias)
   {
   int index = lookup.Integer(string);

   if (index < 0)
      {
      aliases.Push(alias);
      lookup.SetInteger(string, aliases.Length() - 1);
      }
   else
      aliases[index] = alias;
   }

const String & StringAlias::GetAlias(const String & string) const
   {
   int index = lookup.Integer(string);

   if (index < 0)
      return string;
   else
      return aliases[index];
   }

bool StringAlias::ReadFromFile(const char * filename)
   {
   IFILE input = ifopen(filename, "rt");

   if (input == NULL)
      return false;

   ReadFromFile(input);

   ifclose(input);

   return true;
   }

bool StringAlias::ReadFromFile(IFILE & input)
   {
   StringArray lines, tokens;
   lines.Read(input);

   for (int j = 0; j < lines.Length(); j++)
      {
      tokens.ReplaceTokens(lines[j]);

      if (tokens.Length() != 2) continue;

      SetAlias(tokens[0], tokens[1]);
      }

   return true;
   }


#include "FastaLike.h"

bool FastaLike::OpenArchive(const char * filename)
   {
   CloseArchive();

   input = ifopen(filename, "rb");

   return input != NULL;
   }

void FastaLike::CloseArchive()
   {
   if (input != NULL) ifclose(input);

   buffer.Clear();
   input = NULL;
   }

int FastaLike::Load(IFILE & input, String & shortLabel, String & label, String & sequence)
   {
   sequence.Clear();

   if (buffer[0] == '>')
      {
      label = buffer;
      buffer.Clear();
      }
   else
      label.ReadLine(input);

   int delim = label.FindChar(' ');
   if (delim > 0)
      shortLabel = label.SubStr(1, delim - 1);
   else
      shortLabel = label.SubStr(1);

   if (input == NULL || ifeof(input))
      return -1;

   if (label[0] != '>' && label[0] != '@')
      printf("WARNING: Invalid FAST[A/Q] sequence header (expecting line to begin with '>' or '@')\n");

   while (!ifeof(input))
      {
      buffer.ReadLine(input);
      // printf("BASES: %s\n", (const char *) buffer);

      if (buffer[0] == '>') break;

      int length = sequence.Length();
      for (int i = 0; i < buffer.Length(); i++)
         if (buffer[i] != ' ' || buffer[i] != '\t')
         {
         sequence += (char) buffer[i];
         length++;
         }

      sequence.SetLength(length);
      }

   return sequence.Length();
   }

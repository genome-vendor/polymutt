#include "MemoryInfo.h"

String & IntegerInfo(double integer)
   {
   static String info;

   if (integer < 10000)
      {
      info.printf("%.0f", integer);
      return info;
      }

   if (integer < 1000. * 1000.)
      info.printf("%.1f K", (integer + 999) / 1000.);
   else if (integer < 1024. * 1024. * 1024.)
      info.printf("%.1f M", (integer + 1000. * 1000. - 1) / (1000. * 1000.));
   else if (integer < 1000. * 1000. * 1000. * 1000.)
      info.printf("%.1f G", integer / (1000. * 1000. * 1000.));
   else
      info.printf("%.1f T", integer / (1000. * 1000. * 1000. * 1000.));

   return info;
   }

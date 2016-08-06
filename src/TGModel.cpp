#include "TGlib.h"
int main(int argc, char* arg[])
{
  BigTherm Therm(arg[1]);
  Therm.TakataDamage();
  return 0;
}

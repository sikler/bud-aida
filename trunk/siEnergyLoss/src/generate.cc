#include "../interface/ClusterGenerator.h"

#include <cstdlib>
#include <cstring>

/*****************************************************************************/
int main(int arg, char *arc[])
{
  ClusterGenerator theGenerator;

  if(strcmp(arc[3],"-pixel") == 0)
    theGenerator.run(atoi(arc[2]), 0);

  if(strcmp(arc[3],"-strip") == 0)
    theGenerator.run(atoi(arc[2]), 1);

  return 0;
}

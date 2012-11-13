#include <algorithm>
#include "ZShape/HFZeeVBTF/interface/SelectElectron.h"

void selectElectron::add(int v)
{ allowed.push_back(v);}

bool selectElectron::doesElePass(int v)
{ 
  if ( find(allowed.begin(),allowed.end(),v) != allowed.end() ) return true;
  else return false;
}

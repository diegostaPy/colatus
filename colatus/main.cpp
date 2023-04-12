
//
//  main.cpp
//  NBodySim
//
//  Created by José Ricardo da Silva  Júnior on 06/02/12.
//  Modified by Diego H. Stalder Diaz   18/09/12
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "Direct.h"

int main (int argc,  char * argv[])
{
		NBody::Simulators::Direct direct(argv[1]);
		direct.Go();

	return 0;
}


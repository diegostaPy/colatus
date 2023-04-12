/*
 * main.cpp
 *
 *  Created on: 26/01/2011
 *      Author: josericardodasilvajunior
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <libxml/xmlwriter.h>
#include <libxml/xmlstring.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <sstream>
#include "inc/AppMain.h"
#include <helper_cuda.h>
#include "jrfxgl.h"
int fromChar(char* value){
	int i;
	std::istringstream is((char*)value);
	is >> i;
	return i;
}
double fromCharF(char* value){
	double i;
	std::istringstream is((char*)value);
	is >> i;

	return i;
}

params read_parameters(const char *filename){
// Ler dados do xml
	params P;std::string parameter;
	std::string _xmlFilename = std::string(filename);
	std::string _mDataFilename;
	int flag;
	xmlDocPtr doc;
	xmlNodePtr cur;
	doc = xmlParseFile(_xmlFilename.c_str());

	if (doc == NULL){
		fprintf(stderr, "Error reading file!");
	}

	cur = xmlDocGetRootElement(doc);
	if (cur == NULL){
		fprintf(stderr, "File is empty!");
		xmlFreeDoc(doc);
	}

	// Read data
	if (!xmlStrcmp(cur->name, (const xmlChar*) "NBody_Sim_out")){
		xmlChar* r_buf;

	    r_buf= xmlGetProp(cur, (const xmlChar*) "NumBodies");
	    P.numbodies= fromChar((char*) r_buf);
		xmlFree(r_buf);

	   	r_buf= xmlGetProp(cur, (const xmlChar*) "Box");
	   	P.L= (float)fromChar((char*) r_buf);
		xmlFree(r_buf);
		r_buf= xmlGetProp(cur, (const xmlChar*) "Output_base_filename");
		P.filename=std::string((char*)r_buf);
		xmlFree(r_buf);
		r_buf= xmlGetProp(cur, (const xmlChar*) "Snapshots");
		P.snapshots= fromChar((char*) r_buf);
		xmlFree(r_buf);
		r_buf=xmlGetProp(cur, (const xmlChar*) "File_binary");
		parameter=std::string((char*)r_buf);
		P.BinInput=((parameter.compare("Yes"))== 0);
		parameter.clear();
		xmlFree(r_buf);
	}

	xmlFree(doc);
	return P;
	}
int main(int argc, const  char *argv[])
{
printf("\n Filename %s Read ouput files ...",argv[1]);
	params P=read_parameters(argv[1]);
printf("...done");
                AppMain app(P);
                app.go();


    

    return 0;
}


#ifndef _OPTIONS_H_
#define _OPTIONS_H_

#include <stdlib.h>
#include <string>
#include <cstring>
#include <stdio.h>

using namespace std;

extern bool ILSA;
extern bool NOISOLATE;
<<<<<<< .merge_file_N36WcI
extern bool USERDATA;
extern bool PARAMS;
=======
//extern bool USERDATA;
//extern bool PARAMS;
extern bool LIMIT_DISTANCE;
>>>>>>> .merge_file_qTispI
extern bool BPP_ENABLED;
extern bool SUBOPT_ENABLED;
extern bool CONS_ENABLED;
extern bool VERBOSE;
<<<<<<< .merge_file_N36WcI
extern bool SHAPE_ENABLED;
=======
extern bool PARAM_DIR;
>>>>>>> .merge_file_qTispI

extern string seqfile;
extern string constraintsFile;
extern string shapeFile;
extern string outputFile;
extern string paramDir;

extern int suboptDelta;
extern int nThreads;

extern bool LIMIT_DISTANCE;
extern int contactDistance;


void help();
void parse_options(int argc, char** argv);
void printRunConfiguration(string seq);


#endif

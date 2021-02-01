#ifndef MODELSELECTOR_H
#define MODELSELECTOR_H

#include "prefixspan.h"
#include "learnerspp.h"
#include "database.h"

#include <vector>

using namespace std;

class ModelSelector{
public:
  virtual ~ModelSelector() {}
  virtual uint select(const vector<double> &aSolutionPath,const Database &aDatabase,PrefixSpan &aPrefix,LearnerSPP &aLearner,const vector<uint> &aOptions, PrefixSpan &tPrefix, double &answer) = 0;
  virtual uint select_with_parameter(const vector<double> &aSolutionPath,const Database &aDatabase,vector<PrefixSpan> &aPrefixs,LearnerSPP &aLearner,const vector<uint> &aOptions, PrefixSpan &tPrefix,uint &select_prefixs_i) = 0;
};

#endif

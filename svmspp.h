#ifndef SVMSPP_H
#define SVMSPP_H

#include "prefixspan.h"
#include "database.h"
#include "learnerspp.h"
#include <chrono>
#include <vector>
#include <unordered_map>

using namespace std;

//L1-reguralized L2-SVM
class SVMSPP : public LearnerSPP{
private:

  uint mN;

  double mBias;
  vector<double> mR;
  uint mT;
  uint mMaxIter;
  uint mFreq;
  double mEps;
  double mRatio;
  double clac_sup(const vector<Event> &aSequence,const vector<Event> &aPattern,const uint aSupMode);

  /*
  struct keyPattern{
    vector<Event> Transac;
    vector<Event> pat;
    bool operator < (const keyPattern& k) const{
      if(sizeof(Transac) < sizeof(k.Transac)){
        return true;
      }else if(sizeof(Transac) > sizeof(k.Transac)){
        return false;
      }else if(sizeof(pat) < sizeof(k.pat)){
        return true;
      }else if(sizeof(pat) > sizeof(k.pat)){
        return false;
      }
      return false;
    }
  };


  //pattern と系列に応じたsupportを保存するmapを作成
  map<keyPattern, double> mSupportMap;
  */
  map<pair<vector<Event>, vector<Event>>, double > mSupportMap;
  //unordered_map <pair<vector<Event>, vector<Event>>, double > mSupportMap;

public:
  SVMSPP(uint aMaxIter, uint aFreq, double aEps,double aRatio){
  mMaxIter = aMaxIter;
  mFreq = aFreq;
  mEps = aEps;
  mRatio = aRatio;
  };
  //aOption[0]に1~aLambdas.size()までの数字をi入れることでaLambdas[i-1]での学習を行う
  virtual void learn(PrefixSpan &aPrefix,const vector<double> &aLambdas,const vector<uint> &aOptions);
  //λmaxの値を計算して返す
  virtual double get_lambda_max(PrefixSpan &aPrefix);
  //現在学習されているモデルで予測を行う
  virtual vector<double> predict(const PrefixSpan &aPrefix,const vector<vector<Event> > &aTransaction, const vector<uint> &testId,PrefixSpan &tPrefix);
  //solution path全てででの予測値Y^を返す. aOption[0]の設定learnと同じ
  virtual vector<vector<double> > get_all_predict(PrefixSpan &aPrefix,const vector<double> &aLambdas,const vector<vector<Event> > &aTransaction,const vector<uint> &aOptions, const vector<uint> &testId,PrefixSpan &tPrefix);

  virtual double get_bias(void);
};

#endif

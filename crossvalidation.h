#ifndef CROSSVALIDATION_H
#define CROSSVALIDATION_H

#include "modelselector.h"
#include <random>

//learnerspp
#include "svmspp.h"
#include "lassospp.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

class CrossValidation: public ModelSelector{
private:
	vector<uint> mTrainIdx;
	vector<uint> mTestIdx;
	vector<uint> mNumDiv;


public:
	CrossValidation(){};
	void calc_accuracy(const vector<vector<double> > &aYHats, const vector<double> &aY, vector<double> &aCnt,const uint &aLossType);
	void next_train_test();

	//テンプレートでaLearnerの型を抽象化してます． aLearnerの型によって処理が変わるなら具体的な型名をつけたものを実装としてかけばオーバーロードされます．

	template<class T> uint select(const vector<double> &aSolutionPath, const Database &aDatabase, PrefixSpan &aPrefix, T &aLearner, const vector<uint> &aOptions, PrefixSpan &tPrefix, double &answer){
		uint tK = aOptions[0];
		uint tAve = aOptions[1];
		uint tLossType = aOptions[2];

		//全データからテストデータと学習データのk分割をこしらえる
		vector<vector<Event> > tTransaction = aDatabase.get_transaction();
		vector<double> tY = aDatabase.get_y();
		uint tN = tY.size();
		if(tN < tK){
			cout << "error:CV,n<k" << '\n';
			exit(1);
		}

		//各λの正解数を保持
		vector<double> tCnt(aSolutionPath.size() - 1, 0);

		//1.データがk分割すると何個づつになるかを計算する
		//例えば 19 = [4,5,5,5]
		vector<uint> tNumDiv;
		for(uint j = 0; j < tK; ++j){
			tNumDiv.push_back((tN + j) / tK);
		}

		#pragma omp parallel
		{
			//各λの正解数をスレッド毎に保持
			vector<double> tCnt_private(aSolutionPath.size() - 1, 0);

			#pragma omp for
			for(uint i = 0; i < tAve; ++i){
				//inex をシャッフル
				vector<uint> tTrainIdx;
				for(uint j = 0; j < tN; ++j){
					tTrainIdx.push_back(j);
				}

				for(uint j = 0; j < tN; ++j){
					uint k = j + (rand() % (tN - j));
					swap(tTrainIdx[j], tTrainIdx[k]);
				}

				//k-foldの結果がtCntに入る
				k_fold(aSolutionPath, aPrefix, aLearner, tTransaction, tY, tK, tCnt_private, tTrainIdx, tNumDiv,tLossType,tPrefix);
			}

			//tCntへ結果をまとめる．念のため排他的に処理
			#pragma omp critical
			{
				for(uint j = 0; j < tCnt.size(); ++j){
					tCnt[j] += tCnt_private[j];
				}
			}

		}

		vector<double>::iterator tIter;
		cout << '\n';
		switch (tLossType) {
      // L2SVM
			case 1:
			for(uint i = 0; i < tCnt.size(); ++i){
				cout << "λ[" << i + 1 << "]=" << aSolutionPath[i + 1] << " " << tCnt[i] << "/" << aOptions[1] * tN << "=" << tCnt[i] / (double) (aOptions[1] * tN) << '\n';
			}
			//最大値選び
			tIter = max_element(tCnt.begin(), tCnt.end());
			break;
      //Lasso
			case 2:
			for(uint i = 0; i < tCnt.size(); ++i){
				cout << "λ[" << i + 1 << "]=" << aSolutionPath[i + 1] << " RMSE:" <<sqrt(tCnt[i] / (double) (aOptions[1] * tN)) << '\n';
			}
			//最大値選び
			tIter = min_element(tCnt.begin(), tCnt.end());
			break;
			default:
			std::cout << "error:CV output" << '\n';
		}




		return distance(tCnt.begin(), tIter) + 1;
	}

	template<class T> uint select_with_parameter(const vector<double> &aSolutionPath, const Database &aDatabase, vector<PrefixSpan> &aPrefixs, T &aLearner, const vector<uint> &aOptions, PrefixSpan &tPrefix, uint &select_prefixs_i){
		uint tK = aOptions[0];
		uint tAve = aOptions[1];
		uint tLossType = aOptions[2];

		//全データからテストデータと学習データのk分割をこしらえる
		vector<vector<Event> > tTransaction = aDatabase.get_transaction();
		vector<double> tY = aDatabase.get_y();
		uint tN = tY.size();
		if(tN < tK){
			cout << "error:CV,n<k" << '\n';
			exit(1);
		}

		//各λの正解数を保持
		vector<double> tCnt(aSolutionPath.size() - 1, 0);
		vector<vector<double>> tCnts(aPrefixs.size(),vector<double>(aSolutionPath.size() - 1, 0));

		//1.データがk分割すると何個づつになるかを計算する
		//例えば 19 = [4,5,5,5]
		vector<uint> tNumDiv;
		for(uint j = 0; j < tK; ++j){
			tNumDiv.push_back((tN + j) / tK);
		}

		#pragma omp parallel
		{
			//各λの正解数をスレッド毎に保持
			vector<double> tCnt_private(aSolutionPath.size() - 1, 0);
			vector<vector<double>> tCnts_private(aPrefixs.size(),vector<double>(aSolutionPath.size() - 1, 0));

			#pragma omp for
			for(uint i = 0; i < tAve; ++i){
				//inex をシャッフル
				vector<uint> tTrainIdx;
				for(uint j = 0; j < tN; ++j){
					tTrainIdx.push_back(j);
				}

				for(uint j = 0; j < tN; ++j){
					uint k = j + (rand() % (tN - j));
					swap(tTrainIdx[j], tTrainIdx[k]);
				}

				//k-foldの結果がtCntに入る
				for(uint j = 0; j < aPrefixs.size();j++){
					k_fold(aSolutionPath, aPrefixs[j], aLearner, tTransaction, tY, tK, tCnt_private, tTrainIdx, tNumDiv,tLossType,tPrefix);
					tCnts_private[j] = tCnt_private;
			 	}
			}

			//tCntへ結果をまとめる．念のため排他的に処理
			#pragma omp critical
			{
				for(uint j = 0; j < tCnts.size(); ++j){
					for(uint k = 0; k < tCnts[j].size(); ++k){
						tCnts[j][k] += tCnts_private[j][k];
					}
				}
			}

		}

		vector<double>::iterator tIter;
		double Max = 0;
		double Min = 1000;
		uint answer;
		cout << '\n';
		switch (tLossType) {
			// L2SVM
			case 1:
			for(uint i = 0; i < tCnts.size(); ++i){
				cout << "aPrefix["<< i << "]" << endl;
				for(uint j = 0; j < tCnts[i].size(); ++j){
					cout << "λ[" << j + 1 << "]=" << aSolutionPath[j + 1] << " " << tCnts[i][j] << "/" << aOptions[1] * tN << "=" << tCnts[i][j] / (double) (aOptions[1] * tN) << '\n';
				}
			}
			//最大値選び
			for (uint i = 0; i < tCnts.size(); ++i){
				double max_tCnt = max(max_tCnt,(double)*max_element(tCnts[i].begin(), tCnts[i].end()));
				if(Max < max_tCnt){
					Max = max_tCnt;
					select_prefixs_i = i;
					answer = distance(tCnts[i].begin(), max_element(tCnts[i].begin(), tCnts[i].end())) + 1;
				}
			}
			break;
			//Lasso
			case 2:
			for(uint i = 0; i < tCnts.size(); ++i){
				cout << "aPrefix["<< i << "]" << endl;
				for(uint j = 0; j < tCnts[i].size(); ++j){
					cout << "λ[" << j + 1 << "]=" << aSolutionPath[j + 1] << " RMSE:" <<sqrt(tCnts[i][j] / (double) (aOptions[1] * tN)) << '\n';
				}
			}
			//最小値選び
			for (uint i = 0; i < tCnts.size(); ++i){
				double min_tCnt = min(min_tCnt,(double)*min_element(tCnts[i].begin(), tCnts[i].end()));
				if(Min < min_tCnt){
					Min = min_tCnt;
					select_prefixs_i = i;
					answer = distance(tCnts[i].begin(), min_element(tCnts[i].begin(), tCnts[i].end())) + 1;
				}
			}
			break;
			default:
			std::cout << "error:CV output" << '\n';
		}




		return answer;
	}

	template<class T> void k_fold(const vector<double> &aSolutionPath, PrefixSpan &aPrefix, T &aLearner, const vector<vector<Event>> &aTransaction, const vector<double> &aY, const uint &aK, vector<double> &aCnt, vector<uint> aTrainIdx, vector<uint> aNumDiv,const uint &aLossType,PrefixSpan &tPrefix){

		vector<uint> tSVMOption = { (uint) aSolutionPath.size() };

		#pragma omp parallel
		{
			//aPrefix, aLearnerはそれぞれのスレッドでprivateに持つようにする
			PrefixSpan tPrefix_private = aPrefix;

			//スレッドごとにインスタンスのコピーを持ちたいので派生クラスTで初期化
			T tLearner_private = aLearner;

			#pragma omp for
			for(uint i = 0; i < aK; ++i){
				vector<uint> tTestIdx;
				vector<uint> tTrainIdx = aTrainIdx;
				uint tSumIdx = 0;
				for(uint j = 0; j < i; ++j){
					tSumIdx += aNumDiv[j];
				}

				copy(tTrainIdx.begin() + tSumIdx, tTrainIdx.begin() + tSumIdx + aNumDiv[i], back_inserter(tTestIdx));
				tTrainIdx.erase(tTrainIdx.begin() + tSumIdx, tTrainIdx.begin() + tSumIdx + aNumDiv[i]);

				vector<vector<Event> > tTrainTransaction(tTrainIdx.size());
				vector<double> tTrainY(tTrainIdx.size());
				for(uint j = 0; j < tTrainIdx.size(); ++j){
					tTrainTransaction[j] = aTransaction[tTrainIdx[j]];
					tTrainY[j] = aY[tTrainIdx[j]];
				}

				vector<vector<Event>> tTestTransaction(tTestIdx.size());
				vector<double> tTestY(tTestIdx.size());
				for(uint j = 0; j < tTestIdx.size(); ++j){
					tTestTransaction[j] = aTransaction[tTestIdx[j]];
					tTestY[j] = aY[tTestIdx[j]];
				}

				tPrefix_private.init(tTrainTransaction, tTrainY);

				vector<vector<double> > tYHats = tLearner_private.get_all_predict(tPrefix_private, aSolutionPath, tTestTransaction, tSVMOption, tTestIdx,tPrefix);

				calc_accuracy(tYHats, tTestY, aCnt,aLossType);
			}
		}

	}

	void k_fold(const vector<double> &aSolutionPath, PrefixSpan &aPrefix, LearnerSPP &aLearner, const vector<vector<Event> > &aTransaction, const vector<double> &aY, const uint &aK, vector<double> &aCnt,const uint &aLossType,PrefixSpan &tPrefix);

	//option[k-foldのk,何回の平均をとるか,...]
	uint select(const vector<double> &aSolutionPath, const Database &aDatabase, PrefixSpan &aPrefix, LearnerSPP &aLearner, const vector<uint> &aOptionsno,PrefixSpan &tPrefix, double &answer);

	uint select_with_parameter(const vector<double> &aSolutionPath, const Database &aDatabase, vector<PrefixSpan> &aPrefixs, LearnerSPP &aLearner, const vector<uint> &aOptions, PrefixSpan &tPrefix, uint &select_prefixs_i);
};

#endif

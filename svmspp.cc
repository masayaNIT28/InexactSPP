#include "svmspp.h"

double SVMSPP::get_lambda_max(PrefixSpan &aPrefix){
//cout << "get_lambda_max" << endl;
	mN = aPrefix.mY.size();
	mBias = 0;
	for(uint i = 0; i < mN; ++i){
		mBias += aPrefix.mY[i];
	}
	mBias /= mN;
	//cout << "get_lambda_max" << endl;

	mR.shrink_to_fit();
	mR.resize(mN);
	//cout << "get_lambda_max" << endl;
	for(uint i = 0; i < mN; ++i){
		//損失
		mR[i] = 1 - aPrefix.mY[i] * mBias;
	}
//cout << "get_lambda_max" << endl;
	aPrefix.get_maxnode(mR, 1);
	//cout << "get_lambda_max" << endl;

	return aPrefix.mMaxnode.val;
}
//sppは木構造を利用しているのでデータではなくaPrefixを入力する
void SVMSPP::learn(PrefixSpan &aPrefix, const vector<double> &aLambdas, const vector<uint> &aOptions){

	if(aOptions.empty()){
		cout << "error:svm.learn option size is incorrect." << '\n';
		exit(1);
	}
	//λを求めるための処理
	if(aLambdas.empty()){
		cout << "error:lambda is empty. Set lambda." << '\n';
		exit(1);
	}
	mT = aLambdas.size();
	double tTrueLambdaMax = get_lambda_max(aPrefix);

	// compute solution path
	uint activesize = 1;
	double L1norm = 0;
	double maxval = tTrueLambdaMax; //lambda_max
	for(uint t = 1; (t < aOptions[0]) && (t < mT); ++t){

		double lam = aLambdas[t];

		cout << "aLambdas["<< t << "]=" << lam << endl;

		if(lam > tTrueLambdaMax){
			cout << "skip:λ=" << lam << " > λmax=" << tTrueLambdaMax << '\n';
			continue;
		}

		vector<uint> index;
		vector<uint *> x;
		vector<double *> support;
		// vector<uint> norm;
		vector<double> norm;
		vector<uint> xSize;
		vector<double> grad;
		vector<double> w;
		for(uint iter = 0; iter <= mMaxIter; ++iter){
			// calculate dual and safe screening
			if(iter % mFreq == 0){
				double loss = 0;
				double oTr = 0;
				for(uint i = 0; i < mN; ++i){
					if(mR[i] > 0){
						loss += mR[i] * mR[i];
						oTr += mR[i];
					}
				}

				if(iter >= 1){
					maxval = 0;
					for(uint s = 0; s < activesize; ++s){

						uint j = index[s];
						grad[j] = 0;

						// for (uint k = 0; k < norm[j]; ++k) {
						for(uint k = 0; k < xSize[j]; ++k){
							uint id = x[j][k];
							if(mR[id] > 0){
								grad[j] += aPrefix.mY[id] * mR[id] * support[j][k];
								//@
							}
						}
						if(fabs(grad[j]) > maxval) maxval = fabs(grad[j]);
					}
				}

				double alpha = min(max(oTr / (lam * loss), 0.0), 1 / maxval);
				double primal = 0.5 * loss + lam * L1norm;
				double dual = -0.5 * lam * lam * alpha * alpha * loss + lam * alpha * oTr;
				double gap = primal - dual;

				if(gap / primal < mEps){
					if(index.empty()){
						break;
					}
					uint active = 0;
					for(uint s = 0; s < activesize; s++){
						uint j = index[s];
						aPrefix.mActive[j]->w = w[j];
						if(w[j] != 0) active++;
					}

					//cout << "[iter " << iter << "] primal: " << primal << ", dual: " << dual << ", gap: " << gap/primal << ", activeset: " << active << endl;
					break;
				}

				double radius = sqrt(2 * gap) / lam;
				if(iter == 0){
					//各λループのはじめに一回だけscreeningを行う(現在のラムダでの最大の木を構築している)
					//λは徐々に小さくなって行くの木は前回のλより恐らく大きくなる
					//w更新するたびバウンドがきつくなり，ある周期で非アクティブを決め，木が小さくなって行く
					aPrefix.safe_screening(mR, alpha, radius, 1);
					activesize = aPrefix.mActive.size();

					index.resize(activesize);
					x.resize(activesize);
					norm.resize(activesize);
					grad.resize(activesize);
					w.resize(activesize);
					xSize.resize(activesize);
					support.resize(activesize);

					for(uint j = 0; j < activesize; ++j){
						index[j] = j;
						x[j] = aPrefix.mActive[j]->x.data();
						support[j] = aPrefix.mActive[j]->support.data();
						// norm[j] = aPrefix.Active[j]->x.size();
						xSize[j] = aPrefix.mActive[j]->x.size();
						norm[j] = aPrefix.mActive[j]->supportSumsqr;
						grad[j] = 0;
						w[j] = aPrefix.mActive[j]->w;
					}

				}else{
					//最初の１回目以外は非アクティブノード検出を行う．（木の大きさは変わらない）
					for(uint s = 0; s < activesize; ++s){
						uint j = index[s];
						double score = fabs(alpha * grad[j]) + radius * sqrt(norm[j]);

						if(score < 1){
							if(w[j] != 0){
								// for (uint k = 0; k < norm[j]; ++k) {
								for(uint k = 0; k < xSize[j]; ++k){
									uint id = x[j][k];
									mR[id] += aPrefix.mY[id] * w[j] * support[j][k];  				//@
								}
							}
							//学習中にwが0になってもActiveのサイズは変更していない
							aPrefix.mActive[j]->w = 0;
							activesize--;
							swap(index[s], index[activesize]);
							s--;
						}
					}
				}
				// cout << "[iter " << iter << "] primal: " << primal << ", dual: " << dual << ", gap: " << gap/primal << ", activeset: " << activesize << endl;
			}

			//ランダムに混ぜている
			for(uint j = 0; j < activesize; ++j){
				uint i = j + rand() % (activesize - j);
				swap(index[i], index[j]);
			}

			// update w
			//newton法
			L1norm = 0;
			for(uint s = 0; s < activesize; ++s){
				uint j = index[s];

				double G = 0;
				double H = 0;
				double loss_old = lam * fabs(w[j]);  		//@
				// for (uint k = 0; k < norm[j]; ++k) {
				for(uint k = 0; k < xSize[j]; ++k){
					uint idx = x[j][k];
					if(mR[idx] > 0){
						G -= aPrefix.mY[idx] * mR[idx] * support[j][k];
						H += support[j][k] * support[j][k];
						loss_old += 0.5 * mR[idx] * mR[idx];
					}
				}
				H = max(H, 1e-12);

				double Gp = G + lam;
				double Gn = G - lam;
				double d;
				if(Gp <= H * w[j]) d = -Gp / H;
				else if(Gn >= H * w[j]) d = -Gn / H;
				else d = -w[j];

				if(fabs(d) < 1e-12){
					L1norm += fabs(w[j]);
					continue;
				}

				double d_old = 0;

				double bound = 0.5 * (G * d + lam * fabs(w[j] + d) - lam * fabs(w[j]));
				uint linesearch;
				for(linesearch = 0; linesearch < 10; ++linesearch){
					double loss_new = lam * fabs(w[j] + d);
					// for (uint k = 0; k < norm[j]; ++k) {
					for(uint k = 0; k < xSize[j]; ++k){
						uint idx = x[j][k];
						//損失の入れ替え
						mR[idx] += aPrefix.mY[idx] * (d_old - d) * support[j][k];
						if(mR[idx] > 0){
							loss_new += 0.5 * mR[idx] * mR[idx];
						}
					}

					if(loss_new - loss_old <= bound){
						break;
					}else{
						d_old = d;
						d *= 0.5;
						bound *= 0.5;
					}
				}
				w[j] += d;
				L1norm += fabs(w[j]);
			}

			// bias
			uint nn = 0;
			double tmp = 0;
			for(uint i = 0; i < mN; i++){
				if(mR[i] > 0){
					tmp += aPrefix.mY[i] * mR[i] + mBias;
					nn++;
				}
			}
			double bias_old = mBias;
			mBias = tmp / nn;
			for(uint i = 0; i < mN; i++){

				mR[i] += aPrefix.mY[i] * (bias_old - mBias);
			}
		}

		aPrefix.add_history(t);

	}

}

vector<double> SVMSPP::predict(const PrefixSpan &aPrefix, const vector<vector<Event> > &aTransaction, const vector<uint> &testId,PrefixSpan &tPrefix){
	uint tTransactionSize = aTransaction.size();
	vector<double> tYHat(tTransactionSize, mBias);  			//要素を全てBiasで初期化


	//keyPattern tKey;

	//for(uint i : testId){
		//cout << i << endl;
	//}



//分割する前の木からActive patternのNodeを見つけ, testのsupportを見つけてくる

//cout << "aa" << endl;
	for(auto tActiveNode: aPrefix.mActive){
		//cout << "bb" << endl;
		vector<Event> tPattern = tActiveNode->pattern;
		double tW = tActiveNode->w;
		if(tW != 0){
			pair<vector<Event>, vector<Event>> tKey;
			tKey.first = aTransaction[0]; //どれでもいいけどpatternが含まれてればいい
			tKey.second = tPattern;
			if(mSupportMap.find(tKey) == mSupportMap.end()){
			  //cout << "cc" << endl;
				//cout << tActiveNode->patternStr << endl;
				//Node見つけてきてそのsupportのリストを吐き出す
				vector<double> support = tPrefix.calcActivesupport(tPattern, testId);
				//cout << "support size = " << support.size() << endl;
				//mapに入れる
				for(int i = 0; i < testId.size(); i++){
					//cout << "dd" << endl;
					//cout << "testId[" << i << "] = "<< testId[i] << endl;
					//cout << "support[" << i << "] = " << support[i] << endl;
					tKey.first = aTransaction[i];
					mSupportMap.emplace(tKey,support[i]);
				}
			}
		}
	}
	//cout << "aa" << endl;

	cout << mSupportMap.size() << endl;

	for(auto tActiveNode : aPrefix.mActive){
		vector<Event> tPattern = tActiveNode->pattern;
		double tW = tActiveNode->w;
		if(tW != 0){
			pair<vector<Event>, vector<Event>> tKey;
			tKey.second = tPattern;
			for(uint i = 0; i < tTransactionSize; ++i){
				tKey.first = aTransaction[i];
				tYHat[i] += tW * mSupportMap[tKey];
			}
		}
	}

	/*
	for(auto tActiveNode : aPrefix.mActive){
		vector<Event> tPattern = tActiveNode->pattern;
		double tW = tActiveNode->w;

		for(uint i = 0; i < tTransactionSize; ++i){
			//tKey = {aTransaction[i], tPattern};
			pair<vector<Event>, vector<Event>> tKey;
			tKey.first = aTransaction[i];
			tKey.second = tPattern;
			double tmp;
			if(mSupportMap.find(tKey) == mSupportMap.end()){
				if(aPrefix.calcSup(aTransaction[i], tPattern) == 1){
					tmp = 1;
				}else{
					tmp = aPrefix.calcInexactSup(aTransaction[i], tPattern);
				}
				count++;
				//cout << "map" << endl;
				mSupportMap.emplace(tKey, tmp);
			}else{
				//cout << "Discovery!!" << endl;
				tmp = mSupportMap[tKey];
			}
			//cout << tmp << endl;
			tYHat[i] += tW * tmp;
		}
	}
	*/



/*
	for(auto tActiveNode : aPrefix.mActive){
		vector<Event> tPattern = tActiveNode->pattern;
		double tW = tActiveNode->w;
		for(uint i = 0; i < tTransactionSize; ++i){
			tYHat[i] += tW * aPrefix.calcInexactSup(aTransaction[i], tPattern);
		}
	}
	*/


	//for(int i = 0; i <tYHat.size(); i++){
		//cout << "YHat[" << i << "] = " << tYHat[i] << endl;
	//}

	//cout << "mapInsert count = " << count << endl;
	return tYHat;
}

vector<vector<double> > SVMSPP::get_all_predict(PrefixSpan &aPrefix, const vector<double> &aLambdas, const vector<vector<Event> > &aTransaction, const vector<uint> &aOptions, const vector<uint> &testId,PrefixSpan &tPrefix){

	mT = aLambdas.size();
	//とりあえずsolution path すべてで学習を行う
	vector<vector<double> > tYHats;

	double tTrueLambdaMax = get_lambda_max(aPrefix);

	cout << "tTrueLambdaMax" << tTrueLambdaMax << endl;
	get_lambda_max(tPrefix);

	if(aOptions.empty()){
		cout << "error:svm.learn option size is incorrect." << '\n';
		exit(1);
	}
	//λを求めるための処理
	if(aLambdas.empty()){
		cout << "error:lambda is empty. Set lambda." << '\n';
		exit(1);
	}

	// compute solution path
	uint activesize = 1;
	double L1norm = 0;
	double maxval = tTrueLambdaMax; //lambda_max

	for(uint t = 1; (t < aOptions[0]) && (t < mT); ++t){
		clock_t start = clock();
		double lam = aLambdas[t];
		cout << "aLambdas[" << t <<"] = " << lam << endl;

		if(lam > tTrueLambdaMax){
			//cout << "aa" << endl;
			vector<double> tV(aTransaction.size(), mBias);
			tYHats.push_back(tV);
			continue;
		}
		//cout << "bb" << endl;
		vector<uint> index;
		vector<uint *> x;
		vector<double *> support;
		// vector<uint> norm;
		vector<double> norm;
		vector<uint> xSize;
		vector<double> grad;
		vector<double> w;
		for(uint iter = 0; iter <= mMaxIter; ++iter){
			// calculate dual and safe screening
			if(iter % mFreq == 0){
				double loss = 0;
				double oTr = 0;
				for(uint i = 0; i < mN; ++i){
					if(mR[i] > 0){
						loss += mR[i] * mR[i];
						oTr += mR[i];
					}
				}
				if(iter >= 1){
					maxval = 0;
					for(uint s = 0; s < activesize; ++s){

						uint j = index[s];
						grad[j] = 0;

						// for (uint k = 0; k < norm[j]; ++k) {
						for(uint k = 0; k < xSize[j]; ++k){
							uint id = x[j][k];
							if(mR[id] > 0){
								grad[j] += aPrefix.mY[id] * mR[id] * support[j][k];
								//@
							}
						}
						if(fabs(grad[j]) > maxval) maxval = fabs(grad[j]);
					}
				}

				double alpha = min(max(oTr / (lam * loss), 0.0), 1 / maxval);
				double primal = 0.5 * loss + lam * L1norm;
				double dual = -0.5 * lam * lam * alpha * alpha * loss + lam * alpha * oTr;
				double gap = primal - dual;
				if(gap / primal < mEps){
					//最初から収束条件を満たしている時
					if(index.empty()){
						cout << "gap = " << gap << endl;
						cout << "primal = " << primal << endl;
						cout << "empty break!" << endl;
						break;
					}
					uint active = 0;
					for(uint s = 0; s < activesize; s++){
						uint j = index[s];
						aPrefix.mActive[j]->w = w[j];
						if(w[j] != 0) active++;
					}

					//cout << "[iter " << iter << "] primal: " << primal << ", dual: " << dual << ", gap: " << gap/primal << ", activeset: " << active << endl;
					break;
				}
				double radius = sqrt(2 * gap) / lam;
				if(iter == 0){

					//各λループのはじめに一回だけscreeningを行う(現在のラムダでの最大の木を構築している)
					//λは徐々に小さくなって行くの木は前回のλより恐らく大きくなる
					//w更新するたびバウンドがきつくなり，ある周期で非アクティブを決め，木が小さくなって行く
					aPrefix.safe_screening(mR, alpha, radius, 1);
					activesize = aPrefix.mActive.size();

					index.resize(activesize);
					x.resize(activesize);
					norm.resize(activesize);
					grad.resize(activesize);
					w.resize(activesize);
					xSize.resize(activesize);
					support.resize(activesize);

					for(uint j = 0; j < activesize; ++j){
						index[j] = j;
						x[j] = aPrefix.mActive[j]->x.data();
						support[j] = aPrefix.mActive[j]->support.data();
						// norm[j] = aPrefix.Active[j]->x.size();
						xSize[j] = aPrefix.mActive[j]->x.size();
						norm[j] = aPrefix.mActive[j]->supportSumsqr;
						grad[j] = 0;
						w[j] = aPrefix.mActive[j]->w;
					}

				}else{
					//最初の１回目以外は非アクティブノード検出を行う．（木の大きさは変わらない）
					for(uint s = 0; s < activesize; ++s){
						uint j = index[s];
						double score = fabs(alpha * grad[j]) + radius * sqrt(norm[j]);

						if(score < 1){
							if(w[j] != 0){
								// for (uint k = 0; k < norm[j]; ++k) {
								for(uint k = 0; k < xSize[j]; ++k){
									uint id = x[j][k];
									mR[id] += aPrefix.mY[id] * w[j] * support[j][k];          //@
								}
							}
							//学習中にwが0になってもActiveのサイズは変更していない
							aPrefix.mActive[j]->w = 0;
							activesize--;
							swap(index[s], index[activesize]);
							s--;
						}
					}
				}
				//cout << "[iter " << iter << "] primal: " << primal << ", dual: " << dual << ", gap: " << gap/primal << ", activeset: " << activesize << endl;
			}

			//ランダムに混ぜている
			for(uint j = 0; j < activesize; ++j){
				uint i = j + rand() % (activesize - j);
				swap(index[i], index[j]);
			}

			// update w
			//newton法
			L1norm = 0;
			for(uint s = 0; s < activesize; ++s){
				uint j = index[s];

				double G = 0;
				double H = 0;
				double loss_old = lam * fabs(w[j]);      //@
				// for (uint k = 0; k < norm[j]; ++k) {
				for(uint k = 0; k < xSize[j]; ++k){
					uint idx = x[j][k];
					if(mR[idx] > 0){
						G -= aPrefix.mY[idx] * mR[idx] * support[j][k];
						H += support[j][k] * support[j][k];
						loss_old += 0.5 * mR[idx] * mR[idx];
					}
				}
				H = max(H, 1e-12);

				double Gp = G + lam;
				double Gn = G - lam;
				double d;
				if(Gp <= H * w[j]) d = -Gp / H;
				else if(Gn >= H * w[j]) d = -Gn / H;
				else d = -w[j];

				if(fabs(d) < 1e-12){
					L1norm += fabs(w[j]);
					continue;
				}

				double d_old = 0;

				double bound = 0.5 * (G * d + lam * fabs(w[j] + d) - lam * fabs(w[j]));
				uint linesearch;
				for(linesearch = 0; linesearch < 10; ++linesearch){
					double loss_new = lam * fabs(w[j] + d);
					// for (uint k = 0; k < norm[j]; ++k) {
					for(uint k = 0; k < xSize[j]; ++k){
						uint idx = x[j][k];
						//損失の入れ替え
						mR[idx] += aPrefix.mY[idx] * (d_old - d) * support[j][k];
						if(mR[idx] > 0){
							loss_new += 0.5 * mR[idx] * mR[idx];
						}
					}

					if(loss_new - loss_old <= bound){
						break;
					}else{
						d_old = d;
						d *= 0.5;
						bound *= 0.5;
					}
				}
				w[j] += d;
				L1norm += fabs(w[j]);
			}

			// bias
			uint nn = 0;
			double tmp = 0;
			for(uint i = 0; i < mN; i++){
				if(mR[i] > 0){
					tmp += aPrefix.mY[i] * mR[i] + mBias;
					nn++;
				}
			}
			double bias_old = mBias;
			mBias = tmp / nn;
			for(uint i = 0; i < mN; i++){

				mR[i] += aPrefix.mY[i] * (bias_old - mBias);
			}
		}
		//cout << "cc" << endl;

		clock_t end = clock();
		const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000.0;
		//printf("最適化 time %lf[ms]\n", time);

		clock_t start1 = clock();

		//cout << "dd" << endl;

		tYHats.push_back(predict(aPrefix, aTransaction, testId,tPrefix));
		clock_t end1 = clock();
		const double time1 = static_cast<double>(end1 - start1) / CLOCKS_PER_SEC * 1000.0;
		//printf("predict  time %lf[ms]\n", time1);
		//cout << "ee" << endl;
	}

	mSupportMap.clear();
	return tYHats;
}

double SVMSPP::get_bias(void){
	return mBias;
}

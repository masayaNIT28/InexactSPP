#include "prefixspan.h"



/**
 * @fn
 * 数字として管理しているパターンを文字列にする
 * @return あるノードにおけるパターンを文字列で出力
 *mEventSizeとは１；１；１だったら３ということ
 */
string PrefixSpan::pat2str(void){
	stringstream ss;
	for(uint i = 0; i < mPattern.size(); ++i){
		ss << (i ? " " : "");
		if(mPattern[i].first.size() > 0){
			ss << "(";
			for(uint j = 0; j < mPattern[i].first.size(); ++j){
				ss << mPattern[i].first[j];
				ss << (j + 1 < mPattern[i].first.size() ? "_" : "");
			}
			ss << ")";
			if(mPattern[i].second.size() > 0) ss << ":";
		}
		for(uint j = 0; j < mItemSize; j++){
			if(mPattern[i].second[j] == MAXVAL) ss << "*";
			else ss << mPattern[i].second[j];
			ss << (j + 1 < mItemSize ? ":" : "");
		}
	}
	return ss.str();
}

string PrefixSpan::pat2str(const vector<Event> aPattern){
	stringstream ss;
	for(uint i = 0; i < aPattern.size(); ++i){
		ss << (i ? " " : "");
		if(aPattern[i].first.size() > 0){
			ss << "(";
			for(auto it : aPattern[i].first){
				ss << it;
				if(it != aPattern[i].first.back()) ss << "_";
			}
			ss << ")";
			if(aPattern[i].second.size() > 0) ss << ":";
		}
		for(auto it : aPattern[i].second){
			if(it == MAXVAL) ss << "*";
			else ss << it;
			if(it != aPattern[i].second.back()) ss << ":";
		}
	}
	return ss.str();
}

void PrefixSpan::printTree(string aFilename, double bias){
	ofstream tFile(aFilename);
	uint tCnt = 0;
	for(auto it : mActive){
		if(it->w != 0){
			tCnt++;
		}
	}

	tFile << "Size of data :" << mN;
	tFile << ",Bias," << bias << '\n';
	tFile << "Size of mTree :" << mTree.size() << '\n';
	tFile << "Size of Active :" << tCnt << '\n';

	tFile << "mPattern,supportSum,addλ,w,list[id:sup]" << '\n';
	//木丸ごと書き込む時はこっち
	// for(auto itr = mTree.begin();itr != mTree.end();++itr){
	// 		 tFile << itr->patternStr << "," <<itr->supportSum << "," << itr->w << "," << itr->isLeaf << '\n';
	// }
	// アクティブだけ
	for(auto it : mActive){
		if(it->w != 0){
			tFile << it->patternStr << "," << it->supportSum << "," << it->addLambda << "," << it->w << ",[";
			for(uint j = 0; j < it->support.size(); j++)
				tFile << it->x[j] << ":" << it->support[j] << ",";
			tFile << "]" << '\n';
		}

	}
}

void PrefixSpan::printTree(string aFilename){
	ofstream tFile(aFilename);
	uint tCnt = 0;
	for(auto it : mActive){
		if(it->w != 0){
			tCnt++;
		}
	}

	tFile << "Size of data :" << mN << '\n';
	tFile << "Size of mTree :" << mTree.size() << '\n';
	tFile << "Size of Active :" << tCnt << '\n';
	//tFile << "Bias," << bias << '\n';
	tFile << "mPattern,supportSum,addλ,w,list[id:sup]" << '\n';
	//木丸ごと書き込む時はこっち
	// for(auto itr = mTree.begin();itr != mTree.end();++itr){
	// 		 tFile << itr->patternStr << "," <<itr->supportSum << "," << itr->w << "," << itr->isLeaf << '\n';
	// }
	// アクティブだけ
	for(auto it : mActive){
		if(it->w != 0){
			tFile << it->patternStr << "," << it->supportSum << "," << it->addLambda << "," << it->w << ",[";
			for(uint j = 0; j < it->support.size(); j++)
				tFile << it->x[j] << ":" << it->support[j] << ",";
			tFile << "]" << '\n';
		}

	}
}


void PrefixSpan::add_history(uint aLambda){
	for(auto it : mActive){
		if(it->w != 0){
			if(it->addLambda == -1){
				it->addLambda = aLambda;
			}
		}else{
			it->addLambda = -1;
		}
	}
}


/**
* あるデータの中の任意のパターンのサポートを数える
* mode:0,(1,2,3以外の数字が入った時) supportは１レコード1まで
* mode:1 単純なパターンの数え上げ，1レコードにいくらでも可能
* @return supprt数
**/
int PrefixSpan::calcSup(const vector<Event> &aTransaction, const vector<Event> &aPattern) const{
	uint patSize = aPattern.size();
	uint k;
	int num = 0;
	int interval;

	//想定しているパターンよりデータが短い場合は0を返す
	if(aTransaction.size() < patSize) return num;

	for(uint i = 0; i <= aTransaction.size() - patSize; i++){

		interval = 0;

		k = 0;

		for(uint j = i; j < aTransaction.size(); j++){

			if(mInterval >= 0 && interval > mInterval){
				break;
			}
			if(aTransaction[j] == aPattern[k]){
				k++;
				interval = -1;
				if(k == patSize){
					//0/1の場合
					if(mSupMode == 0){
						return 1;
					}
					k = 0;
					num++;
					i = j;
					break;
				}
			}
			interval++;
		}
	}
	return num;
}

/**
 * transactionの中のpattern数を数える
 */
double PrefixSpan::calcSup(uint aId, vector<Event> aPattern){
	switch(mSupMode){
		case 1:
			return calcPat(aId, aPattern);
			break;
		default:
			return 1;
			break;
	}
}

/**
* trainTransaction[aId]であるデータ内からパターンがaPatternであるものを数え上げる
* @return 数え上げた数
*/
int PrefixSpan::calcPat(uint aId, vector<Event> aPattern){
	uint patSize = aPattern.size();
	uint k;
	int num = 0;
	int interval;

	//想定しているパターンよりデータが短い場合は0を返す
	if(mTransaction[aId].size() < patSize) return num;

	for(uint i = 0; i <= mTransaction[aId].size() - patSize; i++){

		interval = 0;

		k = 0;
		for(uint j = i; j < mTransaction[aId].size(); j++){

			if(mInterval >= 0 && interval > mInterval){
				break;
			}
			if(mTransaction[aId][j] == aPattern[k]){
				k++;
				interval = -1;
				if(k == patSize){
					k = 0;
					num++;
					i = j;
					break;
				}
			}
			interval++;
		}
	}
	return num;
}

//このポジションでの類似度

double PrefixSpan::calcInexactSup(projectDB tPdb, projectDB &minpDB, vector<Event> aPattern){
	//cout << "calcInexactSup" << endl;
	//cout << pat2str(aPattern) << endl;
	uint patSize = aPattern.size();
	Position tminpos;
	projectDB tminpDB;
	double max = mMaxGap*mMaxpat;

	double min = mMinitemdifferent;

	uint fmode = 1;
	uint flag = 0;
	//cout << "max = " << max << endl;
	//cout << "min = " << min << endl;
	//同じidのものだけが引数として入力される
	uint id = tPdb[0].sequence_index;
	//cout << id << endl;
	//cout << pat2str(mTransaction[id]) << endl;
	if(mTransaction[id].size() < patSize){
		return 0.0;
	}

	for(uint i = 0, size = tPdb.size(); i < size ;i++){
		//cout << "bb" << endl;
		Position pos = tPdb[i];

		uint event = pos.event_index;
		uint item = pos.itemset_index;

		if(event == 0){
			for(uint j = 0; j < mTransaction[id].size() && j+patSize-1 < mTransaction[id].size(); j++){
				uint tmp = 0;
				for(uint k = 0; k < patSize && j+k < mTransaction[id].size(); k++){
					//cout << "k " << k << endl;

					Event a = mTransaction[id][j+k];
					//vector<Event> aa;
					//aa.push_back(a);
					//cout << pat2str(aa) << endl;
					//cout << "asize  " << a.second.size() << endl;
					Event b = aPattern[k];
					for(uint l = 0; l < b.second.size(); l++){
						//cout << "l " << l << endl;
						int c;
						if(fmode == 1){
							//cout << b.second[l] << endl;

							if(a.second[l]==b.second[l]){
								//cout << "cc" << endl;
								c=0;
							}else{
								//cout << "dd" << endl;
								c=1;
							}
						}else{
							c=a.second[l]-b.second[l];
						}
						//cout << "cc" << endl;
						tmp+= abs(c);
						if(tmp > mMinitemdifferent){
							l = b.second.size();
							k = patSize;
							break;
						}
					}
				}
				if(min == tmp){
					tminpos.sequence_index = id;
					tminpos.event_index = i+patSize-1;
					tminpos.itemset_index = 0;
					tminpDB.push_back(tminpos);
				}else if(min > tmp){
					min = tmp;
					tminpDB = {};
					tminpos.sequence_index = id;
					tminpos.event_index = i+patSize-1;
					tminpos.itemset_index = 0;
					tminpDB.push_back(tminpos);
					//cout << "min" << min << endl;
				}
			}
		}else{
			if(event+1 < mTransaction[id].size() && mTransaction[id][event+1] == aPattern[patSize-1]){
				tminpos.sequence_index = id;
				tminpos.event_index = i+patSize-1;
				tminpos.itemset_index = 0;
				tminpDB.push_back(tminpos);
				flag = 1;
			}else{
				uint tmp=0;
				for(int j = patSize-1;j >= 0 && mTransaction[id].size() > event+j-patSize+2; j--){
					//cout << "j " << j << endl;
					Event a = mTransaction[id][event+j-patSize+2];
					Event b = aPattern[j];
					for(uint k = 0; k < b.second.size(); k++){
						//cout << "k " << k << endl;
						int c;
						if(fmode == 1){
							if(a.second[k]==b.second[k]){
								c=0;
							}else{
								c=1;
							}
						}else{
							c=a.second[k]-b.second[k];
						}
						tmp+= abs(c);
						if(tmp > mMinitemdifferent){
							k = b.second.size();
							j = 0;
							break;
						}
					}
				}
				if(min == tmp){
					tminpos.sequence_index = id;
					tminpos.event_index = i+patSize-1;
					tminpos.itemset_index = 0;
					tminpDB.push_back(tminpos);
				}else if(min > tmp){
					min = tmp;
					tminpDB = {};
					tminpos.sequence_index = id;
					tminpos.event_index = i+patSize-1;
					tminpos.itemset_index = 0;
					tminpDB.push_back(tminpos);
					//cout << "min" << min << endl;
				}
			}
		}
	}

	//cout << "flag" << endl;
	if(flag == 0){
		for(uint i = 0; i < mTransaction[id].size(); i++){
			uint tmp = 0;
			//cout << "i " << i << endl;
			if(mTransaction[id][i] == aPattern[patSize-1] && i >= patSize-1){
				for(int j = patSize-1; j >= 0; j--){
					//cout << "j " << j << endl;
					Event a = mTransaction[id][i+j-patSize+1];
					Event b = aPattern[j];
					for(uint k = 0; k < b.second.size(); k++){
						//cout << "k " << k << endl;
						int c;
						if(fmode == 1){
							if(a.second[k]==b.second[k]){
								c=0;
							}else{
								c=1;
							}
						}else{
							c=a.second[k]-b.second[k];
						}
						tmp+= abs(c);
						if(tmp > mMinitemdifferent){
							k = b.second.size();
							j = 0;
							break;
						}
					}
				}
				if(min == tmp){
					tminpos.sequence_index = id;
					tminpos.event_index = i+patSize-1;
					tminpos.itemset_index = 0;
					tminpDB.push_back(tminpos);
				}else if(min > tmp){
					min = tmp;
					tminpDB = {};
					tminpos.sequence_index = id;
					tminpos.event_index = i+patSize-1;
					tminpos.itemset_index = 0;
					tminpDB.push_back(tminpos);
					//cout << "min" << min << endl;
				}
			}
		}
	}
	//cout << "cc" << endl;
	if(tminpDB.size() == 0){
		//cout << "cc" << endl;
		return 0.0;
	}else{
		for(uint i = 0; i < tminpDB.size(); i++){
			Position p = tminpDB[i];
			minpDB.push_back(p);
		}
			//cout << "cc" << endl;

		return 1.0 - min/max;
	}
	return 0.0;
}

double PrefixSpan::calcInexactSup(const vector<Event> &aTransaction, const vector<Event> aPattern) const{
	//cout << "calcInexactSup" << endl;
	//cout << pat2str(aPattern) << endl;
	uint patSize = aPattern.size();
	double max = mMaxGap*mMaxpat;

	double min = max;

	uint fmode = 1;
	uint flag = 0;
	//cout << "max = " << max << endl;
	//cout << "min = " << min << endl;
	//同じidのものだけが引数として入力される
	uint size = aTransaction.size();
	//cout << id << endl;
	//cout << pat2str(mTransaction[id]) << endl;
	if(size < patSize){
		return 0.0;
	}

	for(uint i = 0; i < size ;i++){
		//cout << "bb" << endl;
			for(uint j = 0; j < size && j+patSize-1 < size; j++){
				uint tmp = 0;
				for(uint k = 0; k < patSize && j+k < size; k++){
					//cout << "k " << k << endl;

					Event a = aTransaction[j+k];
					//vector<Event> aa;
					//aa.push_back(a);
					//cout << pat2str(aa) << endl;
					//cout << "asize  " << a.second.size() << endl;
					Event b = aPattern[k];
					for(uint l = 0; l < b.second.size(); l++){
						//cout << "l " << l << endl;
						int c;
						if(fmode == 1){
							//cout << b.second[l] << endl;

							if(a.second[l]==b.second[l]){
								//cout << "cc" << endl;
								c=0;
							}else{
								//cout << "dd" << endl;
								c=1;
							}
						}else{
							c=a.second[l]-b.second[l];
						}
						//cout << "cc" << endl;
						tmp+= abs(c);
						if(tmp > mMinitemdifferent){
							l = b.second.size();
							k = patSize;
							break;
						}
					}
				}
				if(min > tmp){
					min = tmp;
					//cout << "min" << min << endl;
				}
			}
		}

	if(min <= mMinitemdifferent){
		return 1.0 - min/max;
	}
	return 0.0;
}


//初期の木を構築している
void PrefixSpan::init(const vector<vector<Event> > aTransaction, const vector<double> aY){

	mTransaction = aTransaction;
	mY = aY;
	mItemSize = mTransaction[0][0].second.size();
	mN = mTransaction.size();
	if(mTransaction[0][0].first.size() > 0) mFlagItemsetExist = 1;

	Itemset tItemset;
	vector<uint> tWild(mItemSize, MAXVAL);
	Event tEvent;
	tEvent.first = tItemset;
	tEvent.second = tWild;
	mWildEvent = tEvent;

	//nはデータ数
	if(mTransaction.empty() || mY.empty()){
		cout << "error:Data or label is empty." << '\n';
		exit(1);
	}
	mPattern.clear();
	mTree.clear();
	mActive.clear();
	mR.clear();
	mCheckClosed.clear();
	mCheckPDB.clear();


	//木の根をここで作成している
	//パターン、通常のsupportDB、InexactのsupportDB
	map<Event, projectDB> tCounter;

	map<Event, uint> dupCheck;
	//イベント間の間隔に制限を設けない時は初期ノードの重複を許さない
	if(mInterval < 0){
		//ルートにはロードしたデータのイベントが全て入っている
		for(uint i = 0; i < mN; ++i){
			//イベントごとにマッピングされる
			for(uint j = 0, size = mTransaction[i].size(); j < size; ++j){
				Event tEvent;
				vector<uint> tItemset;
				tEvent.second = mTransaction[i][j].second;
				uint k = 0;
				if(mTransaction[i][j].first.size() > 0){
					for(; k < mTransaction[i][j].first.size(); ++k){
						tItemset.push_back(mTransaction[i][j].first[k]);
						tEvent.first = tItemset;
						if(dupCheck.find(tEvent) == dupCheck.end()){
							dupCheck[tEvent] = 0;
							//i jはi番目の系列のj番目のイベント，アイテムセットがあるならk番目のアイテム
							Position tPosition;
							tPosition.sequence_index = i;
							tPosition.event_index = j;
							tPosition.itemset_index = k;
							tCounter[tEvent].push_back(tPosition);
						}
						tItemset.clear();
					}
				}else{
					tEvent.first = tItemset;
					if(dupCheck.find(tEvent) == dupCheck.end()){
						dupCheck[tEvent] = 0;
						Position tPosition;
						tPosition.sequence_index = i;
						tPosition.event_index = j;
						tPosition.itemset_index = k;
						tCounter[tEvent].push_back(tPosition);
					}
				}
			}
			dupCheck.clear();
		}
	}else{
		for(uint i = 0; i < mN; ++i){
			for(uint j = 0, size = mTransaction[i].size(); j < size; ++j){
				Event tEvent;
				vector<uint> tItemset;
				tEvent.second = mTransaction[i][j].second;
				uint k = 0;
				if(mTransaction[i][j].first.size() > 0){
					for(; k < mTransaction[i][j].first.size(); ++k){
						tItemset.push_back(mTransaction[i][j].first[k]);
						tEvent.first = tItemset;
						//i jはi番目の系列のj番目のイベント，アイテムセットがあるならk番目のアイテム
						Position tPosition;
						tPosition.sequence_index = i;
						tPosition.event_index = j;
						tPosition.itemset_index = k;
						tCounter[tEvent].push_back(tPosition);
						tItemset.clear();
					}
				}else{
					tEvent.first = tItemset;
					Position tPosition;
					tPosition.sequence_index = i;
					tPosition.event_index = j;
					tPosition.itemset_index = k;
					tCounter[tEvent].push_back(tPosition);
				}
			}
		}
	}
	mPattern.clear();

	//root用の空ノード
	Node tRoot;
	mTree.push_back(tRoot);

	//tCounter.begin()でイテレータの最初を取り出して，.end()で終わりを取り出せる
	for(auto it = tCounter.begin(), end = tCounter.end(); it != end; ++it){
		//キーつまりイベントが格納される
		mPattern.push_back(it->first);
		Node tNode;
		tNode.parent = mTree.begin();
		tNode.pattern = mPattern;
		tNode.patternStr = pat2str();
		tNode.supportSum = 0;
		tNode.supportSumsqr = 0;
		tNode.pdb = it->second;
		tNode.w = 0;
		tNode.val = 0;
		tNode.closed = true;
		//ProjectDB tIpdb;


		uint oid = MAXVAL;
		vector<uint> idList;
		for(uint i = 0; i < mTransaction.size(); i++){
			idList.push_back(i);
		}
		//projectDBのサイズ回ループ
		for(uint i = 0, size = tNode.pdb.size(); i < size; ++i){
			//idはレコード番号
			uint id = tNode.pdb[i].sequence_index;
			if(id != oid){
				idList.erase(remove(idList.begin(), idList.end(), id), idList.end());
				// tNode.supportSum++;
				double tSup = calcSup(id, mPattern);
				tNode.supportSum += tSup;
				tNode.supportSumsqr+= tSup*tSup;
				tNode.support.push_back(tSup);
				tNode.x.push_back(id);
			}
			oid = id;
		}
		//support1以外かつ閾値以上のidについて計算
		for(uint i = 0; i < idList.size(); i++){
			uint id = idList[i];
			//for(uint j = 0; j < mTransaction[id].size() ;j++){
				//for(uint k = 0; k < mTransaction[id][j].first.size(); k++){
					//Position tPosition;
					//tPosition.sequence_index = id;
					//tPosition.event_index = j;
					//tPosition.itemset_index = k;
					//tIpd.push_back(tPosition);
				//}
			//}
			//double tSup = calcInexactSup(id, mPattern);

			double tSup = 1.0- 1.0/(mMaxGap*mMaxpat);
			if(tSup >= mMinInexactsup){
				Position pos;
				pos.sequence_index = id;
				pos.event_index = 0;
				pos.itemset_index = 0;
				//tNode.Inexact_idList.push_back(id);
				tNode.inexact_pdb.push_back(pos);
				tNode.supportSum += tSup;
				tNode.supportSumsqr += tSup*tSup;
				tNode.support.push_back(tSup);
				tNode.x.push_back(id);
			}
		}

		bool tAddFlag = true;
		if(mMinsup >= 1 && mMinsup > tNode.supportSum) tAddFlag = false;
		else if(mMinsup < 1 && mMinsup > (double) tNode.supportSum / mN) tAddFlag = false;

		if(tAddFlag){
			tNode.pdb_size = sum_items_pdb(tNode.pdb);
			//tNode.inexact_pdb(tPosition);
			mTree.push_back(tNode);
			mTree.begin()->child.push_back(prev(mTree.end()));
		}
		mPattern.pop_back();
	}

}

/**
 * 入力された系列が包含関係かどうか調べる．
 * 包含関係なら1番目か2番目のどちらが短いか（subsequence）を返す．
 * @return 1 or 2 or 0(包含関係ではない or 全く同一系列)
 * @author takuto
 */
uint PrefixSpan::isSubsequence(const vector<Event> aSeq1, const vector<Event> aSeq2){
	vector<Event> tShort_seq;
	vector<Event> tLong_seq;
	uint tSub_seq;

	if(aSeq1.size() < aSeq2.size()){
		tShort_seq = aSeq1;
		tLong_seq = aSeq2;
		tSub_seq = 1;
	}else if(aSeq1.size() > aSeq2.size()){
		tShort_seq = aSeq2;
		tLong_seq = aSeq1;
		tSub_seq = 2;
	}else if(mFlagItemsetExist){
		uint tSum1 = 0;
		uint tSum2 = 0;
		for(auto it : aSeq1){
			tSum1 += it.first.size();
		}
		for(auto it : aSeq2){
			tSum2 += it.first.size();
		}

		if(tSum1 > tSum2){
			tShort_seq = aSeq2;
			tLong_seq = aSeq1;
			tSub_seq = 2;
		}else if(tSum1 < tSum2){
			tShort_seq = aSeq1;
			tLong_seq = aSeq2;
			tSub_seq = 1;
		}else{
			return 0;
		}

	}else{
		return 0;
	}

	uint tCount = 0;
	uint diff_size = tLong_seq.size() - tShort_seq.size();

	for(uint i = 0; i <= diff_size; ++i){
		for(uint it = 0; it < tShort_seq.size(); ++it){
			if(isInclude(tLong_seq[it + i], tShort_seq[it]) || tShort_seq[it] == mWildEvent || tLong_seq[it + i] == mWildEvent){
				tCount++;
			}else{
				tCount = 0;
				break;
			}

			if(tCount == tShort_seq.size()){
				return tSub_seq;
			}
		}
	}

	return 0;
}

/**
 * Event1がEvent2を包含しているか調べる
 * @return true (Event1は2を包含している) or false
 */
bool PrefixSpan::isInclude(const Event aEvent1, const Event aEvent2){
	if(aEvent1.first.size() < aEvent2.first.size()){
		return false;
	}else if(aEvent1.second != aEvent2.second){
		return false;
	}else{
		uint i = 0;
		bool tExistFlag = false;
		for(auto itr : aEvent2.first){
			tExistFlag = false;
			if(i == aEvent1.first.size()){
				return false;
			}
			while(itr >= aEvent1.first[i]){
				if(itr == aEvent1.first[i]){
					tExistFlag = true;
					i++;
					break;
				}else{
					i++;
				}

				if(i == aEvent1.first.size()){
					return false;
				}
			}

			if(!tExistFlag){
				return false;
			}
		}
		return true;
	}

}


/**
 * aNodeIterがaChildIterの親かどうか辿って調べる
 */
bool PrefixSpan::isParent(const TreeIter aNodeIter, const TreeIter aChildIter){
	auto current = aChildIter->parent;
	while(current->pattern.size() >= aNodeIter->pattern.size()){
		if(current == aNodeIter) return true;

		current = current->parent;
	}
	return false;
}

/**
 * CloSpan用
 * @return PDB内の総アイテム数
 * @author takuto
 */
uint PrefixSpan::sum_items_pdb(const projectDB aPdb){
	uint tSum = 0;
	for(auto itr : aPdb){
		uint id = itr.sequence_index;
		uint j = itr.event_index;
		if(mTransaction[id][j].first.size() > 0){
			tSum += mTransaction[id][j].first.size() - 1 - itr.itemset_index;
		}else{
			tSum++;
		}
		for(++j; j < mTransaction[id].size(); ++j){
			if(mTransaction[id][j].first.size() > 0){
				tSum += mTransaction[id][j].first.size();
			}else{
				tSum++;
			}
		}
	}
	return tSum;

}

/**
 * CloSpan用
 * @return あるパターンを含むIDの和(IDの重複しない)
 */
uint PrefixSpan::sumId(const vector<uint> aIdList){
	uint tSum = 0;
	// setを用いてidの重複を消す
	set<uint> tIdSet;

	for(auto tId: aIdList){
		tIdSet.insert(tId);
	}

	for(auto itr = tIdSet.begin();itr != tIdSet.end();++itr){
    	tSum += *itr;
	}

	return tSum;
}


/**
 * @fn
 * SPPC計算を行っている(type == 2の時)
 * @return 枝が切ることができればtrue
 */
bool PrefixSpan::calculate(Node& aNode){
	//cout << "calculate" << endl;
	if(aNode.supportSum < mMinsup) return true;

	aNode.val = 0;
	double p = 0, m = 0;
	//見ているノードmPatternを持っているid数で回す
	for(uint i = 0; i < aNode.x.size(); ++i){
		uint id = aNode.x[i];
		if(mType == 1 || mType == 2){ // svm
			if(mR[id] > 0){
				//(初めは)vが1-ybとなる，これはλmaxのときの損失の値
				//@
				double val = mAlpha * mR[id] * mY[id] * aNode.support[i];
				aNode.val += val;
				(val > 0) ? p += val : m += val;
			}
		}else if(mType == 3 || mType == 4){ //lasso
			double val = mAlpha * mR[id] * aNode.support[i];
			aNode.val += val;
			(val > 0) ? p += val : m += val;
		}else if(mType == 5 || mType == 6){ //logistic
			double val = (mY[id] > 0) ? 1 / (1 + mR[id]) : 1 / (1 + mR[id]) - 1;
			val *= mAlpha * aNode.support[i];
			aNode.val += val;
			(val > 0) ? p += val : m += val;
		}

	}
	aNode.val = fabs(aNode.val);

	if(mType == 1 || mType == 3 || mType == 5){ // get_mMaxnode
		if(max(p, -m) < mMaxnode.val) return true;
	}else if(mType == 2 || mType == 4 || mType == 6){ // safe_screening
		if(max(p, -m) + mRadius * sqrt(aNode.supportSumsqr) < 1) return true;
	}

	return false;
}

/**
 * @fn
 * project_new内でのみ使用
 * サポートを数えている
 * @return 葉になることができればtrue
 */
bool PrefixSpan::calculate_new(Node &aNode){
	//cout << "calculate_new" << endl;
	//cout << aNode.patternStr << endl;
	uint oid = MAXVAL;
	double p = 0, m = 0;

	vector<uint> idList;
	projectDB oProDB = aNode.parent->pdb;
	projectDB tInexact_pdb = aNode.parent->inexact_pdb;
	for(uint i = 0, size = oProDB.size(); i < size; i++){
		uint id = oProDB[i].sequence_index;
		if(oid != id){
			idList.push_back(id);
		}
		oid = id;
	}
	oid = MAXVAL;
	//見ているノードの持っているデータベースの長さ回まわす
	//pdb.sizeは同じデータ内であっても複数存在する可能性あり
	//init()でやっていることを行っているので新しいノードを展開した時に多分使用する
	for(uint i = 0, size = aNode.pdb.size(); i < size; ++i){
		uint id = aNode.pdb[i].sequence_index;
		if(oid != id){
			idList.erase(remove(idList.begin(), idList.end(), id), idList.end());
			// node.supportSum++;
			double tSup = calcSup(id, aNode.pattern);
			aNode.supportSum += tSup;
			aNode.supportSumsqr += tSup*tSup;
			aNode.support.push_back(tSup);
			aNode.x.push_back(id);
			if(mType == 1 || mType == 2){ // svm
				if(mR[id] > 0){
					double val = mAlpha * mR[id] * mY[id] * tSup;
					aNode.val += val;
					(val > 0) ? p += val : m += val;
				}
			}else if(mType == 3 || mType == 4){ // lasso
				double val = mAlpha * mR[id] * tSup;
				aNode.val += val;
				(val > 0) ? p += val : m += val;
			}else if(mType == 5 || mType == 6){ // logistic
				double val = (mY[id] > 0) ? 1 / (1 + mR[id]) : 1 / (1 + mR[id]) - 1;
				val *= mAlpha * tSup;
				aNode.val += val;
				(val > 0) ? p += val : m += val;
			}
		}
		oid = id;
	}

	//support1出なくなった系列をinexact_pdbに入れる
	if(mMinitemdifferent > 0){
		for(uint i = 0, size = idList.size(); i < size; ++i){
			Position tPosition;
			tPosition.sequence_index = idList[i];
			tPosition.event_index = 0;
			tPosition.itemset_index = 0;
			tInexact_pdb.push_back(tPosition);
		}
	}

	oid = MAXVAL;

	//support1以外かつ閾値以上のidについて計算

	//同じidをまとめて入れてcalcする
	projectDB tpdb;
	projectDB ans;
	//cout << "aa" << endl;
	for(uint i = 0, size = tInexact_pdb.size(); i < size; ++i){
		uint id = tInexact_pdb[i].sequence_index;
		uint next;
		if(i != size-1){
			next = tInexact_pdb[i+1].sequence_index;
		}else{
			next = id;
		}
		tpdb.push_back(tInexact_pdb[i]);
		if(id != next){

		// node.supportSum++;
			double tSup = calcInexactSup(tpdb, ans, aNode.pattern);
			tpdb ={};
			//cout << tSup << endl;
			if(tSup > 0.0){
				//cout << tSup << endl;
				aNode.supportSum += tSup;
				aNode.supportSumsqr += tSup*tSup;
				aNode.support.push_back(tSup);
				aNode.x.push_back(id);
				if(mType == 1 || mType == 2){ // svm
					if(mR[id] > 0){
						double val = mAlpha * mR[id] * mY[id] * tSup;
						aNode.val += val;
						(val > 0) ? p += val : m += val;
					}
				}else if(mType == 3 || mType == 4){ // lasso
					double val = mAlpha * mR[id] * tSup;
					aNode.val += val;
					(val > 0) ? p += val : m += val;
				}else if(mType == 5 || mType == 6){ // logistic
					double val = (mY[id] > 0) ? 1 / (1 + mR[id]) : 1 / (1 + mR[id]) - 1;
					val *= mAlpha * tSup;
					aNode.val += val;
					(val > 0) ? p += val : m += val;
				}
			}
		}
	}

	if(ans.size() > 0){
		aNode.inexact_pdb = ans;
	}
	aNode.val = fabs(aNode.val);


	if(mMinsup > aNode.supportSum) return true;
	else if(mMinsup < 1 && mMinsup > (double) aNode.supportSum / mN) return true;

	if(mType == 1 || mType == 3 || mType == 5){
		if(max(p, -m) < mMaxnode.val) return true;
	}else if(mType == 2 || mType == 4 || mType == 6){
		//SPPC計算 これより下は見る必要がないということ(ここが葉になる)
		if(max(p, -m) + mRadius * sqrt(aNode.supportSumsqr) < 1) return true;
	}

	return false;
}

/**
 * @fn
 *今見ているpdbの要素全ての一つ後ろの要素をcountermapに突っ込んで
 *それら全てをノードとして生やすためにproject_newする
 *
 */
void PrefixSpan::project(const TreeIter aNodeIter){
	//cout << "project" << endl;
	projectDB tPdb = aNodeIter->pdb;
	//porjectDB tInexactPdb = aNodeIter->inexact_pdb;

	// scan projected database

	//ここからI-Extension
	//大前提：Itemsetは重複するアイテムは存在しない．探索は飛びを幾つでも許すもののみ
	if(mFlagItemsetExist){
		map<Event, projectDB> tCounter;

		for(uint i = 0, size = tPdb.size(); i < size; ++i){
			uint id = tPdb[i].sequence_index;
			uint j = tPdb[i].event_index;
			if(mTransaction[id][j].first.size() - 1 > tPdb[i].itemset_index){
				uint k = tPdb[i].itemset_index + 1;
				for(; k < mTransaction[id][j].first.size(); ++k){
					Event tEvent;
					vector<uint> tItemset = mPattern.back().first;
					tItemset.push_back(mTransaction[id][j].first[k]);
					tEvent.first = tItemset;
					tEvent.second = mTransaction[id][j].second;
					Position tPos;
					tPos.sequence_index = id;
					tPos.event_index = j;
					tPos.itemset_index = k;
					tCounter[tEvent].push_back(tPos);
				}
			}
		}

		// project: next event
		Event tSaveEvent = mPattern.back();
		for(auto it = tCounter.begin(), end = tCounter.end(); it != end; ++it){
			mPattern.pop_back();
			mPattern.push_back(it->first);
			project_new(it->second, aNodeIter);
		}
		mPattern.pop_back();
		mPattern.push_back(tSaveEvent);
	}

	if(mPattern.size() < mMaxpat){
		//ここからS-Extension
		map<Event, projectDB> tCounter;

		if(mInterval < 0){
			map<Event, uint> dupCheck;
			for(uint i = 0, size = tPdb.size(); i < size; ++i){
				uint id = tPdb[i].sequence_index;
				uint trsize = mTransaction[id].size();
				//シーケンスの一つ隣のindexを取得
				uint j = tPdb[i].event_index + 1;
				//最初に出てきたものだけを採用
				for(; j < trsize; j++){
					Event tEvent;
					vector<uint> tItemset;
					tEvent.second = mTransaction[id][j].second;
					uint k = 0;

					Position tPosition;
					tPosition.sequence_index = id;
					tPosition.event_index = j;

					if(mTransaction[id][j].first.size() > 0){
						for(; k < mTransaction[id][j].first.size(); ++k){
							tItemset.push_back(mTransaction[id][j].first[k]);
							tEvent.first = tItemset;
							if(dupCheck.find(tEvent) == dupCheck.end()){
								dupCheck[tEvent] = 0;
								tPosition.itemset_index = k;
								tCounter[tEvent].push_back(tPosition);
							}
							tItemset.clear();
						}
					}else{
						tEvent.first = tItemset;
						if(dupCheck.find(tEvent) == dupCheck.end()){
							dupCheck[tEvent] = 0;
							tPosition.itemset_index = k;
							tCounter[tEvent].push_back(tPosition);
						}
					}
				}
				dupCheck.clear();
			}
		}else{
			vector<Position> dupCheck;
			for(uint i = 0, size = tPdb.size(); i < size; ++i){
				uint id = tPdb[i].sequence_index;
				uint trsize = mTransaction[id].size();
				//シーケンスの一つ隣のindexを取得
				uint j = tPdb[i].event_index + 1;
				uint j_tmp = j;

				for(; j < trsize; j++){
					//j-kは最大インターバルの計算
					if((int) j - j_tmp > mInterval){
						break;
					}

					Event tEvent;
					vector<uint> tItemset;
					tEvent.second = mTransaction[id][j].second;
					uint k = 0;

					//pdbに全く同じものが入らないようにする
					Position tPos;
					tPos.sequence_index = id;
					tPos.event_index = j;

					if(mTransaction[id][j].first.size() > 0){
						for(; k < mTransaction[id][j].first.size(); ++k){
							tItemset.push_back(mTransaction[id][j].first[k]);
							tEvent.first = tItemset;
							tPos.itemset_index = k;
							if(find(dupCheck.begin(), dupCheck.end(), tPos) == dupCheck.end()){
								dupCheck.push_back(tPos);
								tCounter[tEvent].push_back(tPos);
							}
							tItemset.clear();
						}
					}else{
						tPos.itemset_index = k;
						tEvent.first = tItemset;
						if(find(dupCheck.begin(), dupCheck.end(), tPos) == dupCheck.end()){
							dupCheck.push_back(tPos);
							tCounter[tEvent].push_back(tPos);
						}
					}
				}
			}
		}

		// project: next event
		for(auto it = tCounter.begin(), end = tCounter.end(); it != end; ++it){
			mPattern.push_back(it->first);
			project_new(it->second, aNodeIter);
			mPattern.pop_back();
		}
	}

}

/**
 * @fn
 *新しいノードを作成する
 *それが枝切りできるかチェックしてできるなら葉とする
 *できないなら木にこのノードを追加
 *UBを計算して展開すべきものかを見極める
 *葉にぶつかるまでprojectを掘る
 */
void PrefixSpan::project_new(projectDB& aPdb, const TreeIter aParent){
	//cout << "project_new" << endl;
	Node tNode;
	tNode.parent = aParent;
	tNode.pattern = mPattern;
	tNode.patternStr = pat2str();
	//cout << tNode.patternStr << endl;
	tNode.supportSum = 0;
	tNode.supportSumsqr = 0;
	tNode.pdb = aPdb;
	tNode.w = 0;
	tNode.val = 0;
	tNode.closed = true;
	tNode.pdb_size = sum_items_pdb(tNode.pdb);

	//リストの最後にNodeを挿入
	TreeIter tCurrent = mTree.insert(mTree.end(), tNode);

	//新しく作ったnodeが条件を満たすかどうかをしらべる
	//サポートのカウントもしている
	bool tFlag = calculate_new(*tCurrent);

	//ここでCloSpanチェック
	//tree全てのnodeと新しいnodeを比べて，包含関係にあってPDBが同じものがないかチェックする
	//詳細はCloSpanの論文参照
	if(mCloSpan == 1){
		//親と同じサポートなら親はclosedではない
		if(tCurrent->parent->supportSum == tCurrent->supportSum){
			tCurrent->parent->closed = false;
		}
		mFlagCScheckEnd = 0;

		checkProjectedDB(tCurrent);
		if(mFlagCScheckEnd == 2){
			mTree.pop_back();
			return;
		}
	}

	//親に子供のイテレータを追加
	aParent->child.push_back(tCurrent);

	if(tFlag){
		return;
	}else{
		//葉になるとき
		// S-Extention
		if(tNode.parent->pattern.size() != tNode.pattern.size()){
			tNode.maxChildLength = tNode.pattern.size();
			updateParentsMaxChildLength(tNode.parent, tNode.maxChildLength);
		}else{
			// I-Extension
			Event tLast = tNode.pattern.back();
			tNode.maxItemsetSize = tLast.first.size();
			updateParentsMaxChildItemsetSize(tNode.parent, tNode.maxItemsetSize);
		}
	}


	if(mType == 1 || mType == 3 || mType == 5){ // get_mMaxnode
		if(tCurrent->val > mMaxnode.val){
			mMaxnode = *tCurrent;
		}
	}else if(mType == 2 || mType == 4 || mType == 6){ // safe_screening
		//SPPC計算
		double score = tCurrent->val + mRadius * sqrt(tCurrent->supportSumsqr);
		if(score >= 1){
			//closed sequential mPatternかどうかのフラグ
			bool tAddFlag = true;
			//ここで同一サポートかつより長いパターンを残すようにするチェックを挟む
			uint tKey = tCurrent->supportSum;
			if(mCloSpan == 1 && mCheckClosed.find(tKey) != mCheckClosed.end()){
				for(auto itr = mCheckClosed[tKey].begin(); itr != mCheckClosed[tKey].end();){
					if(!(*itr)->closed){
						itr++;
						continue;
					}
					if(tCurrent->supportSum == (*itr)->supportSum){
						uint tWhichSub = isSubsequence(tCurrent->pattern, (*itr)->pattern);
						if(tWhichSub == 1){
							tAddFlag = false;
							break;
						}else if(tWhichSub == 2){
							mActive.erase(remove(mActive.begin(), mActive.end(), (*itr)), mActive.end());
							itr = mCheckClosed[tKey].erase(itr);
							continue;
						}else if(tCurrent->patternStr.compare((*itr)->patternStr)== 0){
							// 同一パターンがすでに存在していたらなら追加しない
							tAddFlag = false;
							break;
						}
					}
					itr++;
				}
			}

			if(tAddFlag){
				mCheckClosed[tKey].push_back(tCurrent);
				mActive.push_back(tCurrent);
			}

		}
		//if (score >= 1)で漏れたやつは切れないけどnonアクティブ あるかは不明
	}

	if(mFlagCScheckEnd == 0){
		project(tCurrent);
	}
}


/**
 * CloSpan用
 * Pruning可能性をチェックし可能ならpruningする
 */
void PrefixSpan::checkProjectedDB(const TreeIter aCurrentIter){
	// キーをProjedtedDBのidの和＋サポートにする(論文とは異なる)
	uint tKey = sumId(aCurrentIter->x) + aCurrentIter->supportSum + aCurrentIter->pdb_size;

	if(mCheckPDB.find(tKey) != mCheckPDB.end()){
		for(auto itr = mCheckPDB[tKey].begin(); itr != mCheckPDB[tKey].end();){
			// 比較するパターン間の系列長の差と比較先の最大葉ノード系列長が最大系列長を越える場合はPruning不可.
			if(abs((int)(aCurrentIter->pattern.size() - (*itr)->pattern.size())) + (*itr)->maxChildLength > mMaxpat){
				itr++;
				continue;
			}

			if(aCurrentIter->pdb_size == (*itr)->pdb_size){
				uint tWhichSub = isSubsequence(aCurrentIter->pattern, (*itr)->pattern);
				if(tWhichSub == 1){
					mFlagCScheckEnd = 2;
					return;
				}else if(tWhichSub == 2){
					//子ノードのつけかえ
					aCurrentIter->child = (*itr)->child;
					//不飽和ノードの親が持つ子ノードリストから，この不飽和ノードを削除
					(*itr)->parent->child.erase(find((*itr)->parent->child.begin(), (*itr)->parent->child.end(), *itr));

					//mActiveに存在するなら削除
					//これをしないとmActiveが持っているイテレータが未定義値を指すようになる危険性がある
					mActive.erase(remove(mActive.begin(), mActive.end(), (*itr)), mActive.end());


					// これをやるとSegmentation faultがおきてしまう(理由不明)
					// uint tKeyCheckClosed = sumId(aCurrentIter->x) + aCurrentIter->supportSum;
					// mCheckClosed[tKeyCheckClosed].erase(itr);

					//木から不飽和ノードを削除
					// mTree.erase(*itr);
					// 削除ではなくNonClosedってことにだけする
					(*itr)->closed = false;

					itr = mCheckPDB[tKey].erase(itr);

					//子ノードの持つ親ノードつけかえ
					for(auto it : aCurrentIter->child){
						it->parent = aCurrentIter;
					}

					//子ノード以下のパターンを全て更新
					for(auto it : aCurrentIter->child){
						childPatternUpdate(it);
					}
					mFlagCScheckEnd = 1;
					return;
				}
			}

			itr++;
		}

	}

	mCheckPDB[tKey].push_back(aCurrentIter);
	return;
}

/**
 * 子のパターンを更新
 */
void PrefixSpan::childPatternUpdate(const TreeIter aChildIter){
	Event tChildLast = aChildIter->pattern.back();

	Event tParentLast = mPattern.back();
	//I-Extensionかどうかのチェック
	bool tIE = false;
	if(tChildLast.first.size() > tParentLast.first.size()){
		mPattern.pop_back();
		tIE = true;
	}
	mPattern.push_back(tChildLast);
	aChildIter->pattern = mPattern;
	aChildIter->patternStr = pat2str();
	for(auto it : aChildIter->child){
		childPatternUpdate(it);
	}
	mPattern.pop_back();
	if(tIE) mPattern.push_back(tParentLast);
}
/**
 * 葉ノードの系列長を更新できる限り親，親の親．．．と伝えていく
 * 親の持っているMaxChildLenghtが上回ったら終了
 * 第一引数にはみたいノードを渡す
 */
void PrefixSpan::updateParentsMaxChildLength(const TreeIter aNodeIter,const uint aMaxChildLength){
	if(aNodeIter->maxChildLength >= aMaxChildLength){
		return;
	}

	// 更新
	aNodeIter->maxChildLength = aMaxChildLength;


	// 長さ１まで来たor最小系列長まで来たなら終わり
	if(aNodeIter->pattern.size()==1){
		return;
	}

	updateParentsMaxChildLength(aNodeIter->parent,aMaxChildLength);
}

/**
 * ノードのアイテムセットサイズをI-Extensionである限り親，親の親．．．と伝えていく
 * 親がS-Extensionまたは親の持っているMaxChildItemsetSizeが上回ったら終了
 * 第一引数にはみたいノードを渡す
 */
void PrefixSpan::updateParentsMaxChildItemsetSize(const TreeIter aNodeIter, const uint aMaxChildItemSize){
	if(aNodeIter->maxItemsetSize >= aMaxChildItemSize){
		return;
	}

	// 更新
	aNodeIter->maxItemsetSize = aMaxChildItemSize;

	Event tLastEvent = aNodeIter->pattern.back();
	if(tLastEvent.first.size()== 1){
		return;
	}

	// S-Extensionであったら終わり
	if(aNodeIter->parent->pattern.size() != aNodeIter->pattern.size()){
		return;
	}

	updateParentsMaxChildItemsetSize(aNodeIter->parent, aMaxChildItemSize);
}

/**
 * 木を辿る
 */
void PrefixSpan::searchTree_forMaxnode(const TreeIter aNodeIter){

	mPattern = aNodeIter->pattern;
	bool flag = calculate(*aNodeIter);
	if(flag){return;}

	if(aNodeIter->closed){
		// 最大値更新
		if(aNodeIter->val > mMaxnode.val) mMaxnode = *aNodeIter;
	}

	//子の探索
	if(aNodeIter->child.size() == 0){
		project(aNodeIter);
	}else{
		for(auto it : aNodeIter->child){
			searchTree_forMaxnode(it);
		}
	}

}

//木を辿る
void PrefixSpan::searchTree_forSafeScreening(const TreeIter aNodeIter){
	mPattern = aNodeIter->pattern;
	bool flag = calculate(*aNodeIter);
	if(flag){
		aNodeIter->w = 0;
		return;
	}

	if(aNodeIter->closed){
		//uj部分が和の絶対値で少し小さい
		double score = aNodeIter->val + mRadius * sqrt(aNodeIter->supportSumsqr);
		if(score < 1){
			aNodeIter->w = 0;
		}else{
			//closed sequential mPatternかどうかのフラグ
			bool tAddFlag = true;
			//ここで同一サポートかつより長いパターンを残すようにするチェックを挟む
			uint tKey = sumId(aNodeIter->x) + aNodeIter->supportSum;
			if(mCloSpan == 1 && mCheckClosed.find(tKey) != mCheckClosed.end()){
				for(auto itr = mCheckClosed[tKey].begin(); itr != mCheckClosed[tKey].end();){
					if(!(*itr)->closed){
						itr++;
						continue;
					}
					if(aNodeIter->supportSum == (*itr)->supportSum){
						uint tWhichSub = isSubsequence(aNodeIter->pattern, (*itr)->pattern);
						if(tWhichSub == 1){
							tAddFlag = false;
							break;
						}else if(tWhichSub == 2){
							mActive.erase(remove(mActive.begin(), mActive.end(), (*itr)), mActive.end());
							itr = mCheckClosed[tKey].erase(itr);
							continue;
						}else if(aNodeIter->patternStr.compare((*itr)->patternStr)== 0){
							// 同一パターンがすでに存在していたらなら追加しない
							tAddFlag = false;
							break;
						}
					}
					itr++;
				}
			}

			if(tAddFlag){
				mCheckClosed[tKey].push_back(aNodeIter);
				mActive.push_back(aNodeIter);
			}

		}
	}

	//子の探索
	if(aNodeIter->child.size() == 0){
		project(aNodeIter);
	}
	//}else{
		for(auto it : aNodeIter->child){
			searchTree_forSafeScreening(it);
		}
	//}
}


/**
 * @fn
 *こいつを使うことで根から枝を伸ばし切る
 *
 */
void PrefixSpan::get_maxnode(const vector<double> &_v, int solver){
	switch(solver){
		case 1: // svm
			mType = 1;
			break;		case 2: //lasso
			mType = 3;
			break;
		case 3: //logistic
			mType = 5;
			break;
		default:
			cerr << "Error: unknown solver type" << endl;
			exit(1);
			break;
	}
	mR = _v;
	mAlpha = 1;
	mMaxnode.val = 0;
	mPattern.clear();

	TreeIter root = mTree.begin();

	//アイテム集合の木を作成し，深さ優先探索．着目しているアイテム集合をIとする
	for(auto it : root->child){
		mPattern = (it->pattern);
		searchTree_forMaxnode(it);
		//project(it);
	}

}

//残差，α，半径を用いて枝切りできるかを調べる
void PrefixSpan::safe_screening(vector<double> &_v, double _alpha, double _radius, int solver){
	switch(solver){
		case 1: // svm
			mType = 2;
			break;
		case 2: // lasso
			mType = 4;
			break;
		case 3: //logistic
			mType = 6;
			break;
		default:
			cerr << "Error: unknown solver type" << endl;
			exit(1);
			break;
	}
	mR = _v;
	mAlpha = _alpha;
	mRadius = _radius;
	mPattern.clear();
	mActive.clear();
	mCheckClosed.clear();
	mCheckPDB.clear();

	TreeIter root = mTree.begin();

	//アイテム集合の木を作成し，深さ優先探索．着目しているアイテム集合をIとする
	for(auto it : root->child){
		mPattern = (it->pattern);
		searchTree_forSafeScreening(it);
		//project(it);
	}

}

//active pattern の test系列に対するsupportを計算する
vector<double> PrefixSpan::calcActivesupport(vector<Event>  &tPattern,const vector<uint>  &testId){
  vector<double> ans;

	TreeIter root = mTree.begin();
	TreeIter n = root;
	int flag = 0;
	//cout << "calc" << endl;
	searchTree_Activesupport(tPattern, root, n);
	//cout << "calc" << endl;


	vector<uint> id = n -> x;
	vector<double> tsupport = n -> support;
	//for(int j = 0;j < id.size(); j++){
		//cout << "id[" << j << "] = " << id[j] << endl;
		//cout << "support[" << j << "] = " << tsupport[j] << endl;
	//}

	for(int i = 0 ;i < testId.size();i++){
		flag=0;
		//cout << testId[i] << endl;
		for(int j = 0;j < id.size(); j++){
			if(id[j] == testId[i]){
				//cout << id[j] << " = " <<testId[i] << endl;
				ans.push_back(tsupport[j]);
				//cout << tsupport[i] << endl;
				flag = 1;
				break;
			}
		}
		if(flag == 0){
			ans.push_back(0);
			//cout << "support 0" << endl;
		}
	}
	return ans;
}

//再帰用
int PrefixSpan::searchTree_Activesupport(vector<Event>  &tPattern, TreeIter &root,TreeIter &n){
	//cout << "searchA" << endl;
	if((root->child).size() == 0){
		return 0;
	}
	for(auto it : root->child){
				//cout << it->patternStr << endl;
		//cout << "searchB" << endl;

		if(n != mTree.begin()){
			return 0;
		}
		mPattern = (it->pattern);
		if(tPattern == mPattern){

			n = it;
		}
		searchTree_Activesupport(tPattern, it, n);
		//cout << "searchC" << endl;
	}
	return 0;
}

#include "database.h"

void Database::read(const char *aFilename){
	ifstream tFile(aFilename);
	if(!tFile){
		cerr << "Error: cannot open file" << endl;
		exit(1);
	}
	uint flag=0;
	uint tItemSize = 0;
	double tLabel;
	string tLine;
	vector<Event> tSequence;

	//全データを読み込む
	while(getline(tFile, tLine)){
		tSequence.clear();
		stringstream ss1(tLine);
		ss1 >> tLabel;
		mY.push_back(tLabel);
		string eventstring;

		// ここで１レコードを読む
		while(ss1 >> eventstring){
			stringstream ss2(eventstring);
			Itemset tItemset;
			vector<uint> tItem;
			string itemstring;
			int tmp;
			//ここのループで1:1:1というある時間でのデータを読む
			while(getline(ss2, itemstring, ':')){
				if(contain(itemstring, '(')){
					string tString;
					stringstream ss3(itemstring);
					if(contain(itemstring, '_')){
						while(getline(ss3, tString, '_')){
							if(contain(tString, '(')){
								if(contain(tString, ')')) break;

								tString.erase(tString.begin());
							}else if(contain(tString, ')')){
								tString.pop_back();
							}
							tItemset.push_back(stoi(tString));
						}
					}else {
						itemstring.erase(itemstring.begin());
						itemstring.pop_back();

						tItemset.push_back(stoi(itemstring));
					}

				}else{
					tmp = stoi(itemstring);
					uint tVal = (tmp < 0) ? 0xffffffff : tmp; // wild card
					tItem.push_back(tVal);

				}
			}

			//全要素が同じフォーマットになっているかをチェックしている 1:1:1 と1:1みたいなものは共存できない.
			if(!tItemSize){
				tItemSize = tItem.size();
			}else{
				if(tItemSize != tItem.size()){
					cerr << "Format Error: different Event Size at line: " << mTransaction.size() << ", event: " << tSequence.size() << endl;
					exit(-1);
				}
			}

			Event tEvent;
			tEvent.first = tItemset;
			tEvent.second = tItem;
			if(flag==0){
				min=tItem;
				max=tItem;
				flag =1;
			}else{
				for(int i=0;i<min.size();i++){
					if(tItem[i] < min[i]){
						min[i]=tItem[i];
					}
					if(tItem[i] > max[i]){
						max[i]=tItem[i];
					}
				}
			}
			tSequence.push_back(tEvent);
		}
		mTransaction.push_back(tSequence);
	}
}

vector<vector<Event> > Database::get_transaction() const{
	return mTransaction;
}

vector<double> Database::get_y() const{
	return mY;
}

uint Database::get_maxGap() const{
	uint count=0;
	uint fmode=1;
	if(fmode==0){
		for(int i=0;i<min.size();i++){
			//cout << max[i]<< endl;
			//cout << min[i]<< endl;
			count+=max[i]-min[i];
		}
		return count;
	}else{
		return max.size();
	}
}

/* 入力形式として1次元のみアイテムセットを許す
 * 例　(1,2,3):10 はアイテムセット(1,2,3)と10の2次元イベントとする
 * アイテムセットが存在しない場合はEventの後ろのvectorだけ使用する
 * アイテムセットは必ず丸かっこで囲んでおくこと
 */
#ifndef DATABASE_H
#define DATABASE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

using uint = unsigned int;
using Itemset = vector<uint>;
using Event = pair<Itemset, vector<uint>>;
/*
struct Event{
	Itemset first;
	vector<uint> second;

	bool operator  == (const Event& eve)const{
		if(sizeof(first)==sizeof(eve.first)){
			for(int i =0; i < sizeof(first);i++){
				if(first[i]!=eve.first[i]){
					return false;
				}
			}
			return true;
		}
		return false;
	}
	bool operator  < (const Event& eve) const{
		if(sizeof(first) < sizeof(eve.first)){
			return true;
		}else{
			return false;
		}
		return false;
	}
};
*/

class Database{
private:
	vector<vector<Event> > mTransaction;
	vector<double> mY;
	vector<uint> min;
	vector<uint> max;
public:
	void read(const char *aFilename); //txtデータをlabelとtransactionの形にする．
	vector<vector<Event> > get_transaction() const; //transactionを取り出す.
	vector<double> get_y() const; //yを取り出す
	vector< pair<Event, double> > get_fList(); //fuzzyListを取り出す
  uint get_maxGap() const;
	template<class T, class U>
	//sの中にvが含まれていればTrueを返す
	bool contain(const std::basic_string<T>& s, const U& v){
		return s.find(v) != std::basic_string<T>::npos;
	}
};

#endif

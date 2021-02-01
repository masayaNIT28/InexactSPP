#ifndef PREFIXSPAN_H
#define PREFIXSPAN_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <set>

#include "database.h"

using namespace std;

using uint = unsigned int;
using uchar = unsigned char;


class PrefixSpan{
private:
	//何番目のアイテム集合の前から何番目かを表す
	struct Position{
		uint sequence_index;
		uint event_index;
		uint itemset_index;
		bool operator==(const Position& x){
			return (sequence_index == x.sequence_index && event_index == x.event_index && itemset_index == x.itemset_index);
		}
		bool operator!=(const Position& x){
			return (sequence_index != x.sequence_index || event_index != x.event_index || itemset_index != x.itemset_index);
		}
	};

	using projectDB = vector< Position >;

	struct Node{
		list<Node>::iterator parent;
		list<list<Node>::iterator> child;
		vector<Event> pattern;
		string patternStr;
		double supportSum;
		vector<double> support;
		projectDB pdb;
		vector<uint> x;
		double w;
		double val;
		int addLambda = -1;
		bool closed; //CloSpan用．残すべきならTrue clospanをしてないときは常にTrue
		uint pdb_size; //Clospan用.論文中におけるデータベースサイズ. pdbの数ではない.
		uint maxChildLength; //CloSpan用. 自ノードと子ノードの中の最大シケーンス長
		uint maxItemsetSize; //CloSpan用. 自ノード以下での最大アイテムセットサイズ. I-Extensionが区切りとなる.

		//inexact matching
		vector<uint> Inexact_idList;
		projectDB inexact_pdb;
		double supportSumsqr;
		vector<uint> Inexact_x;
		vector<double> Inexact_support;

	};


	using TreeIter = list<Node>::iterator;

	const uint MAXVAL = 0xffffffff;
	uint count=0;
	uint mItemSize;
	double mMinsup;
	uint mMaxpat;
	int mInterval;
	uint mSupMode;
	uint mCloSpan;
	int mType;
	vector<double> mR;
	double mAlpha;
	double mRadius;

	vector< vector<Event> > mTransaction; //ここにファイルから読み込んだデータが入る
	vector<Event> mPattern;
	list<Node> mTree;
	uint mFlagCScheckEnd = 0; //CloSpanチェック用
	unordered_map<uint, list<TreeIter>> mCheckPDB; //CloseSpanでのpruning確認用
	unordered_map<uint, list<TreeIter>> mCheckClosed; //closedなパターンの確認用
	uint mFlagItemsetExist = 0; //Itemsetチェック用
	Event mWildEvent;

	uint mMaxGap;
	uint mMinInexactsup;
	uint mMinitemdifferent;

	vector< vector<Event> > tTest;

	string pat2str(void);
	string pat2str(const vector<Event> aPattern);
	double calcSup(uint aId, vector<Event> aPattern);
	int calcPat(uint aId,vector<Event> aPattern);
	uint isSubsequence(const vector<Event> aSeq1, const vector<Event> aSeq2);
	bool isInclude(const Event aEvent1, const Event aEvent2);
	void childPatternUpdate(const TreeIter aChildIter);
	void updateParentsMaxChildLength(const TreeIter aNodeIter,const uint aMaxChildLength);
	void updateParentsMaxChildItemsetSize(const TreeIter aNodeIter, const uint aMaxChildItemSize);
	void searchTree_forMaxnode(const TreeIter aNodeIter);
	void searchTree_forSafeScreening(const TreeIter aNodeIter);
	void project(const TreeIter aNodeIter);
	void project_new(projectDB& aPdb, const TreeIter aParent);
	bool calculate(Node& aNode);
	bool calculate_new(Node& aNode);
	bool isParent(const TreeIter aNodeIter, const TreeIter aChildIter);
	void checkProjectedDB(const TreeIter aCurrentIter);
	uint sum_items_pdb(const projectDB aPdb);
	uint sumId(const vector<uint> aIdList);

	//inexact matching
	double calcInexactSup(projectDB tPdb, projectDB &minpDB, vector<Event> aPattern);

public:

	uint mN;
	vector<double> mY;
	Node mMaxnode;
	vector<TreeIter> mActive;

 	//コンストラクタ
	PrefixSpan(double aMinsup, uint aMaxpat, int aInterval,uint aSupMode, uint aCloSpan, uint aMaxGap, uint aMinitemdifferent) {
		mMinsup = aMinsup;
		mMaxpat = aMaxpat;
		mInterval = aInterval;
		mSupMode = aSupMode;
		mCloSpan = aCloSpan;
		mMaxGap = aMaxGap;
		mMinitemdifferent = aMinitemdifferent;
		mMinInexactsup = 1.0 - aMinitemdifferent/(aMaxGap*aMaxpat);
	};

	void init(const vector< vector<Event> > aTransaction,const vector<double> aY);
	void get_maxnode(const vector<double> &_v, int solver);
	void safe_screening(vector<double> &_v, double _alpha, double _radius, int solver);
	void printTree(string aFilename, double bias);
	void printTree(string aFilename);
	void add_history(uint aLambda);
	int calcSup(const vector<Event> &aTransaction,const vector<Event> &aPattern) const;
	double calcInexactSup(const vector<Event> &aTransaction, const vector<Event> aPattern) const;
	vector<double> calcActivesupport(vector<Event>  &tPattern,const vector<uint>  &testId);
	int searchTree_Activesupport(vector<Event>  &tPattern, TreeIter &root,TreeIter &n);

};

#endif

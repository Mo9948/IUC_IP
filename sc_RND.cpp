/*_________________________________________________________________________________________________________
|                                                                                                          |
|                 IP solution methods for the Max IUC problem (RND Instances -- Strong VIs added)          |
|                                                                                                          |
|                 Copyright (c) 2019 Seyedmohammadhossein Hosseinian. All rights reserved.                 |
|                                                                                                          |
|__________________________________________________________________________________________________________|


***  READ ME  ********************************************************************************************

(1) Please delete the header of the instance files before running the code.

***********************************************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <Windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <queue>
#include <stack>
#include <string>
#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <math.h>
#include <ppl.h>

#define eta_fH 0.5
#define eta_S 0

using namespace std;

ILOSTLBEGIN

#define Replace_OT_fH 1
#define addVI_fourHole 0
#define addVI_star 0
#define addVI_wheelFan 1
#define printVI 0
#define autoCut 1
#define dispProcess 1
#define printLP 0
#define timeLimit 72000
#define minCycleSize 7
#define minPathSize 7
#define preprocessTime 1
#define numberWheelFan 50

#pragma region "Time Record"

double get_wall_time() {
	LARGE_INTEGER time, freq;
	if (!QueryPerformanceFrequency(&freq)) {
		return 0;
	}
	if (!QueryPerformanceCounter(&time)) {
		return 0;
	}
	return (double)time.QuadPart / freq.QuadPart;
}

double get_cpu_time() {
	FILETIME a, b, c, d;
	if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0) {
		return
			(double)(d.dwLowDateTime |
			((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
	}
	else {
		return 0;
	}
}

#pragma endregion

#pragma region "Utility"

struct ot {
	int vertex_1;
	int vertex_2;
	int vertex_3;
};

bool isOT(int** Adj, int i, int j, int k) {  //is Open Triangle?
	int temp = Adj[i][j] + Adj[i][k] + Adj[j][k];
	if (temp == 2) return true;
	else return false;
}

bool isfourHole(int** Adj, int i, int j, int k, int l) {
	int tempI = Adj[i][j] + Adj[i][k] + Adj[i][l];
	int tempJ = Adj[j][i] + Adj[j][k] + Adj[j][l];
	int tempK = Adj[k][i] + Adj[k][j] + Adj[k][l];
	int tempL = Adj[l][i] + Adj[l][j] + Adj[l][k];
	if (tempI == 2 && tempJ == 2 && tempK == 2 && tempL == 2) return true;
	else return false;
}

struct ver {
	int n;
	int degCol;
};

bool cmpdegINC(ver const & a, ver const & b) { return (a.degCol < b.degCol); }
bool cmpdegDEC(ver const & a, ver const & b) { return (a.degCol > b.degCol); }

// *** GENERATE 4-HOLES ***

struct fourHole {
	int number;
	int** main;
};

fourHole makeFourHole(int** Adj, int N) {
	vector<vector<int>> AllFourHole;
	vector<int> tempHole;
	int tempI;
	int tempJ;
	int tempK;
	int tempL;
	for (int i = 0; i < N - 3; i++) {
		for (int j = i + 1; j < N - 2; j++) {
			for (int k = j + 1; k < N - 1; k++) {
				for (int l = k + 1; l < N; l++) {
					tempI = Adj[i][j] + Adj[i][k] + Adj[i][l];
					tempJ = Adj[j][i] + Adj[j][k] + Adj[j][l];
					tempK = Adj[k][i] + Adj[k][j] + Adj[k][l];
					tempL = Adj[l][i] + Adj[l][j] + Adj[l][k];
					if (tempI == 2 && tempJ == 2 && tempK == 2 && tempL == 2) {
						tempHole.clear();
						tempHole.push_back(i);
						tempHole.push_back(j);
						tempHole.push_back(k);
						tempHole.push_back(l);
						AllFourHole.push_back(tempHole);
					}
				}
			}
		}
	}
	fourHole AllfH;
	AllfH.number = AllFourHole.size();
	int** temMat = new int*[AllFourHole.size()];
	for (int i = 0; i < AllFourHole.size(); i++) {
		temMat[i] = new int[4];
	}
	for (int i = 0; i < AllFourHole.size(); i++) {
		for (int j = 0; j < 4; j++) {
			temMat[i][j] = AllFourHole[i][j];
		}
	}
	AllfH.main = temMat;
	return AllfH;
}

// *** GENERATE STARS ***

struct star {
	int number;
	int* sizes;
	int** main;
};

void checkStar(int** Adj, const vector<int> & Star) {
	bool isStar = TRUE;
	for (int i = 1; i < Star.size(); i++) {
		if (Adj[Star[0]][Star[i]] != 1) {
			isStar = FALSE;
			break;
		}
	}
	int neigh;
	for (int i = 1; i < Star.size(); i++) {
		neigh = 0;
		for (int j = 1; j < Star.size(); j++) {
			if (Adj[Star[i]][Star[j]] == 1) neigh++;
		}
		if (neigh != 0) isStar = FALSE;
		break;
	}
	if (isStar == FALSE) cout << "ERROR --- NOT a Star!" << endl;
};

star makeStar(int** Adj, int N, int v) {
	vector<vector<int>> vAllStar;
	vector<ver> NEIG;
	ver tempVer;
	for (int i = 0; i < N; i++) {
		if (Adj[v][i] == 1) {
			tempVer.n = i;
			NEIG.push_back(tempVer);
		}
	}
	int d;
	for (int i = 0; i < NEIG.size(); i++) {
		d = 0;
		for (int j = 0; j < NEIG.size(); j++) {
			if (Adj[NEIG[i].n][NEIG[j].n] == 1) d++;
		}
		NEIG[i].degCol = d;
	}
	//sort(NEIG.begin(), NEIG.end(), cmpdegINC);
	sort(NEIG.begin(), NEIG.end(), cmpdegDEC);
	vector<int> list;
	for (int i = 0; i < NEIG.size(); i++) {
		list.push_back(NEIG[i].n);
	}
	vector<int> alaki(1);
	alaki[0] = -1;
	vAllStar.push_back(alaki);
	int maxNo = 0;
	for (vector<int>::iterator itr_1 = list.begin(); itr_1 != list.end(); itr_1++) {
		int k = 0;
		bool has_neighbor = TRUE;
		while (has_neighbor == TRUE) {
			if (vAllStar[k].size() == 1) {
				has_neighbor = FALSE;
			}
			else {
				bool connected = FALSE;
				for (vector<int>::iterator itr_2 = vAllStar[k].begin(); itr_2 != vAllStar[k].end(); itr_2++) {
					if (*itr_2 != -1) {
						if (Adj[*itr_1][*itr_2] == 1) {
							k += 1;
							connected = TRUE;
							if (k > maxNo) {
								maxNo = k;
								vAllStar.push_back(alaki);
							}
							break;
						}
					}
				}
				if (connected == FALSE) {
					has_neighbor = FALSE;
				}
			}
		}
		vAllStar[k].push_back(*itr_1);
	}
	for (int i = 0; i < vAllStar.size(); i++) {
		vAllStar[i][0] = v;
		checkStar(Adj, vAllStar[i]);
	}
	vector<int> TEMP;
	for (int i = 0; i < vAllStar.size(); i++) {
		if (vAllStar[i].size() > 3) {
			TEMP.push_back(i);
		}
	}
	star AllSt;
	AllSt.number = TEMP.size();
	int* tempArray = new int[TEMP.size()];
	for (int i = 0; i < TEMP.size(); i++) {
		tempArray[i] = vAllStar[TEMP[i]].size();
	}
	AllSt.sizes = tempArray;
	int** temMat = new int*[TEMP.size()];
	for (int i = 0; i < TEMP.size(); i++) {
		temMat[i] = new int[vAllStar[TEMP[i]].size()];
	}
	for (int i = 0; i < TEMP.size(); i++) {
		for (int j = 0; j < vAllStar[TEMP[i]].size(); j++) {
			temMat[i][j] = vAllStar[TEMP[i]][j];
		}
	}
	AllSt.main = temMat;
	return AllSt;
}

// *** GENERATE WHEELS AND FANS ***

struct cyclePath {
	int numberC;
	int* sizesC;
	int** mainC;
	int numberP;
	int* sizesP;
	int** mainP;
};

struct wheelFan {
	int numberW;
	int* sizesW;
	int** mainW;
	int numberF;
	int* sizesF;
	int** mainF;
};

int* degreeLabeling(int** Adj, int N) {
	int *degree = new int[N];
	int *color = new int[N];
	int *label = new int[N];
	for (int i = 0; i < N; i++) {
		label[i] = -1;
		color[i] = 0;
		degree[i] = 0;
		for (int j = 0; j < N; j++) {
			if (Adj[i][j] == 1) degree[i]++;
		}
	}
	int min_degree;
	int v;
	for (int t = 0; t < N; t++) {
		v = NULL;
		min_degree = N;
		for (int i = 0; i < N; i++) {
			if (color[i] == 0 && degree[i] < min_degree) {
				v = i;
				min_degree = degree[i];
			}
		}
		label[v] = t;
		color[v] = 1;
		for (int i = 0; i < N; i++) {
			if (Adj[v][i] == 1 && color[i] == 0) {
				degree[i]--;
			}
		}
	}
	delete[] color;
	delete[] degree;
	return label;
}

int** makeNEWAdj(int** Adj, int N, int* label) {
	int** newAdj = new int*[N];
	for (int i = 0; i < N; i++) {
		newAdj[i] = new int[N];
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			newAdj[label[i]][label[j]] = Adj[i][j];
		}
	}
	return newAdj;
}

stack<vector<int>> makeTriplets(int** newAdj, int N) {
	stack<vector<int>> T;
	for (int i = 0; i < N - 2; i++) {
		for (int j = i + 1; j < N - 1; j++) {
			for (int k = j + 1; k < N; k++) {
				if (newAdj[i][j] == 1 && newAdj[i][k] == 1 && newAdj[j][k] != 1) {
					vector<int> temp;
					temp.push_back(j);
					temp.push_back(i);
					temp.push_back(k);
					T.push(temp);
				}
			}
		}
	}
	return T;
}

bool checkSingleCycle(int** Adj, const vector<int> & Cyc) {
	bool isCycle = TRUE;
	if (Adj[Cyc[0]][Cyc.back()] != 1) isCycle = FALSE;
	for (int i = 0; i < Cyc.size() - 1; i++) {
		if (Adj[Cyc[i]][Cyc[i + 1]] != 1) {
			isCycle = FALSE;
			break;
		}
	}
	int neigh;
	for (int i = 0; i < Cyc.size(); i++) {
		neigh = 0;
		for (int j = 0; j < Cyc.size(); j++) {
			if (Adj[Cyc[i]][Cyc[j]] > 0) neigh++;
		}
		if (neigh != 2) {
			isCycle = FALSE;
			break;
		}
	}
	//if (isCycle == FALSE) cout << "ERROR! -- NOT A Cycle!" << endl;
	return isCycle;
}

bool checkSinglePath(int** Adj, const vector<int> & Path) {
	bool isPath = TRUE;
	for (int i = 0; i < Path.size() - 1; i++) {
		if (Adj[Path[i]][Path[i + 1]] != 1) {
			isPath = FALSE;
			break;
		}
	}
	int neigh;
	for (int i = 0; i < Path.size(); i++) {
		neigh = 0;
		for (int j = 0; j < Path.size(); j++) {
			if (Adj[Path[i]][Path[j]] > 0) neigh++;
		}
		if (i == 0 || i == Path.size() - 1) {
			if (neigh != 1) {
				isPath = FALSE;
				break;
			}
		}
		else {
			if (neigh != 2) {
				isPath = FALSE;
				break;
			}
		}
	}
	//if (isPath == FALSE) cout << "ERROR! -- NOT A Path!" << endl;
	return isPath;
}

void checkAllWheelFan(int** Adj, const wheelFan & WF) {
	bool allGood = TRUE;
	for (int r = 0; r < WF.numberW; r++) {
		bool isWheel = TRUE;
		for (int i = 1; i < WF.sizesW[r]; i++) {
			if (Adj[WF.mainW[r][0]][WF.mainW[r][i]] != 1) {
				isWheel = FALSE;
				break;
			}
		}
		if (Adj[WF.mainW[r][1]][WF.mainW[r][(WF.sizesW[r]) - 1]] != 1) isWheel = FALSE;
		for (int i = 1; i < WF.sizesW[r] - 1; i++) {
			if (Adj[WF.mainW[r][i]][WF.mainW[r][i + 1]] != 1) {
				isWheel = FALSE;
				break;
			}
		}
		int deg;
		for (int i = 1; i < WF.sizesW[r]; i++) {
			deg = 0;
			for (int j = 1; j < WF.sizesW[r]; j++) {
				if (Adj[WF.mainW[r][i]][WF.mainW[r][j]] == 1) deg++;
			}
			if (deg != 2) {
				isWheel = FALSE;
				break;
			}
		}
		if (isWheel == FALSE) allGood = FALSE;
	}
	//
	for (int r = 0; r < WF.numberF; r++) {
		bool isFan = TRUE;
		for (int i = 1; i < WF.sizesF[r]; i++) {
			if (Adj[WF.mainF[r][0]][WF.mainF[r][i]] != 1) {
				isFan = FALSE;
				break;
			}
		}
		for (int i = 1; i < WF.sizesF[r] - 1; i++) {
			if (Adj[WF.mainF[r][i]][WF.mainF[r][i + 1]] != 1) {
				isFan = FALSE;
				break;
			}
		}
		int deg;
		for (int i = 1; i < WF.sizesF[r]; i++) {
			deg = 0;
			for (int j = 1; j < WF.sizesF[r]; j++) {
				if (Adj[WF.mainF[r][i]][WF.mainF[r][j]] == 1) deg++;
			}
			if (i == 1 || i == WF.sizesF[r] - 1) {
				if (deg != 1) {
					isFan = FALSE;
					break;
				}
			}
			else {
				if (deg != 2) {
					isFan = FALSE;
					break;
				}
			}
		}
		if (isFan == FALSE) allGood = FALSE;
	}
	if (allGood == TRUE) {
		cout << endl;
		cout << "   *** ALL WHEEL AND FAN SUBGRAPHS VERIFIED! ***   " << endl;
		cout << endl;
	}
	else {
		cout << " ----- ERROR ! SOMETHING WRONG ----- " << endl;
	}
}

cyclePath makeCyclePath(int** vAdj, int vN) {
	int* label = degreeLabeling(vAdj, vN);
	int** newAdj = makeNEWAdj(vAdj, vN, label);
	vector<vector<int>> C;
	vector<vector<int>> P;
	stack<vector<int>> T = makeTriplets(newAdj, vN);
	//cout << "# initial triplets: " << T.size();
	int counter = 0;
	double startTime = get_cpu_time();
	double endTime;
	vector<int> temp;
	while (T.size() != 0 && counter < numberWheelFan) {
		temp.clear();
		for (int i = 0; i < T.top().size(); i++) {
			temp.push_back(T.top()[i]);
		}
		T.pop();
		bool hope = FALSE;
		for (int v = 0; v < vN; v++) {
			if (v > temp[1] && newAdj[v][temp.back()] == 1) {
				bool connected = FALSE;
				for (int j = 1; j < temp.size() - 1; j++) {
					if (newAdj[v][temp[j]] == 1) {
						connected = TRUE;
						break;
					}
				}
				if (connected == FALSE) {
					if (newAdj[v][temp[0]] == 1) {
						temp.push_back(v);
						if ((temp.size() >= minCycleSize) && (temp.size() % 6 != 0)) {
							bool isCycle = checkSingleCycle(newAdj, temp);
							if (isCycle == TRUE) {
								C.push_back(temp);
								hope = TRUE;
								counter++;
								if (counter == numberWheelFan) break;
							}
						}
					}
					else {
						temp.push_back(v);
						T.push(temp);
						hope = TRUE;
					}
				}
			}
		}
		if (hope == FALSE && temp.size() >= minPathSize && temp.size() % 3 != 2) {
			bool isPath = checkSinglePath(newAdj, temp);
			if (isPath == TRUE) {
				P.push_back(temp);
				counter++;
				if (counter == numberWheelFan) break;
			}
		}
		endTime = get_cpu_time();
		if ((endTime - startTime) > preprocessTime) break;
	}
	endTime = get_cpu_time();
	int* backLabel = new int[vN];
	for (int i = 0; i < vN; i++) {
		backLabel[label[i]] = i;
	}
	delete[] label;
	for (int i = 0; i < vN; i++) delete[] newAdj[i];
	delete[] newAdj;
	for (int i = 0; i < C.size(); i++) {
		for (int j = 0; j < C[i].size(); j++) {
			C[i][j] = backLabel[C[i][j]];
		}
	}
	for (int i = 0; i < P.size(); i++) {
		for (int j = 0; j < P[i].size(); j++) {
			P[i][j] = backLabel[P[i][j]];
		}
	}
	delete[] backLabel;
	cyclePath CP;
	CP.numberC = C.size();
	int* tempArrayC = new int[C.size()];
	for (int i = 0; i < C.size(); i++) {
		tempArrayC[i] = C[i].size();
	}
	CP.sizesC = tempArrayC;
	int** temMatC = new int*[C.size()];
	for (int i = 0; i < C.size(); i++) {
		temMatC[i] = new int[C[i].size()];
	}
	for (int i = 0; i < C.size(); i++) {
		for (int j = 0; j < C[i].size(); j++) {
			temMatC[i][j] = C[i][j];
		}
	}
	CP.mainC = temMatC;
	CP.numberP = P.size();
	int* tempArrayP = new int[P.size()];
	for (int i = 0; i < P.size(); i++) {
		tempArrayP[i] = P[i].size();
	}
	CP.sizesP = tempArrayP;
	int** temMatP = new int*[P.size()];
	for (int i = 0; i < P.size(); i++) {
		temMatP[i] = new int[P[i].size()];
	}
	for (int i = 0; i < P.size(); i++) {
		for (int j = 0; j < P[i].size(); j++) {
			temMatP[i][j] = P[i][j];
		}
	}
	CP.mainP = temMatP;
	return CP;
}

wheelFan makeAllWheelFan(int** Adj, int N) {
	vector<vector<int>> AllW;
	vector<vector<int>> AllF;
	vector<int> temp;
	for (int v = 0; v < N; v++) {
		vector<int> NEIG;
		for (int i = 0; i < N; i++) {
			if (Adj[v][i] == 1) {
				NEIG.push_back(i);
			}
		}
		int** vAdj = new int*[NEIG.size()];
		for (int i = 0; i < NEIG.size(); i++) {
			vAdj[i] = new int[NEIG.size()];
		}
		for (int i = 0; i < NEIG.size(); i++) {
			for (int j = 0; j < NEIG.size(); j++) {
				vAdj[i][j] = Adj[NEIG[i]][NEIG[j]];
			}
		}
		cyclePath CP = makeCyclePath(vAdj, NEIG.size());
		for (int i = 0; i < CP.numberC; i++) {
			temp.clear();
			temp.push_back(v);
			for (int j = 0; j < CP.sizesC[i]; j++) {
				temp.push_back(NEIG[CP.mainC[i][j]]);
			}
			AllW.push_back(temp);
		}
		for (int i = 0; i < CP.numberP; i++) {
			temp.clear();
			temp.push_back(v);
			for (int j = 0; j < CP.sizesP[i]; j++) {
				temp.push_back(NEIG[CP.mainP[i][j]]);
			}
			AllF.push_back(temp);
		}
		for (int i = 0; i < NEIG.size(); i++) delete[] vAdj[i];
		delete[] vAdj;
		delete[] CP.sizesC;
		for (int i = 0; i < CP.numberC; i++) delete[] CP.mainC[i];
		delete[] CP.mainC;
		delete[] CP.sizesP;
		for (int i = 0; i < CP.numberP; i++) delete[] CP.mainP[i];
		delete[] CP.mainP;
	}
	wheelFan WF;
	WF.numberW = AllW.size();
	int* tempArrayW = new int[AllW.size()];
	for (int i = 0; i < AllW.size(); i++) {
		tempArrayW[i] = AllW[i].size();
	}
	WF.sizesW = tempArrayW;
	int** temMatW = new int*[AllW.size()];
	for (int i = 0; i < AllW.size(); i++) {
		temMatW[i] = new int[AllW[i].size()];
	}
	for (int i = 0; i < AllW.size(); i++) {
		for (int j = 0; j < AllW[i].size(); j++) {
			temMatW[i][j] = AllW[i][j];
		}
	}
	WF.mainW = temMatW;
	WF.numberF = AllF.size();
	int* tempArrayF = new int[AllF.size()];
	for (int i = 0; i < AllF.size(); i++) {
		tempArrayF[i] = AllF[i].size();
	}
	WF.sizesF = tempArrayF;
	int** temMatF = new int*[AllF.size()];
	for (int i = 0; i < AllF.size(); i++) {
		temMatF[i] = new int[AllF[i].size()];
	}
	for (int i = 0; i < AllF.size(); i++) {
		for (int j = 0; j < AllF[i].size(); j++) {
			temMatF[i][j] = AllF[i][j];
		}
	}
	WF.mainF = temMatF;
	return WF;
}

void checkSingleWheel(int** Adj, int* wheel, int sizeW) {
	bool isWheel = TRUE;
	for (int i = 1; i < sizeW; i++) {
		if (Adj[wheel[0]][wheel[i]] != 1) {
			isWheel = FALSE;
			break;
		}
	}
	if (isWheel == TRUE) {
		if (Adj[wheel[1]][wheel[sizeW - 1]] != 1) isWheel = FALSE;
		for (int i = 1; i < sizeW - 1; i++) {
			if (Adj[wheel[i]][wheel[i + 1]] != 1) {
				isWheel = FALSE;
				break;
			}
		}
		int neigh;
		for (int i = 1; i < sizeW; i++) {
			neigh = 0;
			for (int j = 1; j < sizeW; j++) {
				if (Adj[wheel[i]][wheel[j]] == 1) neigh++;
			}
			if (neigh != 2) {
				isWheel = FALSE;
				break;
			}
		}
	}
	if (isWheel == FALSE) cout << "ERROR! -- NOT A WHEEL!" << endl;
};

void checkSingleFan(int** Adj, int* fan, int sizeF) {
	bool isFan = TRUE;
	for (int i = 1; i < sizeF; i++) {
		if (Adj[fan[0]][fan[i]] != 1) {
			isFan = FALSE;
			break;
		}
	}
	if (isFan == TRUE) {
		for (int i = 1; i < sizeF - 1; i++) {
			if (Adj[fan[i]][fan[i + 1]] != 1) {
				isFan = FALSE;
				break;
			}
		}
		int neigh;
		for (int i = 1; i < sizeF; i++) {
			neigh = 0;
			for (int j = 1; j < sizeF; j++) {
				if (Adj[fan[i]][fan[j]] == 1) neigh++;
			}
			if (i == 1 || i == sizeF - 1) {
				if (neigh != 1) {
					isFan = FALSE;
					break;
				}
			}
			else {
				if (neigh != 2) {
					isFan = FALSE;
					break;
				}
			}
		}
	}
	if (isFan == FALSE) cout << "ERROR! -- NOT A FAN!" << endl;
};

void checkAllWheelFanAgain(int** Adj, wheelFan WF) {
	for (int i = 0; i < WF.numberW; i++) {
		checkSingleWheel(Adj, WF.mainW[i], WF.sizesW[i]);
	}
	for (int i = 0; i < WF.numberF; i++) {
		checkSingleFan(Adj, WF.mainF[i], WF.sizesF[i]);
	}
}

// *** REPLACE OTs WITH 4-HOLEs ***

struct rep {
	int orgOT;
	int OTnum;
	int** OTmat;
	int fHnum;
	int** fHmat;
};

rep replace(int** Adj, int N) {
	rep OTfH;
	OTfH.fHnum = 0;
	OTfH.fHmat = NULL;
	vector<vector<int>> OTlist;
	vector<vector<int>> fHlist;
	vector<int> temp;
	for (int i = 0; i < N - 2; i++) {
		for (int j = i + 1; j < N - 1; j++) {
			for (int k = j + 1; k < N; k++) {
				if (isOT(Adj, i, j, k)) {
					temp.clear();
					temp.push_back(i);
					temp.push_back(j);
					temp.push_back(k);
					OTlist.push_back(temp);
				}
			}
		}
	}
	OTfH.orgOT = OTlist.size();
	int* OTcheck = new int[OTlist.size()];
	for (int i = 0; i < OTlist.size(); i++) {
		OTcheck[i] = 1;
	}
	for (int i = 0; i < N - 3; i++) {
		for (int j = i + 1; j < N - 2; j++) {
			for (int k = j + 1; k < N - 1; k++) {
				for (int l = k + 1; l < N; l++) {
					int tempI = Adj[i][j] + Adj[i][k] + Adj[i][l];
					int tempJ = Adj[j][i] + Adj[j][k] + Adj[j][l];
					int tempK = Adj[k][i] + Adj[k][j] + Adj[k][l];
					int tempL = Adj[l][i] + Adj[l][j] + Adj[l][k];
					if (tempI == 2 && tempJ == 2 && tempK == 2 && tempL == 2) {
						temp.clear();
						temp.push_back(i);
						temp.push_back(j);
						temp.push_back(k);
						temp.push_back(l);
						fHlist.push_back(temp);
					}
				}
			}
		}
	}
	if (fHlist.size() > 0) {
		int* fHcheck = new int[fHlist.size()];
		for (int i = 0; i < fHlist.size(); i++) {
			fHcheck[i] = 0;
		}
		vector<vector<int>> OT_Adj;		//Adjacency list for OTs (i.e., for each OT, it gives the list of 4-holes involving that OT)
		vector<int> alaki;
		alaki.push_back(-1);
		for (int i = 0; i < OTlist.size(); i++) {
			OT_Adj.push_back(alaki);
		}
		vector<vector<int>> fH_Adj;		//Adjacency list for 4-holes (i.e., for each 4-hole, it gives the list of OTs that it contains)
		for (int i = 0; i < fHlist.size(); i++) {
			fH_Adj.push_back(alaki);
		}
		for (int r = 0; r < fHlist.size(); r++) {
			int i, j, k, l;
			i = fHlist[r][0];
			j = fHlist[r][1];
			k = fHlist[r][2];
			l = fHlist[r][3];
			//i
			int iSTART = 0;
			int iEND = NULL;
			while (OTlist[iSTART][0] != i) {
				iSTART++;
			}
			iEND = iSTART;
			while (OTlist[iEND][0] == i) {
				iEND++;
				if (iEND == OTlist.size()) break;
			}
			iEND--;
			//ij
			int jSTART = iSTART;
			int jEND = NULL;
			while (OTlist[jSTART][1] != j) {
				jSTART++;
			}
			jEND = jSTART;
			while (OTlist[jEND][1] == j) {
				jEND++;
				if (jEND == OTlist.size()) break;
			}
			jEND--;
			//ijk
			for (int row = jSTART; row <= jEND; row++) {
				if (OTlist[row][2] == k) {
					fH_Adj[r].push_back(row);
					OT_Adj[row].push_back(r);
					break;
				}
			}
			//ijl
			for (int row = jSTART; row <= jEND; row++) {
				if (OTlist[row][2] == l) {
					fH_Adj[r].push_back(row);
					OT_Adj[row].push_back(r);
					break;
				}
			}
			//ik
			int kSTART = iSTART;
			int kEND = NULL;
			while (OTlist[kSTART][1] != k) {
				kSTART++;
			}
			kEND = kSTART;
			while (OTlist[kEND][1] == k) {
				kEND++;
				if (kEND == OTlist.size()) break;
			}
			kEND--;
			//ikl
			for (int row = kSTART; row <= kEND; row++) {
				if (OTlist[row][2] == l) {
					fH_Adj[r].push_back(row);
					OT_Adj[row].push_back(r);
					break;
				}
			}
			//j
			jSTART = 0;
			jEND = NULL;
			while (OTlist[jSTART][0] != j) {
				jSTART++;
			}
			jEND = jSTART;
			while (OTlist[jEND][0] == j) {
				jEND++;
				if (jEND == OTlist.size()) break;
			}
			jEND--;
			//jk
			kSTART = jSTART;
			kEND = NULL;
			while (OTlist[kSTART][1] != k) {
				kSTART++;
			}
			kEND = kSTART;
			while (OTlist[kEND][1] == k) {
				kEND++;
				if (kEND == OTlist.size()) break;
			}
			kEND--;
			//jkl
			for (int row = kSTART; row <= kEND; row++) {
				if (OTlist[row][2] == l) {
					fH_Adj[r].push_back(row);
					OT_Adj[row].push_back(r);
					break;
				}
			}
		}
		// The replacement part starts here
		vector<int> notCovered;
		for (int r = 0; r < fHlist.size(); r++) {
			notCovered.push_back(fH_Adj[r].size() - 1);
			if (notCovered[r] != 4) cout << " *** ERROR in making fH_Adj ! *** " << endl;
		}
		vector<ver> forSort;
		ver verTemp;
		for (int r = 0; r < fHlist.size(); r++) {
			verTemp.n = r;
			verTemp.degCol = notCovered[r];
			forSort.push_back(verTemp);
		}
		int maxNotCoverted = 4;
		int currentfH;
		while (maxNotCoverted > 0) {
			currentfH = forSort.front().n;
			fHcheck[currentfH] = 1;
			//for (int p = 1; p < fH_Adj[currentfH].size(); p++) {
			for (int p = 1; p < 5; p++) {		// fH_Adj[currentfH].size() is 5 for all
				for (int q = 1; q < OT_Adj[fH_Adj[currentfH][p]].size(); q++) {
					if (OTcheck[fH_Adj[currentfH][p]] == 1 && notCovered[OT_Adj[fH_Adj[currentfH][p]][q]] > 0) {
						notCovered[OT_Adj[fH_Adj[currentfH][p]][q]]--;
					}
				}
				OTcheck[fH_Adj[currentfH][p]] = 0;
			}
			forSort.erase(forSort.begin());
			for (int s = 0; s < forSort.size(); s++) {
				forSort[s].degCol = notCovered[forSort[s].n];
			}
			sort(forSort.begin(), forSort.end(), cmpdegDEC);
			maxNotCoverted = forSort.front().degCol;
		}
		int used_fH = 0;
		for (int i = 0; i < fHlist.size(); i++) {
			if (fHcheck[i] == 1) used_fH++;
		}
		OTfH.fHnum = used_fH;
		int** final_fH = new int*[used_fH];
		for (int i = 0; i < used_fH; i++) {
			final_fH[i] = new int[4];
		}
		int z = 0;
		for (int i = 0; i < fHlist.size(); i++) {
			if (fHcheck[i] == 1) {
				final_fH[z][0] = fHlist[i][0];
				final_fH[z][1] = fHlist[i][1];
				final_fH[z][2] = fHlist[i][2];
				final_fH[z][3] = fHlist[i][3];
				z++;
			}
		}
		OTfH.fHmat = final_fH;

		delete[] fHcheck;

	}	//if fHlist.size()>0 ends here
	int remained_OT = 0;
	for (int i = 0; i < OTlist.size(); i++) {
		if (OTcheck[i] == 1) remained_OT++;
	}
	OTfH.OTnum = remained_OT;
	int** final_OT = new int*[remained_OT];
	for (int i = 0; i < remained_OT; i++) {
		final_OT[i] = new int[3];
	}
	int z = 0;
	for (int i = 0; i < OTlist.size(); i++) {
		if (OTcheck[i] == 1) {
			final_OT[z][0] = OTlist[i][0];
			final_OT[z][1] = OTlist[i][1];
			final_OT[z][2] = OTlist[i][2];
			z++;
		}
	}
	OTfH.OTmat = final_OT;

	delete[] OTcheck;

	return OTfH;
}

#pragma endregion

int main() {

	IloEnv env;

	filebuf fb;
	fb.open("#output.txt", ios::out);
	ostream os(&fb);

	int N;
	string GRAPH;
	ifstream inst("#instances.txt");
	while (inst >> N >> GRAPH) {

		int** Adj = new int*[N];
		for (int i = 0; i < N; i++) {
			Adj[i] = new int[N];
		}
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				Adj[i][j] = 0;
			}
		}
		int numEdge = 0;
		int s, t;
		ifstream data(GRAPH);
		while (data >> s >> t) {
			Adj[s - 1][t - 1] = 1;
			Adj[t - 1][s - 1] = 1;
			numEdge++;
		}
		int m = 0;
		srand(time(NULL));

		// ================ MODEL ================

		// *** VARIABLES ***

		IloIntVarArray x(env);
		for (int i = 0; i < N; i++) {
			x.add(IloIntVar(env, 0, 1));
		}

		// *** OBJECTIVE ***

		IloModel modelIP_org(env);

		IloExpr objFuncIP(env);
		for (int i = 0; i < N; i++) {
			objFuncIP += x[i];
		}
		modelIP_org.add(IloMaximize(env, objFuncIP));
		objFuncIP.end();

		// *** CONSTRAINTS ***

		// (1) OT inequalities

		IloRangeArray constraintsIP_org(env);

		double cpu_9, cpu_10;
		int remOT, addfH;

		if (Replace_OT_fH == 0) {
			for (int i = 0; i < N - 2; i++) {
				for (int j = i + 1; j < N - 1; j++) {
					for (int k = j + 1; k < N; k++) {
						if (isOT(Adj, i, j, k)) {
							constraintsIP_org.add(IloRange(env, -IloInfinity, x[i] + x[j] + x[k], 2.0));
							m++;
						}
					}
				}
			}
		}
		else {
			cpu_9 = get_cpu_time();
			rep OTfH = replace(Adj, N);
			m = OTfH.orgOT;
			remOT = OTfH.OTnum;
			addfH = OTfH.fHnum;
			for (int i = 0; i < OTfH.fHnum; i++) {
				constraintsIP_org.add(IloRange(env, -IloInfinity, x[OTfH.fHmat[i][0]] + x[OTfH.fHmat[i][1]] + x[OTfH.fHmat[i][2]] + x[OTfH.fHmat[i][3]], 2));
			}
			for (int i = 0; i < OTfH.OTnum; i++) {
				constraintsIP_org.add(IloRange(env, -IloInfinity, x[OTfH.OTmat[i][0]] + x[OTfH.OTmat[i][1]] + x[OTfH.OTmat[i][2]], 2));
			}
			cpu_10 = get_cpu_time();
			if (OTfH.OTnum != 0) {
				for (int i = 0; i < OTfH.OTnum; i++)  delete[] OTfH.OTmat[i];
				delete[] OTfH.OTmat;
			}
			if (OTfH.fHnum != 0) {
				for (int i = 0; i < OTfH.fHnum; i++)  delete[] OTfH.fHmat[i];
				delete[] OTfH.fHmat;
			}
		}
		modelIP_org.add(constraintsIP_org);

		cout << "==============================================================" << endl;
		cout << "     IP solution methods of the Max IUC problem (RND set)      " << endl;
		cout << "==============================================================" << endl;
		cout << endl;
		cout << "Instance :              " << GRAPH << endl;
		cout << "# Vertices :            " << N << endl;
		cout << "# Edges :               " << numEdge << endl;
		cout << "# Open Triangles :      " << m << endl;
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;
		cout << endl;

		os << "==============================================================" << endl;
		os << "     IP solution methods of the Max IUC problem (RND set)      " << endl;
		os << "==============================================================" << endl;
		os << endl;
		os << "Instance :              " << GRAPH << endl;
		os << "# Vertices :            " << N << endl;
		os << "# Edges :               " << numEdge << endl;
		os << "# Open Triangles :      " << m << endl;
		os << endl;
		os << "--------------------------------------------------------------" << endl;
		os << endl;

		if (Replace_OT_fH == 1) {
			cout << " (1) Replacing OTs with 4-hole inequalities" << endl;
			cout << endl;
			cout << "Remaining # of OT inequalities :    " << remOT << endl;
			cout << "Added # of 4-Hole inequalities :    " << addfH << endl;
			cout << "Replacing OT/4-Hole CPU (sec) :     " << (cpu_10 - cpu_9) << endl;
			cout << endl;
			os << " (1) Replacing OTs with 4-hole inequalities" << endl;
			os << endl;
			os << "Remaining # of OT inequalities :    " << remOT << endl;
			os << "Added # of 4-Hole inequalities :    " << addfH << endl;
			os << "Replacing OT/4-Hole CPU (sec) :     " << (cpu_10 - cpu_9) << endl;
			os << endl;
		}

		IloModel modelLP(env);
		modelLP.add(modelIP_org);
		for (int i = 0; i < N; i++) {
			modelLP.add(IloConversion(env, x[i], IloNumVar::Float));
		}
		IloCplex cplexLP(modelLP);
		cplexLP.setOut(env.getNullStream());

		// (2) Adding four-hole inequalities

		if (addVI_fourHole == 1) {

			cplexLP.solve();

			int added_fH_cut = 0;
			IloRangeArray constraintsIP_fH(env);

			cout << " (2) 4-hole inequalities" << endl;
			cout << endl;
			cout << "LP rlx obj value (before cuts) :   " << cplexLP.getObjValue() << endl;
			os << " (2) 4-hole inequalities" << endl;
			os << endl;
			os << "LP rlx obj value (before cuts) :   " << cplexLP.getObjValue() << endl;

			double cpu_3 = get_cpu_time();
			fourHole fH = makeFourHole(Adj, N);
			for (int i = 0; i < fH.number; i++) {
				double temp = 0;
				for (int j = 0; j < 4; j++) {
					temp += cplexLP.getValue(x[fH.main[i][j]]);
				}
				if (temp > 2 + eta_fH) {
					added_fH_cut++;
					IloExpr new_cut(env);
					for (int j = 0; j < 4; j++) {
						new_cut += x[fH.main[i][j]];
					}
					constraintsIP_fH.add(IloRange(env, -IloInfinity, new_cut, 2));
					modelLP.add(new_cut <= 2);
					new_cut.end();
					cplexLP.solve();
				}
			}
			double cpu_4 = get_cpu_time();
			cout << "LP rlx obj value (after cuts) :    " << cplexLP.getObjValue() << endl;
			cout << "Total # 4-holes :                  " << fH.number << endl;
			cout << "# Added 4-hole ineq. :             " << added_fH_cut << endl;
			cout << "Generating 4-hole cuts CPU (sec) : " << (cpu_4 - cpu_3) << endl;
			cout << endl;
			os << "LP rlx obj value (after cuts) :    " << cplexLP.getObjValue() << endl;
			os << "Total # 4-holes :                  " << fH.number << endl;
			os << "# Added 4-hole ineq. :             " << added_fH_cut << endl;
			os << "Generating 4-hole cuts CPU (sec) : " << (cpu_4 - cpu_3) << endl;
			os << endl;

			modelIP_org.add(constraintsIP_fH);

			for (int i = 0; i < fH.number; i++)  delete[] fH.main[i];
			delete[] fH.main;
		}

		// (3) Adding Star inequalities --- USELESS AFTER ADDING ENOUGH NUMBER OF 4-HOLE INEQUALITIES

		if (addVI_star == 1) {

			cplexLP.solve();

			int totalStar = 0;
			int added_Star_cut = 0;
			IloRangeArray constraintsIP_star(env);

			cout << " (3) Star inequalities" << endl;
			cout << endl;
			cout << "LP rlx obj value (before cuts) :   " << cplexLP.getObjValue() << endl;
			os << " (3) Star inequalities" << endl;
			os << endl;
			os << "LP rlx obj value (before cuts) :   " << cplexLP.getObjValue() << endl;

			double cpu_5 = get_cpu_time();
			for (int v = 0; v < N; v++) {
				star vS = makeStar(Adj, N, v);
				int mini = (vS.number < 5 ? vS.number : 5);
				for (int i = 0; i < mini; i++) {
					totalStar += vS.number;
					if (vS.main[i][0] != v) cout << " *** Error! ***" << endl;
					int Icard = vS.sizes[i] - 1;
					int RHS = Icard;
					double temp = (Icard - 1)*cplexLP.getValue(x[vS.main[i][0]]);
					for (int j = 1; j < vS.sizes[i]; j++) {
						temp += cplexLP.getValue(x[vS.main[i][j]]);
					}
					if (temp > Icard + eta_S) {
						added_Star_cut++;
						IloExpr expStar(env);
						if (printVI == 1) os << " Star --- ";
						for (int j = 1; j < vS.sizes[i]; j++) {
							expStar += x[vS.main[i][j]];
							if (printVI == 1) os << "x" << vS.main[i][j] + 1 << " + ";
						}
						expStar += (Icard - 1)*x[vS.main[i][0]];
						if (printVI == 1) os << "(" << (Icard - 1) << ")x" << vS.main[i][0] + 1 << " <= " << RHS << endl;
						constraintsIP_star.add(IloRange(env, -IloInfinity, expStar, RHS));
						modelLP.add(expStar <= RHS);
						expStar.end();
						cplexLP.solve();
					}
				}
				for (int i = 0; i < vS.number; i++)  delete[] vS.main[i];
				delete[] vS.main;
				delete[] vS.sizes;
			}
			double cpu_6 = get_cpu_time();
			cout << "LP rlx obj value (after cuts) :    " << cplexLP.getObjValue() << endl;
			cout << "Total # Stars found :              " << totalStar << endl;
			cout << "# Added Star ineq. :               " << added_Star_cut << endl;
			cout << "Generating Star cuts CPU (sec) :   " << (cpu_6 - cpu_5) << endl;
			cout << endl;
			os << "LP rlx obj value (after cuts) :    " << cplexLP.getObjValue() << endl;
			os << "Total # Stars found :              " << totalStar << endl;
			os << "# Added Star ineq. :               " << added_Star_cut << endl;
			os << "Generating Star cuts CPU (sec) :   " << (cpu_6 - cpu_5) << endl;
			os << endl;

			modelIP_org.add(constraintsIP_star);
		}

		// (4) Adding Wheel & Fan inequalities

		if (addVI_wheelFan == 1) {

			IloRangeArray constraintsIP_WF(env);

			cout << " (4) Wheel and Fan inequalities" << endl;
			cout << endl;
			os << " (4) Wheel and Fan inequalities" << endl;
			os << endl;

			double cpu_7 = get_cpu_time();
			wheelFan WF = makeAllWheelFan(Adj, N);
			//checkAllWheelFan(Adj, WF);
			for (int i = 0; i < WF.numberW; i++) {
				int Hcard = WF.sizesW[i] - 1;
				int q = Hcard / 3;
				int r = Hcard % 3;
				int RHS = 2 * q + floor(2 * r / 3);
				IloExpr expWheel(env);
				if (printVI == 1) os << "Wheel --- ";
				for (int j = 1; j < WF.sizesW[i]; j++) {
					expWheel += x[WF.mainW[i][j]];
					if (printVI == 1) os << "x" << WF.mainW[i][j] + 1 << " + ";
				}
				int coeff = 2 * (q - 1) + floor(2 * r / 3);
				expWheel += coeff*x[WF.mainW[i][0]];
				if (printVI == 1) os << "(" << coeff << ")x" << WF.mainW[i][0] + 1 << " <= " << RHS << endl;
				constraintsIP_WF.add(IloRange(env, -IloInfinity, expWheel, RHS));
				expWheel.end();
			}
			for (int i = 0; i < WF.numberF; i++) {
				int Pcard = WF.sizesF[i] - 1;
				int q = Pcard / 3;
				int r = Pcard % 3;
				int RHS = 2 * q + floor(2 * (r + 1) / 3);
				IloExpr expFan(env);
				if (printVI == 1) os << "  Fan --- ";
				for (int j = 1; j < WF.sizesF[i]; j++) {
					expFan += x[WF.mainF[i][j]];
					if (printVI == 1) os << "x" << WF.mainF[i][j] + 1 << " + ";
				}
				int coeff = 2 * (q - 1) + floor(2 * (r + 1) / 3);
				expFan += coeff*x[WF.mainF[i][0]];
				if (printVI == 1) os << "(" << coeff << ")x" << WF.mainF[i][0] + 1 << " <= " << RHS << endl;
				constraintsIP_WF.add(IloRange(env, -IloInfinity, expFan, RHS));
				expFan.end();
			}
			double cpu_8 = get_cpu_time();
			cout << "min size of chordless cycle :       " << minCycleSize << endl;
			cout << "min size of chordless path :        " << minPathSize << endl;
			cout << "# Added Wheel ineq. :               " << WF.numberW << endl;
			cout << "# Added Fan ineq. :                 " << WF.numberF << endl;
			cout << "Generating W&F cuts CPU (sec) :     " << (cpu_8 - cpu_7) << endl;
			cout << endl;
			cout << "--------------------------------------------------------------" << endl;
			os << "min size of chordless cycle :       " << minCycleSize << endl;
			os << "min size of chordless path :        " << minPathSize << endl;
			os << "# Added Wheel ineq. :               " << WF.numberW << endl;
			os << "# Added Fan ineq. :                 " << WF.numberF << endl;
			os << "Generating W&F cuts CPU (sec) :     " << (cpu_8 - cpu_7) << endl;
			os << endl;
			os << "--------------------------------------------------------------" << endl;

			modelIP_org.add(constraintsIP_WF);

			for (int i = 0; i < WF.numberW; i++)  delete[] WF.mainW[i];
			delete[] WF.mainW;
			delete[] WF.sizesW;
			for (int i = 0; i < WF.numberF; i++)  delete[] WF.mainF[i];
			delete[] WF.mainF;
			delete[] WF.sizesF;
		}

		// *** SOLVE ***

		for (int i = 0; i < N; i++)  delete[] Adj[i];
		delete[] Adj;

		IloCplex cplexIP_org(modelIP_org);
		//cplexIP_org.addLazyConstraints(constraintsIP_org);
		if (printLP == 1) cplexIP_org.exportModel("IUC_VI.lp");
		if (dispProcess == 0) cplexIP_org.setOut(env.getNullStream());
		cplexIP_org.setParam(IloCplex::ClockType, 1);     //option 1: CPU time 2: wall time
		cplexIP_org.setParam(IloCplex::TiLim, timeLimit);
		if (autoCut == 0) cplexIP_org.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 1);    // Turns automatic cuts off

		double cpu_1 = get_cpu_time();
		cplexIP_org.solve();
		double cpu_2 = get_cpu_time();

		cout << "--------------------------------------------------------------" << endl;
		cout << endl;
		cout << " ***  Solution time after VIs added   *** " << endl;
		cout << endl;
		cout << "Exit flag :             " << cplexIP_org.getCplexStatus() << endl;
		cout << "Obj value :             " << cplexIP_org.getObjValue() << endl;
		cout << "CPU time (sec) :        " << (cpu_2 - cpu_1) << endl;
		cout << "#BnB nodes :            " << cplexIP_org.getNnodes() << endl;
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;

		os << "--------------------------------------------------------------" << endl;
		os << endl;
		os << " ***  Solution time after VIs added   *** " << endl;
		os << endl;
		os << "Exit flag :             " << cplexIP_org.getCplexStatus() << endl;
		os << "Obj value :             " << cplexIP_org.getObjValue() << endl;
		os << "CPU time (sec) :        " << (cpu_2 - cpu_1) << endl;
		os << "#BnB nodes :            " << cplexIP_org.getNnodes() << endl;
		os << endl;
		os << "--------------------------------------------------------------" << endl;
		//os << "Optimal Solution of IP: " << endl;
		//for (int i = 0; i < N; i++) {
		//os << "x" << i + 1 << "=" << cplexIP_org.getValue(x[i]) << ",  ";
		//}
		//os << endl;
		//os << "--------------------------------------------------------------" << endl;

		cout << "==============================================================" << endl;
		cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
		cout << "==============================================================" << endl;

		os << "==============================================================" << endl;
		os << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
		os << "==============================================================" << endl;
		//system("pause");
	}
	fb.close();
	return 0;
}

/*_________________________________________________________________________________________________________
|                                                                                                          |
|           IP solution methods for the Max IUC problem (application of HA / SD / FW inequalities)         |
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
#include <string>
#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <math.h>

using namespace std;

ILOSTLBEGIN

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

#pragma endregion

int main() {

	IloEnv env;

	filebuf fb;
	fb.open("#output.txt", ios::out);
	ostream os(&fb);

	int N;
	string GRAPH, PARTITION;
	ifstream inst("#instances.txt");
	while (inst >> N >> GRAPH >> PARTITION) {

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
		IloModel modelIP_holeAnti(env);
		IloModel modelIP_starDouble(env);
		IloModel modelIP_fanWheel(env);

		IloExpr objFuncIP(env);
		for (int i = 0; i < N; i++) {
			objFuncIP += x[i];
		}
		modelIP_org.add(IloMaximize(env, objFuncIP));
		modelIP_holeAnti.add(IloMaximize(env, objFuncIP));
		modelIP_starDouble.add(IloMaximize(env, objFuncIP));
		modelIP_fanWheel.add(IloMaximize(env, objFuncIP));
		objFuncIP.end();

		// *** CONSTRAINTS ***

		// (1) OT inequalities

		IloRangeArray constraintsIP_org(env);
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
		modelIP_org.add(constraintsIP_org);
		modelIP_holeAnti.add(constraintsIP_org);
		modelIP_starDouble.add(constraintsIP_org);
		modelIP_fanWheel.add(constraintsIP_org);

		for (int i = 0; i < N; i++)  delete[] Adj[i];
		delete[] Adj;

		cout << "==============================================================" << endl;
		cout << "         IP solution methods of the Max IUC problem           " << endl;
		cout << "==============================================================" << endl;
		cout << endl;
		cout << "Instance :              " << GRAPH << endl;
		cout << "# Vertices :            " << N << endl;
		cout << "# Edges :               " << numEdge << endl;
		cout << "# Open Triangles :      " << m << endl;
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;
		cout << endl;
		cout << " ***  Additional Constraints   *** " << endl;
		cout << endl;

		os << "==============================================================" << endl;
		os << "         IP solution methods of the Max IUC problem           " << endl;
		os << "==============================================================" << endl;
		os << endl;
		os << "Instance :              " << GRAPH << endl;
		os << "# Vertices :            " << N << endl;
		os << "# Edges :               " << numEdge << endl;
		os << "# Open Triangles :      " << m << endl;
		os << endl;
		os << "--------------------------------------------------------------" << endl;
		os << endl;
		os << " ***  Additional Constraints   *** " << endl;
		os << endl;

		// (2) Substructure inequalities

		IloRangeArray constraintsIP_holeAnti(env);
		IloRangeArray constraintsIP_starDouble(env);
		IloRangeArray constraintsIP_fanWheel(env);
		int numHole = 0;
		int numAnti = 0;
		int numStar = 0;
		int numDouble = 0;
		int numFan = 0;
		int numWheel = 0;
		char c;
		int type, sVer, eVer;
		ifstream subs(PARTITION);
		while (subs >> c >> type >> sVer >> eVer) {
			if (type == 1) {
				cout << " Hole --- ";
				os << " Hole --- ";
				numHole++;
				int card = eVer - sVer + 1;
				int q = card / 3;
				int r = card % 3;
				int RHS = 2 * q + floor(2 * r / 3);
				IloExpr expHole(env);
				cout << "x: ";
				os << "x: ";
				for (int i = sVer; i < eVer + 1; i++) {
					expHole += x[i - 1];
					cout << i << " ";
					os << i << " ";
				}
				constraintsIP_holeAnti.add(IloRange(env, -IloInfinity, expHole, RHS));
				expHole.end();
				cout << "  <=  " << RHS << endl;
				os << "  <=  " << RHS << endl;
			}
			else if (type == 2) {
				cout << " Anti --- ";
				os << " Anti --- ";
				numAnti++;
				int card = eVer - sVer + 1;
				int RHS = floor(card / 2);
				IloExpr expAnti(env);
				cout << "x: ";
				os << "x: ";
				for (int i = sVer; i < eVer + 1; i++) {
					expAnti += x[i - 1];
					cout << i << " ";
					os << i << " ";
				}
				constraintsIP_holeAnti.add(IloRange(env, -IloInfinity, expAnti, RHS));
				expAnti.end();
				cout << "  <=  " << RHS << endl;
				os << "  <=  " << RHS << endl;
			}
			else if (type == 3) {
				cout << " Star --- ";
				os << " Star --- ";
				numStar++;
				int Icard = eVer - sVer;
				int RHS = Icard;
				IloExpr expStar(env);
				cout << "x: ";
				os << "x: ";
				for (int i = sVer + 1; i < eVer + 1; i++) {
					expStar += x[i - 1];
					cout << i << " ";
					os << i << " ";
				}
				expStar += (Icard - 1)*x[sVer - 1];
				cout << "(" << (Icard - 1) << ")*" << sVer << " ";
				os << "(" << (Icard - 1) << ")*" << sVer << " ";
				constraintsIP_starDouble.add(IloRange(env, -IloInfinity, expStar, RHS));
				expStar.end();
				cout << "  <=  " << RHS << endl;
				os << "  <=  " << RHS << endl;
			}
			else if (type == 4) {
				cout << "Doubl --- ";
				os << "Doubl --- ";
				numDouble++;
				int Icard = eVer - sVer - 1;
				int RHS = Icard;
				IloExpr expDoub_1(env);
				IloExpr expDoub_2(env);
				cout << "x: ";
				os << "x: ";
				for (int i = sVer + 2; i < eVer + 1; i++) {
					expDoub_1 += x[i - 1];
					expDoub_2 += x[i - 1];
					cout << i << " ";
					os << i << " ";
				}
				expDoub_1 += (Icard - 1)*x[sVer - 1];
				expDoub_1 += x[sVer];
				cout << "(" << (Icard - 1) << ")*" << sVer << " ";
				cout << "(" << (1) << ")*" << sVer + 1 << " ";
				os << "(" << (Icard - 1) << ")*" << sVer << " ";
				os << "(" << (1) << ")*" << sVer + 1 << " ";
				expDoub_2 += x[sVer - 1];
				expDoub_2 += (Icard - 1)*x[sVer];
				constraintsIP_starDouble.add(IloRange(env, -IloInfinity, expDoub_1, RHS));
				constraintsIP_starDouble.add(IloRange(env, -IloInfinity, expDoub_2, RHS));
				expDoub_1.end();
				expDoub_2.end();
				cout << "  <=  " << RHS << endl;
				cout << "Doubl --- ";
				cout << "x: IS same as above and ";
				cout << "(" << (1) << ")*" << sVer << " ";
				cout << "(" << (Icard - 1) << ")*" << sVer + 1 << " ";
				cout << "  <=  " << RHS << endl;
				os << "  <=  " << RHS << endl;
				os << "Doubl --- ";
				os << "x: IS same as above and ";
				os << "(" << (1) << ")*" << sVer << " ";
				os << "(" << (Icard - 1) << ")*" << sVer + 1 << " ";
				os << "  <=  " << RHS << endl;
			}
			else if (type == 5) {
				cout << "  Fan --- ";
				os << "  Fan --- ";
				numFan++;
				int Pcard = eVer - sVer;
				int q = Pcard / 3;
				int r = Pcard % 3;
				int RHS = 2 * q + floor(2 * (r + 1) / 3);
				IloExpr expFan(env);
				cout << "x: ";
				os << "x: ";
				for (int i = sVer + 1; i < eVer + 1; i++) {
					expFan += x[i - 1];
					cout << i << " ";
					os << i << " ";
				}
				int coeff = 2 * (q - 1) + floor(2 * (r + 1) / 3);
				expFan += coeff*x[sVer - 1];
				cout << "(" << coeff << ")*" << sVer << " ";
				os << "(" << coeff << ")*" << sVer << " ";
				constraintsIP_fanWheel.add(IloRange(env, -IloInfinity, expFan, RHS));
				expFan.end();
				cout << "  <=  " << RHS << endl;
				os << "  <=  " << RHS << endl;
			}
			else if (type == 6) {
				cout << "Wheel --- ";
				os << "Wheel --- ";
				numWheel++;
				int Hcard = eVer - sVer;
				int q = Hcard / 3;
				int r = Hcard % 3;
				int RHS = 2 * q + floor(2 * r / 3);
				IloExpr expWheel(env);
				cout << "x: ";
				os << "x: ";
				for (int i = sVer + 1; i < eVer + 1; i++) {
					expWheel += x[i - 1];
					cout << i << " ";
					os << i << " ";
				}
				int coeff = 2 * (q - 1) + floor(2 * r / 3);
				expWheel += coeff*x[sVer - 1];
				cout << "(" << coeff << ")*" << sVer << " ";
				os << "(" << coeff << ")*" << sVer << " ";
				constraintsIP_fanWheel.add(IloRange(env, -IloInfinity, expWheel, RHS));
				expWheel.end();
				cout << "  <=  " << RHS << endl;
				os << "  <=  " << RHS << endl;
			}
		}
		modelIP_holeAnti.add(constraintsIP_holeAnti);
		modelIP_starDouble.add(constraintsIP_starDouble);
		modelIP_fanWheel.add(constraintsIP_fanWheel);

		cout << endl << endl;
		cout << "Total # additional inequalities:  " << numHole + numAnti + numStar + (2 * numDouble) + numFan + numWheel << endl;
		cout << "# Hole :                          " << numHole << endl;
		cout << "# Anti-hole :                     " << numAnti << endl;
		cout << "# Star :                          " << numStar << endl;
		cout << "# Double-star :                   " << numDouble << " x 2" << endl;
		cout << "# Fan :                           " << numFan << endl;
		cout << "# Wheel :                         " << numWheel << endl;
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;

		os << endl << endl;
		os << "Total # additional inequalities:  " << numHole + numAnti + numStar + (2 * numDouble) + numFan + numWheel << endl;
		os << "# Hole :                          " << numHole << endl;
		os << "# Anti-hole :                     " << numAnti << endl;
		os << "# Star :                          " << numStar << endl;
		os << "# Double-star :                   " << numDouble << " x 2" << endl;
		os << "# Fan :                           " << numFan << endl;
		os << "# Wheel :                         " << numWheel << endl;
		os << endl;
		os << "--------------------------------------------------------------" << endl;

		// *** SOLVE ***

		IloCplex cplexIP_org(modelIP_org);
		cplexIP_org.setOut(env.getNullStream());
		cplexIP_org.setParam(IloCplex::ClockType, 1);
		cplexIP_org.setParam(IloCplex::TiLim, 36000.000);
		cplexIP_org.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 1);    // Turns automatic cuts off

		double cpu_1 = get_cpu_time();
		cplexIP_org.solve();
		double cpu_2 = get_cpu_time();

		cout << endl;
		cout << " ***  Original Formulation   *** " << endl;
		cout << endl;
		cout << "Exit flag :             " << cplexIP_org.getCplexStatus() << endl;
		cout << "Obj value :             " << cplexIP_org.getObjValue() << endl;
		cout << "CPU time (sec) :        " << (cpu_2 - cpu_1) << endl;
		cout << "#BnB nodes :            " << cplexIP_org.getNnodes() << endl;
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;

		os << endl;
		os << " ***  Original Formulation   *** " << endl;
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

		IloCplex cplexIP_holeAnti(modelIP_holeAnti);
		cplexIP_holeAnti.setOut(env.getNullStream());
		cplexIP_holeAnti.setParam(IloCplex::ClockType, 1);
		cplexIP_holeAnti.setParam(IloCplex::TiLim, 36000.000);
		cplexIP_holeAnti.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 1);    // Turns automatic cuts off

		double cpu_3 = get_cpu_time();
		cplexIP_holeAnti.solve();
		double cpu_4 = get_cpu_time();

		cout << endl;
		cout << " ***  Hole and Antihole inequalities added   *** " << endl;
		cout << endl;
		cout << "Exit flag :             " << cplexIP_holeAnti.getCplexStatus() << endl;
		cout << "Obj value :             " << cplexIP_holeAnti.getObjValue() << endl;
		cout << "CPU time (sec) :        " << (cpu_4 - cpu_3) << endl;
		cout << "#BnB nodes :            " << cplexIP_holeAnti.getNnodes() << endl;
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;

		os << endl;
		os << " ***  Hole and Antihole inequalities added   *** " << endl;
		os << endl;
		os << "Exit flag :             " << cplexIP_holeAnti.getCplexStatus() << endl;
		os << "Obj value :             " << cplexIP_holeAnti.getObjValue() << endl;
		os << "CPU time (sec) :        " << (cpu_4 - cpu_3) << endl;
		os << "#BnB nodes :            " << cplexIP_holeAnti.getNnodes() << endl;
		os << endl;
		os << "--------------------------------------------------------------" << endl;
		//os << "Optimal Solution of IP: " << endl;
		//for (int i = 0; i < N; i++) {
		//os << "x" << i + 1 << "=" << cplexIP_holeAnti.getValue(x[i]) << ",  ";
		//}
		//os << endl;
		//os << "--------------------------------------------------------------" << endl;

		IloCplex cplexIP_starDouble(modelIP_starDouble);
		cplexIP_starDouble.setOut(env.getNullStream());
		cplexIP_starDouble.setParam(IloCplex::ClockType, 1);
		cplexIP_starDouble.setParam(IloCplex::TiLim, 36000.000);
		cplexIP_starDouble.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 1);    // Turns automatic cuts off

		double cpu_5 = get_cpu_time();
		cplexIP_starDouble.solve();
		double cpu_6 = get_cpu_time();

		cout << endl;
		cout << " ***  Star and Double-star inequalities added   *** " << endl;
		cout << endl;
		cout << "Exit flag :             " << cplexIP_starDouble.getCplexStatus() << endl;
		cout << "Obj value :             " << cplexIP_starDouble.getObjValue() << endl;
		cout << "CPU time (sec) :        " << (cpu_6 - cpu_5) << endl;
		cout << "#BnB nodes :            " << cplexIP_starDouble.getNnodes() << endl;
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;

		os << endl;
		os << " ***  Star and Double-star inequalities added   *** " << endl;
		os << endl;
		os << "Exit flag :             " << cplexIP_starDouble.getCplexStatus() << endl;
		os << "Obj value :             " << cplexIP_starDouble.getObjValue() << endl;
		os << "CPU time (sec) :        " << (cpu_6 - cpu_5) << endl;
		os << "#BnB nodes :            " << cplexIP_starDouble.getNnodes() << endl;
		os << endl;
		os << "--------------------------------------------------------------" << endl;
		//os << "Optimal Solution of IP: " << endl;
		//for (int i = 0; i < N; i++) {
		//os << "x" << i + 1 << "=" << cplexIP_starDouble.getValue(x[i]) << ",  ";
		//}
		//os << endl;
		//os << "--------------------------------------------------------------" << endl;

		IloCplex cplexIP_fanWheel(modelIP_fanWheel);
		cplexIP_fanWheel.setOut(env.getNullStream());
		cplexIP_fanWheel.setParam(IloCplex::ClockType, 1);
		cplexIP_fanWheel.setParam(IloCplex::TiLim, 36000.000);
		cplexIP_fanWheel.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 1);    // Turns automatic cuts off

		double cpu_7 = get_cpu_time();
		cplexIP_fanWheel.solve();
		double cpu_8 = get_cpu_time();

		cout << endl;
		cout << " ***  Fan and Wheel inequalities added   *** " << endl;
		cout << endl;
		cout << "Exit flag :             " << cplexIP_fanWheel.getCplexStatus() << endl;
		cout << "Obj value :             " << cplexIP_fanWheel.getObjValue() << endl;
		cout << "CPU time (sec) :        " << (cpu_8 - cpu_7) << endl;
		cout << "#BnB nodes :            " << cplexIP_fanWheel.getNnodes() << endl;
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;
		cout << "==============================================================" << endl;
		cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
		cout << "==============================================================" << endl;

		os << endl;
		os << " ***  Fan and Wheel inequalities added   *** " << endl;
		os << endl;
		os << "Exit flag :             " << cplexIP_fanWheel.getCplexStatus() << endl;
		os << "Obj value :             " << cplexIP_fanWheel.getObjValue() << endl;
		os << "CPU time (sec) :        " << (cpu_8 - cpu_7) << endl;
		os << "#BnB nodes :            " << cplexIP_fanWheel.getNnodes() << endl;
		os << endl;
		os << "--------------------------------------------------------------" << endl;
		//os << "Optimal Solution of IP: " << endl;
		//for (int i = 0; i < N; i++) {
		//os << "x" << i + 1 << "=" << cplexIP_fanWheel.getValue(x[i]) << ",  ";
		//}
		//os << endl;
		//os << "--------------------------------------------------------------" << endl;
		os << "==============================================================" << endl;
		os << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
		os << "==============================================================" << endl;
	}
	fb.close();
	return 0;
}

/*_________________________________________________________________________________________________________
|                                                                                                          |
|           IP solution methods for the Max IUC problem (combination of HA & SD & FW inequalities)         |
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
		IloModel modelIP_HASD(env);
		IloModel modelIP_HAFW(env);
		IloModel modelIP_SDFW(env);
		IloModel modelIP_ALL(env);

		IloExpr objFuncIP(env);
		for (int i = 0; i < N; i++) {
			objFuncIP += x[i];
		}
		modelIP_org.add(IloMaximize(env, objFuncIP));
		modelIP_HASD.add(IloMaximize(env, objFuncIP));
		modelIP_HAFW.add(IloMaximize(env, objFuncIP));
		modelIP_SDFW.add(IloMaximize(env, objFuncIP));
		modelIP_ALL.add(IloMaximize(env, objFuncIP));
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
		modelIP_HASD.add(constraintsIP_org);
		modelIP_HAFW.add(constraintsIP_org);
		modelIP_SDFW.add(constraintsIP_org);
		modelIP_ALL.add(constraintsIP_org);

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
				numHole++;
				int card = eVer - sVer + 1;
				int q = card / 3;
				int r = card % 3;
				int RHS = 2 * q + floor(2 * r / 3);
				IloExpr expHole(env);
				for (int i = sVer; i < eVer + 1; i++) {
					expHole += x[i - 1];
				}
				constraintsIP_holeAnti.add(IloRange(env, -IloInfinity, expHole, RHS));
				expHole.end();
			}
			else if (type == 2) {
				numAnti++;
				int card = eVer - sVer + 1;
				int RHS = floor(card / 2);
				IloExpr expAnti(env);
				for (int i = sVer; i < eVer + 1; i++) {
					expAnti += x[i - 1];
				}
				constraintsIP_holeAnti.add(IloRange(env, -IloInfinity, expAnti, RHS));
				expAnti.end();
			}
			else if (type == 3) {
				numStar++;
				int Icard = eVer - sVer;
				int RHS = Icard;
				IloExpr expStar(env);
				for (int i = sVer + 1; i < eVer + 1; i++) {
					expStar += x[i - 1];
				}
				expStar += (Icard - 1)*x[sVer - 1];
				constraintsIP_starDouble.add(IloRange(env, -IloInfinity, expStar, RHS));
				expStar.end();
			}
			else if (type == 4) {
				numDouble++;
				int Icard = eVer - sVer - 1;
				int RHS = Icard;
				IloExpr expDoub_1(env);
				IloExpr expDoub_2(env);
				for (int i = sVer + 2; i < eVer + 1; i++) {
					expDoub_1 += x[i - 1];
					expDoub_2 += x[i - 1];
				}
				expDoub_1 += (Icard - 1)*x[sVer - 1];
				expDoub_1 += x[sVer];
				expDoub_2 += x[sVer - 1];
				expDoub_2 += (Icard - 1)*x[sVer];
				constraintsIP_starDouble.add(IloRange(env, -IloInfinity, expDoub_1, RHS));
				constraintsIP_starDouble.add(IloRange(env, -IloInfinity, expDoub_2, RHS));
				expDoub_1.end();
				expDoub_2.end();
			}
			else if (type == 5) {
				numFan++;
				int Pcard = eVer - sVer;
				int q = Pcard / 3;
				int r = Pcard % 3;
				int RHS = 2 * q + floor(2 * (r + 1) / 3);
				IloExpr expFan(env);
				for (int i = sVer + 1; i < eVer + 1; i++) {
					expFan += x[i - 1];
				}
				int coeff = 2 * (q - 1) + floor(2 * (r + 1) / 3);
				expFan += coeff*x[sVer - 1];
				constraintsIP_fanWheel.add(IloRange(env, -IloInfinity, expFan, RHS));
				expFan.end();
			}
			else if (type == 6) {
				numWheel++;
				int Hcard = eVer - sVer;
				int q = Hcard / 3;
				int r = Hcard % 3;
				int RHS = 2 * q + floor(2 * r / 3);
				IloExpr expWheel(env);
				for (int i = sVer + 1; i < eVer + 1; i++) {
					expWheel += x[i - 1];
				}
				int coeff = 2 * (q - 1) + floor(2 * r / 3);
				expWheel += coeff*x[sVer - 1];
				constraintsIP_fanWheel.add(IloRange(env, -IloInfinity, expWheel, RHS));
				expWheel.end();
			}
		}
		modelIP_HASD.add(constraintsIP_holeAnti);
		modelIP_HASD.add(constraintsIP_starDouble);
		modelIP_HAFW.add(constraintsIP_holeAnti);
		modelIP_HAFW.add(constraintsIP_fanWheel);
		modelIP_SDFW.add(constraintsIP_starDouble);
		modelIP_SDFW.add(constraintsIP_fanWheel);
		modelIP_ALL.add(constraintsIP_holeAnti);
		modelIP_ALL.add(constraintsIP_starDouble);
		modelIP_ALL.add(constraintsIP_fanWheel);

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
		cout << "CPLEX TIME :            " << cplexIP_org.getCplexTime() << endl;
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

		IloCplex cplexIP_HASD(modelIP_HASD);
		cplexIP_HASD.setOut(env.getNullStream());
		cplexIP_HASD.setParam(IloCplex::ClockType, 1);
		cplexIP_HASD.setParam(IloCplex::TiLim, 36000.000);
		cplexIP_HASD.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 1);    // Turns automatic cuts off

		double cpu_3 = get_cpu_time();
		cplexIP_HASD.solve();
		double cpu_4 = get_cpu_time();

		cout << endl;
		cout << " ***  HA & SD inequalities added   *** " << endl;
		cout << endl;
		cout << "CPLEX TIME :            " << cplexIP_HASD.getCplexTime() << endl;
		cout << "Exit flag :             " << cplexIP_HASD.getCplexStatus() << endl;
		cout << "Obj value :             " << cplexIP_HASD.getObjValue() << endl;
		cout << "CPU time (sec) :        " << (cpu_4 - cpu_3) << endl;
		cout << "#BnB nodes :            " << cplexIP_HASD.getNnodes() << endl;
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;

		os << endl;
		os << " ***  HA & SD inequalities added   *** " << endl;
		os << endl;
		os << "Exit flag :             " << cplexIP_HASD.getCplexStatus() << endl;
		os << "Obj value :             " << cplexIP_HASD.getObjValue() << endl;
		os << "CPU time (sec) :        " << (cpu_4 - cpu_3) << endl;
		os << "#BnB nodes :            " << cplexIP_HASD.getNnodes() << endl;
		os << endl;
		os << "--------------------------------------------------------------" << endl;
		//os << "Optimal Solution of IP: " << endl;
		//for (int i = 0; i < N; i++) {
		//os << "x" << i + 1 << "=" << cplexIP_HASD.getValue(x[i]) << ",  ";
		//}
		//os << endl;
		//os << "--------------------------------------------------------------" << endl;

		IloCplex cplexIP_HAFW(modelIP_HAFW);
		cplexIP_HAFW.setOut(env.getNullStream());
		cplexIP_HAFW.setParam(IloCplex::ClockType, 1);
		cplexIP_HAFW.setParam(IloCplex::TiLim, 36000.000);
		cplexIP_HAFW.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 1);    // Turns automatic cuts off

		double cpu_5 = get_cpu_time();
		cplexIP_HAFW.solve();
		double cpu_6 = get_cpu_time();

		cout << endl;
		cout << " ***  HA & FW inequalities added   *** " << endl;
		cout << endl;
		cout << "CPLEX TIME :            " << cplexIP_HAFW.getCplexTime() << endl;
		cout << "Exit flag :             " << cplexIP_HAFW.getCplexStatus() << endl;
		cout << "Obj value :             " << cplexIP_HAFW.getObjValue() << endl;
		cout << "CPU time (sec) :        " << (cpu_6 - cpu_5) << endl;
		cout << "#BnB nodes :            " << cplexIP_HAFW.getNnodes() << endl;
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;

		os << endl;
		os << " ***  HA & FW inequalities added   *** " << endl;
		os << endl;
		os << "Exit flag :             " << cplexIP_HAFW.getCplexStatus() << endl;
		os << "Obj value :             " << cplexIP_HAFW.getObjValue() << endl;
		os << "CPU time (sec) :        " << (cpu_6 - cpu_5) << endl;
		os << "#BnB nodes :            " << cplexIP_HAFW.getNnodes() << endl;
		os << endl;
		os << "--------------------------------------------------------------" << endl;
		//os << "Optimal Solution of IP: " << endl;
		//for (int i = 0; i < N; i++) {
		//os << "x" << i + 1 << "=" << cplexIP_HAFW.getValue(x[i]) << ",  ";
		//}
		//os << endl;
		//os << "--------------------------------------------------------------" << endl;

		IloCplex cplexIP_SDFW(modelIP_SDFW);
		cplexIP_SDFW.setOut(env.getNullStream());
		cplexIP_SDFW.setParam(IloCplex::ClockType, 1);
		cplexIP_SDFW.setParam(IloCplex::TiLim, 36000.000);
		cplexIP_SDFW.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 1);    // Turns automatic cuts off

		double cpu_7 = get_cpu_time();
		cplexIP_SDFW.solve();
		double cpu_8 = get_cpu_time();

		cout << endl;
		cout << " ***  SD & FW inequalities added   *** " << endl;
		cout << endl;
		cout << "CPLEX TIME :            " << cplexIP_SDFW.getCplexTime() << endl;
		cout << "Exit flag :             " << cplexIP_SDFW.getCplexStatus() << endl;
		cout << "Obj value :             " << cplexIP_SDFW.getObjValue() << endl;
		cout << "CPU time (sec) :        " << (cpu_8 - cpu_7) << endl;
		cout << "#BnB nodes :            " << cplexIP_SDFW.getNnodes() << endl;
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;

		os << endl;
		os << " ***  SD & FW inequalities added   *** " << endl;
		os << endl;
		os << "Exit flag :             " << cplexIP_SDFW.getCplexStatus() << endl;
		os << "Obj value :             " << cplexIP_SDFW.getObjValue() << endl;
		os << "CPU time (sec) :        " << (cpu_8 - cpu_7) << endl;
		os << "#BnB nodes :            " << cplexIP_SDFW.getNnodes() << endl;
		os << endl;
		os << "--------------------------------------------------------------" << endl;
		//os << "Optimal Solution of IP: " << endl;
		//for (int i = 0; i < N; i++) {
		//os << "x" << i + 1 << "=" << cplexIP_SDFW.getValue(x[i]) << ",  ";
		//}
		//os << endl;
		//os << "--------------------------------------------------------------" << endl;

		IloCplex cplexIP_ALL(modelIP_ALL);
		//cplexIP_ALL.exportModel("IUC_all.lp");
		cplexIP_ALL.setOut(env.getNullStream());
		cplexIP_ALL.setParam(IloCplex::ClockType, 1);
		cplexIP_ALL.setParam(IloCplex::TiLim, 36000.000);
		cplexIP_ALL.setParam(IloCplex::Param::MIP::Limits::CutsFactor, 1);    // Turns automatic cuts off

		double cpu_9 = get_cpu_time();
		cplexIP_ALL.solve();
		double cpu_10 = get_cpu_time();

		cout << endl;
		cout << " ***  ALL inequalities added   *** " << endl;
		cout << endl;
		cout << "CPLEX TIME :            " << cplexIP_ALL.getCplexTime() << endl;
		cout << "Exit flag :             " << cplexIP_ALL.getCplexStatus() << endl;
		cout << "Obj value :             " << cplexIP_ALL.getObjValue() << endl;
		cout << "CPU time (sec) :        " << (cpu_10 - cpu_9) << endl;
		cout << "#BnB nodes :            " << cplexIP_ALL.getNnodes() << endl;
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;
		cout << "==============================================================" << endl;
		cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
		cout << "==============================================================" << endl;

		os << endl;
		os << " ***  ALL inequalities added   *** " << endl;
		os << endl;
		os << "Exit flag :             " << cplexIP_ALL.getCplexStatus() << endl;
		os << "Obj value :             " << cplexIP_ALL.getObjValue() << endl;
		os << "CPU time (sec) :        " << (cpu_10 - cpu_9) << endl;
		os << "#BnB nodes :            " << cplexIP_ALL.getNnodes() << endl;
		os << endl;
		os << "--------------------------------------------------------------" << endl;
		//os << "Optimal Solution of IP: " << endl;
		//for (int i = 0; i < N; i++) {
		//os << "x" << i + 1 << "=" << cplexIP_ALL.getValue(x[i]) << ",  ";
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

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <sstream>
#include "HomotopyContinuation.h"
#include "BethePolynomialHomotopy.h"
#include "BetheRootCache.h"
using namespace std;

void test2() {
    UnityEquation start(1, 4);
    std::vector<var_t> coefs({-20.0,0,1.0,0,1.0});
    PolynomialFunction target(coefs);
    SimpleHomotopy sh(&start, &target);
    sh.setSteps(10);
    SimpleHomotopyContinuation hc;
    
    for (int i = 0; i < start.numberOfRoots(); i++) {
        Vector startRoot = start.getRoot(i);
        Solution sol(startRoot);
        hc.solve(sh, sol);
        std::cout << sol.root() << std::endl;
//        PrintVector(root);
//        std::cout << std::endl;
    }
}

void helpInfo() {
    std::cout << "Usage:\n\tbs -L <spin chain length> -M <magnon number> [-b <value of beta> -r <index of root> -l <log level> -s]" << std::endl;
}

void outputVector(const Vector& v) {
    std::cout << '[';
    for (int i = 0; i < v.size(); i++) {
        if (i) std::cout << ' ';
        std::cout << v[i];
    }
    std::cout << ']';
}

int main(int argc, char* argv[])
{
    bool displayHelp = false;
    int L = -1;
    int M = -1;
    int r = -1;
    int nr = 30;
    int maxDepth = 10;
    LogLevel level = Error;
    elem_t beta = 2 * Pi;
    bool saveResult = false;

    for (int i = 1; i < argc && !displayHelp; i++) {
        std::string cmd = argv[i];
        if (cmd == "-L") {
            // 
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            L = ToExpression(val, -1);
        } else if (cmd == "-M") {
            // 
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            M = ToExpression(val, -1);
        } else if (cmd == "-r") {
            // 
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            r = ToExpression(val, -1);
        } else if (cmd == "-b") {
            // 
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            std::size_t pos = val.find("Pi");
            if (pos == 0) {
                beta = Pi;
            } else if (pos != std::string::npos) {
                beta = ToExpression(val.substr(0, pos), 2.0) * Pi;
            } else {
                beta = ToExpression(val, 2.0 * Pi);
            }
        } else if (cmd == "-nr") {
            // 
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            nr = ToExpression(val, nr);
        } else if (cmd == "-s") {
            saveResult = true;
        } else if (cmd == "-l") {
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            level = (LogLevel)ToExpression(val, 2);
        } else if (cmd == "-d") {
            i++;
            if (i >= argc) {
                displayHelp = true;
                break;
            }
            string val = argv[i];
            maxDepth = ToExpression(val, maxDepth);
            if (maxDepth <= 0) {
                displayHelp = true;
            }
        }
    }
    
    if (L <= 0 || M <= 0 || displayHelp) {
        helpInfo();
        return -1;
    }
    // run homotopy continuation
    BethePolynomialHomotopy hom(L, M, beta);
    BetheRootCache cache(L, M);
    BetheHomotopyContinuation hc;
    hc.setMaxNumberOfNewtonRaphonEachStep(nr);
    hc.setMaxDepth(maxDepth);
    
    std::vector<int> rootIndex;
    if (r != -1) {
        rootIndex.push_back(r);
    } else {
        for (int i = 0; i < cache.numberOfRoots(); i++) {
            rootIndex.push_back(i);
            r = i;
        }
    }
    
    SET_LOG_LEVEL(level);
    for (int i = 0; i < rootIndex.size(); i++) {
        int id = rootIndex[i];
        Vector root = cache.getRoot(id);
        std::cout << "processing root #" << id << ":\t" << green;
        outputVector(root);
        std::cout << def << std::endl;
        InvertableSolution sol(root, L);
        if (sol.indexOfSingularRoot1() >= 0) { 
            std::cout << "it's a singular root: eps=" << sol.epsilon() << ", c=" << sol.c() << std::endl;
        }
        
        if (saveResult) {
            std::string file = "BS_L" + ToString(L)+ "M" + ToString(M)+ "r" + ToString(id) + ".txt";
            hc.setTraceFile(file);
            std::cout << "file: " << file << std::endl;
        }
        hc.solve(hom, sol);
        std::cout << "root at beta=" << beta << ":\t" << yellow;
        Vector finalRoot = sol.normalRoot();
        chop(finalRoot);
        outputVector(finalRoot);
        std::cout << def << std::endl;
        std::cout << "errors:\t" << hc.errors() << std::endl << "gaps:\t" << hc.gaps() << std::endl << std::endl;
    }
    
	return 0;
}


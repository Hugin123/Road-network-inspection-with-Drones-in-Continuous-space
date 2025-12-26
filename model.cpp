//
// Created by Hugin on 2025.12.22.
//

#include "model.h"


void modelSolver(const std::string& path, const std::string& logname){
    // PROCEDURE: 建立环境
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 1);
    env.start();

    // 建立空模型
    GRBModel model = GRBModel(env);

    // 输入变量
    // 建立弧-点对照表
    std::map<std::pair<int, int>, int> EdgeIndex;
    std::map<int, std::pair<int, int>> RevEdgeIndex;
    for(int e = 0; e < E_Demand.size(); e++){
        int i = E_Demand[e].first;
        int j = E_Demand[e].second;
        EdgeIndex[{i, j}] = e;
        RevEdgeIndex[e] = {i, j};
    }
    std::map<std::pair<int, int>, int> ArcIndex;
    std::map<int, std::pair<int, int>> RevArcIndex;
    for(int a = 0; a < E_DouLDemand.size(); a++){
        int i = E_DouLDemand[a].first;
        int j = E_DouLDemand[a].second;
        ArcIndex[{i, j}] = a;
        RevArcIndex[a] = {i,j};
    }

    // PROCEDURE: 正式建立决策变量
    GRBCube x(E_DouLDemand.size(), GRBMatrix(E_DouLDemand.size(), GRBVector(k)));
    GRBCube y(E_DouLDemand.size(), GRBMatrix(mathbb_B.size(), GRBVector(k)));
    GRBVector zeta(k);
    GRBCube d(E_DouLDemand.size(), GRBMatrix(E_DouLDemand.size(), GRBVector(k)));
    GRBCube f(E_DouLDemand.size(), GRBMatrix(E_DouLDemand.size(), GRBVector(k)));
    GRBMatrix z(E_Demand.size(), GRBVector(sigma.size()));
    GRBVector xi(E_Demand.size());
    GRBVector c_x(E_DouLDemand.size());
    GRBVector c_y(E_DouLDemand.size());
    GRBCube sqrt_term1(E_DouLDemand.size(), GRBMatrix(E_DouLDemand.size(), GRBVector(k)));
    GRBCube sqrt_term2(E_DouLDemand.size(), GRBMatrix(E_DouLDemand.size(), GRBVector(k)));
    GRBCube sqrt_term3(E_DouLDemand.size(), GRBMatrix(E_DouLDemand.size(), GRBVector(k)));
    GRBCube sqrt_term4(E_DouLDemand.size(), GRBMatrix(E_DouLDemand.size(), GRBVector(k)));
    GRBCube sqrt_term5(E_DouLDemand.size(), GRBMatrix(E_DouLDemand.size(), GRBVector(k)));
    GRBCube sqrt_term6(E_DouLDemand.size(), GRBMatrix(E_DouLDemand.size(), GRBVector(k)));
    GRBCube sqrt_term7(E_DouLDemand.size(), GRBMatrix(E_DouLDemand.size(), GRBVector(k)));
    GRBTensor4d alpha(2, GRBCube(E_DouLDemand.size(),
                                 GRBMatrix(E_DouLDemand.size(), GRBVector(k))));
    GRBTensor4d beta(2, GRBCube(E_DouLDemand.size(),
                                 GRBMatrix(E_DouLDemand.size(), GRBVector(k))));


    // 为变量赋予具体类型
    for(int i = 0; i < E_DouLDemand.size(); i++){
        c_x[i] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                              "c_x_(" + std::to_string(RevArcIndex[i].first) + "," +
                              std::to_string(RevArcIndex[i].second) + ")");
        c_y[i] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                              "c_y_(" + std::to_string(RevArcIndex[i].first) + "," +
                              std::to_string(RevArcIndex[i].second) + ")");

        for(int j = 0; j < E_DouLDemand.size(); j++){
            for(int k1 = 0; k1 < k; k1++){
                x[i][j][k1] = model.addVar(0.0, 1.0, 0, GRB_BINARY,
                                           "x_(" + std::to_string(RevArcIndex[i].first) + "," +
                                           std::to_string(RevArcIndex[i].second) + ")(" +
                                           std::to_string(RevArcIndex[j].first) + "," +
                                           std::to_string(RevArcIndex[j].second) + ")" +
                                           std::to_string(k1));
                d[i][j][k1] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                                           "d_(" + std::to_string(RevArcIndex[i].first) + "," +
                                           std::to_string(RevArcIndex[i].second) + ")(" +
                                           std::to_string(RevArcIndex[j].first) + "," +
                                           std::to_string(RevArcIndex[j].second) + ")" +
                                           std::to_string(k1));
                f[i][j][k1] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                                           "f_(" + std::to_string(RevArcIndex[i].first) + "," +
                                           std::to_string(RevArcIndex[i].second) + ")(" +
                                           std::to_string(RevArcIndex[j].first) + "," +
                                           std::to_string(RevArcIndex[j].second) + ")" +
                                           std::to_string(k1));
                sqrt_term1[i][j][k1] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                                                    "sqrt-term1_(" + std::to_string(RevArcIndex[i].first)
                                                    + "," + std::to_string(RevArcIndex[i].second) + ")(" +
                                                    std::to_string(RevArcIndex[j].first) + "," +
                                                    std::to_string(RevArcIndex[j].second) + ")" +
                                                    std::to_string(k1));
                sqrt_term2[i][j][k1] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                                                    "sqrt-term2_(" + std::to_string(RevArcIndex[i].first)
                                                    + "," + std::to_string(RevArcIndex[i].second) + ")(" +
                                                    std::to_string(RevArcIndex[j].first) + "," +
                                                    std::to_string(RevArcIndex[j].second) + ")" +
                                                    std::to_string(k1));
                sqrt_term3[i][j][k1] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                                                    "sqrt-term3_(" + std::to_string(RevArcIndex[i].first)
                                                    + "," + std::to_string(RevArcIndex[i].second) + ")(" +
                                                    std::to_string(RevArcIndex[j].first) + "," +
                                                    std::to_string(RevArcIndex[j].second) + ")" +
                                                    std::to_string(k1));
                sqrt_term4[i][j][k1] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                                                    "sqrt-term4_(" + std::to_string(RevArcIndex[i].first)
                                                    + "," + std::to_string(RevArcIndex[i].second) + ")(" +
                                                    std::to_string(RevArcIndex[j].first) + "," +
                                                    std::to_string(RevArcIndex[j].second) + ")" +
                                                    std::to_string(k1));
                sqrt_term5[i][j][k1] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                                                    "sqrt-term5_(" + std::to_string(RevArcIndex[i].first)
                                                    + "," + std::to_string(RevArcIndex[i].second) + ")(" +
                                                    std::to_string(RevArcIndex[j].first) + "," +
                                                    std::to_string(RevArcIndex[j].second) + ")" +
                                                    std::to_string(k1));
                sqrt_term6[i][j][k1] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                                                    "sqrt-term6_(" + std::to_string(RevArcIndex[i].first)
                                                    + "," + std::to_string(RevArcIndex[i].second) + ")(" +
                                                    std::to_string(RevArcIndex[j].first) + "," +
                                                    std::to_string(RevArcIndex[j].second) + ")" +
                                                    std::to_string(k1));
                sqrt_term7[i][j][k1] = model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS,
                                                    "sqrt-term7_(" + std::to_string(RevArcIndex[i].first)
                                                    + "," + std::to_string(RevArcIndex[i].second) + ")(" +
                                                    std::to_string(RevArcIndex[j].first) + "," +
                                                    std::to_string(RevArcIndex[j].second) + ")" +
                                                    std::to_string(k1));
                for(int b = 0; b < 2; b++){
                    alpha[b][i][j][k1] = model.addVar(0.0, 1.0, 0, GRB_BINARY,
                                                      "alpha-" + std::to_string(b) + ",(" +
                                                      std::to_string(RevArcIndex[i].first) + "," +
                                                      std::to_string(RevArcIndex[i].second) + ")(" +
                                                      std::to_string(RevArcIndex[j].first) + "," +
                                                      std::to_string(RevArcIndex[j].second) + ")" +
                                                      std::to_string(k1));
                    beta[b][i][j][k1] = model.addVar(0.0, 1.0, 0, GRB_BINARY,
                                                      "beta-" + std::to_string(b) + ",(" +
                                                      std::to_string(RevArcIndex[i].first) + "," +
                                                      std::to_string(RevArcIndex[i].second) + ")(" +
                                                      std::to_string(RevArcIndex[j].first) + "," +
                                                      std::to_string(RevArcIndex[j].second) + ")" +
                                                      std::to_string(k1));
                }
            }
        }
        for(int j = 0; j < mathbb_B.size(); j++){
            for(int k1 = 0; k1 < k; k1++){
                y[i][j][k1] = model.addVar(0.0, 1.0, 0, GRB_BINARY,
                                           "y_(" + std::to_string(RevArcIndex[i].first) + "," +
                                           std::to_string(RevArcIndex[i].second) + ")(" +
                                           std::to_string(mathbb_B[j].first) + "," +
                                           std::to_string(mathbb_B[j].second) + ")" +
                                           std::to_string(k1));
            }
        }
    }

    for(int i = 0; i < E_Demand.size(); i++){
        xi[i] = model.addVar(0.0, 1.0, 0, GRB_BINARY, "xi_(" +
                             std::to_string(RevEdgeIndex[i].first) + "," +
                             std::to_string(RevEdgeIndex[i].second) + ")");
        for(int k1 = 0; k1 < sigma.size(); k1++) {
            z[i][k1] = model.addVar(0.0, 1.0, 0, GRB_BINARY, "z_(" +
                                                             std::to_string(RevEdgeIndex[i].first) + "," +
                                                             std::to_string(RevEdgeIndex[i].second) + ")," +
                                                             std::to_string(sigma[k1]));
        }
    }

    for(int i = 0; i < k; i++){
        zeta[i] = model.addVar(0.0, 1.0, 0, GRB_BINARY, "zeta_" + std::to_string(i));
    }

    // PROCEDURE: 建立目标函数
    GRBLinExpr obj = 0.0;
    for(int i = 0; i < E_DouLDemand.size(); i++){
        for(int j = 0; j < E_DouLDemand.size(); j++){
            for(int k1 = 0; k1 < k; k1++){
//                obj += x[i][j][k1] * 1000 * Theta_E;
                obj += d[i][j][k1] * 1;
            }
        }
    }
    for(int k1 = 0; k1 < k; k1++){
        obj += zeta[k1] * Theta_D;
    }
    model.setObjective(obj, GRB_MINIMIZE);

    // PROCEDURE: 输入约束
    GRBLinExpr c1 = 0.0;
    GRBLinExpr c2 = 0.0;
    GRBLinExpr c3 = 0.0;

    // con2.2
    for(int i = 0; i < E_DouLDemand.size(); i++){
        for (int k1 = 0; k1 < k; k1++){
            c1 = 0.0;
            c2 = 0.0;
            for(int j = 0; j < E_DouLDemand.size(); j++){
                c1 += x[i][j][k1];
                c2 += x[j][i][k1];
            }
            model.addConstr(c1 == c2, "con2.2-(" + std::to_string(RevArcIndex[i].first) + "," +
                            std::to_string(RevArcIndex[i].second) + ")," +
                            std::to_string(k1));
        }
    }

    // con2.3, con2.4, con2.5.1, con2.5.2
    for(int k1 = 0; k1 < k; k1++){
        c1 = 0;
        for(int i = 0; i < E_DouLDemand.size(); i++){
            c1 += x[ArcIndex[{1, 1}]][i][k1];
            model.addConstr(x[i][i][k1] == 0, "con2.5.1-(" + std::to_string(RevArcIndex[i].first) + "," +
                                              std::to_string(RevArcIndex[i].second) + ")," +
                                              std::to_string(k1));
            int j = ArcIndex[{RevArcIndex[i].second, RevArcIndex[i].first}];
            model.addConstr(x[i][j][k1] == 0, "con2.5.2-(" + std::to_string(RevArcIndex[i].first) + "," +
                                              std::to_string(RevArcIndex[i].second) + ")," +
                                              std::to_string(k1));
        }
        model.addConstr(c1 <= zeta[k1], "con2.3-" + std::to_string(k1));
        if (k1 > 0){
            model.addConstr(zeta[k1 - 1] <= zeta[k1], "con2.4-" + std::to_string(k1));
        }
    }

    // con2.6
    for(int i = 0; i < E_Demand.size(); i++){
        c1 = 0.0;
        int fir = RevEdgeIndex[i].first;
        int sec = RevEdgeIndex[i].second;
        for(int k1 = 0; k1 < k; k1++){
            for(int j = 0; j < E_DouLDemand.size(); j++){
                c1 += x[j][ArcIndex[{fir, sec}]][k1] + x[j][ArcIndex[{sec, fir}]][k1];
            }
        }
        model.addConstr(c1 == 1 + xi[i], "con2.6-(" + std::to_string(RevEdgeIndex[i].first) + "," +
                                         std::to_string(RevEdgeIndex[i].second) + ")");
    }

    // con2.7
    for(int i = 0; i < E_DouLDemand.size(); i++){
        for(int k1 = 0; k1 < k; k1++){
            c1 = 0.0;
            for(int j = 0; j < mathbb_B.size(); j++){
                c1 += y[i][j][k1];
            }
            c2 = 0.0;
            for(int j = 0; j < E_DouLDemand.size(); j++){
                c2 += x[j][i][k1];
            }
            model.addConstr(c1 == c2, "con2.7-(" + std::to_string(RevArcIndex[i].first) + "," +
                                                   std::to_string(RevArcIndex[i].second) + ")," +
                                                   std::to_string(k1));
        }
    }

    // con2.8
    for(int k1 = 0; k1 < k; k1++){
        model.addConstr(y[ArcIndex[{0, 0}]][2][k1] == 1, "con2.8-" + std::to_string(k1));
    }

    // con2.9
    for(int i = 1; i < E_Demand.size(); i++){
        int fir = RevEdgeIndex[i].first;
        int sec = RevEdgeIndex[i].second;
        c1 = 0.0;
        for(int k1 = 0; k1 < k; k1++){
            c1 += y[ArcIndex[{fir, sec}]][2][k1];
            c1 += y[ArcIndex[{sec, fir}]][2][k1];
        }
        c1 += xi[i];
        model.addConstr(c1 == 1, "con2.9-(" + std::to_string(RevEdgeIndex[i].first) + "," +
                                              std::to_string(RevEdgeIndex[i].second) + ")");
    }

    // con2.10, con2.11
    for(int i = 0; i < E_Demand.size(); i++){
        c1 = 0.0;
        int fir = RevEdgeIndex[i].first;
        int sec = RevEdgeIndex[i].second;
        for(int k1 = 0; k1 < k; k1++){
            c1 += y[ArcIndex[{fir, sec}]][2][k1] + y[ArcIndex[{sec, fir}]][2][k1];
        }
        model.addConstr(c1 <= 1 + xi[i] + M * (1 - z[i][0]), "con2.10-" +
                                                                std::to_string(RevEdgeIndex[i].first) + "," +
                                                                std::to_string(RevEdgeIndex[i].second) + ")");
        model.addConstr(c1 >= 1 + xi[i] - M * (1 - z[i][0]), "con2.11-" +
                                                             std::to_string(RevEdgeIndex[i].first) + "," +
                                                             std::to_string(RevEdgeIndex[i].second) + ")");
    }

    // con2.12, con2.13
    for(int i = 0; i < E_Demand.size(); i++){
        c1 = 0.0;
        int fir = RevEdgeIndex[i].first;
        int sec = RevEdgeIndex[i].second;
        for(int k1 = 0; k1 < k; k1++){
            c1 += y[ArcIndex[{fir, sec}]][0][k1] + y[ArcIndex[{fir, sec}]][1][k1];
        }
        model.addConstr(c1 <= 1 + xi[i] + M * (1 - z[i][1]), "con2.12-" +
                                                             std::to_string(RevEdgeIndex[i].first) + "," +
                                                             std::to_string(RevEdgeIndex[i].second) + ")");
        model.addConstr(c1 >= 1 + xi[i] - M * (1 - z[i][1]), "con2.13-" +
                                                             std::to_string(RevEdgeIndex[i].first) + "," +
                                                             std::to_string(RevEdgeIndex[i].second) + ")");
    }

    // con2.14, con2.15
    for(int i = 0; i < E_Demand.size(); i++){
        c1 = 0.0;
        int fir = RevEdgeIndex[i].first;
        int sec = RevEdgeIndex[i].second;
        for(int k1 = 0; k1 < k; k1++){
            c1 += y[ArcIndex[{fir, sec}]][0][k1] + y[ArcIndex[{sec, fir}]][0][k1];
        }
        model.addConstr(c1 <= 1 + xi[i] + M * (1 - z[i][2]), "con2.14-" +
                                                             std::to_string(RevEdgeIndex[i].first) + "," +
                                                             std::to_string(RevEdgeIndex[i].second) + ")");
        model.addConstr(c1 >= 1 + xi[i] - M * (1 - z[i][2]), "con2.15-" +
                                                             std::to_string(RevEdgeIndex[i].first) + "," +
                                                             std::to_string(RevEdgeIndex[i].second) + ")");
    }

    // con2.16, con2.17
    for(int i = 0; i < E_Demand.size(); i++){
        c1 = 0.0;
        int fir = RevEdgeIndex[i].first;
        int sec = RevEdgeIndex[i].second;
        for(int k1 = 0; k1 < k; k1++){
            c1 += y[ArcIndex[{fir, sec}]][1][k1] + y[ArcIndex[{sec, fir}]][1][k1];
        }
        model.addConstr(c1 <= 1 + xi[i] + M * (1 - z[i][3]), "con2.16-" +
                                                             std::to_string(RevEdgeIndex[i].first) + "," +
                                                             std::to_string(RevEdgeIndex[i].second) + ")");
        model.addConstr(c1 >= 1 + xi[i] - M * (1 - z[i][3]), "con2.17-" +
                                                             std::to_string(RevEdgeIndex[i].first) + "," +
                                                             std::to_string(RevEdgeIndex[i].second) + ")");
    }

    // con2.18, con2.19
    for(int i = 0; i < E_Demand.size(); i++){
        c1 = 0.0;
        int fir = RevEdgeIndex[i].first;
        int sec = RevEdgeIndex[i].second;
        for(int k1 = 0; k1 < k; k1++){
            c1 += y[ArcIndex[{sec, fir}]][1][k1] + y[ArcIndex[{sec, fir}]][0][k1];
        }
        model.addConstr(c1 <= 1 + xi[i] + M * (1 - z[i][4]), "con2.18-" +
                                                             std::to_string(RevEdgeIndex[i].first) + "," +
                                                             std::to_string(RevEdgeIndex[i].second) + ")");
        model.addConstr(c1 >= 1 + xi[i] - M * (1 - z[i][4]), "con2.19-" +
                                                             std::to_string(RevEdgeIndex[i].first) + "," +
                                                             std::to_string(RevEdgeIndex[i].second) + ")");
    }

    // con2.20
    for(int i = 0; i < E_Demand.size(); i++){
        c1 = 0.0;
        for(int sig = 0; sig < sigma.size(); sig++){
            c1 += z[i][sig];
        }
        model.addConstr(c1 == 1, "con2.20-" +
                                 std::to_string(RevEdgeIndex[i].first) + "," +
                                 std::to_string(RevEdgeIndex[i].second) + ")");
    }

    // con2.21, con2.22
    GRBLinExpr c4 = 0.0;
    GRBLinExpr c5 = 0.0;
    for(int i = 0; i < E_DouLDemand.size(); i++){
        for(int b = 0; b <= 1; b++){
            for(int k1 = 0; k1 < k; k1++){
                c1 = 0.0;
                c2 = 0.0;
                for(int j = 0; j < E_DouLDemand.size(); j++){
                    c1 += alpha[b][i][j][k1];
                    c2 += beta[b][j][i][k1];
                }
                model.addConstr(c1 <= y[i][b][k1], "con2.21-(" + std::to_string(RevArcIndex[i].first) +
                                                       "," + std::to_string(RevArcIndex[i].second) + ")(" +
                                                       std::to_string(mathbb_B[b].first) + "," +
                                                       std::to_string(mathbb_B[b].second) + ")" +
                                                       std::to_string(k1));
                model.addConstr(c2 <= y[i][b][k1], "con2.22-(" + std::to_string(RevArcIndex[i].first) +
                                                       "," + std::to_string(RevArcIndex[i].second) + ")(" +
                                                       std::to_string(mathbb_B[b].first) + "," +
                                                       std::to_string(mathbb_B[b].second) + ")" +
                                                       std::to_string(k1));
            }
        }
    }

    // con2.23, con2.24
    for(int i = 0; i < E_DouLDemand.size(); i++){
        for(int j = 0; j < E_DouLDemand.size(); j++){
            for(int k1 = 0; k1 < k; k1++){
                model.addQConstr(sqrt_term1[i][j][k1] >=
                                (p_x[RevArcIndex[i].second] - p_x[RevArcIndex[i].first]) *
                                (p_x[RevArcIndex[i].second] - p_x[RevArcIndex[i].first]) +
                                (p_y[RevArcIndex[i].second] - p_y[RevArcIndex[i].first]) *
                                (p_y[RevArcIndex[i].second] - p_y[RevArcIndex[i].first]) -
                                M * M * (2 - x[i][j][k1] - y[i][2][k1]),
                                "con2.23-(" + std::to_string(RevArcIndex[i].first) + "," +
                                std::to_string(RevArcIndex[i].second) + ")(" +
                                std::to_string(RevArcIndex[j].first) + "," +
                                std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
                model.addQConstr(sqrt_term1[i][j][k1] <=
                                 (p_x[RevArcIndex[i].second] - p_x[RevArcIndex[i].first]) *
                                 (p_x[RevArcIndex[i].second] - p_x[RevArcIndex[i].first]) +
                                 (p_y[RevArcIndex[i].second] - p_y[RevArcIndex[i].first]) *
                                 (p_y[RevArcIndex[i].second] - p_y[RevArcIndex[i].first]) +
                                 M * M * (2 - x[i][j][k1] - y[i][2][k1]),
                                "con2.24-(" + std::to_string(RevArcIndex[i].first) + "," +
                                std::to_string(RevArcIndex[i].second) + ")(" +
                                std::to_string(RevArcIndex[j].first) + "," +
                                std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
            }
        }
    }

    // con2.23, con2.24
    for(int i = 0; i < E_DouLDemand.size(); i++){
        for(int j = 0; j < E_DouLDemand.size(); j++){
            for(int k1 = 0; k1 < k; k1++){
                model.addQConstr(sqrt_term2[i][j][k1] >=
                                 (c_x[i] - p_x[RevArcIndex[i].first]) *
                                 (c_x[i] - p_x[RevArcIndex[i].first]) +
                                 (c_y[i] - p_y[RevArcIndex[i].first]) *
                                 (c_y[i] - p_y[RevArcIndex[i].first]) -
                                 M * M * (2 - x[i][j][k1] - y[i][1][k1]),
                                 "con2.23-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
                model.addQConstr(sqrt_term2[i][j][k1] <=
                                 (c_x[i] - p_x[RevArcIndex[i].first]) *
                                 (c_x[i] - p_x[RevArcIndex[i].first]) +
                                 (c_y[i] - p_y[RevArcIndex[i].first]) *
                                 (c_y[i] - p_y[RevArcIndex[i].first]) +
                                 M * M * (2 - x[i][j][k1] - y[i][1][k1]),
                                 "con2.24-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
            }
        }
    }

    // con2.25, con2.26
    for(int i = 0; i < E_DouLDemand.size(); i++){
        for(int j = 0; j < E_DouLDemand.size(); j++){
            for(int k1 = 0; k1 < k; k1++){
                model.addQConstr(sqrt_term3[i][j][k1] >=
                                 (c_x[i] - p_x[RevArcIndex[i].second]) *
                                 (c_x[i] - p_x[RevArcIndex[i].second]) +
                                 (c_y[i] - p_y[RevArcIndex[i].second]) *
                                 (c_y[i] - p_y[RevArcIndex[i].second]) -
                                 M * M * (2 - x[i][j][k1] - y[i][0][k1]),
                                 "con2.25-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
                model.addQConstr(sqrt_term3[i][j][k1] <=
                                 (c_x[i] - p_x[RevArcIndex[i].second]) *
                                 (c_x[i] - p_x[RevArcIndex[i].second]) +
                                 (c_y[i] - p_y[RevArcIndex[i].second]) *
                                 (c_y[i] - p_y[RevArcIndex[i].second]) +
                                 M * M * (2 - x[i][j][k1] - y[i][0][k1]),
                                 "con2.26-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
            }
        }
    }

    // con2.27, con2.28
    for(int i = 0; i < E_DouLDemand.size(); i++){
        for(int j = 0; j < E_DouLDemand.size(); j++){
            for(int k1 = 0; k1 < k; k1++){
                model.addQConstr(sqrt_term4[i][j][k1] >=
                                 (c_x[i] - p_x[RevArcIndex[j].first]) *
                                 (c_x[i] - p_x[RevArcIndex[j].first]) +
                                 (c_y[i] - p_y[RevArcIndex[j].first]) *
                                 (c_y[i] - p_y[RevArcIndex[j].first]) -
                                 M * M * (3 - x[i][j][k1] - y[i][1][k1] - (y[j][1][k1] + y[j][2][k1])),
                                 "con2.27-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
                model.addQConstr(sqrt_term4[i][j][k1] <=
                                 (c_x[i] - p_x[RevArcIndex[j].first]) *
                                 (c_x[i] - p_x[RevArcIndex[j].first]) +
                                 (c_y[i] - p_y[RevArcIndex[j].first]) *
                                 (c_y[i] - p_y[RevArcIndex[j].first]) +
                                 M * M * (3 - x[i][j][k1] - y[i][1][k1] - (y[j][1][k1] + y[j][2][k1])),
                                 "con2.28-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
            }
        }
    }

    // con2.29, con2.30
    for(int i = 0; i < E_DouLDemand.size(); i++){
        for(int j = 0; j < E_DouLDemand.size(); j++){
            for(int k1 = 0; k1 < k; k1++){
                model.addQConstr(sqrt_term5[i][j][k1] >=
                                 (c_x[i] - c_x[j]) * (c_x[i] - c_x[j]) +
                                 (c_y[i] - c_y[j]) * (c_y[i] - c_y[j]) -
                                 M * M * (3 - x[i][j][k1] - y[i][1][k1] - y[j][0][k1]),
                                 "con2.29-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
                model.addQConstr(sqrt_term5[i][j][k1] <=
                                 (c_x[i] - c_x[j]) * (c_x[i] - c_x[j]) +
                                 (c_y[i] - c_y[j]) * (c_y[i] - c_y[j]) +
                                 M * M * (3 - x[i][j][k1] - y[i][1][k1] - y[j][0][k1]),
                                 "con2.30-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
            }
        }
    }

    // con2.31, con2.32
    for(int i = 0; i < E_DouLDemand.size(); i++){
        for(int j = 0; j < E_DouLDemand.size(); j++){
            for(int k1 = 0; k1 < k; k1++){
                model.addQConstr(sqrt_term6[i][j][k1] >=
                                 (p_x[RevArcIndex[i].second] - p_x[RevArcIndex[j].first]) *
                                 (p_x[RevArcIndex[i].second] - p_x[RevArcIndex[j].first]) +
                                 (p_y[RevArcIndex[i].second] - p_y[RevArcIndex[j].first]) *
                                 (p_y[RevArcIndex[i].second] - p_y[RevArcIndex[j].first]) -
                                 M * M * (3 - x[i][j][k1] - (y[i][0][k1] + y[i][2][k1]) -
                                            (y[j][1][k1] + y[j][2][k1])),
                                 "con2.31-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
                model.addQConstr(sqrt_term6[i][j][k1] <=
                                 (p_x[RevArcIndex[i].second] - p_x[RevArcIndex[j].first]) *
                                 (p_x[RevArcIndex[i].second] - p_x[RevArcIndex[j].first]) +
                                 (p_y[RevArcIndex[i].second] - p_y[RevArcIndex[j].first]) *
                                 (p_y[RevArcIndex[i].second] - p_y[RevArcIndex[j].first]) +
                                 M * M * (3 - x[i][j][k1] - (y[i][0][k1] + y[i][2][k1]) -
                                          (y[j][1][k1] + y[j][2][k1])),
                                 "con2.32-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
            }
        }
    }

    // con2.33, con2.34
    for(int i = 0; i < E_DouLDemand.size(); i++){
        for(int j = 0; j < E_DouLDemand.size(); j++){
            for(int k1 = 0; k1 < k; k1++){
                model.addQConstr(sqrt_term7[i][j][k1] >=
                                 (p_x[RevArcIndex[i].second] - c_x[j]) *
                                 (p_x[RevArcIndex[i].second] - c_x[j]) +
                                 (p_y[RevArcIndex[i].second] - c_y[j]) *
                                 (p_y[RevArcIndex[i].second] - c_y[j]) -
                                 M * M * (3 - x[i][j][k1] - (y[i][0][k1] + y[i][2][k1]) - y[j][0][k1]),
                                 "con2.33-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
                model.addQConstr(sqrt_term7[i][j][k1] <=
                                 (p_x[RevArcIndex[i].second] - c_x[j]) *
                                 (p_x[RevArcIndex[i].second] - c_x[j]) +
                                 (p_y[RevArcIndex[i].second] - c_y[j]) *
                                 (p_y[RevArcIndex[i].second] - c_y[j]) +
                                 M * M * (3 - x[i][j][k1] - (y[i][0][k1] + y[i][2][k1]) - y[j][0][k1]),
                                 "con2.34-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
            }
        }
    }

    // con2.35
    for(int i = 0; i < E_DouLDemand.size(); i++){
        for(int j = 0; j < E_DouLDemand.size(); j++){
            for(int k1 = 0; k1 < k; k1++){
                model.addQConstr(d[i][j][k1] * d[i][j][k1] >=
                                 sqrt_term1[i][j][k1] + sqrt_term2[i][j][k1] + sqrt_term3[i][j][k1] +
                                 sqrt_term4[i][j][k1] + sqrt_term5[i][j][k1] + sqrt_term6[i][j][k1] +
                                 sqrt_term7[i][j][k1],
                                 "con2.35-(" + std::to_string(RevArcIndex[i].first) + "," +
                                 std::to_string(RevArcIndex[i].second) + ")(" +
                                 std::to_string(RevArcIndex[j].first) + "," +
                                 std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
            }
        }
    }

    // con2.36, con2.37
    for(int i = 0; i < E_DouLDemand.size(); i++){
        int tmp_i = RevArcIndex[i].first;
        int tmp_j = RevArcIndex[i].second;
        model.addConstr((p_y[tmp_j] - p_y[tmp_i]) * (c_x[i] - p_x[tmp_i]) - (p_x[tmp_j] - p_x[tmp_i]) *
                        (c_y[i] - p_y[tmp_i]) == 0,
                        "con2.36-(" + std::to_string(RevArcIndex[i].first) + "," +
                        std::to_string(RevArcIndex[i].second) + ")");
        model.addQConstr((c_x[i] - p_x[tmp_i]) * (c_x[i] - p_x[tmp_j]) +
                        (c_y[i] - p_y[tmp_i]) * (c_y[i] - p_y[tmp_j]) <= 0,
                        "con2.37-(" + std::to_string(RevArcIndex[i].first) + "," +
                        std::to_string(RevArcIndex[i].second) + ")");
    }

    // con2.38
    for(int k1 = 0; k1 < k; k1++){
        c1 = 0.0;
        for(int i = 0; i < E_DouLDemand.size(); i++){
            for(int j = 0; j < E_DouLDemand.size(); j++){
                c1 += d[i][j][k1];
            }
        }
        model.addConstr(c1 <= Q, "con2.38-" + std::to_string(k1));
    }

    // con2.39, con2.40
    for(int i = 1; i < E_DouLDemand.size(); i++){
        int tmp_i = ArcIndex[{RevArcIndex[i].second, RevArcIndex[i].first}];
        model.addConstr(c_x[i] == c_x[tmp_i], "con2.39-(" + std::to_string(RevArcIndex[i].first) +
                                              "," + std::to_string(RevArcIndex[i].second) + ")");
        model.addConstr(c_y[i] == c_y[tmp_i], "con2.40-(" + std::to_string(RevArcIndex[i].first) +
                                              "," + std::to_string(RevArcIndex[i].second) + ")");
    }

    // con2.41, con2.42
    for(int i = 1; i < E_DouLDemand.size(); i++){
        for(int j = 0; j < E_DouLDemand.size(); j++){
            for(int k1 = 0; k1 < k; k1++){
                c1 = 0.0;
                for(int i1 = 0; i1 < E_DouLDemand.size(); i1++){
                    c1 += f[i][i1][k1];
                }
                model.addConstr(c1 >= f[j][i][k1] - M * M * (1-x[j][i][k1]),
                                "con2.41-(" + std::to_string(RevArcIndex[i].first) + "," +
                                std::to_string(RevArcIndex[i].second) + ")(" +
                                std::to_string(RevArcIndex[j].first) + "," +
                                std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
                model.addConstr(c1 <= f[j][i][k1] + M * M * (1-x[j][i][k1]),
                                "con2.42-(" + std::to_string(RevArcIndex[i].first) + "," +
                                std::to_string(RevArcIndex[i].second) + ")(" +
                                std::to_string(RevArcIndex[j].first) + "," +
                                std::to_string(RevArcIndex[j].second) + ")," + std::to_string(k1));
            }
        }
    }

    // con2.43, con2.44
    for(int k1 = 0; k1 < k; k1++){
        c1 = 0.0;
        c2 = 0.0;
        c3 = 0.0;
        for(int i = 0; i < E_DouLDemand.size(); i++){
            c1 += f[0][i][k1];
            c2 += x[0][i][k1];
            c3 += d[0][i][k1];
        }
        model.addConstr(c1 >= Q - c3 - M * M * (1 - c2),
                        "con2.43-" + std::to_string(k1));
        model.addConstr(c1 <= Q - c3 + M * M * (1 - c2),
                        "con2.44-" + std::to_string(k1));
    }

    // con2.45, con2.46
    for(int k1 = 0; k1 < k; k1++){
        c1 = 0.0;
        c2 = 0.0;
        c3 = 0.0;
        for(int i = 0; i < E_DouLDemand.size(); i++){
            c1 += f[i][0][k1];
            c2 += x[i][0][k1];
            c3 += d[i][0][k1];
        }
        model.addConstr(c1 >= Q - c3 - M * M * (1 - c2),
                        "con2.45-" + std::to_string(k1));
        model.addConstr(c1 <= Q - c3 + M * M * (1 - c2),
                        "con2.46-" + std::to_string(k1));
    }

    // 设置调试算例: 随机算例\5-3-1(0).txt
    model.addConstr(x[ArcIndex[{0, 0}]][ArcIndex[{3, 1}]][0] == 1, "solution_con_1");
    model.addConstr(x[ArcIndex[{3, 1}]][ArcIndex[{1, 2}]][0] == 1, "solution_con_2");
    model.addConstr(x[ArcIndex[{1, 2}]][ArcIndex[{4, 3}]][0] == 1, "solution_con_3");
    model.addConstr(x[ArcIndex[{4, 3}]][ArcIndex[{3, 1}]][0] == 1, "solution_con_4");
    model.addConstr(x[ArcIndex[{3, 1}]][ArcIndex[{0, 0}]][0] == 1, "solution_con_6");
    model.addConstr(y[ArcIndex[{3, 1}]][0][0] == 1, "solution_con_7");
    model.addConstr(y[ArcIndex[{3, 1}]][1][0] == 1, "solution_con_8");
    model.addConstr(y[ArcIndex[{1, 2}]][2][0] == 1, "solution_con_9");
    model.addConstr(y[ArcIndex[{4, 3}]][2][0] == 1, "solution_con_10");



    // 设置求解参数
    std::cout << "模型正式运行" << std::endl;
    model.set("OutputFlag", "1");               // 1表示输出日志, 0表示不输出
    model.set("TimeLimit", "14400");            // 限制求解时间最长为14400秒
    model.set("MIPGap", "1e-4");                // 要求上下界差值为0.01%时停止求解
    model.set("FuncNonLinear", "1");            // 要求切换求解算法，使得求解混合整数非线性模型更快、精度更高
    model.set("DisplayInterval", "10");         // 设置输出间隔时长为10秒
    model.set("LogFile", path + logname.substr(0, logname.size()-4) + ".log");     //设置输出命名为logname的求解日志文件

    // 开始求解，如果模型不可行，输出冲突条件
    model.optimize();
    if (model.get(GRB_IntAttr_Status) == 3){
        model.computeIIS();
        model.write(path + logname.substr(0, logname.size()-4) + ".ilp");
    }

    // 求解结果写入文档, 分为两个文档, 一为结果文档包含求解的目标函数、求解时间等; 二为具体信息包含各个决策变量取值
    model.write(path + logname.substr(0, logname.size()-4) + ".lp");  // 输出线性/二次模型
    std::ofstream ResFile;
    std::ofstream VarFile;

    // 后续补足result需要获取的内容
    ResFile.open(path + "result.txt", std::ios::app);

    // 将当前需要的内容输入variable文件中
    VarFile.open(path + logname.substr(0, logname.size()-4) + "-VarFlie.txt");
    // 需要注意
    // 1 - 模型已加载但未求解
    // 2 - 找到了最优解
    // 3 - 模型不可行
    // 5 - 模型无界
    // 9 - 达到时间限制
    // 11 - 求解被中断
    if (model.get(GRB_IntAttr_Status) == 2 || model.get(GRB_IntAttr_Status) == 9){
        VarFile.precision(6);
//        VarFile.width(10);
        VarFile.setf(std::ios::fixed);

        VarFile << "Obj : " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
        VarFile << "Time : " << model.get(GRB_DoubleAttr_Runtime) << std::endl;
        VarFile << "Gap : " << model.get(GRB_DoubleAttr_MIPGap) << std::endl;
        VarFile << std::endl;

        for(int i = 0; i < E_DouLDemand.size(); i++){
            for(int j = 0; j < E_DouLDemand.size(); j++){
                for(int k1 = 0; k1 < k; k1++){
                    if (x[i][j][k1].get(GRB_DoubleAttr_X) > 0){
                        VarFile.width(20);
                        VarFile << x[i][j][k1].get(GRB_StringAttr_VarName) << "="
                                << x[i][j][k1].get(GRB_DoubleAttr_X) << "\n";
                    }
                }
            }
        }

        VarFile << std::endl;
        for(int i = 0; i < E_DouLDemand.size(); i++) {
            for(int j = 0; j < mathbb_B.size(); j++){
                for(int k1 = 0; k1 < k; k1++){
                    if (y[i][j][k1].get(GRB_DoubleAttr_X) > 0){
                        VarFile.width(20);
                        VarFile << y[i][j][k1].get(GRB_StringAttr_VarName) << "="
                                << y[i][j][k1].get(GRB_DoubleAttr_X) << "\n";
                    }
                }
            }
        }

        VarFile << std::endl;
        for(int i = 0; i < E_Demand.size(); i++){
            if (xi[i].get(GRB_DoubleAttr_X) > 0) {
                VarFile.width(20);
                VarFile << xi[i].get(GRB_StringAttr_VarName) << "="
                        << xi[i].get(GRB_DoubleAttr_X) << "\n";
            }
        }

        VarFile << std::endl;
        for(int i = 0; i < E_Demand.size(); i++){
            for(int k1 = 0; k1 < sigma.size(); k1++) {
                if (z[i][k1].get(GRB_DoubleAttr_X) > 0){
                    VarFile.width(20);
                    VarFile << z[i][k1].get(GRB_StringAttr_VarName) << "="
                            << z[i][k1].get(GRB_DoubleAttr_X) << "\n";
                }
            }
        }

        VarFile << std::endl;
        for(int i = 0; i < k; i++){
            if (zeta[i].get(GRB_DoubleAttr_X) > 0){
                VarFile.width(20);
                VarFile << zeta[i].get(GRB_StringAttr_VarName) << "="
                        << zeta[i].get(GRB_DoubleAttr_X) << "\n";
            }
        }

        VarFile << std::endl;
        for(int i = 0; i < E_DouLDemand.size(); i++) {
            for(int j = 0; j < E_DouLDemand.size(); j++){
                for(int k1 = 0; k1 < k; k1++){
                    if (d[i][j][k1].get(GRB_DoubleAttr_X) > 0){
                        VarFile.width(20);
                        VarFile << d[i][j][k1].get(GRB_StringAttr_VarName) << "="
                                << d[i][j][k1].get(GRB_DoubleAttr_X) << "\n";
                    }
                }
            }
        }

        VarFile << std::endl;
        for(int i = 0; i < E_DouLDemand.size(); i++){
            VarFile << c_x[i].get(GRB_StringAttr_VarName) << " = ";
            VarFile << "(" << std::to_string(p_x[RevArcIndex[i].first]) << ","
                    << std::to_string(p_y[RevArcIndex[i].first]) << ")\t\t";
            VarFile << "(" << c_x[i].get(GRB_DoubleAttr_X) << ","
                    << c_y[i].get(GRB_DoubleAttr_X) << ")\t\t";
            VarFile << "(" << std::to_string(p_x[RevArcIndex[i].second]) << ","
                    << std::to_string(p_y[RevArcIndex[i].second]) << ")";
            VarFile << std::endl;
        }

        VarFile << std::endl;
        for(int i = 0; i < E_DouLDemand.size(); i++){
            for(int j = 0; j < E_DouLDemand.size(); j++){
                for(int k1 = 0; k1 < k; k1++){
                    if (sqrt_term1[i][j][k1].get(GRB_DoubleAttr_X) > 0){
                        VarFile.width(20);
                        VarFile << sqrt_term1[i][j][k1].get(GRB_StringAttr_VarName) << "="
                                << sqrt_term1[i][j][k1].get(GRB_DoubleAttr_X) << "\n";
                    }
                    if (sqrt_term2[i][j][k1].get(GRB_DoubleAttr_X) > 0){
                        VarFile.width(20);
                        VarFile << sqrt_term2[i][j][k1].get(GRB_StringAttr_VarName) << "="
                                << sqrt_term2[i][j][k1].get(GRB_DoubleAttr_X) << "\n";
                    }
                    if (sqrt_term3[i][j][k1].get(GRB_DoubleAttr_X) > 0){
                        VarFile.width(20);
                        VarFile << sqrt_term3[i][j][k1].get(GRB_StringAttr_VarName) << "="
                                << sqrt_term3[i][j][k1].get(GRB_DoubleAttr_X) << "\n";
                    }
                    if (sqrt_term4[i][j][k1].get(GRB_DoubleAttr_X) > 0){
                        VarFile.width(20);
                        VarFile << sqrt_term4[i][j][k1].get(GRB_StringAttr_VarName) << "="
                                << sqrt_term4[i][j][k1].get(GRB_DoubleAttr_X) << "\n";
                    }
                    if (sqrt_term5[i][j][k1].get(GRB_DoubleAttr_X) > 0){
                        VarFile.width(20);
                        VarFile << sqrt_term5[i][j][k1].get(GRB_StringAttr_VarName) << "="
                                << sqrt_term5[i][j][k1].get(GRB_DoubleAttr_X) << "\n";
                    }
                    if (sqrt_term6[i][j][k1].get(GRB_DoubleAttr_X) > 0){
                        VarFile.width(20);
                        VarFile << sqrt_term6[i][j][k1].get(GRB_StringAttr_VarName) << "="
                                << sqrt_term6[i][j][k1].get(GRB_DoubleAttr_X) << "\n";
                    }
                    if (sqrt_term7[i][j][k1].get(GRB_DoubleAttr_X) > 0){
                        VarFile.width(20);
                        VarFile << sqrt_term7[i][j][k1].get(GRB_StringAttr_VarName) << "="
                                << sqrt_term7[i][j][k1].get(GRB_DoubleAttr_X) << "\n";
                    }
                }
            }
        }

        VarFile << std::endl;
        for(int i = 0; i < E_DouLDemand.size(); i++){
            for(int j = 0; j < E_DouLDemand.size(); j++){
                for(int k1 = 0; k1 < k; k1++){
                    if (f[i][j][k1].get(GRB_DoubleAttr_X) > 0) {
                        VarFile.width(20);
                        VarFile << f[i][j][k1].get(GRB_StringAttr_VarName) << "="
                                << f[i][j][k1].get(GRB_DoubleAttr_X) << "\n";
                    }
                }
            }
        }

        VarFile << std::endl;
        VarFile.close();
    }
}


//
// Created by Hugin on 2025.12.22.
//

#ifndef BISHE_INIT_DATA_H
#define BISHE_INIT_DATA_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <cmath>
#include <algorithm>


/****        算例基础参数         ****/
extern int n_B, n_pie;          // 分别表示基站点，路网端点
extern int n_nodes, n_edges;    // 分别表示节点数量，边数量
extern int k;                   // 表示无人机数量
extern double Q;                // 无人机电池容量
extern double v;                // 无人机速度
extern double Theta_E, Theta_D; // 无人机能源成本系数（元/kwh），无人机单次调用成本
extern double e_D, e_N;         // 无人机执行任务的能耗，无人机转移的能耗（kwh/m）
extern int M;                   // 足够大的常数

/****        坐标列表            ****/
extern std::vector<double> p_x;
extern std::vector<double> p_y;

/****        需求边集合          ****/
extern std::vector<std::pair<int, int>> E_Demand;       // 需求边集合
extern std::vector<std::pair<int, int>> E_DouLDemand;   // 双向弧集合
extern std::vector<std::pair<int, int>> mathbb_B;       // 需求边的访问状态集合
extern std::vector<int> sigma;                          // 辅助变量z覆盖的条件数量

void TXTDatareader(const std::string& filename);                    // 读取TXT文档参数
std::vector<double> parseCoordinateList(const std::string& str);    // 将str转换为向量

#endif //BISHE_INIT_DATA_H

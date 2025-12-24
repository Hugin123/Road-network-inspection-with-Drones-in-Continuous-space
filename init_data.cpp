//
// Created by Hugin on 2025.12.22.
//

#include "init_data.h"


/****        算例基础参数         ****/
int n_B, n_pie;          // 分别表示基站点，路网端点
int n_nodes, n_edges;    // 分别表示节点数量，边数量
int k;                   // 表示无人机数量
double Q;                // 无人机电池容量
double v;                // 无人机速度
double Theta_E, Theta_D; // 无人机能源成本系数（元/kwh），无人机单次调用成本
double e_D, e_N;         // 无人机执行任务的能耗，无人机转移的能耗（kwh/m）
int M;                   // 足够大的常数

/****        坐标列表            ****/
std::vector<double> p_x;
std::vector<double> p_y;

/****        需求边集合          ****/
std::vector<std::pair<int, int>> E_Demand;       // 需求边集合
std::vector<std::pair<int, int>> E_DouLDemand;   // 双向弧集合
std::vector<std::pair<int, int>> mathbb_B;       // 需求边的访问状态集合
std::vector<int> sigma;                          // 辅助变量z覆盖的条件数量



void TXTDatareader(const std::string& filename){
    std::ifstream file(filename);

    if (!file.is_open()){
        std::cerr << "无法打开文件: " << filename << std::endl;
        return ;
    }

    std::string line;

    // 读取前五行
    std::getline(file, line);
    n_B = std::stoi(line);
    std::getline(file, line);
    n_pie = std::stoi(line);
    std::getline(file, line);
    n_nodes = std::stoi(line);
    std::getline(file, line);
    n_edges = std::stoi(line);
    std::getline(file, line);
    k = std::stoi(line);
    std::getline(file, line);

    // 读取算例相关参数
    std::getline(file, line);
    Q = std::stod(line);
    std::getline(file, line);
    v = std::stod(line);
    std::getline(file, line);
    Theta_E = std::stod(line);
    std::getline(file, line);
    Theta_D = std::stod(line);
    std::getline(file, line);
    e_D = std::stod(line);
    std::getline(file, line);
    e_N = std::stod(line);
    std::getline(file, line);
    M = std::stoi(line);
    std::getline(file, line);

    // 读取坐标列表
    std::getline(file, line);
    p_x = parseCoordinateList(line);
    std::getline(file, line);
    p_y = parseCoordinateList(line);
    std::getline(file, line);

    // 读取需求边连线
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        // 移除可能的制表符分隔符
        size_t tabPos = line.find('\t');
        if (tabPos != std::string::npos) {
            line = line.substr(0, tabPos);
        }

        // 解析括号中的数字对
        if (line[0] == '(') {
            // 移除括号和可能的空格
            line = line.substr(1, line.find(')') - 1);

            // 分割逗号分隔的两个数字
            size_t commaPos = line.find(',');
            if (commaPos != std::string::npos) {
                int first = stoi(line.substr(0, commaPos));
                int second = stoi(line.substr(commaPos + 1));
                E_Demand.emplace_back(first, second);
                E_DouLDemand.emplace_back(first, second);
                if (first == second) continue;
                E_DouLDemand.emplace_back(second, first);
            }
        }
    }

    file.close();

    // 其余集合填充
    mathbb_B.emplace_back(0, 1);
    mathbb_B.emplace_back(1, 0);
    mathbb_B.emplace_back(1, 1);

    for(int i = 1; i <= 5; i++) {
        sigma.push_back(i);
    }
}


std::vector<double> parseCoordinateList(const std::string& str) {
    std::vector<double> coordinates;
    std::stringstream ss(str);
    std::string token;

    while (getline(ss, token, ',')) {
        // 移除空格
        token.erase(remove(token.begin(), token.end(), ' '), token.end());
        if (!token.empty()) {
            coordinates.push_back(stod(token));
        }
    }

    return coordinates;
}


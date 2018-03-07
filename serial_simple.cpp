#include <iostream>
#include <fstream>
#include <deque>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <set>

#define __memory

using namespace std;

typedef pair<vector<bool>, vector<pair<int, int> > > playarea;


struct board_fitness_compare{
    bool operator()(const pair<int, playarea>& lhs, const pair<int, playarea>& rhs){
      return lhs.first > rhs.first;
    }
} board_fitness_compare;

string make_representation(const playarea &from){
    string repre;
    repre += to_string(from.second.back().first);
    repre += to_string(from.second.back().second);
    repre += ':';
    repre += to_string(from.second.size());
    repre += ':';
    for (bool i : from.first) {
        repre += to_string(i);
    }
    return repre;
}

void moveHorse(vector<pair<int, playarea> > &next_boards, const playarea &old, const vector<pair<int, int> > &peons, const int x, const int y, const int p, const unsigned int upperLimit, const int size, set<string>& visited){
    auto old_second_back = old.second.back();
    if(old.second.size()+1 > upperLimit){
        //cout << "   reached limit" << endl;
        return;
    }
    if(old_second_back.first+x < 0 || old_second_back.first+x >= size){
        //cout << "    move out of board" << endl;
        return;
    }
    if(old_second_back.second+y < 0 || old_second_back.second+y >= size){
        //cout << "    move out of board" << endl;
        return;
    }

    int count = 0;
    for (auto &&it : old.first) {
        if(it)count++;
    }
    auto found = find(peons.begin(), peons.end(),make_pair(old_second_back.first+x,old_second_back.second+y));
    if(found != peons.end())
        count++;

    if(old.second.size()+1 + (p-count) >= upperLimit){
      //cout << "    move doesnt lead to solution" << endl;
        return;
    }

    playarea next = old;
    //cout << "    adding move " << old.second.back().first+x << " " << old.second.back().second+y;
    next.second.emplace_back(old_second_back.first+x,old_second_back.second+y);
    auto next_second_back = next.second.back();

    int best_distance = INT32_MAX;
    for (const auto &peon : peons) {
        int distance = abs(peon.first - next_second_back.first) + abs(peon.second - next_second_back.second);
        if(distance > 0 && distance < best_distance) best_distance = distance;
    }

    //cout << " peons: ";
    int taken = 0;
    for(unsigned int i = 0; i < peons.size(); ++i ){
        if(peons[i].first == next_second_back.first && peons[i].second == next_second_back.second){
            next.first[i] = true;
            taken = 1;
        }
        //cout << next.first[i];
    }

    int fitness = 8*taken-best_distance;
    //cout << " fitness of this move is " << fitness << endl;
    next_boards.emplace_back(fitness, next);
}

int getMoves(const playarea &old, const vector<pair<int, int> > &peons, deque<playarea> &space, const int p, const unsigned int upperLimit, int size, set<string>& visited){
    //cout << "  generating moves from " << old.second.back().first << " " << old.second.back().second << endl;
    vector<pair<int, playarea> > next_boards;
    moveHorse(next_boards, old, peons, 2, -1, p, upperLimit, size, visited);
    moveHorse(next_boards, old, peons, 1, -2, p, upperLimit, size, visited);

    moveHorse(next_boards, old, peons, -1,-2, p, upperLimit, size, visited);
    moveHorse(next_boards, old, peons, -2,-1, p, upperLimit, size, visited);

    moveHorse(next_boards, old, peons, -2,1, p, upperLimit, size, visited);
    moveHorse(next_boards, old, peons, -1,2, p, upperLimit, size, visited);

    moveHorse(next_boards, old, peons, 1,2, p, upperLimit, size, visited);
    moveHorse(next_boards, old, peons, 2,1, p, upperLimit, size, visited);
    sort(next_boards.begin(), next_boards.end(), board_fitness_compare);
    for (auto &next_board : next_boards) {
#ifdef _memory
      visited.insert(make_representation(next_board.second));
#endif
      space.push_back(next_board.second);
    }
    //cout << "  added: " << static_cast<int>(next_boards.size()) << " moves" << endl;
    return static_cast<int>(next_boards.size() - 1);
}

int main(int argc, char* argv[]){
    if (argc < 2 || argc > 2) {
        cerr << "Usage:" << argv[0] << " FILE " << endl;
        return 1;
    }

    fstream file;
    file.open(argv[1]);

    int K;
    unsigned int upperLimit;
    file >> K >> upperLimit;

    vector<pair<int, int> > peons;

    playarea play;

    int p = 0;

    for(int row = 0; row < K; row++) {
        for (int col = 0; col < K; col++) {
            auto c = static_cast<char>(file.get());
            while(c == '\r' || c == '\n')c = static_cast<char>(file.get());
            if(c == '3'){
              play.second.emplace_back(row,col);
            }
            if(c == '1'){
              peons.emplace_back(row, col);
              play.first.push_back(false);
              p++;
            };
        }
    }

    deque<playarea> space;
    set<string> visited;
#ifdef _memory
    visited.insert(make_representation(play));
#endif
    space.push_back(play);

    playarea tmp;

    tmp = space.front();
    getMoves(tmp, peons, space, p, upperLimit, K, visited);

    space.pop_front();

    int depth = 3;
    double iteration = ((pow(8,depth)-1)/7)-1;

    //cout << "starting BFS" << endl;

    clock_t begin = clock();
    for(int i = 0; i < iteration; i++){
        tmp = space.front();
        getMoves(tmp, peons, space, p, upperLimit, K, visited);
        space.pop_front();
    }

    //cout << "BFS finished with " << space.size() << " states generated" << endl;

    iteration = ((pow(8,upperLimit)-depth)/7)-1;
    unsigned long best = upperLimit+1;
    //cout << "starting DFS" << endl;
    while(!space.empty()){
        tmp = space.back();
        space.pop_back();
        //cout << space.size() << ":" << tmp.second.size();

        getMoves(tmp, peons, space, p, best, K, visited);

        int trues = 0;
        for (auto &&it : tmp.first) {
          if(it)trues++;
        }

        if(trues == p && tmp.second.size() < best){

            best = tmp.second.size();
            cout << "  Found solution: " << tmp.second.size()-1 << " with moves: ";
            for (auto &it : tmp.second) {
                auto found = find(peons.begin(), peons.end(), make_pair(it.first ,it.second));
                if(found != peons.end()) cout << " *";
                else cout << ' ';
                cout << "(" << it.first << "," << it.second << ")";
            }
            cout << endl;
        }
    }
    clock_t end = clock();
    space.clear();
    visited.clear();
    peons.clear();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //cout << "Best solution found with: " << best-1 << " moves" << endl;
    cout << "Calculated in: " << elapsed_secs << "s" << endl;
    return 0;
}
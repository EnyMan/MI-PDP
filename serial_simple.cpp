#include <iostream>
#include <fstream>
#include <deque>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <bitset>

using namespace std;

typedef pair<vector<bool>, vector<pair<int, int> > > playarea;

/*
struct board_fitness_compare{
    bool operator()(const pair<int, pair>& lhs, const pair<int, auto>& rhs){
      return lhs.first < rhs.first;
    }
} board_fitness_compare;
*/
void moveHorse(vector<pair<int, playarea> > &next_boards, const playarea &old, const vector<pair<int, int> > &peons, const int x, const int y, const int p, const unsigned int upperLimit, const int size){
    if(old.second.size()+1 > upperLimit){
      return;
    }
    if(old.second.back().first+x < 0 || old.second.back().first+x >= size){
      return;
    }
    if(old.second.back().second+y < 0 || old.second.back().second+y >= size){
      return;
    }
    int count = 0;
    for (auto &&it : old.first) {
      if(it)count++;
    }
    if(old.second.size()+1 + (p-count) >= upperLimit){
      return;
    }

    int fitness = 0;
    playarea next = old;
    next.second.emplace_back(old.second.back().first+x,old.second.back().second+y);

    for(unsigned int i = 0; i < peons.size(); ++i ){
      if(peons[i].first == old.second.back().first && peons[i].second == old.second.back().second){
        next.first[i] = true;
      }
    }

    next_boards.emplace_back(fitness, next);
}

int getMoves(const playarea &old, const vector<pair<int, int> > &peons, deque<playarea> &space, const int p, const unsigned int upperLimit, int size){
    vector<pair<int, playarea> > next_boards;
    moveHorse(next_boards, old, peons, 2, -1, p, upperLimit, size);
    moveHorse(next_boards, old, peons, 1, -2, p, upperLimit, size);

    moveHorse(next_boards, old, peons, -1,-2, p, upperLimit, size);
    moveHorse(next_boards, old, peons, -2,-1, p, upperLimit, size);

    moveHorse(next_boards, old, peons, -2,1, p, upperLimit, size);
    moveHorse(next_boards, old, peons, -1,2, p, upperLimit, size);

    moveHorse(next_boards, old, peons, 1,2, p, upperLimit, size);
    moveHorse(next_boards, old, peons, 2,1, p, upperLimit, size);
    //sort(next_boards.begin(), next_boards.end(), board_fitness_compare);
    for (auto &next_board : next_boards) {
      space.push_back(next_board.second);
    }
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
              p++;
            };
        }
    }



    deque<playarea> space;

    space.push_back(play);

    playarea tmp = space.front();
    getMoves(tmp, peons, space, p, upperLimit, K);

    space.pop_front();

    int depth = 3;
    double iteration = ((pow(8,depth)-1)/7)-1;

    clock_t begin = clock();
    for(int i = 0; i < iteration; i++){
        tmp = space.front();
        getMoves(tmp, peons, space, p, upperLimit, K);
        space.pop_front();
    }

    cout << space.size() << endl;

    iteration = ((pow(8,upperLimit)-depth)/7)-1;
    unsigned long best = upperLimit;
    while(!space.empty()){
    //for(int i = 0; i < iteration-1; i++){
        tmp = space.back();
        getMoves(tmp, peons, space, p, best, K);

        int trues = 0;
        for (auto &&it : tmp.first) {
          if(it)trues++;
        }

        if(trues == p && tmp.second.size() < best){
            cout << space.size() << endl;
            best = tmp.second.size();
            cout << "Found solution: " << tmp.second.size() << " with moves: " << endl;
            for (auto &it : tmp.second) {
              cout << "(" << it.first << "," << it.second << ")" << endl;
            }
        }
        for(auto it = space.end(); it != space.end()-9; --it){
          if(tmp == *it){
            space.erase(it);
          }
        }

    }
    space.clear();
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Best solution found with: " << best << " moves" << endl << "Calculated in: " << elapsed_secs << "s" << endl;
    return 0;
}
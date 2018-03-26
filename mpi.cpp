#include <iostream>
#include <fstream>
#include <deque>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <set>
#include <bitset>
#include "mpi.h"

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

void moveHorse(vector<pair<int, playarea> > &next_boards, const playarea &old, const vector<pair<int, int> > &peons, const int x, const int y, const int p, const unsigned int upperLimit, const int size){
    auto old_second_back = old.second.back();
    if(old.second.size()+1 > upperLimit){
        return;
    }
    if(old_second_back.first+x < 0 || old_second_back.first+x >= size){
        return;
    }
    if(old_second_back.second+y < 0 || old_second_back.second+y >= size){
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
        return;
    }

    playarea next = old;
    next.second.emplace_back(old_second_back.first+x,old_second_back.second+y);
    auto next_second_back = next.second.back();

    int best_distance = INT32_MAX;
    for (const auto &peon : peons) {
        int distance = abs(peon.first - next_second_back.first) + abs(peon.second - next_second_back.second);
        if(distance > 0 && distance < best_distance) best_distance = distance;
    }

    int taken = 0;
    for(unsigned int i = 0; i < peons.size(); ++i ){
        if(peons[i].first == next_second_back.first && peons[i].second == next_second_back.second){
            next.first[i] = true;
            taken = 1;
        }
    }

    int fitness = 8*taken-best_distance;
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
    sort(next_boards.begin(), next_boards.end(), board_fitness_compare);
    for (auto &next_board : next_boards) {
      space.push_back(next_board.second);
    }
    return static_cast<int>(next_boards.size() - 1);
}

int main(int argc, char* argv[]){

    int my_rank;
    int p;

    /* start up MPI */
    MPI_Init( &argc, &argv );

    /* find out process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if (my_rank != 0) {

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

      for (int row = 0; row < K; row++) {
        for (int col = 0; col < K; col++) {
          auto c = static_cast<char>(file.get());
          while (c == '\r' || c == '\n')c = static_cast<char>(file.get());
          if (c == '3') {
            play.second.emplace_back(row, col);
          }
          if (c == '1') {
            peons.emplace_back(row, col);
            play.first.push_back(false);
            p++;
          };
        }
      }

      deque<playarea> space;

      space.push_back(play);

      playarea tmp;

      tmp = space.front();
      getMoves(tmp, peons, space, p, upperLimit, K);

      space.pop_front();

      int depth = 3;
      double iteration = ((pow(8, depth) - 1) / 7) - 1;


      /* think of parallelling this */
      while(space.size() < pow(8,p)){
      //for (int i = 0; i < iteration; i++) {
        tmp = space.front();
        getMoves(tmp, peons, space, p, upperLimit, K);
        space.pop_front();
      }

      string best_solution;

      vector<bool> working;
      working.reserve(p);
      int flag;
      MPI_Status status;
      while(!space.empty()) {
        for (int i = 0; i < p; i++) {
          if(!working[i]) {
            /*
             * serialize state+best_solution
             * (x,y)(x,y)(x,y)...:upper_limit
             * send serialized state+best_solution
             */
            MPI_Send(MPI_COMM_WORLD);
          } else {
            /* when not sending tasks receive best_solution s from slaves */
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
              int msg_size;
              MPI_Get_count(&status, MPI_CHAR, &msg_size);
              char *message = new char[msg_size+1];
              MPI_Recv(&message, msg_size, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
              working[status.MPI_SOURCE] = false;
              cout << message << endl;
            }
          }
        }
      }

      for (int i = 0; i < p; i++) {
        /* send END msg */
      }
      for (int i = 0; i < p; i++) {
        /* receive last best_solutions */
      }
      cout << best_solution << endl;
    } else {
      /* == SLAVES == */
      bool work = true;
      int flag;
      MPI_Status status;
      while(work){
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
        int msg_size;
        MPI_Get_count(&status, MPI_CHAR, &msg_size);
        auto *message = new char[msg_size+1];
        MPI_Recv(&message, msg_size, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if(message == "END"){
          work = false;
        } else {
          /* DO THE ACTUAL DFS TASK IN PARALLEL */
        }
      }


      int iteration = ((pow(8, upperLimit) - depth) / 7) - 1;
      unsigned long best = upperLimit + 1;
      string best_moves;
      while (!space.empty()) {
        tmp = space.back();
        space.pop_back();

        getMoves(tmp, peons, space, p, best, K);

        int trues = 0;
        for (auto &&it : tmp.first) {
          if (it)trues++;
        }

        if (trues == p && tmp.second.size() < best) {
          best_moves.clear();
          best = tmp.second.size();
          for (auto &it : tmp.second) {
            auto found = find(peons.begin(), peons.end(), make_pair(it.first, it.second));
            if (it != *tmp.second.begin()) best_moves += ' ';
            if (found != peons.end()) best_moves += "*";
            best_moves += "(";
            best_moves += to_string(it.first);
            best_moves += ",";
            best_moves += to_string(it.second);
            best_moves += ")";
          }
        }
      }
    }
    MPI_Finalize();
    return 0;
}
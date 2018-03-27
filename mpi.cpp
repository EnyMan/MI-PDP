#include <iostream>
#include <fstream>
#include <deque>
#include <vector>
#include <cmath>
#include <ctime>
#include <sstream>
#include <algorithm>
#include <set>
#include <unistd.h>
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
    int cpus;

    /* start up MPI */
    MPI_Init( &argc, &argv );

    /* find out process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &cpus);

    if (my_rank == 0) {
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

      string init_msg = to_string(K);
      init_msg.append(":");
      for(auto peon : peons){
        init_msg.append(to_string(peon.first));
        init_msg.append(",");
        init_msg.append(to_string(peon.second));
        init_msg.append(";");
      }
      for (int i = 1; i < cpus; i++) {
        MPI_Send(init_msg.c_str(), init_msg.size(), MPI_CHAR, i, 3, MPI_COMM_WORLD);
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
      while(space.size() < pow(8,cpus-1)){
      //for (int i = 0; i < iteration; i++) {
        tmp = space.front();
        getMoves(tmp, peons, space, p, upperLimit, K);
        space.pop_front();
      }
      cout << "generated " << pow(8,cpus-1) << " states" << endl;
      string best_solution;
      int best = upperLimit;
      vector<bool> working(static_cast<unsigned long>(cpus-1), false);
      int flag;
      MPI_Status status;
      while(!space.empty()) {
        for (int i = 1; i < cpus; i++) {
          if(!working[i]) {
            tmp = space.back();
            space.pop_back();
            string state_msg = to_string(best);
            state_msg.append(":");
            for(auto move : tmp.second){
              state_msg.append(to_string(move.first));
              state_msg.append(",");
              state_msg.append(to_string(move.second));
              state_msg.append(";");
            }
            MPI_Send(state_msg.c_str(), state_msg.size(), MPI_CHAR, i, 0, MPI_COMM_WORLD);
            working[i] = true;
          } else {
            /* when not sending tasks receive best_solution s from slaves */
            MPI_Iprobe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &flag, &status);
            if (flag) {
              int msg_size;
              MPI_Get_count(&status, MPI_CHAR, &msg_size);
              char *message = new char[msg_size+1];
              MPI_Recv(message, msg_size, MPI_CHAR, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
              working[status.MPI_SOURCE] = false;
              cout << message << endl;
            }
          }
        }
      }

      for (int i = 1; i < p; i++) {
        MPI_Send("END", 3, MPI_CHAR, i, 0, MPI_COMM_WORLD);
      }
      for (int i = 1; i < p; i++) {
        /* receive last best_solutions */
      }
      cout << best_solution << endl;
    } else {
      /* == SLAVES == */
      bool work = true;
      int flag;
      int msg_size;
      MPI_Status status;

      MPI_Probe(MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_CHAR, &msg_size);

      auto message = new char[msg_size+1];
      MPI_Recv(message, msg_size, MPI_CHAR, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
      message[msg_size] = '\0';
      string str(message);

      int K;
      vector<pair<int, int> > peons;
      bool split = false;
      string str_tmp;
      pair<int, int> pair_tmp;
      for(char c : str){
        if(!split){
          if(c == ':'){
            split = true;
            K = stoi(str_tmp);
            str_tmp.clear();
          } else {
            str_tmp.push_back(c);
          }
        } else {
          if(c == ';'){
            pair_tmp.second = stoi(str_tmp);
            str_tmp.clear();
            peons.push_back(pair_tmp);
          }else if(c == ','){
            pair_tmp.first = stoi(str_tmp);
            str_tmp.clear();
          } else{
            str_tmp.push_back(c);
          }
        }
      }
      free(message);
      while(work){
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        msg_size;
        MPI_Get_count(&status, MPI_CHAR, &msg_size);
        message = new char[msg_size+1];
        MPI_Recv(message, msg_size, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        message[msg_size] = '\0';
        string msg_str(message);
        cout << msg_str << endl;
        if(msg_str == "END"){
          work = false;
          /**/
        } else {
          int best;
          playarea state;
          split = false;
          str_tmp.clear();
          for(char c : msg_str){
            if(!split){
              if(c == ':'){
                split = true;
                best = stoi(str_tmp);
                str_tmp.clear();
              } else {
                str_tmp.push_back(c);
              }
            } else {
              if(c == ';'){
                pair_tmp.second = stoi(str_tmp);
                auto is_peon = find(peons.begin(),peons.end(),pair_tmp);

                state.first.push_back(!(is_peon == peons.end()));
                state.second.push_back(pair_tmp);
                str_tmp.clear();
              }else if(c == ','){
                pair_tmp.first = stoi(str_tmp);
                str_tmp.clear();
              } else{
                str_tmp.push_back(c);
              }
            }
          }
          /* DO THE ACTUAL DFS TASK IN PARALLEL */
        }
      }


      /*int iteration = ((pow(8, upperLimit) - depth) / 7) - 1;
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
      }*/
    }
    MPI_Finalize();
    return 0;
}
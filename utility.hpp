#include <string>
#include <vector>
#include <utility>
#include <iostream>
#ifndef UTILITY_H
#define UTILITY_H

class BotAction {
public:
    BotAction(bool shooting, double rotation, double acceleration);

    bool shooting;
    double rotation;
    double acceleration;
};

class GameState {
public:
    GameState(std::string line);
    std::ostream& submit_move(std::ostream& stream, const std::vector<BotAction>& actions);

    //int iteration;
    //Parameters parameters;
    //std::vector<std::vector<Agent>> agents;
    //std::vector<GameObject> asteroids;
    //std::vector<GameObject> bullets;
};


#endif

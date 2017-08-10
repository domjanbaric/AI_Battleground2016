#include<iostream>
#include <string>
#include <fstream>
#include "utility.hpp"

#include "json.hpp"
using nlohmann::json;

// docs: https://github.com/nlohmann/json

BotAction::BotAction(bool _shooting, double _rotation, double _acceleration) :
    shooting(_shooting), rotation(_rotation), acceleration(_acceleration) {}

GameState::GameState(std::string line) {
    auto state = json::parse(line);
    // extract fields
}

std::ostream& GameState::submit_move(std::ostream& stream, const std::vector<BotAction>& actions) {
    json j = json::array();
    for (auto a: actions) {
        j.push_back({{"shooting", a.shooting}, {"rotation", a.rotation}, {"acceleration", a.acceleration}});
    }
    stream << j << std::endl;
    return stream;
}

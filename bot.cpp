#include <iostream>
#include <string>
#include"complex"
#include <vector>
#include <fstream>
#include<ctime>
#include <thread>
#include<complex>
#include "utility.hpp"
using namespace std;
#include "json.hpp"
using nlohmann::json;

/*struct defaultt
{
int width=1920, height=1080;
double maxroatation, maxspeed, bulletspeed, bulletradius;
};*/

//defaultt df;

class player {

public: double x, y;
public: double vx, vy;
public: double radius;
public: int team;
public: bool alive;
public: int timetoshoot = -1;

public:
	player() {};
	player(double _x, double _y, int _team, double _radius, double _vx, double _vy, bool _alive) {
		x = _x;
		y = _y;
		team = _team;
		radius = _radius;
		vx = _vx;
		vy = _vy;
		alive = _alive;
	};


	//GET
	double getRadius() {
		return radius;
	};
	pair<double, double> getSpeed() {
		pair<double, double> temp;
		temp.first = vx;
		temp.second = vy;
		return temp;
	};
	pair<double, double> getCord() {
		pair<double, double> temp;
		temp.first = x;
		temp.second = y;
		return temp;
	};

	int getTeam() {
		return team;
	}


	//SET
	double setRadius(double R) {
		radius = R;
	};
	double setSpeed(double _vx, double _vy) {
		vx = _vx;
		vy = _vy;
	};
	pair<double, double> setCord(pair<double, double> temp) {
		x = temp.first;
		y = temp.second;
	};
	int setTeam(int temp) {
		team = temp;
	};

	//f
	int isCollision(player &player2) {
		pair<double, double> temp = player2.getCord();
		return pow(x - temp.first, 2) + pow(y - temp.second, 2) - (radius + player2.getRadius())*(radius + player2.getRadius()) <= 0;
	};

	void collision(player *player2) {

		pair<double, double> v2;
		v2 = player2->getSpeed();

		double speed1 = sqrt(vx*vx + vy*vy);
		double speed2 = sqrt(v2.first*v2.first + v2.second*v2.second);

		double effect1 = sqrt(player2->getRadius() / radius) / 3;
		double effect2 = sqrt(radius / player2->getRadius()) / 3;

		pair<double, double> unit1;
		unit1.first = (vx / speed1)*effect1*speed2;
		unit1.second = (vy / speed1)*effect1*speed2;

		pair<double, double> unit2;
		unit2.first = (v2.first / speed2)*effect2*speed1;
		unit2.second = (v2.second / speed2)*effect2*speed1;

		vx += unit1.first;
		vy += unit1.second;

		player2->vx += unit2.first;
		player2->vy += unit2.second;
	};

	int isCollisionEdge() {
		if (x + radius >= 1920 || x - radius <= 0 || y + radius >= 1080 || y - radius<0)
			return 1;
		else
			return 0;
	};

	void collisionEdge() {
		if (x + radius >= 1920)
			vx = -vx;
		else if (x - radius<0)
			vx = -vx;
		else if (y + radius >= 1080)
			vy = -vy;
		else if (y - radius <= 0)
			vy = -vy;
	};

	void motion() {
		x += 2 * vx;
		y += 2 * vy;
	};

	void acceleration(double a) {
		double phi = atan(vy / vx);
		vx += cos(phi);
		vy += sin(phi);
	}

	void rotation(double phi) {
		vx = vx*cos(phi) + vy*sin(phi);
		vy = -vx*sin(phi) + vy*cos(phi);
	};


};


void dfs(vector<player> dfsplayer, double rotation, int dubina, int *rez)
{
	int vsize = dfsplayer.size();
	dfsplayer[0].rotation(rotation);

	if (dubina < 3)
	for (int i = 0; i < vsize; ++i)
	{
		dfsplayer[i].x += dfsplayer[i].vx;
		dfsplayer[i].y += dfsplayer[i].vy;
	}
	else
	for (int i = 0; i < vsize; ++i)
		dfsplayer[i].motion();
	int i;
	for (i = 0; i < vsize&& dfsplayer[i].team != 1; ++i)
		;
	for (; i < vsize; ++i)
	if (dfsplayer[0].isCollision(dfsplayer[i]))
	{
		*rez = 0;
		return;
	}
	for (i = 0; i < vsize&& dfsplayer[i].team != -1; ++i)
		;
	for (; i < vsize; i++)
	for (int j = i + 1; j < vsize; j++)
	if (dfsplayer[i].isCollision(dfsplayer[j]))
		dfsplayer[i].collision(&dfsplayer[j]);
	if (dfsplayer[0].isCollisionEdge()) {
		{
			*rez = 0;
			return;
		}
	}
	for (i = 0; i < vsize&& dfsplayer[i].team != -1; ++i)
		;
	for (; i < vsize&&dfsplayer[i].team == -1; ++i)
	if (dfsplayer[i].isCollisionEdge())
		dfsplayer[i].collisionEdge();
	for (; i < dfsplayer.size(); i++)
	if (dfsplayer[i].isCollisionEdge())
		dfsplayer.erase(dfsplayer.begin() + i);
	if (dubina == 6)
	{
		*rez = 1;
		return;
	}
	int f1, f2, f3;
	dfs(dfsplayer, -0.2, dubina + 1, &f1);
	dfs(dfsplayer, 0, dubina + 1, &f2);
	dfs(dfsplayer, 0.2, dubina + 1, &f3);
	*rez = f1 + f2 + f3;
}

void dfs2(vector<player> dfsplayer, double rotation, int dubina, int *rez)
{
	int vsize = dfsplayer.size();
	dfsplayer[1].rotation(rotation);

	if (dubina < 3)
	for (int i = 0; i < vsize; ++i)
	{
		dfsplayer[i].x += dfsplayer[i].vx;
		dfsplayer[i].y += dfsplayer[i].vy;
	}
	else
	for (int i = 0; i < vsize; ++i)
		dfsplayer[i].motion();
	int i;
	for (i = 0; i < vsize&& dfsplayer[i].team != 1; ++i)
		;
	for (; i < vsize; ++i)
	if (dfsplayer[1].isCollision(dfsplayer[i]))
	{
		*rez = 0;
		return;
	}
	for (i = 0; i < vsize&& dfsplayer[i].team != -1; ++i)
		;
	for (; i < vsize; i++)
	for (int j = i + 1; j < vsize; j++)
	if (dfsplayer[i].isCollision(dfsplayer[j]))
		dfsplayer[i].collision(&dfsplayer[j]);
	if (dfsplayer[1].isCollisionEdge()) {
		{
			*rez = 0;
			return;
		}
	}
	for (i = 0; i < vsize&& dfsplayer[i].team != -1; ++i)
		;
	for (; i < vsize&&dfsplayer[i].team == -1; ++i)
	if (dfsplayer[i].isCollisionEdge())
		dfsplayer[i].collisionEdge();
	for (; i < dfsplayer.size(); i++)
	if (dfsplayer[i].isCollisionEdge())
		dfsplayer.erase(dfsplayer.begin() + i);
	if (dubina == 5)
	{
		*rez = 1;
		return;
	}
	int f1, f2, f3;
	dfs2(dfsplayer, -0.2, dubina + 1, &f1);
	dfs2(dfsplayer, 0, dubina + 1, &f2);
	dfs2(dfsplayer, 0.2, dubina + 1, &f3);
	*rez = f1 + f2 + f3;
}

/*void runrun(vector<player> players,int p,int *opasnost)
{
int i,vsize=players.size();
for (i = 0;players[i].team != -1;++i)
;
*opasnost = i;
++i;
for (;i < vsize;++i)
if ((players[i].x - players[p].x)*(players[i].x - players[p].x) + (players[i].y - players[p].y)*(players[i].y - players[p].y) < (players[*opasnost].x - players[p].x)*(players[*opasnost].x - players[p].x) + (players[*opasnost].y - players[p].y)*(players[*opasnost].y - players[p].y))
*opasnost = i;
double ud = (players[*opasnost].x - players[p].x)*(players[*opasnost].x - players[p].x) + (players[*opasnost].y - players[p].y)*(players[*opasnost].y - players[p].y);
if (players[p].x < ud || players[p].y < ud || 1920 - players[p].x < ud || 1080 - players[p].x < ud)
*opasnost = -1;
return;
}*/

void toPolar(double *x, double *y, double mix, double miy)
{
	complex<double> mycomplex(*x - mix, *y - miy);
	double temp = abs(mycomplex);
	*y = atan2(-*y + miy, *x - mix);
	*x = temp;
}
void shoot(vector<player> game, bool *rj)
{
	*rj = false;

	vector<player> on, asteroid, mi;
	int i = 0;

	while (i<game.size())
	{
		if (game.at(i).team == 0)
			mi.push_back(game.at(i));
		if (game.at(i).team == 1)
			on.push_back(game.at(i));
		if (game.at(i).team == -1)
			asteroid.push_back(game.at(i));
		i++;
	}

	for (int i = 0; i < on.size(); i++)
	{
		toPolar(&on.at(i).x, &on.at(i).y, mi.at(0).x, mi.at(0).y);

	}
	double mini = 10000;
	for (int i = 0; i < asteroid.size(); i++)
	{
		toPolar(&asteroid.at(i).x, &asteroid.at(i).y, mi.at(0).x, mi.at(0).y);

	}
	if (mi.size() == 2)
	{
		toPolar(&mi.at(1).x, &mi.at(1).y, mi.at(0).x, mi.at(0).y);
		mi.at(0).y = atan2(mi.at(0).vx, mi.at(0).vy) - 3.1415 / 2;
		mi.at(0).x = 0;
		if (mi.at(0).y + 3.14 / 6 > mi.at(1).y and mi.at(0).y - 3.14 / 6 < mi.at(1).y)
		{
			*rj = false;
			return;
		}
	}
	else
	{
		mi.at(0).y = atan2(mi.at(0).vx, mi.at(0).vy) - 3.1415 / 2;
		mi.at(0).x = 0;
	}
	for (int i = 0; i < asteroid.size(); i++)
	{
		if (mi.at(0).y + 3.14 / 6 > asteroid.at(i).y and mi.at(0).y - 3.14 / 6 < asteroid.at(i).y)
		{
			double pomoc = atan2(asteroid.at(i).vx, asteroid.at(i).vy) - 3.1415 / 2 + 3.1415;
			double temp2 = abs(atan(mi.at(0).radius / asteroid.at(i).x));
			if (pomoc - mi.at(0).y<1.1*temp2 and pomoc - mi.at(0).y>-1.1*temp2)
			{
				*rj = false;
				return;
			}
			else if (pomoc - 3.1415 - mi.at(0).y<1.1*temp2 and pomoc - 3.1415 - mi.at(0).y>-1.1*temp2)
			{
				*rj = false;
				return;
			}
			else
			{
				*rj = true;
				return;
			}

			if (asteroid.at(i).x < mini)
			{
				mini = asteroid.at(i).x;
			}
		}
	}
	for (int i = 0; i < on.size(); i++)
	{
		if (mi.at(0).y + 3.14 / 6 > on.at(i).y and mi.at(0).y - 3.14 / 6 < on.at(i).y and on.at(i).x<mini)
		{
			*rj = true;
			return;
		}
	}
	return;

}
void shoot2(vector<player> game, bool *rj)
{
	*rj = false;
	vector<player> on, asteroid, mi;
	int i = 0;
	while (i<game.size())
	{
		if (game.at(i).team == 0)
			mi.push_back(game.at(i));
		if (game.at(i).team == 1)
			on.push_back(game.at(i));
		if (game.at(i).team == -1)
			asteroid.push_back(game.at(i));
		i++;
	}
	for (int i = 0; i < on.size(); i++)
	{
		toPolar(&on.at(i).x, &on.at(i).y, mi.at(1).x, mi.at(1).y);
	}
	double mini = 10000;
	for (int i = 0; i < asteroid.size(); i++)
	{
		toPolar(&asteroid.at(i).x, &asteroid.at(i).y, mi.at(1).x, mi.at(1).y);

	}
	toPolar(&mi.at(0).x, &mi.at(0).y, mi.at(1).x, mi.at(1).y);
	mi.at(1).y = atan2(mi.at(1).vx, mi.at(1).vy) - 3.1415 / 2;
	mi.at(1).x = 0;
	for (int i = 0; i < asteroid.size(); i++)
	{
		if (mi.at(1).y + 3.14 / 6 > asteroid.at(i).y and mi.at(1).y - 3.14 / 6 < asteroid.at(i).y) {
			double pomoc = atan2(asteroid.at(i).vx, asteroid.at(i).vy) - 3.1415 / 2 + 3.1415;
			double temp2 = abs(atan(mi.at(1).radius / asteroid.at(i).x));
			if (pomoc - mi.at(1).y<1.1*temp2 and pomoc - mi.at(1).y>-1.1*temp2)
			{
				*rj = false;
				return;
			}
			else if (pomoc - 3.1415 - mi.at(1).y<1.1*temp2 and pomoc - 3.1415 - mi.at(1).y>-1.1*temp2)
			{
				*rj = false;
				return;
			}
			else
			{
				*rj = true;
				return;
			}
			if (asteroid.at(i).x < mini)
			{
				mini = asteroid.at(i).x;
			}
		}
	}
	for (int i = 0; i < on.size(); i++)
	{
		if (mi.at(1).y + 3.14 / 6 > on.at(i).y and mi.at(1).y - 3.14 / 6 < on.at(i).y and on.at(i).x<mini)
		{
			*rj = true;
			return;
		}
	}
	if (mi.at(1).y + 3.14 / 6 > mi.at(0).y and mi.at(1).y - 3.14 / 6 < mi.at(0).y)
	{
		*rj = false;
		return;
	}
	return;
}

int main() {
	for (string line; getline(cin, line);){
		clock_t t;
		t = clock();
		auto state = json::parse(line);
		vector<player> players;
		player temp;
		if (state["agents"].at(0).at(0)["alive"]) {
			temp.x = state["agents"].at(0).at(0)["object"]["x"];
			temp.y = state["agents"].at(0).at(0)["object"]["y"];
			temp.vx = state["agents"].at(0).at(0)["object"]["vector"]["i"];
			temp.vy = state["agents"].at(0).at(0)["object"]["vector"]["j"];
			temp.radius = state["agents"].at(0).at(0)["object"]["radius"];
			temp.timetoshoot = state["agents"].at(0).at(0)["turnsUntilShooting"];
			temp.team = 0;
			players.push_back(temp);
		}
		if (state["agents"].at(0).size() == 2 && state["agents"].at(0).at(1)["alive"])
		{
			temp.x = state["agents"].at(0).at(1)["object"]["x"];
			temp.y = state["agents"].at(0).at(1)["object"]["y"];
			temp.vx = state["agents"].at(0).at(1)["object"]["vector"]["i"];
			temp.vy = state["agents"].at(0).at(1)["object"]["vector"]["j"];
			temp.radius = state["agents"].at(0).at(1)["object"]["radius"];
			temp.timetoshoot = state["agents"].at(0).at(1)["turnsUntilShooting"];
			temp.team = 0;
			players.push_back(temp);
		}
		if (state["agents"].at(1).at(0)["alive"]) {
			temp.x = state["agents"].at(1).at(0)["object"]["x"];
			temp.y = state["agents"].at(1).at(0)["object"]["y"];
			temp.vx = state["agents"].at(1).at(0)["object"]["vector"]["i"];
			temp.vy = state["agents"].at(1).at(0)["object"]["vector"]["j"];
			temp.radius = state["agents"].at(1).at(0)["object"]["radius"];
			temp.timetoshoot = state["agents"].at(1).at(0)["turnsUntilShooting"];
			temp.team = 1;
			players.push_back(temp);
		}
		if (state["agents"].at(1).size() == 2 && state["agents"].at(1).at(1)["alive"])
		{
			temp.x = state["agents"].at(1).at(1)["object"]["x"];
			temp.y = state["agents"].at(1).at(1)["object"]["y"];
			temp.vx = state["agents"].at(1).at(1)["object"]["vector"]["i"];
			temp.vy = state["agents"].at(1).at(1)["object"]["vector"]["j"];
			temp.radius = state["agents"].at(1).at(1)["object"]["radius"];
			temp.timetoshoot = state["agents"].at(1).at(1)["turnsUntilShooting"];
			temp.team = 1;
			players.push_back(temp);
		}
		for (int i = 0; i < state["asteroids"].size(); ++i)
		{
			temp.x = state["asteroids"].at(i)["x"];
			temp.y = state["asteroids"].at(i)["y"];
			temp.vx = state["asteroids"].at(i)["vector"]["i"];
			temp.vy = state["asteroids"].at(i)["vector"]["j"];
			temp.radius = state["asteroids"].at(i)["radius"];
			temp.timetoshoot = -1;
			temp.team = -1;
			players.push_back(temp);
		}
		for (int i = 0; i < state["bullets"].size(); ++i)
		{
			temp.x = state["bullets"].at(i)["x"];
			temp.y = state["bullets"].at(i)["y"];
			temp.vx = state["bullets"].at(i)["vector"]["i"];
			temp.vy = state["bullets"].at(i)["vector"]["j"];
			temp.radius = state["bullets"].at(i)["radius"];
			temp.timetoshoot = -1;
			temp.team = -2;
			players.push_back(temp);
		}
		t = clock() - t;
		cerr << "time1  " << ((float)t) / CLOCKS_PER_SEC << endl;
		t = clock();
		bool pucanje;
		vector<BotAction> v;
		/*int z = 1;
		if (players[z].team == 0)
		++z;
		if (players[z + 1].team == 1)
		{
		if ((players[z].x - players[0].x)*(players[z].x - players[0].x) + (players[z].y - players[0].y)*(players[z].y - players[0].y) > (players[z + 1].x - players[0].x)*(players[z + 1].x - players[0].x) + (players[z + 1].y - players[0].y)*(players[z + 1].y - players[0].y))
		++z;
		}
		player enemy = players[z], mi = players[0];

		toPolar(&enemy.x, &enemy.y, mi.x, mi.y);
		double kut = atan2(mi.vx, mi.vy) - 3.1415 / 2;
		kut = enemy.y - kut;
		cerr << "TEST " << kut << endl;*/
		int d1, d2, d3, op;
		thread t1(dfs, players, 0.2, 0, &d1);
		thread t2(dfs, players, 0, 0, &d2);
		thread t3(dfs, players, -0.2, 0, &d3);
		thread t4(shoot, players, &pucanje);
		/*int d4;
		dfs(players, kut, 0, &d4);*/
		t1.join();
		t2.join();
		t3.join();
		t4.join();
		t = clock() - t;
		cerr << "time2  " << ((float)t) / CLOCKS_PER_SEC << endl;
		t = clock();
		if (d1 == 0 && d1 == 0 && d1 == 0)
			v.push_back(BotAction(pucanje, 0, -3));
		if (d1 == d2 && d2 == d3) {
			if (players[0].x < 1920 - players[0].x && players[0].y < 1080 - players[0].y && players[0].x < players[0].y) {
				if (players[0].vy > 0)
					v.push_back(BotAction(pucanje, -0.2, 3));
				else if (players[0].vy < 0)
					v.push_back(BotAction(pucanje, 0.2, 3));
			}
			else if (players[0].x > 1920 - players[0].x && 1920 - players[0].x < 1080 - players[0].y && 1920 - players[0].x < players[0].y) {
				if (players[0].vy > 0)
					v.push_back(BotAction(pucanje, 0.2, 3));
				else if (players[0].vy < 0)
					v.push_back(BotAction(pucanje, -0.2, 3));
			}
			else if (players[0].y < 1920 - players[0].x && players[0].y < 1080 - players[0].y && players[0].y < players[0].x) {
				if (players[0].vx < 0)
					v.push_back(BotAction(pucanje, -0.2, 3));
				else if (players[0].vx > 0)
					v.push_back(BotAction(pucanje, 0.2, 3));
			}
			else {
				if (players[0].vx < 0)
					v.push_back(BotAction(pucanje, 0.2, 3));
				else if (players[0].vx > 0)
					v.push_back(BotAction(pucanje, -0.2, 3));
			}
			/*int z = 1;
			if (players[z].team == 0)
			++z;
			if (players[z + 1].team == 1)
			{
			if ((players[z].x - players[0].x)*(players[z].x - players[0].x) + (players[z].y - players[0].y)*(players[z].y - players[0].y) < (players[z + 1].x - players[0].x)*(players[z + 1].x - players[0].x) + (players[z + 1].y - players[0].y)*(players[z + 1].y - players[0].y))
			++z;
			}
			player enemy=players[z], mi=players[0];

			toPolar(&enemy.x, &enemy.y, mi.x, mi.y);
			double kut = atan2(mi.vx, mi.vy) - 3.1415 / 2;
			kut = enemy.y - kut;
			v.push_back(BotAction(pucanje, -kut, 3));*/
		}
		else if (max(d1, max(d2, d3)) == d2)
			v.push_back(BotAction(pucanje, 0, 3));
		else if (max(d1, max(d2, d3)) == d1)
			v.push_back(BotAction(pucanje, -0.2, 3));
		else if (max(d1, max(d2, d3)) == d3)
			v.push_back(BotAction(pucanje, 0.2, 3));
		if (players[1].team == 0) {
			thread t1(dfs2, players, 0.2, 0, &d1);
			thread t2(dfs2, players, 0, 0, &d2);
			thread t3(dfs2, players, -0.2, 0, &d3);
			thread t4(shoot2, players, &pucanje);
			t1.join();
			t2.join();
			t3.join();
			t4.join();
			if (d1 == 0 && d1 == 0 && d1 == 0)
				v.push_back(BotAction(pucanje, 0, -3));
			else if (d1 == d2&&d2 == d3)
			{
				if (players[0].x < 1920 - players[0].x && players[0].y < 1080 - players[0].y && players[0].x < players[0].y) {
					if (players[0].vy > 0)
						v.push_back(BotAction(pucanje, -0.2, 3));
					else if (players[0].vy < 0)
						v.push_back(BotAction(pucanje, 0.2, 3));
				}
				else if (players[0].x > 1920 - players[0].x && 1920 - players[0].x < 1080 - players[0].y && 1920 - players[0].x < players[0].y) {
					if (players[0].vy > 0)
						v.push_back(BotAction(pucanje, 0.2, 3));
					else if (players[0].vy < 0)
						v.push_back(BotAction(pucanje, -0.2, 3));
				}
				else if (players[0].y < 1920 - players[0].x && players[0].y < 1080 - players[0].y && players[0].y < players[0].x) {
					if (players[0].vx < 0)
						v.push_back(BotAction(pucanje, -0.2, 3));
					else if (players[0].vx > 0)
						v.push_back(BotAction(pucanje, 0.2, 3));
				}
				else {
					if (players[0].vx < 0)
						v.push_back(BotAction(pucanje, 0.2, 3));
					else if (players[0].vx > 0)
						v.push_back(BotAction(pucanje, -0.2, 3));
				}
				/*int z = 2;
				if (players[z + 1].team == 1)
				{
				if ((players[z].x - players[1].x)*(players[z].x - players[1].x) + (players[z].y - players[1].y)*(players[z].y - players[1].y) < (players[z + 1].x - players[1].x)*(players[z + 1].x - players[1].x) + (players[z + 1].y - players[1].y)*(players[z + 1].y - players[1].y))
				++z;
				}
				player enemy = players[z], mi = players[0];

				toPolar(&enemy.x, &enemy.y, mi.x, mi.y);
				double kut = atan2(mi.vx, mi.vy) - 3.1415 / 2;
				kut = enemy.y - kut;
				v.push_back(BotAction(pucanje, -kut, 3));*/

			}
			else if (max(d1, max(d2, d3)) == d2)
				v.push_back(BotAction(pucanje, 0, 3));
			else if (max(d1, max(d2, d3)) == d1)
				v.push_back(BotAction(pucanje, -0.2, 3));
			else if (max(d1, max(d2, d3)) == d3)
				v.push_back(BotAction(pucanje, 0, 3));
		}
		json j = json::array();
		for (auto a : v) {
			j.push_back({ { "shooting", a.shooting }, { "rotation", a.rotation }, { "acceleration", a.acceleration } });
		}
		cout << j << std::endl;
		t = clock() - t;
		cerr << "time1  " << ((float)t) / CLOCKS_PER_SEC << endl;
	}
}

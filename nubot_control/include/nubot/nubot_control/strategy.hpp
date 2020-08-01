#ifndef STRATEGY_H
#define STRATEGY_H
#include <cmath>
#include"core.hpp"
//#include "nubot/nubot_control/goaliestrategy.h"
#include "nubot/nubot_control/role_assignment.h"
#include "nubot/nubot_control/world_model_info.h"
#include "nubot/nubot_control/goaliestrategy.h"
#include "nubot/nubot_control/plan.h"
#include "nubot/nubot_control/activerole.h"
#include "nubot/nubot_control/midfieldrole.h"
#include "nubot/nubot_control/passiverole.h"
#include "nubot/nubot_control/assistrole.h"
namespace nubot{

class Strategy
{

public:
    Strategy();
    Strategy(World_Model_Info & _world_model, Plan & _plan);
    ~Strategy();
   void selectRole(DPoint ball_pos_);
   void selectAction();
   void process();
   bool passStrategy();
   void fuzzy_evaluation();
   //守门员位置1
   void Market_act();
   void get_line(double*para, DPoint2d pose, DPoint pos_ball);
public:
    RoleAssignment    RoleAssignment_;
    World_Model_Info * world_model_;
    int selected_role_;
    int selected_action_;
    Plan * m_plan_;
    ActiveRole        ActiveRole_;
    AssistRole        AssistRole_;
    PassiveRole       PassiveRole_;
    MidfieldRole      MidfieldRole_;
    GoalieStrategy    goalie_strategy_;
    bool  auto_competition;

    //群体策略
    double weight[6] = {0.4112, 0.1761, 0.1761, 0.1108, 0.0518, 0.0839};
    double pose_us[5][3], pose_oppo[5][2], score[2], price[5][2];
    DPoint position_ball, Obstacle[8];
    int status;


};

}
#endif // STRATEGY_H

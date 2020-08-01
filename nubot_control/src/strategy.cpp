#include "nubot/nubot_control/strategy.hpp"
#include "algorithm"

using namespace nubot;
Strategy::Strategy(World_Model_Info & _world_model, Plan & _plan){
    world_model_ = & _world_model;
    m_plan_ = & _plan;
    RoleAssignment_.world_model_ = world_model_;

    ActiveRole_.m_plan_ = m_plan_;
    ActiveRole_.world_model_ = world_model_;

    AssistRole_.world_model_ = world_model_;
    AssistRole_.plan_ = m_plan_;

    PassiveRole_.world_model_ = world_model_;
    PassiveRole_.plan_ = m_plan_;

    MidfieldRole_.world_model_ = world_model_;
    MidfieldRole_.plan_ = m_plan_;
    selected_role_ = NOROLE;
    selected_action_ = Positioned;

}
Strategy::Strategy(){
    RoleAssignment_.world_model_ = world_model_;
    selected_role_ = NOROLE;
}
Strategy::~Strategy(){
}

//unable to update score & pose of oppos
//unkown ID of goalkeeper, so as position
void  Strategy::selectRole(DPoint ball_pos_)
    {
        for(std::size_t i = 0 ; i < OUR_TEAM ; i++)
        {
            pose_us[i][0] = world_model_->RobotInfo_[i].getLocation().x_;
            pose_us[i][1] = world_model_->RobotInfo_[i].getLocation().y_;
            pose_us[i][2] = world_model_->RobotInfo_[i].getHead().radian_;
        }
        for(std::size_t i = 0 ; i < world_model_->Opponents_.size() ; i++)
        {
            pose_oppo[i][0] = world_model_->Opponents_[i].x_;
            pose_oppo[i][1] = world_model_->Opponents_[i].y_;
        }
        position_ball = ball_pos_;
    }
void Strategy::fuzzy_evaluation()
    {
        double B1[3],B2[3],B3[3],B4[3],B5[3],B6[3],B7[3],B8[3];
        //B1 球的位置
        //球在Ａ１范围
        if(position_ball.x_>=-1100 && position_ball.x_<=-550 && position_ball.y_>=-500 && position_ball.y_<=500)
        {B1[0] = 0.9; B1[1]=0.1; B1[2]=0.0;}   
        //球在Ａ3范围
        else if(position_ball.x_>=0 && position_ball.x_<550 && position_ball.y_>=-700 && position_ball.y_<=700)
        {B1[0]=0.2; B1[1]=0.3; B1[2]=0.5;}
        //球在Ａ4范围
        else if(position_ball.x_>=550 && position_ball.x_<=1100 && position_ball.y_>=-700 && position_ball.y_<=700)
        {B1[0]=0.0; B1[1]=0.1; B1[2]=0.9;}        
        //球在Ａ2范围
        else if((position_ball.x_>=-1100 && position_ball.x_<0 && position_ball.y_>500 && position_ball.y_<=700) || (position_ball.x_>=-1100 && position_ball.x_<0 && position_ball.y_>=-700 && position_ball.y_<-500) || (position_ball.x_>-550 && position_ball.x_<0 && position_ball.y_>=-500 && position_ball.y_<=500))
        {B1[0]=0.5; B1[1]=0.3; B1[2]=0.2;}
        //球出界
        else   {B1[0]=-1; B1[1]=-1; B1[2]=-1;}
        
        //B2 己方x轴
        //unkown ID of goalkeeper, so as position
        //mean_x = sum(pose_us(:,1))/5.0;  gap_x = 11.0 - mean_x;
        double gap_x = 1100, L=2200;
        for(std::size_t i = 0 ; i < OUR_TEAM ; i++)
            gap_x+=pose_us[i][0]/5.0;

        if(gap_x>=0 && gap_x<L/2)
            B2[0] = 0.0;
        else if(L/2<=gap_x && gap_x<2.0*L/3.0)
            B2[0] = (gap_x - 1100)/(L/6.0);
        else if(2.0*L/3.0<=gap_x && gap_x<=L)
            B2[0] = 1.0;
        else
            std::cout<<"gap_x is error";
        
        if(0.0<=gap_x && gap_x<L/3.0)
            B2[1] = 0.0;
        else if(L/3.0<=gap_x && gap_x<L/2.0)
            B2[1] = 1 - (L/2.0 - gap_x)/(L/6.0);
        else if(L/2.0<=gap_x && gap_x<2.0*L/3.0)
            B2[1] = 1 - (gap_x - L/2.0)/(L/6.0);
        else if(2.0*L/3.0<=gap_x && gap_x<=L)
            B2[1] = 0.0;
        else
            std::cout<<"gap_x is error";
    
        if(0<=gap_x && gap_x<L/3.0)
            B2[2] = 1.0;
        else if(L/3.0<=gap_x && gap_x<L/2.0)
            B2[2] = (L/2.0 - gap_x)/(L/6.0);
        else if(L/2.0<=gap_x && gap_x <= L)
            B2[2] = 0.0;
        else 
            std::cout<<"gap_x is error";
        
        //B3 敌方x轴
        gap_x = 1100, L=2200;
        for(std::size_t i = 0 ; i < world_model_->Opponents_.size() ; i++)
            gap_x+=pose_oppo[i][0]/5.0;

        if(0.0<=gap_x && gap_x<L/3.0)
            B3[0] = 1.0;
        else if(L/3.0<=gap_x && gap_x<L/2.0)
            B3[0] = (L/2.0 - gap_x)/(L/6.0);
        else if(2.0*L/3.0<=gap_x && gap_x<=L)
            B3[0] = 0.0;
        else
            std::cout<<"gap_x is error";
        
        if(0.0<=gap_x && gap_x<L/3.0)
            B3[1] = 0.0;
        else if(L/3.0<=gap_x && gap_x<L/2.0)
            B3[1] = 1 - (L/2.0 - gap_x)/(L/6.0);
        else if(L/2.0<=gap_x && gap_x<2.0*L/3.0)
            B3[1] = 1 - (gap_x - L/2.0)/(L/6.0);
        else if(2.0*L/3.0<=gap_x && gap_x<=L)
            B3[1] = 0.0;
        else
            std::cout<<"gap_x is error";
    
        if(0<=gap_x && gap_x<L/2.0)
            B3[2] = 0.0;
        else if(L/2.0<=gap_x && gap_x<2.0*L/3.0)
            B3[2] = (gap_x - L/2.0)/(L/6.0);
        else if(2.0*L/3.0<=gap_x && gap_x <= L)
            B3[2] = 1.0;
        else 
            std::cout<<"gap_x is error";

        // B4 己方与球的距离
        double distance_sum = 0, mean_dis;
        for(std::size_t i = 0 ; i < OUR_TEAM ; i++)
            distance_sum += DPoint(pose_us[i][0],pose_us[i][1]).distance(position_ball);
        mean_dis = distance_sum/OUR_TEAM;

        if(0.0<=mean_dis && mean_dis<5.0*L/8.0)
            B4[0] = 0.0;
        else if(5.0*L/8.0<=mean_dis && mean_dis<7.0*L/8.0)
            B4[0] = (mean_dis - 5.0*L/8.0)/(L/4.0);
        else if(7.0*L/8.0<=mean_dis && mean_dis<=L)
            B4[0] = 1.0;
        else
            std::cout<<"mean_dis is error!";
        
        if(0.0<=mean_dis && mean_dis<L/8.0)
            B4[1] = 0.0;
        else if(L/8.0<=mean_dis && mean_dis<3.0*L/8.0)
            B4[1] = 1 - (3.0*L/8.0 - mean_dis)/(L/4.0);
        else if(3.0*L/8.0<=mean_dis && mean_dis<5.0*L/8.0)
            B4[1] = 1 ;
        else if(5.0*L/8.0<=mean_dis && mean_dis<7.0*L/8.0)
            B4[1] = 1.0 - (mean_dis - 5.0*L/8.0)/(L/4.0);
        else if(7.0*L/8.0<=mean_dis && mean_dis <=L)
            B4[1] = 0.0;
        else
            std::cout<<"mean_dis is error!";
    
        if(0<=mean_dis && mean_dis<L/8.0)
            B4[2] = 1.0;
        else if(L/8.0<=mean_dis && mean_dis<3.0*L/8.0)
            B4[2] = (3.0*L/8.0 - mean_dis)/(L/4.0);
        else if(3.0*L/8.0<=mean_dis && mean_dis <= L)
            B4[2] = 0.0;
        else 
            std::cout<<"mean_dis is error!";
        
        //B5 己方朝向
        double deta, mean_angle;
        for(std::size_t i = 0 ; i < OUR_TEAM ; i++)
        {
            deta = abs(DPoint(pose_us[i][0],pose_us[i][1]).angle(position_ball).radian_);
            if(deta>180) deta=360-deta;
            mean_angle += deta/OUR_TEAM;
        }
        
        if(0.0<=mean_angle && mean_angle<45.0)
            B5[0] = 0.0;
        else if(45.0<=mean_angle && mean_angle<60.0)
            B5[0] = (mean_angle - 45.0)/15.0;
        else if(60<=mean_angle && mean_angle<=180)
            B5[0] = 1.0;
        else
            std::cout<<"mean_angle is error!";
        
        if(0.0<=mean_angle && mean_angle<15.0)
            B5[1] = 0.0;
        else if(15.0<=mean_angle && mean_angle<45.0)
            B5[1] = 1 - (45.0 - mean_angle)/30.0;
        else if(45.0<=mean_angle && mean_angle<60.0)
            B5[1] = 1 - (mean_angle - 45)/15.0;
        else if(60.0<=mean_angle && mean_angle<=180.0)
            B5[1] = 0.0;
        else
            std::cout<<"mean_angle is error!";
    
        if(0<=mean_angle && mean_angle<15.0)
            B5[2] = 1.0;
        else if(15.0<=mean_angle && mean_angle<45.0)
            B5[2] = (45.0 - mean_angle)/30.0;
        else if(45.0<=mean_angle && mean_angle <= 180.0)
            B5[2] = 0.0;
        else 
            std::cout<<"mean_angle is error!";

        // B6 敌方与球的距离
        distance_sum = 0.0;
        for(std::size_t i = 0 ; i < world_model_->Opponents_.size() ; i++)
            distance_sum += DPoint(pose_oppo[i][0],pose_oppo[i][1]).distance(position_ball);
        mean_dis = distance_sum/OUR_TEAM;

        if(0.0<=mean_dis && mean_dis<L/8.0)
            B6[0] = 1.0;
        else if(L/8.0<=mean_dis && mean_dis<3.0*L/8.0)
            B6[0] = (3.0*L/8.0 - mean_dis)/(L/4.0);
        else if(3.0*L/8.0<=mean_dis && mean_dis<=L)
            B6[0] = 0.0;
        else
            std::cout<<"mean_dis is error!";
        
        if(0.0<=mean_dis && mean_dis<L/8.0)
            B6[1] = 0.0;
        else if(L/8.0<=mean_dis && mean_dis<3.0*L/8.0)
            B6[1] = (3.0*L/8.0 - mean_dis)/(L/4.0);
        else if(3.0*L/8.0<=mean_dis && mean_dis<5.0*L/8.0)
            B6[1] = 1.0 ;
        else if(5.0*L/8.0<=mean_dis && mean_dis<7.0*L/8.0)
            B6[1] = (mean_dis - 5.0*L/8.0)/(L/4.0);
        else if(7.0*L/8.0<=mean_dis && mean_dis <=L)
            B6[1] = 0.0;
        else
            std::cout<<"mean_dis is error!";
    
        if(0<=mean_dis && mean_dis<5.0*L/8.0)
            B6[2] = 0.0;
        else if(5.0*L/8.0<=mean_dis && mean_dis<7.0*L/8.0)
            B6[2] = (mean_dis - 5.0*L/8.0)/(L/4.0);
        else if(7.0*L/8.0<=mean_dis && mean_dis <= L)
            B6[2] = 1.0;
        else 
            std::cout<<"mean_dis is error!";
        
        // 计分
        double max = 0;
        for(std::size_t i = 0 ; i < 3 ; i++)
        {
            score[i]=weight[0]*B1[i]+weight[1]*B2[i]+weight[2]*B1[3]+weight[3]*B4[i]+weight[4]*B5[i]+weight[5]*B6[i];
            if(score[i]>max)
                {max=score[i];status=i+1;}
        }
    }


void Strategy::Market_act()
{
    double L=2200, H=1400, k1=0.6, k2=0.4, L_us, angle_us_ball;
    price[0][0]=1000;
    for(std::size_t i = 1 ; i < 5 ; i++)
    {
        L_us = abs(-L/2.0 - pose_us[i][0]);
        angle_us_ball = abs(DPoint(pose_us[i][0],pose_us[i][1]).angle(position_ball).radian_);
        price[i][0] = 1.0/(k1*(L_us/L) + k2*(angle_us_ball/(180/2)));
    }
    price[0][1]=0;
    k1=0.2875, k2=0.2125;
    double k3=0.5, k4=1, para[3], beta=0, sita, L_ball_robot, L_diag, angle_yaw_ball, theta_us_ball_goal;
    // ID of goalkeeper is 1, serial is 0
    for(std::size_t i = 1 ; i < 5 ; i++)
    {
        // k1
        L_ball_robot = DPoint(pose_us[i][0],pose_us[i][1]).distance(position_ball);
        L_diag = sqrt(L*L + H*H);
        // k2
        angle_yaw_ball = abs(pose_us[i][2] - DPoint(pose_us[i][0],pose_us[i][1]).angle(position_ball).radian_);
        if(angle_yaw_ball >180.0)
            angle_yaw_ball = 360.0 - angle_yaw_ball;
        // k3
        theta_us_ball_goal = abs(DPoint(1100,0).angle(position_ball).radian_ - 
                                DPoint(pose_us[i][0],pose_us[i][1]).angle(position_ball).radian_);
        if(theta_us_ball_goal > 180.0)
            theta_us_ball_goal = 360.0 - theta_us_ball_goal;
        // beta
        beta=0;
        for(std::size_t j = 1 ; j < 5 ; j++)
        {
            Obstacle[j-1]=DPoint(pose_us[j][0],pose_us[j][1]);
            Obstacle[j+3]=DPoint(pose_oppo[j][0],pose_oppo[j][1]);
        }
        get_line( para, DPoint(pose_us[i][0],pose_us[i][1]), position_ball);
        for(std::size_t j = 0 ; j < world_model_->Opponents_.size()+OUR_TEAM-2 ; j++)
        {
            if(j!=i-1)continue;
            if(para[0]==1)
            {
                sita = atan(para[1]);
                if(sita<0)sita = sita + 3.1415926;
                double y1 = para[1] * Obstacle[j].x_ + para[2] + 0.5/cos(sita);
                double y2 = para[1] * Obstacle[j].x_ + para[2] - 0.5/cos(sita);
                if(((y1 - Obstacle[j].y_)*(y2 - Obstacle[j].y_) <= 0) && Obstacle[j].x_ >= min(pose_us[i][0], position_ball.x_) && Obstacle[j].x_ <= max(pose_us[i][0], position_ball.x_) && Obstacle[j].y_ >= min(pose_us[i][1], position_ball.y_) && Obstacle[j].y_ <= max(pose_us[i][1], position_ball.y_))
                {
                    beta = 1.0*k4;
                    break;
                }
            }
            else if(para[0]==2)
            {
                if((para[1] - 0.5) <= Obstacle[j].x_ && Obstacle[j].x_ <=(para[1] + 0.5) &&  Obstacle[j].y_ >= min(pose_us[1][1], position_ball.y_) && Obstacle[j].y_ <= max(pose_us[i][1], position_ball.y_))
                {
                    beta = 1.0*k4;
                    break;
                }
            }
            else
            {
                beta = 1.0*k4;
            }
            price[i][1] = 1.0/(k1*(L_ball_robot/L_diag) + k2*(abs(angle_yaw_ball)/180.0) + k3*(theta_us_ball_goal)/180.0 + beta);
        }
    }
}
void Strategy::get_line(double*para, DPoint2d pose, DPoint pos_ball)
{
    double dx, dy;
    dx = pos_ball.x_-pose.x_;
    dy = pos_ball.y_-pose.y_;
    if(dx!=0)
    {
        para[0]=1;
        para[1]=tan(dy/dx);
        para[2]=pose.y_-para[1]*pose.x_;
    }
    else if(dy==0)
    {
        para[0]=2;
        para[1]=pose.x_;}
    else
        para[0]=-1;    
}

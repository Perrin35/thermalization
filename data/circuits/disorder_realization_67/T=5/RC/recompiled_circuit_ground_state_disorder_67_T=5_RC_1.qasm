OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9031653) q[0];
sx q[0];
rz(-0.16024661) q[0];
sx q[0];
rz(1.8538374) q[0];
rz(-0.84914485) q[1];
sx q[1];
rz(-1.1916817) q[1];
sx q[1];
rz(-2.3828659) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37739524) q[0];
sx q[0];
rz(-1.8599417) q[0];
sx q[0];
rz(1.8125066) q[0];
x q[1];
rz(-2.0082194) q[2];
sx q[2];
rz(-2.6133399) q[2];
sx q[2];
rz(0.6483486) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.030144) q[1];
sx q[1];
rz(-1.2093102) q[1];
sx q[1];
rz(-2.3356781) q[1];
rz(-pi) q[2];
rz(2.4843744) q[3];
sx q[3];
rz(-2.2003678) q[3];
sx q[3];
rz(-2.8384125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8410926) q[2];
sx q[2];
rz(-0.79885834) q[2];
sx q[2];
rz(2.63499) q[2];
rz(1.5233585) q[3];
sx q[3];
rz(-0.62464276) q[3];
sx q[3];
rz(-3.0860331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8697206) q[0];
sx q[0];
rz(-2.6756918) q[0];
sx q[0];
rz(3.077935) q[0];
rz(-1.1977389) q[1];
sx q[1];
rz(-1.5143737) q[1];
sx q[1];
rz(1.0187842) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5658582) q[0];
sx q[0];
rz(-1.3745148) q[0];
sx q[0];
rz(-1.422276) q[0];
x q[1];
rz(2.8009861) q[2];
sx q[2];
rz(-1.34769) q[2];
sx q[2];
rz(2.4702213) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9460822) q[1];
sx q[1];
rz(-1.2942593) q[1];
sx q[1];
rz(0.72317381) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0080454) q[3];
sx q[3];
rz(-2.3192647) q[3];
sx q[3];
rz(-2.453367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7055052) q[2];
sx q[2];
rz(-2.5477396) q[2];
sx q[2];
rz(2.2994821) q[2];
rz(-1.0574794) q[3];
sx q[3];
rz(-0.92856854) q[3];
sx q[3];
rz(1.9717982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8975526) q[0];
sx q[0];
rz(-1.0201447) q[0];
sx q[0];
rz(-1.9496339) q[0];
rz(-0.92487088) q[1];
sx q[1];
rz(-2.4332739) q[1];
sx q[1];
rz(-1.301773) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90067139) q[0];
sx q[0];
rz(-2.0427454) q[0];
sx q[0];
rz(-0.94743418) q[0];
rz(-pi) q[1];
rz(1.5315751) q[2];
sx q[2];
rz(-1.3717781) q[2];
sx q[2];
rz(-1.7165556) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9684549) q[1];
sx q[1];
rz(-2.0357642) q[1];
sx q[1];
rz(-1.370621) q[1];
rz(-0.93940063) q[3];
sx q[3];
rz(-0.41920788) q[3];
sx q[3];
rz(1.7260176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3785582) q[2];
sx q[2];
rz(-0.083746567) q[2];
sx q[2];
rz(-3.0540826) q[2];
rz(-1.3148426) q[3];
sx q[3];
rz(-2.0851236) q[3];
sx q[3];
rz(0.99353138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73115504) q[0];
sx q[0];
rz(-1.1719828) q[0];
sx q[0];
rz(-2.3325969) q[0];
rz(-1.0357098) q[1];
sx q[1];
rz(-0.46329841) q[1];
sx q[1];
rz(-2.1083924) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3952156) q[0];
sx q[0];
rz(-2.4762634) q[0];
sx q[0];
rz(1.6341342) q[0];
x q[1];
rz(0.90649267) q[2];
sx q[2];
rz(-0.64721738) q[2];
sx q[2];
rz(2.3075019) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3174008) q[1];
sx q[1];
rz(-0.5017952) q[1];
sx q[1];
rz(-2.7609894) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5979126) q[3];
sx q[3];
rz(-1.9900232) q[3];
sx q[3];
rz(2.4329684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0674627) q[2];
sx q[2];
rz(-2.9758657) q[2];
sx q[2];
rz(1.6039675) q[2];
rz(1.4175203) q[3];
sx q[3];
rz(-1.1105024) q[3];
sx q[3];
rz(-2.0541151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23987016) q[0];
sx q[0];
rz(-0.19342315) q[0];
sx q[0];
rz(-1.1892009) q[0];
rz(2.4173648) q[1];
sx q[1];
rz(-1.8502356) q[1];
sx q[1];
rz(-1.6667746) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3132202) q[0];
sx q[0];
rz(-0.68183866) q[0];
sx q[0];
rz(0.71558715) q[0];
rz(-pi) q[1];
rz(1.7314265) q[2];
sx q[2];
rz(-0.77996563) q[2];
sx q[2];
rz(0.51973625) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2853299) q[1];
sx q[1];
rz(-1.26795) q[1];
sx q[1];
rz(-1.9252127) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.311077) q[3];
sx q[3];
rz(-0.6336824) q[3];
sx q[3];
rz(0.049402852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7597947) q[2];
sx q[2];
rz(-2.2319824) q[2];
sx q[2];
rz(2.6749715) q[2];
rz(-1.7032547) q[3];
sx q[3];
rz(-1.2983026) q[3];
sx q[3];
rz(2.2013825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4284215) q[0];
sx q[0];
rz(-2.3155825) q[0];
sx q[0];
rz(3.1312422) q[0];
rz(1.6185282) q[1];
sx q[1];
rz(-1.7489988) q[1];
sx q[1];
rz(-1.7656322) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.812777) q[0];
sx q[0];
rz(-1.3294018) q[0];
sx q[0];
rz(-1.9184979) q[0];
rz(2.9508123) q[2];
sx q[2];
rz(-1.6990802) q[2];
sx q[2];
rz(0.14036638) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.048307) q[1];
sx q[1];
rz(-1.3897589) q[1];
sx q[1];
rz(0.50706086) q[1];
rz(0.74919219) q[3];
sx q[3];
rz(-1.0674879) q[3];
sx q[3];
rz(-1.0096514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7213664) q[2];
sx q[2];
rz(-2.0927636) q[2];
sx q[2];
rz(0.72727195) q[2];
rz(-0.51491245) q[3];
sx q[3];
rz(-1.3313096) q[3];
sx q[3];
rz(1.4179199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10776831) q[0];
sx q[0];
rz(-1.533778) q[0];
sx q[0];
rz(-2.7346101) q[0];
rz(2.175323) q[1];
sx q[1];
rz(-0.85710183) q[1];
sx q[1];
rz(0.13872096) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4086558) q[0];
sx q[0];
rz(-1.4545024) q[0];
sx q[0];
rz(3.1058806) q[0];
rz(-pi) q[1];
rz(-0.27757711) q[2];
sx q[2];
rz(-1.7834366) q[2];
sx q[2];
rz(2.9523498) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3073716) q[1];
sx q[1];
rz(-1.2117996) q[1];
sx q[1];
rz(-1.9407843) q[1];
x q[2];
rz(1.1769663) q[3];
sx q[3];
rz(-2.7714249) q[3];
sx q[3];
rz(0.98370508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2762642) q[2];
sx q[2];
rz(-1.3778957) q[2];
sx q[2];
rz(-0.31065568) q[2];
rz(1.8336512) q[3];
sx q[3];
rz(-1.6303948) q[3];
sx q[3];
rz(2.7697897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0528316) q[0];
sx q[0];
rz(-0.3449057) q[0];
sx q[0];
rz(0.088223591) q[0];
rz(-0.010559646) q[1];
sx q[1];
rz(-1.5225531) q[1];
sx q[1];
rz(0.47058502) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4751401) q[0];
sx q[0];
rz(-1.0658403) q[0];
sx q[0];
rz(2.3078437) q[0];
rz(-0.80704524) q[2];
sx q[2];
rz(-1.7084439) q[2];
sx q[2];
rz(2.0375843) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7631665) q[1];
sx q[1];
rz(-1.4861123) q[1];
sx q[1];
rz(1.3099109) q[1];
rz(-pi) q[2];
rz(-1.9235189) q[3];
sx q[3];
rz(-1.7021057) q[3];
sx q[3];
rz(2.2985947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18870246) q[2];
sx q[2];
rz(-2.1850093) q[2];
sx q[2];
rz(1.0799705) q[2];
rz(-0.88932577) q[3];
sx q[3];
rz(-1.5417128) q[3];
sx q[3];
rz(1.5842452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66460669) q[0];
sx q[0];
rz(-2.7598858) q[0];
sx q[0];
rz(-0.81740022) q[0];
rz(-0.74642247) q[1];
sx q[1];
rz(-2.4669929) q[1];
sx q[1];
rz(-2.2019763) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8864) q[0];
sx q[0];
rz(-1.5879244) q[0];
sx q[0];
rz(3.1411132) q[0];
rz(-pi) q[1];
rz(-2.4822761) q[2];
sx q[2];
rz(-2.7020279) q[2];
sx q[2];
rz(-2.2826113) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0683953) q[1];
sx q[1];
rz(-1.4540352) q[1];
sx q[1];
rz(1.8224656) q[1];
rz(-pi) q[2];
rz(1.1577179) q[3];
sx q[3];
rz(-1.7677757) q[3];
sx q[3];
rz(0.22991163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8115936) q[2];
sx q[2];
rz(-0.96572319) q[2];
sx q[2];
rz(2.7962371) q[2];
rz(-1.6164448) q[3];
sx q[3];
rz(-1.350178) q[3];
sx q[3];
rz(-2.7143872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78373194) q[0];
sx q[0];
rz(-1.9851728) q[0];
sx q[0];
rz(3.0586437) q[0];
rz(0.16231617) q[1];
sx q[1];
rz(-1.4444618) q[1];
sx q[1];
rz(-2.4457795) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9385949) q[0];
sx q[0];
rz(-1.8502529) q[0];
sx q[0];
rz(1.6142815) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4415353) q[2];
sx q[2];
rz(-1.5678355) q[2];
sx q[2];
rz(-1.2128613) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3491427) q[1];
sx q[1];
rz(-0.44466838) q[1];
sx q[1];
rz(0.044388958) q[1];
rz(-pi) q[2];
rz(1.8553461) q[3];
sx q[3];
rz(-0.88594809) q[3];
sx q[3];
rz(-1.7189792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7210228) q[2];
sx q[2];
rz(-0.15779725) q[2];
sx q[2];
rz(-1.2201307) q[2];
rz(-2.5681791) q[3];
sx q[3];
rz(-1.8107332) q[3];
sx q[3];
rz(2.8964608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86380105) q[0];
sx q[0];
rz(-1.6391123) q[0];
sx q[0];
rz(3.0085051) q[0];
rz(-0.0028903891) q[1];
sx q[1];
rz(-2.7255701) q[1];
sx q[1];
rz(2.4400673) q[1];
rz(-2.1720282) q[2];
sx q[2];
rz(-2.8429016) q[2];
sx q[2];
rz(0.071879172) q[2];
rz(0.057249476) q[3];
sx q[3];
rz(-0.80116873) q[3];
sx q[3];
rz(2.7309668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

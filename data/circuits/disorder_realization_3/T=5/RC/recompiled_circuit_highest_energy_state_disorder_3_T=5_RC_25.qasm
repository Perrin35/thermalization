OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1138104) q[0];
sx q[0];
rz(-1.8960928) q[0];
sx q[0];
rz(-1.0863289) q[0];
rz(2.0431986) q[1];
sx q[1];
rz(-1.9280704) q[1];
sx q[1];
rz(-1.1462125) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1725464) q[0];
sx q[0];
rz(-1.2451225) q[0];
sx q[0];
rz(-2.9708751) q[0];
x q[1];
rz(-2.5511247) q[2];
sx q[2];
rz(-1.9158949) q[2];
sx q[2];
rz(-2.9546628) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0806345) q[1];
sx q[1];
rz(-0.74392156) q[1];
sx q[1];
rz(-0.77772909) q[1];
rz(-1.0686773) q[3];
sx q[3];
rz(-1.7699336) q[3];
sx q[3];
rz(-0.56875436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9132797) q[2];
sx q[2];
rz(-1.5756807) q[2];
sx q[2];
rz(0.3621462) q[2];
rz(2.1892138) q[3];
sx q[3];
rz(-2.5824472) q[3];
sx q[3];
rz(0.70498103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1700965) q[0];
sx q[0];
rz(-2.9269452) q[0];
sx q[0];
rz(1.5317408) q[0];
rz(1.8929298) q[1];
sx q[1];
rz(-1.5208533) q[1];
sx q[1];
rz(0.85266399) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0097875) q[0];
sx q[0];
rz(-0.62776792) q[0];
sx q[0];
rz(0.96012299) q[0];
rz(1.8478099) q[2];
sx q[2];
rz(-1.4528034) q[2];
sx q[2];
rz(-1.5548776) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.98140796) q[1];
sx q[1];
rz(-1.7716737) q[1];
sx q[1];
rz(-0.9717047) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1308822) q[3];
sx q[3];
rz(-1.2398071) q[3];
sx q[3];
rz(2.6382642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2109005) q[2];
sx q[2];
rz(-0.96052581) q[2];
sx q[2];
rz(1.9453913) q[2];
rz(0.026737468) q[3];
sx q[3];
rz(-2.6452711) q[3];
sx q[3];
rz(-1.668821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6579599) q[0];
sx q[0];
rz(-2.8909029) q[0];
sx q[0];
rz(0.87580097) q[0];
rz(-0.35975131) q[1];
sx q[1];
rz(-1.7774372) q[1];
sx q[1];
rz(-2.13805) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10570603) q[0];
sx q[0];
rz(-1.0125475) q[0];
sx q[0];
rz(1.4869177) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0613173) q[2];
sx q[2];
rz(-0.70827863) q[2];
sx q[2];
rz(0.88513155) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4565125) q[1];
sx q[1];
rz(-2.4488423) q[1];
sx q[1];
rz(2.9167006) q[1];
x q[2];
rz(2.4677213) q[3];
sx q[3];
rz(-1.1316852) q[3];
sx q[3];
rz(-1.278745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2063736) q[2];
sx q[2];
rz(-2.424365) q[2];
sx q[2];
rz(0.75500542) q[2];
rz(-3.0435009) q[3];
sx q[3];
rz(-0.57049975) q[3];
sx q[3];
rz(-0.3977972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0841242) q[0];
sx q[0];
rz(-1.5115154) q[0];
sx q[0];
rz(-2.496642) q[0];
rz(-1.0954674) q[1];
sx q[1];
rz(-2.7332833) q[1];
sx q[1];
rz(1.1114978) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0065466) q[0];
sx q[0];
rz(-1.5113976) q[0];
sx q[0];
rz(-3.0778372) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6531282) q[2];
sx q[2];
rz(-1.7373184) q[2];
sx q[2];
rz(-2.4660267) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9330628) q[1];
sx q[1];
rz(-2.8840319) q[1];
sx q[1];
rz(2.527586) q[1];
rz(-pi) q[2];
rz(0.74097775) q[3];
sx q[3];
rz(-1.40687) q[3];
sx q[3];
rz(-1.0865097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5081386) q[2];
sx q[2];
rz(-2.3544669) q[2];
sx q[2];
rz(0.48466361) q[2];
rz(2.700452) q[3];
sx q[3];
rz(-1.7287247) q[3];
sx q[3];
rz(-1.887656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50705528) q[0];
sx q[0];
rz(-0.84369722) q[0];
sx q[0];
rz(1.7150568) q[0];
rz(2.2818395) q[1];
sx q[1];
rz(-0.25329241) q[1];
sx q[1];
rz(-0.8304798) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5317212) q[0];
sx q[0];
rz(-0.97519704) q[0];
sx q[0];
rz(-2.2723212) q[0];
x q[1];
rz(2.203622) q[2];
sx q[2];
rz(-1.5245078) q[2];
sx q[2];
rz(-0.37424311) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1089835) q[1];
sx q[1];
rz(-2.2936037) q[1];
sx q[1];
rz(1.384114) q[1];
x q[2];
rz(0.40404218) q[3];
sx q[3];
rz(-2.7528304) q[3];
sx q[3];
rz(-2.5498912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6972203) q[2];
sx q[2];
rz(-1.0151007) q[2];
sx q[2];
rz(0.74697906) q[2];
rz(0.30484453) q[3];
sx q[3];
rz(-1.3957142) q[3];
sx q[3];
rz(-1.2386809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9656669) q[0];
sx q[0];
rz(-1.8400064) q[0];
sx q[0];
rz(2.9081705) q[0];
rz(-0.4625136) q[1];
sx q[1];
rz(-1.9513444) q[1];
sx q[1];
rz(-2.7375284) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44495917) q[0];
sx q[0];
rz(-0.6650266) q[0];
sx q[0];
rz(1.2103266) q[0];
x q[1];
rz(-1.839371) q[2];
sx q[2];
rz(-1.982455) q[2];
sx q[2];
rz(1.1546763) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5517496) q[1];
sx q[1];
rz(-2.1300364) q[1];
sx q[1];
rz(-1.4435121) q[1];
rz(-pi) q[2];
rz(2.5109291) q[3];
sx q[3];
rz(-1.3913097) q[3];
sx q[3];
rz(1.6398182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.81249753) q[2];
sx q[2];
rz(-1.3500682) q[2];
sx q[2];
rz(-2.4883032) q[2];
rz(-1.6434044) q[3];
sx q[3];
rz(-2.523962) q[3];
sx q[3];
rz(-2.5352855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46635982) q[0];
sx q[0];
rz(-1.2132069) q[0];
sx q[0];
rz(-2.8208222) q[0];
rz(-0.65387154) q[1];
sx q[1];
rz(-0.32952285) q[1];
sx q[1];
rz(-3.0810862) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084235926) q[0];
sx q[0];
rz(-0.31937803) q[0];
sx q[0];
rz(1.0127629) q[0];
rz(-pi) q[1];
rz(2.4662911) q[2];
sx q[2];
rz(-0.85682288) q[2];
sx q[2];
rz(1.7973822) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6459319) q[1];
sx q[1];
rz(-1.0905318) q[1];
sx q[1];
rz(-2.2731794) q[1];
rz(-pi) q[2];
rz(-0.860608) q[3];
sx q[3];
rz(-2.6250556) q[3];
sx q[3];
rz(2.9252657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3939646) q[2];
sx q[2];
rz(-1.4614033) q[2];
sx q[2];
rz(0.72236577) q[2];
rz(-2.8584621) q[3];
sx q[3];
rz(-0.88909283) q[3];
sx q[3];
rz(0.50962555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79001456) q[0];
sx q[0];
rz(-0.55792648) q[0];
sx q[0];
rz(-0.18490069) q[0];
rz(-2.6376873) q[1];
sx q[1];
rz(-1.3495219) q[1];
sx q[1];
rz(-0.46461836) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4330209) q[0];
sx q[0];
rz(-1.4056217) q[0];
sx q[0];
rz(-0.1689706) q[0];
rz(-pi) q[1];
rz(-2.9275168) q[2];
sx q[2];
rz(-0.42595562) q[2];
sx q[2];
rz(-1.730012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2459348) q[1];
sx q[1];
rz(-0.57809752) q[1];
sx q[1];
rz(-0.74191414) q[1];
rz(2.9161386) q[3];
sx q[3];
rz(-2.4092393) q[3];
sx q[3];
rz(-1.454753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9087312) q[2];
sx q[2];
rz(-1.6249012) q[2];
sx q[2];
rz(-2.0880584) q[2];
rz(-0.72707027) q[3];
sx q[3];
rz(-1.9126242) q[3];
sx q[3];
rz(-0.54522902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9497249) q[0];
sx q[0];
rz(-1.1911012) q[0];
sx q[0];
rz(0.18549347) q[0];
rz(-2.5718451) q[1];
sx q[1];
rz(-2.2173917) q[1];
sx q[1];
rz(1.1572908) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3678903) q[0];
sx q[0];
rz(-1.7130525) q[0];
sx q[0];
rz(2.7223047) q[0];
x q[1];
rz(1.6389388) q[2];
sx q[2];
rz(-0.63387094) q[2];
sx q[2];
rz(-2.9958452) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9461575) q[1];
sx q[1];
rz(-1.3570254) q[1];
sx q[1];
rz(0.56758387) q[1];
rz(2.7140541) q[3];
sx q[3];
rz(-2.2898009) q[3];
sx q[3];
rz(-2.2710272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40413228) q[2];
sx q[2];
rz(-2.9170333) q[2];
sx q[2];
rz(1.7411211) q[2];
rz(2.776966) q[3];
sx q[3];
rz(-2.3783763) q[3];
sx q[3];
rz(2.3188685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51065651) q[0];
sx q[0];
rz(-0.83844227) q[0];
sx q[0];
rz(-3.1070218) q[0];
rz(0.336054) q[1];
sx q[1];
rz(-2.2269109) q[1];
sx q[1];
rz(-2.5539982) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60025763) q[0];
sx q[0];
rz(-2.6980639) q[0];
sx q[0];
rz(-2.8044273) q[0];
rz(-pi) q[1];
rz(0.35878344) q[2];
sx q[2];
rz(-2.3780895) q[2];
sx q[2];
rz(0.43649188) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42660824) q[1];
sx q[1];
rz(-2.9054925) q[1];
sx q[1];
rz(-1.3900422) q[1];
rz(-pi) q[2];
rz(2.0301271) q[3];
sx q[3];
rz(-0.94375347) q[3];
sx q[3];
rz(0.88739838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8079638) q[2];
sx q[2];
rz(-0.89666349) q[2];
sx q[2];
rz(0.72820747) q[2];
rz(0.41094574) q[3];
sx q[3];
rz(-2.9185037) q[3];
sx q[3];
rz(1.2724916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67644607) q[0];
sx q[0];
rz(-1.5560173) q[0];
sx q[0];
rz(-1.7898855) q[0];
rz(0.43988718) q[1];
sx q[1];
rz(-0.37275795) q[1];
sx q[1];
rz(2.8511924) q[1];
rz(-1.8710623) q[2];
sx q[2];
rz(-1.9004442) q[2];
sx q[2];
rz(1.6729542) q[2];
rz(1.6174118) q[3];
sx q[3];
rz(-1.4641932) q[3];
sx q[3];
rz(-2.1144397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

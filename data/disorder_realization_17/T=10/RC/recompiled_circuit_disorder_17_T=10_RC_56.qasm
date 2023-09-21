OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2330115) q[0];
sx q[0];
rz(-1.1865948) q[0];
sx q[0];
rz(-1.9957805) q[0];
rz(0.97283483) q[1];
sx q[1];
rz(-1.4714779) q[1];
sx q[1];
rz(0.29247984) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2746409) q[0];
sx q[0];
rz(-0.39817087) q[0];
sx q[0];
rz(-1.7142332) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5195261) q[2];
sx q[2];
rz(-0.79203696) q[2];
sx q[2];
rz(-0.37026065) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7293207) q[1];
sx q[1];
rz(-2.1643157) q[1];
sx q[1];
rz(-3.0862942) q[1];
rz(-pi) q[2];
rz(-0.0019449751) q[3];
sx q[3];
rz(-2.0446152) q[3];
sx q[3];
rz(-0.67809826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4108489) q[2];
sx q[2];
rz(-1.0935254) q[2];
sx q[2];
rz(2.6921819) q[2];
rz(2.4959026) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(1.8507563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9280424) q[0];
sx q[0];
rz(-0.88590652) q[0];
sx q[0];
rz(2.399562) q[0];
rz(1.4713326) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(-1.7785633) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3516386) q[0];
sx q[0];
rz(-0.75766701) q[0];
sx q[0];
rz(2.555048) q[0];
x q[1];
rz(-1.7861373) q[2];
sx q[2];
rz(-1.1110348) q[2];
sx q[2];
rz(0.41972566) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.14582536) q[1];
sx q[1];
rz(-1.9103423) q[1];
sx q[1];
rz(-1.5911566) q[1];
rz(-pi) q[2];
rz(1.8369254) q[3];
sx q[3];
rz(-1.3849272) q[3];
sx q[3];
rz(0.64418018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.30119511) q[2];
sx q[2];
rz(-2.6448554) q[2];
sx q[2];
rz(-1.8691241) q[2];
rz(0.33723351) q[3];
sx q[3];
rz(-1.971258) q[3];
sx q[3];
rz(-0.01005323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0740046) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(0.2628251) q[0];
rz(2.7858531) q[1];
sx q[1];
rz(-1.6825312) q[1];
sx q[1];
rz(-1.9062818) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10997406) q[0];
sx q[0];
rz(-1.585916) q[0];
sx q[0];
rz(1.5605687) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4624865) q[2];
sx q[2];
rz(-0.24178594) q[2];
sx q[2];
rz(2.0568648) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3706627) q[1];
sx q[1];
rz(-1.8924894) q[1];
sx q[1];
rz(1.6817723) q[1];
rz(-0.87576207) q[3];
sx q[3];
rz(-1.2998298) q[3];
sx q[3];
rz(2.392644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0170903) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(0.18033218) q[2];
rz(0.63594937) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(-2.7699871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76698774) q[0];
sx q[0];
rz(-2.7174482) q[0];
sx q[0];
rz(-0.71907991) q[0];
rz(0.51302296) q[1];
sx q[1];
rz(-1.2027556) q[1];
sx q[1];
rz(-2.1069353) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661515) q[0];
sx q[0];
rz(-0.062128566) q[0];
sx q[0];
rz(-2.8176869) q[0];
x q[1];
rz(3.0940975) q[2];
sx q[2];
rz(-1.4357899) q[2];
sx q[2];
rz(-2.6035655) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82995854) q[1];
sx q[1];
rz(-1.8533851) q[1];
sx q[1];
rz(0.71210536) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47311584) q[3];
sx q[3];
rz(-1.72662) q[3];
sx q[3];
rz(1.6254049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2085312) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(-0.2362403) q[2];
rz(1.158372) q[3];
sx q[3];
rz(-2.1393496) q[3];
sx q[3];
rz(-0.85197824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029595705) q[0];
sx q[0];
rz(-2.1313666) q[0];
sx q[0];
rz(0.90233666) q[0];
rz(-2.2205655) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(0.62228084) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3233933) q[0];
sx q[0];
rz(-0.33938956) q[0];
sx q[0];
rz(-2.1470977) q[0];
rz(1.4786517) q[2];
sx q[2];
rz(-0.92458692) q[2];
sx q[2];
rz(0.19942936) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.2052853) q[1];
sx q[1];
rz(-0.19166066) q[1];
sx q[1];
rz(-1.081341) q[1];
rz(-pi) q[2];
rz(2.4683687) q[3];
sx q[3];
rz(-2.554318) q[3];
sx q[3];
rz(-2.4622963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2174012) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(3.1203111) q[2];
rz(1.4422669) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(-0.91059476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845881) q[0];
sx q[0];
rz(-0.5963043) q[0];
sx q[0];
rz(-3.0969627) q[0];
rz(2.1394829) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(-0.034428509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32193004) q[0];
sx q[0];
rz(-0.18580431) q[0];
sx q[0];
rz(1.4968605) q[0];
rz(-pi) q[1];
rz(-3.0197633) q[2];
sx q[2];
rz(-1.0708772) q[2];
sx q[2];
rz(0.027539754) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1457445) q[1];
sx q[1];
rz(-1.8954344) q[1];
sx q[1];
rz(-2.1374628) q[1];
rz(-2.9488726) q[3];
sx q[3];
rz(-0.7145213) q[3];
sx q[3];
rz(-2.2477828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.78952152) q[2];
sx q[2];
rz(-1.8692724) q[2];
sx q[2];
rz(0.96674031) q[2];
rz(0.012185193) q[3];
sx q[3];
rz(-2.2120078) q[3];
sx q[3];
rz(-3.0325082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0063342) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(-0.99408856) q[0];
rz(-0.055123568) q[1];
sx q[1];
rz(-1.8693285) q[1];
sx q[1];
rz(-0.73928839) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3294322) q[0];
sx q[0];
rz(-1.433916) q[0];
sx q[0];
rz(3.0457004) q[0];
rz(1.0871068) q[2];
sx q[2];
rz(-1.0488044) q[2];
sx q[2];
rz(-2.1453478) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.90034396) q[1];
sx q[1];
rz(-0.56300357) q[1];
sx q[1];
rz(2.7808933) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7474653) q[3];
sx q[3];
rz(-0.65791241) q[3];
sx q[3];
rz(3.1042299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58549515) q[2];
sx q[2];
rz(-1.408351) q[2];
sx q[2];
rz(0.56376702) q[2];
rz(-2.901315) q[3];
sx q[3];
rz(-0.53301817) q[3];
sx q[3];
rz(-1.7668004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4450842) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(0.58404303) q[0];
rz(2.7208327) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(2.8631794) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95094341) q[0];
sx q[0];
rz(-1.6961432) q[0];
sx q[0];
rz(-2.1423253) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3827219) q[2];
sx q[2];
rz(-1.9577868) q[2];
sx q[2];
rz(0.92162161) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.091374) q[1];
sx q[1];
rz(-0.41793567) q[1];
sx q[1];
rz(1.0068847) q[1];
rz(2.897103) q[3];
sx q[3];
rz(-1.036662) q[3];
sx q[3];
rz(2.8492846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2723096) q[2];
sx q[2];
rz(-0.81988207) q[2];
sx q[2];
rz(0.38796866) q[2];
rz(-1.7913473) q[3];
sx q[3];
rz(-0.46348059) q[3];
sx q[3];
rz(2.8665682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824654) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(2.4097089) q[0];
rz(-0.16648079) q[1];
sx q[1];
rz(-2.6961168) q[1];
sx q[1];
rz(-3.0678715) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8548944) q[0];
sx q[0];
rz(-2.694283) q[0];
sx q[0];
rz(1.8064503) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6905626) q[2];
sx q[2];
rz(-1.6511859) q[2];
sx q[2];
rz(2.5188841) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2697061) q[1];
sx q[1];
rz(-1.6082113) q[1];
sx q[1];
rz(2.5579631) q[1];
rz(-pi) q[2];
rz(-2.5684486) q[3];
sx q[3];
rz(-1.4775652) q[3];
sx q[3];
rz(1.7236818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5471197) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(2.7129042) q[2];
rz(1.8244913) q[3];
sx q[3];
rz(-1.7842112) q[3];
sx q[3];
rz(-1.930442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4093032) q[0];
sx q[0];
rz(-1.4855054) q[0];
sx q[0];
rz(-1.0428585) q[0];
rz(1.5111142) q[1];
sx q[1];
rz(-1.7621367) q[1];
sx q[1];
rz(1.7932549) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52085984) q[0];
sx q[0];
rz(-3.0254786) q[0];
sx q[0];
rz(1.7945047) q[0];
rz(-pi) q[1];
rz(-2.3073763) q[2];
sx q[2];
rz(-1.5110821) q[2];
sx q[2];
rz(-0.63876736) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.41801449) q[1];
sx q[1];
rz(-1.2472767) q[1];
sx q[1];
rz(0.07467204) q[1];
x q[2];
rz(1.3799558) q[3];
sx q[3];
rz(-1.3340501) q[3];
sx q[3];
rz(-0.56425795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3412791) q[2];
sx q[2];
rz(-2.6276402) q[2];
sx q[2];
rz(0.094816118) q[2];
rz(1.173165) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(0.15869424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.512758) q[0];
sx q[0];
rz(-0.47369581) q[0];
sx q[0];
rz(1.0439903) q[0];
rz(-1.5402773) q[1];
sx q[1];
rz(-1.5834783) q[1];
sx q[1];
rz(-1.6316354) q[1];
rz(0.082967233) q[2];
sx q[2];
rz(-1.7832179) q[2];
sx q[2];
rz(-0.28945343) q[2];
rz(0.81090609) q[3];
sx q[3];
rz(-1.1699642) q[3];
sx q[3];
rz(-1.5311833) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

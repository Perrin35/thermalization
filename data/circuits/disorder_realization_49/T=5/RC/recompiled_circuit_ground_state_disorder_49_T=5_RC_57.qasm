OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9743118) q[0];
sx q[0];
rz(-0.33742961) q[0];
sx q[0];
rz(-2.9396102) q[0];
rz(2.1972411) q[1];
sx q[1];
rz(-0.75466067) q[1];
sx q[1];
rz(1.4442297) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9044735) q[0];
sx q[0];
rz(-2.5982862) q[0];
sx q[0];
rz(2.7271557) q[0];
rz(-pi) q[1];
rz(-3.0336122) q[2];
sx q[2];
rz(-1.7758581) q[2];
sx q[2];
rz(-0.16879665) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.18563823) q[1];
sx q[1];
rz(-2.1954838) q[1];
sx q[1];
rz(-2.5491979) q[1];
rz(-pi) q[2];
rz(-1.2060479) q[3];
sx q[3];
rz(-1.5284131) q[3];
sx q[3];
rz(-0.11755951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.28295383) q[2];
sx q[2];
rz(-0.084246548) q[2];
sx q[2];
rz(-1.6072744) q[2];
rz(2.9547847) q[3];
sx q[3];
rz(-2.6818633) q[3];
sx q[3];
rz(1.648858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9635791) q[0];
sx q[0];
rz(-0.78129617) q[0];
sx q[0];
rz(-0.70469967) q[0];
rz(-2.1251382) q[1];
sx q[1];
rz(-1.1926788) q[1];
sx q[1];
rz(1.9889529) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3160313) q[0];
sx q[0];
rz(-1.7411147) q[0];
sx q[0];
rz(2.1209149) q[0];
rz(-pi) q[1];
rz(-2.5637758) q[2];
sx q[2];
rz(-1.0989099) q[2];
sx q[2];
rz(-0.058171169) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1543442) q[1];
sx q[1];
rz(-2.1385405) q[1];
sx q[1];
rz(-2.5310234) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74781466) q[3];
sx q[3];
rz(-1.4749817) q[3];
sx q[3];
rz(-0.15296061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6013201) q[2];
sx q[2];
rz(-2.2728964) q[2];
sx q[2];
rz(-1.3199838) q[2];
rz(-3.0705304) q[3];
sx q[3];
rz(-3.0607405) q[3];
sx q[3];
rz(-2.0474056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0367301) q[0];
sx q[0];
rz(-2.7920089) q[0];
sx q[0];
rz(0.87376755) q[0];
rz(-1.8050487) q[1];
sx q[1];
rz(-0.68160325) q[1];
sx q[1];
rz(0.85493404) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23531518) q[0];
sx q[0];
rz(-1.7830669) q[0];
sx q[0];
rz(2.2215823) q[0];
rz(-2.1941691) q[2];
sx q[2];
rz(-0.50259841) q[2];
sx q[2];
rz(-1.1029152) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6715321) q[1];
sx q[1];
rz(-1.4619267) q[1];
sx q[1];
rz(-0.30372365) q[1];
x q[2];
rz(0.35936648) q[3];
sx q[3];
rz(-0.99726652) q[3];
sx q[3];
rz(1.9952967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9285589) q[2];
sx q[2];
rz(-1.1308257) q[2];
sx q[2];
rz(-2.0862759) q[2];
rz(-0.94163752) q[3];
sx q[3];
rz(-1.0712653) q[3];
sx q[3];
rz(-0.94282237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4554491) q[0];
sx q[0];
rz(-0.84857714) q[0];
sx q[0];
rz(2.4461179) q[0];
rz(-1.4675568) q[1];
sx q[1];
rz(-1.8753884) q[1];
sx q[1];
rz(0.094559018) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7692663) q[0];
sx q[0];
rz(-2.8594198) q[0];
sx q[0];
rz(1.7442987) q[0];
rz(-pi) q[1];
rz(2.2945061) q[2];
sx q[2];
rz(-1.0930335) q[2];
sx q[2];
rz(-1.0003045) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.27457224) q[1];
sx q[1];
rz(-2.7944075) q[1];
sx q[1];
rz(0.95496655) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26043268) q[3];
sx q[3];
rz(-2.2420001) q[3];
sx q[3];
rz(0.77653054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.62119421) q[2];
sx q[2];
rz(-0.79402557) q[2];
sx q[2];
rz(-2.3801079) q[2];
rz(0.17101184) q[3];
sx q[3];
rz(-1.2820425) q[3];
sx q[3];
rz(-3.0448992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8347725) q[0];
sx q[0];
rz(-2.4691395) q[0];
sx q[0];
rz(0.49224725) q[0];
rz(1.8222088) q[1];
sx q[1];
rz(-1.7183813) q[1];
sx q[1];
rz(2.6548903) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0224441) q[0];
sx q[0];
rz(-0.41765019) q[0];
sx q[0];
rz(0.29830427) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38551824) q[2];
sx q[2];
rz(-2.3340014) q[2];
sx q[2];
rz(2.453192) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4449205) q[1];
sx q[1];
rz(-1.2081523) q[1];
sx q[1];
rz(1.1128694) q[1];
x q[2];
rz(1.1718122) q[3];
sx q[3];
rz(-2.5907142) q[3];
sx q[3];
rz(-0.18607947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1897366) q[2];
sx q[2];
rz(-2.2978013) q[2];
sx q[2];
rz(-2.0437415) q[2];
rz(2.3073933) q[3];
sx q[3];
rz(-2.4662377) q[3];
sx q[3];
rz(1.7043017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51682353) q[0];
sx q[0];
rz(-0.56684816) q[0];
sx q[0];
rz(-0.045850642) q[0];
rz(2.397873) q[1];
sx q[1];
rz(-0.36879483) q[1];
sx q[1];
rz(-0.060294453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9550388) q[0];
sx q[0];
rz(-0.56772029) q[0];
sx q[0];
rz(2.4283152) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1072537) q[2];
sx q[2];
rz(-0.75159479) q[2];
sx q[2];
rz(1.3857402) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0452779) q[1];
sx q[1];
rz(-1.244925) q[1];
sx q[1];
rz(-2.9067944) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0457959) q[3];
sx q[3];
rz(-2.2471273) q[3];
sx q[3];
rz(2.6696882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9149949) q[2];
sx q[2];
rz(-1.8151585) q[2];
sx q[2];
rz(2.162852) q[2];
rz(1.8799051) q[3];
sx q[3];
rz(-1.0190957) q[3];
sx q[3];
rz(0.88920465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9642445) q[0];
sx q[0];
rz(-1.526399) q[0];
sx q[0];
rz(-2.2571795) q[0];
rz(-2.2824967) q[1];
sx q[1];
rz(-2.0180549) q[1];
sx q[1];
rz(-0.20800796) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73586845) q[0];
sx q[0];
rz(-0.62851539) q[0];
sx q[0];
rz(2.5121538) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9875097) q[2];
sx q[2];
rz(-2.1325958) q[2];
sx q[2];
rz(2.3903008) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0178725) q[1];
sx q[1];
rz(-1.8353288) q[1];
sx q[1];
rz(-2.2295647) q[1];
rz(1.5670942) q[3];
sx q[3];
rz(-0.50338689) q[3];
sx q[3];
rz(-2.9682669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38006833) q[2];
sx q[2];
rz(-0.91975776) q[2];
sx q[2];
rz(-2.5999787) q[2];
rz(-1.8875061) q[3];
sx q[3];
rz(-1.5897635) q[3];
sx q[3];
rz(0.90887535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33614531) q[0];
sx q[0];
rz(-0.31444028) q[0];
sx q[0];
rz(-1.0910777) q[0];
rz(-2.3167141) q[1];
sx q[1];
rz(-2.7179317) q[1];
sx q[1];
rz(1.0728015) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049041852) q[0];
sx q[0];
rz(-1.9922087) q[0];
sx q[0];
rz(-1.7803935) q[0];
rz(-1.0457951) q[2];
sx q[2];
rz(-3.1154251) q[2];
sx q[2];
rz(1.8963842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.271928) q[1];
sx q[1];
rz(-2.8638693) q[1];
sx q[1];
rz(0.50599338) q[1];
x q[2];
rz(1.1680702) q[3];
sx q[3];
rz(-1.793705) q[3];
sx q[3];
rz(0.18833489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23952809) q[2];
sx q[2];
rz(-0.54486474) q[2];
sx q[2];
rz(2.876335) q[2];
rz(2.9910763) q[3];
sx q[3];
rz(-2.033332) q[3];
sx q[3];
rz(-2.7786541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4758258) q[0];
sx q[0];
rz(-0.095223991) q[0];
sx q[0];
rz(1.4775403) q[0];
rz(-2.3382969) q[1];
sx q[1];
rz(-1.7684312) q[1];
sx q[1];
rz(-1.0172179) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4519891) q[0];
sx q[0];
rz(-1.7740806) q[0];
sx q[0];
rz(0.13138055) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3147589) q[2];
sx q[2];
rz(-2.0190329) q[2];
sx q[2];
rz(-2.2729276) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2443631) q[1];
sx q[1];
rz(-1.9112226) q[1];
sx q[1];
rz(-1.8805518) q[1];
rz(1.8735879) q[3];
sx q[3];
rz(-1.2583411) q[3];
sx q[3];
rz(-2.2175077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9251755) q[2];
sx q[2];
rz(-1.567013) q[2];
sx q[2];
rz(-2.4851921) q[2];
rz(0.21550719) q[3];
sx q[3];
rz(-1.7301205) q[3];
sx q[3];
rz(1.6925252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086275252) q[0];
sx q[0];
rz(-1.7973987) q[0];
sx q[0];
rz(-2.2917746) q[0];
rz(-0.96018106) q[1];
sx q[1];
rz(-0.4147059) q[1];
sx q[1];
rz(0.93920952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67341833) q[0];
sx q[0];
rz(-2.3395754) q[0];
sx q[0];
rz(-1.2002263) q[0];
x q[1];
rz(-3.1338336) q[2];
sx q[2];
rz(-1.9853704) q[2];
sx q[2];
rz(-1.3778292) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.997949) q[1];
sx q[1];
rz(-2.392329) q[1];
sx q[1];
rz(-0.56203385) q[1];
rz(-0.84374238) q[3];
sx q[3];
rz(-2.4795723) q[3];
sx q[3];
rz(0.5211422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5822997) q[2];
sx q[2];
rz(-1.4037932) q[2];
sx q[2];
rz(2.5913008) q[2];
rz(-0.127921) q[3];
sx q[3];
rz(-3.0030799) q[3];
sx q[3];
rz(0.96854717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.48548231) q[0];
sx q[0];
rz(-0.67018296) q[0];
sx q[0];
rz(-0.67076587) q[0];
rz(2.7863964) q[1];
sx q[1];
rz(-1.4757481) q[1];
sx q[1];
rz(-0.3955985) q[1];
rz(0.14329362) q[2];
sx q[2];
rz(-1.7096277) q[2];
sx q[2];
rz(-1.1451677) q[2];
rz(1.1200503) q[3];
sx q[3];
rz(-1.7786296) q[3];
sx q[3];
rz(2.8854388) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

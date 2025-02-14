OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.98915339) q[0];
sx q[0];
rz(-1.5664772) q[0];
sx q[0];
rz(-1.1251261) q[0];
rz(-5.2611051) q[1];
sx q[1];
rz(2.4740969) q[1];
sx q[1];
rz(11.752887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8047377) q[0];
sx q[0];
rz(-1.2053524) q[0];
sx q[0];
rz(2.9904537) q[0];
rz(-pi) q[1];
rz(-0.44806077) q[2];
sx q[2];
rz(-1.2431113) q[2];
sx q[2];
rz(-2.0304012) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2250925) q[1];
sx q[1];
rz(-1.719702) q[1];
sx q[1];
rz(2.6970106) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7343821) q[3];
sx q[3];
rz(-1.4814875) q[3];
sx q[3];
rz(2.341604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.4890613) q[2];
sx q[2];
rz(-1.6901313) q[2];
sx q[2];
rz(-2.2129464) q[2];
rz(-1.5422025) q[3];
sx q[3];
rz(-1.8079115) q[3];
sx q[3];
rz(-1.1716051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9376675) q[0];
sx q[0];
rz(-1.7610022) q[0];
sx q[0];
rz(-3.0145338) q[0];
rz(-0.98310414) q[1];
sx q[1];
rz(-1.3652722) q[1];
sx q[1];
rz(-0.7712706) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8255804) q[0];
sx q[0];
rz(-0.86992747) q[0];
sx q[0];
rz(-2.6722145) q[0];
rz(-pi) q[1];
rz(-3.0808582) q[2];
sx q[2];
rz(-1.9078476) q[2];
sx q[2];
rz(0.13523808) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.123536) q[1];
sx q[1];
rz(-0.50790826) q[1];
sx q[1];
rz(0.77202173) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5653789) q[3];
sx q[3];
rz(-2.0348843) q[3];
sx q[3];
rz(-0.6512385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1473006) q[2];
sx q[2];
rz(-1.2039801) q[2];
sx q[2];
rz(0.0083943923) q[2];
rz(-0.66347915) q[3];
sx q[3];
rz(-1.9050262) q[3];
sx q[3];
rz(2.8765163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0405149) q[0];
sx q[0];
rz(-0.82656693) q[0];
sx q[0];
rz(0.44152942) q[0];
rz(-2.1532374) q[1];
sx q[1];
rz(-2.007273) q[1];
sx q[1];
rz(-0.13557869) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7429333) q[0];
sx q[0];
rz(-3.1236095) q[0];
sx q[0];
rz(-2.6989486) q[0];
x q[1];
rz(-1.1889815) q[2];
sx q[2];
rz(-0.890632) q[2];
sx q[2];
rz(1.176468) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5244658) q[1];
sx q[1];
rz(-1.2167261) q[1];
sx q[1];
rz(-2.0984142) q[1];
rz(-pi) q[2];
rz(0.19702487) q[3];
sx q[3];
rz(-0.84614119) q[3];
sx q[3];
rz(-1.0850832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2077937) q[2];
sx q[2];
rz(-1.7094882) q[2];
sx q[2];
rz(-3.1033893) q[2];
rz(0.52538747) q[3];
sx q[3];
rz(-0.62756687) q[3];
sx q[3];
rz(0.6853404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.4811089) q[0];
sx q[0];
rz(-0.79367343) q[0];
sx q[0];
rz(2.5307181) q[0];
rz(1.3350217) q[1];
sx q[1];
rz(-1.7673312) q[1];
sx q[1];
rz(2.9023721) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0611899) q[0];
sx q[0];
rz(-2.1227269) q[0];
sx q[0];
rz(-0.30776382) q[0];
rz(-pi) q[1];
rz(2.6692713) q[2];
sx q[2];
rz(-1.9237674) q[2];
sx q[2];
rz(-1.8582839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.661631) q[1];
sx q[1];
rz(-1.1203655) q[1];
sx q[1];
rz(1.2110787) q[1];
x q[2];
rz(1.3705105) q[3];
sx q[3];
rz(-2.0242175) q[3];
sx q[3];
rz(0.1399006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7188344) q[2];
sx q[2];
rz(-1.6966635) q[2];
sx q[2];
rz(-3.139843) q[2];
rz(-2.8478029) q[3];
sx q[3];
rz(-1.9521451) q[3];
sx q[3];
rz(2.8747115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9002429) q[0];
sx q[0];
rz(-1.7647864) q[0];
sx q[0];
rz(-2.2985261) q[0];
rz(-0.33755606) q[1];
sx q[1];
rz(-1.866021) q[1];
sx q[1];
rz(-1.5379803) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46066901) q[0];
sx q[0];
rz(-2.3516555) q[0];
sx q[0];
rz(1.3000751) q[0];
x q[1];
rz(-1.5174687) q[2];
sx q[2];
rz(-1.3966718) q[2];
sx q[2];
rz(1.2542031) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.027315779) q[1];
sx q[1];
rz(-1.7661347) q[1];
sx q[1];
rz(2.8802425) q[1];
rz(-pi) q[2];
rz(-0.36253039) q[3];
sx q[3];
rz(-2.9177319) q[3];
sx q[3];
rz(2.8739342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7097077) q[2];
sx q[2];
rz(-1.0872492) q[2];
sx q[2];
rz(-2.0951927) q[2];
rz(0.31878582) q[3];
sx q[3];
rz(-0.71109486) q[3];
sx q[3];
rz(-0.54767245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043561291) q[0];
sx q[0];
rz(-1.8836319) q[0];
sx q[0];
rz(2.4936254) q[0];
rz(1.3994392) q[1];
sx q[1];
rz(-1.49767) q[1];
sx q[1];
rz(-1.8125777) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9861222) q[0];
sx q[0];
rz(-2.8948445) q[0];
sx q[0];
rz(2.0339436) q[0];
x q[1];
rz(0.29000303) q[2];
sx q[2];
rz(-1.776374) q[2];
sx q[2];
rz(-2.9221591) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15550286) q[1];
sx q[1];
rz(-0.8900607) q[1];
sx q[1];
rz(2.7377241) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59698128) q[3];
sx q[3];
rz(-2.5425445) q[3];
sx q[3];
rz(2.0839276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38137388) q[2];
sx q[2];
rz(-1.9332644) q[2];
sx q[2];
rz(-1.5488497) q[2];
rz(-0.069843944) q[3];
sx q[3];
rz(-1.8826238) q[3];
sx q[3];
rz(0.86863345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1182564) q[0];
sx q[0];
rz(-1.2160439) q[0];
sx q[0];
rz(2.1997531) q[0];
rz(0.16009227) q[1];
sx q[1];
rz(-1.6611049) q[1];
sx q[1];
rz(-2.9494185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3775786) q[0];
sx q[0];
rz(-1.6068216) q[0];
sx q[0];
rz(1.3467623) q[0];
rz(0.040080796) q[2];
sx q[2];
rz(-1.049384) q[2];
sx q[2];
rz(1.557795) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0380437) q[1];
sx q[1];
rz(-0.86665857) q[1];
sx q[1];
rz(1.1188933) q[1];
x q[2];
rz(-1.6113847) q[3];
sx q[3];
rz(-0.49792624) q[3];
sx q[3];
rz(-2.9663309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.75752246) q[2];
sx q[2];
rz(-1.9033868) q[2];
sx q[2];
rz(-0.43357098) q[2];
rz(-1.5571669) q[3];
sx q[3];
rz(-1.6470563) q[3];
sx q[3];
rz(0.43237329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50345355) q[0];
sx q[0];
rz(-1.764955) q[0];
sx q[0];
rz(-1.7973416) q[0];
rz(1.9380219) q[1];
sx q[1];
rz(-1.9529587) q[1];
sx q[1];
rz(-1.7787836) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8289611) q[0];
sx q[0];
rz(-1.5488708) q[0];
sx q[0];
rz(1.5410822) q[0];
x q[1];
rz(-2.3154738) q[2];
sx q[2];
rz(-1.4419793) q[2];
sx q[2];
rz(1.6507738) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4688022) q[1];
sx q[1];
rz(-1.4458477) q[1];
sx q[1];
rz(-0.23016696) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31764389) q[3];
sx q[3];
rz(-0.90932019) q[3];
sx q[3];
rz(0.7484439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5726996) q[2];
sx q[2];
rz(-0.77304825) q[2];
sx q[2];
rz(2.1577238) q[2];
rz(-0.24108663) q[3];
sx q[3];
rz(-2.2086996) q[3];
sx q[3];
rz(2.4510395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47138658) q[0];
sx q[0];
rz(-0.79600483) q[0];
sx q[0];
rz(0.27467003) q[0];
rz(1.341691) q[1];
sx q[1];
rz(-0.62961737) q[1];
sx q[1];
rz(1.0460269) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5513902) q[0];
sx q[0];
rz(-0.98197637) q[0];
sx q[0];
rz(1.9639652) q[0];
rz(-pi) q[1];
rz(1.7685031) q[2];
sx q[2];
rz(-1.6501763) q[2];
sx q[2];
rz(-2.6153836) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1362125) q[1];
sx q[1];
rz(-0.58717218) q[1];
sx q[1];
rz(3.0239952) q[1];
rz(-pi) q[2];
rz(1.4946497) q[3];
sx q[3];
rz(-0.5042432) q[3];
sx q[3];
rz(-1.7011736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2785953) q[2];
sx q[2];
rz(-1.7509165) q[2];
sx q[2];
rz(-2.7421303) q[2];
rz(0.56600371) q[3];
sx q[3];
rz(-2.1750906) q[3];
sx q[3];
rz(-0.81542265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.852916) q[0];
sx q[0];
rz(-1.4194019) q[0];
sx q[0];
rz(-2.5592819) q[0];
rz(0.18320006) q[1];
sx q[1];
rz(-2.2661426) q[1];
sx q[1];
rz(2.3351672) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2516625) q[0];
sx q[0];
rz(-1.4668873) q[0];
sx q[0];
rz(-3.0534239) q[0];
x q[1];
rz(2.8060032) q[2];
sx q[2];
rz(-1.3177455) q[2];
sx q[2];
rz(-2.0584681) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.90102531) q[1];
sx q[1];
rz(-1.7738924) q[1];
sx q[1];
rz(-2.3939449) q[1];
x q[2];
rz(-2.8424758) q[3];
sx q[3];
rz(-2.7795305) q[3];
sx q[3];
rz(0.55940926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9091984) q[2];
sx q[2];
rz(-0.71467233) q[2];
sx q[2];
rz(-0.8052899) q[2];
rz(0.14690873) q[3];
sx q[3];
rz(-0.78607905) q[3];
sx q[3];
rz(-2.4033191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.88937) q[0];
sx q[0];
rz(-1.3852373) q[0];
sx q[0];
rz(-1.2607384) q[0];
rz(-1.5421142) q[1];
sx q[1];
rz(-1.5972932) q[1];
sx q[1];
rz(1.6398026) q[1];
rz(-0.84757126) q[2];
sx q[2];
rz(-1.4847652) q[2];
sx q[2];
rz(2.8105856) q[2];
rz(2.4965548) q[3];
sx q[3];
rz(-1.8492263) q[3];
sx q[3];
rz(3.1254569) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1155137) q[0];
sx q[0];
rz(-1.4839412) q[0];
sx q[0];
rz(-0.32615647) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(1.3500554) q[1];
sx q[1];
rz(4.6842484) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.125995) q[0];
sx q[0];
rz(-1.7863569) q[0];
sx q[0];
rz(0.21685812) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3165663) q[2];
sx q[2];
rz(-2.1581274) q[2];
sx q[2];
rz(1.4884782) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3781158) q[1];
sx q[1];
rz(-0.33707481) q[1];
sx q[1];
rz(0.76428767) q[1];
rz(-pi) q[2];
rz(0.013734038) q[3];
sx q[3];
rz(-2.3547958) q[3];
sx q[3];
rz(-2.9247583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8618384) q[2];
sx q[2];
rz(-0.97936169) q[2];
sx q[2];
rz(-2.2564783) q[2];
rz(-2.4195813) q[3];
sx q[3];
rz(-1.4530028) q[3];
sx q[3];
rz(-0.0074145934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9279813) q[0];
sx q[0];
rz(-0.95887029) q[0];
sx q[0];
rz(-2.0425178) q[0];
rz(2.4765769) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(-0.87759334) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3646334) q[0];
sx q[0];
rz(-1.3898802) q[0];
sx q[0];
rz(-2.4448256) q[0];
rz(-pi) q[1];
rz(-0.34061956) q[2];
sx q[2];
rz(-2.7567731) q[2];
sx q[2];
rz(2.91586) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21786015) q[1];
sx q[1];
rz(-1.6777615) q[1];
sx q[1];
rz(-0.033543368) q[1];
rz(-0.86032805) q[3];
sx q[3];
rz(-2.6521157) q[3];
sx q[3];
rz(-1.5599172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.39891222) q[2];
sx q[2];
rz(-1.8360527) q[2];
sx q[2];
rz(1.1304643) q[2];
rz(1.2997262) q[3];
sx q[3];
rz(-1.2239417) q[3];
sx q[3];
rz(-1.4484423) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.995342) q[0];
sx q[0];
rz(-1.2820219) q[0];
sx q[0];
rz(0.28999844) q[0];
rz(0.6668123) q[1];
sx q[1];
rz(-2.1077483) q[1];
sx q[1];
rz(-0.07382948) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2436115) q[0];
sx q[0];
rz(-1.6323166) q[0];
sx q[0];
rz(0.089540066) q[0];
x q[1];
rz(-0.29088144) q[2];
sx q[2];
rz(-1.8996432) q[2];
sx q[2];
rz(1.8287303) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2245582) q[1];
sx q[1];
rz(-2.9552166) q[1];
sx q[1];
rz(1.264155) q[1];
rz(0.35145268) q[3];
sx q[3];
rz(-1.6822527) q[3];
sx q[3];
rz(-0.66672882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4905711) q[2];
sx q[2];
rz(-1.9475503) q[2];
sx q[2];
rz(0.64669615) q[2];
rz(-1.1086639) q[3];
sx q[3];
rz(-2.3587148) q[3];
sx q[3];
rz(1.1289319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9639503) q[0];
sx q[0];
rz(-2.9681866) q[0];
sx q[0];
rz(-1.1886764) q[0];
rz(2.1229318) q[1];
sx q[1];
rz(-0.97266346) q[1];
sx q[1];
rz(-1.4368988) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68543363) q[0];
sx q[0];
rz(-1.6628656) q[0];
sx q[0];
rz(0.7325367) q[0];
rz(-pi) q[1];
rz(-0.75886274) q[2];
sx q[2];
rz(-0.43735158) q[2];
sx q[2];
rz(-0.26668374) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1632299) q[1];
sx q[1];
rz(-2.0320738) q[1];
sx q[1];
rz(-0.39050885) q[1];
rz(-pi) q[2];
rz(2.339429) q[3];
sx q[3];
rz(-2.9608316) q[3];
sx q[3];
rz(0.67644955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.556095) q[2];
sx q[2];
rz(-1.4336339) q[2];
sx q[2];
rz(1.1222703) q[2];
rz(2.1155817) q[3];
sx q[3];
rz(-0.75338537) q[3];
sx q[3];
rz(2.1508353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.68200237) q[0];
sx q[0];
rz(-2.2408709) q[0];
sx q[0];
rz(-2.7686152) q[0];
rz(2.9176118) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(1.8251098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7861159) q[0];
sx q[0];
rz(-1.9516203) q[0];
sx q[0];
rz(-2.064408) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8226132) q[2];
sx q[2];
rz(-1.8251112) q[2];
sx q[2];
rz(-3.0248259) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5455294) q[1];
sx q[1];
rz(-1.5830333) q[1];
sx q[1];
rz(1.6227325) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54721197) q[3];
sx q[3];
rz(-1.1037165) q[3];
sx q[3];
rz(-0.78391778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0405154) q[2];
sx q[2];
rz(-0.83156362) q[2];
sx q[2];
rz(-0.80580795) q[2];
rz(-0.53330437) q[3];
sx q[3];
rz(-1.1321944) q[3];
sx q[3];
rz(1.0114975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39078113) q[0];
sx q[0];
rz(-1.823714) q[0];
sx q[0];
rz(-3.0506296) q[0];
rz(-0.85995752) q[1];
sx q[1];
rz(-1.1227612) q[1];
sx q[1];
rz(1.3202753) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0108171) q[0];
sx q[0];
rz(-2.5512619) q[0];
sx q[0];
rz(-0.70461313) q[0];
rz(-0.57029057) q[2];
sx q[2];
rz(-1.7208793) q[2];
sx q[2];
rz(-2.1652086) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49530416) q[1];
sx q[1];
rz(-1.0083645) q[1];
sx q[1];
rz(-0.56519392) q[1];
rz(-0.91464197) q[3];
sx q[3];
rz(-1.3387036) q[3];
sx q[3];
rz(2.7616012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0992574) q[2];
sx q[2];
rz(-2.1741185) q[2];
sx q[2];
rz(-0.60097224) q[2];
rz(-2.6565334) q[3];
sx q[3];
rz(-0.22189134) q[3];
sx q[3];
rz(1.4453567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8619974) q[0];
sx q[0];
rz(-1.9626564) q[0];
sx q[0];
rz(0.55554187) q[0];
rz(0.034596054) q[1];
sx q[1];
rz(-0.75841537) q[1];
sx q[1];
rz(1.7506036) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44156528) q[0];
sx q[0];
rz(-1.0374984) q[0];
sx q[0];
rz(-1.1839068) q[0];
rz(-pi) q[1];
rz(0.36632914) q[2];
sx q[2];
rz(-1.804367) q[2];
sx q[2];
rz(0.69586588) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6185551) q[1];
sx q[1];
rz(-0.68261519) q[1];
sx q[1];
rz(1.7726266) q[1];
rz(-pi) q[2];
x q[2];
rz(0.075918003) q[3];
sx q[3];
rz(-2.3644991) q[3];
sx q[3];
rz(-2.0457552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3283219) q[2];
sx q[2];
rz(-2.6997824) q[2];
sx q[2];
rz(-2.7461046) q[2];
rz(-1.8528806) q[3];
sx q[3];
rz(-1.6059395) q[3];
sx q[3];
rz(0.66974631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.678858) q[0];
sx q[0];
rz(-2.8023219) q[0];
sx q[0];
rz(-1.6495552) q[0];
rz(2.18816) q[1];
sx q[1];
rz(-1.1089193) q[1];
sx q[1];
rz(-1.7038201) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8673082) q[0];
sx q[0];
rz(-2.0831997) q[0];
sx q[0];
rz(-2.9318277) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4687924) q[2];
sx q[2];
rz(-2.9252508) q[2];
sx q[2];
rz(2.1397482) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7319665) q[1];
sx q[1];
rz(-2.1831174) q[1];
sx q[1];
rz(2.1698976) q[1];
rz(3.1215454) q[3];
sx q[3];
rz(-1.123395) q[3];
sx q[3];
rz(0.1964387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4961204) q[2];
sx q[2];
rz(-0.43473736) q[2];
sx q[2];
rz(0.17871857) q[2];
rz(2.2802165) q[3];
sx q[3];
rz(-1.9390315) q[3];
sx q[3];
rz(-2.7698959) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9361967) q[0];
sx q[0];
rz(-1.9367171) q[0];
sx q[0];
rz(-2.0478915) q[0];
rz(0.73668346) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(-1.0587943) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5537542) q[0];
sx q[0];
rz(-0.61843473) q[0];
sx q[0];
rz(2.0536325) q[0];
rz(-pi) q[1];
rz(2.2718614) q[2];
sx q[2];
rz(-2.5517002) q[2];
sx q[2];
rz(1.215495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.141507) q[1];
sx q[1];
rz(-0.39751378) q[1];
sx q[1];
rz(1.2500398) q[1];
x q[2];
rz(2.4478854) q[3];
sx q[3];
rz(-2.5543164) q[3];
sx q[3];
rz(1.3035843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.27292353) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(1.6607364) q[2];
rz(-2.7311834) q[3];
sx q[3];
rz(-1.4199665) q[3];
sx q[3];
rz(2.9836392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2720298) q[0];
sx q[0];
rz(-0.27150387) q[0];
sx q[0];
rz(2.8503382) q[0];
rz(0.60925305) q[1];
sx q[1];
rz(-1.4657425) q[1];
sx q[1];
rz(1.4321009) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8490484) q[0];
sx q[0];
rz(-1.8869072) q[0];
sx q[0];
rz(-2.6228117) q[0];
x q[1];
rz(0.45639313) q[2];
sx q[2];
rz(-0.49955873) q[2];
sx q[2];
rz(1.511614) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4869625) q[1];
sx q[1];
rz(-1.6143867) q[1];
sx q[1];
rz(2.589588) q[1];
x q[2];
rz(-1.1389144) q[3];
sx q[3];
rz(-1.6500041) q[3];
sx q[3];
rz(-1.6457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46135205) q[2];
sx q[2];
rz(-1.1169008) q[2];
sx q[2];
rz(-2.0001901) q[2];
rz(-1.5348148) q[3];
sx q[3];
rz(-1.1780058) q[3];
sx q[3];
rz(-2.6721568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5466945) q[0];
sx q[0];
rz(-1.5190769) q[0];
sx q[0];
rz(1.4357823) q[0];
rz(0.36874157) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(-0.22696162) q[2];
sx q[2];
rz(-1.3703176) q[2];
sx q[2];
rz(-0.3581518) q[2];
rz(-2.7183919) q[3];
sx q[3];
rz(-1.7508218) q[3];
sx q[3];
rz(-1.4765061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

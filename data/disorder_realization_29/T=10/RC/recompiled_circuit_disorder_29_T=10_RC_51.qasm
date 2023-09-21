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
rz(2.8154362) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(-1.7915373) q[1];
sx q[1];
rz(-1.5426558) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.125995) q[0];
sx q[0];
rz(-1.7863569) q[0];
sx q[0];
rz(-0.21685812) q[0];
x q[1];
rz(-1.3165663) q[2];
sx q[2];
rz(-0.98346522) q[2];
sx q[2];
rz(-1.6531144) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.029763873) q[1];
sx q[1];
rz(-1.811869) q[1];
sx q[1];
rz(-1.8087216) q[1];
rz(-1.5845675) q[3];
sx q[3];
rz(-2.3574986) q[3];
sx q[3];
rz(0.23628326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8618384) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(2.2564783) q[2];
rz(-2.4195813) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(-3.1341781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2136114) q[0];
sx q[0];
rz(-2.1827224) q[0];
sx q[0];
rz(-1.0990748) q[0];
rz(0.66501578) q[1];
sx q[1];
rz(-1.4140833) q[1];
sx q[1];
rz(-0.87759334) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7769593) q[0];
sx q[0];
rz(-1.7517125) q[0];
sx q[0];
rz(0.69676708) q[0];
rz(-pi) q[1];
rz(2.8009731) q[2];
sx q[2];
rz(-2.7567731) q[2];
sx q[2];
rz(-0.22573267) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6192012) q[1];
sx q[1];
rz(-3.0295105) q[1];
sx q[1];
rz(1.8735318) q[1];
x q[2];
rz(1.1869933) q[3];
sx q[3];
rz(-1.882453) q[3];
sx q[3];
rz(-0.66031885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7426804) q[2];
sx q[2];
rz(-1.8360527) q[2];
sx q[2];
rz(-2.0111283) q[2];
rz(-1.2997262) q[3];
sx q[3];
rz(-1.9176509) q[3];
sx q[3];
rz(1.6931504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.995342) q[0];
sx q[0];
rz(-1.2820219) q[0];
sx q[0];
rz(2.8515942) q[0];
rz(-2.4747804) q[1];
sx q[1];
rz(-2.1077483) q[1];
sx q[1];
rz(3.0677632) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66729546) q[0];
sx q[0];
rz(-1.4814261) q[0];
sx q[0];
rz(-1.6325634) q[0];
rz(-1.2286085) q[2];
sx q[2];
rz(-1.8456736) q[2];
sx q[2];
rz(0.16155044) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2245582) q[1];
sx q[1];
rz(-2.9552166) q[1];
sx q[1];
rz(-1.264155) q[1];
rz(-1.6894433) q[3];
sx q[3];
rz(-1.9199748) q[3];
sx q[3];
rz(-2.196764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4905711) q[2];
sx q[2];
rz(-1.1940424) q[2];
sx q[2];
rz(2.4948965) q[2];
rz(-1.1086639) q[3];
sx q[3];
rz(-2.3587148) q[3];
sx q[3];
rz(-2.0126608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9639503) q[0];
sx q[0];
rz(-0.17340604) q[0];
sx q[0];
rz(1.9529163) q[0];
rz(-2.1229318) q[1];
sx q[1];
rz(-2.1689292) q[1];
sx q[1];
rz(-1.4368988) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68543363) q[0];
sx q[0];
rz(-1.6628656) q[0];
sx q[0];
rz(-2.409056) q[0];
rz(2.3827299) q[2];
sx q[2];
rz(-0.43735158) q[2];
sx q[2];
rz(2.8749089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1632299) q[1];
sx q[1];
rz(-2.0320738) q[1];
sx q[1];
rz(2.7510838) q[1];
rz(-0.12636633) q[3];
sx q[3];
rz(-1.7003945) q[3];
sx q[3];
rz(-1.6882997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.556095) q[2];
sx q[2];
rz(-1.4336339) q[2];
sx q[2];
rz(-2.0193224) q[2];
rz(-1.026011) q[3];
sx q[3];
rz(-2.3882073) q[3];
sx q[3];
rz(-2.1508353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68200237) q[0];
sx q[0];
rz(-0.90072173) q[0];
sx q[0];
rz(2.7686152) q[0];
rz(0.22398082) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(-1.8251098) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017942863) q[0];
sx q[0];
rz(-2.0262449) q[0];
sx q[0];
rz(0.42670593) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3189795) q[2];
sx q[2];
rz(-1.8251112) q[2];
sx q[2];
rz(0.11676678) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5455294) q[1];
sx q[1];
rz(-1.5830333) q[1];
sx q[1];
rz(-1.6227325) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0373711) q[3];
sx q[3];
rz(-1.08764) q[3];
sx q[3];
rz(2.0869568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.10107723) q[2];
sx q[2];
rz(-0.83156362) q[2];
sx q[2];
rz(0.80580795) q[2];
rz(-2.6082883) q[3];
sx q[3];
rz(-1.1321944) q[3];
sx q[3];
rz(2.1300952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39078113) q[0];
sx q[0];
rz(-1.823714) q[0];
sx q[0];
rz(-0.090963013) q[0];
rz(2.2816351) q[1];
sx q[1];
rz(-1.1227612) q[1];
sx q[1];
rz(1.3202753) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1307756) q[0];
sx q[0];
rz(-2.5512619) q[0];
sx q[0];
rz(0.70461313) q[0];
rz(1.3930416) q[2];
sx q[2];
rz(-1.0076992) q[2];
sx q[2];
rz(2.4515738) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7399866) q[1];
sx q[1];
rz(-2.0410129) q[1];
sx q[1];
rz(-2.2120038) q[1];
rz(1.201186) q[3];
sx q[3];
rz(-0.69023057) q[3];
sx q[3];
rz(-1.481101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0992574) q[2];
sx q[2];
rz(-2.1741185) q[2];
sx q[2];
rz(-2.5406204) q[2];
rz(2.6565334) q[3];
sx q[3];
rz(-0.22189134) q[3];
sx q[3];
rz(-1.4453567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8619974) q[0];
sx q[0];
rz(-1.1789362) q[0];
sx q[0];
rz(-2.5860508) q[0];
rz(-0.034596054) q[1];
sx q[1];
rz(-2.3831773) q[1];
sx q[1];
rz(1.7506036) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9075508) q[0];
sx q[0];
rz(-0.64767917) q[0];
sx q[0];
rz(2.5729022) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5552942) q[2];
sx q[2];
rz(-0.43160298) q[2];
sx q[2];
rz(1.4177711) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2512867) q[1];
sx q[1];
rz(-1.444) q[1];
sx q[1];
rz(2.2433953) q[1];
x q[2];
rz(-1.6452541) q[3];
sx q[3];
rz(-2.3450608) q[3];
sx q[3];
rz(-1.9394685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3283219) q[2];
sx q[2];
rz(-0.44181028) q[2];
sx q[2];
rz(0.39548809) q[2];
rz(-1.8528806) q[3];
sx q[3];
rz(-1.5356531) q[3];
sx q[3];
rz(-0.66974631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.678858) q[0];
sx q[0];
rz(-0.33927074) q[0];
sx q[0];
rz(1.6495552) q[0];
rz(2.18816) q[1];
sx q[1];
rz(-1.1089193) q[1];
sx q[1];
rz(1.4377726) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8673082) q[0];
sx q[0];
rz(-1.058393) q[0];
sx q[0];
rz(-2.9318277) q[0];
rz(-pi) q[1];
rz(2.9713694) q[2];
sx q[2];
rz(-1.7049689) q[2];
sx q[2];
rz(1.9113049) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9287712) q[1];
sx q[1];
rz(-2.050424) q[1];
sx q[1];
rz(2.4368083) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5290456) q[3];
sx q[3];
rz(-0.44781993) q[3];
sx q[3];
rz(-2.8988422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.64547223) q[2];
sx q[2];
rz(-0.43473736) q[2];
sx q[2];
rz(0.17871857) q[2];
rz(-0.86137613) q[3];
sx q[3];
rz(-1.9390315) q[3];
sx q[3];
rz(-2.7698959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9361967) q[0];
sx q[0];
rz(-1.2048756) q[0];
sx q[0];
rz(2.0478915) q[0];
rz(2.4049092) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(-2.0827983) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5878384) q[0];
sx q[0];
rz(-2.5231579) q[0];
sx q[0];
rz(-2.0536325) q[0];
rz(-pi) q[1];
rz(2.043622) q[2];
sx q[2];
rz(-1.2038004) q[2];
sx q[2];
rz(-0.25640139) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4150548) q[1];
sx q[1];
rz(-1.6931567) q[1];
sx q[1];
rz(1.9499669) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1684253) q[3];
sx q[3];
rz(-2.0109004) q[3];
sx q[3];
rz(2.6228867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8686691) q[2];
sx q[2];
rz(-0.69512853) q[2];
sx q[2];
rz(1.4808562) q[2];
rz(-0.41040928) q[3];
sx q[3];
rz(-1.7216262) q[3];
sx q[3];
rz(-0.15795344) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8695628) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(2.8503382) q[0];
rz(-2.5323396) q[1];
sx q[1];
rz(-1.4657425) q[1];
sx q[1];
rz(-1.7094918) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0390022) q[0];
sx q[0];
rz(-2.0615091) q[0];
sx q[0];
rz(-1.9309994) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6860793) q[2];
sx q[2];
rz(-1.7835155) q[2];
sx q[2];
rz(0.46609512) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1283975) q[1];
sx q[1];
rz(-2.5880475) q[1];
sx q[1];
rz(-3.058606) q[1];
rz(-pi) q[2];
rz(2.0026783) q[3];
sx q[3];
rz(-1.4915885) q[3];
sx q[3];
rz(1.6457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6802406) q[2];
sx q[2];
rz(-1.1169008) q[2];
sx q[2];
rz(2.0001901) q[2];
rz(1.5348148) q[3];
sx q[3];
rz(-1.9635868) q[3];
sx q[3];
rz(0.46943584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5948982) q[0];
sx q[0];
rz(-1.5190769) q[0];
sx q[0];
rz(1.4357823) q[0];
rz(0.36874157) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(-0.73451191) q[2];
sx q[2];
rz(-0.30167689) q[2];
sx q[2];
rz(-2.6405356) q[2];
rz(2.7183919) q[3];
sx q[3];
rz(-1.3907708) q[3];
sx q[3];
rz(1.6650865) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

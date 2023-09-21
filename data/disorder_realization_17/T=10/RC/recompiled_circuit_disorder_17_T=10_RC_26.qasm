OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9085812) q[0];
sx q[0];
rz(-1.9549978) q[0];
sx q[0];
rz(-1.1458122) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(2.8491128) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777893) q[0];
sx q[0];
rz(-1.6262494) q[0];
sx q[0];
rz(1.9652912) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77941676) q[2];
sx q[2];
rz(-1.5343108) q[2];
sx q[2];
rz(-1.905029) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9521192) q[1];
sx q[1];
rz(-1.5249624) q[1];
sx q[1];
rz(-2.1650251) q[1];
rz(-pi) q[2];
rz(2.0446159) q[3];
sx q[3];
rz(-1.5690656) q[3];
sx q[3];
rz(-0.8918106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73074377) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(2.6921819) q[2];
rz(0.64569008) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(-1.8507563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21355024) q[0];
sx q[0];
rz(-2.2556861) q[0];
sx q[0];
rz(-2.399562) q[0];
rz(-1.6702601) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(-1.7785633) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.610696) q[0];
sx q[0];
rz(-0.96141059) q[0];
sx q[0];
rz(1.0884398) q[0];
rz(1.3554553) q[2];
sx q[2];
rz(-2.0305579) q[2];
sx q[2];
rz(2.721867) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.14582536) q[1];
sx q[1];
rz(-1.9103423) q[1];
sx q[1];
rz(1.550436) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19248776) q[3];
sx q[3];
rz(-1.3093595) q[3];
sx q[3];
rz(2.1646433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.30119511) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(1.8691241) q[2];
rz(-2.8043591) q[3];
sx q[3];
rz(-1.1703346) q[3];
sx q[3];
rz(-3.1315394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067588016) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(2.8787676) q[0];
rz(-0.35573959) q[1];
sx q[1];
rz(-1.6825312) q[1];
sx q[1];
rz(1.2353108) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6806157) q[0];
sx q[0];
rz(-1.5810228) q[0];
sx q[0];
rz(-0.015120487) q[0];
x q[1];
rz(1.6791061) q[2];
sx q[2];
rz(-0.24178594) q[2];
sx q[2];
rz(2.0568648) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23535145) q[1];
sx q[1];
rz(-1.6760577) q[1];
sx q[1];
rz(2.8180442) q[1];
rz(-pi) q[2];
rz(-1.1615109) q[3];
sx q[3];
rz(-0.73771362) q[3];
sx q[3];
rz(1.1324594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1245023) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(0.18033218) q[2];
rz(-2.5056433) q[3];
sx q[3];
rz(-1.9790117) q[3];
sx q[3];
rz(-0.37160555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76698774) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(2.4225127) q[0];
rz(2.6285697) q[1];
sx q[1];
rz(-1.2027556) q[1];
sx q[1];
rz(-1.0346574) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7906404) q[0];
sx q[0];
rz(-1.5119023) q[0];
sx q[0];
rz(1.590593) q[0];
rz(-pi) q[1];
rz(3.0940975) q[2];
sx q[2];
rz(-1.4357899) q[2];
sx q[2];
rz(0.53802711) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50463488) q[1];
sx q[1];
rz(-2.2491498) q[1];
sx q[1];
rz(-1.9370609) q[1];
rz(-2.6684768) q[3];
sx q[3];
rz(-1.72662) q[3];
sx q[3];
rz(-1.5161878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2085312) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(-2.9053524) q[2];
rz(1.158372) q[3];
sx q[3];
rz(-2.1393496) q[3];
sx q[3];
rz(2.2896144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1119969) q[0];
sx q[0];
rz(-1.010226) q[0];
sx q[0];
rz(-0.90233666) q[0];
rz(-0.92102712) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(2.5193118) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8392004) q[0];
sx q[0];
rz(-1.3883739) q[0];
sx q[0];
rz(-1.8586041) q[0];
rz(-pi) q[1];
rz(-0.12139608) q[2];
sx q[2];
rz(-0.65181323) q[2];
sx q[2];
rz(2.7898942) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8473377) q[1];
sx q[1];
rz(-1.6604742) q[1];
sx q[1];
rz(1.4011821) q[1];
rz(-2.6617674) q[3];
sx q[3];
rz(-1.92355) q[3];
sx q[3];
rz(1.477581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2174012) q[2];
sx q[2];
rz(-1.8175586) q[2];
sx q[2];
rz(0.021281555) q[2];
rz(1.6993258) q[3];
sx q[3];
rz(-2.4271624) q[3];
sx q[3];
rz(-0.91059476) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4845881) q[0];
sx q[0];
rz(-0.5963043) q[0];
sx q[0];
rz(3.0969627) q[0];
rz(2.1394829) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(-0.034428509) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39715595) q[0];
sx q[0];
rz(-1.7560871) q[0];
sx q[0];
rz(0.013884355) q[0];
rz(-2.0738515) q[2];
sx q[2];
rz(-1.6776553) q[2];
sx q[2];
rz(1.6018794) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.039636314) q[1];
sx q[1];
rz(-2.4974681) q[1];
sx q[1];
rz(2.1307751) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9488726) q[3];
sx q[3];
rz(-2.4270714) q[3];
sx q[3];
rz(0.89380985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3520711) q[2];
sx q[2];
rz(-1.8692724) q[2];
sx q[2];
rz(-2.1748523) q[2];
rz(0.012185193) q[3];
sx q[3];
rz(-0.92958486) q[3];
sx q[3];
rz(3.0325082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1352585) q[0];
sx q[0];
rz(-0.2671347) q[0];
sx q[0];
rz(-2.1475041) q[0];
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
rz(1.3294322) q[0];
sx q[0];
rz(-1.433916) q[0];
sx q[0];
rz(-3.0457004) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5653966) q[2];
sx q[2];
rz(-1.155876) q[2];
sx q[2];
rz(2.8232099) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.48078295) q[1];
sx q[1];
rz(-1.0479095) q[1];
sx q[1];
rz(1.3516264) q[1];
x q[2];
rz(2.5218511) q[3];
sx q[3];
rz(-1.3337787) q[3];
sx q[3];
rz(1.851351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.58549515) q[2];
sx q[2];
rz(-1.408351) q[2];
sx q[2];
rz(0.56376702) q[2];
rz(-2.901315) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(-1.3747922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4450842) q[0];
sx q[0];
rz(-2.962528) q[0];
sx q[0];
rz(2.5575496) q[0];
rz(2.7208327) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(-0.27841321) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70008343) q[0];
sx q[0];
rz(-2.1372876) q[0];
sx q[0];
rz(2.9928815) q[0];
rz(-pi) q[1];
rz(1.7588708) q[2];
sx q[2];
rz(-1.9577868) q[2];
sx q[2];
rz(-0.92162161) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1380996) q[1];
sx q[1];
rz(-1.7894735) q[1];
sx q[1];
rz(1.9298645) q[1];
rz(-pi) q[2];
rz(2.897103) q[3];
sx q[3];
rz(-1.036662) q[3];
sx q[3];
rz(-0.29230803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.86928308) q[2];
sx q[2];
rz(-0.81988207) q[2];
sx q[2];
rz(2.753624) q[2];
rz(-1.7913473) q[3];
sx q[3];
rz(-2.6781121) q[3];
sx q[3];
rz(-2.8665682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824654) q[0];
sx q[0];
rz(-2.9282741) q[0];
sx q[0];
rz(0.7318837) q[0];
rz(2.9751119) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(3.0678715) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28669824) q[0];
sx q[0];
rz(-2.694283) q[0];
sx q[0];
rz(1.8064503) q[0];
x q[1];
rz(1.45103) q[2];
sx q[2];
rz(-1.4904067) q[2];
sx q[2];
rz(0.62270852) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8970866) q[1];
sx q[1];
rz(-2.5569041) q[1];
sx q[1];
rz(3.0737682) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1707465) q[3];
sx q[3];
rz(-2.5617544) q[3];
sx q[3];
rz(0.29614007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5471197) q[2];
sx q[2];
rz(-2.9056845) q[2];
sx q[2];
rz(-2.7129042) q[2];
rz(1.8244913) q[3];
sx q[3];
rz(-1.3573815) q[3];
sx q[3];
rz(-1.2111506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73228943) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(1.0428585) q[0];
rz(1.5111142) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(-1.7932549) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6207328) q[0];
sx q[0];
rz(-0.11611406) q[0];
sx q[0];
rz(1.3470879) q[0];
rz(-pi) q[1];
rz(2.3073763) q[2];
sx q[2];
rz(-1.5110821) q[2];
sx q[2];
rz(-2.5028253) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.18689449) q[1];
sx q[1];
rz(-2.8098626) q[1];
sx q[1];
rz(1.7897254) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24095778) q[3];
sx q[3];
rz(-1.7562508) q[3];
sx q[3];
rz(-2.0897739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.80031359) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(0.094816118) q[2];
rz(1.173165) q[3];
sx q[3];
rz(-1.4011819) q[3];
sx q[3];
rz(2.9828984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(1.3576635) q[2];
sx q[2];
rz(-1.6518946) q[2];
sx q[2];
rz(1.2988731) q[2];
rz(-2.6125828) q[3];
sx q[3];
rz(-0.88376868) q[3];
sx q[3];
rz(0.39467011) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
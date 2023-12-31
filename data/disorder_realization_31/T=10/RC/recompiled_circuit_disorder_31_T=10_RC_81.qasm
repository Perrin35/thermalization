OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(0.84258643) q[0];
rz(-1.0358345) q[1];
sx q[1];
rz(-2.0422715) q[1];
sx q[1];
rz(1.6834747) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5320839) q[0];
sx q[0];
rz(-2.3713787) q[0];
sx q[0];
rz(0.30755933) q[0];
rz(1.5904434) q[2];
sx q[2];
rz(-2.1845449) q[2];
sx q[2];
rz(0.27054271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.018651389) q[1];
sx q[1];
rz(-2.7416347) q[1];
sx q[1];
rz(-2.8040228) q[1];
rz(-pi) q[2];
rz(1.6928715) q[3];
sx q[3];
rz(-2.1616518) q[3];
sx q[3];
rz(-1.6216244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6136916) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(0.17949417) q[2];
rz(1.2256631) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.74801385) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(-2.8161312) q[0];
rz(1.7851967) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(-1.1546086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5529454) q[0];
sx q[0];
rz(-1.5536904) q[0];
sx q[0];
rz(3.1228035) q[0];
x q[1];
rz(-2.1968368) q[2];
sx q[2];
rz(-1.895004) q[2];
sx q[2];
rz(-2.7446483) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0305811) q[1];
sx q[1];
rz(-1.4861373) q[1];
sx q[1];
rz(-0.76534033) q[1];
x q[2];
rz(-1.8335908) q[3];
sx q[3];
rz(-1.3920708) q[3];
sx q[3];
rz(-0.59584111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4521728) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(-0.88341218) q[2];
rz(2.6702821) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31323355) q[0];
sx q[0];
rz(-1.4947083) q[0];
sx q[0];
rz(-1.6261684) q[0];
rz(2.5405163) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(-2.0498958) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1356782) q[0];
sx q[0];
rz(-2.1164829) q[0];
sx q[0];
rz(2.2360327) q[0];
rz(-pi) q[1];
rz(0.41468427) q[2];
sx q[2];
rz(-0.95699246) q[2];
sx q[2];
rz(2.2874122) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8852383) q[1];
sx q[1];
rz(-1.2043722) q[1];
sx q[1];
rz(0.79856915) q[1];
rz(2.6767119) q[3];
sx q[3];
rz(-2.998623) q[3];
sx q[3];
rz(-2.3425992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(0.88095218) q[2];
rz(1.3736003) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(-1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3110733) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(2.7048892) q[0];
rz(-0.23315915) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(2.8312347) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4219907) q[0];
sx q[0];
rz(-1.4388226) q[0];
sx q[0];
rz(-0.3814401) q[0];
x q[1];
rz(2.9878346) q[2];
sx q[2];
rz(-2.4506844) q[2];
sx q[2];
rz(-1.549987) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51018184) q[1];
sx q[1];
rz(-1.5522172) q[1];
sx q[1];
rz(-1.2127962) q[1];
x q[2];
rz(-0.035590812) q[3];
sx q[3];
rz(-1.6944052) q[3];
sx q[3];
rz(1.3328758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0115396) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(2.0641573) q[2];
rz(3.0854026) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(-1.594054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.189165) q[0];
sx q[0];
rz(-2.0963033) q[0];
sx q[0];
rz(-0.24965832) q[0];
rz(1.5646308) q[1];
sx q[1];
rz(-2.3639634) q[1];
sx q[1];
rz(0.87019428) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10660431) q[0];
sx q[0];
rz(-2.7634794) q[0];
sx q[0];
rz(2.5614221) q[0];
rz(-pi) q[1];
rz(1.382803) q[2];
sx q[2];
rz(-2.2586939) q[2];
sx q[2];
rz(-0.57304136) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1340027) q[1];
sx q[1];
rz(-0.78467272) q[1];
sx q[1];
rz(0.44962928) q[1];
rz(-pi) q[2];
rz(-1.8364041) q[3];
sx q[3];
rz(-2.5634273) q[3];
sx q[3];
rz(0.8997013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8683118) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(2.4678521) q[2];
rz(-0.30361787) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7917787) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(-0.2579903) q[0];
rz(-2.7164283) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.4917096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30116044) q[0];
sx q[0];
rz(-1.668881) q[0];
sx q[0];
rz(0.43099404) q[0];
rz(-0.67955534) q[2];
sx q[2];
rz(-1.9050042) q[2];
sx q[2];
rz(-2.213775) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7219833) q[1];
sx q[1];
rz(-2.1188542) q[1];
sx q[1];
rz(2.4649058) q[1];
x q[2];
rz(1.2061938) q[3];
sx q[3];
rz(-1.4758665) q[3];
sx q[3];
rz(-1.7942384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.012718) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(3.0498665) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(2.2475524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82350746) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(0.41123018) q[0];
rz(0.86589083) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(0.033989865) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995178) q[0];
sx q[0];
rz(-2.3410428) q[0];
sx q[0];
rz(-2.9151025) q[0];
rz(-0.54079536) q[2];
sx q[2];
rz(-2.135709) q[2];
sx q[2];
rz(-0.62513798) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.36044869) q[1];
sx q[1];
rz(-0.89372674) q[1];
sx q[1];
rz(3.0767246) q[1];
rz(-pi) q[2];
rz(-1.4230698) q[3];
sx q[3];
rz(-0.96102321) q[3];
sx q[3];
rz(-2.8019398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6035446) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(-0.87654385) q[2];
rz(2.792568) q[3];
sx q[3];
rz(-1.9411496) q[3];
sx q[3];
rz(-2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.2440764) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(-0.39392719) q[0];
rz(-0.36755964) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(-1.4454909) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1720393) q[0];
sx q[0];
rz(-1.5656099) q[0];
sx q[0];
rz(1.1374377) q[0];
rz(-pi) q[1];
rz(-2.5567899) q[2];
sx q[2];
rz(-1.0324761) q[2];
sx q[2];
rz(0.98758299) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1576924) q[1];
sx q[1];
rz(-1.0219814) q[1];
sx q[1];
rz(1.7556612) q[1];
rz(-pi) q[2];
rz(0.74650713) q[3];
sx q[3];
rz(-2.3254447) q[3];
sx q[3];
rz(1.5619123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8470856) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(0.40714804) q[2];
rz(1.6242705) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.80609926) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(-1.2517713) q[0];
rz(-2.4720526) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(-2.8318185) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1447434) q[0];
sx q[0];
rz(-2.092917) q[0];
sx q[0];
rz(-1.7168619) q[0];
rz(1.3938815) q[2];
sx q[2];
rz(-1.5614911) q[2];
sx q[2];
rz(2.6975346) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5619547) q[1];
sx q[1];
rz(-1.5127752) q[1];
sx q[1];
rz(1.3647563) q[1];
rz(-pi) q[2];
rz(-2.7731032) q[3];
sx q[3];
rz(-0.62928761) q[3];
sx q[3];
rz(3.1143509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4391675) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(-1.9343728) q[2];
rz(1.0369982) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(-0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0913775) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(-1.2058831) q[0];
rz(0.58569113) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(1.6419798) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3033894) q[0];
sx q[0];
rz(-1.8614385) q[0];
sx q[0];
rz(-3.035726) q[0];
rz(3.0707804) q[2];
sx q[2];
rz(-2.376308) q[2];
sx q[2];
rz(0.27809696) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8780898) q[1];
sx q[1];
rz(-0.51551688) q[1];
sx q[1];
rz(1.131119) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5851192) q[3];
sx q[3];
rz(-2.9687772) q[3];
sx q[3];
rz(-0.9848136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2575834) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(0.94669) q[2];
rz(0.36869129) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(-2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7832227) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(-3.070667) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(0.46438607) q[2];
sx q[2];
rz(-1.1190363) q[2];
sx q[2];
rz(-2.6418532) q[2];
rz(2.5881913) q[3];
sx q[3];
rz(-0.80080606) q[3];
sx q[3];
rz(-1.3047119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

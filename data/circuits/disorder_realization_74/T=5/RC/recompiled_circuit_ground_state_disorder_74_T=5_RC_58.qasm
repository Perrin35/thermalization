OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5689019) q[0];
sx q[0];
rz(-0.98839086) q[0];
sx q[0];
rz(2.6925777) q[0];
rz(-0.16945893) q[1];
sx q[1];
rz(-3.0100477) q[1];
sx q[1];
rz(1.1313255) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5558452) q[0];
sx q[0];
rz(-2.1679584) q[0];
sx q[0];
rz(0.35564977) q[0];
rz(-0.65931084) q[2];
sx q[2];
rz(-2.3917195) q[2];
sx q[2];
rz(-1.4091968) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1234963) q[1];
sx q[1];
rz(-0.61401788) q[1];
sx q[1];
rz(-1.2662925) q[1];
rz(-pi) q[2];
rz(1.2901957) q[3];
sx q[3];
rz(-2.5803356) q[3];
sx q[3];
rz(-0.11229501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.955287) q[2];
sx q[2];
rz(-0.39262843) q[2];
sx q[2];
rz(0.10360959) q[2];
rz(-2.7747532) q[3];
sx q[3];
rz(-1.6678526) q[3];
sx q[3];
rz(-1.8240671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9570626) q[0];
sx q[0];
rz(-0.4758895) q[0];
sx q[0];
rz(-2.6153508) q[0];
rz(-1.2435675) q[1];
sx q[1];
rz(-1.4105816) q[1];
sx q[1];
rz(-0.19613656) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58914069) q[0];
sx q[0];
rz(-1.5040483) q[0];
sx q[0];
rz(-1.1132973) q[0];
rz(-pi) q[1];
rz(-0.7095269) q[2];
sx q[2];
rz(-1.6975743) q[2];
sx q[2];
rz(1.4122054) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89333234) q[1];
sx q[1];
rz(-0.47139097) q[1];
sx q[1];
rz(1.8880647) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1297471) q[3];
sx q[3];
rz(-1.5603754) q[3];
sx q[3];
rz(1.9886147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5474995) q[2];
sx q[2];
rz(-2.8671691) q[2];
sx q[2];
rz(2.7653232) q[2];
rz(2.2606692) q[3];
sx q[3];
rz(-1.8062402) q[3];
sx q[3];
rz(2.3295565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45102099) q[0];
sx q[0];
rz(-0.32830992) q[0];
sx q[0];
rz(-1.0115393) q[0];
rz(-0.66863376) q[1];
sx q[1];
rz(-1.9790383) q[1];
sx q[1];
rz(-0.54723251) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.512151) q[0];
sx q[0];
rz(-2.9392588) q[0];
sx q[0];
rz(0.095937177) q[0];
rz(2.3509432) q[2];
sx q[2];
rz(-1.8230652) q[2];
sx q[2];
rz(-2.6961435) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1909132) q[1];
sx q[1];
rz(-0.70531323) q[1];
sx q[1];
rz(-1.1937792) q[1];
x q[2];
rz(1.3667447) q[3];
sx q[3];
rz(-2.3746852) q[3];
sx q[3];
rz(-0.48108654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68613595) q[2];
sx q[2];
rz(-0.28527173) q[2];
sx q[2];
rz(1.5020465) q[2];
rz(2.5655668) q[3];
sx q[3];
rz(-1.4833769) q[3];
sx q[3];
rz(2.7801133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6308924) q[0];
sx q[0];
rz(-2.3402813) q[0];
sx q[0];
rz(-0.33962387) q[0];
rz(-2.6009808) q[1];
sx q[1];
rz(-0.70053354) q[1];
sx q[1];
rz(2.2115754) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4444816) q[0];
sx q[0];
rz(-2.9822095) q[0];
sx q[0];
rz(-1.1109933) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9564863) q[2];
sx q[2];
rz(-1.0218595) q[2];
sx q[2];
rz(-1.5654636) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8254357) q[1];
sx q[1];
rz(-0.47940578) q[1];
sx q[1];
rz(-1.5569219) q[1];
x q[2];
rz(0.50765462) q[3];
sx q[3];
rz(-2.1498907) q[3];
sx q[3];
rz(0.41989014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.12863079) q[2];
sx q[2];
rz(-1.4350812) q[2];
sx q[2];
rz(1.5677412) q[2];
rz(0.74357998) q[3];
sx q[3];
rz(-0.79135528) q[3];
sx q[3];
rz(2.1775406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-1.4670694) q[0];
sx q[0];
rz(-0.85551298) q[0];
sx q[0];
rz(-3.0849482) q[0];
rz(1.4757587) q[1];
sx q[1];
rz(-2.00878) q[1];
sx q[1];
rz(1.7830361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76641335) q[0];
sx q[0];
rz(-2.4046899) q[0];
sx q[0];
rz(2.3794258) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65746324) q[2];
sx q[2];
rz(-1.7320447) q[2];
sx q[2];
rz(-0.25875124) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8090933) q[1];
sx q[1];
rz(-2.1488071) q[1];
sx q[1];
rz(0.9851458) q[1];
rz(1.452946) q[3];
sx q[3];
rz(-2.1292344) q[3];
sx q[3];
rz(1.8359566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8654827) q[2];
sx q[2];
rz(-1.7296187) q[2];
sx q[2];
rz(2.0152246) q[2];
rz(-1.8051091) q[3];
sx q[3];
rz(-2.0650605) q[3];
sx q[3];
rz(1.0895464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69018501) q[0];
sx q[0];
rz(-1.0253588) q[0];
sx q[0];
rz(0.13352808) q[0];
rz(-2.1639157) q[1];
sx q[1];
rz(-2.1656499) q[1];
sx q[1];
rz(0.28688637) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8650019) q[0];
sx q[0];
rz(-1.6185068) q[0];
sx q[0];
rz(-0.86930958) q[0];
rz(-pi) q[1];
x q[1];
rz(3.098549) q[2];
sx q[2];
rz(-1.2428987) q[2];
sx q[2];
rz(-0.69059935) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.62090767) q[1];
sx q[1];
rz(-2.1322827) q[1];
sx q[1];
rz(-2.8771299) q[1];
x q[2];
rz(1.828015) q[3];
sx q[3];
rz(-0.99204274) q[3];
sx q[3];
rz(-2.677315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7202683) q[2];
sx q[2];
rz(-1.4807533) q[2];
sx q[2];
rz(1.6555017) q[2];
rz(2.7739575) q[3];
sx q[3];
rz(-1.0272107) q[3];
sx q[3];
rz(1.0953085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7169645) q[0];
sx q[0];
rz(-2.3724738) q[0];
sx q[0];
rz(2.6677483) q[0];
rz(-2.4773856) q[1];
sx q[1];
rz(-1.217548) q[1];
sx q[1];
rz(-1.3833822) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63603386) q[0];
sx q[0];
rz(-2.2657388) q[0];
sx q[0];
rz(0.64993422) q[0];
rz(-pi) q[1];
rz(1.6331312) q[2];
sx q[2];
rz(-2.7228055) q[2];
sx q[2];
rz(2.5588148) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.72402945) q[1];
sx q[1];
rz(-1.8010407) q[1];
sx q[1];
rz(-1.3413221) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17857213) q[3];
sx q[3];
rz(-1.7925486) q[3];
sx q[3];
rz(0.54820433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.15709269) q[2];
sx q[2];
rz(-1.4089156) q[2];
sx q[2];
rz(0.3248997) q[2];
rz(2.7609008) q[3];
sx q[3];
rz(-0.84563962) q[3];
sx q[3];
rz(2.0438173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.0710881) q[0];
sx q[0];
rz(-2.4779713) q[0];
sx q[0];
rz(2.2027503) q[0];
rz(-0.70010575) q[1];
sx q[1];
rz(-0.63799262) q[1];
sx q[1];
rz(2.8525888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8696006) q[0];
sx q[0];
rz(-1.4933407) q[0];
sx q[0];
rz(1.3415706) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6436725) q[2];
sx q[2];
rz(-2.1564031) q[2];
sx q[2];
rz(2.917054) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.9146198) q[1];
sx q[1];
rz(-1.4767854) q[1];
sx q[1];
rz(-2.0932122) q[1];
rz(-pi) q[2];
rz(-1.4988171) q[3];
sx q[3];
rz(-1.1907309) q[3];
sx q[3];
rz(-2.7615508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8812022) q[2];
sx q[2];
rz(-1.5217047) q[2];
sx q[2];
rz(0.41268665) q[2];
rz(0.9203426) q[3];
sx q[3];
rz(-2.6350382) q[3];
sx q[3];
rz(1.2296366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6270139) q[0];
sx q[0];
rz(-2.161442) q[0];
sx q[0];
rz(0.29943109) q[0];
rz(-2.0191655) q[1];
sx q[1];
rz(-1.2010682) q[1];
sx q[1];
rz(2.1103512) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.497488) q[0];
sx q[0];
rz(-1.2046763) q[0];
sx q[0];
rz(-3.1300504) q[0];
rz(-pi) q[1];
rz(-1.3716566) q[2];
sx q[2];
rz(-1.1165035) q[2];
sx q[2];
rz(0.21541883) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.448062) q[1];
sx q[1];
rz(-1.8603357) q[1];
sx q[1];
rz(2.2953643) q[1];
x q[2];
rz(2.2935872) q[3];
sx q[3];
rz(-2.3591745) q[3];
sx q[3];
rz(2.5944124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1624182) q[2];
sx q[2];
rz(-1.3434429) q[2];
sx q[2];
rz(-0.24270414) q[2];
rz(1.8414712) q[3];
sx q[3];
rz(-2.0827115) q[3];
sx q[3];
rz(2.9259031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31371394) q[0];
sx q[0];
rz(-0.21610459) q[0];
sx q[0];
rz(-1.4208273) q[0];
rz(2.8315262) q[1];
sx q[1];
rz(-2.2646751) q[1];
sx q[1];
rz(-1.5440595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2472947) q[0];
sx q[0];
rz(-1.159512) q[0];
sx q[0];
rz(0.572834) q[0];
x q[1];
rz(3.0528131) q[2];
sx q[2];
rz(-1.9114466) q[2];
sx q[2];
rz(1.3442049) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.76031715) q[1];
sx q[1];
rz(-2.7997428) q[1];
sx q[1];
rz(2.443497) q[1];
rz(-pi) q[2];
rz(-2.4214823) q[3];
sx q[3];
rz(-2.930899) q[3];
sx q[3];
rz(1.9062476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1342423) q[2];
sx q[2];
rz(-1.4126567) q[2];
sx q[2];
rz(-1.3664112) q[2];
rz(-0.8775231) q[3];
sx q[3];
rz(-1.6759422) q[3];
sx q[3];
rz(-2.2519978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9919745) q[0];
sx q[0];
rz(-1.5852954) q[0];
sx q[0];
rz(-1.4854767) q[0];
rz(0.86434518) q[1];
sx q[1];
rz(-1.0135916) q[1];
sx q[1];
rz(2.7085173) q[1];
rz(-2.6411459) q[2];
sx q[2];
rz(-1.7519578) q[2];
sx q[2];
rz(2.0448207) q[2];
rz(-1.0917615) q[3];
sx q[3];
rz(-1.619966) q[3];
sx q[3];
rz(1.0109284) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

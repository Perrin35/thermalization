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
rz(2.9721337) q[1];
sx q[1];
rz(-0.13154498) q[1];
sx q[1];
rz(2.0102672) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1695263) q[0];
sx q[0];
rz(-0.6837877) q[0];
sx q[0];
rz(1.0975361) q[0];
rz(-0.65931084) q[2];
sx q[2];
rz(-0.74987312) q[2];
sx q[2];
rz(1.4091968) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1234963) q[1];
sx q[1];
rz(-2.5275748) q[1];
sx q[1];
rz(1.2662925) q[1];
rz(-1.027358) q[3];
sx q[3];
rz(-1.7187331) q[3];
sx q[3];
rz(1.9224059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1863056) q[2];
sx q[2];
rz(-2.7489642) q[2];
sx q[2];
rz(-0.10360959) q[2];
rz(-0.3668395) q[3];
sx q[3];
rz(-1.47374) q[3];
sx q[3];
rz(-1.8240671) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1845301) q[0];
sx q[0];
rz(-0.4758895) q[0];
sx q[0];
rz(-2.6153508) q[0];
rz(1.8980252) q[1];
sx q[1];
rz(-1.7310111) q[1];
sx q[1];
rz(-2.9454561) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2945902) q[0];
sx q[0];
rz(-2.6795912) q[0];
sx q[0];
rz(1.7209956) q[0];
x q[1];
rz(1.7372484) q[2];
sx q[2];
rz(-0.86814082) q[2];
sx q[2];
rz(0.26674092) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.89333234) q[1];
sx q[1];
rz(-2.6702017) q[1];
sx q[1];
rz(-1.253528) q[1];
rz(1.581218) q[3];
sx q[3];
rz(-1.5589514) q[3];
sx q[3];
rz(-0.41769496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5474995) q[2];
sx q[2];
rz(-2.8671691) q[2];
sx q[2];
rz(0.37626949) q[2];
rz(-2.2606692) q[3];
sx q[3];
rz(-1.3353525) q[3];
sx q[3];
rz(-0.81203619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45102099) q[0];
sx q[0];
rz(-0.32830992) q[0];
sx q[0];
rz(2.1300533) q[0];
rz(-0.66863376) q[1];
sx q[1];
rz(-1.1625544) q[1];
sx q[1];
rz(0.54723251) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035346383) q[0];
sx q[0];
rz(-1.5515455) q[0];
sx q[0];
rz(-2.9401642) q[0];
x q[1];
rz(1.9220566) q[2];
sx q[2];
rz(-0.81163663) q[2];
sx q[2];
rz(1.3724302) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8142952) q[1];
sx q[1];
rz(-1.8117827) q[1];
sx q[1];
rz(2.2404033) q[1];
rz(1.3667447) q[3];
sx q[3];
rz(-2.3746852) q[3];
sx q[3];
rz(2.6605061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.68613595) q[2];
sx q[2];
rz(-0.28527173) q[2];
sx q[2];
rz(1.5020465) q[2];
rz(-0.57602588) q[3];
sx q[3];
rz(-1.4833769) q[3];
sx q[3];
rz(2.7801133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5107002) q[0];
sx q[0];
rz(-2.3402813) q[0];
sx q[0];
rz(-0.33962387) q[0];
rz(2.6009808) q[1];
sx q[1];
rz(-0.70053354) q[1];
sx q[1];
rz(0.9300173) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7226573) q[0];
sx q[0];
rz(-1.6412853) q[0];
sx q[0];
rz(1.427729) q[0];
x q[1];
rz(-1.9564863) q[2];
sx q[2];
rz(-2.1197332) q[2];
sx q[2];
rz(-1.5654636) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31615694) q[1];
sx q[1];
rz(-2.6621869) q[1];
sx q[1];
rz(-1.5569219) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93149473) q[3];
sx q[3];
rz(-2.3911797) q[3];
sx q[3];
rz(1.9281333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0129619) q[2];
sx q[2];
rz(-1.7065115) q[2];
sx q[2];
rz(-1.5738515) q[2];
rz(-0.74357998) q[3];
sx q[3];
rz(-2.3502374) q[3];
sx q[3];
rz(-0.96405205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6745233) q[0];
sx q[0];
rz(-0.85551298) q[0];
sx q[0];
rz(3.0849482) q[0];
rz(-1.665834) q[1];
sx q[1];
rz(-2.00878) q[1];
sx q[1];
rz(-1.3585565) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76641335) q[0];
sx q[0];
rz(-2.4046899) q[0];
sx q[0];
rz(0.76216682) q[0];
x q[1];
rz(2.4841294) q[2];
sx q[2];
rz(-1.4095479) q[2];
sx q[2];
rz(-0.25875124) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.92734071) q[1];
sx q[1];
rz(-2.3433629) q[1];
sx q[1];
rz(0.70294522) q[1];
rz(2.580022) q[3];
sx q[3];
rz(-1.4709146) q[3];
sx q[3];
rz(-2.9390854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.27611) q[2];
sx q[2];
rz(-1.411974) q[2];
sx q[2];
rz(-2.0152246) q[2];
rz(1.3364835) q[3];
sx q[3];
rz(-1.0765321) q[3];
sx q[3];
rz(-1.0895464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69018501) q[0];
sx q[0];
rz(-2.1162338) q[0];
sx q[0];
rz(0.13352808) q[0];
rz(2.1639157) q[1];
sx q[1];
rz(-0.97594273) q[1];
sx q[1];
rz(0.28688637) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27659076) q[0];
sx q[0];
rz(-1.6185068) q[0];
sx q[0];
rz(2.2722831) q[0];
rz(-pi) q[1];
rz(1.2426161) q[2];
sx q[2];
rz(-1.6115453) q[2];
sx q[2];
rz(0.894067) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.520685) q[1];
sx q[1];
rz(-2.1322827) q[1];
sx q[1];
rz(2.8771299) q[1];
x q[2];
rz(-1.3135776) q[3];
sx q[3];
rz(-0.99204274) q[3];
sx q[3];
rz(0.46427765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7202683) q[2];
sx q[2];
rz(-1.6608394) q[2];
sx q[2];
rz(-1.6555017) q[2];
rz(-0.36763516) q[3];
sx q[3];
rz(-2.1143819) q[3];
sx q[3];
rz(-1.0953085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7169645) q[0];
sx q[0];
rz(-2.3724738) q[0];
sx q[0];
rz(-2.6677483) q[0];
rz(-2.4773856) q[1];
sx q[1];
rz(-1.217548) q[1];
sx q[1];
rz(-1.3833822) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5067399) q[0];
sx q[0];
rz(-0.91284767) q[0];
sx q[0];
rz(0.94292504) q[0];
rz(-pi) q[1];
rz(-0.02772133) q[2];
sx q[2];
rz(-1.1528735) q[2];
sx q[2];
rz(0.51456416) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4175632) q[1];
sx q[1];
rz(-1.340552) q[1];
sx q[1];
rz(-1.3413221) q[1];
x q[2];
rz(1.7960127) q[3];
sx q[3];
rz(-1.3966421) q[3];
sx q[3];
rz(2.1586777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15709269) q[2];
sx q[2];
rz(-1.7326771) q[2];
sx q[2];
rz(0.3248997) q[2];
rz(-0.38069185) q[3];
sx q[3];
rz(-0.84563962) q[3];
sx q[3];
rz(-1.0977753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.070504524) q[0];
sx q[0];
rz(-2.4779713) q[0];
sx q[0];
rz(0.93884236) q[0];
rz(0.70010575) q[1];
sx q[1];
rz(-2.5036) q[1];
sx q[1];
rz(2.8525888) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61902133) q[0];
sx q[0];
rz(-0.24174015) q[0];
sx q[0];
rz(-1.2416583) q[0];
rz(-pi) q[1];
rz(2.1949057) q[2];
sx q[2];
rz(-2.3923229) q[2];
sx q[2];
rz(2.5890337) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2269729) q[1];
sx q[1];
rz(-1.4767854) q[1];
sx q[1];
rz(-2.0932122) q[1];
rz(-pi) q[2];
rz(-2.9634776) q[3];
sx q[3];
rz(-0.38649544) q[3];
sx q[3];
rz(-0.18807553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2603904) q[2];
sx q[2];
rz(-1.5217047) q[2];
sx q[2];
rz(-2.728906) q[2];
rz(-0.9203426) q[3];
sx q[3];
rz(-0.50655443) q[3];
sx q[3];
rz(1.2296366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6270139) q[0];
sx q[0];
rz(-2.161442) q[0];
sx q[0];
rz(-0.29943109) q[0];
rz(2.0191655) q[1];
sx q[1];
rz(-1.2010682) q[1];
sx q[1];
rz(-2.1103512) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64410463) q[0];
sx q[0];
rz(-1.2046763) q[0];
sx q[0];
rz(-3.1300504) q[0];
x q[1];
rz(1.769936) q[2];
sx q[2];
rz(-1.1165035) q[2];
sx q[2];
rz(-2.9261738) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.448062) q[1];
sx q[1];
rz(-1.281257) q[1];
sx q[1];
rz(-2.2953643) q[1];
rz(-2.5599285) q[3];
sx q[3];
rz(-1.0137034) q[3];
sx q[3];
rz(0.34639369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1624182) q[2];
sx q[2];
rz(-1.7981497) q[2];
sx q[2];
rz(-0.24270414) q[2];
rz(-1.3001214) q[3];
sx q[3];
rz(-2.0827115) q[3];
sx q[3];
rz(-0.21568957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8278787) q[0];
sx q[0];
rz(-2.9254881) q[0];
sx q[0];
rz(1.4208273) q[0];
rz(0.31006649) q[1];
sx q[1];
rz(-0.87691751) q[1];
sx q[1];
rz(-1.5440595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9102073) q[0];
sx q[0];
rz(-2.4500896) q[0];
sx q[0];
rz(-2.4639936) q[0];
x q[1];
rz(0.088779595) q[2];
sx q[2];
rz(-1.9114466) q[2];
sx q[2];
rz(-1.3442049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6621993) q[1];
sx q[1];
rz(-1.3536204) q[1];
sx q[1];
rz(2.875476) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.710911) q[3];
sx q[3];
rz(-1.4129253) q[3];
sx q[3];
rz(0.50413361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1342423) q[2];
sx q[2];
rz(-1.728936) q[2];
sx q[2];
rz(-1.3664112) q[2];
rz(0.8775231) q[3];
sx q[3];
rz(-1.4656504) q[3];
sx q[3];
rz(0.88959488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9919745) q[0];
sx q[0];
rz(-1.5852954) q[0];
sx q[0];
rz(-1.4854767) q[0];
rz(2.2772475) q[1];
sx q[1];
rz(-2.1280011) q[1];
sx q[1];
rz(-0.43307532) q[1];
rz(1.7766118) q[2];
sx q[2];
rz(-1.0792776) q[2];
sx q[2];
rz(0.57224032) q[2];
rz(-1.4644365) q[3];
sx q[3];
rz(-2.660236) q[3];
sx q[3];
rz(-0.65421792) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

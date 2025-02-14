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
rz(2.514578) q[0];
sx q[0];
rz(-0.63627565) q[0];
sx q[0];
rz(-1.0778435) q[0];
rz(1.0765422) q[1];
sx q[1];
rz(-0.62907469) q[1];
sx q[1];
rz(2.10973) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92562308) q[0];
sx q[0];
rz(-1.5703989) q[0];
sx q[0];
rz(-3.1374323) q[0];
rz(2.7912895) q[2];
sx q[2];
rz(-1.7997348) q[2];
sx q[2];
rz(0.13267429) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.061432211) q[1];
sx q[1];
rz(-0.54418635) q[1];
sx q[1];
rz(-0.70744343) q[1];
rz(-pi) q[2];
rz(-1.4937819) q[3];
sx q[3];
rz(-1.9056013) q[3];
sx q[3];
rz(0.96559262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91118139) q[2];
sx q[2];
rz(-1.958467) q[2];
sx q[2];
rz(-2.6739056) q[2];
rz(2.6905401) q[3];
sx q[3];
rz(-0.39877287) q[3];
sx q[3];
rz(-2.8635645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4462047) q[0];
sx q[0];
rz(-2.9187293) q[0];
sx q[0];
rz(-0.15431246) q[0];
rz(-0.34126869) q[1];
sx q[1];
rz(-0.35366615) q[1];
sx q[1];
rz(-2.7214859) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4898191) q[0];
sx q[0];
rz(-2.4742527) q[0];
sx q[0];
rz(-2.8267008) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6252363) q[2];
sx q[2];
rz(-2.4307361) q[2];
sx q[2];
rz(-2.9756851) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.005689335) q[1];
sx q[1];
rz(-1.5721675) q[1];
sx q[1];
rz(-2.6708009) q[1];
rz(2.3942727) q[3];
sx q[3];
rz(-1.6856717) q[3];
sx q[3];
rz(-2.9385027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61098617) q[2];
sx q[2];
rz(-0.89908081) q[2];
sx q[2];
rz(-2.3685624) q[2];
rz(-0.84345877) q[3];
sx q[3];
rz(-1.5638899) q[3];
sx q[3];
rz(-0.57190603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.968349) q[0];
sx q[0];
rz(-1.0272212) q[0];
sx q[0];
rz(-2.9395043) q[0];
rz(0.48187065) q[1];
sx q[1];
rz(-0.20129573) q[1];
sx q[1];
rz(-0.906382) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4545091) q[0];
sx q[0];
rz(-0.89620279) q[0];
sx q[0];
rz(-0.97954155) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34445539) q[2];
sx q[2];
rz(-1.1397425) q[2];
sx q[2];
rz(-1.5004917) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4587332) q[1];
sx q[1];
rz(-2.5762229) q[1];
sx q[1];
rz(0.302178) q[1];
x q[2];
rz(0.1679095) q[3];
sx q[3];
rz(-1.862226) q[3];
sx q[3];
rz(-1.341602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.071094461) q[2];
sx q[2];
rz(-2.0281894) q[2];
sx q[2];
rz(-0.58007288) q[2];
rz(1.5813367) q[3];
sx q[3];
rz(-1.7507078) q[3];
sx q[3];
rz(-1.4238547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.4105014) q[0];
sx q[0];
rz(-3.0803362) q[0];
sx q[0];
rz(0.047792338) q[0];
rz(1.2740678) q[1];
sx q[1];
rz(-1.85227) q[1];
sx q[1];
rz(-3.0963669) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52818283) q[0];
sx q[0];
rz(-0.65887132) q[0];
sx q[0];
rz(1.3263561) q[0];
rz(-pi) q[1];
rz(2.2061646) q[2];
sx q[2];
rz(-1.3029939) q[2];
sx q[2];
rz(1.2768465) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3056335) q[1];
sx q[1];
rz(-3.1352512) q[1];
sx q[1];
rz(-2.0672227) q[1];
rz(2.1158695) q[3];
sx q[3];
rz(-1.9506321) q[3];
sx q[3];
rz(1.7120509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.73146003) q[2];
sx q[2];
rz(-1.2612017) q[2];
sx q[2];
rz(0.88978466) q[2];
rz(-2.5951923) q[3];
sx q[3];
rz(-2.140464) q[3];
sx q[3];
rz(-2.8304097) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50362098) q[0];
sx q[0];
rz(-1.6872971) q[0];
sx q[0];
rz(1.4171492) q[0];
rz(-1.1196989) q[1];
sx q[1];
rz(-1.3816625) q[1];
sx q[1];
rz(-1.959257) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.946163) q[0];
sx q[0];
rz(-2.7481005) q[0];
sx q[0];
rz(-0.18184967) q[0];
x q[1];
rz(-2.6835713) q[2];
sx q[2];
rz(-1.1586896) q[2];
sx q[2];
rz(-2.1864708) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1685698) q[1];
sx q[1];
rz(-0.72271361) q[1];
sx q[1];
rz(0.074345868) q[1];
rz(-pi) q[2];
rz(-2.0145217) q[3];
sx q[3];
rz(-0.87165735) q[3];
sx q[3];
rz(-0.8389896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.22415367) q[2];
sx q[2];
rz(-2.3771693) q[2];
sx q[2];
rz(2.7177641) q[2];
rz(-2.0297) q[3];
sx q[3];
rz(-0.49911505) q[3];
sx q[3];
rz(-2.4901938) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3424585) q[0];
sx q[0];
rz(-3.0693711) q[0];
sx q[0];
rz(0.95241958) q[0];
rz(-0.77596387) q[1];
sx q[1];
rz(-2.2064078) q[1];
sx q[1];
rz(1.5663358) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090388894) q[0];
sx q[0];
rz(-3.0984833) q[0];
sx q[0];
rz(0.48844047) q[0];
x q[1];
rz(-2.1774954) q[2];
sx q[2];
rz(-1.5390054) q[2];
sx q[2];
rz(0.83039415) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4500931) q[1];
sx q[1];
rz(-0.97292346) q[1];
sx q[1];
rz(1.8370502) q[1];
x q[2];
rz(2.0616343) q[3];
sx q[3];
rz(-1.4075402) q[3];
sx q[3];
rz(1.3364244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7427407) q[2];
sx q[2];
rz(-2.3336918) q[2];
sx q[2];
rz(0.74637949) q[2];
rz(-1.0910723) q[3];
sx q[3];
rz(-1.7079791) q[3];
sx q[3];
rz(-2.5892042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93254507) q[0];
sx q[0];
rz(-0.046367558) q[0];
sx q[0];
rz(-0.5433425) q[0];
rz(-0.53687334) q[1];
sx q[1];
rz(-0.32506341) q[1];
sx q[1];
rz(3.0184025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42483271) q[0];
sx q[0];
rz(-2.2062629) q[0];
sx q[0];
rz(-0.038095564) q[0];
x q[1];
rz(-0.51080475) q[2];
sx q[2];
rz(-1.9419365) q[2];
sx q[2];
rz(-2.5696511) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.035288485) q[1];
sx q[1];
rz(-2.114944) q[1];
sx q[1];
rz(-2.418251) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76456548) q[3];
sx q[3];
rz(-0.95029059) q[3];
sx q[3];
rz(-0.71128774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3049551) q[2];
sx q[2];
rz(-1.387549) q[2];
sx q[2];
rz(0.48509625) q[2];
rz(1.9890316) q[3];
sx q[3];
rz(-1.7771114) q[3];
sx q[3];
rz(-0.047957234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3845859) q[0];
sx q[0];
rz(-2.204019) q[0];
sx q[0];
rz(-0.49569976) q[0];
rz(0.063830201) q[1];
sx q[1];
rz(-2.3921236) q[1];
sx q[1];
rz(-3.0705423) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.166693) q[0];
sx q[0];
rz(-1.3314795) q[0];
sx q[0];
rz(0.47310754) q[0];
rz(-2.536473) q[2];
sx q[2];
rz(-1.7150884) q[2];
sx q[2];
rz(2.6168106) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7531705) q[1];
sx q[1];
rz(-2.8533397) q[1];
sx q[1];
rz(2.1436195) q[1];
x q[2];
rz(2.6943915) q[3];
sx q[3];
rz(-1.4942823) q[3];
sx q[3];
rz(-2.2390128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0157328) q[2];
sx q[2];
rz(-1.4733529) q[2];
sx q[2];
rz(2.2802172) q[2];
rz(-1.8371948) q[3];
sx q[3];
rz(-2.5671037) q[3];
sx q[3];
rz(-1.5484352) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17643377) q[0];
sx q[0];
rz(-0.49089828) q[0];
sx q[0];
rz(-1.7023671) q[0];
rz(0.08055117) q[1];
sx q[1];
rz(-1.4713902) q[1];
sx q[1];
rz(1.390994) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0223402) q[0];
sx q[0];
rz(-2.0827522) q[0];
sx q[0];
rz(2.6520024) q[0];
rz(-1.5059696) q[2];
sx q[2];
rz(-0.30538426) q[2];
sx q[2];
rz(0.3665273) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30069253) q[1];
sx q[1];
rz(-2.024902) q[1];
sx q[1];
rz(-0.91030376) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3426314) q[3];
sx q[3];
rz(-0.50833251) q[3];
sx q[3];
rz(1.5647183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3535658) q[2];
sx q[2];
rz(-1.2949508) q[2];
sx q[2];
rz(3.0143152) q[2];
rz(-0.92653972) q[3];
sx q[3];
rz(-2.8997731) q[3];
sx q[3];
rz(-1.658879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7552898) q[0];
sx q[0];
rz(-0.64978623) q[0];
sx q[0];
rz(1.6408828) q[0];
rz(-2.3312148) q[1];
sx q[1];
rz(-2.2893298) q[1];
sx q[1];
rz(-2.3694029) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15543518) q[0];
sx q[0];
rz(-0.86830189) q[0];
sx q[0];
rz(-2.2566354) q[0];
rz(-pi) q[1];
rz(1.7842845) q[2];
sx q[2];
rz(-1.8597892) q[2];
sx q[2];
rz(1.9910781) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4160847) q[1];
sx q[1];
rz(-0.5495175) q[1];
sx q[1];
rz(-1.548442) q[1];
rz(-1.5072892) q[3];
sx q[3];
rz(-2.4606703) q[3];
sx q[3];
rz(-1.6977967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7605674) q[2];
sx q[2];
rz(-0.80485359) q[2];
sx q[2];
rz(-2.5123361) q[2];
rz(0.73090807) q[3];
sx q[3];
rz(-0.84354246) q[3];
sx q[3];
rz(1.4430911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6998941) q[0];
sx q[0];
rz(-0.48342539) q[0];
sx q[0];
rz(-0.40406686) q[0];
rz(0.40456698) q[1];
sx q[1];
rz(-1.4063944) q[1];
sx q[1];
rz(-0.688596) q[1];
rz(-2.8597833) q[2];
sx q[2];
rz(-1.6360313) q[2];
sx q[2];
rz(-2.4665063) q[2];
rz(2.6647207) q[3];
sx q[3];
rz(-0.58860368) q[3];
sx q[3];
rz(0.5640201) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

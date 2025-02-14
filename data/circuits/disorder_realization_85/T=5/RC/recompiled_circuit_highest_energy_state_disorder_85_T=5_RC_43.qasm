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
rz(-0.83952251) q[0];
sx q[0];
rz(-1.5877725) q[0];
sx q[0];
rz(0.96437469) q[0];
rz(-2.7388465) q[1];
sx q[1];
rz(-1.4637113) q[1];
sx q[1];
rz(0.95199624) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5247241) q[0];
sx q[0];
rz(-1.6646806) q[0];
sx q[0];
rz(-1.0794207) q[0];
rz(-0.098564221) q[2];
sx q[2];
rz(-1.9910887) q[2];
sx q[2];
rz(-1.9914876) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.18033548) q[1];
sx q[1];
rz(-1.1828268) q[1];
sx q[1];
rz(1.1731824) q[1];
x q[2];
rz(0.75593573) q[3];
sx q[3];
rz(-1.2159676) q[3];
sx q[3];
rz(-2.1906022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6160148) q[2];
sx q[2];
rz(-1.792045) q[2];
sx q[2];
rz(2.7363321) q[2];
rz(-1.1303834) q[3];
sx q[3];
rz(-0.78961343) q[3];
sx q[3];
rz(-1.214504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.90527117) q[0];
sx q[0];
rz(-3.1089678) q[0];
sx q[0];
rz(-2.3765833) q[0];
rz(-2.8969823) q[1];
sx q[1];
rz(-1.2449934) q[1];
sx q[1];
rz(-2.5027221) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8804358) q[0];
sx q[0];
rz(-0.64468274) q[0];
sx q[0];
rz(-3.0993483) q[0];
rz(-2.5193754) q[2];
sx q[2];
rz(-1.2583311) q[2];
sx q[2];
rz(1.086233) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7201405) q[1];
sx q[1];
rz(-1.1126592) q[1];
sx q[1];
rz(1.1790183) q[1];
rz(-pi) q[2];
rz(-1.2876533) q[3];
sx q[3];
rz(-0.75732175) q[3];
sx q[3];
rz(0.18268302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8581533) q[2];
sx q[2];
rz(-2.6946113) q[2];
sx q[2];
rz(-0.24432527) q[2];
rz(-2.0325932) q[3];
sx q[3];
rz(-1.8651483) q[3];
sx q[3];
rz(-0.23787704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37794381) q[0];
sx q[0];
rz(-2.0574594) q[0];
sx q[0];
rz(-1.9042683) q[0];
rz(1.411865) q[1];
sx q[1];
rz(-2.6928163) q[1];
sx q[1];
rz(1.3317187) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4770289) q[0];
sx q[0];
rz(-0.77858276) q[0];
sx q[0];
rz(-0.12523366) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2197415) q[2];
sx q[2];
rz(-1.7760385) q[2];
sx q[2];
rz(2.308941) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6159042) q[1];
sx q[1];
rz(-1.5153023) q[1];
sx q[1];
rz(1.6041808) q[1];
rz(-0.94966763) q[3];
sx q[3];
rz(-2.0916208) q[3];
sx q[3];
rz(2.4609712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1639013) q[2];
sx q[2];
rz(-1.0670263) q[2];
sx q[2];
rz(-2.4816371) q[2];
rz(1.6948949) q[3];
sx q[3];
rz(-1.4283254) q[3];
sx q[3];
rz(-1.1032633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8834943) q[0];
sx q[0];
rz(-2.3907008) q[0];
sx q[0];
rz(2.0080436) q[0];
rz(2.6253888) q[1];
sx q[1];
rz(-1.8323003) q[1];
sx q[1];
rz(-2.5925327) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7507197) q[0];
sx q[0];
rz(-1.4060806) q[0];
sx q[0];
rz(2.9956942) q[0];
x q[1];
rz(-0.98117657) q[2];
sx q[2];
rz(-1.7038888) q[2];
sx q[2];
rz(-0.99118628) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.085103206) q[1];
sx q[1];
rz(-1.6360549) q[1];
sx q[1];
rz(-0.27537055) q[1];
x q[2];
rz(-1.9677343) q[3];
sx q[3];
rz(-1.8324495) q[3];
sx q[3];
rz(1.73571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0893726) q[2];
sx q[2];
rz(-2.1286271) q[2];
sx q[2];
rz(-2.842556) q[2];
rz(-0.61797577) q[3];
sx q[3];
rz(-1.0733913) q[3];
sx q[3];
rz(-1.7834024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1842136) q[0];
sx q[0];
rz(-3.1035677) q[0];
sx q[0];
rz(0.46603611) q[0];
rz(2.274463) q[1];
sx q[1];
rz(-2.6301818) q[1];
sx q[1];
rz(-2.3066511) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9098489) q[0];
sx q[0];
rz(-0.0590894) q[0];
sx q[0];
rz(0.50414796) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.088836121) q[2];
sx q[2];
rz(-1.4918909) q[2];
sx q[2];
rz(-0.14278459) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9795436) q[1];
sx q[1];
rz(-2.0503373) q[1];
sx q[1];
rz(-0.060258106) q[1];
x q[2];
rz(-0.72839177) q[3];
sx q[3];
rz(-1.0755223) q[3];
sx q[3];
rz(0.033897696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9338108) q[2];
sx q[2];
rz(-0.89263478) q[2];
sx q[2];
rz(0.031522838) q[2];
rz(-2.3837714) q[3];
sx q[3];
rz(-1.1008681) q[3];
sx q[3];
rz(-0.56263629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.23444489) q[0];
sx q[0];
rz(-1.7666768) q[0];
sx q[0];
rz(3.1100519) q[0];
rz(0.08983287) q[1];
sx q[1];
rz(-0.49982163) q[1];
sx q[1];
rz(2.8057742) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0425371) q[0];
sx q[0];
rz(-0.49199781) q[0];
sx q[0];
rz(2.3943731) q[0];
x q[1];
rz(-2.8029595) q[2];
sx q[2];
rz(-2.0969318) q[2];
sx q[2];
rz(-2.3735685) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.3396158) q[1];
sx q[1];
rz(-1.812979) q[1];
sx q[1];
rz(2.1276317) q[1];
rz(-pi) q[2];
rz(2.8434136) q[3];
sx q[3];
rz(-2.0070873) q[3];
sx q[3];
rz(-2.0661708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.971656) q[2];
sx q[2];
rz(-2.0270831) q[2];
sx q[2];
rz(-1.2706832) q[2];
rz(0.058567889) q[3];
sx q[3];
rz(-1.6978369) q[3];
sx q[3];
rz(-0.34846714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.612065) q[0];
sx q[0];
rz(-2.6030774) q[0];
sx q[0];
rz(2.8756323) q[0];
rz(-0.82414857) q[1];
sx q[1];
rz(-1.3105086) q[1];
sx q[1];
rz(1.5305653) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6067995) q[0];
sx q[0];
rz(-1.3086119) q[0];
sx q[0];
rz(1.8305837) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0826999) q[2];
sx q[2];
rz(-2.4097171) q[2];
sx q[2];
rz(0.99725396) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2648523) q[1];
sx q[1];
rz(-2.1365668) q[1];
sx q[1];
rz(-1.2575722) q[1];
rz(-pi) q[2];
rz(0.77218382) q[3];
sx q[3];
rz(-1.8265542) q[3];
sx q[3];
rz(-2.135475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1223569) q[2];
sx q[2];
rz(-1.1160206) q[2];
sx q[2];
rz(-0.96987152) q[2];
rz(2.9979749) q[3];
sx q[3];
rz(-0.93846455) q[3];
sx q[3];
rz(-2.5041049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.403991) q[0];
sx q[0];
rz(-2.8856475) q[0];
sx q[0];
rz(3.0373489) q[0];
rz(-0.741611) q[1];
sx q[1];
rz(-0.18709083) q[1];
sx q[1];
rz(0.43824497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61442626) q[0];
sx q[0];
rz(-2.2402555) q[0];
sx q[0];
rz(1.7015083) q[0];
rz(1.9096228) q[2];
sx q[2];
rz(-2.5164587) q[2];
sx q[2];
rz(1.8115316) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91854535) q[1];
sx q[1];
rz(-2.1761736) q[1];
sx q[1];
rz(3.1001841) q[1];
x q[2];
rz(-0.53536587) q[3];
sx q[3];
rz(-1.4450049) q[3];
sx q[3];
rz(-2.4411502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7354342) q[2];
sx q[2];
rz(-0.27190748) q[2];
sx q[2];
rz(-2.820365) q[2];
rz(0.99212232) q[3];
sx q[3];
rz(-1.0935874) q[3];
sx q[3];
rz(-1.7295624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.853249) q[0];
sx q[0];
rz(-3.0217331) q[0];
sx q[0];
rz(2.993592) q[0];
rz(-1.1389698) q[1];
sx q[1];
rz(-2.482246) q[1];
sx q[1];
rz(-2.151087) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55427021) q[0];
sx q[0];
rz(-2.0519216) q[0];
sx q[0];
rz(1.6552129) q[0];
rz(-pi) q[1];
rz(-1.1075145) q[2];
sx q[2];
rz(-0.89416891) q[2];
sx q[2];
rz(-2.7156242) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79117486) q[1];
sx q[1];
rz(-1.5759849) q[1];
sx q[1];
rz(-1.5395916) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0664987) q[3];
sx q[3];
rz(-1.9073581) q[3];
sx q[3];
rz(2.9838205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.972435) q[2];
sx q[2];
rz(-2.3507698) q[2];
sx q[2];
rz(-2.5992744) q[2];
rz(1.4646685) q[3];
sx q[3];
rz(-1.2277536) q[3];
sx q[3];
rz(2.1577468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.9024314) q[0];
sx q[0];
rz(-0.79817525) q[0];
sx q[0];
rz(-2.8620128) q[0];
rz(2.3565049) q[1];
sx q[1];
rz(-2.3468192) q[1];
sx q[1];
rz(2.7395172) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.077686) q[0];
sx q[0];
rz(-2.5878667) q[0];
sx q[0];
rz(-2.5095449) q[0];
x q[1];
rz(-0.78530757) q[2];
sx q[2];
rz(-2.1673137) q[2];
sx q[2];
rz(-1.9526838) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6732503) q[1];
sx q[1];
rz(-1.7627009) q[1];
sx q[1];
rz(0.90190355) q[1];
rz(-pi) q[2];
rz(-2.8567186) q[3];
sx q[3];
rz(-1.1732297) q[3];
sx q[3];
rz(-1.8101102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98099199) q[2];
sx q[2];
rz(-2.4805785) q[2];
sx q[2];
rz(0.92657363) q[2];
rz(-2.5911234) q[3];
sx q[3];
rz(-2.2678943) q[3];
sx q[3];
rz(1.2800823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67113039) q[0];
sx q[0];
rz(-2.0736546) q[0];
sx q[0];
rz(0.14508844) q[0];
rz(-0.078770272) q[1];
sx q[1];
rz(-1.4046184) q[1];
sx q[1];
rz(1.2974993) q[1];
rz(-2.3433122) q[2];
sx q[2];
rz(-1.5599121) q[2];
sx q[2];
rz(-2.6679089) q[2];
rz(-1.79803) q[3];
sx q[3];
rz(-0.5323635) q[3];
sx q[3];
rz(-1.3743286) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

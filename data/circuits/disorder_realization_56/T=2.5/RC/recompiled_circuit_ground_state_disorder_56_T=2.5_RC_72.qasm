OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.174515) q[0];
sx q[0];
rz(1.4580026) q[0];
sx q[0];
rz(9.7249029) q[0];
rz(-1.9441654) q[1];
sx q[1];
rz(-1.5265042) q[1];
sx q[1];
rz(1.5010897) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34632698) q[0];
sx q[0];
rz(-3.1259968) q[0];
sx q[0];
rz(-0.25176425) q[0];
rz(-pi) q[1];
rz(1.1741224) q[2];
sx q[2];
rz(-2.3331917) q[2];
sx q[2];
rz(-1.8044745) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7607494) q[1];
sx q[1];
rz(-0.92807209) q[1];
sx q[1];
rz(1.2973143) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80880161) q[3];
sx q[3];
rz(-1.5743269) q[3];
sx q[3];
rz(0.88909066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40452051) q[2];
sx q[2];
rz(-1.9998735) q[2];
sx q[2];
rz(0.70297757) q[2];
rz(0.56973714) q[3];
sx q[3];
rz(-0.8786141) q[3];
sx q[3];
rz(1.5574633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1341781) q[0];
sx q[0];
rz(-0.61892048) q[0];
sx q[0];
rz(1.9084357) q[0];
rz(-0.59457072) q[1];
sx q[1];
rz(-2.1023127) q[1];
sx q[1];
rz(-2.0436683) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6195589) q[0];
sx q[0];
rz(-0.21169835) q[0];
sx q[0];
rz(1.2489399) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64867257) q[2];
sx q[2];
rz(-0.59174109) q[2];
sx q[2];
rz(1.6206738) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1008901) q[1];
sx q[1];
rz(-1.6482254) q[1];
sx q[1];
rz(3.1162683) q[1];
rz(-pi) q[2];
rz(2.3499346) q[3];
sx q[3];
rz(-1.5053621) q[3];
sx q[3];
rz(1.017788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.42830959) q[2];
sx q[2];
rz(-1.2998591) q[2];
sx q[2];
rz(1.606288) q[2];
rz(-2.4226268) q[3];
sx q[3];
rz(-2.1956367) q[3];
sx q[3];
rz(-0.21397056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83800256) q[0];
sx q[0];
rz(-2.328023) q[0];
sx q[0];
rz(-2.549262) q[0];
rz(-1.8968286) q[1];
sx q[1];
rz(-2.0015621) q[1];
sx q[1];
rz(2.9959784) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94731047) q[0];
sx q[0];
rz(-0.5049476) q[0];
sx q[0];
rz(0.28316747) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4136366) q[2];
sx q[2];
rz(-2.4470235) q[2];
sx q[2];
rz(2.5099011) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56359953) q[1];
sx q[1];
rz(-1.0219928) q[1];
sx q[1];
rz(2.8072458) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26992195) q[3];
sx q[3];
rz(-2.0441491) q[3];
sx q[3];
rz(-0.68655754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.738872) q[2];
sx q[2];
rz(-2.1108184) q[2];
sx q[2];
rz(-1.2325475) q[2];
rz(-0.23972073) q[3];
sx q[3];
rz(-2.0870049) q[3];
sx q[3];
rz(-2.6565552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0204912) q[0];
sx q[0];
rz(-0.70206577) q[0];
sx q[0];
rz(0.030601587) q[0];
rz(-2.5864511) q[1];
sx q[1];
rz(-1.6629013) q[1];
sx q[1];
rz(2.0487002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7858148) q[0];
sx q[0];
rz(-0.91453881) q[0];
sx q[0];
rz(-2.2758765) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.061528607) q[2];
sx q[2];
rz(-2.3752651) q[2];
sx q[2];
rz(-1.7249223) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.38072571) q[1];
sx q[1];
rz(-2.4740701) q[1];
sx q[1];
rz(1.236626) q[1];
x q[2];
rz(2.2956156) q[3];
sx q[3];
rz(-1.032077) q[3];
sx q[3];
rz(1.2332476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0838202) q[2];
sx q[2];
rz(-1.7744935) q[2];
sx q[2];
rz(-0.27285451) q[2];
rz(1.3911635) q[3];
sx q[3];
rz(-2.302156) q[3];
sx q[3];
rz(2.1896037) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6735753) q[0];
sx q[0];
rz(-1.8033569) q[0];
sx q[0];
rz(1.8433628) q[0];
rz(-0.6908373) q[1];
sx q[1];
rz(-1.9419443) q[1];
sx q[1];
rz(-1.615049) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0218791) q[0];
sx q[0];
rz(-1.0391616) q[0];
sx q[0];
rz(0.74703947) q[0];
rz(-pi) q[1];
rz(-1.2882907) q[2];
sx q[2];
rz(-1.9108678) q[2];
sx q[2];
rz(-1.9486537) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5954819) q[1];
sx q[1];
rz(-2.271436) q[1];
sx q[1];
rz(-3.107772) q[1];
x q[2];
rz(0.68256179) q[3];
sx q[3];
rz(-0.45387156) q[3];
sx q[3];
rz(-2.0477432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4096628) q[2];
sx q[2];
rz(-2.3356428) q[2];
sx q[2];
rz(-0.16732495) q[2];
rz(-1.3903728) q[3];
sx q[3];
rz(-0.53283397) q[3];
sx q[3];
rz(-0.65756857) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.14199) q[0];
sx q[0];
rz(-2.7775192) q[0];
sx q[0];
rz(-1.863119) q[0];
rz(1.4356042) q[1];
sx q[1];
rz(-1.6537138) q[1];
sx q[1];
rz(-2.7659168) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4091051) q[0];
sx q[0];
rz(-1.300749) q[0];
sx q[0];
rz(2.9069515) q[0];
rz(-2.9922688) q[2];
sx q[2];
rz(-1.4027601) q[2];
sx q[2];
rz(-0.72650331) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3766119) q[1];
sx q[1];
rz(-1.3479479) q[1];
sx q[1];
rz(1.8972023) q[1];
x q[2];
rz(2.385684) q[3];
sx q[3];
rz(-1.7497853) q[3];
sx q[3];
rz(2.8132015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0756691) q[2];
sx q[2];
rz(-1.064294) q[2];
sx q[2];
rz(-0.86090487) q[2];
rz(1.143645) q[3];
sx q[3];
rz(-2.836561) q[3];
sx q[3];
rz(-2.8633269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.025573108) q[0];
sx q[0];
rz(-1.8208068) q[0];
sx q[0];
rz(0.75468165) q[0];
rz(2.2017551) q[1];
sx q[1];
rz(-2.266423) q[1];
sx q[1];
rz(-1.2831203) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5425576) q[0];
sx q[0];
rz(-1.577566) q[0];
sx q[0];
rz(-1.8775131) q[0];
rz(-pi) q[1];
rz(2.2290984) q[2];
sx q[2];
rz(-2.5844838) q[2];
sx q[2];
rz(-2.0832286) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2312113) q[1];
sx q[1];
rz(-1.8613792) q[1];
sx q[1];
rz(0.39931764) q[1];
rz(-0.8831034) q[3];
sx q[3];
rz(-1.9839483) q[3];
sx q[3];
rz(0.93078155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54520404) q[2];
sx q[2];
rz(-2.0241006) q[2];
sx q[2];
rz(-1.0835353) q[2];
rz(-0.74294535) q[3];
sx q[3];
rz(-1.5321621) q[3];
sx q[3];
rz(1.4325745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98899406) q[0];
sx q[0];
rz(-2.3567663) q[0];
sx q[0];
rz(2.3866744) q[0];
rz(-2.6424291) q[1];
sx q[1];
rz(-0.9639591) q[1];
sx q[1];
rz(2.4430433) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7592418) q[0];
sx q[0];
rz(-2.1672591) q[0];
sx q[0];
rz(-1.1940707) q[0];
rz(-pi) q[1];
rz(-2.6766308) q[2];
sx q[2];
rz(-1.7002281) q[2];
sx q[2];
rz(-1.429582) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0341728) q[1];
sx q[1];
rz(-2.1028958) q[1];
sx q[1];
rz(-0.023253127) q[1];
rz(-pi) q[2];
rz(-1.5683062) q[3];
sx q[3];
rz(-0.41908136) q[3];
sx q[3];
rz(2.1650139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90427202) q[2];
sx q[2];
rz(-0.99232435) q[2];
sx q[2];
rz(0.80238706) q[2];
rz(2.5896416) q[3];
sx q[3];
rz(-2.3410083) q[3];
sx q[3];
rz(-2.5210181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.17860831) q[0];
sx q[0];
rz(-1.4405788) q[0];
sx q[0];
rz(-2.0340023) q[0];
rz(-1.3811318) q[1];
sx q[1];
rz(-1.3637204) q[1];
sx q[1];
rz(1.507087) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3699781) q[0];
sx q[0];
rz(-0.42015892) q[0];
sx q[0];
rz(2.7922947) q[0];
rz(-pi) q[1];
rz(2.8232024) q[2];
sx q[2];
rz(-1.9620775) q[2];
sx q[2];
rz(-1.1423781) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9025165) q[1];
sx q[1];
rz(-2.7488144) q[1];
sx q[1];
rz(2.7536489) q[1];
rz(-pi) q[2];
rz(1.4816557) q[3];
sx q[3];
rz(-0.74588585) q[3];
sx q[3];
rz(-0.28320492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0845118) q[2];
sx q[2];
rz(-2.9456186) q[2];
sx q[2];
rz(-1.8692807) q[2];
rz(-1.6614527) q[3];
sx q[3];
rz(-2.0062165) q[3];
sx q[3];
rz(-2.541466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24499527) q[0];
sx q[0];
rz(-2.5801881) q[0];
sx q[0];
rz(-1.8580612) q[0];
rz(3.1237579) q[1];
sx q[1];
rz(-0.63648883) q[1];
sx q[1];
rz(-3.0160115) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9367919) q[0];
sx q[0];
rz(-1.5849587) q[0];
sx q[0];
rz(2.1250179) q[0];
rz(-pi) q[1];
rz(0.56350817) q[2];
sx q[2];
rz(-0.76125188) q[2];
sx q[2];
rz(0.13680563) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.24841845) q[1];
sx q[1];
rz(-1.1496953) q[1];
sx q[1];
rz(-1.1014654) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8379032) q[3];
sx q[3];
rz(-1.3767164) q[3];
sx q[3];
rz(1.1052856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3632726) q[2];
sx q[2];
rz(-1.2691701) q[2];
sx q[2];
rz(-0.8030836) q[2];
rz(-0.17656365) q[3];
sx q[3];
rz(-1.3253788) q[3];
sx q[3];
rz(-2.363502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9618027) q[0];
sx q[0];
rz(-0.50146865) q[0];
sx q[0];
rz(1.5487221) q[0];
rz(1.8190307) q[1];
sx q[1];
rz(-1.7772728) q[1];
sx q[1];
rz(-1.5900236) q[1];
rz(-1.5079458) q[2];
sx q[2];
rz(-1.3286776) q[2];
sx q[2];
rz(1.5059289) q[2];
rz(0.17038067) q[3];
sx q[3];
rz(-0.34210404) q[3];
sx q[3];
rz(0.31755527) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7425793) q[0];
sx q[0];
rz(-1.9785545) q[0];
sx q[0];
rz(1.1385588) q[0];
rz(-0.63602716) q[1];
sx q[1];
rz(-0.50995246) q[1];
sx q[1];
rz(0.62503254) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.910927) q[0];
sx q[0];
rz(-1.5269484) q[0];
sx q[0];
rz(0.14104985) q[0];
x q[1];
rz(-1.6051859) q[2];
sx q[2];
rz(-0.54044881) q[2];
sx q[2];
rz(-2.2480223) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.6909862) q[1];
sx q[1];
rz(-2.4391973) q[1];
sx q[1];
rz(-0.90249636) q[1];
rz(-pi) q[2];
rz(-0.19842438) q[3];
sx q[3];
rz(-0.82272595) q[3];
sx q[3];
rz(-2.2446475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.92224377) q[2];
sx q[2];
rz(-0.85180989) q[2];
sx q[2];
rz(2.7395978) q[2];
rz(2.3123907) q[3];
sx q[3];
rz(-1.3196557) q[3];
sx q[3];
rz(1.3623665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336693) q[0];
sx q[0];
rz(-0.55183721) q[0];
sx q[0];
rz(2.0358987) q[0];
rz(2.3528174) q[1];
sx q[1];
rz(-0.62216798) q[1];
sx q[1];
rz(2.6355991) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9049022) q[0];
sx q[0];
rz(-2.8941029) q[0];
sx q[0];
rz(-0.31544073) q[0];
x q[1];
rz(2.3621735) q[2];
sx q[2];
rz(-2.1723618) q[2];
sx q[2];
rz(-2.0815408) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6044652) q[1];
sx q[1];
rz(-1.2187276) q[1];
sx q[1];
rz(0.2381937) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6384505) q[3];
sx q[3];
rz(-0.50853679) q[3];
sx q[3];
rz(-0.10756216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2371696) q[2];
sx q[2];
rz(-1.0639031) q[2];
sx q[2];
rz(2.2696631) q[2];
rz(1.0540086) q[3];
sx q[3];
rz(-1.0796615) q[3];
sx q[3];
rz(-1.4409298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0525518) q[0];
sx q[0];
rz(-0.65610039) q[0];
sx q[0];
rz(-1.22714) q[0];
rz(2.6350806) q[1];
sx q[1];
rz(-2.3498693) q[1];
sx q[1];
rz(2.2167749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.310257) q[0];
sx q[0];
rz(-1.4228205) q[0];
sx q[0];
rz(-3.0382115) q[0];
x q[1];
rz(3.001865) q[2];
sx q[2];
rz(-2.681571) q[2];
sx q[2];
rz(-2.7163028) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.23942415) q[1];
sx q[1];
rz(-1.8730867) q[1];
sx q[1];
rz(-0.92642529) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6684202) q[3];
sx q[3];
rz(-0.71801502) q[3];
sx q[3];
rz(-2.022746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0792803) q[2];
sx q[2];
rz(-1.5996876) q[2];
sx q[2];
rz(2.6463032) q[2];
rz(1.7700206) q[3];
sx q[3];
rz(-0.74426952) q[3];
sx q[3];
rz(-1.3220968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1247193) q[0];
sx q[0];
rz(-1.158241) q[0];
sx q[0];
rz(2.25368) q[0];
rz(-1.3814231) q[1];
sx q[1];
rz(-2.1235178) q[1];
sx q[1];
rz(2.6240614) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71947384) q[0];
sx q[0];
rz(-1.5848012) q[0];
sx q[0];
rz(-1.3097358) q[0];
rz(-pi) q[1];
rz(-2.5597455) q[2];
sx q[2];
rz(-0.80815017) q[2];
sx q[2];
rz(0.86607546) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1932839) q[1];
sx q[1];
rz(-1.6689321) q[1];
sx q[1];
rz(-3.1317284) q[1];
rz(-pi) q[2];
rz(-0.66756396) q[3];
sx q[3];
rz(-2.5680411) q[3];
sx q[3];
rz(1.1508326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0722644) q[2];
sx q[2];
rz(-1.1257409) q[2];
sx q[2];
rz(0.98108712) q[2];
rz(-0.4782933) q[3];
sx q[3];
rz(-2.7893453) q[3];
sx q[3];
rz(-0.52811629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9454055) q[0];
sx q[0];
rz(-1.7178752) q[0];
sx q[0];
rz(-3.1086573) q[0];
rz(1.6163274) q[1];
sx q[1];
rz(-0.5739637) q[1];
sx q[1];
rz(0.93564916) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.118555) q[0];
sx q[0];
rz(-1.584125) q[0];
sx q[0];
rz(0.47001919) q[0];
rz(-pi) q[1];
rz(2.6695365) q[2];
sx q[2];
rz(-0.85534125) q[2];
sx q[2];
rz(-3.1394342) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5577573) q[1];
sx q[1];
rz(-1.5084861) q[1];
sx q[1];
rz(2.3642217) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35772459) q[3];
sx q[3];
rz(-2.5265363) q[3];
sx q[3];
rz(-1.1411203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4719438) q[2];
sx q[2];
rz(-0.38413298) q[2];
sx q[2];
rz(-1.6432537) q[2];
rz(2.5390427) q[3];
sx q[3];
rz(-0.82585255) q[3];
sx q[3];
rz(1.2171827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5947386) q[0];
sx q[0];
rz(-0.29404077) q[0];
sx q[0];
rz(0.038851693) q[0];
rz(-0.36000577) q[1];
sx q[1];
rz(-1.729676) q[1];
sx q[1];
rz(0.14883277) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1921944) q[0];
sx q[0];
rz(-1.697505) q[0];
sx q[0];
rz(-1.697798) q[0];
x q[1];
rz(2.5720949) q[2];
sx q[2];
rz(-1.6844771) q[2];
sx q[2];
rz(-1.4194429) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.383713) q[1];
sx q[1];
rz(-2.0603018) q[1];
sx q[1];
rz(0.70546024) q[1];
rz(-pi) q[2];
rz(2.7018188) q[3];
sx q[3];
rz(-2.5086702) q[3];
sx q[3];
rz(1.1732111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0845324) q[2];
sx q[2];
rz(-2.3134505) q[2];
sx q[2];
rz(2.7866936) q[2];
rz(2.0830294) q[3];
sx q[3];
rz(-0.94616977) q[3];
sx q[3];
rz(-2.6088349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8232089) q[0];
sx q[0];
rz(-1.4223149) q[0];
sx q[0];
rz(-0.80867714) q[0];
rz(-3.0478364) q[1];
sx q[1];
rz(-2.3010727) q[1];
sx q[1];
rz(-2.7309928) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4559993) q[0];
sx q[0];
rz(-2.5846722) q[0];
sx q[0];
rz(-0.81057517) q[0];
rz(-pi) q[1];
rz(-1.3099395) q[2];
sx q[2];
rz(-2.2078288) q[2];
sx q[2];
rz(2.2230679) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0714662) q[1];
sx q[1];
rz(-1.9570597) q[1];
sx q[1];
rz(1.9420619) q[1];
x q[2];
rz(-2.9867044) q[3];
sx q[3];
rz(-2.4714176) q[3];
sx q[3];
rz(0.46524276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9984596) q[2];
sx q[2];
rz(-2.748558) q[2];
sx q[2];
rz(3.105799) q[2];
rz(-1.3068457) q[3];
sx q[3];
rz(-1.1467609) q[3];
sx q[3];
rz(0.80865639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3150113) q[0];
sx q[0];
rz(-0.97493521) q[0];
sx q[0];
rz(-2.5052729) q[0];
rz(-2.4782205) q[1];
sx q[1];
rz(-1.1095108) q[1];
sx q[1];
rz(0.62796193) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54701824) q[0];
sx q[0];
rz(-1.5421405) q[0];
sx q[0];
rz(-3.0648607) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3971252) q[2];
sx q[2];
rz(-1.6518804) q[2];
sx q[2];
rz(2.9604767) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9414569) q[1];
sx q[1];
rz(-1.0147328) q[1];
sx q[1];
rz(2.0651613) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9598755) q[3];
sx q[3];
rz(-2.9510806) q[3];
sx q[3];
rz(0.2976892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5593354) q[2];
sx q[2];
rz(-0.11056837) q[2];
sx q[2];
rz(-1.7250693) q[2];
rz(2.0420117) q[3];
sx q[3];
rz(-2.1018335) q[3];
sx q[3];
rz(-0.85272461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.002554) q[0];
sx q[0];
rz(-0.89109963) q[0];
sx q[0];
rz(0.6148327) q[0];
rz(2.6333206) q[1];
sx q[1];
rz(-0.95242396) q[1];
sx q[1];
rz(-2.8654548) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17441882) q[0];
sx q[0];
rz(-1.1456657) q[0];
sx q[0];
rz(-0.14226724) q[0];
rz(-pi) q[1];
rz(0.15900826) q[2];
sx q[2];
rz(-1.0463042) q[2];
sx q[2];
rz(-0.25689143) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.32807531) q[1];
sx q[1];
rz(-2.349252) q[1];
sx q[1];
rz(0.66285073) q[1];
x q[2];
rz(2.2302141) q[3];
sx q[3];
rz(-1.8448534) q[3];
sx q[3];
rz(2.1942133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95163661) q[2];
sx q[2];
rz(-1.0389682) q[2];
sx q[2];
rz(-0.098380066) q[2];
rz(1.240587) q[3];
sx q[3];
rz(-2.0799347) q[3];
sx q[3];
rz(-0.043717535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28546178) q[0];
sx q[0];
rz(-2.0933445) q[0];
sx q[0];
rz(0.35828006) q[0];
rz(-2.1367392) q[1];
sx q[1];
rz(-2.256911) q[1];
sx q[1];
rz(2.7994432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8703259) q[0];
sx q[0];
rz(-1.7881601) q[0];
sx q[0];
rz(-2.5408265) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.762879) q[2];
sx q[2];
rz(-1.0193664) q[2];
sx q[2];
rz(0.086520606) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.085386203) q[1];
sx q[1];
rz(-1.4634988) q[1];
sx q[1];
rz(0.57319586) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10327424) q[3];
sx q[3];
rz(-1.2359796) q[3];
sx q[3];
rz(2.0044199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3442012) q[2];
sx q[2];
rz(-0.48365334) q[2];
sx q[2];
rz(2.7698216) q[2];
rz(-0.19112912) q[3];
sx q[3];
rz(-1.1647977) q[3];
sx q[3];
rz(-2.718486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38901781) q[0];
sx q[0];
rz(-0.49674635) q[0];
sx q[0];
rz(2.4865271) q[0];
rz(0.3599421) q[1];
sx q[1];
rz(-1.5725726) q[1];
sx q[1];
rz(-1.6307065) q[1];
rz(-0.062535738) q[2];
sx q[2];
rz(-1.2690496) q[2];
sx q[2];
rz(-3.0435111) q[2];
rz(1.5866847) q[3];
sx q[3];
rz(-2.3835328) q[3];
sx q[3];
rz(2.6450139) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

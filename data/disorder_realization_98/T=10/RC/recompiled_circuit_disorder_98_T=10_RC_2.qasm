OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37031072) q[0];
sx q[0];
rz(4.1376576) q[0];
sx q[0];
rz(7.1538038) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(-0.14970782) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4275442) q[0];
sx q[0];
rz(-2.2872426) q[0];
sx q[0];
rz(0.9057522) q[0];
x q[1];
rz(2.3157273) q[2];
sx q[2];
rz(-0.68544938) q[2];
sx q[2];
rz(-0.65537383) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.183061) q[1];
sx q[1];
rz(-1.5325938) q[1];
sx q[1];
rz(1.7512291) q[1];
rz(-pi) q[2];
rz(-1.0204131) q[3];
sx q[3];
rz(-1.7341511) q[3];
sx q[3];
rz(-2.6640716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.25847882) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(0.70409888) q[2];
rz(-2.1885833) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(-1.7378418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(-2.5090704) q[0];
rz(-0.44644341) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(0.65223637) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95603847) q[0];
sx q[0];
rz(-1.5421876) q[0];
sx q[0];
rz(1.5831328) q[0];
rz(-pi) q[1];
rz(-1.2781497) q[2];
sx q[2];
rz(-1.3424982) q[2];
sx q[2];
rz(1.4676859) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3995041) q[1];
sx q[1];
rz(-1.794145) q[1];
sx q[1];
rz(0.29466596) q[1];
rz(-pi) q[2];
rz(0.84516256) q[3];
sx q[3];
rz(-2.3361428) q[3];
sx q[3];
rz(-1.0172539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5474881) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(-1.9799505) q[2];
rz(1.9836327) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(-1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050425477) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(-2.2913349) q[0];
rz(-2.6440874) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(-1.7920378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3916546) q[0];
sx q[0];
rz(-0.91021252) q[0];
sx q[0];
rz(-1.2664938) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0565874) q[2];
sx q[2];
rz(-0.075767013) q[2];
sx q[2];
rz(-0.55759341) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68223665) q[1];
sx q[1];
rz(-2.5767234) q[1];
sx q[1];
rz(2.6911246) q[1];
x q[2];
rz(1.549167) q[3];
sx q[3];
rz(-0.4261741) q[3];
sx q[3];
rz(0.91859761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0321908) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(2.2375977) q[2];
rz(-0.30113014) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(-1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49144739) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(1.7310463) q[0];
rz(-2.5097805) q[1];
sx q[1];
rz(-1.3316863) q[1];
sx q[1];
rz(3.1052123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43544337) q[0];
sx q[0];
rz(-0.94067803) q[0];
sx q[0];
rz(-2.8773727) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58156275) q[2];
sx q[2];
rz(-1.8459324) q[2];
sx q[2];
rz(-3.0951701) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.034678) q[1];
sx q[1];
rz(-2.2639096) q[1];
sx q[1];
rz(2.9702859) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4724457) q[3];
sx q[3];
rz(-1.3373168) q[3];
sx q[3];
rz(1.0055621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92695421) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(1.9533763) q[2];
rz(0.67048091) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(2.5454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(-0.16648509) q[0];
rz(2.3855551) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(0.23434815) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7974632) q[0];
sx q[0];
rz(-1.2327694) q[0];
sx q[0];
rz(-2.5812838) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.095586153) q[2];
sx q[2];
rz(-1.0786973) q[2];
sx q[2];
rz(-2.0795859) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.105141) q[1];
sx q[1];
rz(-2.3804133) q[1];
sx q[1];
rz(2.9245604) q[1];
rz(-pi) q[2];
rz(-1.8946394) q[3];
sx q[3];
rz(-1.4959335) q[3];
sx q[3];
rz(0.74136855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4328737) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(-2.9210572) q[2];
rz(-2.7045414) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41546145) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(2.3244526) q[0];
rz(-2.5754886) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(1.1436499) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4139347) q[0];
sx q[0];
rz(-1.9415138) q[0];
sx q[0];
rz(1.5278221) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8191889) q[2];
sx q[2];
rz(-0.69903261) q[2];
sx q[2];
rz(-2.2088745) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7856577) q[1];
sx q[1];
rz(-2.6895084) q[1];
sx q[1];
rz(1.7788586) q[1];
x q[2];
rz(-1.850071) q[3];
sx q[3];
rz(-2.0719299) q[3];
sx q[3];
rz(0.53370332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3879261) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(-0.4894408) q[2];
rz(2.9135381) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(-2.7155546) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.56931) q[0];
sx q[0];
rz(-2.4991878) q[0];
sx q[0];
rz(1.2868767) q[0];
rz(-0.6634179) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(-1.2333262) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42686233) q[0];
sx q[0];
rz(-2.7573702) q[0];
sx q[0];
rz(-1.6003438) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8823422) q[2];
sx q[2];
rz(-0.8493648) q[2];
sx q[2];
rz(1.7989858) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.308097) q[1];
sx q[1];
rz(-1.4596246) q[1];
sx q[1];
rz(-1.584946) q[1];
rz(-0.60621467) q[3];
sx q[3];
rz(-0.57176916) q[3];
sx q[3];
rz(-1.9619463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.104091) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(2.1255778) q[2];
rz(3.0715023) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.12524097) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(-1.4455147) q[0];
rz(2.9267172) q[1];
sx q[1];
rz(-2.3863249) q[1];
sx q[1];
rz(1.8833556) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90342605) q[0];
sx q[0];
rz(-1.4362207) q[0];
sx q[0];
rz(1.1358791) q[0];
x q[1];
rz(-2.4969205) q[2];
sx q[2];
rz(-1.03627) q[2];
sx q[2];
rz(-1.0169741) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5453323) q[1];
sx q[1];
rz(-0.55570554) q[1];
sx q[1];
rz(-0.19821367) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2860252) q[3];
sx q[3];
rz(-1.9868402) q[3];
sx q[3];
rz(2.1895529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4593279) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(-2.5788467) q[2];
rz(2.941926) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(-2.8905706) q[0];
rz(-0.42731467) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(3.1138611) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6513034) q[0];
sx q[0];
rz(-2.2130744) q[0];
sx q[0];
rz(2.519033) q[0];
x q[1];
rz(2.2698195) q[2];
sx q[2];
rz(-1.7773526) q[2];
sx q[2];
rz(-2.2733462) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.30222826) q[1];
sx q[1];
rz(-1.6703969) q[1];
sx q[1];
rz(0.60217963) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9412882) q[3];
sx q[3];
rz(-0.071404608) q[3];
sx q[3];
rz(-2.2927473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.015908265) q[2];
sx q[2];
rz(-0.97970825) q[2];
sx q[2];
rz(-1.998418) q[2];
rz(0.14287359) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(-2.2843602) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3906355) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(-0.45387682) q[0];
rz(2.4699396) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(-2.8840816) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2902381) q[0];
sx q[0];
rz(-2.5332753) q[0];
sx q[0];
rz(0.71382199) q[0];
rz(-2.0673413) q[2];
sx q[2];
rz(-2.0601344) q[2];
sx q[2];
rz(0.41000965) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.820206) q[1];
sx q[1];
rz(-2.4452129) q[1];
sx q[1];
rz(-0.022298261) q[1];
x q[2];
rz(-2.0503644) q[3];
sx q[3];
rz(-1.7864831) q[3];
sx q[3];
rz(-1.2558162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8132849) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(2.5349687) q[2];
rz(-2.666752) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(-2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3474779) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(0.22656245) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(-0.67129927) q[2];
sx q[2];
rz(-2.56834) q[2];
sx q[2];
rz(2.7111531) q[2];
rz(-1.0971309) q[3];
sx q[3];
rz(-0.53285014) q[3];
sx q[3];
rz(0.22900029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

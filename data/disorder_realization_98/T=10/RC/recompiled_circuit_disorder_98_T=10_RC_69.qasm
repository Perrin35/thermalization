OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7712819) q[0];
sx q[0];
rz(-0.9960649) q[0];
sx q[0];
rz(2.2709742) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(-0.14970782) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71404845) q[0];
sx q[0];
rz(-0.85435003) q[0];
sx q[0];
rz(0.9057522) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50617354) q[2];
sx q[2];
rz(-2.0548327) q[2];
sx q[2];
rz(1.6137705) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.183061) q[1];
sx q[1];
rz(-1.5325938) q[1];
sx q[1];
rz(1.7512291) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0204131) q[3];
sx q[3];
rz(-1.7341511) q[3];
sx q[3];
rz(2.6640716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8831138) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(-2.4374938) q[2];
rz(2.1885833) q[3];
sx q[3];
rz(-2.1694031) q[3];
sx q[3];
rz(1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54884058) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(2.5090704) q[0];
rz(-0.44644341) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(-2.4893563) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3632293) q[0];
sx q[0];
rz(-0.031154545) q[0];
sx q[0];
rz(-2.7345783) q[0];
rz(-pi) q[1];
rz(-1.2781497) q[2];
sx q[2];
rz(-1.3424982) q[2];
sx q[2];
rz(-1.6739068) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3398509) q[1];
sx q[1];
rz(-2.7738214) q[1];
sx q[1];
rz(-2.4778609) q[1];
rz(-2.2964301) q[3];
sx q[3];
rz(-0.80544986) q[3];
sx q[3];
rz(-2.1243387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5474881) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(1.9799505) q[2];
rz(-1.9836327) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0911672) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(-2.2913349) q[0];
rz(0.49750528) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(1.3495548) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7499381) q[0];
sx q[0];
rz(-0.91021252) q[0];
sx q[0];
rz(-1.2664938) q[0];
x q[1];
rz(3.1042728) q[2];
sx q[2];
rz(-1.63675) q[2];
sx q[2];
rz(1.0730336) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9393443) q[1];
sx q[1];
rz(-2.073624) q[1];
sx q[1];
rz(-1.8400251) q[1];
rz(-pi) q[2];
rz(3.1317741) q[3];
sx q[3];
rz(-1.9968642) q[3];
sx q[3];
rz(0.94235086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1094018) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(-0.90399495) q[2];
rz(2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(-1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.6501453) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(-1.4105463) q[0];
rz(-2.5097805) q[1];
sx q[1];
rz(-1.3316863) q[1];
sx q[1];
rz(3.1052123) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8481537) q[0];
sx q[0];
rz(-1.7833976) q[0];
sx q[0];
rz(2.2178177) q[0];
x q[1];
rz(-0.47469791) q[2];
sx q[2];
rz(-0.63650741) q[2];
sx q[2];
rz(1.1324901) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8424884) q[1];
sx q[1];
rz(-2.4310388) q[1];
sx q[1];
rz(1.3683661) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7500238) q[3];
sx q[3];
rz(-0.25300004) q[3];
sx q[3];
rz(-2.5391425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92695421) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(1.1882163) q[2];
rz(0.67048091) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.4398414) q[0];
sx q[0];
rz(-1.259946) q[0];
sx q[0];
rz(2.9751076) q[0];
rz(2.3855551) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(-0.23434815) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4011824) q[0];
sx q[0];
rz(-0.64490841) q[0];
sx q[0];
rz(2.5572204) q[0];
rz(-pi) q[1];
rz(-1.3946103) q[2];
sx q[2];
rz(-2.641045) q[2];
sx q[2];
rz(-1.2622152) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4489331) q[1];
sx q[1];
rz(-1.4217136) q[1];
sx q[1];
rz(-0.74933185) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2469532) q[3];
sx q[3];
rz(-1.6456592) q[3];
sx q[3];
rz(0.74136855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4328737) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(-2.9210572) q[2];
rz(0.43705127) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7261312) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(-0.81714001) q[0];
rz(2.5754886) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(-1.9979427) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7276579) q[0];
sx q[0];
rz(-1.9415138) q[0];
sx q[0];
rz(1.5278221) q[0];
x q[1];
rz(1.3104865) q[2];
sx q[2];
rz(-2.227265) q[2];
sx q[2];
rz(1.3442163) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1144048) q[1];
sx q[1];
rz(-1.4804375) q[1];
sx q[1];
rz(1.1272217) q[1];
x q[2];
rz(-1.850071) q[3];
sx q[3];
rz(-2.0719299) q[3];
sx q[3];
rz(-2.6078893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3879261) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(0.4894408) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.56931) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(1.8547159) q[0];
rz(-2.4781748) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(1.9082665) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9702643) q[0];
sx q[0];
rz(-1.5597222) q[0];
sx q[0];
rz(1.1867255) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5160518) q[2];
sx q[2];
rz(-2.1893246) q[2];
sx q[2];
rz(-0.90460888) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.26112939) q[1];
sx q[1];
rz(-1.5848586) q[1];
sx q[1];
rz(-3.0304099) q[1];
rz(-pi) q[2];
rz(2.6551649) q[3];
sx q[3];
rz(-1.8842116) q[3];
sx q[3];
rz(3.0048971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.104091) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(1.0160149) q[2];
rz(-3.0715023) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(-1.0664553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(3.0163517) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(-1.4455147) q[0];
rz(0.21487543) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(1.8833556) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94851516) q[0];
sx q[0];
rz(-0.45398871) q[0];
sx q[0];
rz(1.2598739) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3636742) q[2];
sx q[2];
rz(-0.81233835) q[2];
sx q[2];
rz(-1.9922436) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80553493) q[1];
sx q[1];
rz(-1.4667257) q[1];
sx q[1];
rz(-0.54688262) q[1];
rz(-pi) q[2];
rz(-0.52950852) q[3];
sx q[3];
rz(-2.2141075) q[3];
sx q[3];
rz(0.95638004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4593279) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(-0.56274596) q[2];
rz(-0.19966666) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0257618) q[0];
sx q[0];
rz(-2.0985726) q[0];
sx q[0];
rz(-2.8905706) q[0];
rz(2.714278) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(-3.1138611) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32556191) q[0];
sx q[0];
rz(-2.0566018) q[0];
sx q[0];
rz(-2.3150139) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26720033) q[2];
sx q[2];
rz(-0.88951096) q[2];
sx q[2];
rz(-0.53182488) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.30222826) q[1];
sx q[1];
rz(-1.4711958) q[1];
sx q[1];
rz(-2.539413) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1157007) q[3];
sx q[3];
rz(-1.504244) q[3];
sx q[3];
rz(-2.6641012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(1.998418) q[2];
rz(0.14287359) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(0.85723248) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3906355) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(0.45387682) q[0];
rz(0.67165309) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(-0.25751105) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2430902) q[0];
sx q[0];
rz(-1.1872963) q[0];
sx q[0];
rz(2.6570508) q[0];
x q[1];
rz(1.0742513) q[2];
sx q[2];
rz(-2.0601344) q[2];
sx q[2];
rz(0.41000965) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35044893) q[1];
sx q[1];
rz(-2.2669683) q[1];
sx q[1];
rz(-1.552156) q[1];
x q[2];
rz(-1.1274687) q[3];
sx q[3];
rz(-0.52237288) q[3];
sx q[3];
rz(0.70538196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3283078) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(0.60662398) q[2];
rz(-2.666752) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-1.3474779) q[0];
sx q[0];
rz(-1.62962) q[0];
sx q[0];
rz(-0.99123065) q[0];
rz(0.22656245) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(2.6735641) q[2];
sx q[2];
rz(-1.9149018) q[2];
sx q[2];
rz(0.55185774) q[2];
rz(-0.26279454) q[3];
sx q[3];
rz(-1.1017208) q[3];
sx q[3];
rz(2.8337939) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

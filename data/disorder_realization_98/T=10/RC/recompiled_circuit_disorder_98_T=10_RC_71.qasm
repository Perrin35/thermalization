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
rz(-0.87061849) q[0];
rz(2.119996) q[1];
sx q[1];
rz(-2.8586913) q[1];
sx q[1];
rz(0.14970782) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5538841) q[0];
sx q[0];
rz(-0.93548488) q[0];
sx q[0];
rz(2.5250838) q[0];
x q[1];
rz(2.6354191) q[2];
sx q[2];
rz(-2.0548327) q[2];
sx q[2];
rz(1.5278221) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.61923164) q[1];
sx q[1];
rz(-1.751096) q[1];
sx q[1];
rz(-3.1027604) q[1];
rz(-pi) q[2];
rz(2.1211795) q[3];
sx q[3];
rz(-1.4074416) q[3];
sx q[3];
rz(-0.47752105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8831138) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(0.70409888) q[2];
rz(-0.95300931) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(-1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(2.5090704) q[0];
rz(-2.6951492) q[1];
sx q[1];
rz(-1.4182785) q[1];
sx q[1];
rz(-2.4893563) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1855542) q[0];
sx q[0];
rz(-1.5421876) q[0];
sx q[0];
rz(-1.5831328) q[0];
rz(-2.2488238) q[2];
sx q[2];
rz(-2.7724578) q[2];
sx q[2];
rz(-2.6004651) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7420885) q[1];
sx q[1];
rz(-1.794145) q[1];
sx q[1];
rz(-2.8469267) q[1];
x q[2];
rz(2.2324123) q[3];
sx q[3];
rz(-2.069807) q[3];
sx q[3];
rz(3.1391075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5474881) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(1.9799505) q[2];
rz(1.15796) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(-1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050425477) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(0.85025775) q[0];
rz(2.6440874) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(-1.3495548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3916546) q[0];
sx q[0];
rz(-0.91021252) q[0];
sx q[0];
rz(1.2664938) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1042728) q[2];
sx q[2];
rz(-1.63675) q[2];
sx q[2];
rz(-2.0685591) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.50074358) q[1];
sx q[1];
rz(-1.8060246) q[1];
sx q[1];
rz(2.6231223) q[1];
rz(-pi) q[2];
rz(0.0098185929) q[3];
sx q[3];
rz(-1.1447284) q[3];
sx q[3];
rz(0.94235086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1094018) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(1.3999456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6501453) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(-1.7310463) q[0];
rz(2.5097805) q[1];
sx q[1];
rz(-1.3316863) q[1];
sx q[1];
rz(-3.1052123) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.13658) q[0];
sx q[0];
rz(-0.67626017) q[0];
sx q[0];
rz(-1.9146634) q[0];
rz(-pi) q[1];
x q[1];
rz(1.245001) q[2];
sx q[2];
rz(-1.0137644) q[2];
sx q[2];
rz(1.4404802) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7155647) q[1];
sx q[1];
rz(-1.702311) q[1];
sx q[1];
rz(-2.2711666) q[1];
rz(-pi) q[2];
rz(-2.7500238) q[3];
sx q[3];
rz(-0.25300004) q[3];
sx q[3];
rz(2.5391425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.92695421) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(-1.1882163) q[2];
rz(-0.67048091) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(0.59613434) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4398414) q[0];
sx q[0];
rz(-1.259946) q[0];
sx q[0];
rz(0.16648509) q[0];
rz(-2.3855551) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(-0.23434815) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34412947) q[0];
sx q[0];
rz(-1.9088233) q[0];
sx q[0];
rz(2.5812838) q[0];
rz(-pi) q[1];
rz(0.095586153) q[2];
sx q[2];
rz(-2.0628953) q[2];
sx q[2];
rz(1.0620067) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6926596) q[1];
sx q[1];
rz(-1.7198791) q[1];
sx q[1];
rz(-0.74933185) q[1];
rz(1.8022728) q[3];
sx q[3];
rz(-2.8095062) q[3];
sx q[3];
rz(-1.0486697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
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
rz(2.7261312) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(-0.81714001) q[0];
rz(-2.5754886) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(-1.9979427) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4139347) q[0];
sx q[0];
rz(-1.2000788) q[0];
sx q[0];
rz(-1.5278221) q[0];
x q[1];
rz(2.8191889) q[2];
sx q[2];
rz(-2.44256) q[2];
sx q[2];
rz(-0.93271819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1144048) q[1];
sx q[1];
rz(-1.4804375) q[1];
sx q[1];
rz(-1.1272217) q[1];
rz(1.2915217) q[3];
sx q[3];
rz(-2.0719299) q[3];
sx q[3];
rz(0.53370332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75366655) q[2];
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
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7466) q[0];
sx q[0];
rz(-1.9548423) q[0];
sx q[0];
rz(0.011944255) q[0];
rz(-pi) q[1];
rz(-0.8823422) q[2];
sx q[2];
rz(-2.2922278) q[2];
sx q[2];
rz(-1.7989858) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8804633) q[1];
sx q[1];
rz(-1.556734) q[1];
sx q[1];
rz(0.1111828) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60621467) q[3];
sx q[3];
rz(-0.57176916) q[3];
sx q[3];
rz(-1.9619463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.104091) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(-2.1255778) q[2];
rz(-0.070090381) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(-1.0664553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0163517) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(-1.4455147) q[0];
rz(-2.9267172) q[1];
sx q[1];
rz(-2.3863249) q[1];
sx q[1];
rz(-1.8833556) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94851516) q[0];
sx q[0];
rz(-0.45398871) q[0];
sx q[0];
rz(-1.2598739) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2082981) q[2];
sx q[2];
rz(-1.0273232) q[2];
sx q[2];
rz(0.18804929) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3360577) q[1];
sx q[1];
rz(-1.6748669) q[1];
sx q[1];
rz(-0.54688262) q[1];
x q[2];
rz(-2.2860252) q[3];
sx q[3];
rz(-1.1547525) q[3];
sx q[3];
rz(2.1895529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.68226472) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(-2.5788467) q[2];
rz(-0.19966666) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257618) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(0.2510221) q[0];
rz(2.714278) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(0.02773157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32556191) q[0];
sx q[0];
rz(-2.0566018) q[0];
sx q[0];
rz(0.82657878) q[0];
x q[1];
rz(1.2559782) q[2];
sx q[2];
rz(-2.417649) q[2];
sx q[2];
rz(0.9418504) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8393644) q[1];
sx q[1];
rz(-1.4711958) q[1];
sx q[1];
rz(2.539413) q[1];
rz(1.9412882) q[3];
sx q[3];
rz(-3.070188) q[3];
sx q[3];
rz(-0.84884531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1256844) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(-1.998418) q[2];
rz(-0.14287359) q[3];
sx q[3];
rz(-1.5121127) q[3];
sx q[3];
rz(0.85723248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.7509572) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(0.45387682) q[0];
rz(-0.67165309) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(-0.25751105) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6637708) q[0];
sx q[0];
rz(-1.1241233) q[0];
sx q[0];
rz(-1.1429943) q[0];
rz(-pi) q[1];
rz(-1.0742513) q[2];
sx q[2];
rz(-2.0601344) q[2];
sx q[2];
rz(2.731583) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35044893) q[1];
sx q[1];
rz(-0.87462438) q[1];
sx q[1];
rz(-1.5894366) q[1];
rz(-pi) q[2];
rz(-2.8994843) q[3];
sx q[3];
rz(-1.1032411) q[3];
sx q[3];
rz(2.9374591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3283078) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(-2.5349687) q[2];
rz(-0.47484067) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3474779) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(2.9150302) q[1];
sx q[1];
rz(-1.7270052) q[1];
sx q[1];
rz(-2.5546767) q[1];
rz(-2.6735641) q[2];
sx q[2];
rz(-1.2266908) q[2];
sx q[2];
rz(-2.5897349) q[2];
rz(1.0971309) q[3];
sx q[3];
rz(-2.6087425) q[3];
sx q[3];
rz(-2.9125924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

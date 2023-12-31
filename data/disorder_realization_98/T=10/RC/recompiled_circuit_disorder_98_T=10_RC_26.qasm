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
rz(-2.1455278) q[0];
sx q[0];
rz(-2.2709742) q[0];
rz(-7.3047819) q[1];
sx q[1];
rz(2.8586913) q[1];
sx q[1];
rz(18.999264) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38116954) q[0];
sx q[0];
rz(-2.0548577) q[0];
sx q[0];
rz(-0.83597393) q[0];
rz(2.6354191) q[2];
sx q[2];
rz(-2.0548327) q[2];
sx q[2];
rz(-1.6137705) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7357199) q[1];
sx q[1];
rz(-0.18438965) q[1];
sx q[1];
rz(1.3609481) q[1];
rz(2.1211795) q[3];
sx q[3];
rz(-1.7341511) q[3];
sx q[3];
rz(0.47752105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.25847882) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(-0.70409888) q[2];
rz(0.95300931) q[3];
sx q[3];
rz(-2.1694031) q[3];
sx q[3];
rz(-1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(-0.63252226) q[0];
rz(-2.6951492) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(-0.65223637) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3632293) q[0];
sx q[0];
rz(-0.031154545) q[0];
sx q[0];
rz(2.7345783) q[0];
rz(-pi) q[1];
rz(2.2488238) q[2];
sx q[2];
rz(-2.7724578) q[2];
sx q[2];
rz(2.6004651) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0374239) q[1];
sx q[1];
rz(-1.8579322) q[1];
sx q[1];
rz(1.3377405) q[1];
rz(-pi) q[2];
rz(2.2324123) q[3];
sx q[3];
rz(-1.0717857) q[3];
sx q[3];
rz(-3.1391075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5474881) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(-1.9799505) q[2];
rz(1.9836327) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(-1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.050425477) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(2.2913349) q[0];
rz(-2.6440874) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(1.7920378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7499381) q[0];
sx q[0];
rz(-2.2313801) q[0];
sx q[0];
rz(1.2664938) q[0];
rz(-pi) q[1];
rz(0.037319855) q[2];
sx q[2];
rz(-1.5048426) q[2];
sx q[2];
rz(1.0730336) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.68223665) q[1];
sx q[1];
rz(-0.56486928) q[1];
sx q[1];
rz(2.6911246) q[1];
rz(1.1447103) q[3];
sx q[3];
rz(-1.5797371) q[3];
sx q[3];
rz(-2.5090891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1094018) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(-2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(-1.3999456) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49144739) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(1.4105463) q[0];
rz(0.63181216) q[1];
sx q[1];
rz(-1.3316863) q[1];
sx q[1];
rz(-0.036380336) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2934389) q[0];
sx q[0];
rz(-1.358195) q[0];
sx q[0];
rz(-0.92377499) q[0];
rz(-pi) q[1];
x q[1];
rz(1.245001) q[2];
sx q[2];
rz(-1.0137644) q[2];
sx q[2];
rz(1.4404802) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2991043) q[1];
sx q[1];
rz(-0.71055382) q[1];
sx q[1];
rz(-1.7732265) q[1];
rz(-2.7500238) q[3];
sx q[3];
rz(-0.25300004) q[3];
sx q[3];
rz(-0.60245017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92695421) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(-1.1882163) q[2];
rz(2.4711117) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(-2.5454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7017512) q[0];
sx q[0];
rz(-1.259946) q[0];
sx q[0];
rz(-0.16648509) q[0];
rz(-2.3855551) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(-2.9072445) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4317961) q[0];
sx q[0];
rz(-1.0456107) q[0];
sx q[0];
rz(1.1774506) q[0];
rz(0.095586153) q[2];
sx q[2];
rz(-2.0628953) q[2];
sx q[2];
rz(-2.0795859) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8824132) q[1];
sx q[1];
rz(-2.3098574) q[1];
sx q[1];
rz(1.3684567) q[1];
rz(-3.0626416) q[3];
sx q[3];
rz(-1.8936994) q[3];
sx q[3];
rz(-2.3372646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.70871893) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(2.9210572) q[2];
rz(-0.43705127) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(-0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7261312) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(0.81714001) q[0];
rz(-0.56610402) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(-1.9979427) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1724388) q[0];
sx q[0];
rz(-1.6108496) q[0];
sx q[0];
rz(2.7705631) q[0];
x q[1];
rz(-2.4684858) q[2];
sx q[2];
rz(-1.3654725) q[2];
sx q[2];
rz(-2.7538607) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5864582) q[1];
sx q[1];
rz(-1.1291593) q[1];
sx q[1];
rz(-0.099979062) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.850071) q[3];
sx q[3];
rz(-1.0696628) q[3];
sx q[3];
rz(-0.53370332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3879261) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(-2.6521519) q[2];
rz(-2.9135381) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(2.7155546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.56931) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(-1.8547159) q[0];
rz(2.4781748) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(-1.9082665) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7466) q[0];
sx q[0];
rz(-1.1867503) q[0];
sx q[0];
rz(-3.1296484) q[0];
rz(-0.85031021) q[2];
sx q[2];
rz(-1.0734953) q[2];
sx q[2];
rz(-0.26956272) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4349605) q[1];
sx q[1];
rz(-3.0295277) q[1];
sx q[1];
rz(3.0155165) q[1];
x q[2];
rz(-2.535378) q[3];
sx q[3];
rz(-0.57176916) q[3];
sx q[3];
rz(-1.1796463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.037501637) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(1.0160149) q[2];
rz(0.070090381) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(-2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12524097) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(1.4455147) q[0];
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
rz(2.5364752) q[0];
sx q[0];
rz(-1.1400756) q[0];
sx q[0];
rz(0.14819781) q[0];
x q[1];
rz(2.3636742) q[2];
sx q[2];
rz(-2.3292543) q[2];
sx q[2];
rz(1.1493491) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.59626034) q[1];
sx q[1];
rz(-2.5858871) q[1];
sx q[1];
rz(2.943379) q[1];
rz(-pi) q[2];
rz(2.1636837) q[3];
sx q[3];
rz(-2.3330354) q[3];
sx q[3];
rz(-2.9582994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4593279) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(-2.5788467) q[2];
rz(-0.19966666) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257618) q[0];
sx q[0];
rz(-2.0985726) q[0];
sx q[0];
rz(-0.2510221) q[0];
rz(2.714278) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(-3.1138611) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4902892) q[0];
sx q[0];
rz(-2.2130744) q[0];
sx q[0];
rz(2.519033) q[0];
rz(-0.26720033) q[2];
sx q[2];
rz(-2.2520817) q[2];
sx q[2];
rz(2.6097678) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7293538) q[1];
sx q[1];
rz(-2.5322399) q[1];
sx q[1];
rz(-2.9669697) q[1];
x q[2];
rz(-3.1157007) q[3];
sx q[3];
rz(-1.504244) q[3];
sx q[3];
rz(-0.47749146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1256844) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(1.1431747) q[2];
rz(2.9987191) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(2.2843602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7509572) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(-2.6877158) q[0];
rz(-0.67165309) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(-2.8840816) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2902381) q[0];
sx q[0];
rz(-0.60831735) q[0];
sx q[0];
rz(-2.4277707) q[0];
rz(-pi) q[1];
rz(1.0742513) q[2];
sx q[2];
rz(-1.0814582) q[2];
sx q[2];
rz(-0.41000965) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9092907) q[1];
sx q[1];
rz(-1.5564939) q[1];
sx q[1];
rz(2.4453352) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24210838) q[3];
sx q[3];
rz(-1.1032411) q[3];
sx q[3];
rz(0.20413354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3283078) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(0.60662398) q[2];
rz(-2.666752) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(0.84038466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-1.1889585) q[2];
sx q[2];
rz(-2.0094064) q[2];
sx q[2];
rz(-1.1878427) q[2];
rz(2.0444617) q[3];
sx q[3];
rz(-0.53285014) q[3];
sx q[3];
rz(0.22900029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-0.8166135) q[0];
sx q[0];
rz(-0.54083523) q[0];
sx q[0];
rz(1.1582561) q[0];
rz(0.090016063) q[1];
sx q[1];
rz(-2.6675192) q[1];
sx q[1];
rz(2.3064244) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7258361) q[0];
sx q[0];
rz(-0.8716363) q[0];
sx q[0];
rz(1.1582202) q[0];
x q[1];
rz(-1.8878172) q[2];
sx q[2];
rz(-0.95404139) q[2];
sx q[2];
rz(-2.1932909) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7427828) q[1];
sx q[1];
rz(-1.2025857) q[1];
sx q[1];
rz(1.3081934) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99906875) q[3];
sx q[3];
rz(-0.99410996) q[3];
sx q[3];
rz(3.1126693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0917255) q[2];
sx q[2];
rz(-0.89605248) q[2];
sx q[2];
rz(-0.33207616) q[2];
rz(2.8927228) q[3];
sx q[3];
rz(-1.9460461) q[3];
sx q[3];
rz(-2.2539049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87149549) q[0];
sx q[0];
rz(-0.29093727) q[0];
sx q[0];
rz(2.6089597) q[0];
rz(-0.10781413) q[1];
sx q[1];
rz(-2.0262521) q[1];
sx q[1];
rz(-0.32049387) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6924393) q[0];
sx q[0];
rz(-2.4863613) q[0];
sx q[0];
rz(-1.6754608) q[0];
rz(-pi) q[1];
rz(2.2079289) q[2];
sx q[2];
rz(-0.28544989) q[2];
sx q[2];
rz(-2.1674726) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7254843) q[1];
sx q[1];
rz(-1.4112765) q[1];
sx q[1];
rz(1.4166338) q[1];
rz(-pi) q[2];
rz(1.3555834) q[3];
sx q[3];
rz(-1.2485412) q[3];
sx q[3];
rz(2.2029869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.33114854) q[2];
sx q[2];
rz(-0.30511567) q[2];
sx q[2];
rz(-1.6395052) q[2];
rz(-0.33809996) q[3];
sx q[3];
rz(-0.88518849) q[3];
sx q[3];
rz(-2.6147208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51760393) q[0];
sx q[0];
rz(-1.9094587) q[0];
sx q[0];
rz(2.7841618) q[0];
rz(-0.39237157) q[1];
sx q[1];
rz(-2.3452499) q[1];
sx q[1];
rz(-1.9270814) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93335184) q[0];
sx q[0];
rz(-2.9989971) q[0];
sx q[0];
rz(-1.4487793) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25193314) q[2];
sx q[2];
rz(-2.5802543) q[2];
sx q[2];
rz(2.6643945) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7347265) q[1];
sx q[1];
rz(-0.38741437) q[1];
sx q[1];
rz(1.563036) q[1];
rz(0.037596627) q[3];
sx q[3];
rz(-0.98340935) q[3];
sx q[3];
rz(-1.8855349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7963205) q[2];
sx q[2];
rz(-0.97769633) q[2];
sx q[2];
rz(0.39682445) q[2];
rz(-1.8848298) q[3];
sx q[3];
rz(-1.7233012) q[3];
sx q[3];
rz(2.5313012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9339555) q[0];
sx q[0];
rz(-1.3960681) q[0];
sx q[0];
rz(-1.9130094) q[0];
rz(1.2127016) q[1];
sx q[1];
rz(-2.1804501) q[1];
sx q[1];
rz(0.62087762) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29781326) q[0];
sx q[0];
rz(-0.50619805) q[0];
sx q[0];
rz(-1.9202581) q[0];
rz(0.16571705) q[2];
sx q[2];
rz(-1.4172557) q[2];
sx q[2];
rz(2.4927947) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.5073237) q[1];
sx q[1];
rz(-2.7719927) q[1];
sx q[1];
rz(2.3564767) q[1];
rz(-pi) q[2];
rz(-1.6895164) q[3];
sx q[3];
rz(-1.190589) q[3];
sx q[3];
rz(-2.3480727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2584201) q[2];
sx q[2];
rz(-2.0032538) q[2];
sx q[2];
rz(-0.15667285) q[2];
rz(-2.2251718) q[3];
sx q[3];
rz(-2.1082924) q[3];
sx q[3];
rz(-2.5887183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2857392) q[0];
sx q[0];
rz(-1.1612949) q[0];
sx q[0];
rz(-0.39749417) q[0];
rz(-0.022857895) q[1];
sx q[1];
rz(-0.50067478) q[1];
sx q[1];
rz(-2.4668677) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0894273) q[0];
sx q[0];
rz(-2.8244655) q[0];
sx q[0];
rz(-0.84256147) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4762474) q[2];
sx q[2];
rz(-1.6942021) q[2];
sx q[2];
rz(-2.9837554) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0832842) q[1];
sx q[1];
rz(-1.4183908) q[1];
sx q[1];
rz(-0.17249523) q[1];
rz(-pi) q[2];
rz(2.2915693) q[3];
sx q[3];
rz(-0.40215079) q[3];
sx q[3];
rz(-1.7274203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5270093) q[2];
sx q[2];
rz(-0.60636568) q[2];
sx q[2];
rz(0.69925365) q[2];
rz(1.7953385) q[3];
sx q[3];
rz(-2.5739539) q[3];
sx q[3];
rz(0.7114555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4166477) q[0];
sx q[0];
rz(-0.41794932) q[0];
sx q[0];
rz(0.39837343) q[0];
rz(2.1669972) q[1];
sx q[1];
rz(-1.83788) q[1];
sx q[1];
rz(1.7030565) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4476623) q[0];
sx q[0];
rz(-1.8724144) q[0];
sx q[0];
rz(-1.680879) q[0];
x q[1];
rz(-1.826615) q[2];
sx q[2];
rz(-1.5808357) q[2];
sx q[2];
rz(0.61983392) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.12425646) q[1];
sx q[1];
rz(-1.8030232) q[1];
sx q[1];
rz(-0.88084765) q[1];
x q[2];
rz(1.3138198) q[3];
sx q[3];
rz(-2.0517618) q[3];
sx q[3];
rz(2.0723267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4633816) q[2];
sx q[2];
rz(-1.7794098) q[2];
sx q[2];
rz(-1.8360651) q[2];
rz(-1.3156923) q[3];
sx q[3];
rz(-1.7782327) q[3];
sx q[3];
rz(-2.2731884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0618133) q[0];
sx q[0];
rz(-2.7257305) q[0];
sx q[0];
rz(-1.6449991) q[0];
rz(2.1553701) q[1];
sx q[1];
rz(-1.5602427) q[1];
sx q[1];
rz(-2.644002) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7114899) q[0];
sx q[0];
rz(-0.27207366) q[0];
sx q[0];
rz(-1.9338694) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2642483) q[2];
sx q[2];
rz(-0.75226558) q[2];
sx q[2];
rz(1.330842) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5274682) q[1];
sx q[1];
rz(-1.9640199) q[1];
sx q[1];
rz(-2.9725463) q[1];
rz(-pi) q[2];
rz(-1.4769745) q[3];
sx q[3];
rz(-1.2872496) q[3];
sx q[3];
rz(2.97992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30279532) q[2];
sx q[2];
rz(-1.5668198) q[2];
sx q[2];
rz(0.29067579) q[2];
rz(-2.931328) q[3];
sx q[3];
rz(-0.99431521) q[3];
sx q[3];
rz(-1.0342342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5656089) q[0];
sx q[0];
rz(-1.1701595) q[0];
sx q[0];
rz(-1.1822816) q[0];
rz(2.9871509) q[1];
sx q[1];
rz(-1.4192105) q[1];
sx q[1];
rz(-0.90528893) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2574999) q[0];
sx q[0];
rz(-1.582889) q[0];
sx q[0];
rz(-1.8435249) q[0];
rz(0.45640517) q[2];
sx q[2];
rz(-1.785136) q[2];
sx q[2];
rz(-2.6021007) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61991751) q[1];
sx q[1];
rz(-0.82304919) q[1];
sx q[1];
rz(-1.770439) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21568294) q[3];
sx q[3];
rz(-1.2372011) q[3];
sx q[3];
rz(-1.622596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0901383) q[2];
sx q[2];
rz(-1.5608414) q[2];
sx q[2];
rz(2.1913989) q[2];
rz(-3.0692302) q[3];
sx q[3];
rz(-1.3772929) q[3];
sx q[3];
rz(-0.70934057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-2.7224834) q[0];
sx q[0];
rz(-0.95471946) q[0];
sx q[0];
rz(-1.0954274) q[0];
rz(2.0504045) q[1];
sx q[1];
rz(-2.2406816) q[1];
sx q[1];
rz(2.1464164) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68190116) q[0];
sx q[0];
rz(-2.2349305) q[0];
sx q[0];
rz(1.6524397) q[0];
x q[1];
rz(-1.9068933) q[2];
sx q[2];
rz(-1.8310412) q[2];
sx q[2];
rz(-1.1059831) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0867445) q[1];
sx q[1];
rz(-1.5252171) q[1];
sx q[1];
rz(3.0566825) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18539683) q[3];
sx q[3];
rz(-1.3747921) q[3];
sx q[3];
rz(-0.40962266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6301849) q[2];
sx q[2];
rz(-2.6738622) q[2];
sx q[2];
rz(2.4616145) q[2];
rz(-1.4558815) q[3];
sx q[3];
rz(-0.89434353) q[3];
sx q[3];
rz(-1.1792012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.42530123) q[0];
sx q[0];
rz(-0.93451262) q[0];
sx q[0];
rz(2.5262078) q[0];
rz(-1.3353434) q[1];
sx q[1];
rz(-1.6067303) q[1];
sx q[1];
rz(-1.3478442) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5292717) q[0];
sx q[0];
rz(-2.7833496) q[0];
sx q[0];
rz(-0.1310346) q[0];
rz(0.81735264) q[2];
sx q[2];
rz(-2.0869227) q[2];
sx q[2];
rz(1.942526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.19619689) q[1];
sx q[1];
rz(-2.7295503) q[1];
sx q[1];
rz(-0.41007385) q[1];
rz(-pi) q[2];
rz(2.3854783) q[3];
sx q[3];
rz(-1.4817837) q[3];
sx q[3];
rz(-1.4579888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.24796692) q[2];
sx q[2];
rz(-1.8411571) q[2];
sx q[2];
rz(0.51353961) q[2];
rz(0.14704554) q[3];
sx q[3];
rz(-2.5976318) q[3];
sx q[3];
rz(-1.2646382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2697987) q[0];
sx q[0];
rz(-2.2525621) q[0];
sx q[0];
rz(2.9004108) q[0];
rz(1.3006032) q[1];
sx q[1];
rz(-2.0722957) q[1];
sx q[1];
rz(3.121079) q[1];
rz(0.14221556) q[2];
sx q[2];
rz(-1.9420062) q[2];
sx q[2];
rz(0.26376482) q[2];
rz(2.5592531) q[3];
sx q[3];
rz(-2.4262541) q[3];
sx q[3];
rz(-2.2179009) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(2.9958148) q[0];
sx q[0];
rz(4.4895953) q[0];
sx q[0];
rz(9.0564981) q[0];
rz(-0.085973099) q[1];
sx q[1];
rz(-2.2987125) q[1];
sx q[1];
rz(0.24922961) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50685731) q[0];
sx q[0];
rz(-1.0930499) q[0];
sx q[0];
rz(-1.4450106) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31166844) q[2];
sx q[2];
rz(-1.0578007) q[2];
sx q[2];
rz(2.5530961) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.492384) q[1];
sx q[1];
rz(-1.8064152) q[1];
sx q[1];
rz(-2.431303) q[1];
rz(1.6843819) q[3];
sx q[3];
rz(-1.2143597) q[3];
sx q[3];
rz(2.4939052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2117846) q[2];
sx q[2];
rz(-2.4346508) q[2];
sx q[2];
rz(-0.24016538) q[2];
rz(0.54720488) q[3];
sx q[3];
rz(-1.6271084) q[3];
sx q[3];
rz(-2.1521294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3419679) q[0];
sx q[0];
rz(-2.7662321) q[0];
sx q[0];
rz(-2.4976835) q[0];
rz(2.0945235) q[1];
sx q[1];
rz(-1.964485) q[1];
sx q[1];
rz(0.26328304) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1564668) q[0];
sx q[0];
rz(-1.5432285) q[0];
sx q[0];
rz(-0.03027244) q[0];
rz(-pi) q[1];
rz(-3.075454) q[2];
sx q[2];
rz(-1.6355343) q[2];
sx q[2];
rz(-1.3548702) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.689641) q[1];
sx q[1];
rz(-2.3799689) q[1];
sx q[1];
rz(0.93859251) q[1];
x q[2];
rz(-2.7277522) q[3];
sx q[3];
rz(-2.1713421) q[3];
sx q[3];
rz(0.59970705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6836267) q[2];
sx q[2];
rz(-1.3764952) q[2];
sx q[2];
rz(-2.4158884) q[2];
rz(0.30722412) q[3];
sx q[3];
rz(-1.3064462) q[3];
sx q[3];
rz(-1.0544302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82045186) q[0];
sx q[0];
rz(-1.6571925) q[0];
sx q[0];
rz(2.3725574) q[0];
rz(0.10737315) q[1];
sx q[1];
rz(-1.702405) q[1];
sx q[1];
rz(0.69951397) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0663129) q[0];
sx q[0];
rz(-1.9121721) q[0];
sx q[0];
rz(0.34324788) q[0];
rz(-0.0064445297) q[2];
sx q[2];
rz(-0.61476189) q[2];
sx q[2];
rz(2.4561735) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.87772885) q[1];
sx q[1];
rz(-0.45719621) q[1];
sx q[1];
rz(1.3844107) q[1];
x q[2];
rz(2.2992976) q[3];
sx q[3];
rz(-0.96046042) q[3];
sx q[3];
rz(1.5092602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.37381441) q[2];
sx q[2];
rz(-2.5641597) q[2];
sx q[2];
rz(-0.90816298) q[2];
rz(-2.5670037) q[3];
sx q[3];
rz(-1.2535973) q[3];
sx q[3];
rz(-0.38578924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.813756) q[0];
sx q[0];
rz(-1.3264553) q[0];
sx q[0];
rz(-0.5089708) q[0];
rz(3.0834037) q[1];
sx q[1];
rz(-1.4570844) q[1];
sx q[1];
rz(1.0675272) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021085652) q[0];
sx q[0];
rz(-2.3160546) q[0];
sx q[0];
rz(-2.6322281) q[0];
rz(-pi) q[1];
rz(1.7420962) q[2];
sx q[2];
rz(-2.7242085) q[2];
sx q[2];
rz(2.4462552) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4219027) q[1];
sx q[1];
rz(-2.9585144) q[1];
sx q[1];
rz(0.94193926) q[1];
x q[2];
rz(-0.068693585) q[3];
sx q[3];
rz(-1.8634081) q[3];
sx q[3];
rz(0.13394395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.031614583) q[2];
sx q[2];
rz(-1.5974533) q[2];
sx q[2];
rz(-2.1515089) q[2];
rz(1.7456985) q[3];
sx q[3];
rz(-0.11002222) q[3];
sx q[3];
rz(1.2466189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9083967) q[0];
sx q[0];
rz(-1.3968503) q[0];
sx q[0];
rz(-2.0821849) q[0];
rz(1.1147095) q[1];
sx q[1];
rz(-1.474294) q[1];
sx q[1];
rz(-1.4310736) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1170332) q[0];
sx q[0];
rz(-2.9142671) q[0];
sx q[0];
rz(-0.84421279) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1226019) q[2];
sx q[2];
rz(-0.55010527) q[2];
sx q[2];
rz(-0.12223003) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5268577) q[1];
sx q[1];
rz(-0.6893553) q[1];
sx q[1];
rz(2.7562642) q[1];
rz(-pi) q[2];
rz(-0.21417136) q[3];
sx q[3];
rz(-1.7167544) q[3];
sx q[3];
rz(-2.0623061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7913671) q[2];
sx q[2];
rz(-2.3855049) q[2];
sx q[2];
rz(1.4678601) q[2];
rz(1.0427467) q[3];
sx q[3];
rz(-2.0839033) q[3];
sx q[3];
rz(-2.3500672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70822155) q[0];
sx q[0];
rz(-1.2935761) q[0];
sx q[0];
rz(2.9300387) q[0];
rz(-2.4987706) q[1];
sx q[1];
rz(-1.1245518) q[1];
sx q[1];
rz(-2.0932253) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12473561) q[0];
sx q[0];
rz(-2.2679459) q[0];
sx q[0];
rz(-1.4037031) q[0];
rz(-1.5177814) q[2];
sx q[2];
rz(-1.7496193) q[2];
sx q[2];
rz(-0.12763466) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.590946) q[1];
sx q[1];
rz(-2.5027486) q[1];
sx q[1];
rz(3.1025801) q[1];
x q[2];
rz(-3.0934342) q[3];
sx q[3];
rz(-2.1827896) q[3];
sx q[3];
rz(2.4186866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6362777) q[2];
sx q[2];
rz(-1.2039528) q[2];
sx q[2];
rz(-0.25406507) q[2];
rz(-2.9680179) q[3];
sx q[3];
rz(-2.5901399) q[3];
sx q[3];
rz(-1.180163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59488615) q[0];
sx q[0];
rz(-2.7017024) q[0];
sx q[0];
rz(2.34483) q[0];
rz(-0.88675371) q[1];
sx q[1];
rz(-1.3462857) q[1];
sx q[1];
rz(2.5409882) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4457339) q[0];
sx q[0];
rz(-1.2810871) q[0];
sx q[0];
rz(1.6281307) q[0];
x q[1];
rz(-1.7037292) q[2];
sx q[2];
rz(-0.80604751) q[2];
sx q[2];
rz(1.199276) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4485885) q[1];
sx q[1];
rz(-1.7327762) q[1];
sx q[1];
rz(0.88108351) q[1];
x q[2];
rz(-1.0375627) q[3];
sx q[3];
rz(-1.8775465) q[3];
sx q[3];
rz(-0.83707929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0299006) q[2];
sx q[2];
rz(-1.7039958) q[2];
sx q[2];
rz(-1.0924443) q[2];
rz(-1.7732636) q[3];
sx q[3];
rz(-0.60951257) q[3];
sx q[3];
rz(-2.7294559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71548897) q[0];
sx q[0];
rz(-2.0273209) q[0];
sx q[0];
rz(0.45147595) q[0];
rz(-3.0637947) q[1];
sx q[1];
rz(-2.1092238) q[1];
sx q[1];
rz(0.66158867) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1790893) q[0];
sx q[0];
rz(-2.668619) q[0];
sx q[0];
rz(-2.0684558) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0724474) q[2];
sx q[2];
rz(-1.5094286) q[2];
sx q[2];
rz(-2.5453486) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3319131) q[1];
sx q[1];
rz(-2.6617378) q[1];
sx q[1];
rz(-1.2820246) q[1];
rz(-pi) q[2];
rz(0.57406942) q[3];
sx q[3];
rz(-2.6791253) q[3];
sx q[3];
rz(1.9898349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0563858) q[2];
sx q[2];
rz(-1.2678009) q[2];
sx q[2];
rz(2.3351604) q[2];
rz(1.2422397) q[3];
sx q[3];
rz(-2.9759088) q[3];
sx q[3];
rz(-3.0012259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30271444) q[0];
sx q[0];
rz(-1.1352204) q[0];
sx q[0];
rz(-0.43933991) q[0];
rz(1.9413403) q[1];
sx q[1];
rz(-0.83963436) q[1];
sx q[1];
rz(2.1913948) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0341102) q[0];
sx q[0];
rz(-2.378298) q[0];
sx q[0];
rz(1.0176246) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9885811) q[2];
sx q[2];
rz(-1.2875612) q[2];
sx q[2];
rz(-2.6161414) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31017329) q[1];
sx q[1];
rz(-1.9260287) q[1];
sx q[1];
rz(0.11174496) q[1];
rz(-pi) q[2];
rz(2.1860649) q[3];
sx q[3];
rz(-2.5078109) q[3];
sx q[3];
rz(1.2891226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3085559) q[2];
sx q[2];
rz(-2.6335282) q[2];
sx q[2];
rz(1.9527831) q[2];
rz(3.0360119) q[3];
sx q[3];
rz(-0.9413541) q[3];
sx q[3];
rz(1.0134491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.062716) q[0];
sx q[0];
rz(-0.3083516) q[0];
sx q[0];
rz(2.4421332) q[0];
rz(-1.4106916) q[1];
sx q[1];
rz(-2.3532093) q[1];
sx q[1];
rz(-2.9404822) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3956665) q[0];
sx q[0];
rz(-0.89218436) q[0];
sx q[0];
rz(-2.1525488) q[0];
rz(-0.4847862) q[2];
sx q[2];
rz(-2.1109627) q[2];
sx q[2];
rz(0.44819427) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8751018) q[1];
sx q[1];
rz(-1.3492341) q[1];
sx q[1];
rz(-2.2699577) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6789959) q[3];
sx q[3];
rz(-2.9411081) q[3];
sx q[3];
rz(1.6311262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.19645277) q[2];
sx q[2];
rz(-1.2597193) q[2];
sx q[2];
rz(3.073976) q[2];
rz(-1.128528) q[3];
sx q[3];
rz(-2.0566548) q[3];
sx q[3];
rz(1.5170001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.6832798) q[0];
sx q[0];
rz(-1.0673609) q[0];
sx q[0];
rz(1.3491032) q[0];
rz(1.6364527) q[1];
sx q[1];
rz(-0.22947336) q[1];
sx q[1];
rz(-0.089692399) q[1];
rz(0.2978953) q[2];
sx q[2];
rz(-2.5174601) q[2];
sx q[2];
rz(-2.8909825) q[2];
rz(-2.811089) q[3];
sx q[3];
rz(-1.8310343) q[3];
sx q[3];
rz(-1.9559947) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

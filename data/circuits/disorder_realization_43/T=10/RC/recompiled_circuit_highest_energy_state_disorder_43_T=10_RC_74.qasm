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
rz(-2.4515406) q[0];
sx q[0];
rz(-2.2909988) q[0];
sx q[0];
rz(2.9414862) q[0];
rz(-0.19835681) q[1];
sx q[1];
rz(-0.52302066) q[1];
sx q[1];
rz(1.3318292) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9560741) q[0];
sx q[0];
rz(-1.3393434) q[0];
sx q[0];
rz(-1.1197907) q[0];
rz(-pi) q[1];
rz(1.0806141) q[2];
sx q[2];
rz(-2.1128383) q[2];
sx q[2];
rz(1.3203743) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3568452) q[1];
sx q[1];
rz(-2.4618399) q[1];
sx q[1];
rz(0.076514449) q[1];
rz(-pi) q[2];
rz(1.0947919) q[3];
sx q[3];
rz(-1.4972357) q[3];
sx q[3];
rz(-1.5293763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0693822) q[2];
sx q[2];
rz(-2.78077) q[2];
sx q[2];
rz(-1.3145831) q[2];
rz(2.2383111) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(0.34590736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8697934) q[0];
sx q[0];
rz(-1.6558187) q[0];
sx q[0];
rz(2.2858009) q[0];
rz(-2.3675512) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(-0.78786293) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5408527) q[0];
sx q[0];
rz(-2.8898281) q[0];
sx q[0];
rz(-1.4815848) q[0];
x q[1];
rz(-1.5283952) q[2];
sx q[2];
rz(-1.2268492) q[2];
sx q[2];
rz(-2.6478772) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1545002) q[1];
sx q[1];
rz(-0.56114158) q[1];
sx q[1];
rz(-0.89548703) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8773852) q[3];
sx q[3];
rz(-2.4572861) q[3];
sx q[3];
rz(2.5591889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.065980109) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(-0.39609972) q[2];
rz(-1.570545) q[3];
sx q[3];
rz(-1.6377662) q[3];
sx q[3];
rz(-1.3554696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9539255) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(-3.0809825) q[0];
rz(-0.68471471) q[1];
sx q[1];
rz(-1.4002607) q[1];
sx q[1];
rz(-2.8025467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1364839) q[0];
sx q[0];
rz(-1.5500796) q[0];
sx q[0];
rz(-0.58309515) q[0];
rz(-pi) q[1];
rz(-0.18908638) q[2];
sx q[2];
rz(-1.6918618) q[2];
sx q[2];
rz(-0.4640641) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.10968929) q[1];
sx q[1];
rz(-1.8432143) q[1];
sx q[1];
rz(-1.8930356) q[1];
x q[2];
rz(-2.364758) q[3];
sx q[3];
rz(-1.562444) q[3];
sx q[3];
rz(-1.9296822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6591349) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(0.47046146) q[2];
rz(0.86236924) q[3];
sx q[3];
rz(-0.44895288) q[3];
sx q[3];
rz(-1.5615777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3622793) q[0];
sx q[0];
rz(-1.1184432) q[0];
sx q[0];
rz(0.40801868) q[0];
rz(-1.3511924) q[1];
sx q[1];
rz(-1.8441169) q[1];
sx q[1];
rz(1.321235) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7119638) q[0];
sx q[0];
rz(-2.5416059) q[0];
sx q[0];
rz(-2.5866051) q[0];
x q[1];
rz(0.086060103) q[2];
sx q[2];
rz(-2.2241908) q[2];
sx q[2];
rz(1.6376405) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3776748) q[1];
sx q[1];
rz(-1.6435247) q[1];
sx q[1];
rz(0.12158981) q[1];
rz(-0.77797555) q[3];
sx q[3];
rz(-2.5730592) q[3];
sx q[3];
rz(-1.83873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9185751) q[2];
sx q[2];
rz(-1.7097946) q[2];
sx q[2];
rz(-2.4187386) q[2];
rz(-0.84960788) q[3];
sx q[3];
rz(-1.1725715) q[3];
sx q[3];
rz(1.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1866813) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(-0.55996672) q[0];
rz(0.2050744) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(2.8505039) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3223461) q[0];
sx q[0];
rz(-2.8221606) q[0];
sx q[0];
rz(0.56630212) q[0];
x q[1];
rz(-0.79925691) q[2];
sx q[2];
rz(-1.6148668) q[2];
sx q[2];
rz(1.8792626) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0138058) q[1];
sx q[1];
rz(-1.3723137) q[1];
sx q[1];
rz(-2.2747893) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3897393) q[3];
sx q[3];
rz(-1.6883779) q[3];
sx q[3];
rz(-0.29871179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9186972) q[2];
sx q[2];
rz(-1.4783858) q[2];
sx q[2];
rz(0.29829868) q[2];
rz(1.9722624) q[3];
sx q[3];
rz(-2.1746641) q[3];
sx q[3];
rz(1.6303308) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8339612) q[0];
sx q[0];
rz(-1.6809604) q[0];
sx q[0];
rz(-2.5391915) q[0];
rz(-2.7700453) q[1];
sx q[1];
rz(-1.5448152) q[1];
sx q[1];
rz(-1.0317624) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1790947) q[0];
sx q[0];
rz(-2.966605) q[0];
sx q[0];
rz(1.6618769) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40470064) q[2];
sx q[2];
rz(-1.7621303) q[2];
sx q[2];
rz(-1.2197242) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.25615869) q[1];
sx q[1];
rz(-2.6226461) q[1];
sx q[1];
rz(-0.38796723) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86182819) q[3];
sx q[3];
rz(-0.78360451) q[3];
sx q[3];
rz(-0.48101048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3370257) q[2];
sx q[2];
rz(-1.1013384) q[2];
sx q[2];
rz(1.5186914) q[2];
rz(0.8484146) q[3];
sx q[3];
rz(-2.266326) q[3];
sx q[3];
rz(-3.1019822) q[3];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0254211) q[0];
sx q[0];
rz(-3.0905368) q[0];
sx q[0];
rz(-0.99789944) q[0];
rz(-0.2746703) q[1];
sx q[1];
rz(-2.2614567) q[1];
sx q[1];
rz(-1.3708699) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8664198) q[0];
sx q[0];
rz(-2.201797) q[0];
sx q[0];
rz(-1.5874528) q[0];
rz(-pi) q[1];
rz(-1.9267124) q[2];
sx q[2];
rz(-1.7150262) q[2];
sx q[2];
rz(-3.1145417) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6295885) q[1];
sx q[1];
rz(-1.3301714) q[1];
sx q[1];
rz(-2.5179067) q[1];
rz(-pi) q[2];
rz(0.82084772) q[3];
sx q[3];
rz(-2.2162262) q[3];
sx q[3];
rz(2.5643333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2844598) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(-3.1119463) q[2];
rz(-0.61521411) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(-0.342338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0520332) q[0];
sx q[0];
rz(-2.679306) q[0];
sx q[0];
rz(-0.7830559) q[0];
rz(-1.1198593) q[1];
sx q[1];
rz(-0.66711396) q[1];
sx q[1];
rz(1.7436183) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8496765) q[0];
sx q[0];
rz(-1.9789985) q[0];
sx q[0];
rz(1.085039) q[0];
rz(-pi) q[1];
rz(0.88770788) q[2];
sx q[2];
rz(-2.7811433) q[2];
sx q[2];
rz(2.6101024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82247558) q[1];
sx q[1];
rz(-0.90293316) q[1];
sx q[1];
rz(3.130359) q[1];
rz(-pi) q[2];
rz(-2.3572631) q[3];
sx q[3];
rz(-1.722986) q[3];
sx q[3];
rz(1.0206211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29584259) q[2];
sx q[2];
rz(-2.5551899) q[2];
sx q[2];
rz(1.8208549) q[2];
rz(2.1218421) q[3];
sx q[3];
rz(-1.4239862) q[3];
sx q[3];
rz(-1.4199055) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19972292) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(0.95440188) q[0];
rz(-1.9873387) q[1];
sx q[1];
rz(-0.64607611) q[1];
sx q[1];
rz(-2.0571713) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0823776) q[0];
sx q[0];
rz(-0.70967662) q[0];
sx q[0];
rz(-1.5060471) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83553548) q[2];
sx q[2];
rz(-0.37097574) q[2];
sx q[2];
rz(-2.6717466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.919405) q[1];
sx q[1];
rz(-2.0424941) q[1];
sx q[1];
rz(-0.081187831) q[1];
x q[2];
rz(-2.3929651) q[3];
sx q[3];
rz(-1.7065062) q[3];
sx q[3];
rz(-2.3627797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.083100975) q[2];
sx q[2];
rz(-0.19187555) q[2];
sx q[2];
rz(0.43711883) q[2];
rz(-1.7982177) q[3];
sx q[3];
rz(-1.3513887) q[3];
sx q[3];
rz(0.04537151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042353543) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(-1.7310671) q[0];
rz(2.9208185) q[1];
sx q[1];
rz(-2.6103554) q[1];
sx q[1];
rz(0.9332307) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8932874) q[0];
sx q[0];
rz(-2.3961005) q[0];
sx q[0];
rz(-3.1214691) q[0];
x q[1];
rz(1.2104697) q[2];
sx q[2];
rz(-0.98305741) q[2];
sx q[2];
rz(0.34229842) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9737967) q[1];
sx q[1];
rz(-0.78404155) q[1];
sx q[1];
rz(-1.7586437) q[1];
rz(-pi) q[2];
rz(-1.7228863) q[3];
sx q[3];
rz(-0.63816164) q[3];
sx q[3];
rz(1.9877951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.431388) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(-2.5101275) q[2];
rz(-0.7190052) q[3];
sx q[3];
rz(-1.7944261) q[3];
sx q[3];
rz(2.9001111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63856335) q[0];
sx q[0];
rz(-1.4894435) q[0];
sx q[0];
rz(1.3110934) q[0];
rz(-2.6999264) q[1];
sx q[1];
rz(-1.8001945) q[1];
sx q[1];
rz(1.1460907) q[1];
rz(-1.6497816) q[2];
sx q[2];
rz(-2.7597703) q[2];
sx q[2];
rz(-2.4673354) q[2];
rz(-0.56959116) q[3];
sx q[3];
rz(-1.6574331) q[3];
sx q[3];
rz(2.1059753) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

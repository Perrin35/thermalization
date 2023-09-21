OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8060057) q[0];
sx q[0];
rz(-0.94526362) q[0];
sx q[0];
rz(2.6160016) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(2.2367509) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.605699) q[0];
sx q[0];
rz(-2.4568757) q[0];
sx q[0];
rz(2.2936355) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4279847) q[2];
sx q[2];
rz(-2.8461694) q[2];
sx q[2];
rz(-1.3908536) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9368254) q[1];
sx q[1];
rz(-1.966404) q[1];
sx q[1];
rz(2.8863465) q[1];
rz(-pi) q[2];
rz(-2.9526887) q[3];
sx q[3];
rz(-1.8763262) q[3];
sx q[3];
rz(0.73959914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9378172) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(0.096244372) q[2];
rz(1.0359267) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(2.9878785) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44089833) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(2.3764215) q[0];
rz(1.8493429) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(-0.66295019) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63755858) q[0];
sx q[0];
rz(-0.090888977) q[0];
sx q[0];
rz(-0.81764098) q[0];
rz(-2.2234369) q[2];
sx q[2];
rz(-1.895972) q[2];
sx q[2];
rz(-0.52415372) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.029164974) q[1];
sx q[1];
rz(-1.6532835) q[1];
sx q[1];
rz(-0.44177456) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8421721) q[3];
sx q[3];
rz(-2.6357108) q[3];
sx q[3];
rz(-2.8477856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.092676) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(-0.32901397) q[2];
rz(2.4760903) q[3];
sx q[3];
rz(-2.9232959) q[3];
sx q[3];
rz(-1.8238508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.5927785) q[0];
sx q[0];
rz(-0.22286649) q[0];
sx q[0];
rz(2.9192525) q[0];
rz(-2.1242583) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(-0.51868784) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4043536) q[0];
sx q[0];
rz(-0.81439942) q[0];
sx q[0];
rz(-2.0703719) q[0];
rz(-pi) q[1];
rz(-0.2661744) q[2];
sx q[2];
rz(-1.6632348) q[2];
sx q[2];
rz(0.12649525) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31836244) q[1];
sx q[1];
rz(-1.3981817) q[1];
sx q[1];
rz(-0.77962064) q[1];
rz(1.5699925) q[3];
sx q[3];
rz(-1.5161848) q[3];
sx q[3];
rz(-0.87517232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2805933) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(-3.0730491) q[2];
rz(2.5391501) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(-0.13901916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.65072) q[0];
sx q[0];
rz(-0.77712178) q[0];
sx q[0];
rz(-2.9673476) q[0];
rz(2.6113367) q[1];
sx q[1];
rz(-1.6590051) q[1];
sx q[1];
rz(2.6285016) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4671191) q[0];
sx q[0];
rz(-1.6498483) q[0];
sx q[0];
rz(0.62069686) q[0];
rz(-0.60855234) q[2];
sx q[2];
rz(-0.90494472) q[2];
sx q[2];
rz(-0.37492875) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2345703) q[1];
sx q[1];
rz(-1.7090194) q[1];
sx q[1];
rz(0.45781086) q[1];
rz(-pi) q[2];
rz(3.1107535) q[3];
sx q[3];
rz(-1.3645002) q[3];
sx q[3];
rz(-1.4269958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6667368) q[2];
sx q[2];
rz(-0.36617827) q[2];
sx q[2];
rz(0.22988698) q[2];
rz(2.722548) q[3];
sx q[3];
rz(-1.7925526) q[3];
sx q[3];
rz(-2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6722365) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(-0.67681926) q[0];
rz(2.6485486) q[1];
sx q[1];
rz(-0.9489916) q[1];
sx q[1];
rz(-2.5255323) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84232932) q[0];
sx q[0];
rz(-0.063532524) q[0];
sx q[0];
rz(-0.97139831) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8370503) q[2];
sx q[2];
rz(-0.38123044) q[2];
sx q[2];
rz(0.234137) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1787781) q[1];
sx q[1];
rz(-1.5652834) q[1];
sx q[1];
rz(-0.57761901) q[1];
rz(-pi) q[2];
rz(2.2474399) q[3];
sx q[3];
rz(-1.3291306) q[3];
sx q[3];
rz(-1.579293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4074576) q[2];
sx q[2];
rz(-2.5871758) q[2];
sx q[2];
rz(-0.25536728) q[2];
rz(1.6051965) q[3];
sx q[3];
rz(-0.95241773) q[3];
sx q[3];
rz(-0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42246321) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(-2.6690924) q[0];
rz(-2.6155112) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(-2.2568259) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17381829) q[0];
sx q[0];
rz(-0.38823715) q[0];
sx q[0];
rz(-1.3740747) q[0];
rz(-0.21005819) q[2];
sx q[2];
rz(-1.5569599) q[2];
sx q[2];
rz(0.42523281) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5342418) q[1];
sx q[1];
rz(-2.631175) q[1];
sx q[1];
rz(2.1900602) q[1];
x q[2];
rz(1.0370449) q[3];
sx q[3];
rz(-1.3579206) q[3];
sx q[3];
rz(0.45149976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.74449599) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(2.8302144) q[2];
rz(-1.7729676) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41806528) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(0.061070651) q[0];
rz(-3.1014077) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(2.4087002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6717364) q[0];
sx q[0];
rz(-1.0438907) q[0];
sx q[0];
rz(2.0335474) q[0];
rz(2.6635025) q[2];
sx q[2];
rz(-0.92667898) q[2];
sx q[2];
rz(2.9113876) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1317132) q[1];
sx q[1];
rz(-2.1965346) q[1];
sx q[1];
rz(-0.0080607944) q[1];
rz(-pi) q[2];
rz(-0.041386889) q[3];
sx q[3];
rz(-1.6780516) q[3];
sx q[3];
rz(-2.0859352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0733033) q[2];
sx q[2];
rz(-1.1969593) q[2];
sx q[2];
rz(-0.51458365) q[2];
rz(1.2375281) q[3];
sx q[3];
rz(-2.5585744) q[3];
sx q[3];
rz(2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056203689) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(0.7094267) q[0];
rz(1.5052694) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(-2.8628796) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73662169) q[0];
sx q[0];
rz(-1.4355037) q[0];
sx q[0];
rz(-2.9297329) q[0];
rz(-0.3785554) q[2];
sx q[2];
rz(-0.59213973) q[2];
sx q[2];
rz(-2.2373667) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4019926) q[1];
sx q[1];
rz(-0.38312373) q[1];
sx q[1];
rz(2.1869704) q[1];
rz(-pi) q[2];
rz(1.7940984) q[3];
sx q[3];
rz(-1.6468862) q[3];
sx q[3];
rz(-0.014811024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8401106) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(-0.056079496) q[2];
rz(-2.2864443) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(2.7364031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1466325) q[0];
sx q[0];
rz(-0.17630795) q[0];
sx q[0];
rz(2.1726998) q[0];
rz(-0.47337198) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(1.999058) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56209598) q[0];
sx q[0];
rz(-1.2764494) q[0];
sx q[0];
rz(0.094898183) q[0];
x q[1];
rz(-0.70334401) q[2];
sx q[2];
rz(-1.5726657) q[2];
sx q[2];
rz(-1.4700996) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0513623) q[1];
sx q[1];
rz(-1.6962595) q[1];
sx q[1];
rz(2.5999703) q[1];
rz(-2.9386018) q[3];
sx q[3];
rz(-2.2309125) q[3];
sx q[3];
rz(-0.70860329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.896495) q[2];
sx q[2];
rz(-1.921804) q[2];
sx q[2];
rz(1.0207821) q[2];
rz(-0.3237237) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(-0.011172115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97994119) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(2.4401869) q[0];
rz(-2.2258863) q[1];
sx q[1];
rz(-1.0083102) q[1];
sx q[1];
rz(1.9030301) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68574821) q[0];
sx q[0];
rz(-1.0554753) q[0];
sx q[0];
rz(-0.2483764) q[0];
rz(0.752991) q[2];
sx q[2];
rz(-1.5333311) q[2];
sx q[2];
rz(-2.3245036) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4148256) q[1];
sx q[1];
rz(-1.6884202) q[1];
sx q[1];
rz(2.2284501) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.402461) q[3];
sx q[3];
rz(-0.59909648) q[3];
sx q[3];
rz(-1.7408016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.68676585) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(-0.043838538) q[2];
rz(1.2010126) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(-1.003585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6702406) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(3.1162221) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(0.92924835) q[2];
sx q[2];
rz(-2.6570448) q[2];
sx q[2];
rz(-1.6397283) q[2];
rz(0.79904859) q[3];
sx q[3];
rz(-1.6934762) q[3];
sx q[3];
rz(1.3140524) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

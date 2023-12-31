OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(0.84258643) q[0];
rz(-1.0358345) q[1];
sx q[1];
rz(-2.0422715) q[1];
sx q[1];
rz(1.6834747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9561477) q[0];
sx q[0];
rz(-1.7831793) q[0];
sx q[0];
rz(-2.3953715) q[0];
rz(-pi) q[1];
rz(-0.027878472) q[2];
sx q[2];
rz(-2.5275702) q[2];
sx q[2];
rz(-2.8369454) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.38274256) q[1];
sx q[1];
rz(-1.9470125) q[1];
sx q[1];
rz(1.43169) q[1];
rz(0.59431608) q[3];
sx q[3];
rz(-1.4694957) q[3];
sx q[3];
rz(0.017410226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6136916) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(-2.9620985) q[2];
rz(1.2256631) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(-2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3935788) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(-0.32546145) q[0];
rz(1.7851967) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(-1.9869841) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3852859) q[0];
sx q[0];
rz(-0.025408832) q[0];
sx q[0];
rz(2.4029762) q[0];
rz(-pi) q[1];
rz(2.1968368) q[2];
sx q[2];
rz(-1.2465887) q[2];
sx q[2];
rz(-2.7446483) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7696015) q[1];
sx q[1];
rz(-0.76906119) q[1];
sx q[1];
rz(3.0197057) q[1];
x q[2];
rz(2.9566544) q[3];
sx q[3];
rz(-1.8293081) q[3];
sx q[3];
rz(-2.1188494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4521728) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(2.2581805) q[2];
rz(2.6702821) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(-2.3538891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.31323355) q[0];
sx q[0];
rz(-1.4947083) q[0];
sx q[0];
rz(1.5154243) q[0];
rz(-0.60107636) q[1];
sx q[1];
rz(-2.5939012) q[1];
sx q[1];
rz(-1.0916969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0059144817) q[0];
sx q[0];
rz(-1.0251097) q[0];
sx q[0];
rz(-0.90555993) q[0];
rz(-pi) q[1];
rz(0.91471471) q[2];
sx q[2];
rz(-1.9064184) q[2];
sx q[2];
rz(-0.4682954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6669238) q[1];
sx q[1];
rz(-2.3033934) q[1];
sx q[1];
rz(2.073642) q[1];
x q[2];
rz(-1.5063498) q[3];
sx q[3];
rz(-1.6985053) q[3];
sx q[3];
rz(-2.8116022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(-2.2606405) q[2];
rz(-1.7679924) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(-1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3110733) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(0.4367035) q[0];
rz(2.9084335) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(0.31035796) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4219907) q[0];
sx q[0];
rz(-1.4388226) q[0];
sx q[0];
rz(2.7601526) q[0];
rz(1.6967625) q[2];
sx q[2];
rz(-0.8896041) q[2];
sx q[2];
rz(1.7484401) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6314108) q[1];
sx q[1];
rz(-1.5893755) q[1];
sx q[1];
rz(-1.2127962) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4471099) q[3];
sx q[3];
rz(-1.6061155) q[3];
sx q[3];
rz(-0.24231054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13005304) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(-1.0774353) q[2];
rz(-0.056190101) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(-1.594054) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524277) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(0.24965832) q[0];
rz(-1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(-2.2713984) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2244959) q[0];
sx q[0];
rz(-1.3670237) q[0];
sx q[0];
rz(-0.32075551) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9179847) q[2];
sx q[2];
rz(-0.70906559) q[2];
sx q[2];
rz(-0.86415926) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.732547) q[1];
sx q[1];
rz(-0.88102075) q[1];
sx q[1];
rz(-1.1613261) q[1];
rz(-2.1327444) q[3];
sx q[3];
rz(-1.7147439) q[3];
sx q[3];
rz(-0.89509237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.27328086) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-0.67374054) q[2];
rz(2.8379748) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7917787) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(-2.8836024) q[0];
rz(-2.7164283) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.4917096) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30116044) q[0];
sx q[0];
rz(-1.4727117) q[0];
sx q[0];
rz(0.43099404) q[0];
rz(-pi) q[1];
rz(2.4620373) q[2];
sx q[2];
rz(-1.2365885) q[2];
sx q[2];
rz(2.213775) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7266453) q[1];
sx q[1];
rz(-0.84268314) q[1];
sx q[1];
rz(2.3689518) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10156472) q[3];
sx q[3];
rz(-1.933681) q[3];
sx q[3];
rz(-0.18728072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.012718) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(-0.091726124) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(-2.2475524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3180852) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(-2.7303625) q[0];
rz(-2.2757018) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(0.033989865) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8695553) q[0];
sx q[0];
rz(-0.79622686) q[0];
sx q[0];
rz(1.7982593) q[0];
x q[1];
rz(2.2529644) q[2];
sx q[2];
rz(-0.76105984) q[2];
sx q[2];
rz(-1.6737446) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.781144) q[1];
sx q[1];
rz(-2.2478659) q[1];
sx q[1];
rz(-3.0767246) q[1];
rz(-pi) q[2];
rz(0.61492413) q[3];
sx q[3];
rz(-1.691754) q[3];
sx q[3];
rz(1.8254335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5380481) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(0.87654385) q[2];
rz(2.792568) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2440764) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(-2.7476655) q[0];
rz(2.774033) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(1.4454909) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60364265) q[0];
sx q[0];
rz(-1.1374439) q[0];
sx q[0];
rz(-3.135878) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94930737) q[2];
sx q[2];
rz(-2.0645803) q[2];
sx q[2];
rz(2.2311503) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.48965028) q[1];
sx q[1];
rz(-1.7282657) q[1];
sx q[1];
rz(-2.5850992) q[1];
rz(-pi) q[2];
rz(-0.94533841) q[3];
sx q[3];
rz(-1.006554) q[3];
sx q[3];
rz(2.4953147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8470856) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(-0.40714804) q[2];
rz(-1.5173222) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(-2.8919162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(0.80609926) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(1.8898213) q[0];
rz(0.66954008) q[1];
sx q[1];
rz(-1.957683) q[1];
sx q[1];
rz(-0.30977419) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4315902) q[0];
sx q[0];
rz(-0.54034034) q[0];
sx q[0];
rz(-2.8938328) q[0];
rz(-pi) q[1];
rz(-1.6236213) q[2];
sx q[2];
rz(-2.9644358) q[2];
sx q[2];
rz(-1.0747386) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.0032776912) q[1];
sx q[1];
rz(-1.365108) q[1];
sx q[1];
rz(0.059271952) q[1];
rz(-pi) q[2];
rz(-0.59659776) q[3];
sx q[3];
rz(-1.784424) q[3];
sx q[3];
rz(1.240977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70242515) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(1.2072198) q[2];
rz(-2.1045945) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(-2.4859378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0913775) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(1.2058831) q[0];
rz(2.5559015) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(1.4996128) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23696454) q[0];
sx q[0];
rz(-1.4693854) q[0];
sx q[0];
rz(-1.2786091) q[0];
rz(-1.502938) q[2];
sx q[2];
rz(-0.80791622) q[2];
sx q[2];
rz(2.9615336) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91883509) q[1];
sx q[1];
rz(-1.3593874) q[1];
sx q[1];
rz(-2.0445776) q[1];
x q[2];
rz(0.0025000574) q[3];
sx q[3];
rz(-1.7435939) q[3];
sx q[3];
rz(-2.1713184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2575834) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(-0.94669) q[2];
rz(-0.36869129) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(-2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7832227) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(0.070925698) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(2.3161841) q[2];
sx q[2];
rz(-2.5054629) q[2];
sx q[2];
rz(2.7873743) q[2];
rz(-2.4214217) q[3];
sx q[3];
rz(-1.9577033) q[3];
sx q[3];
rz(-2.4693558) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

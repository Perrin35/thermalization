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
rz(0.49864545) q[0];
sx q[0];
rz(-2.5957624) q[0];
sx q[0];
rz(-0.11830615) q[0];
rz(0.91953295) q[1];
sx q[1];
rz(4.4448648) q[1];
sx q[1];
rz(8.2104609) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81033731) q[0];
sx q[0];
rz(-0.48733586) q[0];
sx q[0];
rz(-2.0361855) q[0];
rz(-pi) q[1];
rz(0.6251752) q[2];
sx q[2];
rz(-1.0708276) q[2];
sx q[2];
rz(1.878405) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7711) q[1];
sx q[1];
rz(-0.95540651) q[1];
sx q[1];
rz(0.16754383) q[1];
rz(-pi) q[2];
rz(0.021539979) q[3];
sx q[3];
rz(-1.2976754) q[3];
sx q[3];
rz(2.5813951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.168657) q[2];
sx q[2];
rz(-2.0106222) q[2];
sx q[2];
rz(2.0449779) q[2];
rz(0.24122572) q[3];
sx q[3];
rz(-1.7719519) q[3];
sx q[3];
rz(-2.296804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089652561) q[0];
sx q[0];
rz(-1.0900499) q[0];
sx q[0];
rz(0.3845149) q[0];
rz(2.8334726) q[1];
sx q[1];
rz(-2.6607951) q[1];
sx q[1];
rz(-0.62517977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168361) q[0];
sx q[0];
rz(-0.28795469) q[0];
sx q[0];
rz(0.11174996) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4246419) q[2];
sx q[2];
rz(-1.4745108) q[2];
sx q[2];
rz(-1.7801628) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75698419) q[1];
sx q[1];
rz(-0.9238657) q[1];
sx q[1];
rz(1.6773083) q[1];
x q[2];
rz(-0.92722882) q[3];
sx q[3];
rz(-2.4146898) q[3];
sx q[3];
rz(-2.79799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.008931) q[2];
sx q[2];
rz(-0.93175685) q[2];
sx q[2];
rz(0.30161944) q[2];
rz(0.53145069) q[3];
sx q[3];
rz(-0.88574946) q[3];
sx q[3];
rz(-2.0187995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5136435) q[0];
sx q[0];
rz(-1.0114089) q[0];
sx q[0];
rz(-1.9568141) q[0];
rz(-2.9966677) q[1];
sx q[1];
rz(-1.2624319) q[1];
sx q[1];
rz(-1.5941934) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3959409) q[0];
sx q[0];
rz(-0.81123039) q[0];
sx q[0];
rz(0.70180362) q[0];
rz(-pi) q[1];
rz(1.1357665) q[2];
sx q[2];
rz(-2.7270881) q[2];
sx q[2];
rz(0.5629102) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.84128296) q[1];
sx q[1];
rz(-1.0307068) q[1];
sx q[1];
rz(2.6909091) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7361467) q[3];
sx q[3];
rz(-1.8687752) q[3];
sx q[3];
rz(-1.0007981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.80921119) q[2];
sx q[2];
rz(-0.96497649) q[2];
sx q[2];
rz(-1.3261718) q[2];
rz(2.2281846) q[3];
sx q[3];
rz(-1.2181506) q[3];
sx q[3];
rz(-0.53328812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42295414) q[0];
sx q[0];
rz(-0.47684968) q[0];
sx q[0];
rz(-1.780321) q[0];
rz(-2.4286229) q[1];
sx q[1];
rz(-1.1402592) q[1];
sx q[1];
rz(0.68351173) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28302971) q[0];
sx q[0];
rz(-0.57335317) q[0];
sx q[0];
rz(1.5893776) q[0];
rz(0.2095378) q[2];
sx q[2];
rz(-1.352013) q[2];
sx q[2];
rz(2.0444586) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.57505137) q[1];
sx q[1];
rz(-1.05541) q[1];
sx q[1];
rz(-2.383286) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24432791) q[3];
sx q[3];
rz(-1.4736946) q[3];
sx q[3];
rz(-1.8116784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6354562) q[2];
sx q[2];
rz(-0.41716245) q[2];
sx q[2];
rz(2.6551841) q[2];
rz(-0.5168612) q[3];
sx q[3];
rz(-1.9603399) q[3];
sx q[3];
rz(2.1080871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0559693) q[0];
sx q[0];
rz(-1.6317246) q[0];
sx q[0];
rz(-1.7686718) q[0];
rz(-1.4068475) q[1];
sx q[1];
rz(-2.5739248) q[1];
sx q[1];
rz(-0.47703201) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041056496) q[0];
sx q[0];
rz(-1.9874128) q[0];
sx q[0];
rz(0.21505298) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1978756) q[2];
sx q[2];
rz(-2.5260128) q[2];
sx q[2];
rz(2.1114388) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.20471599) q[1];
sx q[1];
rz(-0.5936247) q[1];
sx q[1];
rz(-1.0148263) q[1];
x q[2];
rz(0.84609109) q[3];
sx q[3];
rz(-2.2109004) q[3];
sx q[3];
rz(-0.37615955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10739747) q[2];
sx q[2];
rz(-0.98661462) q[2];
sx q[2];
rz(0.46727115) q[2];
rz(0.16921903) q[3];
sx q[3];
rz(-1.7305814) q[3];
sx q[3];
rz(-2.8362078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53448236) q[0];
sx q[0];
rz(-0.87962532) q[0];
sx q[0];
rz(2.95209) q[0];
rz(-0.27319187) q[1];
sx q[1];
rz(-1.0739645) q[1];
sx q[1];
rz(-2.3409519) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9439745) q[0];
sx q[0];
rz(-1.6139784) q[0];
sx q[0];
rz(-1.3373242) q[0];
rz(-0.34363644) q[2];
sx q[2];
rz(-0.26310194) q[2];
sx q[2];
rz(-3.0575036) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8756518) q[1];
sx q[1];
rz(-1.7663301) q[1];
sx q[1];
rz(0.13130782) q[1];
x q[2];
rz(-2.0728821) q[3];
sx q[3];
rz(-1.0866345) q[3];
sx q[3];
rz(-1.3997149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.93929401) q[2];
sx q[2];
rz(-2.9677128) q[2];
sx q[2];
rz(0.536971) q[2];
rz(1.8592853) q[3];
sx q[3];
rz(-1.8351646) q[3];
sx q[3];
rz(1.4793226) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.110431) q[0];
sx q[0];
rz(-0.81542504) q[0];
sx q[0];
rz(-1.331331) q[0];
rz(-1.7000465) q[1];
sx q[1];
rz(-1.4794289) q[1];
sx q[1];
rz(-1.9190681) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44935575) q[0];
sx q[0];
rz(-2.6550482) q[0];
sx q[0];
rz(-1.0445717) q[0];
x q[1];
rz(0.55180727) q[2];
sx q[2];
rz(-1.5071609) q[2];
sx q[2];
rz(-2.3752874) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4536087) q[1];
sx q[1];
rz(-2.2609432) q[1];
sx q[1];
rz(-2.860387) q[1];
rz(-pi) q[2];
rz(3.0223373) q[3];
sx q[3];
rz(-2.0125161) q[3];
sx q[3];
rz(0.91631232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20087251) q[2];
sx q[2];
rz(-0.49786374) q[2];
sx q[2];
rz(2.8087924) q[2];
rz(0.29916549) q[3];
sx q[3];
rz(-1.2406415) q[3];
sx q[3];
rz(-3.1201194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0126295) q[0];
sx q[0];
rz(-0.7681995) q[0];
sx q[0];
rz(-2.8666038) q[0];
rz(3.0601652) q[1];
sx q[1];
rz(-0.80137253) q[1];
sx q[1];
rz(-1.3224695) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24182651) q[0];
sx q[0];
rz(-1.6303008) q[0];
sx q[0];
rz(1.3391979) q[0];
rz(-pi) q[1];
rz(-2.7733388) q[2];
sx q[2];
rz(-0.91742951) q[2];
sx q[2];
rz(-2.7453932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15582196) q[1];
sx q[1];
rz(-2.3233674) q[1];
sx q[1];
rz(0.51795141) q[1];
rz(-0.72663088) q[3];
sx q[3];
rz(-1.0612349) q[3];
sx q[3];
rz(-1.0973615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72121173) q[2];
sx q[2];
rz(-0.12004852) q[2];
sx q[2];
rz(1.4018641) q[2];
rz(-1.8574572) q[3];
sx q[3];
rz(-1.1893585) q[3];
sx q[3];
rz(3.106626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5597124) q[0];
sx q[0];
rz(-1.548883) q[0];
sx q[0];
rz(2.6522719) q[0];
rz(0.020546546) q[1];
sx q[1];
rz(-1.2362044) q[1];
sx q[1];
rz(-0.70125854) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6416949) q[0];
sx q[0];
rz(-2.5144082) q[0];
sx q[0];
rz(1.1526965) q[0];
rz(1.6270774) q[2];
sx q[2];
rz(-0.70437925) q[2];
sx q[2];
rz(0.13408184) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.32966954) q[1];
sx q[1];
rz(-0.81336248) q[1];
sx q[1];
rz(1.2438891) q[1];
rz(-2.4550426) q[3];
sx q[3];
rz(-1.1995763) q[3];
sx q[3];
rz(0.29774259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7265085) q[2];
sx q[2];
rz(-2.5231611) q[2];
sx q[2];
rz(2.6611967) q[2];
rz(1.4091617) q[3];
sx q[3];
rz(-1.7879854) q[3];
sx q[3];
rz(0.33717808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3163863) q[0];
sx q[0];
rz(-1.2403064) q[0];
sx q[0];
rz(0.41148841) q[0];
rz(1.2443789) q[1];
sx q[1];
rz(-2.0020516) q[1];
sx q[1];
rz(0.32434514) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5748904) q[0];
sx q[0];
rz(-1.5181203) q[0];
sx q[0];
rz(-1.1979237) q[0];
rz(-pi) q[1];
rz(-1.2191992) q[2];
sx q[2];
rz(-1.2484387) q[2];
sx q[2];
rz(-0.11478648) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.56411769) q[1];
sx q[1];
rz(-1.0077969) q[1];
sx q[1];
rz(-0.098280829) q[1];
rz(-0.12121157) q[3];
sx q[3];
rz(-1.5577661) q[3];
sx q[3];
rz(2.7708294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1856498) q[2];
sx q[2];
rz(-1.1142542) q[2];
sx q[2];
rz(3.0246217) q[2];
rz(0.95514917) q[3];
sx q[3];
rz(-0.17894608) q[3];
sx q[3];
rz(-0.54168934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87880001) q[0];
sx q[0];
rz(-2.265082) q[0];
sx q[0];
rz(-2.2646917) q[0];
rz(2.0211438) q[1];
sx q[1];
rz(-2.3255377) q[1];
sx q[1];
rz(-1.9768523) q[1];
rz(1.2016313) q[2];
sx q[2];
rz(-0.537048) q[2];
sx q[2];
rz(-0.67759003) q[2];
rz(-1.6616115) q[3];
sx q[3];
rz(-2.869893) q[3];
sx q[3];
rz(0.13520959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

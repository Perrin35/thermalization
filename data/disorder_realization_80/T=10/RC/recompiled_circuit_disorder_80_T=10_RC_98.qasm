OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39188448) q[0];
sx q[0];
rz(-0.19667974) q[0];
sx q[0];
rz(-1.952202) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(-0.84140468) q[1];
sx q[1];
rz(0.37766159) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084378622) q[0];
sx q[0];
rz(-2.9998261) q[0];
sx q[0];
rz(2.111582) q[0];
rz(-pi) q[1];
rz(1.8544191) q[2];
sx q[2];
rz(-1.5986773) q[2];
sx q[2];
rz(1.6698128) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9610112) q[1];
sx q[1];
rz(-0.86997021) q[1];
sx q[1];
rz(0.58971528) q[1];
x q[2];
rz(-1.430106) q[3];
sx q[3];
rz(-1.9860387) q[3];
sx q[3];
rz(-1.5822441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1315786) q[2];
sx q[2];
rz(-2.6476314) q[2];
sx q[2];
rz(0.67260355) q[2];
rz(-0.16942313) q[3];
sx q[3];
rz(-0.38893458) q[3];
sx q[3];
rz(-1.3385564) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23068962) q[0];
sx q[0];
rz(-2.2407273) q[0];
sx q[0];
rz(0.22856523) q[0];
rz(-0.16054343) q[1];
sx q[1];
rz(-1.7030145) q[1];
sx q[1];
rz(2.8536318) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1624958) q[0];
sx q[0];
rz(-1.393035) q[0];
sx q[0];
rz(-1.3773247) q[0];
rz(2.7520913) q[2];
sx q[2];
rz(-1.0132388) q[2];
sx q[2];
rz(-0.26706375) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2634695) q[1];
sx q[1];
rz(-2.1681004) q[1];
sx q[1];
rz(-1.0695446) q[1];
x q[2];
rz(-3.037022) q[3];
sx q[3];
rz(-0.35211709) q[3];
sx q[3];
rz(2.5196688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.59445375) q[2];
sx q[2];
rz(-1.2524266) q[2];
sx q[2];
rz(1.4734369) q[2];
rz(-0.93747059) q[3];
sx q[3];
rz(-0.44527403) q[3];
sx q[3];
rz(0.37500769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8163452) q[0];
sx q[0];
rz(-0.49961093) q[0];
sx q[0];
rz(-0.26741272) q[0];
rz(1.7193517) q[1];
sx q[1];
rz(-1.1317252) q[1];
sx q[1];
rz(-0.95169383) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6355977) q[0];
sx q[0];
rz(-1.6615168) q[0];
sx q[0];
rz(-1.3784301) q[0];
rz(-pi) q[1];
rz(1.9556324) q[2];
sx q[2];
rz(-0.54883146) q[2];
sx q[2];
rz(-1.8897111) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0905076) q[1];
sx q[1];
rz(-1.2625492) q[1];
sx q[1];
rz(-2.3329263) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1943201) q[3];
sx q[3];
rz(-2.0200649) q[3];
sx q[3];
rz(-2.505213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13016985) q[2];
sx q[2];
rz(-1.6458076) q[2];
sx q[2];
rz(-2.7974131) q[2];
rz(-2.2551645) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(-0.56604958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-3.0043871) q[0];
sx q[0];
rz(-1.6442278) q[0];
sx q[0];
rz(-1.8970867) q[0];
rz(3.0124774) q[1];
sx q[1];
rz(-1.3543509) q[1];
sx q[1];
rz(-2.7688162) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6308206) q[0];
sx q[0];
rz(-1.6850867) q[0];
sx q[0];
rz(-0.75829102) q[0];
x q[1];
rz(-2.3635025) q[2];
sx q[2];
rz(-1.8847244) q[2];
sx q[2];
rz(2.2997466) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4200538) q[1];
sx q[1];
rz(-2.5905847) q[1];
sx q[1];
rz(2.0481471) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0471441) q[3];
sx q[3];
rz(-1.2762478) q[3];
sx q[3];
rz(2.0457207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5059775) q[2];
sx q[2];
rz(-1.6515235) q[2];
sx q[2];
rz(0.90488952) q[2];
rz(-2.7010226) q[3];
sx q[3];
rz(-0.4959271) q[3];
sx q[3];
rz(2.1499965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6624517) q[0];
sx q[0];
rz(-2.9859556) q[0];
sx q[0];
rz(1.2878081) q[0];
rz(-1.4783391) q[1];
sx q[1];
rz(-1.9961424) q[1];
sx q[1];
rz(2.6073661) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9984765) q[0];
sx q[0];
rz(-1.8792218) q[0];
sx q[0];
rz(-2.8240859) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.062837) q[2];
sx q[2];
rz(-1.8740219) q[2];
sx q[2];
rz(-2.6360896) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2248762) q[1];
sx q[1];
rz(-1.4955048) q[1];
sx q[1];
rz(-3.0776943) q[1];
rz(-pi) q[2];
rz(-1.4388496) q[3];
sx q[3];
rz(-0.85634106) q[3];
sx q[3];
rz(-0.42657846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6405032) q[2];
sx q[2];
rz(-1.9260294) q[2];
sx q[2];
rz(-0.14349288) q[2];
rz(-1.3714553) q[3];
sx q[3];
rz(-1.2797132) q[3];
sx q[3];
rz(0.21970704) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4390398) q[0];
sx q[0];
rz(-2.2655903) q[0];
sx q[0];
rz(0.23705661) q[0];
rz(1.9006231) q[1];
sx q[1];
rz(-0.82890141) q[1];
sx q[1];
rz(-1.3751078) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.188376) q[0];
sx q[0];
rz(-1.2326476) q[0];
sx q[0];
rz(-1.8494649) q[0];
x q[1];
rz(2.3926922) q[2];
sx q[2];
rz(-1.6237215) q[2];
sx q[2];
rz(-0.63268328) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68723893) q[1];
sx q[1];
rz(-0.41077405) q[1];
sx q[1];
rz(0.06128581) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49075134) q[3];
sx q[3];
rz(-1.9487582) q[3];
sx q[3];
rz(-0.85765391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53283006) q[2];
sx q[2];
rz(-0.71981788) q[2];
sx q[2];
rz(-2.9928845) q[2];
rz(0.016629774) q[3];
sx q[3];
rz(-2.7745268) q[3];
sx q[3];
rz(-0.19255157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49611133) q[0];
sx q[0];
rz(-2.8493024) q[0];
sx q[0];
rz(2.2684229) q[0];
rz(0.9219777) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(0.08392863) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31939313) q[0];
sx q[0];
rz(-2.5376352) q[0];
sx q[0];
rz(-1.970406) q[0];
rz(-pi) q[1];
rz(0.35328816) q[2];
sx q[2];
rz(-1.7751667) q[2];
sx q[2];
rz(2.3452961) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1740239) q[1];
sx q[1];
rz(-2.1597383) q[1];
sx q[1];
rz(-2.3322361) q[1];
x q[2];
rz(1.6478959) q[3];
sx q[3];
rz(-0.73148433) q[3];
sx q[3];
rz(0.25792083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38368791) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(2.596358) q[2];
rz(2.7111354) q[3];
sx q[3];
rz(-1.1377708) q[3];
sx q[3];
rz(0.30495131) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6298744) q[0];
sx q[0];
rz(-3.1261303) q[0];
sx q[0];
rz(1.9301201) q[0];
rz(-2.1684872) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(-2.1957695) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66447542) q[0];
sx q[0];
rz(-0.23453377) q[0];
sx q[0];
rz(2.7995336) q[0];
rz(2.2875167) q[2];
sx q[2];
rz(-1.0208566) q[2];
sx q[2];
rz(2.9727109) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.91499) q[1];
sx q[1];
rz(-1.1922925) q[1];
sx q[1];
rz(1.6016017) q[1];
rz(-pi) q[2];
rz(1.5628912) q[3];
sx q[3];
rz(-2.2369011) q[3];
sx q[3];
rz(2.283309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2716081) q[2];
sx q[2];
rz(-1.4166069) q[2];
sx q[2];
rz(-0.28042173) q[2];
rz(2.5366606) q[3];
sx q[3];
rz(-2.0623902) q[3];
sx q[3];
rz(0.98208565) q[3];
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
rz(-2.9873001) q[0];
sx q[0];
rz(-0.91006088) q[0];
sx q[0];
rz(2.5262685) q[0];
rz(-0.92957169) q[1];
sx q[1];
rz(-0.85837448) q[1];
sx q[1];
rz(0.5756793) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0711813) q[0];
sx q[0];
rz(-1.0072664) q[0];
sx q[0];
rz(-1.0067183) q[0];
rz(-pi) q[1];
rz(1.3747146) q[2];
sx q[2];
rz(-1.6819281) q[2];
sx q[2];
rz(0.17734222) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5513735) q[1];
sx q[1];
rz(-3.0312523) q[1];
sx q[1];
rz(1.8054086) q[1];
x q[2];
rz(0.18746312) q[3];
sx q[3];
rz(-2.1037357) q[3];
sx q[3];
rz(-0.77123469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3118887) q[2];
sx q[2];
rz(-2.3788033) q[2];
sx q[2];
rz(3.1196307) q[2];
rz(-0.17045505) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(-0.26486614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4158674) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(-0.6341933) q[0];
rz(3.1126853) q[1];
sx q[1];
rz(-2.3529265) q[1];
sx q[1];
rz(0.96910563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4983738) q[0];
sx q[0];
rz(-1.6341097) q[0];
sx q[0];
rz(-1.4988585) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74659851) q[2];
sx q[2];
rz(-1.714163) q[2];
sx q[2];
rz(1.5310841) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.22419588) q[1];
sx q[1];
rz(-2.2317413) q[1];
sx q[1];
rz(1.7978653) q[1];
x q[2];
rz(1.85384) q[3];
sx q[3];
rz(-1.1453562) q[3];
sx q[3];
rz(-0.89983672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6897631) q[2];
sx q[2];
rz(-1.4844866) q[2];
sx q[2];
rz(-2.7815212) q[2];
rz(-1.5277956) q[3];
sx q[3];
rz(-0.66643047) q[3];
sx q[3];
rz(0.56267363) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9344899) q[0];
sx q[0];
rz(-1.5705382) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(-0.044152505) q[1];
sx q[1];
rz(-1.4587198) q[1];
sx q[1];
rz(-1.1062467) q[1];
rz(-1.953223) q[2];
sx q[2];
rz(-1.8845176) q[2];
sx q[2];
rz(0.11558576) q[2];
rz(-2.6882761) q[3];
sx q[3];
rz(-1.286187) q[3];
sx q[3];
rz(-1.5766531) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

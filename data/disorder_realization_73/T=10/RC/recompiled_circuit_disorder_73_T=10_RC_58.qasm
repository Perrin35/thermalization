OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1459382) q[0];
sx q[0];
rz(3.6448195) q[0];
sx q[0];
rz(10.148944) q[0];
rz(6.9231482) q[1];
sx q[1];
rz(5.7531113) q[1];
sx q[1];
rz(2.35676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1278348) q[0];
sx q[0];
rz(-2.7779967) q[0];
sx q[0];
rz(2.512393) q[0];
rz(2.9169693) q[2];
sx q[2];
rz(-2.7135239) q[2];
sx q[2];
rz(0.12878865) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8514429) q[1];
sx q[1];
rz(-0.57934299) q[1];
sx q[1];
rz(-1.1555374) q[1];
x q[2];
rz(1.1562528) q[3];
sx q[3];
rz(-1.6016377) q[3];
sx q[3];
rz(-0.24658345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5518387) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(-0.067967728) q[2];
rz(0.12456482) q[3];
sx q[3];
rz(-2.8187276) q[3];
sx q[3];
rz(-1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2215866) q[0];
sx q[0];
rz(-0.13555549) q[0];
sx q[0];
rz(0.24366972) q[0];
rz(-2.5098353) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(-1.7858645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60123721) q[0];
sx q[0];
rz(-1.8260801) q[0];
sx q[0];
rz(3.113494) q[0];
rz(-pi) q[1];
rz(-0.83265702) q[2];
sx q[2];
rz(-2.6359733) q[2];
sx q[2];
rz(-0.51072272) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1214952) q[1];
sx q[1];
rz(-1.1250245) q[1];
sx q[1];
rz(0.14264588) q[1];
rz(-2.5494266) q[3];
sx q[3];
rz(-1.5321931) q[3];
sx q[3];
rz(-1.6174699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0624861) q[2];
sx q[2];
rz(-0.9920384) q[2];
sx q[2];
rz(2.8919354) q[2];
rz(0.50659531) q[3];
sx q[3];
rz(-1.6258312) q[3];
sx q[3];
rz(0.33199582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8963985) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(-2.2431592) q[0];
rz(1.3348745) q[1];
sx q[1];
rz(-1.2355665) q[1];
sx q[1];
rz(1.867884) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1153591) q[0];
sx q[0];
rz(-2.6991978) q[0];
sx q[0];
rz(0.94985234) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31295915) q[2];
sx q[2];
rz(-0.3096146) q[2];
sx q[2];
rz(1.5086053) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0395567) q[1];
sx q[1];
rz(-1.2219056) q[1];
sx q[1];
rz(0.7015014) q[1];
rz(-pi) q[2];
rz(1.0728733) q[3];
sx q[3];
rz(-2.2306799) q[3];
sx q[3];
rz(2.4592196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0818103) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(0.95345062) q[2];
rz(0.034514286) q[3];
sx q[3];
rz(-0.78648609) q[3];
sx q[3];
rz(0.22687337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2801441) q[0];
sx q[0];
rz(-2.9252958) q[0];
sx q[0];
rz(-2.8934073) q[0];
rz(2.10363) q[1];
sx q[1];
rz(-1.1231517) q[1];
sx q[1];
rz(-0.074137069) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.733487) q[0];
sx q[0];
rz(-1.8251849) q[0];
sx q[0];
rz(2.0612201) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2494227) q[2];
sx q[2];
rz(-1.8461508) q[2];
sx q[2];
rz(-0.25472578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6756235) q[1];
sx q[1];
rz(-0.33756653) q[1];
sx q[1];
rz(1.9338495) q[1];
rz(1.3965963) q[3];
sx q[3];
rz(-1.3341691) q[3];
sx q[3];
rz(-0.06121204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8903824) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(3.1029491) q[2];
rz(-0.97366992) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(-0.27004778) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6073109) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(-1.3624396) q[0];
rz(0.81659395) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(-1.978925) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29279941) q[0];
sx q[0];
rz(-1.6822364) q[0];
sx q[0];
rz(0.029990002) q[0];
rz(-pi) q[1];
rz(-2.0432203) q[2];
sx q[2];
rz(-0.35944164) q[2];
sx q[2];
rz(3.0430832) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7971024) q[1];
sx q[1];
rz(-1.7624197) q[1];
sx q[1];
rz(1.2002798) q[1];
x q[2];
rz(-2.0883457) q[3];
sx q[3];
rz(-1.9540817) q[3];
sx q[3];
rz(-2.7305207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6267307) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(-0.14222063) q[2];
rz(-2.2375315) q[3];
sx q[3];
rz(-1.8217434) q[3];
sx q[3];
rz(-2.952125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.11480039) q[0];
sx q[0];
rz(-0.27652201) q[0];
sx q[0];
rz(-1.4676771) q[0];
rz(0.57178512) q[1];
sx q[1];
rz(-2.7829058) q[1];
sx q[1];
rz(2.8335559) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7361569) q[0];
sx q[0];
rz(-1.476256) q[0];
sx q[0];
rz(-3.0148538) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61200895) q[2];
sx q[2];
rz(-1.0527805) q[2];
sx q[2];
rz(-2.3305364) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.249237) q[1];
sx q[1];
rz(-2.4294469) q[1];
sx q[1];
rz(-1.874079) q[1];
rz(1.6623392) q[3];
sx q[3];
rz(-1.2237826) q[3];
sx q[3];
rz(2.7029944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77928153) q[2];
sx q[2];
rz(-1.4004536) q[2];
sx q[2];
rz(1.1479088) q[2];
rz(-2.4273196) q[3];
sx q[3];
rz(-2.2439984) q[3];
sx q[3];
rz(-0.64546293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66184735) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(0.1299783) q[0];
rz(0.030844363) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(-2.470509) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6993461) q[0];
sx q[0];
rz(-1.5171577) q[0];
sx q[0];
rz(-1.6058812) q[0];
rz(-1.7550049) q[2];
sx q[2];
rz(-1.6008018) q[2];
sx q[2];
rz(-1.9629994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8819067) q[1];
sx q[1];
rz(-0.38609186) q[1];
sx q[1];
rz(1.2950456) q[1];
rz(-pi) q[2];
rz(0.047366553) q[3];
sx q[3];
rz(-2.2455375) q[3];
sx q[3];
rz(2.588152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0188296) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(-2.5637131) q[2];
rz(-0.028586483) q[3];
sx q[3];
rz(-1.8609906) q[3];
sx q[3];
rz(1.2602497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1639444) q[0];
sx q[0];
rz(-2.2387235) q[0];
sx q[0];
rz(0.41241616) q[0];
rz(1.6917797) q[1];
sx q[1];
rz(-1.342536) q[1];
sx q[1];
rz(1.1669881) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0359356) q[0];
sx q[0];
rz(-1.3577537) q[0];
sx q[0];
rz(0.94353326) q[0];
x q[1];
rz(-1.251986) q[2];
sx q[2];
rz(-2.406771) q[2];
sx q[2];
rz(-2.1679945) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1046909) q[1];
sx q[1];
rz(-0.19100405) q[1];
sx q[1];
rz(-2.9810993) q[1];
rz(-pi) q[2];
rz(-1.8230209) q[3];
sx q[3];
rz(-0.54402292) q[3];
sx q[3];
rz(-2.9094484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2010487) q[2];
sx q[2];
rz(-0.87934914) q[2];
sx q[2];
rz(0.18903014) q[2];
rz(-0.14686251) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(-1.7485025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015633164) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(2.2145859) q[0];
rz(-1.3828297) q[1];
sx q[1];
rz(-0.60950509) q[1];
sx q[1];
rz(-1.6519201) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7324069) q[0];
sx q[0];
rz(-0.86940765) q[0];
sx q[0];
rz(-0.62774815) q[0];
rz(-pi) q[1];
rz(-1.6696879) q[2];
sx q[2];
rz(-0.47404587) q[2];
sx q[2];
rz(2.3941819) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5786963) q[1];
sx q[1];
rz(-1.7637196) q[1];
sx q[1];
rz(1.7460515) q[1];
rz(2.1065815) q[3];
sx q[3];
rz(-2.0097369) q[3];
sx q[3];
rz(2.5715695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5902517) q[2];
sx q[2];
rz(-0.44289032) q[2];
sx q[2];
rz(-1.7112973) q[2];
rz(-2.5643505) q[3];
sx q[3];
rz(-2.2539299) q[3];
sx q[3];
rz(-1.757471) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5883314) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(-2.8531895) q[0];
rz(0.53238955) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(0.14702252) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08976905) q[0];
sx q[0];
rz(-2.4252486) q[0];
sx q[0];
rz(-0.96869529) q[0];
rz(-pi) q[1];
rz(-0.27980079) q[2];
sx q[2];
rz(-0.24926148) q[2];
sx q[2];
rz(-1.3485497) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.30565572) q[1];
sx q[1];
rz(-2.5218997) q[1];
sx q[1];
rz(1.5020919) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9045194) q[3];
sx q[3];
rz(-0.58963886) q[3];
sx q[3];
rz(1.8001363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0816575) q[2];
sx q[2];
rz(-2.7066878) q[2];
sx q[2];
rz(-2.3975513) q[2];
rz(2.384281) q[3];
sx q[3];
rz(-1.7777187) q[3];
sx q[3];
rz(-1.7448759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.1159146) q[0];
sx q[0];
rz(-1.0703351) q[0];
sx q[0];
rz(-1.0967789) q[0];
rz(-0.81746447) q[1];
sx q[1];
rz(-1.2066963) q[1];
sx q[1];
rz(-0.6304601) q[1];
rz(-0.4521162) q[2];
sx q[2];
rz(-1.5099667) q[2];
sx q[2];
rz(-2.920334) q[2];
rz(0.13027262) q[3];
sx q[3];
rz(-2.1289005) q[3];
sx q[3];
rz(2.0990513) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

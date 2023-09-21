OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9956545) q[0];
sx q[0];
rz(-0.50322682) q[0];
sx q[0];
rz(2.4174262) q[0];
rz(-2.5016298) q[1];
sx q[1];
rz(-2.6115186) q[1];
sx q[1];
rz(-2.35676) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0137579) q[0];
sx q[0];
rz(-0.36359596) q[0];
sx q[0];
rz(0.6291997) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41854026) q[2];
sx q[2];
rz(-1.4782018) q[2];
sx q[2];
rz(1.6469524) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.29014978) q[1];
sx q[1];
rz(-2.5622497) q[1];
sx q[1];
rz(-1.1555374) q[1];
rz(-pi) q[2];
rz(-1.1562528) q[3];
sx q[3];
rz(-1.539955) q[3];
sx q[3];
rz(-0.24658345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5518387) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(-0.067967728) q[2];
rz(0.12456482) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(-1.7547866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92000604) q[0];
sx q[0];
rz(-0.13555549) q[0];
sx q[0];
rz(0.24366972) q[0];
rz(-2.5098353) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(1.3557281) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60123721) q[0];
sx q[0];
rz(-1.8260801) q[0];
sx q[0];
rz(3.113494) q[0];
x q[1];
rz(-0.83265702) q[2];
sx q[2];
rz(-0.50561935) q[2];
sx q[2];
rz(-2.6308699) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.3416459) q[1];
sx q[1];
rz(-0.46657545) q[1];
sx q[1];
rz(1.2816309) q[1];
x q[2];
rz(3.0725067) q[3];
sx q[3];
rz(-0.59327263) q[3];
sx q[3];
rz(3.1309576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0791066) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(-0.24965723) q[2];
rz(-0.50659531) q[3];
sx q[3];
rz(-1.6258312) q[3];
sx q[3];
rz(2.8095968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24519414) q[0];
sx q[0];
rz(-1.9165374) q[0];
sx q[0];
rz(-0.89843345) q[0];
rz(-1.3348745) q[1];
sx q[1];
rz(-1.2355665) q[1];
sx q[1];
rz(1.2737087) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1184517) q[0];
sx q[0];
rz(-1.8225192) q[0];
sx q[0];
rz(1.203042) q[0];
rz(2.8286335) q[2];
sx q[2];
rz(-0.3096146) q[2];
sx q[2];
rz(-1.6329873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.102036) q[1];
sx q[1];
rz(-1.9196871) q[1];
sx q[1];
rz(-0.7015014) q[1];
rz(-pi) q[2];
rz(0.7234296) q[3];
sx q[3];
rz(-1.1838786) q[3];
sx q[3];
rz(2.5748411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0818103) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(2.188142) q[2];
rz(-0.034514286) q[3];
sx q[3];
rz(-2.3551066) q[3];
sx q[3];
rz(-2.9147193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2801441) q[0];
sx q[0];
rz(-2.9252958) q[0];
sx q[0];
rz(2.8934073) q[0];
rz(-1.0379627) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(-3.0674556) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4081057) q[0];
sx q[0];
rz(-1.8251849) q[0];
sx q[0];
rz(-1.0803726) q[0];
x q[1];
rz(0.89216994) q[2];
sx q[2];
rz(-1.8461508) q[2];
sx q[2];
rz(0.25472578) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6756235) q[1];
sx q[1];
rz(-2.8040261) q[1];
sx q[1];
rz(1.2077431) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7449964) q[3];
sx q[3];
rz(-1.3341691) q[3];
sx q[3];
rz(3.0803806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2512102) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(-0.038643535) q[2];
rz(2.1679227) q[3];
sx q[3];
rz(-2.645851) q[3];
sx q[3];
rz(-2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6073109) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(1.3624396) q[0];
rz(-2.3249987) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(1.1626676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29279941) q[0];
sx q[0];
rz(-1.6822364) q[0];
sx q[0];
rz(3.1116027) q[0];
x q[1];
rz(1.0983724) q[2];
sx q[2];
rz(-2.782151) q[2];
sx q[2];
rz(0.098509468) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77046493) q[1];
sx q[1];
rz(-0.41509291) q[1];
sx q[1];
rz(-1.0789372) q[1];
rz(0.43443067) q[3];
sx q[3];
rz(-2.0475004) q[3];
sx q[3];
rz(-1.3694976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5148619) q[2];
sx q[2];
rz(-1.0802439) q[2];
sx q[2];
rz(-2.999372) q[2];
rz(0.90406117) q[3];
sx q[3];
rz(-1.8217434) q[3];
sx q[3];
rz(0.18946762) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0267923) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(1.6739155) q[0];
rz(2.5698075) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(-0.30803672) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.98826) q[0];
sx q[0];
rz(-1.6969661) q[0];
sx q[0];
rz(1.6660965) q[0];
rz(-pi) q[1];
rz(2.1790444) q[2];
sx q[2];
rz(-2.0934009) q[2];
sx q[2];
rz(-2.7163497) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.249237) q[1];
sx q[1];
rz(-2.4294469) q[1];
sx q[1];
rz(1.2675136) q[1];
rz(-pi) q[2];
rz(0.24758731) q[3];
sx q[3];
rz(-2.7831804) q[3];
sx q[3];
rz(-0.7022411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77928153) q[2];
sx q[2];
rz(-1.7411391) q[2];
sx q[2];
rz(1.1479088) q[2];
rz(0.71427304) q[3];
sx q[3];
rz(-2.2439984) q[3];
sx q[3];
rz(-0.64546293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66184735) q[0];
sx q[0];
rz(-0.84091887) q[0];
sx q[0];
rz(0.1299783) q[0];
rz(0.030844363) q[1];
sx q[1];
rz(-1.2896616) q[1];
sx q[1];
rz(-0.67108363) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13043159) q[0];
sx q[0];
rz(-1.6058308) q[0];
sx q[0];
rz(0.053671562) q[0];
rz(1.3865878) q[2];
sx q[2];
rz(-1.5407908) q[2];
sx q[2];
rz(1.9629994) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1783501) q[1];
sx q[1];
rz(-1.9415783) q[1];
sx q[1];
rz(3.0313655) q[1];
x q[2];
rz(1.5116793) q[3];
sx q[3];
rz(-2.4654508) q[3];
sx q[3];
rz(-0.47770559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.122763) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(0.57787952) q[2];
rz(-0.028586483) q[3];
sx q[3];
rz(-1.280602) q[3];
sx q[3];
rz(-1.2602497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9776483) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(-0.41241616) q[0];
rz(1.6917797) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(-1.1669881) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1056571) q[0];
sx q[0];
rz(-1.783839) q[0];
sx q[0];
rz(2.1980594) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27600482) q[2];
sx q[2];
rz(-0.88062421) q[2];
sx q[2];
rz(1.3921757) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.941277) q[1];
sx q[1];
rz(-1.759316) q[1];
sx q[1];
rz(-1.6016866) q[1];
rz(-pi) q[2];
rz(1.8230209) q[3];
sx q[3];
rz(-2.5975697) q[3];
sx q[3];
rz(-2.9094484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.940544) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(-0.18903014) q[2];
rz(2.9947301) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(1.3930901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-0.015633164) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(0.92700672) q[0];
rz(1.758763) q[1];
sx q[1];
rz(-0.60950509) q[1];
sx q[1];
rz(1.4896726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4179174) q[0];
sx q[0];
rz(-1.105504) q[0];
sx q[0];
rz(-2.377541) q[0];
rz(-pi) q[1];
rz(0.050612014) q[2];
sx q[2];
rz(-2.0423371) q[2];
sx q[2];
rz(0.6363578) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97396321) q[1];
sx q[1];
rz(-1.7427674) q[1];
sx q[1];
rz(0.19584882) q[1];
rz(-pi) q[2];
rz(-2.1065815) q[3];
sx q[3];
rz(-2.0097369) q[3];
sx q[3];
rz(0.5700232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5902517) q[2];
sx q[2];
rz(-0.44289032) q[2];
sx q[2];
rz(-1.7112973) q[2];
rz(-0.57724214) q[3];
sx q[3];
rz(-2.2539299) q[3];
sx q[3];
rz(-1.3841217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5532613) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(-0.28840315) q[0];
rz(0.53238955) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(-2.9945701) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1823746) q[0];
sx q[0];
rz(-1.9518513) q[0];
sx q[0];
rz(0.94840886) q[0];
rz(-1.5006127) q[2];
sx q[2];
rz(-1.3314221) q[2];
sx q[2];
rz(2.0812876) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3210996) q[1];
sx q[1];
rz(-1.610678) q[1];
sx q[1];
rz(0.95221968) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5649928) q[3];
sx q[3];
rz(-1.4398265) q[3];
sx q[3];
rz(3.1104345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0599351) q[2];
sx q[2];
rz(-2.7066878) q[2];
sx q[2];
rz(-0.74404136) q[2];
rz(0.75731164) q[3];
sx q[3];
rz(-1.7777187) q[3];
sx q[3];
rz(1.7448759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.025678) q[0];
sx q[0];
rz(-1.0703351) q[0];
sx q[0];
rz(-1.0967789) q[0];
rz(-2.3241282) q[1];
sx q[1];
rz(-1.9348963) q[1];
sx q[1];
rz(2.5111326) q[1];
rz(-0.4521162) q[2];
sx q[2];
rz(-1.5099667) q[2];
sx q[2];
rz(-2.920334) q[2];
rz(-1.7759454) q[3];
sx q[3];
rz(-2.57006) q[3];
sx q[3];
rz(-0.80001696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

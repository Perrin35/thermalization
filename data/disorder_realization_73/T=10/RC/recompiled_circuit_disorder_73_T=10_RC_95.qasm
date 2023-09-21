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
rz(-0.72416645) q[0];
rz(0.63996285) q[1];
sx q[1];
rz(-0.53007403) q[1];
sx q[1];
rz(-0.78483265) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1278348) q[0];
sx q[0];
rz(-2.7779967) q[0];
sx q[0];
rz(-2.512393) q[0];
rz(-pi) q[1];
rz(-1.6720812) q[2];
sx q[2];
rz(-1.9874319) q[2];
sx q[2];
rz(-3.0243304) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8514429) q[1];
sx q[1];
rz(-0.57934299) q[1];
sx q[1];
rz(1.1555374) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1562528) q[3];
sx q[3];
rz(-1.6016377) q[3];
sx q[3];
rz(2.8950092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5518387) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(-3.0736249) q[2];
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
rz(-0.92000604) q[0];
sx q[0];
rz(-3.0060372) q[0];
sx q[0];
rz(-2.8979229) q[0];
rz(2.5098353) q[1];
sx q[1];
rz(-1.4032204) q[1];
sx q[1];
rz(-1.7858645) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49039098) q[0];
sx q[0];
rz(-0.25679195) q[0];
sx q[0];
rz(-1.6780361) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83265702) q[2];
sx q[2];
rz(-2.6359733) q[2];
sx q[2];
rz(-2.6308699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6527378) q[1];
sx q[1];
rz(-1.4421717) q[1];
sx q[1];
rz(-2.0205523) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5494266) q[3];
sx q[3];
rz(-1.5321931) q[3];
sx q[3];
rz(-1.6174699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0624861) q[2];
sx q[2];
rz(-0.9920384) q[2];
sx q[2];
rz(-0.24965723) q[2];
rz(-2.6349973) q[3];
sx q[3];
rz(-1.5157615) q[3];
sx q[3];
rz(2.8095968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1153591) q[0];
sx q[0];
rz(-2.6991978) q[0];
sx q[0];
rz(-2.1917403) q[0];
rz(-pi) q[1];
rz(1.6689698) q[2];
sx q[2];
rz(-1.86491) q[2];
sx q[2];
rz(1.3054747) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.85325235) q[1];
sx q[1];
rz(-2.3715092) q[1];
sx q[1];
rz(2.6283162) q[1];
rz(-pi) q[2];
rz(-2.5898315) q[3];
sx q[3];
rz(-0.80358395) q[3];
sx q[3];
rz(1.4078275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0818103) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(-2.188142) q[2];
rz(-0.034514286) q[3];
sx q[3];
rz(-2.3551066) q[3];
sx q[3];
rz(-2.9147193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8614486) q[0];
sx q[0];
rz(-0.21629688) q[0];
sx q[0];
rz(0.24818534) q[0];
rz(1.0379627) q[1];
sx q[1];
rz(-1.1231517) q[1];
sx q[1];
rz(0.074137069) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0291245) q[0];
sx q[0];
rz(-1.0974786) q[0];
sx q[0];
rz(-0.28664696) q[0];
rz(-pi) q[1];
rz(1.147869) q[2];
sx q[2];
rz(-0.72407702) q[2];
sx q[2];
rz(-2.1507182) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3809507) q[1];
sx q[1];
rz(-1.4529072) q[1];
sx q[1];
rz(-1.2537434) q[1];
x q[2];
rz(-1.3965963) q[3];
sx q[3];
rz(-1.8074236) q[3];
sx q[3];
rz(3.0803806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2512102) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(-3.1029491) q[2];
rz(-0.97366992) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(-0.27004778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.6073109) q[0];
sx q[0];
rz(-1.5357635) q[0];
sx q[0];
rz(-1.3624396) q[0];
rz(0.81659395) q[1];
sx q[1];
rz(-1.8530308) q[1];
sx q[1];
rz(1.978925) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2813331) q[0];
sx q[0];
rz(-1.6006002) q[0];
sx q[0];
rz(1.6822862) q[0];
rz(-pi) q[1];
rz(-0.16935279) q[2];
sx q[2];
rz(-1.2522109) q[2];
sx q[2];
rz(2.7404075) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3711277) q[1];
sx q[1];
rz(-2.7264997) q[1];
sx q[1];
rz(1.0789372) q[1];
rz(-pi) q[2];
rz(-0.43443067) q[3];
sx q[3];
rz(-1.0940922) q[3];
sx q[3];
rz(1.7720951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5148619) q[2];
sx q[2];
rz(-1.0802439) q[2];
sx q[2];
rz(0.14222063) q[2];
rz(-2.2375315) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(2.952125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.0267923) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(-1.6739155) q[0];
rz(0.57178512) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(0.30803672) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15333262) q[0];
sx q[0];
rz(-1.6969661) q[0];
sx q[0];
rz(1.6660965) q[0];
rz(0.96254827) q[2];
sx q[2];
rz(-2.0934009) q[2];
sx q[2];
rz(2.7163497) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8923556) q[1];
sx q[1];
rz(-2.4294469) q[1];
sx q[1];
rz(1.2675136) q[1];
rz(-2.7932348) q[3];
sx q[3];
rz(-1.4847241) q[3];
sx q[3];
rz(-1.100988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3623111) q[2];
sx q[2];
rz(-1.7411391) q[2];
sx q[2];
rz(-1.9936838) q[2];
rz(2.4273196) q[3];
sx q[3];
rz(-0.89759421) q[3];
sx q[3];
rz(2.4961297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66184735) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(3.0116144) q[0];
rz(0.030844363) q[1];
sx q[1];
rz(-1.2896616) q[1];
sx q[1];
rz(2.470509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0111611) q[0];
sx q[0];
rz(-1.5357619) q[0];
sx q[0];
rz(-3.0879211) q[0];
rz(-pi) q[1];
rz(-0.030521557) q[2];
sx q[2];
rz(-1.754921) q[2];
sx q[2];
rz(0.39779278) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.259686) q[1];
sx q[1];
rz(-0.38609186) q[1];
sx q[1];
rz(1.2950456) q[1];
rz(-pi) q[2];
rz(-0.89550771) q[3];
sx q[3];
rz(-1.607778) q[3];
sx q[3];
rz(2.0946338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0188296) q[2];
sx q[2];
rz(-1.5315703) q[2];
sx q[2];
rz(2.5637131) q[2];
rz(0.028586483) q[3];
sx q[3];
rz(-1.280602) q[3];
sx q[3];
rz(-1.8813429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9776483) q[0];
sx q[0];
rz(-2.2387235) q[0];
sx q[0];
rz(-0.41241616) q[0];
rz(1.6917797) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(-1.1669881) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3827688) q[0];
sx q[0];
rz(-0.95982691) q[0];
sx q[0];
rz(0.2610892) q[0];
x q[1];
rz(-0.27600482) q[2];
sx q[2];
rz(-0.88062421) q[2];
sx q[2];
rz(1.7494169) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3762714) q[1];
sx q[1];
rz(-1.6011392) q[1];
sx q[1];
rz(-0.18860753) q[1];
rz(2.9917631) q[3];
sx q[3];
rz(-1.0458046) q[3];
sx q[3];
rz(2.6168952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.940544) q[2];
sx q[2];
rz(-0.87934914) q[2];
sx q[2];
rz(-2.9525625) q[2];
rz(0.14686251) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(1.7485025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015633164) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(2.2145859) q[0];
rz(-1.758763) q[1];
sx q[1];
rz(-2.5320876) q[1];
sx q[1];
rz(-1.6519201) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531567) q[0];
sx q[0];
rz(-0.90421593) q[0];
sx q[0];
rz(-0.96320926) q[0];
rz(-pi) q[1];
rz(2.0428558) q[2];
sx q[2];
rz(-1.5257116) q[2];
sx q[2];
rz(2.2301607) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5786963) q[1];
sx q[1];
rz(-1.7637196) q[1];
sx q[1];
rz(1.3955411) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49976607) q[3];
sx q[3];
rz(-2.0511813) q[3];
sx q[3];
rz(-1.2479316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5902517) q[2];
sx q[2];
rz(-0.44289032) q[2];
sx q[2];
rz(1.7112973) q[2];
rz(-0.57724214) q[3];
sx q[3];
rz(-0.8876628) q[3];
sx q[3];
rz(1.3841217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.5532613) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(0.28840315) q[0];
rz(-0.53238955) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(-0.14702252) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4924016) q[0];
sx q[0];
rz(-2.1426139) q[0];
sx q[0];
rz(-2.6834821) q[0];
x q[1];
rz(0.23994259) q[2];
sx q[2];
rz(-1.6389756) q[2];
sx q[2];
rz(-2.6477674) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3210996) q[1];
sx q[1];
rz(-1.5309146) q[1];
sx q[1];
rz(-0.95221968) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9045194) q[3];
sx q[3];
rz(-2.5519538) q[3];
sx q[3];
rz(-1.3414563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0816575) q[2];
sx q[2];
rz(-0.43490484) q[2];
sx q[2];
rz(2.3975513) q[2];
rz(2.384281) q[3];
sx q[3];
rz(-1.363874) q[3];
sx q[3];
rz(-1.3967167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
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
rz(0.4521162) q[2];
sx q[2];
rz(-1.631626) q[2];
sx q[2];
rz(0.22125868) q[2];
rz(-1.0088624) q[3];
sx q[3];
rz(-1.4603793) q[3];
sx q[3];
rz(0.59752656) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
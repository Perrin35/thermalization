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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1544224) q[0];
sx q[0];
rz(-1.3599456) q[0];
sx q[0];
rz(-0.2984557) q[0];
x q[1];
rz(-2.9169693) q[2];
sx q[2];
rz(-2.7135239) q[2];
sx q[2];
rz(3.012804) q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
x q[2];
rz(1.647244) q[3];
sx q[3];
rz(-0.41562286) q[3];
sx q[3];
rz(-1.7474183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5518387) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(0.067967728) q[2];
rz(0.12456482) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92000604) q[0];
sx q[0];
rz(-0.13555549) q[0];
sx q[0];
rz(0.24366972) q[0];
rz(2.5098353) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(-1.3557281) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6512017) q[0];
sx q[0];
rz(-0.25679195) q[0];
sx q[0];
rz(-1.4635565) q[0];
rz(-pi) q[1];
rz(2.3089356) q[2];
sx q[2];
rz(-2.6359733) q[2];
sx q[2];
rz(2.6308699) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7999468) q[1];
sx q[1];
rz(-0.46657545) q[1];
sx q[1];
rz(1.8599618) q[1];
rz(0.069085932) q[3];
sx q[3];
rz(-0.59327263) q[3];
sx q[3];
rz(-3.1309576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0624861) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(-2.8919354) q[2];
rz(0.50659531) q[3];
sx q[3];
rz(-1.6258312) q[3];
sx q[3];
rz(-2.8095968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(0.24519414) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(-0.89843345) q[0];
rz(1.3348745) q[1];
sx q[1];
rz(-1.9060262) q[1];
sx q[1];
rz(1.2737087) q[1];
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
x q[1];
rz(-1.4726228) q[2];
sx q[2];
rz(-1.2766826) q[2];
sx q[2];
rz(1.836118) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85325235) q[1];
sx q[1];
rz(-0.77008343) q[1];
sx q[1];
rz(-0.5132765) q[1];
rz(0.55176118) q[3];
sx q[3];
rz(-0.80358395) q[3];
sx q[3];
rz(-1.7337652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0818103) q[2];
sx q[2];
rz(-2.5368097) q[2];
sx q[2];
rz(-2.188142) q[2];
rz(-0.034514286) q[3];
sx q[3];
rz(-0.78648609) q[3];
sx q[3];
rz(2.9147193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8614486) q[0];
sx q[0];
rz(-0.21629688) q[0];
sx q[0];
rz(0.24818534) q[0];
rz(2.10363) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(-3.0674556) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5384597) q[0];
sx q[0];
rz(-0.54766253) q[0];
sx q[0];
rz(-1.0663701) q[0];
rz(1.147869) q[2];
sx q[2];
rz(-0.72407702) q[2];
sx q[2];
rz(0.99087447) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46596913) q[1];
sx q[1];
rz(-2.8040261) q[1];
sx q[1];
rz(1.2077431) q[1];
rz(-pi) q[2];
rz(-1.3965963) q[3];
sx q[3];
rz(-1.3341691) q[3];
sx q[3];
rz(-3.0803806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8903824) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(-3.1029491) q[2];
rz(0.97366992) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(-2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(-0.81659395) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(1.978925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8487932) q[0];
sx q[0];
rz(-1.6822364) q[0];
sx q[0];
rz(3.1116027) q[0];
x q[1];
rz(0.16935279) q[2];
sx q[2];
rz(-1.8893818) q[2];
sx q[2];
rz(2.7404075) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.77046493) q[1];
sx q[1];
rz(-0.41509291) q[1];
sx q[1];
rz(1.0789372) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88697042) q[3];
sx q[3];
rz(-0.63347048) q[3];
sx q[3];
rz(0.57852832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6267307) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(-0.14222063) q[2];
rz(2.2375315) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(-2.952125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0267923) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(-1.4676771) q[0];
rz(0.57178512) q[1];
sx q[1];
rz(-2.7829058) q[1];
sx q[1];
rz(-0.30803672) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.98826) q[0];
sx q[0];
rz(-1.4446265) q[0];
sx q[0];
rz(1.4754962) q[0];
x q[1];
rz(-2.3601989) q[2];
sx q[2];
rz(-2.3618744) q[2];
sx q[2];
rz(-1.3741072) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0526272) q[1];
sx q[1];
rz(-1.374377) q[1];
sx q[1];
rz(2.2599225) q[1];
rz(-pi) q[2];
rz(0.34835784) q[3];
sx q[3];
rz(-1.4847241) q[3];
sx q[3];
rz(-1.100988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3623111) q[2];
sx q[2];
rz(-1.4004536) q[2];
sx q[2];
rz(-1.1479088) q[2];
rz(-0.71427304) q[3];
sx q[3];
rz(-0.89759421) q[3];
sx q[3];
rz(2.4961297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4797453) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(0.1299783) q[0];
rz(-0.030844363) q[1];
sx q[1];
rz(-1.2896616) q[1];
sx q[1];
rz(0.67108363) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86258793) q[0];
sx q[0];
rz(-0.064084856) q[0];
sx q[0];
rz(-2.5628753) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7550049) q[2];
sx q[2];
rz(-1.6008018) q[2];
sx q[2];
rz(1.9629994) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56747251) q[1];
sx q[1];
rz(-1.4680871) q[1];
sx q[1];
rz(-1.9436388) q[1];
rz(-pi) q[2];
rz(0.047366553) q[3];
sx q[3];
rz(-0.8960552) q[3];
sx q[3];
rz(0.55344068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0188296) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(0.57787952) q[2];
rz(0.028586483) q[3];
sx q[3];
rz(-1.8609906) q[3];
sx q[3];
rz(1.8813429) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9776483) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(0.41241616) q[0];
rz(1.4498129) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(-1.9746045) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0359356) q[0];
sx q[0];
rz(-1.783839) q[0];
sx q[0];
rz(-0.94353326) q[0];
x q[1];
rz(1.251986) q[2];
sx q[2];
rz(-0.73482162) q[2];
sx q[2];
rz(0.97359818) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.941277) q[1];
sx q[1];
rz(-1.3822767) q[1];
sx q[1];
rz(-1.6016866) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9917631) q[3];
sx q[3];
rz(-1.0458046) q[3];
sx q[3];
rz(0.52469745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.940544) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(-2.9525625) q[2];
rz(2.9947301) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(1.3930901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015633164) q[0];
sx q[0];
rz(-1.3110315) q[0];
sx q[0];
rz(2.2145859) q[0];
rz(1.758763) q[1];
sx q[1];
rz(-0.60950509) q[1];
sx q[1];
rz(1.4896726) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4091858) q[0];
sx q[0];
rz(-0.86940765) q[0];
sx q[0];
rz(2.5138445) q[0];
rz(-pi) q[1];
rz(-2.0428558) q[2];
sx q[2];
rz(-1.5257116) q[2];
sx q[2];
rz(-2.2301607) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5786963) q[1];
sx q[1];
rz(-1.7637196) q[1];
sx q[1];
rz(-1.3955411) q[1];
x q[2];
rz(2.3143523) q[3];
sx q[3];
rz(-0.67875553) q[3];
sx q[3];
rz(2.7620706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5513409) q[2];
sx q[2];
rz(-2.6987023) q[2];
sx q[2];
rz(1.7112973) q[2];
rz(2.5643505) q[3];
sx q[3];
rz(-2.2539299) q[3];
sx q[3];
rz(-1.3841217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.5532613) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(2.8531895) q[0];
rz(-2.6092031) q[1];
sx q[1];
rz(-2.6817697) q[1];
sx q[1];
rz(-0.14702252) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4924016) q[0];
sx q[0];
rz(-0.99897879) q[0];
sx q[0];
rz(-2.6834821) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8617919) q[2];
sx q[2];
rz(-2.8923312) q[2];
sx q[2];
rz(-1.793043) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8359369) q[1];
sx q[1];
rz(-0.61969295) q[1];
sx q[1];
rz(1.5020919) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5649928) q[3];
sx q[3];
rz(-1.4398265) q[3];
sx q[3];
rz(-0.031158202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0816575) q[2];
sx q[2];
rz(-2.7066878) q[2];
sx q[2];
rz(0.74404136) q[2];
rz(-2.384281) q[3];
sx q[3];
rz(-1.7777187) q[3];
sx q[3];
rz(-1.3967167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.025678) q[0];
sx q[0];
rz(-2.0712576) q[0];
sx q[0];
rz(2.0448137) q[0];
rz(0.81746447) q[1];
sx q[1];
rz(-1.9348963) q[1];
sx q[1];
rz(2.5111326) q[1];
rz(3.0030737) q[2];
sx q[2];
rz(-0.45590966) q[2];
sx q[2];
rz(1.6675303) q[2];
rz(-3.01132) q[3];
sx q[3];
rz(-2.1289005) q[3];
sx q[3];
rz(2.0990513) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

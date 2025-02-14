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
rz(-1.350116) q[0];
sx q[0];
rz(3.8173563) q[0];
sx q[0];
rz(9.4326333) q[0];
rz(0.76771843) q[1];
sx q[1];
rz(1.6176728) q[1];
sx q[1];
rz(9.6674506) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3033256) q[0];
sx q[0];
rz(-1.6969862) q[0];
sx q[0];
rz(0.85889205) q[0];
x q[1];
rz(-1.2062293) q[2];
sx q[2];
rz(-1.1375712) q[2];
sx q[2];
rz(1.3490167) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3585904) q[1];
sx q[1];
rz(-1.5944408) q[1];
sx q[1];
rz(-2.7882697) q[1];
rz(-1.3245565) q[3];
sx q[3];
rz(-2.221635) q[3];
sx q[3];
rz(-0.090373978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1050038) q[2];
sx q[2];
rz(-1.236311) q[2];
sx q[2];
rz(-0.71199065) q[2];
rz(2.8365734) q[3];
sx q[3];
rz(-0.22235338) q[3];
sx q[3];
rz(3.0960848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63118339) q[0];
sx q[0];
rz(-1.864186) q[0];
sx q[0];
rz(-1.3674059) q[0];
rz(-0.039904682) q[1];
sx q[1];
rz(-2.3343562) q[1];
sx q[1];
rz(-0.59919277) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0586226) q[0];
sx q[0];
rz(-2.6520067) q[0];
sx q[0];
rz(0.24607308) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56053253) q[2];
sx q[2];
rz(-2.2403702) q[2];
sx q[2];
rz(-2.5472484) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89394648) q[1];
sx q[1];
rz(-2.1987913) q[1];
sx q[1];
rz(1.1537854) q[1];
rz(-2.4559095) q[3];
sx q[3];
rz(-1.7910379) q[3];
sx q[3];
rz(1.9626647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0590608) q[2];
sx q[2];
rz(-2.579687) q[2];
sx q[2];
rz(-1.3085636) q[2];
rz(3.1176873) q[3];
sx q[3];
rz(-1.5289565) q[3];
sx q[3];
rz(-1.5628975) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4557274) q[0];
sx q[0];
rz(-2.4754334) q[0];
sx q[0];
rz(-2.300793) q[0];
rz(2.1127545) q[1];
sx q[1];
rz(-2.7327171) q[1];
sx q[1];
rz(-1.6811446) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2551291) q[0];
sx q[0];
rz(-0.48454075) q[0];
sx q[0];
rz(2.8302829) q[0];
rz(-pi) q[1];
rz(0.65030789) q[2];
sx q[2];
rz(-0.57194369) q[2];
sx q[2];
rz(-0.24947333) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9509995) q[1];
sx q[1];
rz(-1.9036607) q[1];
sx q[1];
rz(-2.4458363) q[1];
rz(0.60734235) q[3];
sx q[3];
rz(-1.6173956) q[3];
sx q[3];
rz(2.6452015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8847522) q[2];
sx q[2];
rz(-1.931087) q[2];
sx q[2];
rz(2.9849226) q[2];
rz(-1.6414075) q[3];
sx q[3];
rz(-1.898396) q[3];
sx q[3];
rz(2.959804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0665862) q[0];
sx q[0];
rz(-1.2970507) q[0];
sx q[0];
rz(-3.060925) q[0];
rz(-1.1219885) q[1];
sx q[1];
rz(-2.5009218) q[1];
sx q[1];
rz(1.5544308) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4951404) q[0];
sx q[0];
rz(-0.82927726) q[0];
sx q[0];
rz(0.089228169) q[0];
rz(-pi) q[1];
rz(2.8055259) q[2];
sx q[2];
rz(-3.0185351) q[2];
sx q[2];
rz(0.50257909) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.69103855) q[1];
sx q[1];
rz(-0.86133707) q[1];
sx q[1];
rz(0.67994173) q[1];
rz(-pi) q[2];
rz(-0.46458475) q[3];
sx q[3];
rz(-2.4813093) q[3];
sx q[3];
rz(2.6990776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2598205) q[2];
sx q[2];
rz(-1.0910923) q[2];
sx q[2];
rz(2.1176977) q[2];
rz(-2.3648868) q[3];
sx q[3];
rz(-1.9909765) q[3];
sx q[3];
rz(2.1385433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7316932) q[0];
sx q[0];
rz(-2.2875146) q[0];
sx q[0];
rz(-0.62752974) q[0];
rz(2.4420786) q[1];
sx q[1];
rz(-1.0001837) q[1];
sx q[1];
rz(-2.2342009) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8555657) q[0];
sx q[0];
rz(-0.81840179) q[0];
sx q[0];
rz(1.6256385) q[0];
rz(1.4110231) q[2];
sx q[2];
rz(-1.7587307) q[2];
sx q[2];
rz(2.2794276) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4470565) q[1];
sx q[1];
rz(-1.9354104) q[1];
sx q[1];
rz(-1.0178119) q[1];
rz(-1.6697407) q[3];
sx q[3];
rz(-0.53831783) q[3];
sx q[3];
rz(2.7145999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13713914) q[2];
sx q[2];
rz(-1.8567825) q[2];
sx q[2];
rz(-0.84490204) q[2];
rz(-2.4431303) q[3];
sx q[3];
rz(-0.74028492) q[3];
sx q[3];
rz(-0.79160488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.38248211) q[0];
sx q[0];
rz(-1.6805205) q[0];
sx q[0];
rz(1.8735029) q[0];
rz(-1.6146487) q[1];
sx q[1];
rz(-2.0497649) q[1];
sx q[1];
rz(2.7472034) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1067692) q[0];
sx q[0];
rz(-0.64578694) q[0];
sx q[0];
rz(1.4950947) q[0];
rz(-0.092124513) q[2];
sx q[2];
rz(-0.73706223) q[2];
sx q[2];
rz(-1.9305522) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.907669) q[1];
sx q[1];
rz(-2.4872997) q[1];
sx q[1];
rz(3.0680502) q[1];
x q[2];
rz(2.9378618) q[3];
sx q[3];
rz(-0.36760783) q[3];
sx q[3];
rz(0.59083168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4322728) q[2];
sx q[2];
rz(-0.061881438) q[2];
sx q[2];
rz(-0.56149948) q[2];
rz(2.0105441) q[3];
sx q[3];
rz(-0.80315042) q[3];
sx q[3];
rz(-2.6530182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5442218) q[0];
sx q[0];
rz(-2.9591296) q[0];
sx q[0];
rz(2.9083948) q[0];
rz(-0.11416642) q[1];
sx q[1];
rz(-2.3515067) q[1];
sx q[1];
rz(2.9679969) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4712119) q[0];
sx q[0];
rz(-2.2224373) q[0];
sx q[0];
rz(0.46335133) q[0];
rz(-pi) q[1];
rz(0.38192449) q[2];
sx q[2];
rz(-1.1537103) q[2];
sx q[2];
rz(0.34427364) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.50362316) q[1];
sx q[1];
rz(-1.733128) q[1];
sx q[1];
rz(1.3470525) q[1];
x q[2];
rz(-0.18867698) q[3];
sx q[3];
rz(-1.5828307) q[3];
sx q[3];
rz(0.4922315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0357828) q[2];
sx q[2];
rz(-0.95869392) q[2];
sx q[2];
rz(-2.2868273) q[2];
rz(-0.0828951) q[3];
sx q[3];
rz(-1.9708743) q[3];
sx q[3];
rz(3.0875201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-1.4891124) q[0];
sx q[0];
rz(-1.880045) q[0];
sx q[0];
rz(-1.1789119) q[0];
rz(1.0014125) q[1];
sx q[1];
rz(-1.8779571) q[1];
sx q[1];
rz(-2.9875535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4074041) q[0];
sx q[0];
rz(-1.896669) q[0];
sx q[0];
rz(2.4209321) q[0];
rz(-pi) q[1];
x q[1];
rz(2.270584) q[2];
sx q[2];
rz(-0.99620512) q[2];
sx q[2];
rz(0.99829295) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6012965) q[1];
sx q[1];
rz(-2.5754693) q[1];
sx q[1];
rz(-2.152312) q[1];
x q[2];
rz(0.73518153) q[3];
sx q[3];
rz(-0.38215853) q[3];
sx q[3];
rz(-0.89565403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.68086326) q[2];
sx q[2];
rz(-2.0622084) q[2];
sx q[2];
rz(-0.66780773) q[2];
rz(1.7364511) q[3];
sx q[3];
rz(-1.7738155) q[3];
sx q[3];
rz(-2.5051129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77968303) q[0];
sx q[0];
rz(-1.6661665) q[0];
sx q[0];
rz(-0.59378004) q[0];
rz(1.8608015) q[1];
sx q[1];
rz(-2.180876) q[1];
sx q[1];
rz(-1.8437754) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7379388) q[0];
sx q[0];
rz(-2.5526921) q[0];
sx q[0];
rz(-3.011601) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62784451) q[2];
sx q[2];
rz(-1.3165917) q[2];
sx q[2];
rz(-1.1140299) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5169277) q[1];
sx q[1];
rz(-1.4447277) q[1];
sx q[1];
rz(0.36051118) q[1];
x q[2];
rz(0.16420096) q[3];
sx q[3];
rz(-2.3275725) q[3];
sx q[3];
rz(-2.9661953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7456776) q[2];
sx q[2];
rz(-3.119645) q[2];
sx q[2];
rz(1.6321261) q[2];
rz(-3.0294026) q[3];
sx q[3];
rz(-2.1144919) q[3];
sx q[3];
rz(-1.3794911) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47141948) q[0];
sx q[0];
rz(-1.0059953) q[0];
sx q[0];
rz(-0.26563409) q[0];
rz(0.98948014) q[1];
sx q[1];
rz(-1.2629291) q[1];
sx q[1];
rz(-0.33214733) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0672537) q[0];
sx q[0];
rz(-1.566972) q[0];
sx q[0];
rz(1.2373562) q[0];
rz(-pi) q[1];
rz(-0.99224706) q[2];
sx q[2];
rz(-1.3253115) q[2];
sx q[2];
rz(-3.0089889) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2492318) q[1];
sx q[1];
rz(-2.7977825) q[1];
sx q[1];
rz(-1.4268141) q[1];
x q[2];
rz(1.8197377) q[3];
sx q[3];
rz(-0.1934847) q[3];
sx q[3];
rz(1.8850808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0418479) q[2];
sx q[2];
rz(-1.0724649) q[2];
sx q[2];
rz(0.15038807) q[2];
rz(2.2930875) q[3];
sx q[3];
rz(-0.34981194) q[3];
sx q[3];
rz(1.8457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083241845) q[0];
sx q[0];
rz(-1.8235089) q[0];
sx q[0];
rz(1.9705809) q[0];
rz(-1.8104443) q[1];
sx q[1];
rz(-0.7970627) q[1];
sx q[1];
rz(2.0373559) q[1];
rz(-1.4990357) q[2];
sx q[2];
rz(-2.4007779) q[2];
sx q[2];
rz(3.130198) q[2];
rz(-0.70677291) q[3];
sx q[3];
rz(-1.5903683) q[3];
sx q[3];
rz(-0.49177468) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

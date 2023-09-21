OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(3.9324023) q[0];
sx q[0];
rz(12.232236) q[0];
rz(2.6842527) q[1];
sx q[1];
rz(-2.1973124) q[1];
sx q[1];
rz(1.9231208) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89361184) q[0];
sx q[0];
rz(-0.75151822) q[0];
sx q[0];
rz(3.0889838) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46484868) q[2];
sx q[2];
rz(-1.5343622) q[2];
sx q[2];
rz(-0.51793098) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7114746) q[1];
sx q[1];
rz(-0.700044) q[1];
sx q[1];
rz(2.3858566) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2741332) q[3];
sx q[3];
rz(-0.28775035) q[3];
sx q[3];
rz(-0.61144637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5228287) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(2.5640326) q[2];
rz(-1.9918359) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(-0.66453385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1537271) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(-0.38744774) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(1.739025) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1419066) q[0];
sx q[0];
rz(-1.5748595) q[0];
sx q[0];
rz(3.1121029) q[0];
x q[1];
rz(-2.5820929) q[2];
sx q[2];
rz(-1.3845978) q[2];
sx q[2];
rz(2.6452243) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.034033876) q[1];
sx q[1];
rz(-2.4817913) q[1];
sx q[1];
rz(0.69114139) q[1];
rz(-3.0121147) q[3];
sx q[3];
rz(-2.5054512) q[3];
sx q[3];
rz(2.5667218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7188321) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(0.31769162) q[2];
rz(0.20673949) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(2.3247705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7725672) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(1.4136219) q[0];
rz(2.6638022) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(-0.40107045) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6221878) q[0];
sx q[0];
rz(-1.4675958) q[0];
sx q[0];
rz(-1.2445356) q[0];
x q[1];
rz(1.3813854) q[2];
sx q[2];
rz(-2.2064798) q[2];
sx q[2];
rz(2.1798101) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2420826) q[1];
sx q[1];
rz(-0.84065719) q[1];
sx q[1];
rz(-1.5793369) q[1];
x q[2];
rz(0.19604711) q[3];
sx q[3];
rz(-1.9571597) q[3];
sx q[3];
rz(2.8783609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(2.5857914) q[2];
rz(2.1650971) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34898409) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(-1.697631) q[0];
rz(1.5199039) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(-2.8881883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4454173) q[0];
sx q[0];
rz(-1.6735958) q[0];
sx q[0];
rz(2.1131383) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6035945) q[2];
sx q[2];
rz(-2.6555853) q[2];
sx q[2];
rz(-0.82644586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6654012) q[1];
sx q[1];
rz(-2.5051077) q[1];
sx q[1];
rz(-2.979216) q[1];
x q[2];
rz(-3.0074189) q[3];
sx q[3];
rz(-2.4878256) q[3];
sx q[3];
rz(-1.3282446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.110934) q[2];
sx q[2];
rz(-1.7548283) q[2];
sx q[2];
rz(-2.3542662) q[2];
rz(-2.2287255) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(2.1319938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.426429) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(0.062967904) q[0];
rz(3.0175623) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(0.45809349) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3315017) q[0];
sx q[0];
rz(-2.0134263) q[0];
sx q[0];
rz(1.5733733) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1797921) q[2];
sx q[2];
rz(-0.70448175) q[2];
sx q[2];
rz(1.5693762) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8372247) q[1];
sx q[1];
rz(-2.3886884) q[1];
sx q[1];
rz(1.3992845) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0466988) q[3];
sx q[3];
rz(-1.0720383) q[3];
sx q[3];
rz(-0.14548485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(-3.138792) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(-1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(2.2221185) q[0];
rz(3.0793076) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(-1.8744291) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29612449) q[0];
sx q[0];
rz(-2.2757029) q[0];
sx q[0];
rz(-2.092093) q[0];
rz(1.7153347) q[2];
sx q[2];
rz(-2.309531) q[2];
sx q[2];
rz(1.2046255) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8333203) q[1];
sx q[1];
rz(-0.38494021) q[1];
sx q[1];
rz(-2.1182548) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0165765) q[3];
sx q[3];
rz(-0.52281724) q[3];
sx q[3];
rz(1.0558053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.9617812) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(-2.5578257) q[2];
rz(-2.4328655) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(-0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07847438) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(-1.0725347) q[0];
rz(0.5468927) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(2.0297208) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5522963) q[0];
sx q[0];
rz(-0.49320212) q[0];
sx q[0];
rz(1.769355) q[0];
rz(-pi) q[1];
x q[1];
rz(2.578031) q[2];
sx q[2];
rz(-2.5346018) q[2];
sx q[2];
rz(-2.3129472) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.062473) q[1];
sx q[1];
rz(-1.8603431) q[1];
sx q[1];
rz(1.6505961) q[1];
rz(-pi) q[2];
rz(-2.7518919) q[3];
sx q[3];
rz(-1.1003564) q[3];
sx q[3];
rz(2.5999992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.303858) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(-1.1676577) q[2];
rz(-1.5363103) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(-1.4060098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(-0.73079601) q[0];
rz(-0.90019512) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(-0.75497595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8687815) q[0];
sx q[0];
rz(-0.73505721) q[0];
sx q[0];
rz(1.3520157) q[0];
rz(-pi) q[1];
rz(2.412699) q[2];
sx q[2];
rz(-1.1155827) q[2];
sx q[2];
rz(-2.560937) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0874487) q[1];
sx q[1];
rz(-2.1463697) q[1];
sx q[1];
rz(1.4364442) q[1];
rz(-0.49236492) q[3];
sx q[3];
rz(-1.0400606) q[3];
sx q[3];
rz(-2.9363971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3770611) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(2.5047452) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(-1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.72717845) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(1.138858) q[0];
rz(-0.75421929) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(0.019502217) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7050539) q[0];
sx q[0];
rz(-2.9199613) q[0];
sx q[0];
rz(-1.9428695) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92658652) q[2];
sx q[2];
rz(-0.24878657) q[2];
sx q[2];
rz(-0.34005806) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6448977) q[1];
sx q[1];
rz(-2.1870107) q[1];
sx q[1];
rz(0.22209054) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0476258) q[3];
sx q[3];
rz(-2.7928824) q[3];
sx q[3];
rz(-2.204012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1372244) q[2];
sx q[2];
rz(-1.7251816) q[2];
sx q[2];
rz(-0.84890378) q[2];
rz(0.38765872) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(-1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3432817) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(2.6570901) q[0];
rz(1.3867406) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(1.1482931) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1185547) q[0];
sx q[0];
rz(-1.7114637) q[0];
sx q[0];
rz(-1.5492357) q[0];
rz(-pi) q[1];
rz(0.63126385) q[2];
sx q[2];
rz(-1.751465) q[2];
sx q[2];
rz(-2.0618912) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.53459586) q[1];
sx q[1];
rz(-1.3998704) q[1];
sx q[1];
rz(2.2439438) q[1];
rz(-pi) q[2];
rz(-0.34928068) q[3];
sx q[3];
rz(-1.5216856) q[3];
sx q[3];
rz(-1.1790566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91223532) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(-0.36515507) q[2];
rz(0.12864104) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(-2.685759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.1098332) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(2.1784492) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(1.1744432) q[2];
sx q[2];
rz(-2.195993) q[2];
sx q[2];
rz(3.1396951) q[2];
rz(0.036988463) q[3];
sx q[3];
rz(-1.0109517) q[3];
sx q[3];
rz(3.0299822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-1.090156) q[0];
sx q[0];
rz(-0.12026726) q[0];
sx q[0];
rz(3.0647478) q[0];
rz(-0.85637561) q[1];
sx q[1];
rz(-2.0937347) q[1];
sx q[1];
rz(0.82077789) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73859204) q[0];
sx q[0];
rz(-1.101491) q[0];
sx q[0];
rz(1.655939) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9772867) q[2];
sx q[2];
rz(-2.2191357) q[2];
sx q[2];
rz(-0.68337464) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.32626881) q[1];
sx q[1];
rz(-2.2827049) q[1];
sx q[1];
rz(0.91628847) q[1];
x q[2];
rz(-1.7480811) q[3];
sx q[3];
rz(-1.0330135) q[3];
sx q[3];
rz(0.017740109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1138175) q[2];
sx q[2];
rz(-2.5706988) q[2];
sx q[2];
rz(-1.5296193) q[2];
rz(-2.2714603) q[3];
sx q[3];
rz(-1.8703929) q[3];
sx q[3];
rz(0.99942708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.40365264) q[0];
sx q[0];
rz(-1.0920748) q[0];
sx q[0];
rz(-1.8988443) q[0];
rz(-2.2620762) q[1];
sx q[1];
rz(-2.1915235) q[1];
sx q[1];
rz(-1.3290728) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14318109) q[0];
sx q[0];
rz(-2.3685026) q[0];
sx q[0];
rz(-2.4804658) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6471001) q[2];
sx q[2];
rz(-1.1323511) q[2];
sx q[2];
rz(-0.019767337) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76960388) q[1];
sx q[1];
rz(-2.0298483) q[1];
sx q[1];
rz(0.52476669) q[1];
x q[2];
rz(-1.7326844) q[3];
sx q[3];
rz(-0.88380764) q[3];
sx q[3];
rz(-1.5702198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8326412) q[2];
sx q[2];
rz(-1.2244886) q[2];
sx q[2];
rz(1.857117) q[2];
rz(-0.66713157) q[3];
sx q[3];
rz(-2.7693373) q[3];
sx q[3];
rz(-0.66481203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.852378) q[0];
sx q[0];
rz(-2.7279655) q[0];
sx q[0];
rz(-1.1997696) q[0];
rz(-1.86357) q[1];
sx q[1];
rz(-2.8646902) q[1];
sx q[1];
rz(2.3866381) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65855712) q[0];
sx q[0];
rz(-0.45679507) q[0];
sx q[0];
rz(1.5617338) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4226025) q[2];
sx q[2];
rz(-0.8575927) q[2];
sx q[2];
rz(2.7813965) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36981407) q[1];
sx q[1];
rz(-1.6057456) q[1];
sx q[1];
rz(-0.82280226) q[1];
rz(-1.0126013) q[3];
sx q[3];
rz(-2.1821488) q[3];
sx q[3];
rz(3.0773036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7167012) q[2];
sx q[2];
rz(-2.5548013) q[2];
sx q[2];
rz(1.971604) q[2];
rz(-0.26015002) q[3];
sx q[3];
rz(-1.6157506) q[3];
sx q[3];
rz(2.9366176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78885704) q[0];
sx q[0];
rz(-1.641888) q[0];
sx q[0];
rz(2.4810171) q[0];
rz(-3.0109516) q[1];
sx q[1];
rz(-2.525593) q[1];
sx q[1];
rz(0.43636838) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35231464) q[0];
sx q[0];
rz(-1.3320062) q[0];
sx q[0];
rz(1.055607) q[0];
rz(-pi) q[1];
rz(-0.58664586) q[2];
sx q[2];
rz(-1.9244951) q[2];
sx q[2];
rz(-2.2207137) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1088686) q[1];
sx q[1];
rz(-2.0744262) q[1];
sx q[1];
rz(-0.41028604) q[1];
x q[2];
rz(-0.34109165) q[3];
sx q[3];
rz(-1.4779096) q[3];
sx q[3];
rz(-0.22985458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.067430647) q[2];
sx q[2];
rz(-0.51965153) q[2];
sx q[2];
rz(2.0036428) q[2];
rz(2.2940995) q[3];
sx q[3];
rz(-1.3769826) q[3];
sx q[3];
rz(-1.725089) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6013019) q[0];
sx q[0];
rz(-1.8922946) q[0];
sx q[0];
rz(2.9345595) q[0];
rz(-1.3218468) q[1];
sx q[1];
rz(-1.9925947) q[1];
sx q[1];
rz(1.0386937) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90090758) q[0];
sx q[0];
rz(-2.0863914) q[0];
sx q[0];
rz(-2.2365966) q[0];
rz(1.7798575) q[2];
sx q[2];
rz(-2.1809177) q[2];
sx q[2];
rz(2.3336712) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0952044) q[1];
sx q[1];
rz(-2.5051053) q[1];
sx q[1];
rz(3.1134136) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8349325) q[3];
sx q[3];
rz(-1.6741317) q[3];
sx q[3];
rz(2.5984026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1553104) q[2];
sx q[2];
rz(-2.0955413) q[2];
sx q[2];
rz(-2.5347064) q[2];
rz(0.33489975) q[3];
sx q[3];
rz(-1.5891985) q[3];
sx q[3];
rz(1.5978395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63810054) q[0];
sx q[0];
rz(-1.7146716) q[0];
sx q[0];
rz(2.959429) q[0];
rz(0.69193524) q[1];
sx q[1];
rz(-1.4357932) q[1];
sx q[1];
rz(-1.2729794) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70304027) q[0];
sx q[0];
rz(-1.5219757) q[0];
sx q[0];
rz(0.010558616) q[0];
rz(-pi) q[1];
rz(0.23006205) q[2];
sx q[2];
rz(-1.7126691) q[2];
sx q[2];
rz(1.0070247) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5598076) q[1];
sx q[1];
rz(-2.323287) q[1];
sx q[1];
rz(2.7831828) q[1];
rz(-0.16830651) q[3];
sx q[3];
rz(-2.7790439) q[3];
sx q[3];
rz(-2.1358638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89340297) q[2];
sx q[2];
rz(-2.4290163) q[2];
sx q[2];
rz(2.2656061) q[2];
rz(-0.28410965) q[3];
sx q[3];
rz(-0.94485372) q[3];
sx q[3];
rz(-0.42377728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9559602) q[0];
sx q[0];
rz(-0.9469339) q[0];
sx q[0];
rz(-2.60485) q[0];
rz(0.56450379) q[1];
sx q[1];
rz(-0.91465488) q[1];
sx q[1];
rz(1.5980665) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48115981) q[0];
sx q[0];
rz(-2.8114722) q[0];
sx q[0];
rz(0.31417234) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41588628) q[2];
sx q[2];
rz(-1.3817817) q[2];
sx q[2];
rz(0.11014665) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.35245124) q[1];
sx q[1];
rz(-1.1256212) q[1];
sx q[1];
rz(0.51207249) q[1];
x q[2];
rz(-2.2143977) q[3];
sx q[3];
rz(-0.72228408) q[3];
sx q[3];
rz(2.3658906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87968612) q[2];
sx q[2];
rz(-1.3085111) q[2];
sx q[2];
rz(1.7730664) q[2];
rz(-0.3174828) q[3];
sx q[3];
rz(-1.2513132) q[3];
sx q[3];
rz(2.4449463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19318652) q[0];
sx q[0];
rz(-2.6115311) q[0];
sx q[0];
rz(1.0411221) q[0];
rz(0.55024534) q[1];
sx q[1];
rz(-0.3717652) q[1];
sx q[1];
rz(0.1951898) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6379162) q[0];
sx q[0];
rz(-1.6244672) q[0];
sx q[0];
rz(-0.1366833) q[0];
rz(-pi) q[1];
rz(-1.0894906) q[2];
sx q[2];
rz(-1.6688108) q[2];
sx q[2];
rz(-2.5070407) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.079551815) q[1];
sx q[1];
rz(-1.5925381) q[1];
sx q[1];
rz(-1.8575536) q[1];
rz(2.3161668) q[3];
sx q[3];
rz(-1.22681) q[3];
sx q[3];
rz(1.0937364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0027085) q[2];
sx q[2];
rz(-1.1382853) q[2];
sx q[2];
rz(1.8769334) q[2];
rz(1.5876611) q[3];
sx q[3];
rz(-1.2041661) q[3];
sx q[3];
rz(1.1519661) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78466648) q[0];
sx q[0];
rz(-0.21360989) q[0];
sx q[0];
rz(2.7333562) q[0];
rz(-1.4534072) q[1];
sx q[1];
rz(-1.6994349) q[1];
sx q[1];
rz(0.85809392) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5018117) q[0];
sx q[0];
rz(-2.1970891) q[0];
sx q[0];
rz(-2.7485952) q[0];
x q[1];
rz(1.6995605) q[2];
sx q[2];
rz(-1.4742645) q[2];
sx q[2];
rz(-1.0573204) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.063805997) q[1];
sx q[1];
rz(-1.1764472) q[1];
sx q[1];
rz(-2.0057682) q[1];
rz(1.7210049) q[3];
sx q[3];
rz(-0.69828639) q[3];
sx q[3];
rz(2.5774389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.475829) q[2];
sx q[2];
rz(-2.2593311) q[2];
sx q[2];
rz(2.8909454) q[2];
rz(-2.2541239) q[3];
sx q[3];
rz(-1.1424501) q[3];
sx q[3];
rz(2.7953776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733521) q[0];
sx q[0];
rz(-2.1526985) q[0];
sx q[0];
rz(-2.2762779) q[0];
rz(-2.6189651) q[1];
sx q[1];
rz(-0.80895439) q[1];
sx q[1];
rz(1.4563947) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53082836) q[0];
sx q[0];
rz(-1.3757154) q[0];
sx q[0];
rz(3.0733956) q[0];
x q[1];
rz(-1.3550884) q[2];
sx q[2];
rz(-2.459639) q[2];
sx q[2];
rz(-0.29302675) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2715312) q[1];
sx q[1];
rz(-1.2200331) q[1];
sx q[1];
rz(-2.6501772) q[1];
x q[2];
rz(-0.56444498) q[3];
sx q[3];
rz(-0.58957499) q[3];
sx q[3];
rz(-1.6605476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7437848) q[2];
sx q[2];
rz(-0.91181552) q[2];
sx q[2];
rz(-1.3249116) q[2];
rz(1.2815453) q[3];
sx q[3];
rz(-2.2478588) q[3];
sx q[3];
rz(0.76461422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.926173) q[0];
sx q[0];
rz(-1.8126491) q[0];
sx q[0];
rz(1.5480315) q[0];
rz(-0.46319766) q[1];
sx q[1];
rz(-1.5402272) q[1];
sx q[1];
rz(2.872749) q[1];
rz(2.7364667) q[2];
sx q[2];
rz(-2.4562757) q[2];
sx q[2];
rz(1.9669878) q[2];
rz(1.6921099) q[3];
sx q[3];
rz(-0.90150812) q[3];
sx q[3];
rz(1.5163803) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

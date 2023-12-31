OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.70513201) q[0];
sx q[0];
rz(6.8350514) q[0];
sx q[0];
rz(9.4466136) q[0];
rz(-0.39437374) q[1];
sx q[1];
rz(-1.6819277) q[1];
sx q[1];
rz(0.2149166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267933) q[0];
sx q[0];
rz(-2.7601295) q[0];
sx q[0];
rz(-0.8154072) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46918842) q[2];
sx q[2];
rz(-1.708963) q[2];
sx q[2];
rz(2.5772622) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5948407) q[1];
sx q[1];
rz(-1.1210124) q[1];
sx q[1];
rz(-2.2307598) q[1];
x q[2];
rz(-2.4258852) q[3];
sx q[3];
rz(-1.189609) q[3];
sx q[3];
rz(-2.3472799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.73137838) q[2];
sx q[2];
rz(-1.4593068) q[2];
sx q[2];
rz(0.56420502) q[2];
rz(-1.7764067) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0974225) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(0.92798293) q[0];
rz(1.1652975) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(-0.89675084) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3551536) q[0];
sx q[0];
rz(-0.73620287) q[0];
sx q[0];
rz(-1.8650706) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9727544) q[2];
sx q[2];
rz(-1.642792) q[2];
sx q[2];
rz(0.26693401) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.36205081) q[1];
sx q[1];
rz(-2.3565787) q[1];
sx q[1];
rz(-1.3951468) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0131301) q[3];
sx q[3];
rz(-2.3395174) q[3];
sx q[3];
rz(-0.80004063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8759878) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(-1.0401475) q[2];
rz(-1.6863719) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002157) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(-2.1133912) q[0];
rz(2.0630515) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(2.7064586) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4774287) q[0];
sx q[0];
rz(-1.6782883) q[0];
sx q[0];
rz(0.35350032) q[0];
x q[1];
rz(-0.52559678) q[2];
sx q[2];
rz(-0.28738775) q[2];
sx q[2];
rz(1.4488066) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9738234) q[1];
sx q[1];
rz(-0.5127738) q[1];
sx q[1];
rz(-1.9085401) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35878351) q[3];
sx q[3];
rz(-0.79625087) q[3];
sx q[3];
rz(0.96804726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6083287) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(0.30291525) q[2];
rz(-1.8164002) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(-0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-2.9451697) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(-0.91745013) q[0];
rz(0.67287412) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(2.8767169) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6661975) q[0];
sx q[0];
rz(-2.4817433) q[0];
sx q[0];
rz(3.092479) q[0];
rz(-pi) q[1];
rz(1.869309) q[2];
sx q[2];
rz(-1.8372756) q[2];
sx q[2];
rz(-0.68599115) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3600196) q[1];
sx q[1];
rz(-2.2847166) q[1];
sx q[1];
rz(0.54846958) q[1];
rz(-pi) q[2];
rz(-2.360965) q[3];
sx q[3];
rz(-0.98140162) q[3];
sx q[3];
rz(-1.2300223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.36310568) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(-1.5765566) q[2];
rz(-1.0270843) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(-1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(-0.47873163) q[0];
rz(-1.0331253) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(-0.95265257) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27272308) q[0];
sx q[0];
rz(-0.55991828) q[0];
sx q[0];
rz(0.22229226) q[0];
rz(2.7258337) q[2];
sx q[2];
rz(-1.2649049) q[2];
sx q[2];
rz(-1.3410459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.88541234) q[1];
sx q[1];
rz(-2.6773239) q[1];
sx q[1];
rz(1.0186362) q[1];
rz(-pi) q[2];
x q[2];
rz(0.080658241) q[3];
sx q[3];
rz(-0.47938743) q[3];
sx q[3];
rz(-0.56848923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4218563) q[2];
sx q[2];
rz(-2.77878) q[2];
sx q[2];
rz(-0.53058132) q[2];
rz(1.7355708) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(-2.3099242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.574061) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(1.5166327) q[0];
rz(-1.8364871) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(0.17257246) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0046376) q[0];
sx q[0];
rz(-1.1543659) q[0];
sx q[0];
rz(0.01491551) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0667604) q[2];
sx q[2];
rz(-2.6325588) q[2];
sx q[2];
rz(2.3376655) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3249358) q[1];
sx q[1];
rz(-1.6759733) q[1];
sx q[1];
rz(1.9865958) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8302912) q[3];
sx q[3];
rz(-0.65400306) q[3];
sx q[3];
rz(-1.9037387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3036348) q[2];
sx q[2];
rz(-1.6875608) q[2];
sx q[2];
rz(1.1266358) q[2];
rz(-0.78222328) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(1.8036028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28850266) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(2.4801168) q[0];
rz(-2.181197) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(2.3849934) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264383) q[0];
sx q[0];
rz(-1.4190136) q[0];
sx q[0];
rz(-3.0477235) q[0];
x q[1];
rz(-2.7301894) q[2];
sx q[2];
rz(-1.7559933) q[2];
sx q[2];
rz(1.7298557) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36180624) q[1];
sx q[1];
rz(-2.8621319) q[1];
sx q[1];
rz(-0.96868412) q[1];
rz(-1.6301304) q[3];
sx q[3];
rz(-1.3795128) q[3];
sx q[3];
rz(0.64185601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1371655) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(0.19443092) q[2];
rz(-0.91313177) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(0.066666691) q[0];
rz(-2.8170259) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(2.1527122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7522404) q[0];
sx q[0];
rz(-1.1712495) q[0];
sx q[0];
rz(0.35895343) q[0];
x q[1];
rz(3.0481911) q[2];
sx q[2];
rz(-1.5548692) q[2];
sx q[2];
rz(-1.8708558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1253742) q[1];
sx q[1];
rz(-1.6375293) q[1];
sx q[1];
rz(-1.1557505) q[1];
rz(-2.3905121) q[3];
sx q[3];
rz(-2.4820648) q[3];
sx q[3];
rz(0.86576033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4618335) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(-2.7339593) q[2];
rz(2.3729825) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(2.2176946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75893629) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(2.5323903) q[0];
rz(3.0464879) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(-0.87337714) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0926889) q[0];
sx q[0];
rz(-2.9498219) q[0];
sx q[0];
rz(1.9321241) q[0];
x q[1];
rz(0.14752578) q[2];
sx q[2];
rz(-1.4359056) q[2];
sx q[2];
rz(2.877176) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.020563515) q[1];
sx q[1];
rz(-2.4596446) q[1];
sx q[1];
rz(2.6418583) q[1];
x q[2];
rz(-0.94536762) q[3];
sx q[3];
rz(-2.4646467) q[3];
sx q[3];
rz(-0.35811801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1197027) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(2.4592887) q[2];
rz(-0.37426379) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69797126) q[0];
sx q[0];
rz(-1.3598096) q[0];
sx q[0];
rz(0.95296729) q[0];
rz(0.8264181) q[1];
sx q[1];
rz(-2.4024139) q[1];
sx q[1];
rz(1.7451161) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5293225) q[0];
sx q[0];
rz(-1.6565874) q[0];
sx q[0];
rz(0.32353185) q[0];
rz(-pi) q[1];
rz(-1.868532) q[2];
sx q[2];
rz(-2.2072788) q[2];
sx q[2];
rz(-2.0993078) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.1435946) q[1];
sx q[1];
rz(-1.8815787) q[1];
sx q[1];
rz(2.8423611) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26225984) q[3];
sx q[3];
rz(-1.5159303) q[3];
sx q[3];
rz(0.39615397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6754127) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(-0.15979016) q[2];
rz(2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(3.0970584) q[0];
sx q[0];
rz(-0.67561588) q[0];
sx q[0];
rz(-1.5560879) q[0];
rz(-0.13327577) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(0.940154) q[2];
sx q[2];
rz(-1.4174145) q[2];
sx q[2];
rz(-2.6819475) q[2];
rz(-3.0575183) q[3];
sx q[3];
rz(-2.6064059) q[3];
sx q[3];
rz(2.4168766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

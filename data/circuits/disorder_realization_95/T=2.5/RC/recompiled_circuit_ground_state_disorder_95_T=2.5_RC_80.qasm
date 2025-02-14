OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9606544) q[0];
sx q[0];
rz(-0.063194312) q[0];
sx q[0];
rz(2.5511281) q[0];
rz(-0.2904627) q[1];
sx q[1];
rz(1.7658748) q[1];
sx q[1];
rz(9.4212846) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9866512) q[0];
sx q[0];
rz(-1.9738324) q[0];
sx q[0];
rz(2.0986845) q[0];
x q[1];
rz(0.89031808) q[2];
sx q[2];
rz(-1.8125475) q[2];
sx q[2];
rz(1.7861451) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.063979186) q[1];
sx q[1];
rz(-0.031647041) q[1];
sx q[1];
rz(-0.74546234) q[1];
rz(-pi) q[2];
rz(-0.95432561) q[3];
sx q[3];
rz(-0.69479695) q[3];
sx q[3];
rz(-2.8692226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8840088) q[2];
sx q[2];
rz(-0.85870063) q[2];
sx q[2];
rz(-2.1250471) q[2];
rz(-0.86333418) q[3];
sx q[3];
rz(-0.038067929) q[3];
sx q[3];
rz(-1.4994924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6739552) q[0];
sx q[0];
rz(-1.4731982) q[0];
sx q[0];
rz(-0.13482811) q[0];
rz(0.22678953) q[1];
sx q[1];
rz(-3.0786381) q[1];
sx q[1];
rz(2.8917868) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1824719) q[0];
sx q[0];
rz(-1.9879436) q[0];
sx q[0];
rz(-2.4802101) q[0];
rz(3.089014) q[2];
sx q[2];
rz(-1.0710395) q[2];
sx q[2];
rz(1.6249958) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6625401) q[1];
sx q[1];
rz(-1.6483278) q[1];
sx q[1];
rz(1.5815602) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13447478) q[3];
sx q[3];
rz(-1.933483) q[3];
sx q[3];
rz(2.4464726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9164385) q[2];
sx q[2];
rz(-0.068000451) q[2];
sx q[2];
rz(-0.98006836) q[2];
rz(2.2465536) q[3];
sx q[3];
rz(-2.3990302) q[3];
sx q[3];
rz(-1.9612954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5880244) q[0];
sx q[0];
rz(-2.6346485) q[0];
sx q[0];
rz(-1.5356327) q[0];
rz(-1.6355248) q[1];
sx q[1];
rz(-2.3160544) q[1];
sx q[1];
rz(-2.1855386) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0900537) q[0];
sx q[0];
rz(-0.23583007) q[0];
sx q[0];
rz(0.43311849) q[0];
rz(-pi) q[1];
rz(-1.5044841) q[2];
sx q[2];
rz(-0.8408747) q[2];
sx q[2];
rz(0.17375565) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3967665) q[1];
sx q[1];
rz(-1.464572) q[1];
sx q[1];
rz(-2.3030119) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9185103) q[3];
sx q[3];
rz(-1.082431) q[3];
sx q[3];
rz(-0.85635105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8354336) q[2];
sx q[2];
rz(-0.71194887) q[2];
sx q[2];
rz(-2.5095059) q[2];
rz(2.2575374) q[3];
sx q[3];
rz(-0.017280936) q[3];
sx q[3];
rz(2.250905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.66568351) q[0];
sx q[0];
rz(-1.2983687) q[0];
sx q[0];
rz(0.16715288) q[0];
rz(-1.2627603) q[1];
sx q[1];
rz(-0.94307584) q[1];
sx q[1];
rz(1.7335588) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3638813) q[0];
sx q[0];
rz(-1.1486619) q[0];
sx q[0];
rz(-1.889983) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7825259) q[2];
sx q[2];
rz(-1.5671726) q[2];
sx q[2];
rz(0.45932367) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9973397) q[1];
sx q[1];
rz(-1.3648811) q[1];
sx q[1];
rz(-1.3763672) q[1];
rz(-pi) q[2];
rz(-2.3711331) q[3];
sx q[3];
rz(-1.3589199) q[3];
sx q[3];
rz(-0.4650863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.35272804) q[2];
sx q[2];
rz(-0.015268607) q[2];
sx q[2];
rz(2.7909279) q[2];
rz(-2.8777425) q[3];
sx q[3];
rz(-0.00084547384) q[3];
sx q[3];
rz(-1.2034169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.8769787) q[0];
sx q[0];
rz(-2.3961841) q[0];
sx q[0];
rz(-1.4234446) q[0];
rz(2.8599332) q[1];
sx q[1];
rz(-1.6335082) q[1];
sx q[1];
rz(1.8219832) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21319178) q[0];
sx q[0];
rz(-1.9921148) q[0];
sx q[0];
rz(-2.1405959) q[0];
rz(-pi) q[1];
rz(0.93792589) q[2];
sx q[2];
rz(-3.1073776) q[2];
sx q[2];
rz(0.93955112) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4180243) q[1];
sx q[1];
rz(-2.4279874) q[1];
sx q[1];
rz(-0.51939555) q[1];
rz(-pi) q[2];
rz(1.237147) q[3];
sx q[3];
rz(-1.9556442) q[3];
sx q[3];
rz(0.24144444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.264297) q[2];
sx q[2];
rz(-0.046961203) q[2];
sx q[2];
rz(-1.7128672) q[2];
rz(-1.9127539) q[3];
sx q[3];
rz(-0.052534025) q[3];
sx q[3];
rz(-0.088168941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0838715) q[0];
sx q[0];
rz(-0.11744048) q[0];
sx q[0];
rz(-0.55026662) q[0];
rz(0.30217198) q[1];
sx q[1];
rz(-1.3928394) q[1];
sx q[1];
rz(0.43513939) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9254375) q[0];
sx q[0];
rz(-1.3595681) q[0];
sx q[0];
rz(-1.5351377) q[0];
x q[1];
rz(-0.0074362292) q[2];
sx q[2];
rz(-1.5717634) q[2];
sx q[2];
rz(1.2449679) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8543431) q[1];
sx q[1];
rz(-1.6222427) q[1];
sx q[1];
rz(1.9958985) q[1];
rz(-0.76558785) q[3];
sx q[3];
rz(-2.2556955) q[3];
sx q[3];
rz(-1.9960038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90007323) q[2];
sx q[2];
rz(-3.1406431) q[2];
sx q[2];
rz(2.1692236) q[2];
rz(0.9995681) q[3];
sx q[3];
rz(-3.1344423) q[3];
sx q[3];
rz(1.0424559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5480176) q[0];
sx q[0];
rz(-1.2418208) q[0];
sx q[0];
rz(-1.0663363) q[0];
rz(-1.4284596) q[1];
sx q[1];
rz(-0.26930535) q[1];
sx q[1];
rz(1.8535293) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36670524) q[0];
sx q[0];
rz(-2.2572491) q[0];
sx q[0];
rz(-1.9839601) q[0];
rz(2.7884862) q[2];
sx q[2];
rz(-0.87279183) q[2];
sx q[2];
rz(1.6153796) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3235229) q[1];
sx q[1];
rz(-1.5574291) q[1];
sx q[1];
rz(-1.531015) q[1];
rz(-1.0405447) q[3];
sx q[3];
rz(-1.3589348) q[3];
sx q[3];
rz(-2.7303641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0224033) q[2];
sx q[2];
rz(-3.0722805) q[2];
sx q[2];
rz(-0.27528396) q[2];
rz(-2.9149808) q[3];
sx q[3];
rz(-0.35609069) q[3];
sx q[3];
rz(0.4376469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76448035) q[0];
sx q[0];
rz(-0.28700101) q[0];
sx q[0];
rz(0.58718938) q[0];
rz(1.5234692) q[1];
sx q[1];
rz(-2.0240929) q[1];
sx q[1];
rz(-1.5709741) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29059288) q[0];
sx q[0];
rz(-2.3208566) q[0];
sx q[0];
rz(0.59159578) q[0];
x q[1];
rz(-1.4373407) q[2];
sx q[2];
rz(-2.674006) q[2];
sx q[2];
rz(-3.140492) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.60698673) q[1];
sx q[1];
rz(-3.0879277) q[1];
sx q[1];
rz(2.1855658) q[1];
rz(0.32125485) q[3];
sx q[3];
rz(-2.1200088) q[3];
sx q[3];
rz(0.7859226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9470584) q[2];
sx q[2];
rz(-2.9583866) q[2];
sx q[2];
rz(1.8397231) q[2];
rz(-0.434508) q[3];
sx q[3];
rz(-0.040381581) q[3];
sx q[3];
rz(-0.44809189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3473564) q[0];
sx q[0];
rz(-3.0257822) q[0];
sx q[0];
rz(-1.7606803) q[0];
rz(1.5877089) q[1];
sx q[1];
rz(-2.1788308) q[1];
sx q[1];
rz(-0.10996058) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1292353) q[0];
sx q[0];
rz(-1.7596272) q[0];
sx q[0];
rz(-0.83602943) q[0];
rz(-0.25935632) q[2];
sx q[2];
rz(-2.4217941) q[2];
sx q[2];
rz(-2.6754745) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.090558) q[1];
sx q[1];
rz(-0.75737774) q[1];
sx q[1];
rz(1.5345598) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8300555) q[3];
sx q[3];
rz(-1.7230125) q[3];
sx q[3];
rz(-2.074034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.79003698) q[2];
sx q[2];
rz(-2.07708) q[2];
sx q[2];
rz(-2.7891187) q[2];
rz(2.500109) q[3];
sx q[3];
rz(-0.04647579) q[3];
sx q[3];
rz(-0.36267734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.27956692) q[0];
sx q[0];
rz(-0.27722219) q[0];
sx q[0];
rz(2.5773881) q[0];
rz(1.6212246) q[1];
sx q[1];
rz(-1.0313326) q[1];
sx q[1];
rz(3.0618111) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02139442) q[0];
sx q[0];
rz(-1.7225838) q[0];
sx q[0];
rz(0.89812507) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78861137) q[2];
sx q[2];
rz(-0.067316003) q[2];
sx q[2];
rz(2.2612417) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.141216) q[1];
sx q[1];
rz(-1.3151919) q[1];
sx q[1];
rz(-2.3325898) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8989927) q[3];
sx q[3];
rz(-1.4897457) q[3];
sx q[3];
rz(-1.9417861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0148049) q[2];
sx q[2];
rz(-0.0099651907) q[2];
sx q[2];
rz(2.0662181) q[2];
rz(2.3672095) q[3];
sx q[3];
rz(-0.024024809) q[3];
sx q[3];
rz(2.8703441) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8423691) q[0];
sx q[0];
rz(-1.5689701) q[0];
sx q[0];
rz(-1.5690621) q[0];
rz(-0.52539274) q[1];
sx q[1];
rz(-0.071594302) q[1];
sx q[1];
rz(-0.20803861) q[1];
rz(0.21389773) q[2];
sx q[2];
rz(-2.1281617) q[2];
sx q[2];
rz(-2.8585363) q[2];
rz(-2.0490859) q[3];
sx q[3];
rz(-0.94259613) q[3];
sx q[3];
rz(-2.5362956) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

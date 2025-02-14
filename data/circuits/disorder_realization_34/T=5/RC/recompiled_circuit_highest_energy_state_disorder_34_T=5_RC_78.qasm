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
rz(-1.7582769) q[0];
sx q[0];
rz(4.5780616) q[0];
sx q[0];
rz(8.4599001) q[0];
rz(2.3896253) q[1];
sx q[1];
rz(-2.7120092) q[1];
sx q[1];
rz(-2.092195) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0239068) q[0];
sx q[0];
rz(-1.2313594) q[0];
sx q[0];
rz(-2.1855559) q[0];
rz(-1.8951178) q[2];
sx q[2];
rz(-1.4359546) q[2];
sx q[2];
rz(-0.91712778) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.2856372) q[1];
sx q[1];
rz(-2.1826943) q[1];
sx q[1];
rz(0.4680856) q[1];
rz(-0.5792867) q[3];
sx q[3];
rz(-1.6484954) q[3];
sx q[3];
rz(0.9723817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6875978) q[2];
sx q[2];
rz(-1.6124085) q[2];
sx q[2];
rz(3.122984) q[2];
rz(-2.5971557) q[3];
sx q[3];
rz(-2.8114522) q[3];
sx q[3];
rz(2.4936567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7740771) q[0];
sx q[0];
rz(-1.6870966) q[0];
sx q[0];
rz(0.53502214) q[0];
rz(0.53994838) q[1];
sx q[1];
rz(-2.5060563) q[1];
sx q[1];
rz(1.3998869) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0947123) q[0];
sx q[0];
rz(-2.0167354) q[0];
sx q[0];
rz(-0.42973862) q[0];
rz(-pi) q[1];
rz(0.56324701) q[2];
sx q[2];
rz(-2.0958825) q[2];
sx q[2];
rz(-1.6890749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78838879) q[1];
sx q[1];
rz(-1.3223159) q[1];
sx q[1];
rz(-2.7491991) q[1];
x q[2];
rz(-0.61576188) q[3];
sx q[3];
rz(-1.0887556) q[3];
sx q[3];
rz(-0.69053135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0672368) q[2];
sx q[2];
rz(-1.3398193) q[2];
sx q[2];
rz(-2.2155217) q[2];
rz(2.1615084) q[3];
sx q[3];
rz(-0.96433774) q[3];
sx q[3];
rz(2.4721691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7922908) q[0];
sx q[0];
rz(-1.2783569) q[0];
sx q[0];
rz(1.6287623) q[0];
rz(0.28383645) q[1];
sx q[1];
rz(-2.2161039) q[1];
sx q[1];
rz(1.1121174) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6215206) q[0];
sx q[0];
rz(-3.1315098) q[0];
sx q[0];
rz(2.3501189) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4623423) q[2];
sx q[2];
rz(-1.4949833) q[2];
sx q[2];
rz(0.52559847) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.49431153) q[1];
sx q[1];
rz(-1.0930645) q[1];
sx q[1];
rz(1.7812438) q[1];
rz(-0.25896163) q[3];
sx q[3];
rz(-1.2185214) q[3];
sx q[3];
rz(1.6083628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6185559) q[2];
sx q[2];
rz(-1.7967537) q[2];
sx q[2];
rz(-1.1019361) q[2];
rz(-2.2526422) q[3];
sx q[3];
rz(-2.6933653) q[3];
sx q[3];
rz(0.86047188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3717644) q[0];
sx q[0];
rz(-0.17426057) q[0];
sx q[0];
rz(2.4523822) q[0];
rz(0.31461942) q[1];
sx q[1];
rz(-0.92051053) q[1];
sx q[1];
rz(1.4124195) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5310881) q[0];
sx q[0];
rz(-0.83685368) q[0];
sx q[0];
rz(0.42829163) q[0];
rz(-pi) q[1];
rz(-1.6666404) q[2];
sx q[2];
rz(-1.3587225) q[2];
sx q[2];
rz(-2.1532358) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7594193) q[1];
sx q[1];
rz(-1.628211) q[1];
sx q[1];
rz(0.38019726) q[1];
rz(-2.054677) q[3];
sx q[3];
rz(-1.3984507) q[3];
sx q[3];
rz(2.2215171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.41669258) q[2];
sx q[2];
rz(-0.56513864) q[2];
sx q[2];
rz(-0.55595428) q[2];
rz(-1.8703095) q[3];
sx q[3];
rz(-1.8794182) q[3];
sx q[3];
rz(-0.2909734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3274662) q[0];
sx q[0];
rz(-0.1411345) q[0];
sx q[0];
rz(-2.5471174) q[0];
rz(-2.5796083) q[1];
sx q[1];
rz(-0.85834208) q[1];
sx q[1];
rz(2.1655703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41244477) q[0];
sx q[0];
rz(-1.9322104) q[0];
sx q[0];
rz(-1.1114208) q[0];
x q[1];
rz(0.75404928) q[2];
sx q[2];
rz(-1.591914) q[2];
sx q[2];
rz(1.4625975) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1072717) q[1];
sx q[1];
rz(-1.8684505) q[1];
sx q[1];
rz(-0.68796449) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61686744) q[3];
sx q[3];
rz(-0.36412334) q[3];
sx q[3];
rz(-2.0193651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6550265) q[2];
sx q[2];
rz(-1.8883294) q[2];
sx q[2];
rz(-0.61069926) q[2];
rz(-0.30580172) q[3];
sx q[3];
rz(-2.1761201) q[3];
sx q[3];
rz(-1.8122199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3139451) q[0];
sx q[0];
rz(-0.018445404) q[0];
sx q[0];
rz(-1.6754643) q[0];
rz(0.77955359) q[1];
sx q[1];
rz(-1.164271) q[1];
sx q[1];
rz(-2.4028042) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8429564) q[0];
sx q[0];
rz(-1.8514575) q[0];
sx q[0];
rz(-0.21224169) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5814085) q[2];
sx q[2];
rz(-1.2640177) q[2];
sx q[2];
rz(-0.58394428) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.14911095) q[1];
sx q[1];
rz(-0.68435366) q[1];
sx q[1];
rz(-0.97196399) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8035369) q[3];
sx q[3];
rz(-0.49852405) q[3];
sx q[3];
rz(-0.62832181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5256727) q[2];
sx q[2];
rz(-1.9024666) q[2];
sx q[2];
rz(2.2789148) q[2];
rz(-0.45977965) q[3];
sx q[3];
rz(-1.6254057) q[3];
sx q[3];
rz(-1.9988683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.90726844) q[0];
sx q[0];
rz(-0.37721226) q[0];
sx q[0];
rz(-1.020485) q[0];
rz(0.057295784) q[1];
sx q[1];
rz(-1.4833996) q[1];
sx q[1];
rz(0.78757706) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60815281) q[0];
sx q[0];
rz(-0.69204563) q[0];
sx q[0];
rz(-0.37779053) q[0];
x q[1];
rz(2.8834613) q[2];
sx q[2];
rz(-1.5016593) q[2];
sx q[2];
rz(-0.11786945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.139442) q[1];
sx q[1];
rz(-2.6137335) q[1];
sx q[1];
rz(1.2107641) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.034463783) q[3];
sx q[3];
rz(-0.98317819) q[3];
sx q[3];
rz(1.8922488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.82137498) q[2];
sx q[2];
rz(-1.3221952) q[2];
sx q[2];
rz(0.58464948) q[2];
rz(-0.65230495) q[3];
sx q[3];
rz(-1.2290596) q[3];
sx q[3];
rz(1.339213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.24340165) q[0];
sx q[0];
rz(-1.3404055) q[0];
sx q[0];
rz(0.32824326) q[0];
rz(-0.39168656) q[1];
sx q[1];
rz(-1.990254) q[1];
sx q[1];
rz(1.4788871) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1680964) q[0];
sx q[0];
rz(-2.7910821) q[0];
sx q[0];
rz(-2.4309733) q[0];
rz(0.72490643) q[2];
sx q[2];
rz(-0.6577684) q[2];
sx q[2];
rz(-0.36054128) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6819671) q[1];
sx q[1];
rz(-0.29469583) q[1];
sx q[1];
rz(1.2380283) q[1];
rz(-pi) q[2];
rz(1.201466) q[3];
sx q[3];
rz(-1.7672667) q[3];
sx q[3];
rz(2.7854491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.11999232) q[2];
sx q[2];
rz(-1.3730717) q[2];
sx q[2];
rz(0.78869406) q[2];
rz(2.9005519) q[3];
sx q[3];
rz(-0.92625109) q[3];
sx q[3];
rz(2.327976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4929844) q[0];
sx q[0];
rz(-0.5383752) q[0];
sx q[0];
rz(0.96187821) q[0];
rz(-0.23421639) q[1];
sx q[1];
rz(-2.0030256) q[1];
sx q[1];
rz(-0.73807565) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86185011) q[0];
sx q[0];
rz(-1.8075917) q[0];
sx q[0];
rz(-0.40646942) q[0];
x q[1];
rz(2.5665087) q[2];
sx q[2];
rz(-0.7204537) q[2];
sx q[2];
rz(1.6937814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8705504) q[1];
sx q[1];
rz(-1.9890474) q[1];
sx q[1];
rz(2.8224432) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35154147) q[3];
sx q[3];
rz(-0.69720399) q[3];
sx q[3];
rz(-0.54920025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6568079) q[2];
sx q[2];
rz(-1.022555) q[2];
sx q[2];
rz(-1.2602932) q[2];
rz(-1.8079181) q[3];
sx q[3];
rz(-2.1402054) q[3];
sx q[3];
rz(-0.26237747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8986847) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(0.86724487) q[0];
rz(-1.0137089) q[1];
sx q[1];
rz(-1.8127245) q[1];
sx q[1];
rz(1.0702466) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27425048) q[0];
sx q[0];
rz(-0.19495067) q[0];
sx q[0];
rz(-2.250953) q[0];
rz(-0.81999166) q[2];
sx q[2];
rz(-2.242616) q[2];
sx q[2];
rz(-1.2206248) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6079191) q[1];
sx q[1];
rz(-1.9551139) q[1];
sx q[1];
rz(-0.55748765) q[1];
rz(3.002043) q[3];
sx q[3];
rz(-1.680239) q[3];
sx q[3];
rz(2.8982996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16537198) q[2];
sx q[2];
rz(-2.3207211) q[2];
sx q[2];
rz(1.9099859) q[2];
rz(1.5008789) q[3];
sx q[3];
rz(-0.21017635) q[3];
sx q[3];
rz(0.73178449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3043542) q[0];
sx q[0];
rz(-1.1501034) q[0];
sx q[0];
rz(1.3788086) q[0];
rz(-2.7267743) q[1];
sx q[1];
rz(-1.7638313) q[1];
sx q[1];
rz(1.5324963) q[1];
rz(-1.6258705) q[2];
sx q[2];
rz(-2.016042) q[2];
sx q[2];
rz(2.945937) q[2];
rz(-2.1555156) q[3];
sx q[3];
rz(-0.90860962) q[3];
sx q[3];
rz(-1.2398401) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

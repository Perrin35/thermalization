OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1202886) q[0];
sx q[0];
rz(-2.6074183) q[0];
sx q[0];
rz(2.0913273) q[0];
rz(-1.2644816) q[1];
sx q[1];
rz(-2.2883132) q[1];
sx q[1];
rz(2.5929911) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.825332) q[0];
sx q[0];
rz(-1.863593) q[0];
sx q[0];
rz(1.8789852) q[0];
rz(0.031096641) q[2];
sx q[2];
rz(-0.59760909) q[2];
sx q[2];
rz(0.55127599) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7214515) q[1];
sx q[1];
rz(-2.4597617) q[1];
sx q[1];
rz(-2.2504266) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2554507) q[3];
sx q[3];
rz(-1.7996801) q[3];
sx q[3];
rz(-1.6381519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41142622) q[2];
sx q[2];
rz(-0.41894087) q[2];
sx q[2];
rz(2.1298998) q[2];
rz(-0.88459477) q[3];
sx q[3];
rz(-1.1583637) q[3];
sx q[3];
rz(1.8959034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0737792) q[0];
sx q[0];
rz(-2.1429006) q[0];
sx q[0];
rz(-2.3538537) q[0];
rz(2.1630321) q[1];
sx q[1];
rz(-1.9906882) q[1];
sx q[1];
rz(1.8914793) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.387991) q[0];
sx q[0];
rz(-1.6785673) q[0];
sx q[0];
rz(-2.0600256) q[0];
rz(-0.97888005) q[2];
sx q[2];
rz(-0.47431163) q[2];
sx q[2];
rz(1.7537376) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3652894) q[1];
sx q[1];
rz(-2.3765916) q[1];
sx q[1];
rz(1.6056152) q[1];
x q[2];
rz(2.3552966) q[3];
sx q[3];
rz(-2.8350787) q[3];
sx q[3];
rz(1.802352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1193739) q[2];
sx q[2];
rz(-1.4639857) q[2];
sx q[2];
rz(-3.019943) q[2];
rz(0.71074784) q[3];
sx q[3];
rz(-0.23356479) q[3];
sx q[3];
rz(-0.60230437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0372593) q[0];
sx q[0];
rz(-2.4132044) q[0];
sx q[0];
rz(-2.4816568) q[0];
rz(0.51689369) q[1];
sx q[1];
rz(-2.4362502) q[1];
sx q[1];
rz(-1.0924115) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79688841) q[0];
sx q[0];
rz(-1.832334) q[0];
sx q[0];
rz(-1.198223) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2159816) q[2];
sx q[2];
rz(-0.72941226) q[2];
sx q[2];
rz(1.6075921) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1221083) q[1];
sx q[1];
rz(-2.1017764) q[1];
sx q[1];
rz(2.2506258) q[1];
rz(-pi) q[2];
rz(-1.3140244) q[3];
sx q[3];
rz(-2.2035172) q[3];
sx q[3];
rz(2.4420121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2163781) q[2];
sx q[2];
rz(-0.34319147) q[2];
sx q[2];
rz(1.421831) q[2];
rz(0.4839932) q[3];
sx q[3];
rz(-1.4839987) q[3];
sx q[3];
rz(1.7234195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5487109) q[0];
sx q[0];
rz(-0.084267862) q[0];
sx q[0];
rz(-1.3053869) q[0];
rz(-0.0072172324) q[1];
sx q[1];
rz(-2.9199298) q[1];
sx q[1];
rz(-0.73297393) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78295499) q[0];
sx q[0];
rz(-1.6894537) q[0];
sx q[0];
rz(-1.4374742) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4568366) q[2];
sx q[2];
rz(-1.2780398) q[2];
sx q[2];
rz(1.2639015) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.32506714) q[1];
sx q[1];
rz(-0.96615929) q[1];
sx q[1];
rz(-0.52765347) q[1];
x q[2];
rz(0.72317883) q[3];
sx q[3];
rz(-1.6537602) q[3];
sx q[3];
rz(-1.7570329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8492154) q[2];
sx q[2];
rz(-1.7376309) q[2];
sx q[2];
rz(2.2461069) q[2];
rz(2.605947) q[3];
sx q[3];
rz(-1.5629385) q[3];
sx q[3];
rz(0.061804684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2484922) q[0];
sx q[0];
rz(-0.19344261) q[0];
sx q[0];
rz(2.6336811) q[0];
rz(1.4503362) q[1];
sx q[1];
rz(-1.0117057) q[1];
sx q[1];
rz(1.7844261) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0944509) q[0];
sx q[0];
rz(-3.019382) q[0];
sx q[0];
rz(-0.41441865) q[0];
rz(-1.9903899) q[2];
sx q[2];
rz(-1.9997678) q[2];
sx q[2];
rz(1.8544514) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0101945) q[1];
sx q[1];
rz(-1.6784188) q[1];
sx q[1];
rz(-1.9626928) q[1];
rz(-pi) q[2];
rz(-2.761854) q[3];
sx q[3];
rz(-1.0061227) q[3];
sx q[3];
rz(-2.5329451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7148529) q[2];
sx q[2];
rz(-1.667495) q[2];
sx q[2];
rz(-2.0392141) q[2];
rz(-1.1357931) q[3];
sx q[3];
rz(-1.9250684) q[3];
sx q[3];
rz(-0.77146012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8847467) q[0];
sx q[0];
rz(-0.9698292) q[0];
sx q[0];
rz(2.2127175) q[0];
rz(-0.38201395) q[1];
sx q[1];
rz(-0.99645749) q[1];
sx q[1];
rz(-1.4487723) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8162874) q[0];
sx q[0];
rz(-1.0950118) q[0];
sx q[0];
rz(1.2289117) q[0];
rz(1.5097627) q[2];
sx q[2];
rz(-1.0866797) q[2];
sx q[2];
rz(-2.8060437) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12940059) q[1];
sx q[1];
rz(-1.4621648) q[1];
sx q[1];
rz(-1.9700178) q[1];
rz(-pi) q[2];
rz(-2.9690456) q[3];
sx q[3];
rz(-1.2013519) q[3];
sx q[3];
rz(-1.2404694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.53943071) q[2];
sx q[2];
rz(-2.7950037) q[2];
sx q[2];
rz(-0.76751417) q[2];
rz(-1.2933939) q[3];
sx q[3];
rz(-0.69197217) q[3];
sx q[3];
rz(-2.9366711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9354189) q[0];
sx q[0];
rz(-2.967301) q[0];
sx q[0];
rz(0.25099227) q[0];
rz(0.17768606) q[1];
sx q[1];
rz(-1.6975941) q[1];
sx q[1];
rz(0.65690717) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0278695) q[0];
sx q[0];
rz(-1.3355458) q[0];
sx q[0];
rz(2.9293438) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20505242) q[2];
sx q[2];
rz(-0.66407138) q[2];
sx q[2];
rz(2.1560706) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.71437876) q[1];
sx q[1];
rz(-1.8698815) q[1];
sx q[1];
rz(-2.7370791) q[1];
x q[2];
rz(-2.4961996) q[3];
sx q[3];
rz(-1.925996) q[3];
sx q[3];
rz(-1.3076289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80829197) q[2];
sx q[2];
rz(-0.57922816) q[2];
sx q[2];
rz(-0.91326886) q[2];
rz(2.463786) q[3];
sx q[3];
rz(-1.6401688) q[3];
sx q[3];
rz(2.138413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0388357) q[0];
sx q[0];
rz(-1.1433733) q[0];
sx q[0];
rz(-0.58297771) q[0];
rz(-1.673117) q[1];
sx q[1];
rz(-2.1491094) q[1];
sx q[1];
rz(-2.800422) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4214246) q[0];
sx q[0];
rz(-0.20381323) q[0];
sx q[0];
rz(-1.782062) q[0];
rz(-2.274373) q[2];
sx q[2];
rz(-0.3647621) q[2];
sx q[2];
rz(0.64389578) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9205089) q[1];
sx q[1];
rz(-2.608641) q[1];
sx q[1];
rz(0.78842649) q[1];
rz(-2.736015) q[3];
sx q[3];
rz(-2.2080407) q[3];
sx q[3];
rz(3.0831856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1950281) q[2];
sx q[2];
rz(-0.82411689) q[2];
sx q[2];
rz(0.80671802) q[2];
rz(0.45754704) q[3];
sx q[3];
rz(-2.1541607) q[3];
sx q[3];
rz(-2.257982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.44613999) q[0];
sx q[0];
rz(-2.0841632) q[0];
sx q[0];
rz(2.9651508) q[0];
rz(0.31013075) q[1];
sx q[1];
rz(-1.4896432) q[1];
sx q[1];
rz(-2.2055221) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98956489) q[0];
sx q[0];
rz(-1.7195722) q[0];
sx q[0];
rz(-0.29971931) q[0];
x q[1];
rz(-2.0240729) q[2];
sx q[2];
rz(-1.0220811) q[2];
sx q[2];
rz(2.0224151) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1470799) q[1];
sx q[1];
rz(-1.4777061) q[1];
sx q[1];
rz(-1.9398111) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0210439) q[3];
sx q[3];
rz(-2.0417739) q[3];
sx q[3];
rz(-0.69191832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0940242) q[2];
sx q[2];
rz(-1.3022283) q[2];
sx q[2];
rz(0.0085208323) q[2];
rz(0.73355567) q[3];
sx q[3];
rz(-2.3365884) q[3];
sx q[3];
rz(3.0146397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9200639) q[0];
sx q[0];
rz(-2.7286752) q[0];
sx q[0];
rz(2.371696) q[0];
rz(0.13433111) q[1];
sx q[1];
rz(-0.7901935) q[1];
sx q[1];
rz(-1.92164) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3111356) q[0];
sx q[0];
rz(-0.38549462) q[0];
sx q[0];
rz(2.4804082) q[0];
x q[1];
rz(2.7292615) q[2];
sx q[2];
rz(-2.4930232) q[2];
sx q[2];
rz(-2.1438053) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.091948) q[1];
sx q[1];
rz(-1.1941027) q[1];
sx q[1];
rz(2.6833181) q[1];
rz(3.0705922) q[3];
sx q[3];
rz(-1.7216847) q[3];
sx q[3];
rz(1.7689266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.28025815) q[2];
sx q[2];
rz(-0.65797776) q[2];
sx q[2];
rz(-2.3560143) q[2];
rz(-2.806459) q[3];
sx q[3];
rz(-0.98888713) q[3];
sx q[3];
rz(-2.3194763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32421865) q[0];
sx q[0];
rz(-0.86295177) q[0];
sx q[0];
rz(-1.2865768) q[0];
rz(2.4556976) q[1];
sx q[1];
rz(-2.5527725) q[1];
sx q[1];
rz(1.7935161) q[1];
rz(1.5625207) q[2];
sx q[2];
rz(-2.5604421) q[2];
sx q[2];
rz(1.5887518) q[2];
rz(-1.8635345) q[3];
sx q[3];
rz(-1.5354029) q[3];
sx q[3];
rz(-1.7913747) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

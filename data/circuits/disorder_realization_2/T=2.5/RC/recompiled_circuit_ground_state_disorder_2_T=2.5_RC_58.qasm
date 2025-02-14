OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.34831369) q[0];
sx q[0];
rz(-1.4077185) q[0];
sx q[0];
rz(-0.044535927) q[0];
rz(2.0308004) q[1];
sx q[1];
rz(-1.9171311) q[1];
sx q[1];
rz(0.78723025) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6605741) q[0];
sx q[0];
rz(-2.3957167) q[0];
sx q[0];
rz(-0.77156271) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2455315) q[2];
sx q[2];
rz(-1.6084387) q[2];
sx q[2];
rz(-1.4311439) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8833784) q[1];
sx q[1];
rz(-1.5014663) q[1];
sx q[1];
rz(-0.99154559) q[1];
rz(0.23870339) q[3];
sx q[3];
rz(-1.9145962) q[3];
sx q[3];
rz(-1.3581004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3945776) q[2];
sx q[2];
rz(-2.0614823) q[2];
sx q[2];
rz(2.3874095) q[2];
rz(1.7265823) q[3];
sx q[3];
rz(-0.49608803) q[3];
sx q[3];
rz(2.8126341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0542145) q[0];
sx q[0];
rz(-2.7890451) q[0];
sx q[0];
rz(1.8074328) q[0];
rz(1.2913903) q[1];
sx q[1];
rz(-1.038237) q[1];
sx q[1];
rz(-2.2659567) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9459045) q[0];
sx q[0];
rz(-2.2155846) q[0];
sx q[0];
rz(0.88627215) q[0];
x q[1];
rz(0.062280999) q[2];
sx q[2];
rz(-1.0750293) q[2];
sx q[2];
rz(2.371126) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.27442447) q[1];
sx q[1];
rz(-0.53571415) q[1];
sx q[1];
rz(0.97656004) q[1];
rz(-pi) q[2];
rz(-0.61296566) q[3];
sx q[3];
rz(-2.7904982) q[3];
sx q[3];
rz(-2.8610817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4448173) q[2];
sx q[2];
rz(-2.4397218) q[2];
sx q[2];
rz(0.87265054) q[2];
rz(-0.78898346) q[3];
sx q[3];
rz(-0.9011457) q[3];
sx q[3];
rz(2.9130329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2149684) q[0];
sx q[0];
rz(-1.3525532) q[0];
sx q[0];
rz(-0.0080000814) q[0];
rz(-2.3232715) q[1];
sx q[1];
rz(-0.74940959) q[1];
sx q[1];
rz(1.2368894) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0950355) q[0];
sx q[0];
rz(-2.3716784) q[0];
sx q[0];
rz(-2.0410935) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90571442) q[2];
sx q[2];
rz(-0.53230282) q[2];
sx q[2];
rz(1.5168845) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1874967) q[1];
sx q[1];
rz(-1.8774697) q[1];
sx q[1];
rz(-0.56879725) q[1];
x q[2];
rz(-2.9048484) q[3];
sx q[3];
rz(-1.4567353) q[3];
sx q[3];
rz(-0.1027861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.41039738) q[2];
sx q[2];
rz(-3.0072913) q[2];
sx q[2];
rz(-0.94089874) q[2];
rz(1.1040374) q[3];
sx q[3];
rz(-1.7887812) q[3];
sx q[3];
rz(1.1436536) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10821548) q[0];
sx q[0];
rz(-0.5431076) q[0];
sx q[0];
rz(-0.23750842) q[0];
rz(-0.6595276) q[1];
sx q[1];
rz(-0.57412761) q[1];
sx q[1];
rz(3.1245756) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7369923) q[0];
sx q[0];
rz(-1.7395818) q[0];
sx q[0];
rz(-0.048860839) q[0];
rz(0.44397644) q[2];
sx q[2];
rz(-1.7594565) q[2];
sx q[2];
rz(0.89934811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.58910546) q[1];
sx q[1];
rz(-2.1334339) q[1];
sx q[1];
rz(-1.6943114) q[1];
rz(-pi) q[2];
rz(-0.83570273) q[3];
sx q[3];
rz(-0.91360211) q[3];
sx q[3];
rz(-0.14086715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.90618769) q[2];
sx q[2];
rz(-1.9402639) q[2];
sx q[2];
rz(-0.91442937) q[2];
rz(0.36214456) q[3];
sx q[3];
rz(-2.3320964) q[3];
sx q[3];
rz(-2.5563498) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2528766) q[0];
sx q[0];
rz(-1.7676366) q[0];
sx q[0];
rz(3.0294898) q[0];
rz(1.0653227) q[1];
sx q[1];
rz(-1.0708829) q[1];
sx q[1];
rz(-2.0723453) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1493268) q[0];
sx q[0];
rz(-0.42851028) q[0];
sx q[0];
rz(-0.70121588) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4222071) q[2];
sx q[2];
rz(-2.4771538) q[2];
sx q[2];
rz(1.7094572) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9578609) q[1];
sx q[1];
rz(-1.2535447) q[1];
sx q[1];
rz(-0.31362335) q[1];
x q[2];
rz(-1.3773243) q[3];
sx q[3];
rz(-2.5335823) q[3];
sx q[3];
rz(1.6093169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7142882) q[2];
sx q[2];
rz(-0.71791831) q[2];
sx q[2];
rz(2.60738) q[2];
rz(0.30321768) q[3];
sx q[3];
rz(-1.5568308) q[3];
sx q[3];
rz(1.5155972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41820207) q[0];
sx q[0];
rz(-1.8573107) q[0];
sx q[0];
rz(-2.2499625) q[0];
rz(-1.3806237) q[1];
sx q[1];
rz(-1.9386407) q[1];
sx q[1];
rz(-2.2183529) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6398481) q[0];
sx q[0];
rz(-3.0303133) q[0];
sx q[0];
rz(-2.6529045) q[0];
x q[1];
rz(2.5371952) q[2];
sx q[2];
rz(-2.6816419) q[2];
sx q[2];
rz(-0.15397554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.90369836) q[1];
sx q[1];
rz(-0.35939068) q[1];
sx q[1];
rz(-0.87025799) q[1];
x q[2];
rz(2.6232016) q[3];
sx q[3];
rz(-0.71037358) q[3];
sx q[3];
rz(2.667676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0455857) q[2];
sx q[2];
rz(-0.43232375) q[2];
sx q[2];
rz(-1.2621137) q[2];
rz(-1.1612085) q[3];
sx q[3];
rz(-1.592417) q[3];
sx q[3];
rz(1.8956634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68814174) q[0];
sx q[0];
rz(-0.83204404) q[0];
sx q[0];
rz(-1.6081109) q[0];
rz(0.43117943) q[1];
sx q[1];
rz(-2.2240413) q[1];
sx q[1];
rz(1.4088438) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20074318) q[0];
sx q[0];
rz(-0.701776) q[0];
sx q[0];
rz(0.45325847) q[0];
rz(-pi) q[1];
rz(-1.9555075) q[2];
sx q[2];
rz(-1.6530347) q[2];
sx q[2];
rz(-2.4599894) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3235493) q[1];
sx q[1];
rz(-2.9406639) q[1];
sx q[1];
rz(-2.4646321) q[1];
rz(-1.7944144) q[3];
sx q[3];
rz(-2.6287492) q[3];
sx q[3];
rz(-1.4033069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2048753) q[2];
sx q[2];
rz(-1.1254346) q[2];
sx q[2];
rz(1.8219061) q[2];
rz(0.034505757) q[3];
sx q[3];
rz(-1.6104108) q[3];
sx q[3];
rz(1.3694793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7337604) q[0];
sx q[0];
rz(-0.5540846) q[0];
sx q[0];
rz(1.9246509) q[0];
rz(0.92562428) q[1];
sx q[1];
rz(-1.3652912) q[1];
sx q[1];
rz(-1.9409174) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.417765) q[0];
sx q[0];
rz(-2.1795131) q[0];
sx q[0];
rz(0.62792741) q[0];
rz(-pi) q[1];
rz(-0.56270069) q[2];
sx q[2];
rz(-0.68018736) q[2];
sx q[2];
rz(-2.0923751) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6457121) q[1];
sx q[1];
rz(-2.3420942) q[1];
sx q[1];
rz(1.6281566) q[1];
x q[2];
rz(0.36603277) q[3];
sx q[3];
rz(-1.9683553) q[3];
sx q[3];
rz(1.9983069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.22244464) q[2];
sx q[2];
rz(-1.6067952) q[2];
sx q[2];
rz(1.0486802) q[2];
rz(2.9564296) q[3];
sx q[3];
rz(-1.984963) q[3];
sx q[3];
rz(0.10704253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1222526) q[0];
sx q[0];
rz(-0.66146079) q[0];
sx q[0];
rz(-0.81073236) q[0];
rz(-0.73075378) q[1];
sx q[1];
rz(-2.0850875) q[1];
sx q[1];
rz(2.1810541) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4385242) q[0];
sx q[0];
rz(-1.5674066) q[0];
sx q[0];
rz(1.5717713) q[0];
rz(-pi) q[1];
rz(-1.0658747) q[2];
sx q[2];
rz(-1.9532579) q[2];
sx q[2];
rz(1.7590211) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5170578) q[1];
sx q[1];
rz(-1.0427999) q[1];
sx q[1];
rz(-2.8370884) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1436362) q[3];
sx q[3];
rz(-1.8064326) q[3];
sx q[3];
rz(2.4164661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1824128) q[2];
sx q[2];
rz(-1.7824087) q[2];
sx q[2];
rz(-1.6020927) q[2];
rz(1.0473853) q[3];
sx q[3];
rz(-1.7644707) q[3];
sx q[3];
rz(-1.876095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78028107) q[0];
sx q[0];
rz(-2.2762716) q[0];
sx q[0];
rz(2.5541232) q[0];
rz(2.7777708) q[1];
sx q[1];
rz(-1.1452585) q[1];
sx q[1];
rz(2.2344373) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5840549) q[0];
sx q[0];
rz(-1.4901925) q[0];
sx q[0];
rz(0.35146468) q[0];
rz(2.6689879) q[2];
sx q[2];
rz(-1.9197004) q[2];
sx q[2];
rz(-2.7965703) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2550019) q[1];
sx q[1];
rz(-1.1215409) q[1];
sx q[1];
rz(1.3372861) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.206622) q[3];
sx q[3];
rz(-1.6505906) q[3];
sx q[3];
rz(-1.7062812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5017447) q[2];
sx q[2];
rz(-3.0526243) q[2];
sx q[2];
rz(-2.7196344) q[2];
rz(0.0095327775) q[3];
sx q[3];
rz(-1.5744753) q[3];
sx q[3];
rz(0.52614051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.52787732) q[0];
sx q[0];
rz(-1.2953225) q[0];
sx q[0];
rz(-3.0184826) q[0];
rz(-0.70115024) q[1];
sx q[1];
rz(-1.2260561) q[1];
sx q[1];
rz(-1.4969926) q[1];
rz(0.43177615) q[2];
sx q[2];
rz(-2.3974621) q[2];
sx q[2];
rz(1.5324788) q[2];
rz(-0.27746874) q[3];
sx q[3];
rz(-1.1744842) q[3];
sx q[3];
rz(-0.13406772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

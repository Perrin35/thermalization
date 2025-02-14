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
rz(2.9313791) q[0];
sx q[0];
rz(5.9271521) q[0];
sx q[0];
rz(7.9071101) q[0];
rz(-1.3104982) q[1];
sx q[1];
rz(-0.85135353) q[1];
sx q[1];
rz(0.89827615) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7018902) q[0];
sx q[0];
rz(-2.1680346) q[0];
sx q[0];
rz(2.7053614) q[0];
x q[1];
rz(2.8509628) q[2];
sx q[2];
rz(-2.0465474) q[2];
sx q[2];
rz(-1.0889183) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5772669) q[1];
sx q[1];
rz(-2.942406) q[1];
sx q[1];
rz(0.41173068) q[1];
rz(1.6226852) q[3];
sx q[3];
rz(-0.16434419) q[3];
sx q[3];
rz(-0.60620327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5400759) q[2];
sx q[2];
rz(-2.8105152) q[2];
sx q[2];
rz(2.4430821) q[2];
rz(-1.4861594) q[3];
sx q[3];
rz(-2.555002) q[3];
sx q[3];
rz(-1.9894039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0965213) q[0];
sx q[0];
rz(-0.69913816) q[0];
sx q[0];
rz(-0.29749468) q[0];
rz(-2.7104764) q[1];
sx q[1];
rz(-2.2747206) q[1];
sx q[1];
rz(-1.8284304) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24475141) q[0];
sx q[0];
rz(-2.2479821) q[0];
sx q[0];
rz(-2.0170101) q[0];
rz(-0.044017893) q[2];
sx q[2];
rz(-1.9761397) q[2];
sx q[2];
rz(2.8716033) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0521691) q[1];
sx q[1];
rz(-2.3681297) q[1];
sx q[1];
rz(-1.5816342) q[1];
x q[2];
rz(1.6974988) q[3];
sx q[3];
rz(-1.8167348) q[3];
sx q[3];
rz(-1.9501231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7737274) q[2];
sx q[2];
rz(-0.44168681) q[2];
sx q[2];
rz(2.3750677) q[2];
rz(0.68486989) q[3];
sx q[3];
rz(-1.6133512) q[3];
sx q[3];
rz(-0.49055704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3016475) q[0];
sx q[0];
rz(-1.8122346) q[0];
sx q[0];
rz(0.30461052) q[0];
rz(1.0003164) q[1];
sx q[1];
rz(-1.1023003) q[1];
sx q[1];
rz(-0.89835483) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4173664) q[0];
sx q[0];
rz(-2.6018916) q[0];
sx q[0];
rz(1.4090572) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9581306) q[2];
sx q[2];
rz(-2.4011321) q[2];
sx q[2];
rz(-1.0473521) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.5872941) q[1];
sx q[1];
rz(-1.0815433) q[1];
sx q[1];
rz(2.4913408) q[1];
rz(-pi) q[2];
x q[2];
rz(0.091354185) q[3];
sx q[3];
rz(-2.0915789) q[3];
sx q[3];
rz(-2.6798673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1343214) q[2];
sx q[2];
rz(-1.4192702) q[2];
sx q[2];
rz(0.32709861) q[2];
rz(1.5056115) q[3];
sx q[3];
rz(-1.4666731) q[3];
sx q[3];
rz(1.3350284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026767749) q[0];
sx q[0];
rz(-0.46634316) q[0];
sx q[0];
rz(2.603671) q[0];
rz(2.3665358) q[1];
sx q[1];
rz(-1.8102976) q[1];
sx q[1];
rz(-2.7937826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8224768) q[0];
sx q[0];
rz(-2.4317784) q[0];
sx q[0];
rz(-2.3869724) q[0];
rz(-pi) q[1];
rz(1.6762907) q[2];
sx q[2];
rz(-1.4660133) q[2];
sx q[2];
rz(1.8035165) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.28827846) q[1];
sx q[1];
rz(-2.2396259) q[1];
sx q[1];
rz(2.1200256) q[1];
x q[2];
rz(-1.1775753) q[3];
sx q[3];
rz(-1.4054417) q[3];
sx q[3];
rz(2.4158784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9378918) q[2];
sx q[2];
rz(-0.72526473) q[2];
sx q[2];
rz(0.66246486) q[2];
rz(2.0857816) q[3];
sx q[3];
rz(-2.8231088) q[3];
sx q[3];
rz(1.2663579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98589677) q[0];
sx q[0];
rz(-0.53845423) q[0];
sx q[0];
rz(2.1180617) q[0];
rz(-1.0515949) q[1];
sx q[1];
rz(-1.6513499) q[1];
sx q[1];
rz(0.016062707) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15323205) q[0];
sx q[0];
rz(-1.4356336) q[0];
sx q[0];
rz(-1.022781) q[0];
rz(-pi) q[1];
rz(0.73170029) q[2];
sx q[2];
rz(-0.83544399) q[2];
sx q[2];
rz(1.4994061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2019883) q[1];
sx q[1];
rz(-1.3496205) q[1];
sx q[1];
rz(-0.53316922) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2993126) q[3];
sx q[3];
rz(-1.5010271) q[3];
sx q[3];
rz(1.8673563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2992799) q[2];
sx q[2];
rz(-0.88025847) q[2];
sx q[2];
rz(2.9006145) q[2];
rz(2.0321417) q[3];
sx q[3];
rz(-1.8532608) q[3];
sx q[3];
rz(-0.39458767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0146765) q[0];
sx q[0];
rz(-2.3052445) q[0];
sx q[0];
rz(0.29689223) q[0];
rz(-1.1709921) q[1];
sx q[1];
rz(-1.3208656) q[1];
sx q[1];
rz(-2.6062633) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4142617) q[0];
sx q[0];
rz(-2.8777661) q[0];
sx q[0];
rz(0.78049941) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2601603) q[2];
sx q[2];
rz(-0.77858965) q[2];
sx q[2];
rz(-2.1354577) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8187021) q[1];
sx q[1];
rz(-1.358958) q[1];
sx q[1];
rz(2.0634275) q[1];
rz(-pi) q[2];
rz(2.2933741) q[3];
sx q[3];
rz(-2.9322538) q[3];
sx q[3];
rz(2.6959553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.31763306) q[2];
sx q[2];
rz(-2.3667658) q[2];
sx q[2];
rz(-0.63014692) q[2];
rz(1.0050425) q[3];
sx q[3];
rz(-2.9229087) q[3];
sx q[3];
rz(-2.4439243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.180535) q[0];
sx q[0];
rz(-3.0735885) q[0];
sx q[0];
rz(1.2766174) q[0];
rz(2.0969773) q[1];
sx q[1];
rz(-1.4199384) q[1];
sx q[1];
rz(0.23342625) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40204321) q[0];
sx q[0];
rz(-1.528571) q[0];
sx q[0];
rz(0.84997886) q[0];
x q[1];
rz(-2.4996148) q[2];
sx q[2];
rz(-1.3228088) q[2];
sx q[2];
rz(-0.61680142) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0835159) q[1];
sx q[1];
rz(-0.74878557) q[1];
sx q[1];
rz(-3.0850436) q[1];
rz(-pi) q[2];
rz(-0.57817187) q[3];
sx q[3];
rz(-1.626818) q[3];
sx q[3];
rz(-0.094464555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.6650247) q[2];
sx q[2];
rz(-1.8835521) q[2];
sx q[2];
rz(0.27102077) q[2];
rz(-2.792231) q[3];
sx q[3];
rz(-2.2237399) q[3];
sx q[3];
rz(-0.57749403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2006328) q[0];
sx q[0];
rz(-0.9372434) q[0];
sx q[0];
rz(-1.3295133) q[0];
rz(0.26893523) q[1];
sx q[1];
rz(-2.4766998) q[1];
sx q[1];
rz(-1.3806794) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5912) q[0];
sx q[0];
rz(-1.368528) q[0];
sx q[0];
rz(1.8254542) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5538832) q[2];
sx q[2];
rz(-2.1667987) q[2];
sx q[2];
rz(2.2600095) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0448398) q[1];
sx q[1];
rz(-2.1659454) q[1];
sx q[1];
rz(-2.9446141) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21694048) q[3];
sx q[3];
rz(-0.57121459) q[3];
sx q[3];
rz(1.9344714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5844476) q[2];
sx q[2];
rz(-0.86981213) q[2];
sx q[2];
rz(-2.4073041) q[2];
rz(1.0555438) q[3];
sx q[3];
rz(-1.5493834) q[3];
sx q[3];
rz(-2.1671104) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4771117) q[0];
sx q[0];
rz(-2.8900914) q[0];
sx q[0];
rz(-0.86790458) q[0];
rz(1.3033298) q[1];
sx q[1];
rz(-1.5918599) q[1];
sx q[1];
rz(2.910639) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.445914) q[0];
sx q[0];
rz(-1.1584846) q[0];
sx q[0];
rz(-0.59230174) q[0];
rz(-pi) q[1];
rz(-2.852298) q[2];
sx q[2];
rz(-2.691011) q[2];
sx q[2];
rz(0.10403004) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8764774) q[1];
sx q[1];
rz(-2.6314221) q[1];
sx q[1];
rz(-2.6835346) q[1];
rz(2.4695314) q[3];
sx q[3];
rz(-1.5717744) q[3];
sx q[3];
rz(0.031571183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2086198) q[2];
sx q[2];
rz(-1.4459556) q[2];
sx q[2];
rz(-2.7033973) q[2];
rz(0.70133251) q[3];
sx q[3];
rz(-1.067679) q[3];
sx q[3];
rz(-0.35593885) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0375131) q[0];
sx q[0];
rz(-2.6689745) q[0];
sx q[0];
rz(-1.0802826) q[0];
rz(-2.089962) q[1];
sx q[1];
rz(-2.0104505) q[1];
sx q[1];
rz(2.0463321) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8488377) q[0];
sx q[0];
rz(-2.0803323) q[0];
sx q[0];
rz(1.2408907) q[0];
rz(3.0313086) q[2];
sx q[2];
rz(-0.37026893) q[2];
sx q[2];
rz(-0.25843378) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7113774) q[1];
sx q[1];
rz(-1.681455) q[1];
sx q[1];
rz(-2.2994413) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5445537) q[3];
sx q[3];
rz(-1.9970915) q[3];
sx q[3];
rz(0.36620127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0206535) q[2];
sx q[2];
rz(-1.8755269) q[2];
sx q[2];
rz(-2.2361163) q[2];
rz(-2.4261273) q[3];
sx q[3];
rz(-0.72532907) q[3];
sx q[3];
rz(-2.4874617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41659551) q[0];
sx q[0];
rz(-1.7323957) q[0];
sx q[0];
rz(0.61687627) q[0];
rz(-2.0116518) q[1];
sx q[1];
rz(-1.2810974) q[1];
sx q[1];
rz(1.7249736) q[1];
rz(0.62192179) q[2];
sx q[2];
rz(-2.1991232) q[2];
sx q[2];
rz(-0.029673619) q[2];
rz(0.15177095) q[3];
sx q[3];
rz(-1.3625154) q[3];
sx q[3];
rz(-1.2073928) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

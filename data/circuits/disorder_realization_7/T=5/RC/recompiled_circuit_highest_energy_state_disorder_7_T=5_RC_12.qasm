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
rz(-3.0149674) q[0];
sx q[0];
rz(-1.5746483) q[0];
sx q[0];
rz(-2.6259165) q[0];
rz(0.22817837) q[1];
sx q[1];
rz(-0.74695865) q[1];
sx q[1];
rz(0.42626122) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3328457) q[0];
sx q[0];
rz(-1.3984183) q[0];
sx q[0];
rz(-2.3461949) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49246712) q[2];
sx q[2];
rz(-2.5829007) q[2];
sx q[2];
rz(-0.44670263) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30657739) q[1];
sx q[1];
rz(-1.5131498) q[1];
sx q[1];
rz(-0.60204864) q[1];
rz(-pi) q[2];
rz(0.32483806) q[3];
sx q[3];
rz(-1.0110613) q[3];
sx q[3];
rz(-1.4046275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8969741) q[2];
sx q[2];
rz(-1.3099058) q[2];
sx q[2];
rz(-2.8986325) q[2];
rz(-2.7729559) q[3];
sx q[3];
rz(-0.60550767) q[3];
sx q[3];
rz(-1.2322371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31545562) q[0];
sx q[0];
rz(-0.22268) q[0];
sx q[0];
rz(-2.7742703) q[0];
rz(-0.79633725) q[1];
sx q[1];
rz(-2.0834736) q[1];
sx q[1];
rz(-0.64250362) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87602654) q[0];
sx q[0];
rz(-2.4867704) q[0];
sx q[0];
rz(-0.73237082) q[0];
rz(-2.5902469) q[2];
sx q[2];
rz(-1.2228106) q[2];
sx q[2];
rz(-2.6642193) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5629329) q[1];
sx q[1];
rz(-2.0173434) q[1];
sx q[1];
rz(0.62942735) q[1];
rz(-2.0525371) q[3];
sx q[3];
rz(-1.6828487) q[3];
sx q[3];
rz(0.58603906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.6392886) q[2];
sx q[2];
rz(-2.2059811) q[2];
sx q[2];
rz(-1.3703692) q[2];
rz(0.227452) q[3];
sx q[3];
rz(-1.8919614) q[3];
sx q[3];
rz(1.1915709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3962536) q[0];
sx q[0];
rz(-0.17530137) q[0];
sx q[0];
rz(-1.9434209) q[0];
rz(-2.0612969) q[1];
sx q[1];
rz(-2.9273169) q[1];
sx q[1];
rz(0.094873039) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0563597) q[0];
sx q[0];
rz(-1.0329536) q[0];
sx q[0];
rz(0.30098029) q[0];
x q[1];
rz(0.40887654) q[2];
sx q[2];
rz(-0.95252242) q[2];
sx q[2];
rz(-0.28291075) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2641622) q[1];
sx q[1];
rz(-2.1070711) q[1];
sx q[1];
rz(3.0939878) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1867939) q[3];
sx q[3];
rz(-1.3130762) q[3];
sx q[3];
rz(-2.6741762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0692811) q[2];
sx q[2];
rz(-2.5371964) q[2];
sx q[2];
rz(-2.7351232) q[2];
rz(-1.9629924) q[3];
sx q[3];
rz(-1.7013763) q[3];
sx q[3];
rz(1.1627722) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56226319) q[0];
sx q[0];
rz(-3.0070906) q[0];
sx q[0];
rz(-0.33600268) q[0];
rz(-2.621189) q[1];
sx q[1];
rz(-0.86417472) q[1];
sx q[1];
rz(-0.10428183) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9606303) q[0];
sx q[0];
rz(-1.917836) q[0];
sx q[0];
rz(-0.35023679) q[0];
rz(1.5849131) q[2];
sx q[2];
rz(-1.7272564) q[2];
sx q[2];
rz(-0.9597646) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9469918) q[1];
sx q[1];
rz(-0.90819383) q[1];
sx q[1];
rz(-0.47269447) q[1];
rz(1.453425) q[3];
sx q[3];
rz(-2.0214098) q[3];
sx q[3];
rz(2.4291458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5229554) q[2];
sx q[2];
rz(-2.0433661) q[2];
sx q[2];
rz(-1.9970419) q[2];
rz(-2.5942904) q[3];
sx q[3];
rz(-1.2283044) q[3];
sx q[3];
rz(-1.087629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2523786) q[0];
sx q[0];
rz(-0.19110282) q[0];
sx q[0];
rz(0.55437535) q[0];
rz(0.91122183) q[1];
sx q[1];
rz(-1.2634042) q[1];
sx q[1];
rz(2.8401781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2269991) q[0];
sx q[0];
rz(-2.1560139) q[0];
sx q[0];
rz(2.8666158) q[0];
rz(0.47573702) q[2];
sx q[2];
rz(-0.43125501) q[2];
sx q[2];
rz(0.83275821) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0423066) q[1];
sx q[1];
rz(-1.1908682) q[1];
sx q[1];
rz(-2.9126588) q[1];
rz(-pi) q[2];
rz(-0.88121342) q[3];
sx q[3];
rz(-0.47166079) q[3];
sx q[3];
rz(2.415433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1201799) q[2];
sx q[2];
rz(-2.9605949) q[2];
sx q[2];
rz(-0.0069228355) q[2];
rz(2.210468) q[3];
sx q[3];
rz(-2.0102863) q[3];
sx q[3];
rz(-1.1091703) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50865006) q[0];
sx q[0];
rz(-0.83704346) q[0];
sx q[0];
rz(1.5337926) q[0];
rz(-2.4389229) q[1];
sx q[1];
rz(-2.6103795) q[1];
sx q[1];
rz(-0.30219561) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93389171) q[0];
sx q[0];
rz(-0.91312983) q[0];
sx q[0];
rz(-3.012863) q[0];
x q[1];
rz(-1.2185368) q[2];
sx q[2];
rz(-2.2552236) q[2];
sx q[2];
rz(-0.5796488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8540878) q[1];
sx q[1];
rz(-0.64784986) q[1];
sx q[1];
rz(-1.2811529) q[1];
rz(-0.052107776) q[3];
sx q[3];
rz(-1.4192389) q[3];
sx q[3];
rz(2.1009171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.89221421) q[2];
sx q[2];
rz(-2.0912781) q[2];
sx q[2];
rz(-2.2155679) q[2];
rz(-1.2146436) q[3];
sx q[3];
rz(-0.49321431) q[3];
sx q[3];
rz(2.8652969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72325426) q[0];
sx q[0];
rz(-2.3746018) q[0];
sx q[0];
rz(-1.1085229) q[0];
rz(-2.8077937) q[1];
sx q[1];
rz(-1.8148986) q[1];
sx q[1];
rz(-2.0557859) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5949149) q[0];
sx q[0];
rz(-1.4913861) q[0];
sx q[0];
rz(0.20554849) q[0];
rz(-pi) q[1];
rz(-1.9778697) q[2];
sx q[2];
rz(-2.7086639) q[2];
sx q[2];
rz(-2.0100398) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.27583308) q[1];
sx q[1];
rz(-2.1001022) q[1];
sx q[1];
rz(-2.5342026) q[1];
rz(0.9341888) q[3];
sx q[3];
rz(-1.459834) q[3];
sx q[3];
rz(-2.8257089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6256025) q[2];
sx q[2];
rz(-2.7770999) q[2];
sx q[2];
rz(-1.1642574) q[2];
rz(-0.26816756) q[3];
sx q[3];
rz(-1.8987013) q[3];
sx q[3];
rz(0.75781649) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5371573) q[0];
sx q[0];
rz(-0.51979655) q[0];
sx q[0];
rz(1.9926158) q[0];
rz(-2.8474498) q[1];
sx q[1];
rz(-1.5584757) q[1];
sx q[1];
rz(-2.099096) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54759083) q[0];
sx q[0];
rz(-2.2896575) q[0];
sx q[0];
rz(-0.44525624) q[0];
rz(-0.78977079) q[2];
sx q[2];
rz(-1.4068479) q[2];
sx q[2];
rz(0.81610926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.52400541) q[1];
sx q[1];
rz(-1.397314) q[1];
sx q[1];
rz(-0.37131997) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5339666) q[3];
sx q[3];
rz(-0.66434089) q[3];
sx q[3];
rz(3.1337332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.59757549) q[2];
sx q[2];
rz(-0.67108265) q[2];
sx q[2];
rz(-1.3647122) q[2];
rz(-2.4890066) q[3];
sx q[3];
rz(-0.79129523) q[3];
sx q[3];
rz(-2.4932388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.515601) q[0];
sx q[0];
rz(-0.37726548) q[0];
sx q[0];
rz(-0.01734497) q[0];
rz(3.1248202) q[1];
sx q[1];
rz(-0.78712946) q[1];
sx q[1];
rz(2.9877072) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2543593) q[0];
sx q[0];
rz(-1.6463929) q[0];
sx q[0];
rz(-2.8235954) q[0];
rz(-2.1915386) q[2];
sx q[2];
rz(-0.95648208) q[2];
sx q[2];
rz(-1.4221734) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9964136) q[1];
sx q[1];
rz(-2.6657678) q[1];
sx q[1];
rz(2.4563172) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8762053) q[3];
sx q[3];
rz(-0.079566728) q[3];
sx q[3];
rz(0.1480535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.953557) q[2];
sx q[2];
rz(-2.3253658) q[2];
sx q[2];
rz(-0.44450644) q[2];
rz(-2.9514173) q[3];
sx q[3];
rz(-1.5580274) q[3];
sx q[3];
rz(-1.3727413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1403777) q[0];
sx q[0];
rz(-1.8920521) q[0];
sx q[0];
rz(2.0709399) q[0];
rz(3.095678) q[1];
sx q[1];
rz(-1.4713947) q[1];
sx q[1];
rz(0.79107034) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96595461) q[0];
sx q[0];
rz(-1.6583867) q[0];
sx q[0];
rz(1.9071294) q[0];
rz(-pi) q[1];
rz(2.7439762) q[2];
sx q[2];
rz(-1.8194345) q[2];
sx q[2];
rz(1.719081) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2998878) q[1];
sx q[1];
rz(-1.1763402) q[1];
sx q[1];
rz(1.1839371) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0088455) q[3];
sx q[3];
rz(-1.069396) q[3];
sx q[3];
rz(2.8341334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4967686) q[2];
sx q[2];
rz(-0.18459979) q[2];
sx q[2];
rz(1.4727288) q[2];
rz(1.6139) q[3];
sx q[3];
rz(-2.0852641) q[3];
sx q[3];
rz(1.9014026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.938217) q[0];
sx q[0];
rz(-1.6435517) q[0];
sx q[0];
rz(2.4284651) q[0];
rz(2.6279502) q[1];
sx q[1];
rz(-1.4331663) q[1];
sx q[1];
rz(-1.6930361) q[1];
rz(-1.3327053) q[2];
sx q[2];
rz(-1.0793964) q[2];
sx q[2];
rz(2.4894077) q[2];
rz(-0.41889965) q[3];
sx q[3];
rz(-2.2875026) q[3];
sx q[3];
rz(0.84411375) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

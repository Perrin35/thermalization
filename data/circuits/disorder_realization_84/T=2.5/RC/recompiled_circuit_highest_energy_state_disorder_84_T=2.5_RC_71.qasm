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
rz(-0.62491971) q[0];
sx q[0];
rz(-1.3345557) q[0];
sx q[0];
rz(-0.053243756) q[0];
rz(0.48625311) q[1];
sx q[1];
rz(6.1558131) q[1];
sx q[1];
rz(11.131412) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20986461) q[0];
sx q[0];
rz(-1.7057944) q[0];
sx q[0];
rz(-0.28688669) q[0];
x q[1];
rz(-0.26726079) q[2];
sx q[2];
rz(-1.4556985) q[2];
sx q[2];
rz(-2.5989344) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8813144) q[1];
sx q[1];
rz(-1.6242172) q[1];
sx q[1];
rz(0.39398663) q[1];
rz(-pi) q[2];
rz(2.6463935) q[3];
sx q[3];
rz(-2.1504068) q[3];
sx q[3];
rz(-0.51900348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3367553) q[2];
sx q[2];
rz(-1.3433604) q[2];
sx q[2];
rz(0.27377823) q[2];
rz(1.9233507) q[3];
sx q[3];
rz(-0.67353407) q[3];
sx q[3];
rz(-2.8555253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022920595) q[0];
sx q[0];
rz(-0.76347041) q[0];
sx q[0];
rz(0.77962312) q[0];
rz(0.66462213) q[1];
sx q[1];
rz(-1.9326991) q[1];
sx q[1];
rz(-1.8962616) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0459291) q[0];
sx q[0];
rz(-1.6324658) q[0];
sx q[0];
rz(-1.685919) q[0];
rz(-pi) q[1];
rz(0.092463569) q[2];
sx q[2];
rz(-0.86148724) q[2];
sx q[2];
rz(-2.6576328) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9960963) q[1];
sx q[1];
rz(-0.22095535) q[1];
sx q[1];
rz(-1.2902106) q[1];
rz(-pi) q[2];
rz(0.20540463) q[3];
sx q[3];
rz(-1.0213791) q[3];
sx q[3];
rz(-2.7618264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7094884) q[2];
sx q[2];
rz(-2.4435142) q[2];
sx q[2];
rz(3.091605) q[2];
rz(-1.6677808) q[3];
sx q[3];
rz(-1.7260909) q[3];
sx q[3];
rz(2.2059691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98899984) q[0];
sx q[0];
rz(-2.279156) q[0];
sx q[0];
rz(-1.9631901) q[0];
rz(1.1475457) q[1];
sx q[1];
rz(-0.18745628) q[1];
sx q[1];
rz(0.60290927) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5459203) q[0];
sx q[0];
rz(-1.4688558) q[0];
sx q[0];
rz(-1.7319182) q[0];
rz(0.96430594) q[2];
sx q[2];
rz(-0.78769257) q[2];
sx q[2];
rz(-0.76157969) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5192831) q[1];
sx q[1];
rz(-1.1596679) q[1];
sx q[1];
rz(-2.4101188) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52899482) q[3];
sx q[3];
rz(-2.1784665) q[3];
sx q[3];
rz(-1.120847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9511562) q[2];
sx q[2];
rz(-2.3556605) q[2];
sx q[2];
rz(-2.7583165) q[2];
rz(-1.8303309) q[3];
sx q[3];
rz(-1.0878599) q[3];
sx q[3];
rz(-0.12152984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7168147) q[0];
sx q[0];
rz(-0.99262339) q[0];
sx q[0];
rz(-0.92856652) q[0];
rz(0.83944744) q[1];
sx q[1];
rz(-1.8239559) q[1];
sx q[1];
rz(-2.3323257) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2977375) q[0];
sx q[0];
rz(-2.1240196) q[0];
sx q[0];
rz(-2.7064302) q[0];
x q[1];
rz(1.4014612) q[2];
sx q[2];
rz(-1.9953097) q[2];
sx q[2];
rz(0.2178387) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.41701554) q[1];
sx q[1];
rz(-2.8863686) q[1];
sx q[1];
rz(1.2947114) q[1];
rz(-pi) q[2];
rz(2.7326037) q[3];
sx q[3];
rz(-2.1738767) q[3];
sx q[3];
rz(2.0355253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3889918) q[2];
sx q[2];
rz(-1.5316803) q[2];
sx q[2];
rz(1.6947702) q[2];
rz(1.1270479) q[3];
sx q[3];
rz(-0.73486745) q[3];
sx q[3];
rz(-2.5126422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9160354) q[0];
sx q[0];
rz(-1.9434403) q[0];
sx q[0];
rz(1.4300562) q[0];
rz(0.23722181) q[1];
sx q[1];
rz(-2.1784541) q[1];
sx q[1];
rz(2.846834) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65564102) q[0];
sx q[0];
rz(-1.7136586) q[0];
sx q[0];
rz(0.44012196) q[0];
rz(0.90845256) q[2];
sx q[2];
rz(-1.2725432) q[2];
sx q[2];
rz(0.82240381) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.46550653) q[1];
sx q[1];
rz(-0.63356235) q[1];
sx q[1];
rz(-0.8224727) q[1];
x q[2];
rz(1.421991) q[3];
sx q[3];
rz(-2.1396881) q[3];
sx q[3];
rz(-2.1619145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1416867) q[2];
sx q[2];
rz(-0.20076951) q[2];
sx q[2];
rz(3.0086009) q[2];
rz(0.76987949) q[3];
sx q[3];
rz(-1.8498288) q[3];
sx q[3];
rz(0.31908116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.8517476) q[0];
sx q[0];
rz(-0.3388437) q[0];
sx q[0];
rz(-1.4122562) q[0];
rz(-3.0899561) q[1];
sx q[1];
rz(-2.5187571) q[1];
sx q[1];
rz(2.9488865) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7190921) q[0];
sx q[0];
rz(-1.0570881) q[0];
sx q[0];
rz(1.5490319) q[0];
x q[1];
rz(-2.4102845) q[2];
sx q[2];
rz(-2.0530862) q[2];
sx q[2];
rz(-0.77836365) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.35152516) q[1];
sx q[1];
rz(-1.1173101) q[1];
sx q[1];
rz(2.044911) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1901549) q[3];
sx q[3];
rz(-2.1771113) q[3];
sx q[3];
rz(2.9553138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.073033832) q[2];
sx q[2];
rz(-1.8517588) q[2];
sx q[2];
rz(1.3812836) q[2];
rz(-2.6265788) q[3];
sx q[3];
rz(-0.88740715) q[3];
sx q[3];
rz(-1.380434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.16159049) q[0];
sx q[0];
rz(-1.1450293) q[0];
sx q[0];
rz(2.7951151) q[0];
rz(-2.1222291) q[1];
sx q[1];
rz(-2.4788224) q[1];
sx q[1];
rz(-1.4541218) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2812735) q[0];
sx q[0];
rz(-2.3695282) q[0];
sx q[0];
rz(2.4738929) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90482651) q[2];
sx q[2];
rz(-0.77205333) q[2];
sx q[2];
rz(2.0029298) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.84216181) q[1];
sx q[1];
rz(-2.5379532) q[1];
sx q[1];
rz(-2.0718715) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.725127) q[3];
sx q[3];
rz(-1.6701506) q[3];
sx q[3];
rz(1.1442483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.60160294) q[2];
sx q[2];
rz(-1.8113965) q[2];
sx q[2];
rz(0.23923624) q[2];
rz(-2.7573977) q[3];
sx q[3];
rz(-2.4592063) q[3];
sx q[3];
rz(-2.189883) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9397028) q[0];
sx q[0];
rz(-1.6738142) q[0];
sx q[0];
rz(1.3772759) q[0];
rz(1.2205623) q[1];
sx q[1];
rz(-2.2790597) q[1];
sx q[1];
rz(-0.16622226) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21815366) q[0];
sx q[0];
rz(-1.3138608) q[0];
sx q[0];
rz(1.0983125) q[0];
rz(-pi) q[1];
rz(-1.8883838) q[2];
sx q[2];
rz(-2.2796405) q[2];
sx q[2];
rz(2.6203757) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6320914) q[1];
sx q[1];
rz(-0.21344122) q[1];
sx q[1];
rz(2.9963993) q[1];
rz(-pi) q[2];
rz(1.1508302) q[3];
sx q[3];
rz(-1.8059429) q[3];
sx q[3];
rz(-1.8311892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.347747) q[2];
sx q[2];
rz(-0.85024992) q[2];
sx q[2];
rz(-1.0469077) q[2];
rz(-1.8868123) q[3];
sx q[3];
rz(-1.0898277) q[3];
sx q[3];
rz(1.2966398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60085249) q[0];
sx q[0];
rz(-0.37373251) q[0];
sx q[0];
rz(-2.0284213) q[0];
rz(0.81740776) q[1];
sx q[1];
rz(-1.4459041) q[1];
sx q[1];
rz(-0.36849749) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8857074) q[0];
sx q[0];
rz(-1.8428546) q[0];
sx q[0];
rz(1.2124802) q[0];
rz(-pi) q[1];
rz(-1.4442367) q[2];
sx q[2];
rz(-2.309185) q[2];
sx q[2];
rz(1.0711993) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1296325) q[1];
sx q[1];
rz(-0.76597491) q[1];
sx q[1];
rz(-2.7523044) q[1];
rz(0.59898563) q[3];
sx q[3];
rz(-2.3040207) q[3];
sx q[3];
rz(-2.4323104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9937399) q[2];
sx q[2];
rz(-0.61823121) q[2];
sx q[2];
rz(-1.3573793) q[2];
rz(-2.3944858) q[3];
sx q[3];
rz(-0.96304572) q[3];
sx q[3];
rz(-2.4660306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6152182) q[0];
sx q[0];
rz(-2.6850061) q[0];
sx q[0];
rz(1.0427465) q[0];
rz(2.3710947) q[1];
sx q[1];
rz(-2.1310525) q[1];
sx q[1];
rz(-2.3811293) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7211518) q[0];
sx q[0];
rz(-1.8095224) q[0];
sx q[0];
rz(-1.2308285) q[0];
x q[1];
rz(-1.7508239) q[2];
sx q[2];
rz(-1.6086726) q[2];
sx q[2];
rz(-1.2779209) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.14899602) q[1];
sx q[1];
rz(-1.9933874) q[1];
sx q[1];
rz(-1.484297) q[1];
rz(-2.6728739) q[3];
sx q[3];
rz(-2.4303474) q[3];
sx q[3];
rz(-0.48392344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9383135) q[2];
sx q[2];
rz(-2.5310897) q[2];
sx q[2];
rz(-2.7321613) q[2];
rz(1.5583386) q[3];
sx q[3];
rz(-0.46054545) q[3];
sx q[3];
rz(-0.62622768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4938477) q[0];
sx q[0];
rz(-2.0997601) q[0];
sx q[0];
rz(0.44152015) q[0];
rz(-1.5589177) q[1];
sx q[1];
rz(-1.1987004) q[1];
sx q[1];
rz(-0.62335062) q[1];
rz(2.6263993) q[2];
sx q[2];
rz(-3.0233848) q[2];
sx q[2];
rz(1.530627) q[2];
rz(1.7033475) q[3];
sx q[3];
rz(-1.3302531) q[3];
sx q[3];
rz(2.1147685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

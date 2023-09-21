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
rz(-2.5897265) q[0];
sx q[0];
rz(-0.021835672) q[0];
rz(2.7472189) q[1];
sx q[1];
rz(-1.4596649) q[1];
sx q[1];
rz(2.9266761) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41479933) q[0];
sx q[0];
rz(-0.3814632) q[0];
sx q[0];
rz(0.8154072) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29834892) q[2];
sx q[2];
rz(-2.6539408) q[2];
sx q[2];
rz(-1.8698486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.51337459) q[1];
sx q[1];
rz(-0.77925032) q[1];
sx q[1];
rz(-0.90374225) q[1];
rz(-pi) q[2];
rz(-2.5932556) q[3];
sx q[3];
rz(-0.79474802) q[3];
sx q[3];
rz(-0.37219513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(0.56420502) q[2];
rz(1.365186) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(-1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0441701) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(0.92798293) q[0];
rz(-1.1652975) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(-0.89675084) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7433886) q[0];
sx q[0];
rz(-2.2687015) q[0];
sx q[0];
rz(-0.25701216) q[0];
x q[1];
rz(-1.4977658) q[2];
sx q[2];
rz(-1.739193) q[2];
sx q[2];
rz(1.8254691) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0577382) q[1];
sx q[1];
rz(-1.6946304) q[1];
sx q[1];
rz(-0.79353516) q[1];
rz(-pi) q[2];
rz(0.41665839) q[3];
sx q[3];
rz(-2.2778802) q[3];
sx q[3];
rz(1.7435031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.26560489) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(-1.0401475) q[2];
rz(1.4552207) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.74137694) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(-1.0282015) q[0];
rz(2.0630515) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(-2.7064586) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9520156) q[0];
sx q[0];
rz(-2.7727685) q[0];
sx q[0];
rz(-2.8394305) q[0];
rz(-pi) q[1];
rz(1.4235731) q[2];
sx q[2];
rz(-1.3230811) q[2];
sx q[2];
rz(1.148828) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1060461) q[1];
sx q[1];
rz(-1.4075081) q[1];
sx q[1];
rz(-1.0825023) q[1];
rz(2.7828091) q[3];
sx q[3];
rz(-2.3453418) q[3];
sx q[3];
rz(2.1735454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6083287) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(2.8386774) q[2];
rz(1.8164002) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(-0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9451697) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(2.2241425) q[0];
rz(-0.67287412) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(2.8767169) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47539513) q[0];
sx q[0];
rz(-2.4817433) q[0];
sx q[0];
rz(3.092479) q[0];
rz(0.27819602) q[2];
sx q[2];
rz(-1.2831266) q[2];
sx q[2];
rz(-2.1759335) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.781573) q[1];
sx q[1];
rz(-2.2847166) q[1];
sx q[1];
rz(-2.5931231) q[1];
x q[2];
rz(-0.78062765) q[3];
sx q[3];
rz(-0.98140162) q[3];
sx q[3];
rz(-1.9115703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.36310568) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(-1.5765566) q[2];
rz(-2.1145084) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89001369) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(-2.662861) q[0];
rz(1.0331253) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(0.95265257) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53341502) q[0];
sx q[0];
rz(-2.1153643) q[0];
sx q[0];
rz(1.7081225) q[0];
rz(-pi) q[1];
rz(-1.2383934) q[2];
sx q[2];
rz(-1.1754416) q[2];
sx q[2];
rz(-2.7796641) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6528656) q[1];
sx q[1];
rz(-1.1796724) q[1];
sx q[1];
rz(0.25686849) q[1];
x q[2];
rz(-2.6635366) q[3];
sx q[3];
rz(-1.607967) q[3];
sx q[3];
rz(2.2108848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4218563) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(-0.53058132) q[2];
rz(-1.7355708) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(-2.3099242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.574061) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(1.6249599) q[0];
rz(1.8364871) q[1];
sx q[1];
rz(-1.3508947) q[1];
sx q[1];
rz(0.17257246) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57219244) q[0];
sx q[0];
rz(-1.5571556) q[0];
sx q[0];
rz(1.9872679) q[0];
rz(-pi) q[1];
rz(1.1144981) q[2];
sx q[2];
rz(-1.3367532) q[2];
sx q[2];
rz(-0.32548387) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1291618) q[1];
sx q[1];
rz(-2.7134502) q[1];
sx q[1];
rz(1.8264324) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5112126) q[3];
sx q[3];
rz(-1.7582338) q[3];
sx q[3];
rz(2.558625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3036348) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(1.1266358) q[2];
rz(2.3593694) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(-1.3379898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28850266) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(2.4801168) q[0];
rz(0.96039564) q[1];
sx q[1];
rz(-1.7405225) q[1];
sx q[1];
rz(0.75659928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.341757) q[0];
sx q[0];
rz(-2.9633187) q[0];
sx q[0];
rz(1.0210277) q[0];
x q[1];
rz(-1.7724178) q[2];
sx q[2];
rz(-1.9747509) q[2];
sx q[2];
rz(0.23922761) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7797864) q[1];
sx q[1];
rz(-0.27946073) q[1];
sx q[1];
rz(-0.96868412) q[1];
x q[2];
rz(-2.9499801) q[3];
sx q[3];
rz(-1.629047) q[3];
sx q[3];
rz(-2.2239457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0044272) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(2.9471617) q[2];
rz(-2.2284609) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(0.98541361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(3.074926) q[0];
rz(-2.8170259) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(-2.1527122) q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(1.5867932) q[2];
sx q[2];
rz(-1.664186) q[2];
sx q[2];
rz(0.29856759) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5576396) q[1];
sx q[1];
rz(-1.9848616) q[1];
sx q[1];
rz(-0.072903452) q[1];
x q[2];
rz(0.51560651) q[3];
sx q[3];
rz(-2.0022087) q[3];
sx q[3];
rz(-0.069375667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4618335) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(-2.7339593) q[2];
rz(-0.76861012) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(-0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
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
rz(-1.2520049) q[1];
sx q[1];
rz(-2.2682155) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0489037) q[0];
sx q[0];
rz(-2.9498219) q[0];
sx q[0];
rz(1.2094686) q[0];
rz(1.7071502) q[2];
sx q[2];
rz(-1.4246203) q[2];
sx q[2];
rz(1.855195) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5492591) q[1];
sx q[1];
rz(-2.1570286) q[1];
sx q[1];
rz(1.1997644) q[1];
rz(-pi) q[2];
rz(2.14823) q[3];
sx q[3];
rz(-1.9462898) q[3];
sx q[3];
rz(1.7253699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1197027) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(2.4592887) q[2];
rz(2.7673289) q[3];
sx q[3];
rz(-1.572861) q[3];
sx q[3];
rz(2.9746829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
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
rz(0.69797126) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(-2.1886254) q[0];
rz(0.8264181) q[1];
sx q[1];
rz(-2.4024139) q[1];
sx q[1];
rz(-1.3964765) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29165927) q[0];
sx q[0];
rz(-2.8072661) q[0];
sx q[0];
rz(0.26419421) q[0];
rz(-pi) q[1];
rz(-1.868532) q[2];
sx q[2];
rz(-2.2072788) q[2];
sx q[2];
rz(-2.0993078) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6203306) q[1];
sx q[1];
rz(-1.8552823) q[1];
sx q[1];
rz(1.8950589) q[1];
rz(-pi) q[2];
rz(2.8793328) q[3];
sx q[3];
rz(-1.6256623) q[3];
sx q[3];
rz(-2.7454387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6754127) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(-2.9818025) q[2];
rz(-0.30188489) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(-2.8543499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044534279) q[0];
sx q[0];
rz(-0.67561588) q[0];
sx q[0];
rz(-1.5560879) q[0];
rz(0.13327577) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(0.940154) q[2];
sx q[2];
rz(-1.4174145) q[2];
sx q[2];
rz(-2.6819475) q[2];
rz(-2.6079569) q[3];
sx q[3];
rz(-1.613637) q[3];
sx q[3];
rz(-2.2231495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

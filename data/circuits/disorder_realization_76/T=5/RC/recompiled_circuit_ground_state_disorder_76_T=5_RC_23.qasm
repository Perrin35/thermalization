OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6948833) q[0];
sx q[0];
rz(-1.0816242) q[0];
sx q[0];
rz(2.4021436) q[0];
rz(2.3820355) q[1];
sx q[1];
rz(-1.3146725) q[1];
sx q[1];
rz(-1.7159599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6524871) q[0];
sx q[0];
rz(-1.28363) q[0];
sx q[0];
rz(1.0190359) q[0];
x q[1];
rz(-2.6447634) q[2];
sx q[2];
rz(-1.8112) q[2];
sx q[2];
rz(2.5508326) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6253852) q[1];
sx q[1];
rz(-0.16730669) q[1];
sx q[1];
rz(-1.0870666) q[1];
rz(-pi) q[2];
rz(2.2896122) q[3];
sx q[3];
rz(-1.5379099) q[3];
sx q[3];
rz(-2.6926958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2322959) q[2];
sx q[2];
rz(-2.2984419) q[2];
sx q[2];
rz(0.24093957) q[2];
rz(3.1124034) q[3];
sx q[3];
rz(-1.8029282) q[3];
sx q[3];
rz(-1.1791112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643243) q[0];
sx q[0];
rz(-1.203953) q[0];
sx q[0];
rz(2.6569195) q[0];
rz(-0.67972216) q[1];
sx q[1];
rz(-1.8599963) q[1];
sx q[1];
rz(1.9814804) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.401448) q[0];
sx q[0];
rz(-0.30456802) q[0];
sx q[0];
rz(-0.061621678) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4258575) q[2];
sx q[2];
rz(-1.01075) q[2];
sx q[2];
rz(-0.83362388) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1801123) q[1];
sx q[1];
rz(-2.4870076) q[1];
sx q[1];
rz(-0.41427362) q[1];
x q[2];
rz(2.6061222) q[3];
sx q[3];
rz(-2.5512085) q[3];
sx q[3];
rz(1.930742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29740563) q[2];
sx q[2];
rz(-2.83941) q[2];
sx q[2];
rz(2.6837132) q[2];
rz(-1.1229905) q[3];
sx q[3];
rz(-0.99959683) q[3];
sx q[3];
rz(-0.56639731) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8726525) q[0];
sx q[0];
rz(-0.70916969) q[0];
sx q[0];
rz(-3.1100682) q[0];
rz(-2.8541376) q[1];
sx q[1];
rz(-0.87702409) q[1];
sx q[1];
rz(-1.215975) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8376802) q[0];
sx q[0];
rz(-2.0995576) q[0];
sx q[0];
rz(2.3882475) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7599665) q[2];
sx q[2];
rz(-2.6256621) q[2];
sx q[2];
rz(2.6037773) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85092227) q[1];
sx q[1];
rz(-1.2953723) q[1];
sx q[1];
rz(1.2864134) q[1];
x q[2];
rz(0.6612571) q[3];
sx q[3];
rz(-2.1016069) q[3];
sx q[3];
rz(-2.9531933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.03881255) q[2];
sx q[2];
rz(-1.5422042) q[2];
sx q[2];
rz(0.83596027) q[2];
rz(-1.3577667) q[3];
sx q[3];
rz(-2.416555) q[3];
sx q[3];
rz(0.74762216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2826071) q[0];
sx q[0];
rz(-1.7361807) q[0];
sx q[0];
rz(-3.0614241) q[0];
rz(-2.0591586) q[1];
sx q[1];
rz(-2.9083462) q[1];
sx q[1];
rz(-1.3287883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69705582) q[0];
sx q[0];
rz(-1.2460684) q[0];
sx q[0];
rz(2.6022807) q[0];
x q[1];
rz(1.7432418) q[2];
sx q[2];
rz(-1.8891462) q[2];
sx q[2];
rz(-0.47296745) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7693682) q[1];
sx q[1];
rz(-1.9980248) q[1];
sx q[1];
rz(2.9823496) q[1];
rz(-1.7747709) q[3];
sx q[3];
rz(-1.2414724) q[3];
sx q[3];
rz(2.9993111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3749915) q[2];
sx q[2];
rz(-2.1583755) q[2];
sx q[2];
rz(0.90744606) q[2];
rz(-2.4713016) q[3];
sx q[3];
rz(-0.82796103) q[3];
sx q[3];
rz(-1.7214187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2403253) q[0];
sx q[0];
rz(-2.2690052) q[0];
sx q[0];
rz(-2.7101044) q[0];
rz(-2.0101428) q[1];
sx q[1];
rz(-1.1791469) q[1];
sx q[1];
rz(-0.44050899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90135306) q[0];
sx q[0];
rz(-2.6640737) q[0];
sx q[0];
rz(1.3130472) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9856057) q[2];
sx q[2];
rz(-1.3813263) q[2];
sx q[2];
rz(-1.0111601) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4104328) q[1];
sx q[1];
rz(-1.9010547) q[1];
sx q[1];
rz(1.6124658) q[1];
x q[2];
rz(-2.2045361) q[3];
sx q[3];
rz(-1.0770633) q[3];
sx q[3];
rz(-2.1449094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.12209192) q[2];
sx q[2];
rz(-2.2152405) q[2];
sx q[2];
rz(-2.7276373) q[2];
rz(1.7533938) q[3];
sx q[3];
rz(-1.5475169) q[3];
sx q[3];
rz(0.12510124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51467657) q[0];
sx q[0];
rz(-1.4056982) q[0];
sx q[0];
rz(0.93389121) q[0];
rz(0.67032188) q[1];
sx q[1];
rz(-1.8972081) q[1];
sx q[1];
rz(1.3353039) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7372197) q[0];
sx q[0];
rz(-1.6707194) q[0];
sx q[0];
rz(2.1735682) q[0];
rz(-pi) q[1];
rz(2.9839758) q[2];
sx q[2];
rz(-2.981957) q[2];
sx q[2];
rz(3.1363413) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.851593) q[1];
sx q[1];
rz(-1.4581175) q[1];
sx q[1];
rz(3.1188909) q[1];
rz(1.5794831) q[3];
sx q[3];
rz(-2.207805) q[3];
sx q[3];
rz(-2.6125233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2495217) q[2];
sx q[2];
rz(-0.27270174) q[2];
sx q[2];
rz(1.2707233) q[2];
rz(1.3696085) q[3];
sx q[3];
rz(-1.9000051) q[3];
sx q[3];
rz(-2.6197267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524566) q[0];
sx q[0];
rz(-0.58681762) q[0];
sx q[0];
rz(-0.15368803) q[0];
rz(-2.9391089) q[1];
sx q[1];
rz(-1.9053562) q[1];
sx q[1];
rz(0.94857803) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62371333) q[0];
sx q[0];
rz(-0.9615295) q[0];
sx q[0];
rz(1.3246956) q[0];
x q[1];
rz(-0.3237299) q[2];
sx q[2];
rz(-0.95697953) q[2];
sx q[2];
rz(2.987189) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2413509) q[1];
sx q[1];
rz(-2.7903922) q[1];
sx q[1];
rz(-2.3022149) q[1];
x q[2];
rz(0.92102401) q[3];
sx q[3];
rz(-1.463784) q[3];
sx q[3];
rz(-1.8821074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1367246) q[2];
sx q[2];
rz(-1.0978881) q[2];
sx q[2];
rz(3.1401805) q[2];
rz(3.0794365) q[3];
sx q[3];
rz(-1.7319738) q[3];
sx q[3];
rz(-2.1803161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.38178) q[0];
sx q[0];
rz(-0.75730046) q[0];
sx q[0];
rz(-1.2148452) q[0];
rz(2.3836721) q[1];
sx q[1];
rz(-0.52752033) q[1];
sx q[1];
rz(2.5835999) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76282952) q[0];
sx q[0];
rz(-2.7192594) q[0];
sx q[0];
rz(0.93393737) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4253769) q[2];
sx q[2];
rz(-2.1176257) q[2];
sx q[2];
rz(-2.0273923) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3000189) q[1];
sx q[1];
rz(-1.8487367) q[1];
sx q[1];
rz(1.0701847) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8479061) q[3];
sx q[3];
rz(-1.5449817) q[3];
sx q[3];
rz(1.9896979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5743635) q[2];
sx q[2];
rz(-1.1939253) q[2];
sx q[2];
rz(2.8723259) q[2];
rz(1.1632129) q[3];
sx q[3];
rz(-2.5389157) q[3];
sx q[3];
rz(0.24063024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.6729386) q[0];
sx q[0];
rz(-1.0373632) q[0];
sx q[0];
rz(-2.0682251) q[0];
rz(1.1277554) q[1];
sx q[1];
rz(-0.22722166) q[1];
sx q[1];
rz(-2.5849297) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7391602) q[0];
sx q[0];
rz(-0.35050979) q[0];
sx q[0];
rz(-1.1465766) q[0];
rz(-2.1416592) q[2];
sx q[2];
rz(-1.0733255) q[2];
sx q[2];
rz(1.2116878) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3240668) q[1];
sx q[1];
rz(-2.3285667) q[1];
sx q[1];
rz(0.50285411) q[1];
rz(2.9070517) q[3];
sx q[3];
rz(-1.6681156) q[3];
sx q[3];
rz(-0.095712599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90159455) q[2];
sx q[2];
rz(-1.0866714) q[2];
sx q[2];
rz(-1.979801) q[2];
rz(-2.6084172) q[3];
sx q[3];
rz(-2.6170001) q[3];
sx q[3];
rz(2.9840792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48851442) q[0];
sx q[0];
rz(-0.97452679) q[0];
sx q[0];
rz(2.5467806) q[0];
rz(-2.5626903) q[1];
sx q[1];
rz(-1.6659104) q[1];
sx q[1];
rz(-0.65931177) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4730277) q[0];
sx q[0];
rz(-2.1921232) q[0];
sx q[0];
rz(-1.4061808) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60229199) q[2];
sx q[2];
rz(-2.4458439) q[2];
sx q[2];
rz(0.93590036) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0486006) q[1];
sx q[1];
rz(-1.326418) q[1];
sx q[1];
rz(0.14631573) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80307428) q[3];
sx q[3];
rz(-2.7702906) q[3];
sx q[3];
rz(0.40483958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.79002964) q[2];
sx q[2];
rz(-1.1976676) q[2];
sx q[2];
rz(3.0885546) q[2];
rz(1.3430345) q[3];
sx q[3];
rz(-1.8648632) q[3];
sx q[3];
rz(3.067335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2863083) q[0];
sx q[0];
rz(-1.3636148) q[0];
sx q[0];
rz(-1.6028945) q[0];
rz(3.042649) q[1];
sx q[1];
rz(-1.7715441) q[1];
sx q[1];
rz(-3.0861707) q[1];
rz(-0.40545736) q[2];
sx q[2];
rz(-1.5787072) q[2];
sx q[2];
rz(-2.2070259) q[2];
rz(0.17114279) q[3];
sx q[3];
rz(-2.3928693) q[3];
sx q[3];
rz(0.00015043845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

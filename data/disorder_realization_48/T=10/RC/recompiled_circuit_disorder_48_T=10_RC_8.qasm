OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(-2.350783) q[0];
sx q[0];
rz(2.8074582) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(-0.9442803) q[1];
sx q[1];
rz(1.2184719) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1760362) q[0];
sx q[0];
rz(-0.82057014) q[0];
sx q[0];
rz(-1.6198938) q[0];
rz(-pi) q[1];
rz(-0.081131447) q[2];
sx q[2];
rz(-2.6754224) q[2];
sx q[2];
rz(2.1612338) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4301181) q[1];
sx q[1];
rz(-0.700044) q[1];
sx q[1];
rz(-2.3858566) q[1];
rz(-1.3487885) q[3];
sx q[3];
rz(-1.7553925) q[3];
sx q[3];
rz(2.8649462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.618764) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(-0.5775601) q[2];
rz(1.9918359) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(0.66453385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98786551) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(-2.7541449) q[0];
rz(-2.2024343) q[1];
sx q[1];
rz(-2.1444131) q[1];
sx q[1];
rz(-1.4025677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99968602) q[0];
sx q[0];
rz(-1.5748595) q[0];
sx q[0];
rz(-3.1121029) q[0];
rz(1.7895133) q[2];
sx q[2];
rz(-1.0220851) q[2];
sx q[2];
rz(2.1825841) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.84239292) q[1];
sx q[1];
rz(-1.0789011) q[1];
sx q[1];
rz(-1.1115587) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6658695) q[3];
sx q[3];
rz(-0.94082309) q[3];
sx q[3];
rz(-2.7271987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.42276057) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(0.31769162) q[2];
rz(-0.20673949) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(-0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3690255) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(-1.7279708) q[0];
rz(0.47779045) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(0.40107045) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7558407) q[0];
sx q[0];
rz(-2.7999561) q[0];
sx q[0];
rz(-1.2582448) q[0];
rz(2.4972649) q[2];
sx q[2];
rz(-1.7228848) q[2];
sx q[2];
rz(2.6459141) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8071825) q[1];
sx q[1];
rz(-1.5771598) q[1];
sx q[1];
rz(-2.4114354) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1776351) q[3];
sx q[3];
rz(-1.3893681) q[3];
sx q[3];
rz(1.3822615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5473189) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(-2.5857914) q[2];
rz(-0.9764955) q[3];
sx q[3];
rz(-0.54978168) q[3];
sx q[3];
rz(-2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(-1.4439616) q[0];
rz(1.5199039) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(0.25340432) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4356416) q[0];
sx q[0];
rz(-2.5905529) q[0];
sx q[0];
rz(-1.7680697) q[0];
rz(-pi) q[1];
rz(-1.0850111) q[2];
sx q[2];
rz(-1.5861142) q[2];
sx q[2];
rz(0.77335301) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.036382347) q[1];
sx q[1];
rz(-1.6670334) q[1];
sx q[1];
rz(2.5114245) q[1];
x q[2];
rz(-3.0074189) q[3];
sx q[3];
rz(-2.4878256) q[3];
sx q[3];
rz(1.8133481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0306586) q[2];
sx q[2];
rz(-1.7548283) q[2];
sx q[2];
rz(-2.3542662) q[2];
rz(-2.2287255) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(2.1319938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.426429) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(3.0786247) q[0];
rz(0.12403034) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(0.45809349) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8100909) q[0];
sx q[0];
rz(-1.1281663) q[0];
sx q[0];
rz(-1.5682194) q[0];
rz(-2.1797921) q[2];
sx q[2];
rz(-0.70448175) q[2];
sx q[2];
rz(1.5722164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0702857) q[1];
sx q[1];
rz(-2.3100393) q[1];
sx q[1];
rz(0.15858312) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0701418) q[3];
sx q[3];
rz(-1.6541012) q[3];
sx q[3];
rz(-1.7617776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1725585) q[2];
sx q[2];
rz(-2.2183552) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(-3.138792) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.3389448) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(2.2221185) q[0];
rz(0.062285034) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(-1.2671635) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6309109) q[0];
sx q[0];
rz(-1.1817389) q[0];
sx q[0];
rz(-2.3657777) q[0];
rz(1.426258) q[2];
sx q[2];
rz(-2.309531) q[2];
sx q[2];
rz(-1.2046255) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8900745) q[1];
sx q[1];
rz(-1.2444082) q[1];
sx q[1];
rz(0.20784394) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0503067) q[3];
sx q[3];
rz(-1.7877842) q[3];
sx q[3];
rz(0.12245164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.9617812) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(-0.58376694) q[2];
rz(2.4328655) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07847438) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(-2.069058) q[0];
rz(-0.5468927) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(1.1118719) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3277153) q[0];
sx q[0];
rz(-2.0534671) q[0];
sx q[0];
rz(0.10563235) q[0];
rz(0.56356168) q[2];
sx q[2];
rz(-0.60699082) q[2];
sx q[2];
rz(0.8286455) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3355616) q[1];
sx q[1];
rz(-2.8415488) q[1];
sx q[1];
rz(2.8801444) q[1];
x q[2];
rz(-0.92915793) q[3];
sx q[3];
rz(-0.60141364) q[3];
sx q[3];
rz(-1.8638368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.83773461) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(1.1676577) q[2];
rz(1.5363103) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(-1.4060098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.4090356) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(0.73079601) q[0];
rz(0.90019512) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(2.3866167) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6801493) q[0];
sx q[0];
rz(-1.4247243) q[0];
sx q[0];
rz(-0.84772528) q[0];
x q[1];
rz(-2.412699) q[2];
sx q[2];
rz(-1.1155827) q[2];
sx q[2];
rz(2.560937) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29754408) q[1];
sx q[1];
rz(-2.5522759) q[1];
sx q[1];
rz(0.20357666) q[1];
rz(0.49236492) q[3];
sx q[3];
rz(-2.1015321) q[3];
sx q[3];
rz(-2.9363971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3770611) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(-2.5047452) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(-1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72717845) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(-2.0027347) q[0];
rz(-0.75421929) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(0.019502217) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50197983) q[0];
sx q[0];
rz(-1.6507971) q[0];
sx q[0];
rz(-1.7777068) q[0];
x q[1];
rz(-1.7712014) q[2];
sx q[2];
rz(-1.7192171) q[2];
sx q[2];
rz(0.60147775) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6448977) q[1];
sx q[1];
rz(-2.1870107) q[1];
sx q[1];
rz(-2.9195021) q[1];
rz(-pi) q[2];
rz(1.2583624) q[3];
sx q[3];
rz(-1.4133246) q[3];
sx q[3];
rz(1.0851932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0043682178) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(0.84890378) q[2];
rz(-0.38765872) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7983109) q[0];
sx q[0];
rz(-2.9738975) q[0];
sx q[0];
rz(-0.48450255) q[0];
rz(1.3867406) q[1];
sx q[1];
rz(-1.4258899) q[1];
sx q[1];
rz(-1.1482931) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1185547) q[0];
sx q[0];
rz(-1.430129) q[0];
sx q[0];
rz(-1.5492357) q[0];
rz(-1.3482773) q[2];
sx q[2];
rz(-0.95138022) q[2];
sx q[2];
rz(-0.36048181) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3155568) q[1];
sx q[1];
rz(-0.6912187) q[1];
sx q[1];
rz(-1.8408937) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34928068) q[3];
sx q[3];
rz(-1.619907) q[3];
sx q[3];
rz(1.1790566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2293573) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(-2.7764376) q[2];
rz(3.0129516) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(2.685759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.1098332) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(0.96314349) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(0.49114901) q[2];
sx q[2];
rz(-2.4158203) q[2];
sx q[2];
rz(-0.62266785) q[2];
rz(-3.1046042) q[3];
sx q[3];
rz(-1.0109517) q[3];
sx q[3];
rz(3.0299822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

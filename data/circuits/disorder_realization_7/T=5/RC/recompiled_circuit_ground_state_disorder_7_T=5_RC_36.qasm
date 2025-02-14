OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6586128) q[0];
sx q[0];
rz(-0.40402544) q[0];
sx q[0];
rz(-2.7513096) q[0];
rz(-0.033493869) q[1];
sx q[1];
rz(3.6727603) q[1];
sx q[1];
rz(9.6277278) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34863198) q[0];
sx q[0];
rz(-1.9237776) q[0];
sx q[0];
rz(-2.1499632) q[0];
rz(-pi) q[1];
rz(-1.6338946) q[2];
sx q[2];
rz(-0.089091688) q[2];
sx q[2];
rz(-2.2697946) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.4642849) q[1];
sx q[1];
rz(-2.5618636) q[1];
sx q[1];
rz(-0.34601684) q[1];
rz(-pi) q[2];
rz(0.84439028) q[3];
sx q[3];
rz(-1.5914306) q[3];
sx q[3];
rz(-0.72935361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.94886327) q[2];
sx q[2];
rz(-0.84550965) q[2];
sx q[2];
rz(1.3884937) q[2];
rz(0.071831547) q[3];
sx q[3];
rz(-0.51126945) q[3];
sx q[3];
rz(2.3141919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0440867) q[0];
sx q[0];
rz(-0.16479099) q[0];
sx q[0];
rz(-2.6440115) q[0];
rz(-1.7454106) q[1];
sx q[1];
rz(-1.0284245) q[1];
sx q[1];
rz(2.5713249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.790417) q[0];
sx q[0];
rz(-0.47016682) q[0];
sx q[0];
rz(-0.29542653) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79945081) q[2];
sx q[2];
rz(-2.852147) q[2];
sx q[2];
rz(-1.4674526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5757338) q[1];
sx q[1];
rz(-1.5354146) q[1];
sx q[1];
rz(-2.4206764) q[1];
x q[2];
rz(-2.2482613) q[3];
sx q[3];
rz(-2.3260197) q[3];
sx q[3];
rz(1.4295242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0294864) q[2];
sx q[2];
rz(-1.4542397) q[2];
sx q[2];
rz(-0.71375978) q[2];
rz(-1.7589689) q[3];
sx q[3];
rz(-2.6191923) q[3];
sx q[3];
rz(2.6944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56938982) q[0];
sx q[0];
rz(-2.1694006) q[0];
sx q[0];
rz(-2.8045281) q[0];
rz(-0.49346787) q[1];
sx q[1];
rz(-2.4512873) q[1];
sx q[1];
rz(2.2148671) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0686058) q[0];
sx q[0];
rz(-1.301501) q[0];
sx q[0];
rz(0.26871839) q[0];
rz(-0.83374597) q[2];
sx q[2];
rz(-1.0930702) q[2];
sx q[2];
rz(-1.1857978) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87624967) q[1];
sx q[1];
rz(-2.690965) q[1];
sx q[1];
rz(-1.1875115) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51982359) q[3];
sx q[3];
rz(-0.9614203) q[3];
sx q[3];
rz(-1.5070311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9530764) q[2];
sx q[2];
rz(-2.2761554) q[2];
sx q[2];
rz(2.4600929) q[2];
rz(1.513688) q[3];
sx q[3];
rz(-1.3368006) q[3];
sx q[3];
rz(-0.17076913) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3156768) q[0];
sx q[0];
rz(-3.0382394) q[0];
sx q[0];
rz(-1.0507677) q[0];
rz(-0.61327618) q[1];
sx q[1];
rz(-0.79137099) q[1];
sx q[1];
rz(1.223986) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67658778) q[0];
sx q[0];
rz(-0.95116827) q[0];
sx q[0];
rz(1.9787491) q[0];
rz(-pi) q[1];
rz(-0.87665571) q[2];
sx q[2];
rz(-2.3842065) q[2];
sx q[2];
rz(-0.11981431) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5526841) q[1];
sx q[1];
rz(-0.97637227) q[1];
sx q[1];
rz(1.5556704) q[1];
rz(-pi) q[2];
rz(-1.7324034) q[3];
sx q[3];
rz(-2.1568885) q[3];
sx q[3];
rz(1.3368034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0918538) q[2];
sx q[2];
rz(-0.78511304) q[2];
sx q[2];
rz(-0.66748691) q[2];
rz(-1.5751754) q[3];
sx q[3];
rz(-2.9508041) q[3];
sx q[3];
rz(3.0070087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5176373) q[0];
sx q[0];
rz(-1.0409545) q[0];
sx q[0];
rz(-2.5868296) q[0];
rz(1.293921) q[1];
sx q[1];
rz(-0.41176739) q[1];
sx q[1];
rz(0.88465869) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22065565) q[0];
sx q[0];
rz(-1.0405128) q[0];
sx q[0];
rz(0.028593731) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0366465) q[2];
sx q[2];
rz(-2.2171582) q[2];
sx q[2];
rz(2.4798648) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5003839) q[1];
sx q[1];
rz(-0.51189089) q[1];
sx q[1];
rz(1.1275395) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8507666) q[3];
sx q[3];
rz(-1.9970702) q[3];
sx q[3];
rz(-0.72457641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.027792949) q[2];
sx q[2];
rz(-0.5641368) q[2];
sx q[2];
rz(-1.0998868) q[2];
rz(2.9187628) q[3];
sx q[3];
rz(-1.9375216) q[3];
sx q[3];
rz(-0.13404624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1564388) q[0];
sx q[0];
rz(-2.8414861) q[0];
sx q[0];
rz(1.9578178) q[0];
rz(1.8249594) q[1];
sx q[1];
rz(-0.96375179) q[1];
sx q[1];
rz(-2.1060941) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0619297) q[0];
sx q[0];
rz(-1.890914) q[0];
sx q[0];
rz(-1.0589664) q[0];
rz(-pi) q[1];
rz(-2.7155128) q[2];
sx q[2];
rz(-2.1392864) q[2];
sx q[2];
rz(0.3884494) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70724132) q[1];
sx q[1];
rz(-1.3907258) q[1];
sx q[1];
rz(-1.7456013) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4984958) q[3];
sx q[3];
rz(-0.80567718) q[3];
sx q[3];
rz(-1.8206545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3731132) q[2];
sx q[2];
rz(-0.36102411) q[2];
sx q[2];
rz(-0.38144544) q[2];
rz(-1.2162195) q[3];
sx q[3];
rz(-0.65665025) q[3];
sx q[3];
rz(-0.60429627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65378791) q[0];
sx q[0];
rz(-0.31823802) q[0];
sx q[0];
rz(2.8741264) q[0];
rz(-1.5221315) q[1];
sx q[1];
rz(-0.61264241) q[1];
sx q[1];
rz(-2.6681275) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7601022) q[0];
sx q[0];
rz(-2.6250429) q[0];
sx q[0];
rz(2.3960956) q[0];
rz(0.34469338) q[2];
sx q[2];
rz(-0.50486165) q[2];
sx q[2];
rz(2.8125151) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.59478851) q[1];
sx q[1];
rz(-1.686953) q[1];
sx q[1];
rz(-0.22211566) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4074989) q[3];
sx q[3];
rz(-1.9400404) q[3];
sx q[3];
rz(2.6881517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9240616) q[2];
sx q[2];
rz(-1.6465829) q[2];
sx q[2];
rz(-1.3727429) q[2];
rz(0.30630201) q[3];
sx q[3];
rz(-1.0270216) q[3];
sx q[3];
rz(-2.8021804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31996763) q[0];
sx q[0];
rz(-2.8447633) q[0];
sx q[0];
rz(2.6676275) q[0];
rz(-2.3560246) q[1];
sx q[1];
rz(-0.57890099) q[1];
sx q[1];
rz(2.2483291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5441204) q[0];
sx q[0];
rz(-1.5878146) q[0];
sx q[0];
rz(-1.3025137) q[0];
x q[1];
rz(-0.39126663) q[2];
sx q[2];
rz(-0.80384582) q[2];
sx q[2];
rz(-2.2101457) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.230669) q[1];
sx q[1];
rz(-1.603447) q[1];
sx q[1];
rz(-0.20882512) q[1];
x q[2];
rz(1.8056554) q[3];
sx q[3];
rz(-2.7357172) q[3];
sx q[3];
rz(-2.0879951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6713509) q[2];
sx q[2];
rz(-1.6157776) q[2];
sx q[2];
rz(0.5980171) q[2];
rz(2.9371069) q[3];
sx q[3];
rz(-3.0308767) q[3];
sx q[3];
rz(-0.71271768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0996899) q[0];
sx q[0];
rz(-2.1429017) q[0];
sx q[0];
rz(2.4424851) q[0];
rz(0.39012575) q[1];
sx q[1];
rz(-0.68123078) q[1];
sx q[1];
rz(-0.9763388) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5792865) q[0];
sx q[0];
rz(-1.7111943) q[0];
sx q[0];
rz(-1.3112351) q[0];
rz(-2.7039912) q[2];
sx q[2];
rz(-1.8695306) q[2];
sx q[2];
rz(1.7127812) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7568552) q[1];
sx q[1];
rz(-0.51402521) q[1];
sx q[1];
rz(-1.7208748) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5849708) q[3];
sx q[3];
rz(-1.8020013) q[3];
sx q[3];
rz(1.4299666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1583027) q[2];
sx q[2];
rz(-2.172564) q[2];
sx q[2];
rz(0.14652531) q[2];
rz(-0.26081416) q[3];
sx q[3];
rz(-2.0050037) q[3];
sx q[3];
rz(-2.8543616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2448267) q[0];
sx q[0];
rz(-0.34938669) q[0];
sx q[0];
rz(2.8283327) q[0];
rz(2.3406155) q[1];
sx q[1];
rz(-1.5104834) q[1];
sx q[1];
rz(2.7371178) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.411878) q[0];
sx q[0];
rz(-2.9094271) q[0];
sx q[0];
rz(2.4947255) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5953535) q[2];
sx q[2];
rz(-1.7405542) q[2];
sx q[2];
rz(1.7129829) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8449895) q[1];
sx q[1];
rz(-0.93520404) q[1];
sx q[1];
rz(-0.92439009) q[1];
rz(-pi) q[2];
rz(1.9451109) q[3];
sx q[3];
rz(-2.7817543) q[3];
sx q[3];
rz(-1.3364178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6910088) q[2];
sx q[2];
rz(-2.6477224) q[2];
sx q[2];
rz(2.3502926) q[2];
rz(-2.6386236) q[3];
sx q[3];
rz(-1.0478323) q[3];
sx q[3];
rz(0.80016971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41225152) q[0];
sx q[0];
rz(-1.3164192) q[0];
sx q[0];
rz(-1.3715716) q[0];
rz(-0.62660632) q[1];
sx q[1];
rz(-1.4817487) q[1];
sx q[1];
rz(2.1201835) q[1];
rz(1.7328429) q[2];
sx q[2];
rz(-1.2037983) q[2];
sx q[2];
rz(-1.498602) q[2];
rz(0.23343352) q[3];
sx q[3];
rz(-1.558565) q[3];
sx q[3];
rz(-0.19812921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

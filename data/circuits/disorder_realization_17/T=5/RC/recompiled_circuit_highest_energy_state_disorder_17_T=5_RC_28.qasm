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
rz(2.0789335) q[0];
sx q[0];
rz(-2.3261429) q[0];
sx q[0];
rz(2.4315779) q[0];
rz(-0.31399909) q[1];
sx q[1];
rz(-0.93500885) q[1];
sx q[1];
rz(11.234565) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3478617) q[0];
sx q[0];
rz(-1.8457185) q[0];
sx q[0];
rz(1.8040276) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5153036) q[2];
sx q[2];
rz(-2.7262027) q[2];
sx q[2];
rz(-0.20520575) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1430681) q[1];
sx q[1];
rz(-1.3245163) q[1];
sx q[1];
rz(2.1889436) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7053733) q[3];
sx q[3];
rz(-2.1582104) q[3];
sx q[3];
rz(2.9105216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4189202) q[2];
sx q[2];
rz(-1.2166497) q[2];
sx q[2];
rz(2.326272) q[2];
rz(0.80960649) q[3];
sx q[3];
rz(-1.6428734) q[3];
sx q[3];
rz(-2.0577551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.0631436) q[0];
sx q[0];
rz(-1.0493295) q[0];
sx q[0];
rz(0.11058841) q[0];
rz(0.94509205) q[1];
sx q[1];
rz(-2.7022305) q[1];
sx q[1];
rz(-1.5623215) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3695345) q[0];
sx q[0];
rz(-0.86449776) q[0];
sx q[0];
rz(2.1824271) q[0];
x q[1];
rz(-2.0192102) q[2];
sx q[2];
rz(-1.1875936) q[2];
sx q[2];
rz(-2.2926083) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.64792127) q[1];
sx q[1];
rz(-2.646138) q[1];
sx q[1];
rz(-2.9205771) q[1];
x q[2];
rz(-2.0840629) q[3];
sx q[3];
rz(-2.3239414) q[3];
sx q[3];
rz(3.0000837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9940146) q[2];
sx q[2];
rz(-0.14573228) q[2];
sx q[2];
rz(-0.4501403) q[2];
rz(1.133793) q[3];
sx q[3];
rz(-1.4530051) q[3];
sx q[3];
rz(-0.97761893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5728773) q[0];
sx q[0];
rz(-0.94803888) q[0];
sx q[0];
rz(-0.28133389) q[0];
rz(1.3866792) q[1];
sx q[1];
rz(-1.4298871) q[1];
sx q[1];
rz(-2.5414355) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2413797) q[0];
sx q[0];
rz(-2.8474387) q[0];
sx q[0];
rz(0.16071975) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1776206) q[2];
sx q[2];
rz(-1.5247048) q[2];
sx q[2];
rz(-0.54296903) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9439796) q[1];
sx q[1];
rz(-2.463732) q[1];
sx q[1];
rz(-2.6070836) q[1];
x q[2];
rz(-1.7417817) q[3];
sx q[3];
rz(-2.1372843) q[3];
sx q[3];
rz(-1.8734219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0448138) q[2];
sx q[2];
rz(-2.009095) q[2];
sx q[2];
rz(2.0951994) q[2];
rz(2.7438296) q[3];
sx q[3];
rz(-0.64266509) q[3];
sx q[3];
rz(-3.0090295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1716877) q[0];
sx q[0];
rz(-0.02709087) q[0];
sx q[0];
rz(2.0890253) q[0];
rz(1.196208) q[1];
sx q[1];
rz(-1.2095249) q[1];
sx q[1];
rz(0.41950163) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7561121) q[0];
sx q[0];
rz(-2.1364373) q[0];
sx q[0];
rz(1.8844946) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93650903) q[2];
sx q[2];
rz(-1.7542931) q[2];
sx q[2];
rz(1.2681792) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4070523) q[1];
sx q[1];
rz(-2.190935) q[1];
sx q[1];
rz(-0.46874115) q[1];
rz(0.50856789) q[3];
sx q[3];
rz(-1.3150104) q[3];
sx q[3];
rz(-0.05318197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9516051) q[2];
sx q[2];
rz(-2.0230468) q[2];
sx q[2];
rz(1.0901964) q[2];
rz(-0.53127855) q[3];
sx q[3];
rz(-0.71610206) q[3];
sx q[3];
rz(-1.6147015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4215609) q[0];
sx q[0];
rz(-1.5694542) q[0];
sx q[0];
rz(-1.74362) q[0];
rz(3.0086503) q[1];
sx q[1];
rz(-1.3016737) q[1];
sx q[1];
rz(-0.001210777) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0069790445) q[0];
sx q[0];
rz(-1.9513357) q[0];
sx q[0];
rz(-1.5534205) q[0];
rz(-pi) q[1];
rz(-2.9820061) q[2];
sx q[2];
rz(-1.288646) q[2];
sx q[2];
rz(2.7196376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.397314) q[1];
sx q[1];
rz(-0.71956735) q[1];
sx q[1];
rz(-2.6707763) q[1];
rz(-pi) q[2];
rz(0.026592606) q[3];
sx q[3];
rz(-1.7204434) q[3];
sx q[3];
rz(-2.1872471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91263897) q[2];
sx q[2];
rz(-1.3372083) q[2];
sx q[2];
rz(0.21052989) q[2];
rz(-1.0362961) q[3];
sx q[3];
rz(-0.784289) q[3];
sx q[3];
rz(-1.9265296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9482816) q[0];
sx q[0];
rz(-1.9183777) q[0];
sx q[0];
rz(2.7358828) q[0];
rz(1.7871208) q[1];
sx q[1];
rz(-2.3190119) q[1];
sx q[1];
rz(-0.92232651) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94033891) q[0];
sx q[0];
rz(-2.8329599) q[0];
sx q[0];
rz(1.2135394) q[0];
x q[1];
rz(1.4402661) q[2];
sx q[2];
rz(-2.4707327) q[2];
sx q[2];
rz(-2.6595569) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.74444425) q[1];
sx q[1];
rz(-1.4566629) q[1];
sx q[1];
rz(2.1836917) q[1];
rz(-1.9590553) q[3];
sx q[3];
rz(-1.045533) q[3];
sx q[3];
rz(2.1059259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4313844) q[2];
sx q[2];
rz(-1.2146543) q[2];
sx q[2];
rz(2.7396438) q[2];
rz(-1.8534144) q[3];
sx q[3];
rz(-2.2737019) q[3];
sx q[3];
rz(1.8624381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61889082) q[0];
sx q[0];
rz(-2.0382477) q[0];
sx q[0];
rz(3.0772305) q[0];
rz(2.3035658) q[1];
sx q[1];
rz(-2.3152654) q[1];
sx q[1];
rz(1.3305957) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42444077) q[0];
sx q[0];
rz(-2.1409875) q[0];
sx q[0];
rz(1.7198404) q[0];
rz(-pi) q[1];
rz(2.5352201) q[2];
sx q[2];
rz(-1.9909711) q[2];
sx q[2];
rz(2.5647031) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6846398) q[1];
sx q[1];
rz(-1.7906396) q[1];
sx q[1];
rz(2.6067642) q[1];
rz(2.9693611) q[3];
sx q[3];
rz(-0.96559286) q[3];
sx q[3];
rz(1.8021405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7360721) q[2];
sx q[2];
rz(-1.5084927) q[2];
sx q[2];
rz(2.4580477) q[2];
rz(-2.4077967) q[3];
sx q[3];
rz(-1.7283231) q[3];
sx q[3];
rz(2.2939513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34506327) q[0];
sx q[0];
rz(-1.0836443) q[0];
sx q[0];
rz(1.9684568) q[0];
rz(-2.7660811) q[1];
sx q[1];
rz(-2.651732) q[1];
sx q[1];
rz(1.0367905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4597367) q[0];
sx q[0];
rz(-1.5899183) q[0];
sx q[0];
rz(3.1179423) q[0];
x q[1];
rz(1.8470959) q[2];
sx q[2];
rz(-1.7299998) q[2];
sx q[2];
rz(-2.4197848) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4040096) q[1];
sx q[1];
rz(-2.3669457) q[1];
sx q[1];
rz(2.9402016) q[1];
x q[2];
rz(0.49077197) q[3];
sx q[3];
rz(-2.0476855) q[3];
sx q[3];
rz(-2.8618438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.57572395) q[2];
sx q[2];
rz(-0.91369358) q[2];
sx q[2];
rz(0.73053989) q[2];
rz(0.20914397) q[3];
sx q[3];
rz(-2.6181965) q[3];
sx q[3];
rz(1.8714347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65649477) q[0];
sx q[0];
rz(-2.7182343) q[0];
sx q[0];
rz(-3.094161) q[0];
rz(-0.221953) q[1];
sx q[1];
rz(-1.0739948) q[1];
sx q[1];
rz(-2.2304631) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.297356) q[0];
sx q[0];
rz(-0.043023303) q[0];
sx q[0];
rz(2.4639315) q[0];
rz(-pi) q[1];
rz(-1.5510727) q[2];
sx q[2];
rz(-0.64177401) q[2];
sx q[2];
rz(-0.80427792) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.641727) q[1];
sx q[1];
rz(-1.2033101) q[1];
sx q[1];
rz(-1.5572474) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0044596) q[3];
sx q[3];
rz(-1.3977504) q[3];
sx q[3];
rz(-2.6030948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.86965108) q[2];
sx q[2];
rz(-0.38001529) q[2];
sx q[2];
rz(-0.16858777) q[2];
rz(-2.9224959) q[3];
sx q[3];
rz(-0.85182652) q[3];
sx q[3];
rz(-1.3271837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1472226) q[0];
sx q[0];
rz(-0.34039012) q[0];
sx q[0];
rz(0.45743531) q[0];
rz(1.2262729) q[1];
sx q[1];
rz(-2.8044082) q[1];
sx q[1];
rz(1.7785243) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6206034) q[0];
sx q[0];
rz(-0.61043569) q[0];
sx q[0];
rz(-3.0362398) q[0];
rz(-1.326637) q[2];
sx q[2];
rz(-2.1863822) q[2];
sx q[2];
rz(2.2509839) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0911965) q[1];
sx q[1];
rz(-1.1938057) q[1];
sx q[1];
rz(1.3891549) q[1];
rz(-0.8882723) q[3];
sx q[3];
rz(-1.5049045) q[3];
sx q[3];
rz(2.1988188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7676131) q[2];
sx q[2];
rz(-0.55459905) q[2];
sx q[2];
rz(-0.45200959) q[2];
rz(-0.40397817) q[3];
sx q[3];
rz(-1.890506) q[3];
sx q[3];
rz(2.0154791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44611888) q[0];
sx q[0];
rz(-1.5344545) q[0];
sx q[0];
rz(-2.9547966) q[0];
rz(2.6192464) q[1];
sx q[1];
rz(-0.53032395) q[1];
sx q[1];
rz(1.2293336) q[1];
rz(-1.8277373) q[2];
sx q[2];
rz(-1.2972472) q[2];
sx q[2];
rz(-2.2157584) q[2];
rz(0.26351874) q[3];
sx q[3];
rz(-1.4265887) q[3];
sx q[3];
rz(-1.9909457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

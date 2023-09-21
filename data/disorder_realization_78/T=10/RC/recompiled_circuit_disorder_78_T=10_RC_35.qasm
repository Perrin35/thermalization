OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5821563) q[0];
sx q[0];
rz(-0.59098935) q[0];
sx q[0];
rz(-2.5581869) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(2.2489927) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4757724) q[0];
sx q[0];
rz(-1.8679529) q[0];
sx q[0];
rz(2.1509403) q[0];
rz(-2.4721488) q[2];
sx q[2];
rz(-1.1602243) q[2];
sx q[2];
rz(0.66561156) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8436369) q[1];
sx q[1];
rz(-1.8814109) q[1];
sx q[1];
rz(-2.288504) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6139908) q[3];
sx q[3];
rz(-2.3213263) q[3];
sx q[3];
rz(0.84212069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2261752) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(1.1606476) q[2];
rz(-2.9246269) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(-2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0579257) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(2.5657186) q[0];
rz(1.8946164) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(-1.974568) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8562247) q[0];
sx q[0];
rz(-1.8282187) q[0];
sx q[0];
rz(-3.0811937) q[0];
rz(-pi) q[1];
rz(-1.0300693) q[2];
sx q[2];
rz(-1.3582555) q[2];
sx q[2];
rz(2.9589047) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0868623) q[1];
sx q[1];
rz(-2.3332101) q[1];
sx q[1];
rz(-3.0562835) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59252177) q[3];
sx q[3];
rz(-2.3791109) q[3];
sx q[3];
rz(3.1069063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1077659) q[2];
sx q[2];
rz(-1.1788538) q[2];
sx q[2];
rz(-2.9272184) q[2];
rz(3.068148) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.6630994) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(1.2269155) q[0];
rz(-2.7413209) q[1];
sx q[1];
rz(-1.8881256) q[1];
sx q[1];
rz(-1.0148369) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0600216) q[0];
sx q[0];
rz(-0.9284174) q[0];
sx q[0];
rz(1.5006256) q[0];
x q[1];
rz(2.8995908) q[2];
sx q[2];
rz(-1.9654462) q[2];
sx q[2];
rz(-1.637527) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.240757) q[1];
sx q[1];
rz(-1.490331) q[1];
sx q[1];
rz(2.1686173) q[1];
x q[2];
rz(0.26168163) q[3];
sx q[3];
rz(-0.90173429) q[3];
sx q[3];
rz(0.79479587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0959452) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(-0.8992368) q[2];
rz(-2.4441161) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(0.60788679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65524453) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(0.43193257) q[0];
rz(2.5090384) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(0.63582173) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0441372) q[0];
sx q[0];
rz(-0.53253981) q[0];
sx q[0];
rz(-1.9299279) q[0];
rz(-pi) q[1];
rz(2.2976774) q[2];
sx q[2];
rz(-1.4471495) q[2];
sx q[2];
rz(-2.4713949) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7446049) q[1];
sx q[1];
rz(-0.79345353) q[1];
sx q[1];
rz(2.8810487) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63316791) q[3];
sx q[3];
rz(-0.35053262) q[3];
sx q[3];
rz(-1.7115953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2137961) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(2.4528743) q[2];
rz(-2.8074746) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(0.14373246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14389811) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(1.0744263) q[0];
rz(-0.74514666) q[1];
sx q[1];
rz(-1.4829758) q[1];
sx q[1];
rz(0.27854663) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2715831) q[0];
sx q[0];
rz(-0.9013114) q[0];
sx q[0];
rz(0.55218009) q[0];
rz(0.61168806) q[2];
sx q[2];
rz(-1.6606765) q[2];
sx q[2];
rz(2.9181366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6812233) q[1];
sx q[1];
rz(-1.9378098) q[1];
sx q[1];
rz(2.2296434) q[1];
rz(-2.3260818) q[3];
sx q[3];
rz(-1.1482571) q[3];
sx q[3];
rz(2.5583207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4429861) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(-2.8473575) q[2];
rz(-0.081929835) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(-0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0329523) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(-0.0090573514) q[0];
rz(-0.63502216) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(0.10805282) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8979643) q[0];
sx q[0];
rz(-2.0192696) q[0];
sx q[0];
rz(-1.2178221) q[0];
rz(2.75974) q[2];
sx q[2];
rz(-1.6994886) q[2];
sx q[2];
rz(0.51634386) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.15247004) q[1];
sx q[1];
rz(-1.9580541) q[1];
sx q[1];
rz(2.6168004) q[1];
rz(-pi) q[2];
rz(2.8764358) q[3];
sx q[3];
rz(-1.0401871) q[3];
sx q[3];
rz(1.6924072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99284995) q[2];
sx q[2];
rz(-1.7603346) q[2];
sx q[2];
rz(-1.2711058) q[2];
rz(-0.078401119) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(2.0097849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909797) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(-2.4839731) q[0];
rz(-1.744386) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(-2.2479642) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4511787) q[0];
sx q[0];
rz(-0.72890857) q[0];
sx q[0];
rz(-0.50509392) q[0];
rz(2.0486084) q[2];
sx q[2];
rz(-1.3078948) q[2];
sx q[2];
rz(-0.36793567) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0221755) q[1];
sx q[1];
rz(-2.4319473) q[1];
sx q[1];
rz(0.31125734) q[1];
x q[2];
rz(1.4736389) q[3];
sx q[3];
rz(-0.81340862) q[3];
sx q[3];
rz(-2.017445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.39067337) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(-0.52948362) q[2];
rz(2.6654065) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(-2.7990394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(0.7628409) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(2.0513127) q[0];
rz(3.0293363) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(-1.1539248) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9268122) q[0];
sx q[0];
rz(-1.6048389) q[0];
sx q[0];
rz(2.4916617) q[0];
x q[1];
rz(2.9032193) q[2];
sx q[2];
rz(-0.56406883) q[2];
sx q[2];
rz(-0.43018815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9824144) q[1];
sx q[1];
rz(-2.3783861) q[1];
sx q[1];
rz(-0.089341954) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89973255) q[3];
sx q[3];
rz(-1.2601488) q[3];
sx q[3];
rz(-1.8048546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5899137) q[2];
sx q[2];
rz(-0.38107291) q[2];
sx q[2];
rz(0.38044688) q[2];
rz(-2.0137265) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(1.5765223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8001051) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(0.63968101) q[0];
rz(-1.2387964) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(1.9715086) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3568048) q[0];
sx q[0];
rz(-2.6492282) q[0];
sx q[0];
rz(-2.3193293) q[0];
rz(-2.8662888) q[2];
sx q[2];
rz(-2.346056) q[2];
sx q[2];
rz(0.77992935) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3635892) q[1];
sx q[1];
rz(-0.52007857) q[1];
sx q[1];
rz(2.7705454) q[1];
rz(-pi) q[2];
rz(-2.2421078) q[3];
sx q[3];
rz(-1.6791108) q[3];
sx q[3];
rz(1.2182106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3141979) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(-2.7588552) q[2];
rz(2.2132204) q[3];
sx q[3];
rz(-1.174077) q[3];
sx q[3];
rz(-0.66463566) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9160354) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(-0.25892192) q[0];
rz(2.4312773) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(2.6616667) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40892884) q[0];
sx q[0];
rz(-0.73584475) q[0];
sx q[0];
rz(2.0793414) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2357113) q[2];
sx q[2];
rz(-1.3833628) q[2];
sx q[2];
rz(0.71042176) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.568798) q[1];
sx q[1];
rz(-2.5017782) q[1];
sx q[1];
rz(0.6154284) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6142526) q[3];
sx q[3];
rz(-0.96685997) q[3];
sx q[3];
rz(-2.7690943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.84247983) q[2];
sx q[2];
rz(-2.2128236) q[2];
sx q[2];
rz(1.8245565) q[2];
rz(-1.8995829) q[3];
sx q[3];
rz(-2.1879991) q[3];
sx q[3];
rz(0.73808134) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823572) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(2.1622529) q[1];
sx q[1];
rz(-2.0420488) q[1];
sx q[1];
rz(2.2534823) q[1];
rz(0.32817763) q[2];
sx q[2];
rz(-0.50445088) q[2];
sx q[2];
rz(2.9403461) q[2];
rz(0.27771523) q[3];
sx q[3];
rz(-1.82901) q[3];
sx q[3];
rz(1.3808586) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

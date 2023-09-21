OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10387575) q[0];
sx q[0];
rz(-1.9394983) q[0];
sx q[0];
rz(-1.1480968) q[0];
rz(-1.8885053) q[1];
sx q[1];
rz(-0.94068599) q[1];
sx q[1];
rz(-1.3936477) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049126547) q[0];
sx q[0];
rz(-1.8804212) q[0];
sx q[0];
rz(1.5686839) q[0];
rz(-pi) q[1];
rz(0.68140985) q[2];
sx q[2];
rz(-2.133956) q[2];
sx q[2];
rz(-1.6207221) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.71478903) q[1];
sx q[1];
rz(-1.1183294) q[1];
sx q[1];
rz(-1.2697551) q[1];
rz(-2.1288793) q[3];
sx q[3];
rz(-1.8126243) q[3];
sx q[3];
rz(-1.946132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6538438) q[2];
sx q[2];
rz(-1.2922492) q[2];
sx q[2];
rz(-0.1208819) q[2];
rz(-2.9623048) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(-2.9860935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0497465) q[0];
sx q[0];
rz(-0.76773983) q[0];
sx q[0];
rz(3.0088186) q[0];
rz(-1.4615387) q[1];
sx q[1];
rz(-1.5613873) q[1];
sx q[1];
rz(-2.9002088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7624843) q[0];
sx q[0];
rz(-2.603841) q[0];
sx q[0];
rz(2.3178029) q[0];
rz(-pi) q[1];
rz(-2.3953305) q[2];
sx q[2];
rz(-2.4733739) q[2];
sx q[2];
rz(1.8909188) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6893377) q[1];
sx q[1];
rz(-0.66879767) q[1];
sx q[1];
rz(1.4237088) q[1];
rz(-pi) q[2];
rz(1.4199735) q[3];
sx q[3];
rz(-1.5449459) q[3];
sx q[3];
rz(-2.6005656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1277348) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(2.0347118) q[2];
rz(1.3876623) q[3];
sx q[3];
rz(-0.51968402) q[3];
sx q[3];
rz(-1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36104193) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(2.4011491) q[0];
rz(-2.6904147) q[1];
sx q[1];
rz(-1.2606882) q[1];
sx q[1];
rz(-2.0887451) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8448062) q[0];
sx q[0];
rz(-2.1150981) q[0];
sx q[0];
rz(-1.143572) q[0];
rz(1.2028678) q[2];
sx q[2];
rz(-0.70764467) q[2];
sx q[2];
rz(-0.55756535) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4098674) q[1];
sx q[1];
rz(-1.4508529) q[1];
sx q[1];
rz(2.9883595) q[1];
rz(1.5490407) q[3];
sx q[3];
rz(-2.7113911) q[3];
sx q[3];
rz(2.3019626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8911002) q[2];
sx q[2];
rz(-1.5414457) q[2];
sx q[2];
rz(-2.5946674) q[2];
rz(-0.28918239) q[3];
sx q[3];
rz(-2.6047891) q[3];
sx q[3];
rz(-0.012399013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.1716487) q[0];
sx q[0];
rz(-0.26370731) q[0];
sx q[0];
rz(-1.3522211) q[0];
rz(0.2098473) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(-2.9052177) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.198092) q[0];
sx q[0];
rz(-0.31790942) q[0];
sx q[0];
rz(2.148118) q[0];
x q[1];
rz(-2.0868446) q[2];
sx q[2];
rz(-1.2630672) q[2];
sx q[2];
rz(-2.9204521) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35509767) q[1];
sx q[1];
rz(-1.7384643) q[1];
sx q[1];
rz(-0.0066890072) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27487367) q[3];
sx q[3];
rz(-0.50516869) q[3];
sx q[3];
rz(1.2884017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2410879) q[2];
sx q[2];
rz(-0.97854486) q[2];
sx q[2];
rz(-2.5881055) q[2];
rz(2.2262946) q[3];
sx q[3];
rz(-1.2585636) q[3];
sx q[3];
rz(-2.4966911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47167641) q[0];
sx q[0];
rz(-1.8138509) q[0];
sx q[0];
rz(-2.4374403) q[0];
rz(1.0559121) q[1];
sx q[1];
rz(-2.6864955) q[1];
sx q[1];
rz(2.5456837) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6248557) q[0];
sx q[0];
rz(-1.7863818) q[0];
sx q[0];
rz(-3.0481824) q[0];
rz(2.4679568) q[2];
sx q[2];
rz(-1.4623702) q[2];
sx q[2];
rz(-0.078439586) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6223645) q[1];
sx q[1];
rz(-1.9740826) q[1];
sx q[1];
rz(1.3695903) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3866053) q[3];
sx q[3];
rz(-1.2053688) q[3];
sx q[3];
rz(-0.75152498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.492505) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(2.5904783) q[2];
rz(0.20714949) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9838487) q[0];
sx q[0];
rz(-0.85726964) q[0];
sx q[0];
rz(-2.1898848) q[0];
rz(2.6668008) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(-2.8463083) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8808522) q[0];
sx q[0];
rz(-1.4319001) q[0];
sx q[0];
rz(-1.9689346) q[0];
x q[1];
rz(-1.4163383) q[2];
sx q[2];
rz(-1.8784349) q[2];
sx q[2];
rz(-2.6442106) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8413137) q[1];
sx q[1];
rz(-1.0611371) q[1];
sx q[1];
rz(2.120244) q[1];
rz(-1.4548364) q[3];
sx q[3];
rz(-0.51350683) q[3];
sx q[3];
rz(2.6503369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.39020145) q[2];
sx q[2];
rz(-1.3662806) q[2];
sx q[2];
rz(-2.2078216) q[2];
rz(-1.4592524) q[3];
sx q[3];
rz(-1.3542342) q[3];
sx q[3];
rz(-1.4565844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075832531) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(-0.24205762) q[0];
rz(-0.66479713) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(2.738293) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6751854) q[0];
sx q[0];
rz(-1.9403606) q[0];
sx q[0];
rz(2.7778366) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7935739) q[2];
sx q[2];
rz(-2.4325779) q[2];
sx q[2];
rz(-1.6955171) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.51822) q[1];
sx q[1];
rz(-1.5726591) q[1];
sx q[1];
rz(-1.3606447) q[1];
x q[2];
rz(1.7607479) q[3];
sx q[3];
rz(-1.5459832) q[3];
sx q[3];
rz(-2.1638889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2287801) q[2];
sx q[2];
rz(-2.0704724) q[2];
sx q[2];
rz(-2.7071803) q[2];
rz(2.1408391) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(0.51030695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.0077165724) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(-0.49466053) q[0];
rz(1.6330632) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(2.5411434) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531882) q[0];
sx q[0];
rz(-1.1115523) q[0];
sx q[0];
rz(-1.123239) q[0];
x q[1];
rz(2.0958488) q[2];
sx q[2];
rz(-1.5480435) q[2];
sx q[2];
rz(2.3843228) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8200092) q[1];
sx q[1];
rz(-1.4516186) q[1];
sx q[1];
rz(-2.9585341) q[1];
rz(-2.8377257) q[3];
sx q[3];
rz(-1.7700723) q[3];
sx q[3];
rz(2.5602333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1373458) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(0.11631575) q[2];
rz(-0.4256734) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(-1*pi/12) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9882934) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(1.4784038) q[0];
rz(0.93961811) q[1];
sx q[1];
rz(-1.3213108) q[1];
sx q[1];
rz(-0.41752648) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5696213) q[0];
sx q[0];
rz(-2.6082391) q[0];
sx q[0];
rz(-1.6589952) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0350424) q[2];
sx q[2];
rz(-1.3082192) q[2];
sx q[2];
rz(-3.0343461) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9081887) q[1];
sx q[1];
rz(-1.3797803) q[1];
sx q[1];
rz(1.7407655) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9478491) q[3];
sx q[3];
rz(-1.994641) q[3];
sx q[3];
rz(-1.2731874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4801165) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(-2.647906) q[2];
rz(0.72475973) q[3];
sx q[3];
rz(-1.1586435) q[3];
sx q[3];
rz(0.31989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9676554) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(2.4560112) q[0];
rz(0.29742345) q[1];
sx q[1];
rz(-2.90459) q[1];
sx q[1];
rz(-1.1313653) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2162227) q[0];
sx q[0];
rz(-2.4664481) q[0];
sx q[0];
rz(-2.569909) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93090242) q[2];
sx q[2];
rz(-2.1702592) q[2];
sx q[2];
rz(2.714389) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.79216204) q[1];
sx q[1];
rz(-1.3467448) q[1];
sx q[1];
rz(-0.38778023) q[1];
rz(-pi) q[2];
rz(-2.4990436) q[3];
sx q[3];
rz(-2.4276639) q[3];
sx q[3];
rz(-0.16845265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.82548213) q[2];
sx q[2];
rz(-1.9448514) q[2];
sx q[2];
rz(2.4342009) q[2];
rz(-2.3317544) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(-2.9530318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512882) q[0];
sx q[0];
rz(-1.123883) q[0];
sx q[0];
rz(-0.71258769) q[0];
rz(-0.21223016) q[1];
sx q[1];
rz(-1.6925015) q[1];
sx q[1];
rz(-0.5136516) q[1];
rz(0.60317294) q[2];
sx q[2];
rz(-1.0839673) q[2];
sx q[2];
rz(2.0187335) q[2];
rz(-2.5645732) q[3];
sx q[3];
rz(-0.95303017) q[3];
sx q[3];
rz(-0.87458761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
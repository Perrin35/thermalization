OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.360541) q[0];
sx q[0];
rz(-0.6147576) q[0];
sx q[0];
rz(-1.0714666) q[0];
rz(-0.40845025) q[1];
sx q[1];
rz(4.0790494) q[1];
sx q[1];
rz(11.117878) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14872257) q[0];
sx q[0];
rz(-0.62294423) q[0];
sx q[0];
rz(0.8533661) q[0];
x q[1];
rz(1.229847) q[2];
sx q[2];
rz(-2.1181137) q[2];
sx q[2];
rz(-0.42423074) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9706124) q[1];
sx q[1];
rz(-2.9677848) q[1];
sx q[1];
rz(2.5998678) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2812649) q[3];
sx q[3];
rz(-1.2324047) q[3];
sx q[3];
rz(0.30358728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.19746) q[2];
sx q[2];
rz(-0.84386533) q[2];
sx q[2];
rz(0.79305631) q[2];
rz(-2.872725) q[3];
sx q[3];
rz(-1.8157418) q[3];
sx q[3];
rz(-0.9602921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.8307513) q[0];
sx q[0];
rz(-0.38250592) q[0];
sx q[0];
rz(-2.261396) q[0];
rz(2.510732) q[1];
sx q[1];
rz(-2.3289101) q[1];
sx q[1];
rz(2.0426483) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0202206) q[0];
sx q[0];
rz(-1.3929875) q[0];
sx q[0];
rz(-1.5759379) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3395202) q[2];
sx q[2];
rz(-0.63098365) q[2];
sx q[2];
rz(0.6501261) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0365043) q[1];
sx q[1];
rz(-1.1394115) q[1];
sx q[1];
rz(-0.17891592) q[1];
rz(-pi) q[2];
rz(0.93370943) q[3];
sx q[3];
rz(-1.1843345) q[3];
sx q[3];
rz(-0.704788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7546996) q[2];
sx q[2];
rz(-1.4852445) q[2];
sx q[2];
rz(-2.5999542) q[2];
rz(0.70332876) q[3];
sx q[3];
rz(-2.7024305) q[3];
sx q[3];
rz(0.37731236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2331053) q[0];
sx q[0];
rz(-2.2623514) q[0];
sx q[0];
rz(3.044686) q[0];
rz(-2.5540409) q[1];
sx q[1];
rz(-2.2511626) q[1];
sx q[1];
rz(-2.5403835) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9438697) q[0];
sx q[0];
rz(-1.7515148) q[0];
sx q[0];
rz(0.69036071) q[0];
x q[1];
rz(-0.16232441) q[2];
sx q[2];
rz(-1.5078203) q[2];
sx q[2];
rz(1.8520825) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2612348) q[1];
sx q[1];
rz(-2.6812892) q[1];
sx q[1];
rz(-3.1241547) q[1];
x q[2];
rz(0.91368586) q[3];
sx q[3];
rz(-1.8538215) q[3];
sx q[3];
rz(1.0237657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84843695) q[2];
sx q[2];
rz(-1.4508672) q[2];
sx q[2];
rz(2.4270774) q[2];
rz(1.6589288) q[3];
sx q[3];
rz(-2.8667993) q[3];
sx q[3];
rz(2.7627435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42245427) q[0];
sx q[0];
rz(-1.2353354) q[0];
sx q[0];
rz(1.5041014) q[0];
rz(1.0927041) q[1];
sx q[1];
rz(-1.2810818) q[1];
sx q[1];
rz(0.5853931) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2566168) q[0];
sx q[0];
rz(-2.4451849) q[0];
sx q[0];
rz(2.710706) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4665514) q[2];
sx q[2];
rz(-1.5603258) q[2];
sx q[2];
rz(0.99828966) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5820739) q[1];
sx q[1];
rz(-0.93739707) q[1];
sx q[1];
rz(0.60057171) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0799418) q[3];
sx q[3];
rz(-0.98831359) q[3];
sx q[3];
rz(2.7459308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66712159) q[2];
sx q[2];
rz(-2.4690364) q[2];
sx q[2];
rz(0.39620623) q[2];
rz(-2.951156) q[3];
sx q[3];
rz(-0.93947828) q[3];
sx q[3];
rz(0.52085352) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025295479) q[0];
sx q[0];
rz(-0.46023661) q[0];
sx q[0];
rz(0.39475557) q[0];
rz(-1.0074298) q[1];
sx q[1];
rz(-1.9837244) q[1];
sx q[1];
rz(-1.5347068) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7459864) q[0];
sx q[0];
rz(-1.2885546) q[0];
sx q[0];
rz(-2.372118) q[0];
x q[1];
rz(2.7274969) q[2];
sx q[2];
rz(-2.3813435) q[2];
sx q[2];
rz(2.9000226) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0959209) q[1];
sx q[1];
rz(-1.878698) q[1];
sx q[1];
rz(-0.20521693) q[1];
rz(2.3716912) q[3];
sx q[3];
rz(-2.5320029) q[3];
sx q[3];
rz(1.23884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4452303) q[2];
sx q[2];
rz(-2.4424489) q[2];
sx q[2];
rz(0.17808476) q[2];
rz(-2.1961424) q[3];
sx q[3];
rz(-2.8036696) q[3];
sx q[3];
rz(2.7054355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4431385) q[0];
sx q[0];
rz(-0.66265023) q[0];
sx q[0];
rz(2.9887001) q[0];
rz(2.1235509) q[1];
sx q[1];
rz(-2.1986304) q[1];
sx q[1];
rz(-0.61000383) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4831938) q[0];
sx q[0];
rz(-0.81419277) q[0];
sx q[0];
rz(1.7479595) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0081530054) q[2];
sx q[2];
rz(-0.9969396) q[2];
sx q[2];
rz(-1.8620086) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2392366) q[1];
sx q[1];
rz(-0.78241759) q[1];
sx q[1];
rz(-2.8066325) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69645564) q[3];
sx q[3];
rz(-3.0142733) q[3];
sx q[3];
rz(-0.49139532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4623922) q[2];
sx q[2];
rz(-1.8588763) q[2];
sx q[2];
rz(0.25897762) q[2];
rz(-0.14608832) q[3];
sx q[3];
rz(-0.77272213) q[3];
sx q[3];
rz(-2.5025388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4884969) q[0];
sx q[0];
rz(-0.86576068) q[0];
sx q[0];
rz(-0.75188941) q[0];
rz(-2.409626) q[1];
sx q[1];
rz(-1.3804133) q[1];
sx q[1];
rz(2.6123349) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89443356) q[0];
sx q[0];
rz(-2.7934847) q[0];
sx q[0];
rz(1.0791808) q[0];
rz(-3.0050982) q[2];
sx q[2];
rz(-2.0220827) q[2];
sx q[2];
rz(0.1642483) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6700892) q[1];
sx q[1];
rz(-1.5836599) q[1];
sx q[1];
rz(1.6520581) q[1];
x q[2];
rz(-0.31792171) q[3];
sx q[3];
rz(-2.0057851) q[3];
sx q[3];
rz(0.11162139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.25972128) q[2];
sx q[2];
rz(-2.112381) q[2];
sx q[2];
rz(0.79505801) q[2];
rz(-1.2231539) q[3];
sx q[3];
rz(-1.7984248) q[3];
sx q[3];
rz(2.3651626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.935598) q[0];
sx q[0];
rz(-1.5933651) q[0];
sx q[0];
rz(0.11496168) q[0];
rz(3.0523172) q[1];
sx q[1];
rz(-0.99901366) q[1];
sx q[1];
rz(-2.9866536) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6892825) q[0];
sx q[0];
rz(-1.4703106) q[0];
sx q[0];
rz(1.7013676) q[0];
rz(-0.55122113) q[2];
sx q[2];
rz(-2.1586426) q[2];
sx q[2];
rz(1.2767222) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8304883) q[1];
sx q[1];
rz(-1.957867) q[1];
sx q[1];
rz(1.2933267) q[1];
rz(-pi) q[2];
rz(-1.4890758) q[3];
sx q[3];
rz(-2.4148259) q[3];
sx q[3];
rz(-0.86870199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0838919) q[2];
sx q[2];
rz(-1.8996779) q[2];
sx q[2];
rz(-2.4532301) q[2];
rz(-0.56811959) q[3];
sx q[3];
rz(-2.1391684) q[3];
sx q[3];
rz(-0.68305558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.6883009) q[0];
sx q[0];
rz(-2.330307) q[0];
sx q[0];
rz(-1.9360833) q[0];
rz(2.0506809) q[1];
sx q[1];
rz(-2.5542407) q[1];
sx q[1];
rz(-1.6039414) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6610049) q[0];
sx q[0];
rz(-1.93297) q[0];
sx q[0];
rz(-2.037338) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5700266) q[2];
sx q[2];
rz(-0.56790295) q[2];
sx q[2];
rz(-0.22837328) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6867169) q[1];
sx q[1];
rz(-1.4595965) q[1];
sx q[1];
rz(-0.79373549) q[1];
rz(-pi) q[2];
rz(-1.0991715) q[3];
sx q[3];
rz(-1.0256919) q[3];
sx q[3];
rz(0.68600149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.18498147) q[2];
sx q[2];
rz(-1.1143755) q[2];
sx q[2];
rz(-1.1448917) q[2];
rz(1.4780809) q[3];
sx q[3];
rz(-2.3749115) q[3];
sx q[3];
rz(-2.0293106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0612653) q[0];
sx q[0];
rz(-2.5054131) q[0];
sx q[0];
rz(-1.435745) q[0];
rz(-1.1371293) q[1];
sx q[1];
rz(-2.4500193) q[1];
sx q[1];
rz(2.2045076) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6244753) q[0];
sx q[0];
rz(-1.269425) q[0];
sx q[0];
rz(1.4791489) q[0];
x q[1];
rz(-1.1446196) q[2];
sx q[2];
rz(-1.1010896) q[2];
sx q[2];
rz(-1.5091648) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.67928606) q[1];
sx q[1];
rz(-2.3704154) q[1];
sx q[1];
rz(0.036840082) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61599515) q[3];
sx q[3];
rz(-1.8503891) q[3];
sx q[3];
rz(3.0699025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.81637853) q[2];
sx q[2];
rz(-2.4983695) q[2];
sx q[2];
rz(-0.31036672) q[2];
rz(-2.5988633) q[3];
sx q[3];
rz(-2.2118745) q[3];
sx q[3];
rz(-0.31699666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.19001374) q[0];
sx q[0];
rz(-1.0931451) q[0];
sx q[0];
rz(-2.6059294) q[0];
rz(0.71470064) q[1];
sx q[1];
rz(-1.2477881) q[1];
sx q[1];
rz(1.4664149) q[1];
rz(-1.9541478) q[2];
sx q[2];
rz(-1.9786096) q[2];
sx q[2];
rz(-0.65828029) q[2];
rz(-1.5101931) q[3];
sx q[3];
rz(-2.2987859) q[3];
sx q[3];
rz(1.4684341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

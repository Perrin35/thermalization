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
rz(-0.37762168) q[0];
sx q[0];
rz(-2.7130337) q[0];
sx q[0];
rz(-0.23634401) q[0];
rz(-1.6726681) q[1];
sx q[1];
rz(-1.5863215) q[1];
sx q[1];
rz(-0.16170391) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9028714) q[0];
sx q[0];
rz(-1.6797025) q[0];
sx q[0];
rz(0.98734537) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6290083) q[2];
sx q[2];
rz(-1.3752642) q[2];
sx q[2];
rz(2.9511676) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0031944) q[1];
sx q[1];
rz(-1.5439543) q[1];
sx q[1];
rz(0.0069199847) q[1];
rz(-pi) q[2];
rz(2.7667505) q[3];
sx q[3];
rz(-1.3656022) q[3];
sx q[3];
rz(-2.7480882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8990495) q[2];
sx q[2];
rz(-0.87729064) q[2];
sx q[2];
rz(2.7522932) q[2];
rz(0.23046514) q[3];
sx q[3];
rz(-3.1233628) q[3];
sx q[3];
rz(0.61483312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5742351) q[0];
sx q[0];
rz(-2.1901972) q[0];
sx q[0];
rz(-1.6472598) q[0];
rz(1.5556473) q[1];
sx q[1];
rz(-2.9289398) q[1];
sx q[1];
rz(-1.1344604) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0527549) q[0];
sx q[0];
rz(-2.2688534) q[0];
sx q[0];
rz(1.237117) q[0];
x q[1];
rz(2.0087162) q[2];
sx q[2];
rz(-0.88405245) q[2];
sx q[2];
rz(-1.1816292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5205343) q[1];
sx q[1];
rz(-2.1466594) q[1];
sx q[1];
rz(-0.49382468) q[1];
x q[2];
rz(-0.5364356) q[3];
sx q[3];
rz(-0.89563771) q[3];
sx q[3];
rz(1.9534115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.050934164) q[2];
sx q[2];
rz(-2.2218573) q[2];
sx q[2];
rz(1.8402137) q[2];
rz(-2.0672412) q[3];
sx q[3];
rz(-0.31000546) q[3];
sx q[3];
rz(-1.5955135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12304561) q[0];
sx q[0];
rz(-0.30236852) q[0];
sx q[0];
rz(0.61755919) q[0];
rz(1.1005719) q[1];
sx q[1];
rz(-3.1220084) q[1];
sx q[1];
rz(0.40357959) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8698533) q[0];
sx q[0];
rz(-1.6855441) q[0];
sx q[0];
rz(0.010557584) q[0];
rz(-pi) q[1];
rz(-0.52578747) q[2];
sx q[2];
rz(-1.4830768) q[2];
sx q[2];
rz(1.45873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8051648) q[1];
sx q[1];
rz(-1.4138828) q[1];
sx q[1];
rz(3.0057231) q[1];
rz(-pi) q[2];
rz(-1.8125435) q[3];
sx q[3];
rz(-1.0261593) q[3];
sx q[3];
rz(0.83078362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6214211) q[2];
sx q[2];
rz(-1.4474063) q[2];
sx q[2];
rz(-2.3604895) q[2];
rz(0.33629867) q[3];
sx q[3];
rz(-1.7016442) q[3];
sx q[3];
rz(2.4380016) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4562562) q[0];
sx q[0];
rz(-2.6254613) q[0];
sx q[0];
rz(1.2180895) q[0];
rz(1.3457899) q[1];
sx q[1];
rz(-3.1333874) q[1];
sx q[1];
rz(-1.2340612) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5033213) q[0];
sx q[0];
rz(-1.3809805) q[0];
sx q[0];
rz(1.580761) q[0];
rz(-pi) q[1];
rz(-1.8629563) q[2];
sx q[2];
rz(-2.1641762) q[2];
sx q[2];
rz(1.1012384) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.300894) q[1];
sx q[1];
rz(-1.3188682) q[1];
sx q[1];
rz(-1.9416652) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59438057) q[3];
sx q[3];
rz(-3.1295589) q[3];
sx q[3];
rz(1.2846701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1951083) q[2];
sx q[2];
rz(-0.4230963) q[2];
sx q[2];
rz(-1.1484324) q[2];
rz(-2.3871683) q[3];
sx q[3];
rz(-1.9781338) q[3];
sx q[3];
rz(-0.26328009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36963439) q[0];
sx q[0];
rz(-0.088005528) q[0];
sx q[0];
rz(2.5391286) q[0];
rz(0.98214904) q[1];
sx q[1];
rz(-0.0063449675) q[1];
sx q[1];
rz(-0.51012653) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0786297) q[0];
sx q[0];
rz(-1.7562683) q[0];
sx q[0];
rz(-1.3736891) q[0];
x q[1];
rz(-0.087156239) q[2];
sx q[2];
rz(-1.4498324) q[2];
sx q[2];
rz(0.5145413) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.229631) q[1];
sx q[1];
rz(-2.103064) q[1];
sx q[1];
rz(2.0144573) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81437935) q[3];
sx q[3];
rz(-1.5703716) q[3];
sx q[3];
rz(-0.18691508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0977352) q[2];
sx q[2];
rz(-1.1172224) q[2];
sx q[2];
rz(2.1065693) q[2];
rz(-1.467661) q[3];
sx q[3];
rz(-0.36164713) q[3];
sx q[3];
rz(-0.58290946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3007091) q[0];
sx q[0];
rz(-0.98475921) q[0];
sx q[0];
rz(-2.050052) q[0];
rz(0.40065271) q[1];
sx q[1];
rz(-0.0023829208) q[1];
sx q[1];
rz(-0.56652743) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7687205) q[0];
sx q[0];
rz(-2.3864288) q[0];
sx q[0];
rz(-1.7674957) q[0];
rz(-pi) q[1];
rz(-2.285901) q[2];
sx q[2];
rz(-1.0101748) q[2];
sx q[2];
rz(0.66455807) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.57921806) q[1];
sx q[1];
rz(-1.7953582) q[1];
sx q[1];
rz(-0.57266219) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6831538) q[3];
sx q[3];
rz(-2.1405725) q[3];
sx q[3];
rz(0.36142413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1953676) q[2];
sx q[2];
rz(-2.3905498) q[2];
sx q[2];
rz(1.5283778) q[2];
rz(-2.1402806) q[3];
sx q[3];
rz(-0.87037218) q[3];
sx q[3];
rz(-0.24054578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85750759) q[0];
sx q[0];
rz(-2.3473098) q[0];
sx q[0];
rz(1.1450144) q[0];
rz(-1.6192294) q[1];
sx q[1];
rz(-0.018773627) q[1];
sx q[1];
rz(1.1633263) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0534819) q[0];
sx q[0];
rz(-1.6170752) q[0];
sx q[0];
rz(1.2812216) q[0];
rz(0.14919632) q[2];
sx q[2];
rz(-2.2393804) q[2];
sx q[2];
rz(0.72982349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7419297) q[1];
sx q[1];
rz(-0.96782875) q[1];
sx q[1];
rz(-0.16027995) q[1];
rz(-2.9807253) q[3];
sx q[3];
rz(-2.5013574) q[3];
sx q[3];
rz(1.36659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6254977) q[2];
sx q[2];
rz(-0.98573804) q[2];
sx q[2];
rz(0.37334785) q[2];
rz(-0.18800023) q[3];
sx q[3];
rz(-1.5550273) q[3];
sx q[3];
rz(1.9787623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9173376) q[0];
sx q[0];
rz(-2.6118216) q[0];
sx q[0];
rz(-0.24714558) q[0];
rz(2.9180134) q[1];
sx q[1];
rz(-0.0037184628) q[1];
sx q[1];
rz(-1.5162969) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.926439) q[0];
sx q[0];
rz(-1.6351388) q[0];
sx q[0];
rz(1.6872726) q[0];
rz(-pi) q[1];
rz(-0.4319576) q[2];
sx q[2];
rz(-2.326845) q[2];
sx q[2];
rz(-0.88569631) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5943971) q[1];
sx q[1];
rz(-1.5867579) q[1];
sx q[1];
rz(0.39401024) q[1];
rz(0.070930158) q[3];
sx q[3];
rz(-0.59535691) q[3];
sx q[3];
rz(1.2687253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.197) q[2];
sx q[2];
rz(-0.71114117) q[2];
sx q[2];
rz(1.5468583) q[2];
rz(-2.366015) q[3];
sx q[3];
rz(-1.8577134) q[3];
sx q[3];
rz(-1.4625134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.813886) q[0];
sx q[0];
rz(-2.1729108) q[0];
sx q[0];
rz(-2.0342597) q[0];
rz(-2.2164717) q[1];
sx q[1];
rz(-0.0019625891) q[1];
sx q[1];
rz(-2.3855239) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30413142) q[0];
sx q[0];
rz(-0.63429773) q[0];
sx q[0];
rz(-1.0360121) q[0];
rz(1.5599361) q[2];
sx q[2];
rz(-2.6107222) q[2];
sx q[2];
rz(0.52299352) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.392627) q[1];
sx q[1];
rz(-2.4973329) q[1];
sx q[1];
rz(-0.49796748) q[1];
x q[2];
rz(1.6923201) q[3];
sx q[3];
rz(-1.0022707) q[3];
sx q[3];
rz(-2.4820676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5209311) q[2];
sx q[2];
rz(-0.88520092) q[2];
sx q[2];
rz(1.0010285) q[2];
rz(-1.3641317) q[3];
sx q[3];
rz(-2.2083211) q[3];
sx q[3];
rz(-0.53396839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.119568) q[0];
sx q[0];
rz(-1.3725932) q[0];
sx q[0];
rz(-2.6764828) q[0];
rz(-1.8067092) q[1];
sx q[1];
rz(-0.3723793) q[1];
sx q[1];
rz(-1.575527) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4441285) q[0];
sx q[0];
rz(-0.23374548) q[0];
sx q[0];
rz(1.759619) q[0];
rz(-pi) q[1];
x q[1];
rz(2.100824) q[2];
sx q[2];
rz(-1.841507) q[2];
sx q[2];
rz(1.9418429) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3153957) q[1];
sx q[1];
rz(-0.0028571833) q[1];
sx q[1];
rz(-0.14551659) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7622558) q[3];
sx q[3];
rz(-1.6054244) q[3];
sx q[3];
rz(1.4148764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2472725) q[2];
sx q[2];
rz(-0.044581052) q[2];
sx q[2];
rz(2.0559922) q[2];
rz(1.2224489) q[3];
sx q[3];
rz(-0.51210755) q[3];
sx q[3];
rz(1.2781757) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4492252) q[0];
sx q[0];
rz(-1.7074371) q[0];
sx q[0];
rz(-1.3771124) q[0];
rz(-1.5832681) q[1];
sx q[1];
rz(-0.91455864) q[1];
sx q[1];
rz(0.22462489) q[1];
rz(3.0947826) q[2];
sx q[2];
rz(-3.0407314) q[2];
sx q[2];
rz(-2.8032816) q[2];
rz(-1.1034154) q[3];
sx q[3];
rz(-2.6513908) q[3];
sx q[3];
rz(-1.8309616) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

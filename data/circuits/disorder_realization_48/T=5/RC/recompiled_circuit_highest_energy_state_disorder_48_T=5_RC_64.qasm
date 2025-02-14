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
rz(0.91520619) q[0];
sx q[0];
rz(-0.45867607) q[0];
sx q[0];
rz(0.31804481) q[0];
rz(1.3660499) q[1];
sx q[1];
rz(-2.4429758) q[1];
sx q[1];
rz(0.06279343) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2554277) q[0];
sx q[0];
rz(-2.0768099) q[0];
sx q[0];
rz(0.54765986) q[0];
rz(-pi) q[1];
rz(-0.89102052) q[2];
sx q[2];
rz(-1.9203548) q[2];
sx q[2];
rz(-2.120979) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5972728) q[1];
sx q[1];
rz(-1.8493127) q[1];
sx q[1];
rz(2.0803948) q[1];
rz(-2.1476521) q[3];
sx q[3];
rz(-1.2341043) q[3];
sx q[3];
rz(-1.5722881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.723145) q[2];
sx q[2];
rz(-1.3157088) q[2];
sx q[2];
rz(1.0203699) q[2];
rz(2.4127035) q[3];
sx q[3];
rz(-0.43292361) q[3];
sx q[3];
rz(-1.6093904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6272524) q[0];
sx q[0];
rz(-1.1662551) q[0];
sx q[0];
rz(-3.1033206) q[0];
rz(-3.017784) q[1];
sx q[1];
rz(-2.6688711) q[1];
sx q[1];
rz(-1.570787) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7109315) q[0];
sx q[0];
rz(-1.5607906) q[0];
sx q[0];
rz(-1.5902993) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6699617) q[2];
sx q[2];
rz(-2.3395774) q[2];
sx q[2];
rz(-2.2834509) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99421895) q[1];
sx q[1];
rz(-3.1404853) q[1];
sx q[1];
rz(1.1785517) q[1];
rz(1.4498715) q[3];
sx q[3];
rz(-1.0578511) q[3];
sx q[3];
rz(-0.67739048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6191285) q[2];
sx q[2];
rz(-1.8121441) q[2];
sx q[2];
rz(2.5627356) q[2];
rz(1.045643) q[3];
sx q[3];
rz(-0.78543109) q[3];
sx q[3];
rz(0.75700179) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4724562) q[0];
sx q[0];
rz(-1.0772912) q[0];
sx q[0];
rz(2.6876167) q[0];
rz(-0.16481915) q[1];
sx q[1];
rz(-2.3655128) q[1];
sx q[1];
rz(1.6710612) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9605345) q[0];
sx q[0];
rz(-1.1832602) q[0];
sx q[0];
rz(-2.9563532) q[0];
rz(-pi) q[1];
rz(-1.0452095) q[2];
sx q[2];
rz(-1.8586577) q[2];
sx q[2];
rz(0.20341104) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.20023055) q[1];
sx q[1];
rz(-2.388607) q[1];
sx q[1];
rz(0.45238564) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0461058) q[3];
sx q[3];
rz(-0.98011049) q[3];
sx q[3];
rz(2.5928796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7620324) q[2];
sx q[2];
rz(-1.4780059) q[2];
sx q[2];
rz(-1.7851768) q[2];
rz(2.6421269) q[3];
sx q[3];
rz(-1.6127337) q[3];
sx q[3];
rz(1.0511901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22719638) q[0];
sx q[0];
rz(-0.90407404) q[0];
sx q[0];
rz(-2.0830578) q[0];
rz(1.9805485) q[1];
sx q[1];
rz(-2.5807022) q[1];
sx q[1];
rz(3.0243691) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0974221) q[0];
sx q[0];
rz(-1.5316085) q[0];
sx q[0];
rz(0.54255658) q[0];
rz(-pi) q[1];
rz(2.2160596) q[2];
sx q[2];
rz(-1.4129854) q[2];
sx q[2];
rz(-2.7881546) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99771755) q[1];
sx q[1];
rz(-1.2923601) q[1];
sx q[1];
rz(0.77634676) q[1];
rz(2.6113308) q[3];
sx q[3];
rz(-0.70928364) q[3];
sx q[3];
rz(-1.5269296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3566572) q[2];
sx q[2];
rz(-1.6931809) q[2];
sx q[2];
rz(-1.3680722) q[2];
rz(-1.8562227) q[3];
sx q[3];
rz(-2.0338438) q[3];
sx q[3];
rz(2.4135597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042498978) q[0];
sx q[0];
rz(-0.92785257) q[0];
sx q[0];
rz(-1.7294783) q[0];
rz(2.1389351) q[1];
sx q[1];
rz(-2.899677) q[1];
sx q[1];
rz(1.5826506) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6928685) q[0];
sx q[0];
rz(-0.72843116) q[0];
sx q[0];
rz(-0.65653519) q[0];
rz(-pi) q[1];
rz(-1.4876258) q[2];
sx q[2];
rz(-2.2874333) q[2];
sx q[2];
rz(2.4019474) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4422426) q[1];
sx q[1];
rz(-1.0205871) q[1];
sx q[1];
rz(0.027632874) q[1];
rz(-3.047278) q[3];
sx q[3];
rz(-1.8911165) q[3];
sx q[3];
rz(-1.8220779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58854181) q[2];
sx q[2];
rz(-1.7091227) q[2];
sx q[2];
rz(-2.7480965) q[2];
rz(-0.95997512) q[3];
sx q[3];
rz(-0.67592755) q[3];
sx q[3];
rz(-0.967832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8680962) q[0];
sx q[0];
rz(-3.0143026) q[0];
sx q[0];
rz(-0.18336503) q[0];
rz(0.49194899) q[1];
sx q[1];
rz(-2.349647) q[1];
sx q[1];
rz(-1.2194182) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3122092) q[0];
sx q[0];
rz(-1.6512733) q[0];
sx q[0];
rz(3.0521554) q[0];
rz(-pi) q[1];
rz(3.1111818) q[2];
sx q[2];
rz(-1.4923499) q[2];
sx q[2];
rz(-0.65785656) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9522493) q[1];
sx q[1];
rz(-0.5533411) q[1];
sx q[1];
rz(2.4207741) q[1];
rz(2.0424834) q[3];
sx q[3];
rz(-1.6460643) q[3];
sx q[3];
rz(-1.690563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1313974) q[2];
sx q[2];
rz(-1.7051899) q[2];
sx q[2];
rz(3.0592226) q[2];
rz(0.92665893) q[3];
sx q[3];
rz(-0.72946531) q[3];
sx q[3];
rz(1.6389219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6019186) q[0];
sx q[0];
rz(-2.813756) q[0];
sx q[0];
rz(2.4251921) q[0];
rz(2.7377103) q[1];
sx q[1];
rz(-1.5733746) q[1];
sx q[1];
rz(-1.4366038) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88750171) q[0];
sx q[0];
rz(-1.1932826) q[0];
sx q[0];
rz(-2.98682) q[0];
rz(-pi) q[1];
x q[1];
rz(2.613343) q[2];
sx q[2];
rz(-1.6852323) q[2];
sx q[2];
rz(-0.50811646) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.77715141) q[1];
sx q[1];
rz(-1.6702515) q[1];
sx q[1];
rz(2.742139) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3167123) q[3];
sx q[3];
rz(-0.91478225) q[3];
sx q[3];
rz(-3.0716621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7776103) q[2];
sx q[2];
rz(-2.1668375) q[2];
sx q[2];
rz(0.065356143) q[2];
rz(0.86723793) q[3];
sx q[3];
rz(-1.5012274) q[3];
sx q[3];
rz(1.3245827) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48677483) q[0];
sx q[0];
rz(-1.4210533) q[0];
sx q[0];
rz(-0.73100334) q[0];
rz(-2.3664318) q[1];
sx q[1];
rz(-2.8072) q[1];
sx q[1];
rz(-2.7269272) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1461244) q[0];
sx q[0];
rz(-1.0421317) q[0];
sx q[0];
rz(0.84923262) q[0];
rz(0.53040217) q[2];
sx q[2];
rz(-1.7549522) q[2];
sx q[2];
rz(3.0006486) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.4463628) q[1];
sx q[1];
rz(-1.8703504) q[1];
sx q[1];
rz(1.3079554) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13402588) q[3];
sx q[3];
rz(-1.408159) q[3];
sx q[3];
rz(-2.7664879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.23184648) q[2];
sx q[2];
rz(-2.6555588) q[2];
sx q[2];
rz(3.0492142) q[2];
rz(0.77629027) q[3];
sx q[3];
rz(-1.679136) q[3];
sx q[3];
rz(2.9379454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6542776) q[0];
sx q[0];
rz(-1.4947816) q[0];
sx q[0];
rz(-3.1238632) q[0];
rz(1.0461294) q[1];
sx q[1];
rz(-1.7283864) q[1];
sx q[1];
rz(-2.0808992) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8413139) q[0];
sx q[0];
rz(-1.5848241) q[0];
sx q[0];
rz(3.0728691) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9477884) q[2];
sx q[2];
rz(-2.5193938) q[2];
sx q[2];
rz(-2.5925328) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6010161) q[1];
sx q[1];
rz(-2.3257044) q[1];
sx q[1];
rz(-0.23758446) q[1];
rz(1.6298196) q[3];
sx q[3];
rz(-1.0477598) q[3];
sx q[3];
rz(-0.28429261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9866508) q[2];
sx q[2];
rz(-0.81601802) q[2];
sx q[2];
rz(2.4972534) q[2];
rz(-2.2996357) q[3];
sx q[3];
rz(-1.569845) q[3];
sx q[3];
rz(0.90698609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3802721) q[0];
sx q[0];
rz(-2.705882) q[0];
sx q[0];
rz(-0.45183387) q[0];
rz(-2.1322346) q[1];
sx q[1];
rz(-1.6296856) q[1];
sx q[1];
rz(-1.5712646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0015948) q[0];
sx q[0];
rz(-2.9712935) q[0];
sx q[0];
rz(3.1055121) q[0];
rz(-pi) q[1];
rz(-1.3469996) q[2];
sx q[2];
rz(-0.48529709) q[2];
sx q[2];
rz(-1.8388302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6856076) q[1];
sx q[1];
rz(-2.7311014) q[1];
sx q[1];
rz(1.2430771) q[1];
rz(-0.0088457942) q[3];
sx q[3];
rz(-2.4191751) q[3];
sx q[3];
rz(2.9904108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5181804) q[2];
sx q[2];
rz(-2.0698915) q[2];
sx q[2];
rz(-1.5709467) q[2];
rz(0.10562854) q[3];
sx q[3];
rz(-2.195334) q[3];
sx q[3];
rz(-2.6251729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8138206) q[0];
sx q[0];
rz(-1.0186503) q[0];
sx q[0];
rz(-3.0042197) q[0];
rz(0.2188006) q[1];
sx q[1];
rz(-1.4004424) q[1];
sx q[1];
rz(-0.91011824) q[1];
rz(-1.995718) q[2];
sx q[2];
rz(-0.47987249) q[2];
sx q[2];
rz(2.1940827) q[2];
rz(1.2050592) q[3];
sx q[3];
rz(-2.1954721) q[3];
sx q[3];
rz(1.1385067) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

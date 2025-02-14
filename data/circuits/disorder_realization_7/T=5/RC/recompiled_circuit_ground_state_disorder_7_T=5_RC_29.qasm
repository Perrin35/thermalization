OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.4829798) q[0];
sx q[0];
rz(3.5456181) q[0];
sx q[0];
rz(9.0344949) q[0];
rz(-0.033493869) q[1];
sx q[1];
rz(3.6727603) q[1];
sx q[1];
rz(9.6277278) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73589486) q[0];
sx q[0];
rz(-2.4740088) q[0];
sx q[0];
rz(2.1632458) q[0];
rz(3.13596) q[2];
sx q[2];
rz(-1.4818824) q[2];
sx q[2];
rz(2.2064457) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.4642849) q[1];
sx q[1];
rz(-2.5618636) q[1];
sx q[1];
rz(-2.7955758) q[1];
x q[2];
rz(-0.027598782) q[3];
sx q[3];
rz(-2.2970133) q[3];
sx q[3];
rz(2.3184794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.94886327) q[2];
sx q[2];
rz(-2.296083) q[2];
sx q[2];
rz(-1.753099) q[2];
rz(3.0697611) q[3];
sx q[3];
rz(-2.6303232) q[3];
sx q[3];
rz(2.3141919) q[3];
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
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0440867) q[0];
sx q[0];
rz(-0.16479099) q[0];
sx q[0];
rz(-2.6440115) q[0];
rz(1.3961821) q[1];
sx q[1];
rz(-2.1131682) q[1];
sx q[1];
rz(-2.5713249) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48455301) q[0];
sx q[0];
rz(-1.4385106) q[0];
sx q[0];
rz(0.45251493) q[0];
x q[1];
rz(-1.7811586) q[2];
sx q[2];
rz(-1.3704925) q[2];
sx q[2];
rz(-2.4948073) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.03601053) q[1];
sx q[1];
rz(-0.85042989) q[1];
sx q[1];
rz(-1.6178814) q[1];
rz(0.87941283) q[3];
sx q[3];
rz(-2.0447404) q[3];
sx q[3];
rz(-0.64521257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1121062) q[2];
sx q[2];
rz(-1.4542397) q[2];
sx q[2];
rz(2.4278329) q[2];
rz(-1.7589689) q[3];
sx q[3];
rz(-0.52240038) q[3];
sx q[3];
rz(0.4471603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5722028) q[0];
sx q[0];
rz(-0.97219205) q[0];
sx q[0];
rz(0.33706459) q[0];
rz(0.49346787) q[1];
sx q[1];
rz(-2.4512873) q[1];
sx q[1];
rz(0.92672551) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7125366) q[0];
sx q[0];
rz(-1.829603) q[0];
sx q[0];
rz(-1.8496129) q[0];
rz(-0.83374597) q[2];
sx q[2];
rz(-2.0485224) q[2];
sx q[2];
rz(1.1857978) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7952319) q[1];
sx q[1];
rz(-1.7344001) q[1];
sx q[1];
rz(-1.1490046) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1893004) q[3];
sx q[3];
rz(-2.3626948) q[3];
sx q[3];
rz(2.2924045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1885163) q[2];
sx q[2];
rz(-0.86543721) q[2];
sx q[2];
rz(-2.4600929) q[2];
rz(1.6279047) q[3];
sx q[3];
rz(-1.8047921) q[3];
sx q[3];
rz(2.9708235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3156768) q[0];
sx q[0];
rz(-3.0382394) q[0];
sx q[0];
rz(-1.0507677) q[0];
rz(2.5283165) q[1];
sx q[1];
rz(-0.79137099) q[1];
sx q[1];
rz(-1.9176066) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67658778) q[0];
sx q[0];
rz(-0.95116827) q[0];
sx q[0];
rz(-1.1628435) q[0];
rz(-0.54398016) q[2];
sx q[2];
rz(-1.0144941) q[2];
sx q[2];
rz(-0.73327366) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.026583662) q[1];
sx q[1];
rz(-1.5833276) q[1];
sx q[1];
rz(0.59447713) q[1];
x q[2];
rz(-0.59215109) q[3];
sx q[3];
rz(-1.4363406) q[3];
sx q[3];
rz(-0.32392247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0497389) q[2];
sx q[2];
rz(-2.3564796) q[2];
sx q[2];
rz(-2.4741057) q[2];
rz(-1.5751754) q[3];
sx q[3];
rz(-0.19078855) q[3];
sx q[3];
rz(-3.0070087) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5176373) q[0];
sx q[0];
rz(-2.1006382) q[0];
sx q[0];
rz(2.5868296) q[0];
rz(-1.293921) q[1];
sx q[1];
rz(-0.41176739) q[1];
sx q[1];
rz(-0.88465869) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3356757) q[0];
sx q[0];
rz(-1.5954622) q[0];
sx q[0];
rz(-1.0403344) q[0];
rz(-pi) q[1];
rz(-0.53699643) q[2];
sx q[2];
rz(-2.3648713) q[2];
sx q[2];
rz(-3.1077488) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.139872) q[1];
sx q[1];
rz(-2.0291693) q[1];
sx q[1];
rz(2.9051347) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8507666) q[3];
sx q[3];
rz(-1.9970702) q[3];
sx q[3];
rz(0.72457641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1137997) q[2];
sx q[2];
rz(-0.5641368) q[2];
sx q[2];
rz(1.0998868) q[2];
rz(-2.9187628) q[3];
sx q[3];
rz(-1.9375216) q[3];
sx q[3];
rz(0.13404624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9851538) q[0];
sx q[0];
rz(-2.8414861) q[0];
sx q[0];
rz(-1.1837748) q[0];
rz(1.8249594) q[1];
sx q[1];
rz(-2.1778409) q[1];
sx q[1];
rz(2.1060941) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6661006) q[0];
sx q[0];
rz(-1.0872835) q[0];
sx q[0];
rz(2.7782281) q[0];
x q[1];
rz(-2.7155128) q[2];
sx q[2];
rz(-2.1392864) q[2];
sx q[2];
rz(-2.7531433) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.071152657) q[1];
sx q[1];
rz(-2.8912918) q[1];
sx q[1];
rz(2.3790199) q[1];
x q[2];
rz(-0.76642613) q[3];
sx q[3];
rz(-1.5186678) q[3];
sx q[3];
rz(-2.9418569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7684795) q[2];
sx q[2];
rz(-0.36102411) q[2];
sx q[2];
rz(-2.7601472) q[2];
rz(-1.9253731) q[3];
sx q[3];
rz(-0.65665025) q[3];
sx q[3];
rz(-2.5372964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65378791) q[0];
sx q[0];
rz(-2.8233546) q[0];
sx q[0];
rz(2.8741264) q[0];
rz(1.5221315) q[1];
sx q[1];
rz(-0.61264241) q[1];
sx q[1];
rz(-0.47346514) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43384837) q[0];
sx q[0];
rz(-1.9421541) q[0];
sx q[0];
rz(1.2030364) q[0];
rz(2.7968993) q[2];
sx q[2];
rz(-2.636731) q[2];
sx q[2];
rz(2.8125151) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.639714) q[1];
sx q[1];
rz(-2.8913829) q[1];
sx q[1];
rz(-2.6544957) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7678185) q[3];
sx q[3];
rz(-1.7229986) q[3];
sx q[3];
rz(-1.0579601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21753103) q[2];
sx q[2];
rz(-1.6465829) q[2];
sx q[2];
rz(1.3727429) q[2];
rz(2.8352906) q[3];
sx q[3];
rz(-2.114571) q[3];
sx q[3];
rz(0.33941227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.821625) q[0];
sx q[0];
rz(-2.8447633) q[0];
sx q[0];
rz(2.6676275) q[0];
rz(0.78556806) q[1];
sx q[1];
rz(-2.5626917) q[1];
sx q[1];
rz(-2.2483291) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97800228) q[0];
sx q[0];
rz(-1.8390391) q[0];
sx q[0];
rz(0.017649529) q[0];
rz(-pi) q[1];
rz(-1.1940057) q[2];
sx q[2];
rz(-2.2991141) q[2];
sx q[2];
rz(1.6737936) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.230669) q[1];
sx q[1];
rz(-1.603447) q[1];
sx q[1];
rz(0.20882512) q[1];
rz(-pi) q[2];
rz(-1.3359372) q[3];
sx q[3];
rz(-0.40587546) q[3];
sx q[3];
rz(-1.0535976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.47024176) q[2];
sx q[2];
rz(-1.525815) q[2];
sx q[2];
rz(-2.5435756) q[2];
rz(0.20448576) q[3];
sx q[3];
rz(-3.0308767) q[3];
sx q[3];
rz(-2.428875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04190271) q[0];
sx q[0];
rz(-2.1429017) q[0];
sx q[0];
rz(0.69910753) q[0];
rz(2.7514669) q[1];
sx q[1];
rz(-0.68123078) q[1];
sx q[1];
rz(-2.1652538) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6481144) q[0];
sx q[0];
rz(-2.8472487) q[0];
sx q[0];
rz(-1.0674547) q[0];
rz(-pi) q[1];
rz(-0.43760145) q[2];
sx q[2];
rz(-1.8695306) q[2];
sx q[2];
rz(-1.7127812) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7568552) q[1];
sx q[1];
rz(-2.6275674) q[1];
sx q[1];
rz(-1.4207178) q[1];
rz(-pi) q[2];
rz(1.3003208) q[3];
sx q[3];
rz(-2.1109443) q[3];
sx q[3];
rz(0.00080527472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1583027) q[2];
sx q[2];
rz(-2.172564) q[2];
sx q[2];
rz(0.14652531) q[2];
rz(-0.26081416) q[3];
sx q[3];
rz(-1.1365889) q[3];
sx q[3];
rz(2.8543616) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89676595) q[0];
sx q[0];
rz(-2.792206) q[0];
sx q[0];
rz(2.8283327) q[0];
rz(-2.3406155) q[1];
sx q[1];
rz(-1.5104834) q[1];
sx q[1];
rz(-2.7371178) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0718719) q[0];
sx q[0];
rz(-1.3861462) q[0];
sx q[0];
rz(1.7123366) q[0];
rz(-pi) q[1];
rz(2.8228797) q[2];
sx q[2];
rz(-0.56944263) q[2];
sx q[2];
rz(-0.12886831) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.39242649) q[1];
sx q[1];
rz(-0.87320864) q[1];
sx q[1];
rz(2.4568899) q[1];
rz(-1.1964818) q[3];
sx q[3];
rz(-0.35983837) q[3];
sx q[3];
rz(1.3364178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.45058382) q[2];
sx q[2];
rz(-2.6477224) q[2];
sx q[2];
rz(-2.3502926) q[2];
rz(2.6386236) q[3];
sx q[3];
rz(-2.0937604) q[3];
sx q[3];
rz(-2.3414229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7293411) q[0];
sx q[0];
rz(-1.3164192) q[0];
sx q[0];
rz(-1.3715716) q[0];
rz(-2.5149863) q[1];
sx q[1];
rz(-1.659844) q[1];
sx q[1];
rz(-1.0214092) q[1];
rz(1.7328429) q[2];
sx q[2];
rz(-1.2037983) q[2];
sx q[2];
rz(-1.498602) q[2];
rz(1.5833686) q[3];
sx q[3];
rz(-1.8042121) q[3];
sx q[3];
rz(1.369759) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

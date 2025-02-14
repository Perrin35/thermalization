OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98439944) q[0];
sx q[0];
rz(-1.1097044) q[0];
sx q[0];
rz(1.00534) q[0];
rz(0.73468626) q[1];
sx q[1];
rz(-2.7856196) q[1];
sx q[1];
rz(2.2495143) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9933424) q[0];
sx q[0];
rz(-0.82776946) q[0];
sx q[0];
rz(-1.7048852) q[0];
rz(-pi) q[1];
rz(2.9929225) q[2];
sx q[2];
rz(-0.97236247) q[2];
sx q[2];
rz(-0.27575803) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2764276) q[1];
sx q[1];
rz(-0.49404198) q[1];
sx q[1];
rz(-2.9345291) q[1];
x q[2];
rz(-0.69697505) q[3];
sx q[3];
rz(-2.4066952) q[3];
sx q[3];
rz(-1.01222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4165108) q[2];
sx q[2];
rz(-1.9981013) q[2];
sx q[2];
rz(-1.4408646) q[2];
rz(0.95669389) q[3];
sx q[3];
rz(-2.0243093) q[3];
sx q[3];
rz(-2.8982437) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44593921) q[0];
sx q[0];
rz(-0.39458269) q[0];
sx q[0];
rz(1.12895) q[0];
rz(0.24540643) q[1];
sx q[1];
rz(-1.8281728) q[1];
sx q[1];
rz(-2.479877) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0267994) q[0];
sx q[0];
rz(-0.68959348) q[0];
sx q[0];
rz(2.6845776) q[0];
x q[1];
rz(0.51742571) q[2];
sx q[2];
rz(-0.4330857) q[2];
sx q[2];
rz(0.171207) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9671038) q[1];
sx q[1];
rz(-1.2994819) q[1];
sx q[1];
rz(-1.5378465) q[1];
x q[2];
rz(3.0445736) q[3];
sx q[3];
rz(-1.2811913) q[3];
sx q[3];
rz(-1.9594909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5976065) q[2];
sx q[2];
rz(-0.49763766) q[2];
sx q[2];
rz(-1.0428766) q[2];
rz(-1.4526224) q[3];
sx q[3];
rz(-0.85113227) q[3];
sx q[3];
rz(2.0378621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76688981) q[0];
sx q[0];
rz(-0.68793982) q[0];
sx q[0];
rz(1.3091298) q[0];
rz(-2.6922928) q[1];
sx q[1];
rz(-2.4520912) q[1];
sx q[1];
rz(-1.7151394) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9096722) q[0];
sx q[0];
rz(-1.8350775) q[0];
sx q[0];
rz(2.8876165) q[0];
x q[1];
rz(-3.1235891) q[2];
sx q[2];
rz(-1.1483542) q[2];
sx q[2];
rz(1.2382335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0589185) q[1];
sx q[1];
rz(-0.86784309) q[1];
sx q[1];
rz(-1.4692131) q[1];
x q[2];
rz(-1.8168114) q[3];
sx q[3];
rz(-1.6520471) q[3];
sx q[3];
rz(1.8859091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.36491498) q[2];
sx q[2];
rz(-1.1740351) q[2];
sx q[2];
rz(0.066369973) q[2];
rz(-0.55473173) q[3];
sx q[3];
rz(-1.4812255) q[3];
sx q[3];
rz(-1.5167282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1944815) q[0];
sx q[0];
rz(-2.8567061) q[0];
sx q[0];
rz(2.3696005) q[0];
rz(2.4783065) q[1];
sx q[1];
rz(-2.3978077) q[1];
sx q[1];
rz(-0.1276806) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6948815) q[0];
sx q[0];
rz(-0.90411964) q[0];
sx q[0];
rz(1.5770017) q[0];
rz(-1.9468609) q[2];
sx q[2];
rz(-2.0916307) q[2];
sx q[2];
rz(-1.8686881) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7230715) q[1];
sx q[1];
rz(-1.1643895) q[1];
sx q[1];
rz(1.0999098) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33416602) q[3];
sx q[3];
rz(-2.6892956) q[3];
sx q[3];
rz(-0.75899198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6105303) q[2];
sx q[2];
rz(-2.1805111) q[2];
sx q[2];
rz(0.5985716) q[2];
rz(0.5851723) q[3];
sx q[3];
rz(-1.7681311) q[3];
sx q[3];
rz(-0.055805512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35036206) q[0];
sx q[0];
rz(-1.6329916) q[0];
sx q[0];
rz(2.1685261) q[0];
rz(0.19634253) q[1];
sx q[1];
rz(-1.1390431) q[1];
sx q[1];
rz(-1.6729209) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7356883) q[0];
sx q[0];
rz(-0.87474147) q[0];
sx q[0];
rz(2.619834) q[0];
rz(-2.8873047) q[2];
sx q[2];
rz(-2.0627705) q[2];
sx q[2];
rz(-2.6300583) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.10188516) q[1];
sx q[1];
rz(-2.7927924) q[1];
sx q[1];
rz(-2.0964016) q[1];
rz(-pi) q[2];
rz(1.0319388) q[3];
sx q[3];
rz(-1.6183637) q[3];
sx q[3];
rz(2.4470234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4520182) q[2];
sx q[2];
rz(-1.3268027) q[2];
sx q[2];
rz(-2.9791974) q[2];
rz(1.1299805) q[3];
sx q[3];
rz(-0.88310784) q[3];
sx q[3];
rz(2.0067154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3910386) q[0];
sx q[0];
rz(-2.6226608) q[0];
sx q[0];
rz(0.050405141) q[0];
rz(0.16818908) q[1];
sx q[1];
rz(-1.5610118) q[1];
sx q[1];
rz(0.87356299) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6872183) q[0];
sx q[0];
rz(-2.4719387) q[0];
sx q[0];
rz(-2.7424916) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3375521) q[2];
sx q[2];
rz(-0.14546673) q[2];
sx q[2];
rz(1.625753) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.37994036) q[1];
sx q[1];
rz(-2.4177175) q[1];
sx q[1];
rz(-2.4340802) q[1];
rz(2.843552) q[3];
sx q[3];
rz(-1.3750769) q[3];
sx q[3];
rz(-2.6258371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1217338) q[2];
sx q[2];
rz(-2.9254318) q[2];
sx q[2];
rz(-2.7609694) q[2];
rz(-2.5753283) q[3];
sx q[3];
rz(-1.7411313) q[3];
sx q[3];
rz(2.3316021) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5009163) q[0];
sx q[0];
rz(-0.2114507) q[0];
sx q[0];
rz(-0.40263116) q[0];
rz(1.7598049) q[1];
sx q[1];
rz(-0.80843061) q[1];
sx q[1];
rz(-0.056338739) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1409956) q[0];
sx q[0];
rz(-0.5282481) q[0];
sx q[0];
rz(0.88185593) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1001141) q[2];
sx q[2];
rz(-1.5342481) q[2];
sx q[2];
rz(1.4135897) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90625396) q[1];
sx q[1];
rz(-1.5829931) q[1];
sx q[1];
rz(2.1349299) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0545441) q[3];
sx q[3];
rz(-0.96025459) q[3];
sx q[3];
rz(-2.1757389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9075883) q[2];
sx q[2];
rz(-1.3434709) q[2];
sx q[2];
rz(-0.50055707) q[2];
rz(0.060700011) q[3];
sx q[3];
rz(-1.7874291) q[3];
sx q[3];
rz(2.6589656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.56977) q[0];
sx q[0];
rz(-1.7800542) q[0];
sx q[0];
rz(-1.7806336) q[0];
rz(-0.743615) q[1];
sx q[1];
rz(-1.4185602) q[1];
sx q[1];
rz(-0.13882151) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69037914) q[0];
sx q[0];
rz(-1.0745966) q[0];
sx q[0];
rz(2.1540533) q[0];
x q[1];
rz(2.2108364) q[2];
sx q[2];
rz(-0.70022196) q[2];
sx q[2];
rz(-0.28565684) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.946859) q[1];
sx q[1];
rz(-1.3689965) q[1];
sx q[1];
rz(-1.9029593) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0431792) q[3];
sx q[3];
rz(-2.7076027) q[3];
sx q[3];
rz(-1.4066309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4678141) q[2];
sx q[2];
rz(-2.0970586) q[2];
sx q[2];
rz(0.71411258) q[2];
rz(0.95947391) q[3];
sx q[3];
rz(-1.6474612) q[3];
sx q[3];
rz(-2.1625471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7239083) q[0];
sx q[0];
rz(-0.67643106) q[0];
sx q[0];
rz(-0.12829256) q[0];
rz(2.7543606) q[1];
sx q[1];
rz(-0.59806824) q[1];
sx q[1];
rz(-1.3234352) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5538226) q[0];
sx q[0];
rz(-2.9284796) q[0];
sx q[0];
rz(2.0487259) q[0];
rz(-pi) q[1];
rz(-2.2585715) q[2];
sx q[2];
rz(-0.81420126) q[2];
sx q[2];
rz(-1.3629066) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2880931) q[1];
sx q[1];
rz(-0.14460957) q[1];
sx q[1];
rz(1.4734731) q[1];
rz(-pi) q[2];
rz(-2.9186958) q[3];
sx q[3];
rz(-2.4011297) q[3];
sx q[3];
rz(0.84860138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6748176) q[2];
sx q[2];
rz(-1.520949) q[2];
sx q[2];
rz(-2.9368994) q[2];
rz(-1.8681059) q[3];
sx q[3];
rz(-0.73015648) q[3];
sx q[3];
rz(-1.5654806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8650763) q[0];
sx q[0];
rz(-0.62224329) q[0];
sx q[0];
rz(-2.9283071) q[0];
rz(-0.22008303) q[1];
sx q[1];
rz(-1.8890231) q[1];
sx q[1];
rz(1.782104) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8934568) q[0];
sx q[0];
rz(-1.5611959) q[0];
sx q[0];
rz(-1.5867771) q[0];
rz(0.47410902) q[2];
sx q[2];
rz(-1.342799) q[2];
sx q[2];
rz(1.0158418) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.90725856) q[1];
sx q[1];
rz(-1.9149019) q[1];
sx q[1];
rz(-2.1109778) q[1];
rz(-pi) q[2];
rz(0.90057365) q[3];
sx q[3];
rz(-1.2069205) q[3];
sx q[3];
rz(0.28387851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6243817) q[2];
sx q[2];
rz(-1.5456079) q[2];
sx q[2];
rz(-2.8037996) q[2];
rz(1.4126011) q[3];
sx q[3];
rz(-1.9905041) q[3];
sx q[3];
rz(-0.82144773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.4008041) q[0];
sx q[0];
rz(-0.82397006) q[0];
sx q[0];
rz(-1.9116221) q[0];
rz(-2.7578655) q[1];
sx q[1];
rz(-1.4839254) q[1];
sx q[1];
rz(0.76212777) q[1];
rz(-2.7376851) q[2];
sx q[2];
rz(-1.2474893) q[2];
sx q[2];
rz(2.9328109) q[2];
rz(-2.8970412) q[3];
sx q[3];
rz(-2.0691732) q[3];
sx q[3];
rz(0.3680784) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

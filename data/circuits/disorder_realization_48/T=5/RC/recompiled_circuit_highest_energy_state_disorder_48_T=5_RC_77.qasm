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
rz(-2.2263865) q[0];
sx q[0];
rz(-2.6829166) q[0];
sx q[0];
rz(-0.31804481) q[0];
rz(-1.7755427) q[1];
sx q[1];
rz(-0.69861689) q[1];
sx q[1];
rz(3.0787992) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97205596) q[0];
sx q[0];
rz(-1.0979303) q[0];
sx q[0];
rz(2.1465143) q[0];
rz(0.43831011) q[2];
sx q[2];
rz(-0.93899124) q[2];
sx q[2];
rz(2.8614728) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5972728) q[1];
sx q[1];
rz(-1.2922799) q[1];
sx q[1];
rz(2.0803948) q[1];
rz(-2.1413689) q[3];
sx q[3];
rz(-0.65815845) q[3];
sx q[3];
rz(-0.47129813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.723145) q[2];
sx q[2];
rz(-1.3157088) q[2];
sx q[2];
rz(1.0203699) q[2];
rz(2.4127035) q[3];
sx q[3];
rz(-2.708669) q[3];
sx q[3];
rz(-1.5322022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6272524) q[0];
sx q[0];
rz(-1.9753375) q[0];
sx q[0];
rz(0.038272055) q[0];
rz(0.12380869) q[1];
sx q[1];
rz(-0.47272155) q[1];
sx q[1];
rz(1.570787) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1403303) q[0];
sx q[0];
rz(-1.5512943) q[0];
sx q[0];
rz(-0.01000761) q[0];
x q[1];
rz(-1.471631) q[2];
sx q[2];
rz(-0.80201521) q[2];
sx q[2];
rz(0.85814171) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5396186) q[1];
sx q[1];
rz(-1.5718196) q[1];
sx q[1];
rz(-0.00042331528) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6255076) q[3];
sx q[3];
rz(-1.6760964) q[3];
sx q[3];
rz(2.1886231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6191285) q[2];
sx q[2];
rz(-1.8121441) q[2];
sx q[2];
rz(-0.5788571) q[2];
rz(-1.045643) q[3];
sx q[3];
rz(-0.78543109) q[3];
sx q[3];
rz(2.3845909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66913644) q[0];
sx q[0];
rz(-2.0643015) q[0];
sx q[0];
rz(2.6876167) q[0];
rz(-2.9767735) q[1];
sx q[1];
rz(-2.3655128) q[1];
sx q[1];
rz(-1.6710612) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8225518) q[0];
sx q[0];
rz(-1.7421573) q[0];
sx q[0];
rz(1.1771855) q[0];
rz(-pi) q[1];
rz(-2.8118089) q[2];
sx q[2];
rz(-1.068914) q[2];
sx q[2];
rz(-1.2041853) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1118033) q[1];
sx q[1];
rz(-1.2672499) q[1];
sx q[1];
rz(0.70036594) q[1];
x q[2];
rz(-0.64162125) q[3];
sx q[3];
rz(-0.76867662) q[3];
sx q[3];
rz(-0.25594433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7620324) q[2];
sx q[2];
rz(-1.4780059) q[2];
sx q[2];
rz(1.7851768) q[2];
rz(-2.6421269) q[3];
sx q[3];
rz(-1.5288589) q[3];
sx q[3];
rz(1.0511901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22719638) q[0];
sx q[0];
rz(-2.2375186) q[0];
sx q[0];
rz(2.0830578) q[0];
rz(-1.1610441) q[1];
sx q[1];
rz(-0.5608905) q[1];
sx q[1];
rz(0.11722359) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0974221) q[0];
sx q[0];
rz(-1.5316085) q[0];
sx q[0];
rz(-2.5990361) q[0];
rz(-2.2160596) q[2];
sx q[2];
rz(-1.7286073) q[2];
sx q[2];
rz(-2.7881546) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.99771755) q[1];
sx q[1];
rz(-1.8492325) q[1];
sx q[1];
rz(0.77634676) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53026188) q[3];
sx q[3];
rz(-2.432309) q[3];
sx q[3];
rz(1.6146631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3566572) q[2];
sx q[2];
rz(-1.4484118) q[2];
sx q[2];
rz(1.7735205) q[2];
rz(-1.8562227) q[3];
sx q[3];
rz(-2.0338438) q[3];
sx q[3];
rz(-0.72803298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0990937) q[0];
sx q[0];
rz(-2.2137401) q[0];
sx q[0];
rz(-1.7294783) q[0];
rz(-2.1389351) q[1];
sx q[1];
rz(-0.24191562) q[1];
sx q[1];
rz(-1.5589421) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8914192) q[0];
sx q[0];
rz(-1.0153664) q[0];
sx q[0];
rz(-1.0721747) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6539668) q[2];
sx q[2];
rz(-2.2874333) q[2];
sx q[2];
rz(-0.73964529) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0274899) q[1];
sx q[1];
rz(-1.5472425) q[1];
sx q[1];
rz(2.1211758) q[1];
x q[2];
rz(1.847397) q[3];
sx q[3];
rz(-0.33345727) q[3];
sx q[3];
rz(1.6113623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.58854181) q[2];
sx q[2];
rz(-1.43247) q[2];
sx q[2];
rz(-2.7480965) q[2];
rz(-0.95997512) q[3];
sx q[3];
rz(-2.4656651) q[3];
sx q[3];
rz(0.967832) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2734964) q[0];
sx q[0];
rz(-0.12729004) q[0];
sx q[0];
rz(0.18336503) q[0];
rz(-0.49194899) q[1];
sx q[1];
rz(-2.349647) q[1];
sx q[1];
rz(-1.9221745) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3122092) q[0];
sx q[0];
rz(-1.6512733) q[0];
sx q[0];
rz(-0.089437251) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9398795) q[2];
sx q[2];
rz(-3.0574692) q[2];
sx q[2];
rz(0.28757986) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1893434) q[1];
sx q[1];
rz(-0.5533411) q[1];
sx q[1];
rz(0.72081852) q[1];
rz(-pi) q[2];
rz(1.4063356) q[3];
sx q[3];
rz(-2.6643848) q[3];
sx q[3];
rz(-2.8754614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1313974) q[2];
sx q[2];
rz(-1.7051899) q[2];
sx q[2];
rz(-0.08237002) q[2];
rz(-0.92665893) q[3];
sx q[3];
rz(-0.72946531) q[3];
sx q[3];
rz(-1.6389219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53967404) q[0];
sx q[0];
rz(-2.813756) q[0];
sx q[0];
rz(2.4251921) q[0];
rz(0.40388233) q[1];
sx q[1];
rz(-1.5733746) q[1];
sx q[1];
rz(1.4366038) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5157455) q[0];
sx q[0];
rz(-1.7145918) q[0];
sx q[0];
rz(1.189144) q[0];
x q[1];
rz(1.4384957) q[2];
sx q[2];
rz(-2.0952333) q[2];
sx q[2];
rz(2.14545) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0246504) q[1];
sx q[1];
rz(-2.7305909) q[1];
sx q[1];
rz(-0.25115369) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8248803) q[3];
sx q[3];
rz(-2.2268104) q[3];
sx q[3];
rz(-0.069930596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3639823) q[2];
sx q[2];
rz(-2.1668375) q[2];
sx q[2];
rz(0.065356143) q[2];
rz(2.2743547) q[3];
sx q[3];
rz(-1.5012274) q[3];
sx q[3];
rz(1.8170099) q[3];
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
rz(-0.48677483) q[0];
sx q[0];
rz(-1.7205394) q[0];
sx q[0];
rz(0.73100334) q[0];
rz(-0.77516088) q[1];
sx q[1];
rz(-2.8072) q[1];
sx q[1];
rz(2.7269272) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99546826) q[0];
sx q[0];
rz(-2.099461) q[0];
sx q[0];
rz(2.29236) q[0];
x q[1];
rz(1.3581267) q[2];
sx q[2];
rz(-1.050282) q[2];
sx q[2];
rz(-1.5368324) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6952299) q[1];
sx q[1];
rz(-1.2712423) q[1];
sx q[1];
rz(-1.8336373) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.25423) q[3];
sx q[3];
rz(-0.21036869) q[3];
sx q[3];
rz(-0.31926814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9097462) q[2];
sx q[2];
rz(-2.6555588) q[2];
sx q[2];
rz(-3.0492142) q[2];
rz(-2.3653024) q[3];
sx q[3];
rz(-1.4624566) q[3];
sx q[3];
rz(-2.9379454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6542776) q[0];
sx q[0];
rz(-1.646811) q[0];
sx q[0];
rz(0.017729433) q[0];
rz(-1.0461294) q[1];
sx q[1];
rz(-1.4132063) q[1];
sx q[1];
rz(-2.0808992) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3002787) q[0];
sx q[0];
rz(-1.5567686) q[0];
sx q[0];
rz(-3.0728691) q[0];
rz(0.25814105) q[2];
sx q[2];
rz(-2.1434868) q[2];
sx q[2];
rz(-1.0024459) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.13421397) q[1];
sx q[1];
rz(-1.7430647) q[1];
sx q[1];
rz(-0.80165792) q[1];
rz(-1.5117731) q[3];
sx q[3];
rz(-2.0938329) q[3];
sx q[3];
rz(-2.8573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1549418) q[2];
sx q[2];
rz(-0.81601802) q[2];
sx q[2];
rz(2.4972534) q[2];
rz(2.2996357) q[3];
sx q[3];
rz(-1.5717477) q[3];
sx q[3];
rz(-2.2346066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7613206) q[0];
sx q[0];
rz(-2.705882) q[0];
sx q[0];
rz(-2.6897588) q[0];
rz(-1.0093581) q[1];
sx q[1];
rz(-1.511907) q[1];
sx q[1];
rz(-1.5712646) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0382045) q[0];
sx q[0];
rz(-1.7409835) q[0];
sx q[0];
rz(1.5645932) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0457479) q[2];
sx q[2];
rz(-1.4670851) q[2];
sx q[2];
rz(2.6748859) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.81105197) q[1];
sx q[1];
rz(-1.9582385) q[1];
sx q[1];
rz(-0.13918332) q[1];
rz(-pi) q[2];
rz(-2.4191946) q[3];
sx q[3];
rz(-1.5766451) q[3];
sx q[3];
rz(1.7286144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62341225) q[2];
sx q[2];
rz(-2.0698915) q[2];
sx q[2];
rz(-1.5706459) q[2];
rz(3.0359641) q[3];
sx q[3];
rz(-2.195334) q[3];
sx q[3];
rz(-0.51641974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8138206) q[0];
sx q[0];
rz(-2.1229424) q[0];
sx q[0];
rz(0.13737296) q[0];
rz(0.2188006) q[1];
sx q[1];
rz(-1.4004424) q[1];
sx q[1];
rz(-0.91011824) q[1];
rz(-2.9302421) q[2];
sx q[2];
rz(-1.1366781) q[2];
sx q[2];
rz(-0.47581262) q[2];
rz(-2.6811213) q[3];
sx q[3];
rz(-2.4302767) q[3];
sx q[3];
rz(-2.5828491) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

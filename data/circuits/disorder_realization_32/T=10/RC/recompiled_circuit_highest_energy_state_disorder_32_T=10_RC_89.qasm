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
rz(-1.929317) q[0];
sx q[0];
rz(-1.1317929) q[0];
sx q[0];
rz(-1.3728859) q[0];
rz(2.0139439) q[1];
sx q[1];
rz(-1.7665266) q[1];
sx q[1];
rz(2.3553203) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92106956) q[0];
sx q[0];
rz(-1.2036909) q[0];
sx q[0];
rz(-0.21998946) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.588618) q[2];
sx q[2];
rz(-0.99054256) q[2];
sx q[2];
rz(-1.0667104) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.3226124) q[1];
sx q[1];
rz(-2.8507345) q[1];
sx q[1];
rz(-1.4664654) q[1];
x q[2];
rz(-0.40793307) q[3];
sx q[3];
rz(-2.130748) q[3];
sx q[3];
rz(1.6168646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0218411) q[2];
sx q[2];
rz(-1.7181052) q[2];
sx q[2];
rz(-1.8915668) q[2];
rz(-0.73578468) q[3];
sx q[3];
rz(-2.9671228) q[3];
sx q[3];
rz(-1.5031987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64330548) q[0];
sx q[0];
rz(-1.1730288) q[0];
sx q[0];
rz(3.0702718) q[0];
rz(-1.1530863) q[1];
sx q[1];
rz(-0.76140296) q[1];
sx q[1];
rz(0.52871314) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0604989) q[0];
sx q[0];
rz(-1.1046263) q[0];
sx q[0];
rz(0.3549628) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18888338) q[2];
sx q[2];
rz(-1.3391277) q[2];
sx q[2];
rz(-2.846039) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0866249) q[1];
sx q[1];
rz(-1.0180915) q[1];
sx q[1];
rz(-2.5235559) q[1];
rz(-pi) q[2];
rz(1.0924322) q[3];
sx q[3];
rz(-2.6809664) q[3];
sx q[3];
rz(0.13308976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66136393) q[2];
sx q[2];
rz(-0.39863786) q[2];
sx q[2];
rz(0.028623494) q[2];
rz(2.7875767) q[3];
sx q[3];
rz(-2.386644) q[3];
sx q[3];
rz(2.5825175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70327586) q[0];
sx q[0];
rz(-2.1385312) q[0];
sx q[0];
rz(-0.89135528) q[0];
rz(0.50093961) q[1];
sx q[1];
rz(-1.5653862) q[1];
sx q[1];
rz(3.0827789) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0552154) q[0];
sx q[0];
rz(-0.098345938) q[0];
sx q[0];
rz(-2.6340061) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1040014) q[2];
sx q[2];
rz(-0.57764556) q[2];
sx q[2];
rz(1.4493329) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8938287) q[1];
sx q[1];
rz(-2.2732452) q[1];
sx q[1];
rz(2.48163) q[1];
x q[2];
rz(-2.9887588) q[3];
sx q[3];
rz(-2.9590769) q[3];
sx q[3];
rz(-2.07978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7848876) q[2];
sx q[2];
rz(-2.5829743) q[2];
sx q[2];
rz(2.9050262) q[2];
rz(-2.2208354) q[3];
sx q[3];
rz(-1.0509138) q[3];
sx q[3];
rz(-1.4111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9193566) q[0];
sx q[0];
rz(-1.284143) q[0];
sx q[0];
rz(2.2663569) q[0];
rz(1.5454166) q[1];
sx q[1];
rz(-2.2794006) q[1];
sx q[1];
rz(-0.045305591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6117258) q[0];
sx q[0];
rz(-1.5054724) q[0];
sx q[0];
rz(1.6732814) q[0];
x q[1];
rz(-0.50778054) q[2];
sx q[2];
rz(-0.96755257) q[2];
sx q[2];
rz(-2.809462) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26805763) q[1];
sx q[1];
rz(-2.2378439) q[1];
sx q[1];
rz(0.11063834) q[1];
rz(-pi) q[2];
rz(2.9475582) q[3];
sx q[3];
rz(-0.99655071) q[3];
sx q[3];
rz(1.5057092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8594325) q[2];
sx q[2];
rz(-2.1789447) q[2];
sx q[2];
rz(-1.1939987) q[2];
rz(2.2512839) q[3];
sx q[3];
rz(-1.2229536) q[3];
sx q[3];
rz(2.1931026) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10876656) q[0];
sx q[0];
rz(-0.81619167) q[0];
sx q[0];
rz(-0.68242514) q[0];
rz(-1.8531063) q[1];
sx q[1];
rz(-1.5792184) q[1];
sx q[1];
rz(1.3312181) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85859834) q[0];
sx q[0];
rz(-2.4147075) q[0];
sx q[0];
rz(1.5461393) q[0];
rz(1.9648026) q[2];
sx q[2];
rz(-1.8117684) q[2];
sx q[2];
rz(-2.7888128) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8537123) q[1];
sx q[1];
rz(-1.2128218) q[1];
sx q[1];
rz(1.720721) q[1];
rz(-2.2636637) q[3];
sx q[3];
rz(-0.64833855) q[3];
sx q[3];
rz(0.058365783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0266626) q[2];
sx q[2];
rz(-2.5422091) q[2];
sx q[2];
rz(-0.66693532) q[2];
rz(0.55495787) q[3];
sx q[3];
rz(-0.68735492) q[3];
sx q[3];
rz(-2.8426898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32894593) q[0];
sx q[0];
rz(-1.8829367) q[0];
sx q[0];
rz(-2.4378648) q[0];
rz(0.16920371) q[1];
sx q[1];
rz(-1.5237944) q[1];
sx q[1];
rz(2.9827548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1554398) q[0];
sx q[0];
rz(-1.4503398) q[0];
sx q[0];
rz(-1.4275761) q[0];
rz(1.6936982) q[2];
sx q[2];
rz(-1.4180866) q[2];
sx q[2];
rz(-1.0433152) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4602214) q[1];
sx q[1];
rz(-1.7029783) q[1];
sx q[1];
rz(0.48355196) q[1];
rz(0.63780906) q[3];
sx q[3];
rz(-1.3862228) q[3];
sx q[3];
rz(-2.0710433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3682897) q[2];
sx q[2];
rz(-0.94179073) q[2];
sx q[2];
rz(0.85757315) q[2];
rz(-1.9722021) q[3];
sx q[3];
rz(-1.2811067) q[3];
sx q[3];
rz(0.20839553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9661949) q[0];
sx q[0];
rz(-0.76681391) q[0];
sx q[0];
rz(-2.4229557) q[0];
rz(-0.28193685) q[1];
sx q[1];
rz(-1.7666631) q[1];
sx q[1];
rz(2.4729572) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40067264) q[0];
sx q[0];
rz(-0.95466475) q[0];
sx q[0];
rz(2.7223865) q[0];
rz(-pi) q[1];
x q[1];
rz(1.093542) q[2];
sx q[2];
rz(-0.68533939) q[2];
sx q[2];
rz(-0.46721855) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9033176) q[1];
sx q[1];
rz(-1.897745) q[1];
sx q[1];
rz(-1.9592778) q[1];
rz(-pi) q[2];
rz(0.63124048) q[3];
sx q[3];
rz(-1.5641912) q[3];
sx q[3];
rz(0.94830482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4947027) q[2];
sx q[2];
rz(-1.9560445) q[2];
sx q[2];
rz(2.1051712) q[2];
rz(-0.63452619) q[3];
sx q[3];
rz(-2.1536638) q[3];
sx q[3];
rz(2.3488267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-0.46102872) q[0];
sx q[0];
rz(-0.11756086) q[0];
sx q[0];
rz(0.72599894) q[0];
rz(-1.9793824) q[1];
sx q[1];
rz(-1.1770959) q[1];
sx q[1];
rz(-1.1964218) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0045753) q[0];
sx q[0];
rz(-1.5423623) q[0];
sx q[0];
rz(-1.397055) q[0];
rz(-0.047646626) q[2];
sx q[2];
rz(-1.8777913) q[2];
sx q[2];
rz(-2.7757211) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70564044) q[1];
sx q[1];
rz(-2.0355823) q[1];
sx q[1];
rz(-0.56805261) q[1];
rz(-pi) q[2];
rz(-2.7128588) q[3];
sx q[3];
rz(-2.4247243) q[3];
sx q[3];
rz(2.4459239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90765816) q[2];
sx q[2];
rz(-1.5479167) q[2];
sx q[2];
rz(0.93351239) q[2];
rz(2.524611) q[3];
sx q[3];
rz(-2.139293) q[3];
sx q[3];
rz(2.2945819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9762978) q[0];
sx q[0];
rz(-3.0369861) q[0];
sx q[0];
rz(-3.0883375) q[0];
rz(-2.8540197) q[1];
sx q[1];
rz(-0.92929274) q[1];
sx q[1];
rz(0.25921777) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1179781) q[0];
sx q[0];
rz(-3.0833457) q[0];
sx q[0];
rz(-1.9969166) q[0];
rz(-1.338663) q[2];
sx q[2];
rz(-1.4414424) q[2];
sx q[2];
rz(-2.8393313) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99154918) q[1];
sx q[1];
rz(-0.75356149) q[1];
sx q[1];
rz(1.4808473) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76595347) q[3];
sx q[3];
rz(-1.4608835) q[3];
sx q[3];
rz(3.0725057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5549434) q[2];
sx q[2];
rz(-1.2536851) q[2];
sx q[2];
rz(2.9602236) q[2];
rz(0.3950611) q[3];
sx q[3];
rz(-2.919988) q[3];
sx q[3];
rz(2.1612397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3996537) q[0];
sx q[0];
rz(-1.548865) q[0];
sx q[0];
rz(-2.7594866) q[0];
rz(2.8451879) q[1];
sx q[1];
rz(-1.0045241) q[1];
sx q[1];
rz(1.2688676) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60703245) q[0];
sx q[0];
rz(-1.9968613) q[0];
sx q[0];
rz(-2.0924278) q[0];
x q[1];
rz(3.0090055) q[2];
sx q[2];
rz(-2.2397723) q[2];
sx q[2];
rz(-2.7821409) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.813432) q[1];
sx q[1];
rz(-0.34036553) q[1];
sx q[1];
rz(1.5068547) q[1];
rz(2.6878854) q[3];
sx q[3];
rz(-1.6804916) q[3];
sx q[3];
rz(-0.060893313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2962013) q[2];
sx q[2];
rz(-2.5808344) q[2];
sx q[2];
rz(-2.7517547) q[2];
rz(1.8440394) q[3];
sx q[3];
rz(-1.1417979) q[3];
sx q[3];
rz(-0.84387422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.98904499) q[0];
sx q[0];
rz(-1.7392409) q[0];
sx q[0];
rz(1.637511) q[0];
rz(-1.7851495) q[1];
sx q[1];
rz(-2.103613) q[1];
sx q[1];
rz(1.6115859) q[1];
rz(2.6806954) q[2];
sx q[2];
rz(-0.97958889) q[2];
sx q[2];
rz(-0.19340672) q[2];
rz(2.3572902) q[3];
sx q[3];
rz(-1.3275066) q[3];
sx q[3];
rz(2.2307997) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

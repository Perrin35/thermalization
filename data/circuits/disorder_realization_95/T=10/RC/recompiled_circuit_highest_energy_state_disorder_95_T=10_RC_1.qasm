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
rz(0.79701841) q[0];
sx q[0];
rz(0.92172829) q[0];
sx q[0];
rz(14.447639) q[0];
rz(-2.4951275) q[1];
sx q[1];
rz(-0.87352455) q[1];
sx q[1];
rz(0.33198196) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2811738) q[0];
sx q[0];
rz(-1.5462942) q[0];
sx q[0];
rz(1.9936731) q[0];
x q[1];
rz(1.1951564) q[2];
sx q[2];
rz(-0.90919288) q[2];
sx q[2];
rz(2.313569) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0003854) q[1];
sx q[1];
rz(-1.6024593) q[1];
sx q[1];
rz(-0.75891665) q[1];
rz(-1.6360248) q[3];
sx q[3];
rz(-2.0121775) q[3];
sx q[3];
rz(-2.7595487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6648286) q[2];
sx q[2];
rz(-2.2085184) q[2];
sx q[2];
rz(2.942371) q[2];
rz(-0.18579379) q[3];
sx q[3];
rz(-1.7265065) q[3];
sx q[3];
rz(-1.7004405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.324447) q[0];
sx q[0];
rz(-2.6528093) q[0];
sx q[0];
rz(2.4031438) q[0];
rz(1.6959408) q[1];
sx q[1];
rz(-1.4001458) q[1];
sx q[1];
rz(0.48283985) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8890742) q[0];
sx q[0];
rz(-1.7216428) q[0];
sx q[0];
rz(2.9993527) q[0];
rz(2.966283) q[2];
sx q[2];
rz(-1.2273437) q[2];
sx q[2];
rz(-0.94409787) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0898599) q[1];
sx q[1];
rz(-0.75354939) q[1];
sx q[1];
rz(2.4315351) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0121423) q[3];
sx q[3];
rz(-1.7098655) q[3];
sx q[3];
rz(-2.6731051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.26756755) q[2];
sx q[2];
rz(-1.678062) q[2];
sx q[2];
rz(-1.2919424) q[2];
rz(-2.0180295) q[3];
sx q[3];
rz(-0.70190391) q[3];
sx q[3];
rz(-1.4152214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32347754) q[0];
sx q[0];
rz(-2.1577142) q[0];
sx q[0];
rz(1.4418607) q[0];
rz(-1.2280751) q[1];
sx q[1];
rz(-1.9107198) q[1];
sx q[1];
rz(-2.0679811) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.613205) q[0];
sx q[0];
rz(-2.3484485) q[0];
sx q[0];
rz(2.6499477) q[0];
rz(-pi) q[1];
rz(-0.19811689) q[2];
sx q[2];
rz(-0.92495944) q[2];
sx q[2];
rz(0.4465296) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0328999) q[1];
sx q[1];
rz(-1.5376904) q[1];
sx q[1];
rz(-0.86825235) q[1];
rz(-pi) q[2];
rz(1.4074247) q[3];
sx q[3];
rz(-1.1187807) q[3];
sx q[3];
rz(0.79240914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29232612) q[2];
sx q[2];
rz(-0.30305114) q[2];
sx q[2];
rz(-1.0230052) q[2];
rz(-3.0813713) q[3];
sx q[3];
rz(-1.9950208) q[3];
sx q[3];
rz(2.4165418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4802454) q[0];
sx q[0];
rz(-0.1460954) q[0];
sx q[0];
rz(-2.6222141) q[0];
rz(0.6005148) q[1];
sx q[1];
rz(-2.0527716) q[1];
sx q[1];
rz(-0.75469887) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9924644) q[0];
sx q[0];
rz(-2.0880648) q[0];
sx q[0];
rz(-2.7658092) q[0];
x q[1];
rz(-2.9982996) q[2];
sx q[2];
rz(-0.87277647) q[2];
sx q[2];
rz(-2.3868274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.80513152) q[1];
sx q[1];
rz(-1.0829003) q[1];
sx q[1];
rz(-2.7837903) q[1];
rz(-pi) q[2];
rz(0.58784501) q[3];
sx q[3];
rz(-2.4949772) q[3];
sx q[3];
rz(-2.4481096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34919136) q[2];
sx q[2];
rz(-2.3442522) q[2];
sx q[2];
rz(-0.1758197) q[2];
rz(-2.0154121) q[3];
sx q[3];
rz(-1.7837985) q[3];
sx q[3];
rz(-3.0774097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46614161) q[0];
sx q[0];
rz(-1.4956681) q[0];
sx q[0];
rz(-0.58188907) q[0];
rz(-1.9582845) q[1];
sx q[1];
rz(-0.71154037) q[1];
sx q[1];
rz(-1.8875095) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81800705) q[0];
sx q[0];
rz(-2.3773851) q[0];
sx q[0];
rz(-2.0463405) q[0];
rz(-pi) q[1];
rz(-2.3603201) q[2];
sx q[2];
rz(-0.63129497) q[2];
sx q[2];
rz(-2.6897813) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0662567) q[1];
sx q[1];
rz(-1.3622704) q[1];
sx q[1];
rz(-1.2700524) q[1];
x q[2];
rz(-1.759911) q[3];
sx q[3];
rz(-1.0737891) q[3];
sx q[3];
rz(0.86732098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1648272) q[2];
sx q[2];
rz(-0.81563121) q[2];
sx q[2];
rz(-2.7867479) q[2];
rz(-1.9918848) q[3];
sx q[3];
rz(-1.0747654) q[3];
sx q[3];
rz(0.84158516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0717764) q[0];
sx q[0];
rz(-2.8700097) q[0];
sx q[0];
rz(3.1014882) q[0];
rz(1.4647723) q[1];
sx q[1];
rz(-2.1720839) q[1];
sx q[1];
rz(2.6693595) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036333648) q[0];
sx q[0];
rz(-1.3091631) q[0];
sx q[0];
rz(2.9519269) q[0];
rz(-0.91798616) q[2];
sx q[2];
rz(-1.5734104) q[2];
sx q[2];
rz(1.4823584) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4368867) q[1];
sx q[1];
rz(-0.44389899) q[1];
sx q[1];
rz(2.2870334) q[1];
rz(-1.8566441) q[3];
sx q[3];
rz(-2.3116391) q[3];
sx q[3];
rz(2.2835116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8670696) q[2];
sx q[2];
rz(-2.3752866) q[2];
sx q[2];
rz(1.4806032) q[2];
rz(-0.04145043) q[3];
sx q[3];
rz(-2.4292414) q[3];
sx q[3];
rz(-2.1672772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90254766) q[0];
sx q[0];
rz(-0.79896611) q[0];
sx q[0];
rz(2.3833185) q[0];
rz(2.7039418) q[1];
sx q[1];
rz(-1.2145019) q[1];
sx q[1];
rz(-0.95380107) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8945635) q[0];
sx q[0];
rz(-1.5252542) q[0];
sx q[0];
rz(-1.9532502) q[0];
x q[1];
rz(0.66949797) q[2];
sx q[2];
rz(-0.84247103) q[2];
sx q[2];
rz(2.5491722) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.071324997) q[1];
sx q[1];
rz(-2.4355781) q[1];
sx q[1];
rz(-1.6672177) q[1];
x q[2];
rz(-1.1552951) q[3];
sx q[3];
rz(-1.3127919) q[3];
sx q[3];
rz(0.33694944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7268251) q[2];
sx q[2];
rz(-0.32764062) q[2];
sx q[2];
rz(-0.3978351) q[2];
rz(-1.2658524) q[3];
sx q[3];
rz(-1.7222907) q[3];
sx q[3];
rz(-0.46441594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1606814) q[0];
sx q[0];
rz(-2.2137764) q[0];
sx q[0];
rz(0.30447793) q[0];
rz(-1.132384) q[1];
sx q[1];
rz(-0.67190036) q[1];
sx q[1];
rz(-2.0679881) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17504263) q[0];
sx q[0];
rz(-1.5538408) q[0];
sx q[0];
rz(1.5477563) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4173293) q[2];
sx q[2];
rz(-1.9021735) q[2];
sx q[2];
rz(-1.1294236) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1731092) q[1];
sx q[1];
rz(-2.1595528) q[1];
sx q[1];
rz(-2.8485879) q[1];
rz(0.066929265) q[3];
sx q[3];
rz(-0.64475497) q[3];
sx q[3];
rz(-0.88551846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4094746) q[2];
sx q[2];
rz(-1.2827337) q[2];
sx q[2];
rz(3.0916302) q[2];
rz(-0.91160715) q[3];
sx q[3];
rz(-2.215569) q[3];
sx q[3];
rz(-0.91134206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6066345) q[0];
sx q[0];
rz(-1.7462523) q[0];
sx q[0];
rz(0.362679) q[0];
rz(2.4210988) q[1];
sx q[1];
rz(-0.81169218) q[1];
sx q[1];
rz(1.1788751) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4587817) q[0];
sx q[0];
rz(-1.0014373) q[0];
sx q[0];
rz(-2.1853133) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8209053) q[2];
sx q[2];
rz(-1.4960519) q[2];
sx q[2];
rz(-0.2439258) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7866582) q[1];
sx q[1];
rz(-1.400507) q[1];
sx q[1];
rz(-3.0436467) q[1];
rz(2.5963327) q[3];
sx q[3];
rz(-1.2125748) q[3];
sx q[3];
rz(-2.6521366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1371586) q[2];
sx q[2];
rz(-2.4538071) q[2];
sx q[2];
rz(1.2031817) q[2];
rz(3.1249937) q[3];
sx q[3];
rz(-1.6417475) q[3];
sx q[3];
rz(-1.871292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.5430629) q[0];
sx q[0];
rz(-0.5235343) q[0];
sx q[0];
rz(-1.4434927) q[0];
rz(0.47830018) q[1];
sx q[1];
rz(-2.4595478) q[1];
sx q[1];
rz(1.9313448) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.011199) q[0];
sx q[0];
rz(-0.78576311) q[0];
sx q[0];
rz(-2.4146621) q[0];
x q[1];
rz(-0.6503251) q[2];
sx q[2];
rz(-1.4286094) q[2];
sx q[2];
rz(1.1917758) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.12809556) q[1];
sx q[1];
rz(-1.5186608) q[1];
sx q[1];
rz(-0.32491046) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5860255) q[3];
sx q[3];
rz(-2.7507493) q[3];
sx q[3];
rz(1.9236717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9532507) q[2];
sx q[2];
rz(-2.5276999) q[2];
sx q[2];
rz(0.71315145) q[2];
rz(-0.6330511) q[3];
sx q[3];
rz(-0.96708599) q[3];
sx q[3];
rz(-2.2658394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1322094) q[0];
sx q[0];
rz(-1.3306946) q[0];
sx q[0];
rz(2.7015986) q[0];
rz(0.72883365) q[1];
sx q[1];
rz(-0.76580096) q[1];
sx q[1];
rz(-2.0892807) q[1];
rz(-3.0265774) q[2];
sx q[2];
rz(-1.1830405) q[2];
sx q[2];
rz(-1.9934275) q[2];
rz(-1.1951554) q[3];
sx q[3];
rz(-0.79938625) q[3];
sx q[3];
rz(1.1342794) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

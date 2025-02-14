OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.96707764) q[0];
sx q[0];
rz(-1.4580026) q[0];
sx q[0];
rz(-2.8414677) q[0];
rz(-1.9441654) q[1];
sx q[1];
rz(-1.5265042) q[1];
sx q[1];
rz(-1.6405029) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1688582) q[0];
sx q[0];
rz(-1.5669113) q[0];
sx q[0];
rz(0.015104276) q[0];
x q[1];
rz(2.3387942) q[2];
sx q[2];
rz(-1.8539696) q[2];
sx q[2];
rz(0.51529037) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38084322) q[1];
sx q[1];
rz(-0.92807209) q[1];
sx q[1];
rz(-1.8442783) q[1];
rz(-1.5656823) q[3];
sx q[3];
rz(-0.76200125) q[3];
sx q[3];
rz(0.68540547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40452051) q[2];
sx q[2];
rz(-1.1417192) q[2];
sx q[2];
rz(-2.4386151) q[2];
rz(-0.56973714) q[3];
sx q[3];
rz(-2.2629786) q[3];
sx q[3];
rz(-1.5841293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1341781) q[0];
sx q[0];
rz(-0.61892048) q[0];
sx q[0];
rz(-1.2331569) q[0];
rz(-2.5470219) q[1];
sx q[1];
rz(-1.0392799) q[1];
sx q[1];
rz(-2.0436683) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6195589) q[0];
sx q[0];
rz(-2.9298943) q[0];
sx q[0];
rz(-1.2489399) q[0];
rz(-pi) q[1];
rz(0.64867257) q[2];
sx q[2];
rz(-2.5498516) q[2];
sx q[2];
rz(-1.6206738) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6095396) q[1];
sx q[1];
rz(-1.5960448) q[1];
sx q[1];
rz(1.6482501) q[1];
rz(-pi) q[2];
rz(3.0497562) q[3];
sx q[3];
rz(-0.79376924) q[3];
sx q[3];
rz(2.5240999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.42830959) q[2];
sx q[2];
rz(-1.8417336) q[2];
sx q[2];
rz(1.606288) q[2];
rz(0.71896583) q[3];
sx q[3];
rz(-2.1956367) q[3];
sx q[3];
rz(-0.21397056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83800256) q[0];
sx q[0];
rz(-0.8135697) q[0];
sx q[0];
rz(-2.549262) q[0];
rz(-1.8968286) q[1];
sx q[1];
rz(-1.1400305) q[1];
sx q[1];
rz(-2.9959784) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8728565) q[0];
sx q[0];
rz(-1.4352192) q[0];
sx q[0];
rz(-0.4879293) q[0];
x q[1];
rz(0.8823186) q[2];
sx q[2];
rz(-1.6711418) q[2];
sx q[2];
rz(2.0813297) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1864724) q[1];
sx q[1];
rz(-1.2870409) q[1];
sx q[1];
rz(2.1452745) q[1];
rz(-pi) q[2];
rz(-2.0592732) q[3];
sx q[3];
rz(-1.8104189) q[3];
sx q[3];
rz(2.3828196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40272063) q[2];
sx q[2];
rz(-1.0307743) q[2];
sx q[2];
rz(-1.2325475) q[2];
rz(-2.9018719) q[3];
sx q[3];
rz(-2.0870049) q[3];
sx q[3];
rz(2.6565552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1211014) q[0];
sx q[0];
rz(-2.4395269) q[0];
sx q[0];
rz(-0.030601587) q[0];
rz(2.5864511) q[1];
sx q[1];
rz(-1.4786913) q[1];
sx q[1];
rz(2.0487002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8776833) q[0];
sx q[0];
rz(-1.0315686) q[0];
sx q[0];
rz(-0.79099057) q[0];
rz(-pi) q[1];
rz(1.5116772) q[2];
sx q[2];
rz(-0.80628866) q[2];
sx q[2];
rz(-1.8102243) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.79697414) q[1];
sx q[1];
rz(-0.94616468) q[1];
sx q[1];
rz(0.25298869) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2956156) q[3];
sx q[3];
rz(-1.032077) q[3];
sx q[3];
rz(1.908345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0577724) q[2];
sx q[2];
rz(-1.7744935) q[2];
sx q[2];
rz(0.27285451) q[2];
rz(-1.3911635) q[3];
sx q[3];
rz(-0.83943668) q[3];
sx q[3];
rz(-0.951989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4680173) q[0];
sx q[0];
rz(-1.8033569) q[0];
sx q[0];
rz(-1.2982298) q[0];
rz(-0.6908373) q[1];
sx q[1];
rz(-1.1996484) q[1];
sx q[1];
rz(1.615049) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95193255) q[0];
sx q[0];
rz(-2.2554923) q[0];
sx q[0];
rz(0.71345274) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35299086) q[2];
sx q[2];
rz(-1.8367177) q[2];
sx q[2];
rz(2.860255) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95350515) q[1];
sx q[1];
rz(-1.5966478) q[1];
sx q[1];
rz(0.8698747) q[1];
rz(-0.68256179) q[3];
sx q[3];
rz(-0.45387156) q[3];
sx q[3];
rz(-1.0938494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7319298) q[2];
sx q[2];
rz(-0.80594984) q[2];
sx q[2];
rz(0.16732495) q[2];
rz(1.3903728) q[3];
sx q[3];
rz(-2.6087587) q[3];
sx q[3];
rz(-0.65756857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996027) q[0];
sx q[0];
rz(-2.7775192) q[0];
sx q[0];
rz(1.863119) q[0];
rz(1.4356042) q[1];
sx q[1];
rz(-1.6537138) q[1];
sx q[1];
rz(-2.7659168) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6784793) q[0];
sx q[0];
rz(-2.7857384) q[0];
sx q[0];
rz(-0.87219091) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2907545) q[2];
sx q[2];
rz(-0.22432835) q[2];
sx q[2];
rz(1.6825324) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73114016) q[1];
sx q[1];
rz(-1.8888432) q[1];
sx q[1];
rz(-0.2348301) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75590862) q[3];
sx q[3];
rz(-1.3918073) q[3];
sx q[3];
rz(0.32839113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0756691) q[2];
sx q[2];
rz(-1.064294) q[2];
sx q[2];
rz(-2.2806878) q[2];
rz(-1.143645) q[3];
sx q[3];
rz(-0.30503169) q[3];
sx q[3];
rz(0.27826571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1160195) q[0];
sx q[0];
rz(-1.3207859) q[0];
sx q[0];
rz(-0.75468165) q[0];
rz(0.93983752) q[1];
sx q[1];
rz(-0.87516963) q[1];
sx q[1];
rz(-1.2831203) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96961731) q[0];
sx q[0];
rz(-1.2640868) q[0];
sx q[0];
rz(0.0071010751) q[0];
rz(-pi) q[1];
rz(-1.1129598) q[2];
sx q[2];
rz(-1.2414059) q[2];
sx q[2];
rz(3.0731346) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2312113) q[1];
sx q[1];
rz(-1.8613792) q[1];
sx q[1];
rz(2.742275) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96638443) q[3];
sx q[3];
rz(-2.3570286) q[3];
sx q[3];
rz(-1.094629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54520404) q[2];
sx q[2];
rz(-1.117492) q[2];
sx q[2];
rz(-2.0580573) q[2];
rz(-0.74294535) q[3];
sx q[3];
rz(-1.6094306) q[3];
sx q[3];
rz(-1.4325745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525986) q[0];
sx q[0];
rz(-0.78482634) q[0];
sx q[0];
rz(2.3866744) q[0];
rz(-2.6424291) q[1];
sx q[1];
rz(-0.9639591) q[1];
sx q[1];
rz(2.4430433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7592418) q[0];
sx q[0];
rz(-0.97433358) q[0];
sx q[0];
rz(-1.9475219) q[0];
rz(1.7153984) q[2];
sx q[2];
rz(-1.1100262) q[2];
sx q[2];
rz(0.20587155) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.10741988) q[1];
sx q[1];
rz(-1.0386969) q[1];
sx q[1];
rz(-3.1183395) q[1];
rz(-1.1517161) q[3];
sx q[3];
rz(-1.569783) q[3];
sx q[3];
rz(0.59194293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2373206) q[2];
sx q[2];
rz(-2.1492683) q[2];
sx q[2];
rz(-0.80238706) q[2];
rz(2.5896416) q[3];
sx q[3];
rz(-2.3410083) q[3];
sx q[3];
rz(-2.5210181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9629843) q[0];
sx q[0];
rz(-1.7010138) q[0];
sx q[0];
rz(-2.0340023) q[0];
rz(1.7604609) q[1];
sx q[1];
rz(-1.7778722) q[1];
sx q[1];
rz(1.6345056) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.749562) q[0];
sx q[0];
rz(-1.9641341) q[0];
sx q[0];
rz(-1.419072) q[0];
rz(-2.2198943) q[2];
sx q[2];
rz(-2.6423811) q[2];
sx q[2];
rz(-1.2865024) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6927106) q[1];
sx q[1];
rz(-1.7160985) q[1];
sx q[1];
rz(2.7753745) q[1];
rz(0.82689311) q[3];
sx q[3];
rz(-1.6312459) q[3];
sx q[3];
rz(-1.353144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0570809) q[2];
sx q[2];
rz(-2.9456186) q[2];
sx q[2];
rz(1.2723119) q[2];
rz(-1.48014) q[3];
sx q[3];
rz(-1.1353761) q[3];
sx q[3];
rz(0.60012668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8965974) q[0];
sx q[0];
rz(-0.56140459) q[0];
sx q[0];
rz(-1.2835314) q[0];
rz(0.01783477) q[1];
sx q[1];
rz(-0.63648883) q[1];
sx q[1];
rz(3.0160115) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20480072) q[0];
sx q[0];
rz(-1.5849587) q[0];
sx q[0];
rz(-1.0165748) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1000041) q[2];
sx q[2];
rz(-0.94816899) q[2];
sx q[2];
rz(-2.560844) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8931742) q[1];
sx q[1];
rz(-1.1496953) q[1];
sx q[1];
rz(1.1014654) q[1];
x q[2];
rz(-2.9405648) q[3];
sx q[3];
rz(-1.8327692) q[3];
sx q[3];
rz(2.6233545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77832001) q[2];
sx q[2];
rz(-1.2691701) q[2];
sx q[2];
rz(-2.3385091) q[2];
rz(2.965029) q[3];
sx q[3];
rz(-1.8162138) q[3];
sx q[3];
rz(2.363502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17978996) q[0];
sx q[0];
rz(-2.640124) q[0];
sx q[0];
rz(-1.5928706) q[0];
rz(1.8190307) q[1];
sx q[1];
rz(-1.7772728) q[1];
sx q[1];
rz(-1.5900236) q[1];
rz(-1.6336468) q[2];
sx q[2];
rz(-1.812915) q[2];
sx q[2];
rz(-1.6356638) q[2];
rz(-0.17038067) q[3];
sx q[3];
rz(-2.7994886) q[3];
sx q[3];
rz(-2.8240374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

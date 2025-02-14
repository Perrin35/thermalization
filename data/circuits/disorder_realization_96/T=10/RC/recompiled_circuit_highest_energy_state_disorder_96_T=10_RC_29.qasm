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
rz(0.32956707) q[0];
sx q[0];
rz(-2.1962533) q[0];
sx q[0];
rz(-3.1402631) q[0];
rz(0.20345649) q[1];
sx q[1];
rz(-0.23509547) q[1];
sx q[1];
rz(-2.5703365) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1648096) q[0];
sx q[0];
rz(-2.342988) q[0];
sx q[0];
rz(0.66508663) q[0];
rz(0.6735835) q[2];
sx q[2];
rz(-0.24112186) q[2];
sx q[2];
rz(2.4862289) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9524051) q[1];
sx q[1];
rz(-2.6509977) q[1];
sx q[1];
rz(0.036019939) q[1];
rz(-0.68526973) q[3];
sx q[3];
rz(-2.5766386) q[3];
sx q[3];
rz(-0.23811114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.098238952) q[2];
sx q[2];
rz(-1.9184155) q[2];
sx q[2];
rz(-2.5377048) q[2];
rz(-2.1308925) q[3];
sx q[3];
rz(-0.5637919) q[3];
sx q[3];
rz(-0.92196661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0826223) q[0];
sx q[0];
rz(-2.048546) q[0];
sx q[0];
rz(3.1329204) q[0];
rz(-1.9827093) q[1];
sx q[1];
rz(-0.68717879) q[1];
sx q[1];
rz(3.029356) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85189795) q[0];
sx q[0];
rz(-1.0852975) q[0];
sx q[0];
rz(1.2458891) q[0];
rz(1.6084163) q[2];
sx q[2];
rz(-0.92447399) q[2];
sx q[2];
rz(-0.96137709) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.018377233) q[1];
sx q[1];
rz(-2.0615512) q[1];
sx q[1];
rz(-1.823455) q[1];
rz(-pi) q[2];
rz(2.8772164) q[3];
sx q[3];
rz(-1.2942787) q[3];
sx q[3];
rz(0.087527601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.23942854) q[2];
sx q[2];
rz(-1.3210693) q[2];
sx q[2];
rz(0.93977896) q[2];
rz(0.93722614) q[3];
sx q[3];
rz(-1.3746494) q[3];
sx q[3];
rz(-0.37041131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20060191) q[0];
sx q[0];
rz(-2.2414099) q[0];
sx q[0];
rz(-2.5584333) q[0];
rz(-2.0896301) q[1];
sx q[1];
rz(-2.5556892) q[1];
sx q[1];
rz(-3.1254056) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72618851) q[0];
sx q[0];
rz(-2.4016018) q[0];
sx q[0];
rz(1.9053024) q[0];
rz(-1.2643686) q[2];
sx q[2];
rz(-1.2867905) q[2];
sx q[2];
rz(2.4302182) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5213373) q[1];
sx q[1];
rz(-0.19591051) q[1];
sx q[1];
rz(-0.98448648) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2233769) q[3];
sx q[3];
rz(-0.86091061) q[3];
sx q[3];
rz(0.48274279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6215324) q[2];
sx q[2];
rz(-1.3935057) q[2];
sx q[2];
rz(0.072546093) q[2];
rz(1.0091311) q[3];
sx q[3];
rz(-0.79807177) q[3];
sx q[3];
rz(2.9062241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9811454) q[0];
sx q[0];
rz(-2.7837842) q[0];
sx q[0];
rz(-0.88448802) q[0];
rz(-1.5011939) q[1];
sx q[1];
rz(-1.4720935) q[1];
sx q[1];
rz(-2.5515058) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5948759) q[0];
sx q[0];
rz(-1.1980431) q[0];
sx q[0];
rz(1.0327505) q[0];
rz(-0.8832189) q[2];
sx q[2];
rz(-0.35053634) q[2];
sx q[2];
rz(0.8274629) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6650538) q[1];
sx q[1];
rz(-2.0609984) q[1];
sx q[1];
rz(0.3785822) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6904126) q[3];
sx q[3];
rz(-1.1362459) q[3];
sx q[3];
rz(2.0903319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0044452) q[2];
sx q[2];
rz(-1.659617) q[2];
sx q[2];
rz(1.2443292) q[2];
rz(-1.7826049) q[3];
sx q[3];
rz(-2.6748896) q[3];
sx q[3];
rz(-2.9984503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8878079) q[0];
sx q[0];
rz(-0.65685993) q[0];
sx q[0];
rz(-0.3669056) q[0];
rz(1.7985571) q[1];
sx q[1];
rz(-0.29452205) q[1];
sx q[1];
rz(-1.8273182) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9132694) q[0];
sx q[0];
rz(-2.0073038) q[0];
sx q[0];
rz(-0.69605791) q[0];
rz(-pi) q[1];
rz(-1.8601244) q[2];
sx q[2];
rz(-1.1035894) q[2];
sx q[2];
rz(1.510646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6180743) q[1];
sx q[1];
rz(-2.1879751) q[1];
sx q[1];
rz(2.7568222) q[1];
rz(-pi) q[2];
rz(2.3288127) q[3];
sx q[3];
rz(-2.644499) q[3];
sx q[3];
rz(1.415783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5985976) q[2];
sx q[2];
rz(-2.6948805) q[2];
sx q[2];
rz(2.561595) q[2];
rz(-0.18481208) q[3];
sx q[3];
rz(-1.5671268) q[3];
sx q[3];
rz(-0.68661657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0762894) q[0];
sx q[0];
rz(-0.46173254) q[0];
sx q[0];
rz(0.28934685) q[0];
rz(3.0416327) q[1];
sx q[1];
rz(-0.83671612) q[1];
sx q[1];
rz(2.3438556) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8160523) q[0];
sx q[0];
rz(-2.1959248) q[0];
sx q[0];
rz(2.9751871) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50578588) q[2];
sx q[2];
rz(-2.9120076) q[2];
sx q[2];
rz(-1.2248271) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3926852) q[1];
sx q[1];
rz(-1.8659544) q[1];
sx q[1];
rz(-1.6243851) q[1];
rz(-1.6838668) q[3];
sx q[3];
rz(-0.73323876) q[3];
sx q[3];
rz(-1.7565786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.37819698) q[2];
sx q[2];
rz(-1.244647) q[2];
sx q[2];
rz(0.81983105) q[2];
rz(1.648692) q[3];
sx q[3];
rz(-0.35455743) q[3];
sx q[3];
rz(-2.7259887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019526871) q[0];
sx q[0];
rz(-0.30240914) q[0];
sx q[0];
rz(0.31461883) q[0];
rz(-2.4954691) q[1];
sx q[1];
rz(-1.6982634) q[1];
sx q[1];
rz(-3.019928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2631609) q[0];
sx q[0];
rz(-1.9515349) q[0];
sx q[0];
rz(2.8436766) q[0];
x q[1];
rz(-0.67331566) q[2];
sx q[2];
rz(-0.91860572) q[2];
sx q[2];
rz(-2.7567692) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0371656) q[1];
sx q[1];
rz(-2.011877) q[1];
sx q[1];
rz(1.7165401) q[1];
rz(-2.9470351) q[3];
sx q[3];
rz(-0.97134198) q[3];
sx q[3];
rz(2.636858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6803711) q[2];
sx q[2];
rz(-1.17522) q[2];
sx q[2];
rz(2.1628974) q[2];
rz(2.3693502) q[3];
sx q[3];
rz(-0.52670908) q[3];
sx q[3];
rz(1.0129119) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1145645) q[0];
sx q[0];
rz(-1.9975198) q[0];
sx q[0];
rz(-1.3042599) q[0];
rz(3.0005241) q[1];
sx q[1];
rz(-0.27962676) q[1];
sx q[1];
rz(-1.28349) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6161364) q[0];
sx q[0];
rz(-1.5342664) q[0];
sx q[0];
rz(1.7869048) q[0];
rz(-pi) q[1];
rz(-1.263563) q[2];
sx q[2];
rz(-0.31401411) q[2];
sx q[2];
rz(-2.8548714) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5775958) q[1];
sx q[1];
rz(-1.8297927) q[1];
sx q[1];
rz(-1.6902655) q[1];
rz(-pi) q[2];
x q[2];
rz(0.047824511) q[3];
sx q[3];
rz(-2.3615814) q[3];
sx q[3];
rz(-0.081950233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.86113247) q[2];
sx q[2];
rz(-2.2304163) q[2];
sx q[2];
rz(2.148441) q[2];
rz(-1.447621) q[3];
sx q[3];
rz(-1.2074892) q[3];
sx q[3];
rz(-1.2710458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.3250378) q[0];
sx q[0];
rz(-0.10295454) q[0];
sx q[0];
rz(-0.8763985) q[0];
rz(1.5429629) q[1];
sx q[1];
rz(-1.8537268) q[1];
sx q[1];
rz(-1.1184982) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0792313) q[0];
sx q[0];
rz(-0.97506279) q[0];
sx q[0];
rz(1.1337639) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.411941) q[2];
sx q[2];
rz(-1.083598) q[2];
sx q[2];
rz(-1.2318357) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4913018) q[1];
sx q[1];
rz(-0.87556616) q[1];
sx q[1];
rz(-2.9952666) q[1];
rz(-1.275609) q[3];
sx q[3];
rz(-0.94691197) q[3];
sx q[3];
rz(-2.5518887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8886275) q[2];
sx q[2];
rz(-2.1972563) q[2];
sx q[2];
rz(-1.2523119) q[2];
rz(1.1014994) q[3];
sx q[3];
rz(-1.7606198) q[3];
sx q[3];
rz(-1.8496877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0970704) q[0];
sx q[0];
rz(-0.047053311) q[0];
sx q[0];
rz(-2.4142921) q[0];
rz(2.3749088) q[1];
sx q[1];
rz(-2.2946281) q[1];
sx q[1];
rz(1.9511706) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81556126) q[0];
sx q[0];
rz(-2.4056068) q[0];
sx q[0];
rz(0.18819564) q[0];
x q[1];
rz(1.4486362) q[2];
sx q[2];
rz(-2.4776931) q[2];
sx q[2];
rz(0.66421504) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7720346) q[1];
sx q[1];
rz(-1.2698962) q[1];
sx q[1];
rz(1.8298261) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99440688) q[3];
sx q[3];
rz(-1.1383447) q[3];
sx q[3];
rz(1.2131888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7850354) q[2];
sx q[2];
rz(-1.0871474) q[2];
sx q[2];
rz(-1.9800775) q[2];
rz(1.131743) q[3];
sx q[3];
rz(-0.83860207) q[3];
sx q[3];
rz(2.13818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3840735) q[0];
sx q[0];
rz(-2.2390371) q[0];
sx q[0];
rz(-0.83691103) q[0];
rz(1.8867672) q[1];
sx q[1];
rz(-1.8950987) q[1];
sx q[1];
rz(-1.7991039) q[1];
rz(2.135545) q[2];
sx q[2];
rz(-1.3443832) q[2];
sx q[2];
rz(0.98510712) q[2];
rz(0.014160362) q[3];
sx q[3];
rz(-1.7668528) q[3];
sx q[3];
rz(1.7535221) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

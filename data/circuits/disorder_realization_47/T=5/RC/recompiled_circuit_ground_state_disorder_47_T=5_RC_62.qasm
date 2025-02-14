OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1563675) q[0];
sx q[0];
rz(-1.2824143) q[0];
sx q[0];
rz(-0.13089827) q[0];
rz(-2.6137597) q[1];
sx q[1];
rz(-0.37017828) q[1];
sx q[1];
rz(-2.9642677) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41577121) q[0];
sx q[0];
rz(-0.81233378) q[0];
sx q[0];
rz(1.2052631) q[0];
rz(-pi) q[1];
rz(0.22819569) q[2];
sx q[2];
rz(-1.2773809) q[2];
sx q[2];
rz(0.55390893) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0874112) q[1];
sx q[1];
rz(-2.0435395) q[1];
sx q[1];
rz(2.8768538) q[1];
x q[2];
rz(-2.1471094) q[3];
sx q[3];
rz(-1.3399933) q[3];
sx q[3];
rz(2.1891914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6191972) q[2];
sx q[2];
rz(-1.088524) q[2];
sx q[2];
rz(1.3414475) q[2];
rz(1.3522735) q[3];
sx q[3];
rz(-1.6100581) q[3];
sx q[3];
rz(1.8488098) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.751048) q[0];
sx q[0];
rz(-1.229267) q[0];
sx q[0];
rz(2.6265662) q[0];
rz(0.36188778) q[1];
sx q[1];
rz(-1.7444976) q[1];
sx q[1];
rz(-2.1451758) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5549042) q[0];
sx q[0];
rz(-1.3350272) q[0];
sx q[0];
rz(0.20637189) q[0];
x q[1];
rz(1.9356807) q[2];
sx q[2];
rz(-2.6807086) q[2];
sx q[2];
rz(0.63096607) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39033088) q[1];
sx q[1];
rz(-1.0511025) q[1];
sx q[1];
rz(0.10607509) q[1];
rz(-pi) q[2];
rz(1.0330233) q[3];
sx q[3];
rz(-1.1837848) q[3];
sx q[3];
rz(-1.3082711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4126052) q[2];
sx q[2];
rz(-2.1952486) q[2];
sx q[2];
rz(1.4808572) q[2];
rz(2.6442773) q[3];
sx q[3];
rz(-2.1931084) q[3];
sx q[3];
rz(2.4174387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5363252) q[0];
sx q[0];
rz(-1.8999506) q[0];
sx q[0];
rz(1.9507677) q[0];
rz(2.7438927) q[1];
sx q[1];
rz(-1.5813446) q[1];
sx q[1];
rz(2.2427028) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54972285) q[0];
sx q[0];
rz(-1.6373349) q[0];
sx q[0];
rz(-1.2308685) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4853046) q[2];
sx q[2];
rz(-1.3807266) q[2];
sx q[2];
rz(-0.23641931) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.11375107) q[1];
sx q[1];
rz(-1.6154624) q[1];
sx q[1];
rz(1.0381446) q[1];
rz(-pi) q[2];
rz(-1.8679138) q[3];
sx q[3];
rz(-0.53830244) q[3];
sx q[3];
rz(-0.75002128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1206104) q[2];
sx q[2];
rz(-1.1889428) q[2];
sx q[2];
rz(0.19213842) q[2];
rz(-2.4118679) q[3];
sx q[3];
rz(-0.14484043) q[3];
sx q[3];
rz(-0.63989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9925053) q[0];
sx q[0];
rz(-0.75228107) q[0];
sx q[0];
rz(-1.3080904) q[0];
rz(0.84450841) q[1];
sx q[1];
rz(-1.6788071) q[1];
sx q[1];
rz(2.5194397) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5151241) q[0];
sx q[0];
rz(-0.65931407) q[0];
sx q[0];
rz(-2.4298682) q[0];
rz(-pi) q[1];
rz(2.8573158) q[2];
sx q[2];
rz(-2.948649) q[2];
sx q[2];
rz(0.74954734) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6250861) q[1];
sx q[1];
rz(-1.7343905) q[1];
sx q[1];
rz(1.6229936) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31720576) q[3];
sx q[3];
rz(-1.6807846) q[3];
sx q[3];
rz(0.12817891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.017612351) q[2];
sx q[2];
rz(-1.3717185) q[2];
sx q[2];
rz(2.2464216) q[2];
rz(0.34919843) q[3];
sx q[3];
rz(-0.57217351) q[3];
sx q[3];
rz(-2.9077742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2772086) q[0];
sx q[0];
rz(-1.207749) q[0];
sx q[0];
rz(2.5982017) q[0];
rz(-3.0282989) q[1];
sx q[1];
rz(-1.7430867) q[1];
sx q[1];
rz(1.2921804) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1975511) q[0];
sx q[0];
rz(-2.3359951) q[0];
sx q[0];
rz(1.0260236) q[0];
rz(-1.0903484) q[2];
sx q[2];
rz(-2.4329429) q[2];
sx q[2];
rz(-2.6330269) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1317392) q[1];
sx q[1];
rz(-0.92127548) q[1];
sx q[1];
rz(-2.3765537) q[1];
rz(0.73907799) q[3];
sx q[3];
rz(-2.1391641) q[3];
sx q[3];
rz(3.0773602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.31305227) q[2];
sx q[2];
rz(-1.3581759) q[2];
sx q[2];
rz(-2.6521315) q[2];
rz(2.475907) q[3];
sx q[3];
rz(-0.61167115) q[3];
sx q[3];
rz(2.3556975) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8996443) q[0];
sx q[0];
rz(-2.2388832) q[0];
sx q[0];
rz(-3.0808501) q[0];
rz(1.863106) q[1];
sx q[1];
rz(-0.91438952) q[1];
sx q[1];
rz(1.0260822) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2506977) q[0];
sx q[0];
rz(-1.015626) q[0];
sx q[0];
rz(-2.5722136) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0842956) q[2];
sx q[2];
rz(-0.62556534) q[2];
sx q[2];
rz(-1.4826258) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.68133229) q[1];
sx q[1];
rz(-1.4341913) q[1];
sx q[1];
rz(2.8885319) q[1];
x q[2];
rz(0.52512759) q[3];
sx q[3];
rz(-0.59288247) q[3];
sx q[3];
rz(1.2660932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.25980276) q[2];
sx q[2];
rz(-1.5564352) q[2];
sx q[2];
rz(3.0628824) q[2];
rz(1.4824661) q[3];
sx q[3];
rz(-0.31646287) q[3];
sx q[3];
rz(0.44221529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92597961) q[0];
sx q[0];
rz(-2.7613566) q[0];
sx q[0];
rz(1.444814) q[0];
rz(-0.80698693) q[1];
sx q[1];
rz(-1.746256) q[1];
sx q[1];
rz(-3.1296465) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1769971) q[0];
sx q[0];
rz(-0.93206333) q[0];
sx q[0];
rz(3.1082736) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4623227) q[2];
sx q[2];
rz(-1.0656992) q[2];
sx q[2];
rz(2.6384356) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1216189) q[1];
sx q[1];
rz(-1.6433006) q[1];
sx q[1];
rz(-0.95901476) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9254382) q[3];
sx q[3];
rz(-0.4222479) q[3];
sx q[3];
rz(0.99340445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5936467) q[2];
sx q[2];
rz(-0.30210364) q[2];
sx q[2];
rz(0.83615237) q[2];
rz(-3.0892843) q[3];
sx q[3];
rz(-1.6736504) q[3];
sx q[3];
rz(-0.93594319) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1702105) q[0];
sx q[0];
rz(-0.8571856) q[0];
sx q[0];
rz(-2.6433387) q[0];
rz(-1.0055297) q[1];
sx q[1];
rz(-1.6906831) q[1];
sx q[1];
rz(0.11140579) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83269925) q[0];
sx q[0];
rz(-1.667262) q[0];
sx q[0];
rz(1.1497496) q[0];
rz(1.8209711) q[2];
sx q[2];
rz(-1.9873575) q[2];
sx q[2];
rz(-0.60384679) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3058779) q[1];
sx q[1];
rz(-1.2849766) q[1];
sx q[1];
rz(0.81982433) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4440246) q[3];
sx q[3];
rz(-1.9132735) q[3];
sx q[3];
rz(0.7966744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.19693836) q[2];
sx q[2];
rz(-1.780218) q[2];
sx q[2];
rz(-1.4062175) q[2];
rz(1.1574636) q[3];
sx q[3];
rz(-0.99298733) q[3];
sx q[3];
rz(-1.5520613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.120753) q[0];
sx q[0];
rz(-0.73796219) q[0];
sx q[0];
rz(-1.7027759) q[0];
rz(-0.84398794) q[1];
sx q[1];
rz(-1.3071209) q[1];
sx q[1];
rz(0.2058952) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1263329) q[0];
sx q[0];
rz(-1.5301989) q[0];
sx q[0];
rz(0.6313398) q[0];
rz(3.0768422) q[2];
sx q[2];
rz(-0.68875411) q[2];
sx q[2];
rz(1.9079067) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.89005) q[1];
sx q[1];
rz(-1.6487507) q[1];
sx q[1];
rz(0.13549094) q[1];
rz(0.4959373) q[3];
sx q[3];
rz(-1.933799) q[3];
sx q[3];
rz(-0.41307005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13731185) q[2];
sx q[2];
rz(-1.0575123) q[2];
sx q[2];
rz(-1.5409957) q[2];
rz(0.48504034) q[3];
sx q[3];
rz(-1.3554327) q[3];
sx q[3];
rz(-3.0276827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.6523022) q[0];
sx q[0];
rz(-3.089383) q[0];
sx q[0];
rz(1.8452277) q[0];
rz(2.0908053) q[1];
sx q[1];
rz(-1.4605582) q[1];
sx q[1];
rz(0.60660648) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3171995) q[0];
sx q[0];
rz(-1.6034295) q[0];
sx q[0];
rz(1.4072818) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5818517) q[2];
sx q[2];
rz(-1.4626182) q[2];
sx q[2];
rz(1.6424204) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7672528) q[1];
sx q[1];
rz(-1.2657968) q[1];
sx q[1];
rz(-1.6348331) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38207558) q[3];
sx q[3];
rz(-2.8365072) q[3];
sx q[3];
rz(0.24972734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4198833) q[2];
sx q[2];
rz(-2.0100644) q[2];
sx q[2];
rz(-2.4427872) q[2];
rz(-1.6803668) q[3];
sx q[3];
rz(-2.1278087) q[3];
sx q[3];
rz(-2.6644126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89194311) q[0];
sx q[0];
rz(-1.5105381) q[0];
sx q[0];
rz(1.426209) q[0];
rz(-2.7632948) q[1];
sx q[1];
rz(-2.5143647) q[1];
sx q[1];
rz(-0.73077269) q[1];
rz(-2.8332491) q[2];
sx q[2];
rz(-1.1157805) q[2];
sx q[2];
rz(0.11238712) q[2];
rz(-0.55647464) q[3];
sx q[3];
rz(-1.597098) q[3];
sx q[3];
rz(2.747018) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

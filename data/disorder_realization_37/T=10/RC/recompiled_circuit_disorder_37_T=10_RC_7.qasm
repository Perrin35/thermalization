OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10387575) q[0];
sx q[0];
rz(-1.9394983) q[0];
sx q[0];
rz(1.9934959) q[0];
rz(1.2530874) q[1];
sx q[1];
rz(4.0822786) q[1];
sx q[1];
rz(10.818426) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05605927) q[0];
sx q[0];
rz(-0.30963184) q[0];
sx q[0];
rz(-3.1349896) q[0];
rz(-pi) q[1];
rz(-0.68140985) q[2];
sx q[2];
rz(-1.0076367) q[2];
sx q[2];
rz(1.5208706) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3323101) q[1];
sx q[1];
rz(-0.53762943) q[1];
sx q[1];
rz(-2.5938631) q[1];
rz(-pi) q[2];
rz(-0.2829708) q[3];
sx q[3];
rz(-2.1108147) q[3];
sx q[3];
rz(-2.6178544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6538438) q[2];
sx q[2];
rz(-1.2922492) q[2];
sx q[2];
rz(0.1208819) q[2];
rz(2.9623048) q[3];
sx q[3];
rz(-0.59569734) q[3];
sx q[3];
rz(-2.9860935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0497465) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(3.0088186) q[0];
rz(-1.4615387) q[1];
sx q[1];
rz(-1.5802054) q[1];
sx q[1];
rz(2.9002088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55610181) q[0];
sx q[0];
rz(-1.9560768) q[0];
sx q[0];
rz(0.38498621) q[0];
x q[1];
rz(-0.74626211) q[2];
sx q[2];
rz(-2.4733739) q[2];
sx q[2];
rz(1.2506739) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1387716) q[1];
sx q[1];
rz(-1.4797987) q[1];
sx q[1];
rz(2.2343193) q[1];
rz(-pi) q[2];
rz(-1.4199735) q[3];
sx q[3];
rz(-1.5966468) q[3];
sx q[3];
rz(-2.6005656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0138578) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(1.1068809) q[2];
rz(-1.7539304) q[3];
sx q[3];
rz(-0.51968402) q[3];
sx q[3];
rz(-1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7805507) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(-2.4011491) q[0];
rz(-0.45117798) q[1];
sx q[1];
rz(-1.2606882) q[1];
sx q[1];
rz(2.0887451) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8448062) q[0];
sx q[0];
rz(-1.0264945) q[0];
sx q[0];
rz(1.143572) q[0];
rz(-2.8430976) q[2];
sx q[2];
rz(-0.91909354) q[2];
sx q[2];
rz(3.0534844) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4098674) q[1];
sx q[1];
rz(-1.6907398) q[1];
sx q[1];
rz(-0.15323318) q[1];
rz(-3.1316109) q[3];
sx q[3];
rz(-1.1407033) q[3];
sx q[3];
rz(0.86356589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8911002) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(0.54692522) q[2];
rz(-2.8524103) q[3];
sx q[3];
rz(-2.6047891) q[3];
sx q[3];
rz(-3.1291936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96994394) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(-1.7893715) q[0];
rz(-0.2098473) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(-0.23637493) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5445697) q[0];
sx q[0];
rz(-1.835808) q[0];
sx q[0];
rz(0.17770627) q[0];
rz(0.9985853) q[2];
sx q[2];
rz(-2.5479655) q[2];
sx q[2];
rz(1.8400536) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.826556) q[1];
sx q[1];
rz(-2.9737925) q[1];
sx q[1];
rz(1.5312974) q[1];
rz(-pi) q[2];
rz(-1.4218016) q[3];
sx q[3];
rz(-1.0862724) q[3];
sx q[3];
rz(-1.5414343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9005047) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(-0.55348712) q[2];
rz(-0.91529804) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(2.4966911) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6699162) q[0];
sx q[0];
rz(-1.8138509) q[0];
sx q[0];
rz(-0.70415235) q[0];
rz(-1.0559121) q[1];
sx q[1];
rz(-2.6864955) q[1];
sx q[1];
rz(0.59590894) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51673698) q[0];
sx q[0];
rz(-1.3552109) q[0];
sx q[0];
rz(-0.09341021) q[0];
rz(-pi) q[1];
x q[1];
rz(1.709183) q[2];
sx q[2];
rz(-2.2397537) q[2];
sx q[2];
rz(1.5630747) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0101498) q[1];
sx q[1];
rz(-1.7556659) q[1];
sx q[1];
rz(-0.4106945) q[1];
rz(-pi) q[2];
rz(-0.44645198) q[3];
sx q[3];
rz(-2.7342396) q[3];
sx q[3];
rz(-0.27094597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6490877) q[2];
sx q[2];
rz(-2.0272144) q[2];
sx q[2];
rz(0.55111432) q[2];
rz(0.20714949) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15774396) q[0];
sx q[0];
rz(-0.85726964) q[0];
sx q[0];
rz(0.95170784) q[0];
rz(2.6668008) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(-2.8463083) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5135358) q[0];
sx q[0];
rz(-2.7211468) q[0];
sx q[0];
rz(1.9168617) q[0];
rz(1.7252543) q[2];
sx q[2];
rz(-1.2631577) q[2];
sx q[2];
rz(-0.49738202) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5438248) q[1];
sx q[1];
rz(-0.73111594) q[1];
sx q[1];
rz(2.3901229) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6867562) q[3];
sx q[3];
rz(-0.51350683) q[3];
sx q[3];
rz(-2.6503369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7513912) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(0.93377101) q[2];
rz(1.4592524) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(1.6850083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075832531) q[0];
sx q[0];
rz(-1.9969143) q[0];
sx q[0];
rz(-2.899535) q[0];
rz(-2.4767955) q[1];
sx q[1];
rz(-1.2882065) q[1];
sx q[1];
rz(2.738293) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6751854) q[0];
sx q[0];
rz(-1.9403606) q[0];
sx q[0];
rz(-2.7778366) q[0];
x q[1];
rz(0.34801872) q[2];
sx q[2];
rz(-0.70901477) q[2];
sx q[2];
rz(-1.4460756) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1937716) q[1];
sx q[1];
rz(-1.7809476) q[1];
sx q[1];
rz(0.0019046849) q[1];
rz(1.7014916) q[3];
sx q[3];
rz(-2.9500467) q[3];
sx q[3];
rz(0.46476118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2287801) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(-2.7071803) q[2];
rz(-2.1408391) q[3];
sx q[3];
rz(-2.775511) q[3];
sx q[3];
rz(0.51030695) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0077165724) q[0];
sx q[0];
rz(-3.1354597) q[0];
sx q[0];
rz(-0.49466053) q[0];
rz(-1.6330632) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(-2.5411434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4279815) q[0];
sx q[0];
rz(-2.5117154) q[0];
sx q[0];
rz(2.4226818) q[0];
rz(-1.5254283) q[2];
sx q[2];
rz(-0.52549911) q[2];
sx q[2];
rz(0.85277992) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3214896) q[1];
sx q[1];
rz(-0.21807018) q[1];
sx q[1];
rz(0.58184187) q[1];
rz(-2.5478701) q[3];
sx q[3];
rz(-0.36168081) q[3];
sx q[3];
rz(2.7152674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1373458) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(-0.11631575) q[2];
rz(0.4256734) q[3];
sx q[3];
rz(-0.95723546) q[3];
sx q[3];
rz(11*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.15329926) q[0];
sx q[0];
rz(-0.17803742) q[0];
sx q[0];
rz(-1.4784038) q[0];
rz(-0.93961811) q[1];
sx q[1];
rz(-1.3213108) q[1];
sx q[1];
rz(0.41752648) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5696213) q[0];
sx q[0];
rz(-2.6082391) q[0];
sx q[0];
rz(1.4825975) q[0];
x q[1];
rz(-1.19403) q[2];
sx q[2];
rz(-0.28290877) q[2];
sx q[2];
rz(-2.6434968) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9081887) q[1];
sx q[1];
rz(-1.7618124) q[1];
sx q[1];
rz(1.7407655) q[1];
rz(0.19374356) q[3];
sx q[3];
rz(-1.994641) q[3];
sx q[3];
rz(1.2731874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6614762) q[2];
sx q[2];
rz(-2.5051703) q[2];
sx q[2];
rz(2.647906) q[2];
rz(2.4168329) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(-2.8216968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17393728) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(2.4560112) q[0];
rz(2.8441692) q[1];
sx q[1];
rz(-2.90459) q[1];
sx q[1];
rz(1.1313653) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81998527) q[0];
sx q[0];
rz(-1.2258343) q[0];
sx q[0];
rz(-2.5489775) q[0];
x q[1];
rz(2.4234424) q[2];
sx q[2];
rz(-0.8469204) q[2];
sx q[2];
rz(1.7921599) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2724534) q[1];
sx q[1];
rz(-1.1932045) q[1];
sx q[1];
rz(-1.8121522) q[1];
rz(-0.64254909) q[3];
sx q[3];
rz(-0.71392871) q[3];
sx q[3];
rz(-0.16845265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.82548213) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(0.70739174) q[2];
rz(-2.3317544) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(0.18856089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903044) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(0.21223016) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(2.390776) q[2];
sx q[2];
rz(-2.3859947) q[2];
sx q[2];
rz(1.044556) q[2];
rz(0.57701941) q[3];
sx q[3];
rz(-0.95303017) q[3];
sx q[3];
rz(-0.87458761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

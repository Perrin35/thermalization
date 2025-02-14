OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8706239) q[0];
sx q[0];
rz(-2.5854817) q[0];
sx q[0];
rz(-2.1882353) q[0];
rz(0.019729992) q[1];
sx q[1];
rz(-2.1058197) q[1];
sx q[1];
rz(0.9486202) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9069288) q[0];
sx q[0];
rz(-1.0415823) q[0];
sx q[0];
rz(-2.0509023) q[0];
rz(0.011587338) q[2];
sx q[2];
rz(-0.23749781) q[2];
sx q[2];
rz(1.604014) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.967077) q[1];
sx q[1];
rz(-0.83454718) q[1];
sx q[1];
rz(0.59474919) q[1];
rz(2.0998276) q[3];
sx q[3];
rz(-0.62058016) q[3];
sx q[3];
rz(-2.4247501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3775776) q[2];
sx q[2];
rz(-3.0674051) q[2];
sx q[2];
rz(-2.3589676) q[2];
rz(3.022656) q[3];
sx q[3];
rz(-1.0340034) q[3];
sx q[3];
rz(-2.9940166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9453732) q[0];
sx q[0];
rz(-1.1213028) q[0];
sx q[0];
rz(0.25783208) q[0];
rz(-0.091015426) q[1];
sx q[1];
rz(-1.0992522) q[1];
sx q[1];
rz(-1.6450504) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5071867) q[0];
sx q[0];
rz(-1.0189462) q[0];
sx q[0];
rz(3.0613012) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7690954) q[2];
sx q[2];
rz(-1.2476693) q[2];
sx q[2];
rz(-1.966764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1196668) q[1];
sx q[1];
rz(-0.41321427) q[1];
sx q[1];
rz(-1.6456804) q[1];
rz(-pi) q[2];
rz(2.4739059) q[3];
sx q[3];
rz(-1.8996618) q[3];
sx q[3];
rz(-1.937695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1220793) q[2];
sx q[2];
rz(-1.8668819) q[2];
sx q[2];
rz(-2.7369734) q[2];
rz(2.1520065) q[3];
sx q[3];
rz(-1.8414958) q[3];
sx q[3];
rz(-2.5185481) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0290381) q[0];
sx q[0];
rz(-0.14415388) q[0];
sx q[0];
rz(2.3032904) q[0];
rz(-0.84838947) q[1];
sx q[1];
rz(-1.219607) q[1];
sx q[1];
rz(-2.1489876) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2780571) q[0];
sx q[0];
rz(-0.73672265) q[0];
sx q[0];
rz(-0.20918734) q[0];
rz(-pi) q[1];
rz(2.9783932) q[2];
sx q[2];
rz(-0.33849785) q[2];
sx q[2];
rz(-1.4783354) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4731208) q[1];
sx q[1];
rz(-1.8278012) q[1];
sx q[1];
rz(-2.7196252) q[1];
x q[2];
rz(-2.8512673) q[3];
sx q[3];
rz(-1.6772406) q[3];
sx q[3];
rz(0.31600472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8695716) q[2];
sx q[2];
rz(-1.1676936) q[2];
sx q[2];
rz(-2.0167548) q[2];
rz(0.62075067) q[3];
sx q[3];
rz(-0.97095942) q[3];
sx q[3];
rz(3.0264405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89897412) q[0];
sx q[0];
rz(-1.1599351) q[0];
sx q[0];
rz(0.17380357) q[0];
rz(0.68570343) q[1];
sx q[1];
rz(-1.655429) q[1];
sx q[1];
rz(2.3033843) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4259151) q[0];
sx q[0];
rz(-1.1155778) q[0];
sx q[0];
rz(-2.0865738) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57212215) q[2];
sx q[2];
rz(-1.1231218) q[2];
sx q[2];
rz(-1.3037701) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.77828171) q[1];
sx q[1];
rz(-1.4290591) q[1];
sx q[1];
rz(-2.4635386) q[1];
rz(-pi) q[2];
rz(1.6082786) q[3];
sx q[3];
rz(-2.0333383) q[3];
sx q[3];
rz(-1.9328062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.59914261) q[2];
sx q[2];
rz(-1.4582381) q[2];
sx q[2];
rz(2.7899) q[2];
rz(-1.1951949) q[3];
sx q[3];
rz(-1.9545133) q[3];
sx q[3];
rz(3.1174507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.675932) q[0];
sx q[0];
rz(-2.6213578) q[0];
sx q[0];
rz(-0.35292536) q[0];
rz(-0.30883166) q[1];
sx q[1];
rz(-1.0363204) q[1];
sx q[1];
rz(0.028506361) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7362979) q[0];
sx q[0];
rz(-2.5319244) q[0];
sx q[0];
rz(-1.2627312) q[0];
rz(-pi) q[1];
rz(-1.1813004) q[2];
sx q[2];
rz(-1.2151698) q[2];
sx q[2];
rz(-2.2782922) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9069017) q[1];
sx q[1];
rz(-2.3447835) q[1];
sx q[1];
rz(-1.2654348) q[1];
x q[2];
rz(-2.8978917) q[3];
sx q[3];
rz(-1.9944046) q[3];
sx q[3];
rz(0.085368644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9286524) q[2];
sx q[2];
rz(-3.0781367) q[2];
sx q[2];
rz(2.1707936) q[2];
rz(-2.094723) q[3];
sx q[3];
rz(-1.4528843) q[3];
sx q[3];
rz(0.48986062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8397119) q[0];
sx q[0];
rz(-1.7827001) q[0];
sx q[0];
rz(-3.0506328) q[0];
rz(1.2184527) q[1];
sx q[1];
rz(-0.33088845) q[1];
sx q[1];
rz(2.5023696) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5175745) q[0];
sx q[0];
rz(-1.7057452) q[0];
sx q[0];
rz(-2.8132416) q[0];
x q[1];
rz(-0.12219723) q[2];
sx q[2];
rz(-1.4336042) q[2];
sx q[2];
rz(2.4047763) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51992937) q[1];
sx q[1];
rz(-2.2022793) q[1];
sx q[1];
rz(0.91543053) q[1];
x q[2];
rz(3.1396041) q[3];
sx q[3];
rz(-1.244215) q[3];
sx q[3];
rz(-1.1359362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.062139221) q[2];
sx q[2];
rz(-2.1663351) q[2];
sx q[2];
rz(-0.47387588) q[2];
rz(-1.3300995) q[3];
sx q[3];
rz(-0.88089839) q[3];
sx q[3];
rz(0.37548319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5652931) q[0];
sx q[0];
rz(-1.0813035) q[0];
sx q[0];
rz(2.9190049) q[0];
rz(1.0635771) q[1];
sx q[1];
rz(-1.5779326) q[1];
sx q[1];
rz(1.9482013) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28119606) q[0];
sx q[0];
rz(-1.4544857) q[0];
sx q[0];
rz(0.75006811) q[0];
rz(-pi) q[1];
rz(1.7980099) q[2];
sx q[2];
rz(-1.4526723) q[2];
sx q[2];
rz(0.83641499) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4970253) q[1];
sx q[1];
rz(-0.46176061) q[1];
sx q[1];
rz(0.30739947) q[1];
rz(-pi) q[2];
rz(-0.16321793) q[3];
sx q[3];
rz(-1.8344518) q[3];
sx q[3];
rz(1.4782983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0003164) q[2];
sx q[2];
rz(-2.2681984) q[2];
sx q[2];
rz(0.6401965) q[2];
rz(-0.021942465) q[3];
sx q[3];
rz(-1.4098189) q[3];
sx q[3];
rz(-0.35023165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.491965) q[0];
sx q[0];
rz(-1.3977298) q[0];
sx q[0];
rz(-2.6212027) q[0];
rz(-0.87977663) q[1];
sx q[1];
rz(-1.9089411) q[1];
sx q[1];
rz(2.2487776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3443106) q[0];
sx q[0];
rz(-1.6160674) q[0];
sx q[0];
rz(3.0540953) q[0];
rz(-0.73142902) q[2];
sx q[2];
rz(-2.1831704) q[2];
sx q[2];
rz(0.8548255) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0860111) q[1];
sx q[1];
rz(-2.308929) q[1];
sx q[1];
rz(-2.7559381) q[1];
rz(1.1475943) q[3];
sx q[3];
rz(-2.0484006) q[3];
sx q[3];
rz(1.783361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2616547) q[2];
sx q[2];
rz(-2.0377906) q[2];
sx q[2];
rz(0.086816303) q[2];
rz(0.9451198) q[3];
sx q[3];
rz(-0.34543959) q[3];
sx q[3];
rz(-2.0115578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18411186) q[0];
sx q[0];
rz(-0.090395398) q[0];
sx q[0];
rz(-3.0775253) q[0];
rz(2.0059026) q[1];
sx q[1];
rz(-1.4402729) q[1];
sx q[1];
rz(0.51220977) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24146809) q[0];
sx q[0];
rz(-1.1607803) q[0];
sx q[0];
rz(-1.926169) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86821235) q[2];
sx q[2];
rz(-1.8382475) q[2];
sx q[2];
rz(-2.3600459) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.54359084) q[1];
sx q[1];
rz(-0.27328047) q[1];
sx q[1];
rz(-0.32229801) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.104761) q[3];
sx q[3];
rz(-2.7324711) q[3];
sx q[3];
rz(-0.44884071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8168489) q[2];
sx q[2];
rz(-2.4090359) q[2];
sx q[2];
rz(1.0653488) q[2];
rz(-1.6522853) q[3];
sx q[3];
rz(-0.95732147) q[3];
sx q[3];
rz(3.0491507) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5866933) q[0];
sx q[0];
rz(-0.41189343) q[0];
sx q[0];
rz(-0.33083415) q[0];
rz(1.2384442) q[1];
sx q[1];
rz(-0.60791433) q[1];
sx q[1];
rz(0.27473658) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8023668) q[0];
sx q[0];
rz(-0.5623835) q[0];
sx q[0];
rz(2.2455162) q[0];
rz(-pi) q[1];
rz(-2.8682235) q[2];
sx q[2];
rz(-1.5393352) q[2];
sx q[2];
rz(0.69845573) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3420136) q[1];
sx q[1];
rz(-1.7106317) q[1];
sx q[1];
rz(-2.0796156) q[1];
rz(-3.069271) q[3];
sx q[3];
rz(-2.157271) q[3];
sx q[3];
rz(0.14432913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1642509) q[2];
sx q[2];
rz(-2.3873316) q[2];
sx q[2];
rz(1.784262) q[2];
rz(-2.6409798) q[3];
sx q[3];
rz(-1.5784135) q[3];
sx q[3];
rz(-1.4969426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5951344) q[0];
sx q[0];
rz(-1.559579) q[0];
sx q[0];
rz(2.5478242) q[0];
rz(3.0733227) q[1];
sx q[1];
rz(-1.2208114) q[1];
sx q[1];
rz(0.023718871) q[1];
rz(2.7082607) q[2];
sx q[2];
rz(-2.7637252) q[2];
sx q[2];
rz(-2.6890474) q[2];
rz(-2.3292062) q[3];
sx q[3];
rz(-0.25573041) q[3];
sx q[3];
rz(-2.0144016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

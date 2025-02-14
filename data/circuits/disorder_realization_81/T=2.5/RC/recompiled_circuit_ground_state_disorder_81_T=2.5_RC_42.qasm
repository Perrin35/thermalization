OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9658907) q[0];
sx q[0];
rz(-1.8718636) q[0];
sx q[0];
rz(1.1748535) q[0];
rz(-1.084561) q[1];
sx q[1];
rz(-1.4718055) q[1];
sx q[1];
rz(-2.5785799) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5064237) q[0];
sx q[0];
rz(-2.9384268) q[0];
sx q[0];
rz(0.26655339) q[0];
x q[1];
rz(-2.2038648) q[2];
sx q[2];
rz(-1.9652624) q[2];
sx q[2];
rz(-1.3317837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6463794) q[1];
sx q[1];
rz(-0.36839147) q[1];
sx q[1];
rz(-2.6839593) q[1];
rz(0.75213856) q[3];
sx q[3];
rz(-2.6400551) q[3];
sx q[3];
rz(0.88830321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6072443) q[2];
sx q[2];
rz(-1.4899985) q[2];
sx q[2];
rz(-1.472817) q[2];
rz(-1.8913174) q[3];
sx q[3];
rz(-1.6142802) q[3];
sx q[3];
rz(-2.1691624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8984574) q[0];
sx q[0];
rz(-0.055483015) q[0];
sx q[0];
rz(-2.6836416) q[0];
rz(-0.10174879) q[1];
sx q[1];
rz(-1.190217) q[1];
sx q[1];
rz(-0.32624689) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75263035) q[0];
sx q[0];
rz(-1.9174308) q[0];
sx q[0];
rz(-2.017157) q[0];
rz(-pi) q[1];
rz(-1.7918705) q[2];
sx q[2];
rz(-1.862163) q[2];
sx q[2];
rz(-1.2779026) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.65349039) q[1];
sx q[1];
rz(-2.2970169) q[1];
sx q[1];
rz(2.2745871) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2539034) q[3];
sx q[3];
rz(-2.2884048) q[3];
sx q[3];
rz(-0.16555351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.132167) q[2];
sx q[2];
rz(-2.1475809) q[2];
sx q[2];
rz(-1.2635292) q[2];
rz(3.1215014) q[3];
sx q[3];
rz(-0.76954904) q[3];
sx q[3];
rz(1.2300389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024071368) q[0];
sx q[0];
rz(-0.063952359) q[0];
sx q[0];
rz(2.5653895) q[0];
rz(1.6974712) q[1];
sx q[1];
rz(-1.8918119) q[1];
sx q[1];
rz(-1.261927) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44718868) q[0];
sx q[0];
rz(-1.2617636) q[0];
sx q[0];
rz(0.033559994) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6143965) q[2];
sx q[2];
rz(-1.1644852) q[2];
sx q[2];
rz(0.57002178) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0351959) q[1];
sx q[1];
rz(-0.96265618) q[1];
sx q[1];
rz(1.870283) q[1];
x q[2];
rz(-2.4275196) q[3];
sx q[3];
rz(-2.0104109) q[3];
sx q[3];
rz(-0.022865828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55804092) q[2];
sx q[2];
rz(-2.5463107) q[2];
sx q[2];
rz(-2.9884647) q[2];
rz(-2.2554452) q[3];
sx q[3];
rz(-1.7821507) q[3];
sx q[3];
rz(1.9396293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8837638) q[0];
sx q[0];
rz(-1.184329) q[0];
sx q[0];
rz(-0.15876874) q[0];
rz(1.0348882) q[1];
sx q[1];
rz(-2.1885469) q[1];
sx q[1];
rz(-0.75611702) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0466111) q[0];
sx q[0];
rz(-1.5205407) q[0];
sx q[0];
rz(0.057970164) q[0];
rz(1.6768084) q[2];
sx q[2];
rz(-2.8869625) q[2];
sx q[2];
rz(2.1531596) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8146921) q[1];
sx q[1];
rz(-1.3291825) q[1];
sx q[1];
rz(1.7820226) q[1];
rz(-pi) q[2];
rz(0.71847312) q[3];
sx q[3];
rz(-2.8333896) q[3];
sx q[3];
rz(1.1402477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54255828) q[2];
sx q[2];
rz(-0.90347806) q[2];
sx q[2];
rz(-0.8503882) q[2];
rz(-1.4706069) q[3];
sx q[3];
rz(-1.8149523) q[3];
sx q[3];
rz(0.82408041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65664148) q[0];
sx q[0];
rz(-1.4018207) q[0];
sx q[0];
rz(-2.9630419) q[0];
rz(1.8932331) q[1];
sx q[1];
rz(-1.5310042) q[1];
sx q[1];
rz(-0.53370968) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2399143) q[0];
sx q[0];
rz(-1.9199492) q[0];
sx q[0];
rz(-1.7656183) q[0];
rz(2.1458964) q[2];
sx q[2];
rz(-0.786869) q[2];
sx q[2];
rz(-2.5520419) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2055997) q[1];
sx q[1];
rz(-2.471711) q[1];
sx q[1];
rz(-1.4861466) q[1];
rz(-pi) q[2];
rz(-1.0664135) q[3];
sx q[3];
rz(-1.4850067) q[3];
sx q[3];
rz(-2.498621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8278213) q[2];
sx q[2];
rz(-1.2510108) q[2];
sx q[2];
rz(-1.0531462) q[2];
rz(2.961212) q[3];
sx q[3];
rz(-0.80612055) q[3];
sx q[3];
rz(-2.7789796) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5798222) q[0];
sx q[0];
rz(-1.0012015) q[0];
sx q[0];
rz(-1.9129821) q[0];
rz(2.6749532) q[1];
sx q[1];
rz(-1.2184315) q[1];
sx q[1];
rz(-0.12900464) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056696293) q[0];
sx q[0];
rz(-1.7040623) q[0];
sx q[0];
rz(-0.10198089) q[0];
rz(-0.1618218) q[2];
sx q[2];
rz(-0.84748778) q[2];
sx q[2];
rz(-0.79297334) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1946757) q[1];
sx q[1];
rz(-1.2763775) q[1];
sx q[1];
rz(0.63398449) q[1];
rz(-0.13997649) q[3];
sx q[3];
rz(-1.3050606) q[3];
sx q[3];
rz(-1.1418726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0064319) q[2];
sx q[2];
rz(-0.60144037) q[2];
sx q[2];
rz(0.50430164) q[2];
rz(-1.0935316) q[3];
sx q[3];
rz(-1.3239599) q[3];
sx q[3];
rz(-2.1849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.6495551) q[0];
sx q[0];
rz(-1.2832337) q[0];
sx q[0];
rz(-1.5564224) q[0];
rz(-0.57836142) q[1];
sx q[1];
rz(-1.882694) q[1];
sx q[1];
rz(-3.0392821) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71891448) q[0];
sx q[0];
rz(-2.3862729) q[0];
sx q[0];
rz(-0.97277555) q[0];
rz(-pi) q[1];
rz(-0.39686508) q[2];
sx q[2];
rz(-1.3584653) q[2];
sx q[2];
rz(-3.120829) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7150517) q[1];
sx q[1];
rz(-0.3449966) q[1];
sx q[1];
rz(2.662602) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2874114) q[3];
sx q[3];
rz(-2.183217) q[3];
sx q[3];
rz(1.8068594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9461225) q[2];
sx q[2];
rz(-1.2637694) q[2];
sx q[2];
rz(-0.61402399) q[2];
rz(-2.2139003) q[3];
sx q[3];
rz(-1.5988908) q[3];
sx q[3];
rz(-0.6774925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5880661) q[0];
sx q[0];
rz(-0.29356846) q[0];
sx q[0];
rz(1.9422096) q[0];
rz(-1.9623494) q[1];
sx q[1];
rz(-1.0284547) q[1];
sx q[1];
rz(-0.95338043) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6567568) q[0];
sx q[0];
rz(-1.0547045) q[0];
sx q[0];
rz(2.1950095) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10886635) q[2];
sx q[2];
rz(-2.1061387) q[2];
sx q[2];
rz(-1.1606729) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.39722193) q[1];
sx q[1];
rz(-2.3030641) q[1];
sx q[1];
rz(2.9930755) q[1];
rz(0.48170959) q[3];
sx q[3];
rz(-1.6407971) q[3];
sx q[3];
rz(-0.21177975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.392268) q[2];
sx q[2];
rz(-0.12568036) q[2];
sx q[2];
rz(-2.9607062) q[2];
rz(1.1130029) q[3];
sx q[3];
rz(-1.0180611) q[3];
sx q[3];
rz(-0.38016144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0522633) q[0];
sx q[0];
rz(-2.2398529) q[0];
sx q[0];
rz(0.7641291) q[0];
rz(-2.1185421) q[1];
sx q[1];
rz(-1.5307531) q[1];
sx q[1];
rz(1.6463564) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29158026) q[0];
sx q[0];
rz(-0.33667013) q[0];
sx q[0];
rz(0.48179896) q[0];
rz(-1.1029215) q[2];
sx q[2];
rz(-2.4870424) q[2];
sx q[2];
rz(-2.2783115) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3493986) q[1];
sx q[1];
rz(-2.7490834) q[1];
sx q[1];
rz(-2.9540707) q[1];
rz(1.1782618) q[3];
sx q[3];
rz(-1.1589253) q[3];
sx q[3];
rz(-0.28366006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9493147) q[2];
sx q[2];
rz(-2.1239026) q[2];
sx q[2];
rz(-0.44753543) q[2];
rz(-0.13628515) q[3];
sx q[3];
rz(-0.49301967) q[3];
sx q[3];
rz(2.5435257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78551453) q[0];
sx q[0];
rz(-2.1165753) q[0];
sx q[0];
rz(0.4726952) q[0];
rz(0.39974943) q[1];
sx q[1];
rz(-2.4374297) q[1];
sx q[1];
rz(-0.8185111) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34270014) q[0];
sx q[0];
rz(-1.2620942) q[0];
sx q[0];
rz(-1.3837455) q[0];
rz(-pi) q[1];
rz(0.50869663) q[2];
sx q[2];
rz(-0.74328586) q[2];
sx q[2];
rz(0.41451472) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5732729) q[1];
sx q[1];
rz(-1.5675873) q[1];
sx q[1];
rz(-1.248097) q[1];
x q[2];
rz(-2.6916041) q[3];
sx q[3];
rz(-0.38211029) q[3];
sx q[3];
rz(1.1799605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58273903) q[2];
sx q[2];
rz(-1.0820729) q[2];
sx q[2];
rz(-1.8633128) q[2];
rz(-1.1996972) q[3];
sx q[3];
rz(-1.4408828) q[3];
sx q[3];
rz(0.27967683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77547913) q[0];
sx q[0];
rz(-1.4237325) q[0];
sx q[0];
rz(1.7970418) q[0];
rz(1.7300425) q[1];
sx q[1];
rz(-0.81137864) q[1];
sx q[1];
rz(2.484533) q[1];
rz(-3.0006164) q[2];
sx q[2];
rz(-1.5700414) q[2];
sx q[2];
rz(1.0667917) q[2];
rz(-1.9947407) q[3];
sx q[3];
rz(-0.70953555) q[3];
sx q[3];
rz(-2.4279688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

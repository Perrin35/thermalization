OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2312343) q[0];
sx q[0];
rz(5.4141747) q[0];
sx q[0];
rz(10.5095) q[0];
rz(-1.1551069) q[1];
sx q[1];
rz(-0.81973633) q[1];
sx q[1];
rz(-0.91135946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87093519) q[0];
sx q[0];
rz(-2.3349999) q[0];
sx q[0];
rz(2.8034504) q[0];
rz(-pi) q[1];
rz(-0.16884825) q[2];
sx q[2];
rz(-0.71879866) q[2];
sx q[2];
rz(-2.0618677) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5973917) q[1];
sx q[1];
rz(-2.4185838) q[1];
sx q[1];
rz(0.71031481) q[1];
rz(-0.064685589) q[3];
sx q[3];
rz(-2.0035335) q[3];
sx q[3];
rz(0.81165867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.34065166) q[2];
sx q[2];
rz(-2.030535) q[2];
sx q[2];
rz(-0.079455376) q[2];
rz(0.60845145) q[3];
sx q[3];
rz(-1.2234917) q[3];
sx q[3];
rz(0.89103812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65421739) q[0];
sx q[0];
rz(-2.2779164) q[0];
sx q[0];
rz(-1.7864216) q[0];
rz(1.0379418) q[1];
sx q[1];
rz(-0.67960056) q[1];
sx q[1];
rz(-0.80744809) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6896967) q[0];
sx q[0];
rz(-0.97366714) q[0];
sx q[0];
rz(0.35291617) q[0];
rz(-pi) q[1];
rz(2.7620671) q[2];
sx q[2];
rz(-2.1132601) q[2];
sx q[2];
rz(1.2778953) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.24570736) q[1];
sx q[1];
rz(-1.504532) q[1];
sx q[1];
rz(0.23991983) q[1];
x q[2];
rz(-3.03504) q[3];
sx q[3];
rz(-1.7037183) q[3];
sx q[3];
rz(-0.35655856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9281533) q[2];
sx q[2];
rz(-1.5779053) q[2];
sx q[2];
rz(1.6820924) q[2];
rz(2.6835119) q[3];
sx q[3];
rz(-2.5638678) q[3];
sx q[3];
rz(-2.1817082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9816575) q[0];
sx q[0];
rz(-2.0913048) q[0];
sx q[0];
rz(-2.4666069) q[0];
rz(0.58810294) q[1];
sx q[1];
rz(-2.2876078) q[1];
sx q[1];
rz(-1.3311707) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0818401) q[0];
sx q[0];
rz(-0.65048238) q[0];
sx q[0];
rz(2.4516134) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0841062) q[2];
sx q[2];
rz(-2.170544) q[2];
sx q[2];
rz(-2.5059003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.59440295) q[1];
sx q[1];
rz(-1.8246027) q[1];
sx q[1];
rz(-2.3188792) q[1];
x q[2];
rz(2.0698333) q[3];
sx q[3];
rz(-1.3924358) q[3];
sx q[3];
rz(-0.87770977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5171234) q[2];
sx q[2];
rz(-0.34999592) q[2];
sx q[2];
rz(0.2207174) q[2];
rz(-2.3366426) q[3];
sx q[3];
rz(-1.5419818) q[3];
sx q[3];
rz(1.8055003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6014366) q[0];
sx q[0];
rz(-2.8995081) q[0];
sx q[0];
rz(-2.5228187) q[0];
rz(-3.0586808) q[1];
sx q[1];
rz(-0.34268788) q[1];
sx q[1];
rz(1.6212911) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5280209) q[0];
sx q[0];
rz(-0.53410406) q[0];
sx q[0];
rz(-1.3811124) q[0];
rz(-0.35334584) q[2];
sx q[2];
rz(-2.5328703) q[2];
sx q[2];
rz(2.6578987) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1136696) q[1];
sx q[1];
rz(-1.4884559) q[1];
sx q[1];
rz(-2.0827977) q[1];
rz(2.7079775) q[3];
sx q[3];
rz(-0.40900298) q[3];
sx q[3];
rz(-0.65164372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.80359047) q[2];
sx q[2];
rz(-0.10927304) q[2];
sx q[2];
rz(0.14536157) q[2];
rz(-2.8694782) q[3];
sx q[3];
rz(-1.4067255) q[3];
sx q[3];
rz(1.9947778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014932545) q[0];
sx q[0];
rz(-1.3260051) q[0];
sx q[0];
rz(-3.0354011) q[0];
rz(-2.1821678) q[1];
sx q[1];
rz(-0.54490772) q[1];
sx q[1];
rz(-0.19439654) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17592573) q[0];
sx q[0];
rz(-2.3011294) q[0];
sx q[0];
rz(-1.6863807) q[0];
rz(-pi) q[1];
rz(1.6129877) q[2];
sx q[2];
rz(-0.11494177) q[2];
sx q[2];
rz(1.3745705) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5867956) q[1];
sx q[1];
rz(-1.6117084) q[1];
sx q[1];
rz(0.58920963) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3931469) q[3];
sx q[3];
rz(-0.41209778) q[3];
sx q[3];
rz(0.92890152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83547366) q[2];
sx q[2];
rz(-1.4543616) q[2];
sx q[2];
rz(2.5731738) q[2];
rz(-3.0889555) q[3];
sx q[3];
rz(-1.6391552) q[3];
sx q[3];
rz(3.0047825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0282054) q[0];
sx q[0];
rz(-0.55332342) q[0];
sx q[0];
rz(-2.9521039) q[0];
rz(1.0264617) q[1];
sx q[1];
rz(-0.84954134) q[1];
sx q[1];
rz(2.386327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7758691) q[0];
sx q[0];
rz(-0.9277753) q[0];
sx q[0];
rz(-0.68277208) q[0];
rz(1.009047) q[2];
sx q[2];
rz(-2.0239186) q[2];
sx q[2];
rz(0.31227408) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64769563) q[1];
sx q[1];
rz(-1.4129479) q[1];
sx q[1];
rz(1.3998019) q[1];
x q[2];
rz(-1.9368265) q[3];
sx q[3];
rz(-1.2377973) q[3];
sx q[3];
rz(-0.35375094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.75817627) q[2];
sx q[2];
rz(-1.5864317) q[2];
sx q[2];
rz(-1.7002534) q[2];
rz(1.7656743) q[3];
sx q[3];
rz(-2.223189) q[3];
sx q[3];
rz(2.4413696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43948424) q[0];
sx q[0];
rz(-1.0478042) q[0];
sx q[0];
rz(1.1442319) q[0];
rz(-0.2598091) q[1];
sx q[1];
rz(-2.3016498) q[1];
sx q[1];
rz(-0.98794404) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081962498) q[0];
sx q[0];
rz(-0.10082997) q[0];
sx q[0];
rz(0.16720812) q[0];
x q[1];
rz(-1.1455215) q[2];
sx q[2];
rz(-0.88103154) q[2];
sx q[2];
rz(-0.20343101) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4960031) q[1];
sx q[1];
rz(-0.92964806) q[1];
sx q[1];
rz(1.251613) q[1];
rz(-pi) q[2];
rz(-0.99616237) q[3];
sx q[3];
rz(-1.8271433) q[3];
sx q[3];
rz(-2.4202895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40954956) q[2];
sx q[2];
rz(-1.5012375) q[2];
sx q[2];
rz(-2.5284956) q[2];
rz(0.71550718) q[3];
sx q[3];
rz(-2.3698273) q[3];
sx q[3];
rz(-2.7806921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9382984) q[0];
sx q[0];
rz(-0.84380117) q[0];
sx q[0];
rz(-0.26813689) q[0];
rz(-0.2306436) q[1];
sx q[1];
rz(-1.5985039) q[1];
sx q[1];
rz(-0.014009744) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9833355) q[0];
sx q[0];
rz(-1.7795441) q[0];
sx q[0];
rz(0.76738417) q[0];
x q[1];
rz(2.9815648) q[2];
sx q[2];
rz(-2.0065971) q[2];
sx q[2];
rz(1.7444812) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6790633) q[1];
sx q[1];
rz(-0.53906295) q[1];
sx q[1];
rz(3.0545337) q[1];
rz(-pi) q[2];
rz(-3.1118891) q[3];
sx q[3];
rz(-1.4280768) q[3];
sx q[3];
rz(2.5102455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9739428) q[2];
sx q[2];
rz(-1.3385945) q[2];
sx q[2];
rz(0.30926427) q[2];
rz(-0.15657982) q[3];
sx q[3];
rz(-1.7370109) q[3];
sx q[3];
rz(2.1046765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2579047) q[0];
sx q[0];
rz(-0.64249277) q[0];
sx q[0];
rz(2.8908492) q[0];
rz(-2.5550487) q[1];
sx q[1];
rz(-1.5637014) q[1];
sx q[1];
rz(2.3768545) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4630554) q[0];
sx q[0];
rz(-1.2559811) q[0];
sx q[0];
rz(0.30588715) q[0];
rz(0.7711507) q[2];
sx q[2];
rz(-1.4284819) q[2];
sx q[2];
rz(-2.3754295) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8426399) q[1];
sx q[1];
rz(-3.0184426) q[1];
sx q[1];
rz(-3.0663436) q[1];
x q[2];
rz(-1.6254201) q[3];
sx q[3];
rz(-1.9491073) q[3];
sx q[3];
rz(-0.19669939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5137198) q[2];
sx q[2];
rz(-0.79078117) q[2];
sx q[2];
rz(2.9676843) q[2];
rz(0.74226132) q[3];
sx q[3];
rz(-2.3895013) q[3];
sx q[3];
rz(3.098587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0088418) q[0];
sx q[0];
rz(-0.77021563) q[0];
sx q[0];
rz(-2.8736864) q[0];
rz(1.335089) q[1];
sx q[1];
rz(-0.91064149) q[1];
sx q[1];
rz(1.3722027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5822179) q[0];
sx q[0];
rz(-2.0204883) q[0];
sx q[0];
rz(2.873308) q[0];
rz(-pi) q[1];
rz(1.4698974) q[2];
sx q[2];
rz(-1.865662) q[2];
sx q[2];
rz(-3.1050668) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7160373) q[1];
sx q[1];
rz(-0.73409427) q[1];
sx q[1];
rz(-2.3922763) q[1];
rz(-pi) q[2];
rz(0.31900556) q[3];
sx q[3];
rz(-1.8847163) q[3];
sx q[3];
rz(-0.55578243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.043896) q[2];
sx q[2];
rz(-1.0363657) q[2];
sx q[2];
rz(1.3807266) q[2];
rz(-1.1287639) q[3];
sx q[3];
rz(-2.1916316) q[3];
sx q[3];
rz(-1.5855764) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632521) q[0];
sx q[0];
rz(-1.5338407) q[0];
sx q[0];
rz(-1.1080678) q[0];
rz(-0.68645984) q[1];
sx q[1];
rz(-1.3931128) q[1];
sx q[1];
rz(-1.211094) q[1];
rz(1.8615234) q[2];
sx q[2];
rz(-2.0474993) q[2];
sx q[2];
rz(-1.7119424) q[2];
rz(1.0330647) q[3];
sx q[3];
rz(-1.929639) q[3];
sx q[3];
rz(0.7076984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91035834) q[0];
sx q[0];
rz(-2.2725821) q[0];
sx q[0];
rz(-1.0847217) q[0];
rz(-1.1551069) q[1];
sx q[1];
rz(-0.81973633) q[1];
sx q[1];
rz(-0.91135946) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40075743) q[0];
sx q[0];
rz(-0.82159737) q[0];
sx q[0];
rz(1.9039959) q[0];
x q[1];
rz(-1.7167822) q[2];
sx q[2];
rz(-2.2772191) q[2];
sx q[2];
rz(-1.8391158) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45381308) q[1];
sx q[1];
rz(-1.1247083) q[1];
sx q[1];
rz(0.58961726) q[1];
rz(-pi) q[2];
x q[2];
rz(0.064685589) q[3];
sx q[3];
rz(-1.1380592) q[3];
sx q[3];
rz(0.81165867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.800941) q[2];
sx q[2];
rz(-1.1110577) q[2];
sx q[2];
rz(-0.079455376) q[2];
rz(-2.5331412) q[3];
sx q[3];
rz(-1.918101) q[3];
sx q[3];
rz(-0.89103812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4873753) q[0];
sx q[0];
rz(-2.2779164) q[0];
sx q[0];
rz(-1.7864216) q[0];
rz(1.0379418) q[1];
sx q[1];
rz(-0.67960056) q[1];
sx q[1];
rz(-0.80744809) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8184745) q[0];
sx q[0];
rz(-1.2809296) q[0];
sx q[0];
rz(-2.1978343) q[0];
x q[1];
rz(-1.0196997) q[2];
sx q[2];
rz(-2.4906213) q[2];
sx q[2];
rz(-2.5215182) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.24570736) q[1];
sx q[1];
rz(-1.504532) q[1];
sx q[1];
rz(-0.23991983) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7044675) q[3];
sx q[3];
rz(-1.4651871) q[3];
sx q[3];
rz(-1.228412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9281533) q[2];
sx q[2];
rz(-1.5636874) q[2];
sx q[2];
rz(1.6820924) q[2];
rz(-2.6835119) q[3];
sx q[3];
rz(-0.57772485) q[3];
sx q[3];
rz(0.95988449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.15993519) q[0];
sx q[0];
rz(-1.0502879) q[0];
sx q[0];
rz(-2.4666069) q[0];
rz(0.58810294) q[1];
sx q[1];
rz(-2.2876078) q[1];
sx q[1];
rz(-1.3311707) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0818401) q[0];
sx q[0];
rz(-2.4911103) q[0];
sx q[0];
rz(0.68997927) q[0];
rz(1.0574865) q[2];
sx q[2];
rz(-2.170544) q[2];
sx q[2];
rz(-0.63569234) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2053255) q[1];
sx q[1];
rz(-2.2895439) q[1];
sx q[1];
rz(0.34013744) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0698333) q[3];
sx q[3];
rz(-1.3924358) q[3];
sx q[3];
rz(-0.87770977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5171234) q[2];
sx q[2];
rz(-0.34999592) q[2];
sx q[2];
rz(-0.2207174) q[2];
rz(-2.3366426) q[3];
sx q[3];
rz(-1.5996108) q[3];
sx q[3];
rz(1.3360924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.540156) q[0];
sx q[0];
rz(-0.24208459) q[0];
sx q[0];
rz(2.5228187) q[0];
rz(-0.082911804) q[1];
sx q[1];
rz(-2.7989048) q[1];
sx q[1];
rz(-1.5203016) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6135718) q[0];
sx q[0];
rz(-0.53410406) q[0];
sx q[0];
rz(1.3811124) q[0];
rz(1.8074715) q[2];
sx q[2];
rz(-2.1370721) q[2];
sx q[2];
rz(-0.90619722) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.027923) q[1];
sx q[1];
rz(-1.6531367) q[1];
sx q[1];
rz(2.0827977) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.43361516) q[3];
sx q[3];
rz(-2.7325897) q[3];
sx q[3];
rz(0.65164372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3380022) q[2];
sx q[2];
rz(-0.10927304) q[2];
sx q[2];
rz(-2.9962311) q[2];
rz(-0.27211443) q[3];
sx q[3];
rz(-1.7348671) q[3];
sx q[3];
rz(-1.1468148) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1266601) q[0];
sx q[0];
rz(-1.3260051) q[0];
sx q[0];
rz(3.0354011) q[0];
rz(-0.95942489) q[1];
sx q[1];
rz(-2.5966849) q[1];
sx q[1];
rz(-0.19439654) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9656669) q[0];
sx q[0];
rz(-0.84046326) q[0];
sx q[0];
rz(-1.6863807) q[0];
rz(-pi) q[1];
rz(0.0048695366) q[2];
sx q[2];
rz(-1.4559573) q[2];
sx q[2];
rz(-1.3320992) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5867956) q[1];
sx q[1];
rz(-1.5298843) q[1];
sx q[1];
rz(-2.552383) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9771036) q[3];
sx q[3];
rz(-1.4999564) q[3];
sx q[3];
rz(-2.336647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.306119) q[2];
sx q[2];
rz(-1.6872311) q[2];
sx q[2];
rz(-2.5731738) q[2];
rz(-0.052637188) q[3];
sx q[3];
rz(-1.6391552) q[3];
sx q[3];
rz(-3.0047825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0282054) q[0];
sx q[0];
rz(-0.55332342) q[0];
sx q[0];
rz(0.18948874) q[0];
rz(-2.1151309) q[1];
sx q[1];
rz(-2.2920513) q[1];
sx q[1];
rz(-2.386327) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8403977) q[0];
sx q[0];
rz(-2.2406672) q[0];
sx q[0];
rz(-0.87090454) q[0];
x q[1];
rz(-0.83024518) q[2];
sx q[2];
rz(-2.4355222) q[2];
sx q[2];
rz(-2.490807) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.89596302) q[1];
sx q[1];
rz(-1.7396444) q[1];
sx q[1];
rz(0.16014512) q[1];
rz(0.80249287) q[3];
sx q[3];
rz(-0.48972788) q[3];
sx q[3];
rz(0.51093131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3834164) q[2];
sx q[2];
rz(-1.5551609) q[2];
sx q[2];
rz(-1.4413393) q[2];
rz(-1.3759184) q[3];
sx q[3];
rz(-2.223189) q[3];
sx q[3];
rz(2.4413696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7021084) q[0];
sx q[0];
rz(-1.0478042) q[0];
sx q[0];
rz(1.9973607) q[0];
rz(2.8817835) q[1];
sx q[1];
rz(-0.83994284) q[1];
sx q[1];
rz(-2.1536486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081962498) q[0];
sx q[0];
rz(-0.10082997) q[0];
sx q[0];
rz(0.16720812) q[0];
rz(-pi) q[1];
rz(-1.1455215) q[2];
sx q[2];
rz(-0.88103154) q[2];
sx q[2];
rz(2.9381616) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64558951) q[1];
sx q[1];
rz(-2.2119446) q[1];
sx q[1];
rz(-1.251613) q[1];
x q[2];
rz(-2.8389205) q[3];
sx q[3];
rz(-2.1244123) q[3];
sx q[3];
rz(0.68676567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7320431) q[2];
sx q[2];
rz(-1.6403551) q[2];
sx q[2];
rz(2.5284956) q[2];
rz(2.4260855) q[3];
sx q[3];
rz(-0.77176538) q[3];
sx q[3];
rz(0.36090052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20329423) q[0];
sx q[0];
rz(-2.2977915) q[0];
sx q[0];
rz(-0.26813689) q[0];
rz(0.2306436) q[1];
sx q[1];
rz(-1.5430887) q[1];
sx q[1];
rz(-0.014009744) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5174311) q[0];
sx q[0];
rz(-0.78963477) q[0];
sx q[0];
rz(-0.29615088) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16002782) q[2];
sx q[2];
rz(-1.1349956) q[2];
sx q[2];
rz(-1.3971115) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6790633) q[1];
sx q[1];
rz(-2.6025297) q[1];
sx q[1];
rz(0.087058914) q[1];
rz(-pi) q[2];
rz(3.1118891) q[3];
sx q[3];
rz(-1.4280768) q[3];
sx q[3];
rz(-2.5102455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1676499) q[2];
sx q[2];
rz(-1.8029982) q[2];
sx q[2];
rz(-2.8323284) q[2];
rz(-2.9850128) q[3];
sx q[3];
rz(-1.7370109) q[3];
sx q[3];
rz(1.0369161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2579047) q[0];
sx q[0];
rz(-2.4990999) q[0];
sx q[0];
rz(0.25074348) q[0];
rz(2.5550487) q[1];
sx q[1];
rz(-1.5778912) q[1];
sx q[1];
rz(2.3768545) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79478489) q[0];
sx q[0];
rz(-1.2803923) q[0];
sx q[0];
rz(-1.2417163) q[0];
rz(-2.938835) q[2];
sx q[2];
rz(-2.360095) q[2];
sx q[2];
rz(0.6597214) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7950678) q[1];
sx q[1];
rz(-1.5615614) q[1];
sx q[1];
rz(-3.0187876) q[1];
rz(-pi) q[2];
rz(3.0050817) q[3];
sx q[3];
rz(-0.38204604) q[3];
sx q[3];
rz(-0.34367022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6278729) q[2];
sx q[2];
rz(-0.79078117) q[2];
sx q[2];
rz(2.9676843) q[2];
rz(-2.3993313) q[3];
sx q[3];
rz(-2.3895013) q[3];
sx q[3];
rz(-0.043005634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1327508) q[0];
sx q[0];
rz(-0.77021563) q[0];
sx q[0];
rz(-2.8736864) q[0];
rz(-1.335089) q[1];
sx q[1];
rz(-2.2309512) q[1];
sx q[1];
rz(1.3722027) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2491082) q[0];
sx q[0];
rz(-1.3297446) q[0];
sx q[0];
rz(-2.0349099) q[0];
rz(-pi) q[1];
rz(1.6716953) q[2];
sx q[2];
rz(-1.865662) q[2];
sx q[2];
rz(-0.036525846) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7160373) q[1];
sx q[1];
rz(-0.73409427) q[1];
sx q[1];
rz(2.3922763) q[1];
x q[2];
rz(-0.31900556) q[3];
sx q[3];
rz(-1.8847163) q[3];
sx q[3];
rz(-2.5858102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.043896) q[2];
sx q[2];
rz(-2.1052269) q[2];
sx q[2];
rz(-1.7608661) q[2];
rz(-2.0128287) q[3];
sx q[3];
rz(-2.1916316) q[3];
sx q[3];
rz(-1.5560163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.108528) q[3];
sx q[3];
rz(-1.929639) q[3];
sx q[3];
rz(0.7076984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

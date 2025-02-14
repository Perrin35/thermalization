OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5858894) q[0];
sx q[0];
rz(-1.1893505) q[0];
sx q[0];
rz(1.9777634) q[0];
rz(-3.8883348) q[1];
sx q[1];
rz(0.26489869) q[1];
sx q[1];
rz(12.737552) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25397444) q[0];
sx q[0];
rz(-1.3518392) q[0];
sx q[0];
rz(0.15677126) q[0];
rz(2.2581882) q[2];
sx q[2];
rz(-0.92057487) q[2];
sx q[2];
rz(-1.6028849) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.65003465) q[1];
sx q[1];
rz(-1.2753873) q[1];
sx q[1];
rz(-0.00062382767) q[1];
rz(-pi) q[2];
rz(-1.5762276) q[3];
sx q[3];
rz(-1.5539546) q[3];
sx q[3];
rz(1.668949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.050194) q[2];
sx q[2];
rz(-1.5768496) q[2];
sx q[2];
rz(1.0696627) q[2];
rz(1.4403249) q[3];
sx q[3];
rz(-2.1170719) q[3];
sx q[3];
rz(-1.1721771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9656203) q[0];
sx q[0];
rz(-1.6930641) q[0];
sx q[0];
rz(-2.5536221) q[0];
rz(1.8486456) q[1];
sx q[1];
rz(-0.99423948) q[1];
sx q[1];
rz(-2.9867244) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9706544) q[0];
sx q[0];
rz(-1.0736671) q[0];
sx q[0];
rz(-1.4488104) q[0];
rz(-pi) q[1];
rz(-1.4243519) q[2];
sx q[2];
rz(-1.8093718) q[2];
sx q[2];
rz(-1.6171136) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1701291) q[1];
sx q[1];
rz(-1.2358574) q[1];
sx q[1];
rz(1.0748802) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9294176) q[3];
sx q[3];
rz(-1.7042233) q[3];
sx q[3];
rz(1.9047036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64314848) q[2];
sx q[2];
rz(-2.4109106) q[2];
sx q[2];
rz(-3.0291962) q[2];
rz(-0.82320881) q[3];
sx q[3];
rz(-1.221849) q[3];
sx q[3];
rz(-0.52197758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4560029) q[0];
sx q[0];
rz(-1.0887479) q[0];
sx q[0];
rz(-2.1734557) q[0];
rz(-2.0227506) q[1];
sx q[1];
rz(-1.5348996) q[1];
sx q[1];
rz(-1.8082089) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0007933) q[0];
sx q[0];
rz(-2.5171772) q[0];
sx q[0];
rz(1.2932168) q[0];
rz(-1.1022262) q[2];
sx q[2];
rz(-1.0194687) q[2];
sx q[2];
rz(0.18543772) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0371497) q[1];
sx q[1];
rz(-2.5698476) q[1];
sx q[1];
rz(1.3859831) q[1];
x q[2];
rz(1.795007) q[3];
sx q[3];
rz(-0.36642679) q[3];
sx q[3];
rz(-0.51788766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1916888) q[2];
sx q[2];
rz(-1.0368985) q[2];
sx q[2];
rz(0.26491586) q[2];
rz(2.8252937) q[3];
sx q[3];
rz(-2.8079872) q[3];
sx q[3];
rz(1.3914289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.815627) q[0];
sx q[0];
rz(-2.9209324) q[0];
sx q[0];
rz(-1.6478446) q[0];
rz(1.7296467) q[1];
sx q[1];
rz(-2.1144046) q[1];
sx q[1];
rz(-2.484201) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7058168) q[0];
sx q[0];
rz(-3.0141493) q[0];
sx q[0];
rz(2.5246922) q[0];
rz(-pi) q[1];
rz(2.1543862) q[2];
sx q[2];
rz(-0.87652096) q[2];
sx q[2];
rz(-2.2615446) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3434329) q[1];
sx q[1];
rz(-1.6078336) q[1];
sx q[1];
rz(1.7658893) q[1];
x q[2];
rz(2.3300276) q[3];
sx q[3];
rz(-2.2144284) q[3];
sx q[3];
rz(2.0763458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9221981) q[2];
sx q[2];
rz(-1.9191091) q[2];
sx q[2];
rz(2.738319) q[2];
rz(-1.203631) q[3];
sx q[3];
rz(-1.6324685) q[3];
sx q[3];
rz(0.50886124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68172425) q[0];
sx q[0];
rz(-0.42400703) q[0];
sx q[0];
rz(-0.62393171) q[0];
rz(-0.722018) q[1];
sx q[1];
rz(-1.3474418) q[1];
sx q[1];
rz(-3.0375979) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0777587) q[0];
sx q[0];
rz(-0.88853271) q[0];
sx q[0];
rz(1.3880678) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76461069) q[2];
sx q[2];
rz(-2.6608564) q[2];
sx q[2];
rz(1.3388456) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4482634) q[1];
sx q[1];
rz(-0.86719705) q[1];
sx q[1];
rz(0.87551043) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2448089) q[3];
sx q[3];
rz(-0.95640238) q[3];
sx q[3];
rz(-2.3562285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35088745) q[2];
sx q[2];
rz(-0.71983379) q[2];
sx q[2];
rz(-0.78901115) q[2];
rz(1.2645432) q[3];
sx q[3];
rz(-1.0435373) q[3];
sx q[3];
rz(-2.1544382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3037381) q[0];
sx q[0];
rz(-0.721295) q[0];
sx q[0];
rz(0.20027941) q[0];
rz(-1.6750977) q[1];
sx q[1];
rz(-1.3267582) q[1];
sx q[1];
rz(1.6237578) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87720931) q[0];
sx q[0];
rz(-2.0743999) q[0];
sx q[0];
rz(-1.1491386) q[0];
rz(-pi) q[1];
rz(0.13049651) q[2];
sx q[2];
rz(-2.0708212) q[2];
sx q[2];
rz(-0.44483003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.962783) q[1];
sx q[1];
rz(-0.97985044) q[1];
sx q[1];
rz(1.5750118) q[1];
rz(-2.9041457) q[3];
sx q[3];
rz(-0.81271623) q[3];
sx q[3];
rz(1.0959015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5184861) q[2];
sx q[2];
rz(-0.97272626) q[2];
sx q[2];
rz(2.8002807) q[2];
rz(2.2845279) q[3];
sx q[3];
rz(-1.9144446) q[3];
sx q[3];
rz(0.46522337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.051006) q[0];
sx q[0];
rz(-1.0336646) q[0];
sx q[0];
rz(1.2475913) q[0];
rz(0.11820758) q[1];
sx q[1];
rz(-1.2756313) q[1];
sx q[1];
rz(0.035621312) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1195064) q[0];
sx q[0];
rz(-1.453244) q[0];
sx q[0];
rz(0.022679816) q[0];
x q[1];
rz(0.18892498) q[2];
sx q[2];
rz(-1.9686724) q[2];
sx q[2];
rz(2.973345) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5069711) q[1];
sx q[1];
rz(-1.325532) q[1];
sx q[1];
rz(0.62758201) q[1];
rz(-pi) q[2];
rz(-1.4485967) q[3];
sx q[3];
rz(-2.3482493) q[3];
sx q[3];
rz(2.6650708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.15141618) q[2];
sx q[2];
rz(-1.4750865) q[2];
sx q[2];
rz(1.0756294) q[2];
rz(0.52078024) q[3];
sx q[3];
rz(-2.5213089) q[3];
sx q[3];
rz(-1.552593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.90466475) q[0];
sx q[0];
rz(-1.0228782) q[0];
sx q[0];
rz(-0.9032332) q[0];
rz(-1.5029933) q[1];
sx q[1];
rz(-1.6939949) q[1];
sx q[1];
rz(1.2618056) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0340817) q[0];
sx q[0];
rz(-0.26824646) q[0];
sx q[0];
rz(-2.111042) q[0];
rz(-pi) q[1];
rz(-1.1177344) q[2];
sx q[2];
rz(-1.0894962) q[2];
sx q[2];
rz(2.520176) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8324229) q[1];
sx q[1];
rz(-1.6214402) q[1];
sx q[1];
rz(-3.0502351) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0699959) q[3];
sx q[3];
rz(-2.1163995) q[3];
sx q[3];
rz(0.16014447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5624076) q[2];
sx q[2];
rz(-0.6520485) q[2];
sx q[2];
rz(-3.0273666) q[2];
rz(-1.7646029) q[3];
sx q[3];
rz(-2.0012794) q[3];
sx q[3];
rz(2.637114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0196411) q[0];
sx q[0];
rz(-2.1421102) q[0];
sx q[0];
rz(-1.1720538) q[0];
rz(-1.3837586) q[1];
sx q[1];
rz(-1.9969321) q[1];
sx q[1];
rz(0.11464548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68814502) q[0];
sx q[0];
rz(-0.78399658) q[0];
sx q[0];
rz(-2.408124) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4769606) q[2];
sx q[2];
rz(-2.0619446) q[2];
sx q[2];
rz(2.0744086) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9108565) q[1];
sx q[1];
rz(-1.641526) q[1];
sx q[1];
rz(1.3225698) q[1];
rz(-pi) q[2];
rz(-0.10520868) q[3];
sx q[3];
rz(-2.3423839) q[3];
sx q[3];
rz(-1.5463643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5961479) q[2];
sx q[2];
rz(-0.20716509) q[2];
sx q[2];
rz(-2.2607415) q[2];
rz(-2.7018069) q[3];
sx q[3];
rz(-1.1652911) q[3];
sx q[3];
rz(2.0161276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3038444) q[0];
sx q[0];
rz(-1.4864018) q[0];
sx q[0];
rz(-3.0434171) q[0];
rz(0.57442609) q[1];
sx q[1];
rz(-1.4593294) q[1];
sx q[1];
rz(0.59304768) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0837729) q[0];
sx q[0];
rz(-1.860394) q[0];
sx q[0];
rz(-2.3168091) q[0];
x q[1];
rz(1.0306274) q[2];
sx q[2];
rz(-1.2048651) q[2];
sx q[2];
rz(-0.97971399) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8349466) q[1];
sx q[1];
rz(-3.1241841) q[1];
sx q[1];
rz(-2.2724292) q[1];
rz(2.5092068) q[3];
sx q[3];
rz(-1.2999855) q[3];
sx q[3];
rz(-0.89030823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9888931) q[2];
sx q[2];
rz(-1.1836735) q[2];
sx q[2];
rz(-0.93842554) q[2];
rz(-1.7484131) q[3];
sx q[3];
rz(-1.3104855) q[3];
sx q[3];
rz(2.9032629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68168454) q[0];
sx q[0];
rz(-2.5354698) q[0];
sx q[0];
rz(1.348362) q[0];
rz(2.0943191) q[1];
sx q[1];
rz(-2.2127163) q[1];
sx q[1];
rz(-0.97102078) q[1];
rz(2.1793096) q[2];
sx q[2];
rz(-0.078799993) q[2];
sx q[2];
rz(-2.7344631) q[2];
rz(-0.88225928) q[3];
sx q[3];
rz(-2.8591446) q[3];
sx q[3];
rz(2.1459116) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

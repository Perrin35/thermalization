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
rz(1.4242564) q[0];
sx q[0];
rz(-1.2196701) q[0];
sx q[0];
rz(-0.57504672) q[0];
rz(-2.9991034) q[1];
sx q[1];
rz(-1.2187076) q[1];
sx q[1];
rz(-2.277318) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9271547) q[0];
sx q[0];
rz(-1.3482643) q[0];
sx q[0];
rz(-0.19630918) q[0];
x q[1];
rz(1.7792542) q[2];
sx q[2];
rz(-2.499806) q[2];
sx q[2];
rz(-1.2520977) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4886538) q[1];
sx q[1];
rz(-1.0407742) q[1];
sx q[1];
rz(1.5684134) q[1];
rz(2.4349917) q[3];
sx q[3];
rz(-0.69701435) q[3];
sx q[3];
rz(1.5586588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.076866604) q[2];
sx q[2];
rz(-1.6419502) q[2];
sx q[2];
rz(2.1323252) q[2];
rz(0.81389728) q[3];
sx q[3];
rz(-0.32865694) q[3];
sx q[3];
rz(-2.8024659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0547884) q[0];
sx q[0];
rz(-1.4084933) q[0];
sx q[0];
rz(-2.8993697) q[0];
rz(-1.9107266) q[1];
sx q[1];
rz(-0.63175285) q[1];
sx q[1];
rz(0.79808527) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.674192) q[0];
sx q[0];
rz(-2.5134006) q[0];
sx q[0];
rz(2.9487361) q[0];
rz(-pi) q[1];
rz(1.913979) q[2];
sx q[2];
rz(-1.7084885) q[2];
sx q[2];
rz(-0.80114366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4646513) q[1];
sx q[1];
rz(-1.92917) q[1];
sx q[1];
rz(-2.6531069) q[1];
x q[2];
rz(-1.4105878) q[3];
sx q[3];
rz(-1.9247562) q[3];
sx q[3];
rz(0.33479553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4855087) q[2];
sx q[2];
rz(-1.069243) q[2];
sx q[2];
rz(2.3089224) q[2];
rz(-0.63140702) q[3];
sx q[3];
rz(-2.4000945) q[3];
sx q[3];
rz(-3.1413063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0176508) q[0];
sx q[0];
rz(-2.232382) q[0];
sx q[0];
rz(-0.18774524) q[0];
rz(-0.84367696) q[1];
sx q[1];
rz(-1.2742821) q[1];
sx q[1];
rz(-0.63293308) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7596066) q[0];
sx q[0];
rz(-1.6416024) q[0];
sx q[0];
rz(-1.9340865) q[0];
rz(2.723188) q[2];
sx q[2];
rz(-2.4103571) q[2];
sx q[2];
rz(-0.75218539) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.17281118) q[1];
sx q[1];
rz(-2.079614) q[1];
sx q[1];
rz(-1.225586) q[1];
rz(-pi) q[2];
rz(3.0121481) q[3];
sx q[3];
rz(-1.2264381) q[3];
sx q[3];
rz(-2.5710921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.26744276) q[2];
sx q[2];
rz(-2.1953857) q[2];
sx q[2];
rz(-1.8702742) q[2];
rz(1.4960131) q[3];
sx q[3];
rz(-1.6451903) q[3];
sx q[3];
rz(2.1700844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8895759) q[0];
sx q[0];
rz(-2.3887964) q[0];
sx q[0];
rz(0.65943199) q[0];
rz(1.8239498) q[1];
sx q[1];
rz(-1.8023856) q[1];
sx q[1];
rz(0.79016322) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20445828) q[0];
sx q[0];
rz(-2.7828628) q[0];
sx q[0];
rz(-2.782194) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1217693) q[2];
sx q[2];
rz(-1.6035282) q[2];
sx q[2];
rz(-2.6898877) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15406628) q[1];
sx q[1];
rz(-0.48168698) q[1];
sx q[1];
rz(-0.14194004) q[1];
rz(-pi) q[2];
x q[2];
rz(1.843256) q[3];
sx q[3];
rz(-2.108223) q[3];
sx q[3];
rz(-0.3934653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.059375199) q[2];
sx q[2];
rz(-1.6148115) q[2];
sx q[2];
rz(0.29602948) q[2];
rz(2.5791903) q[3];
sx q[3];
rz(-0.86807576) q[3];
sx q[3];
rz(0.11399046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.247308) q[0];
sx q[0];
rz(-1.9534651) q[0];
sx q[0];
rz(0.2071912) q[0];
rz(1.0317135) q[1];
sx q[1];
rz(-1.1528015) q[1];
sx q[1];
rz(1.656104) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9821859) q[0];
sx q[0];
rz(-0.80538087) q[0];
sx q[0];
rz(-2.5872487) q[0];
rz(-pi) q[1];
rz(1.1274952) q[2];
sx q[2];
rz(-2.5336899) q[2];
sx q[2];
rz(3.0435002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.42669) q[1];
sx q[1];
rz(-1.443092) q[1];
sx q[1];
rz(-0.35121484) q[1];
rz(3.033803) q[3];
sx q[3];
rz(-2.2985085) q[3];
sx q[3];
rz(-2.5586896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1032054) q[2];
sx q[2];
rz(-0.81192553) q[2];
sx q[2];
rz(1.1078328) q[2];
rz(2.2191018) q[3];
sx q[3];
rz(-1.4100217) q[3];
sx q[3];
rz(2.8973575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6237727) q[0];
sx q[0];
rz(-0.51384846) q[0];
sx q[0];
rz(2.9129831) q[0];
rz(0.27433431) q[1];
sx q[1];
rz(-2.3286596) q[1];
sx q[1];
rz(-2.6913604) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0749482) q[0];
sx q[0];
rz(-2.3111176) q[0];
sx q[0];
rz(-1.226783) q[0];
x q[1];
rz(-3.1338056) q[2];
sx q[2];
rz(-1.8000364) q[2];
sx q[2];
rz(2.0321329) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6321559) q[1];
sx q[1];
rz(-2.5773199) q[1];
sx q[1];
rz(0.48308259) q[1];
rz(-pi) q[2];
rz(-2.130533) q[3];
sx q[3];
rz(-2.6626867) q[3];
sx q[3];
rz(2.4541209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36394128) q[2];
sx q[2];
rz(-0.55299091) q[2];
sx q[2];
rz(-0.027776329) q[2];
rz(-0.76644301) q[3];
sx q[3];
rz(-1.1602297) q[3];
sx q[3];
rz(0.69507039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22208333) q[0];
sx q[0];
rz(-1.5064025) q[0];
sx q[0];
rz(-2.5191504) q[0];
rz(0.70603236) q[1];
sx q[1];
rz(-0.89203867) q[1];
sx q[1];
rz(0.58696729) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76510274) q[0];
sx q[0];
rz(-2.4060632) q[0];
sx q[0];
rz(2.2549948) q[0];
rz(2.6640011) q[2];
sx q[2];
rz(-0.53672635) q[2];
sx q[2];
rz(-2.0411185) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.5712292) q[1];
sx q[1];
rz(-2.0763358) q[1];
sx q[1];
rz(1.4897219) q[1];
rz(-pi) q[2];
x q[2];
rz(0.010384287) q[3];
sx q[3];
rz(-2.2828013) q[3];
sx q[3];
rz(0.66085789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.075835) q[2];
sx q[2];
rz(-1.7087874) q[2];
sx q[2];
rz(0.59448057) q[2];
rz(-2.8822656) q[3];
sx q[3];
rz(-1.2649053) q[3];
sx q[3];
rz(-1.9243141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038430564) q[0];
sx q[0];
rz(-0.90534627) q[0];
sx q[0];
rz(2.5627947) q[0];
rz(1.831306) q[1];
sx q[1];
rz(-1.879004) q[1];
sx q[1];
rz(-2.328918) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90234251) q[0];
sx q[0];
rz(-0.50241155) q[0];
sx q[0];
rz(2.8397296) q[0];
rz(-1.9695639) q[2];
sx q[2];
rz(-2.5426455) q[2];
sx q[2];
rz(-0.87397611) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.26749565) q[1];
sx q[1];
rz(-0.90522841) q[1];
sx q[1];
rz(-1.1049306) q[1];
x q[2];
rz(0.50847362) q[3];
sx q[3];
rz(-2.0805854) q[3];
sx q[3];
rz(-0.9289066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7909214) q[2];
sx q[2];
rz(-1.2988043) q[2];
sx q[2];
rz(-0.92362967) q[2];
rz(-2.1212497) q[3];
sx q[3];
rz(-2.2444921) q[3];
sx q[3];
rz(-0.11788192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0608805) q[0];
sx q[0];
rz(-1.7055644) q[0];
sx q[0];
rz(-1.4870148) q[0];
rz(0.05038536) q[1];
sx q[1];
rz(-0.81704187) q[1];
sx q[1];
rz(0.64687669) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4934702) q[0];
sx q[0];
rz(-2.1645344) q[0];
sx q[0];
rz(-2.6208682) q[0];
rz(-pi) q[1];
rz(-2.0503309) q[2];
sx q[2];
rz(-0.96194907) q[2];
sx q[2];
rz(0.11907585) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.17740956) q[1];
sx q[1];
rz(-2.6798267) q[1];
sx q[1];
rz(0.48514556) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7842125) q[3];
sx q[3];
rz(-1.4239256) q[3];
sx q[3];
rz(1.980912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6369624) q[2];
sx q[2];
rz(-1.5550104) q[2];
sx q[2];
rz(2.385425) q[2];
rz(1.470083) q[3];
sx q[3];
rz(-2.3556637) q[3];
sx q[3];
rz(-0.80529958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1038372) q[0];
sx q[0];
rz(-1.8917731) q[0];
sx q[0];
rz(-2.0334429) q[0];
rz(2.8773819) q[1];
sx q[1];
rz(-1.4690396) q[1];
sx q[1];
rz(2.5392551) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36985227) q[0];
sx q[0];
rz(-0.30769545) q[0];
sx q[0];
rz(2.0459396) q[0];
rz(-pi) q[1];
rz(0.17669038) q[2];
sx q[2];
rz(-1.7222705) q[2];
sx q[2];
rz(3.1164411) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.332676) q[1];
sx q[1];
rz(-0.68854587) q[1];
sx q[1];
rz(1.6469514) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4906795) q[3];
sx q[3];
rz(-1.7554566) q[3];
sx q[3];
rz(-2.1391275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5084874) q[2];
sx q[2];
rz(-0.54163951) q[2];
sx q[2];
rz(-0.87149054) q[2];
rz(-1.9320711) q[3];
sx q[3];
rz(-2.8612374) q[3];
sx q[3];
rz(2.5476294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7021983) q[0];
sx q[0];
rz(-1.3855423) q[0];
sx q[0];
rz(2.6227797) q[0];
rz(0.58204542) q[1];
sx q[1];
rz(-2.0379635) q[1];
sx q[1];
rz(1.4720974) q[1];
rz(0.40975801) q[2];
sx q[2];
rz(-0.77959218) q[2];
sx q[2];
rz(-2.6361781) q[2];
rz(-0.68861674) q[3];
sx q[3];
rz(-1.7027161) q[3];
sx q[3];
rz(-0.30218132) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

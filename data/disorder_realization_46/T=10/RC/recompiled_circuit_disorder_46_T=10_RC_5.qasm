OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5912136) q[0];
sx q[0];
rz(-0.0033291078) q[0];
sx q[0];
rz(-0.34531265) q[0];
rz(1.8058864) q[1];
sx q[1];
rz(3.4808573) q[1];
sx q[1];
rz(9.1453287) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0246668) q[0];
sx q[0];
rz(-1.3660087) q[0];
sx q[0];
rz(-1.8860399) q[0];
rz(-pi) q[1];
rz(-1.9136393) q[2];
sx q[2];
rz(-2.3228085) q[2];
sx q[2];
rz(-1.1358827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2742259) q[1];
sx q[1];
rz(-2.9510731) q[1];
sx q[1];
rz(1.5003367) q[1];
x q[2];
rz(-3.0885987) q[3];
sx q[3];
rz(-2.4989359) q[3];
sx q[3];
rz(1.7794123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.44148579) q[2];
sx q[2];
rz(-1.6690648) q[2];
sx q[2];
rz(0.95735615) q[2];
rz(0.29933128) q[3];
sx q[3];
rz(-0.39756164) q[3];
sx q[3];
rz(-2.7295952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9807724) q[0];
sx q[0];
rz(-2.1919724) q[0];
sx q[0];
rz(0.41369307) q[0];
rz(1.7970239) q[1];
sx q[1];
rz(-0.78318703) q[1];
sx q[1];
rz(-2.5059674) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7660852) q[0];
sx q[0];
rz(-1.5264395) q[0];
sx q[0];
rz(-1.4384598) q[0];
rz(-pi) q[1];
rz(-1.7027249) q[2];
sx q[2];
rz(-1.6172098) q[2];
sx q[2];
rz(-0.012243587) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2555274) q[1];
sx q[1];
rz(-0.35554245) q[1];
sx q[1];
rz(-2.38106) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.440968) q[3];
sx q[3];
rz(-1.0430338) q[3];
sx q[3];
rz(-0.81567314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3893434) q[2];
sx q[2];
rz(-1.8698591) q[2];
sx q[2];
rz(0.52865571) q[2];
rz(-1.901249) q[3];
sx q[3];
rz(-0.35959187) q[3];
sx q[3];
rz(2.8958029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5813331) q[0];
sx q[0];
rz(-1.9868877) q[0];
sx q[0];
rz(-2.8833959) q[0];
rz(1.6351581) q[1];
sx q[1];
rz(-2.5841027) q[1];
sx q[1];
rz(-2.7071276) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6601352) q[0];
sx q[0];
rz(-2.0148206) q[0];
sx q[0];
rz(-0.40159479) q[0];
rz(-pi) q[1];
rz(1.4591818) q[2];
sx q[2];
rz(-0.53225213) q[2];
sx q[2];
rz(-0.75512952) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1669958) q[1];
sx q[1];
rz(-1.3612862) q[1];
sx q[1];
rz(2.6796883) q[1];
x q[2];
rz(2.295899) q[3];
sx q[3];
rz(-2.0864668) q[3];
sx q[3];
rz(0.090279467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7210641) q[2];
sx q[2];
rz(-1.3238182) q[2];
sx q[2];
rz(-0.055796441) q[2];
rz(0.93786401) q[3];
sx q[3];
rz(-2.8580229) q[3];
sx q[3];
rz(-2.8985033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0980804) q[0];
sx q[0];
rz(-0.94399095) q[0];
sx q[0];
rz(0.14053024) q[0];
rz(-0.17164104) q[1];
sx q[1];
rz(-1.834603) q[1];
sx q[1];
rz(-2.8780639) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9182501) q[0];
sx q[0];
rz(-1.468987) q[0];
sx q[0];
rz(0.56901594) q[0];
x q[1];
rz(-1.0657004) q[2];
sx q[2];
rz(-2.1720338) q[2];
sx q[2];
rz(-2.2631753) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0861974) q[1];
sx q[1];
rz(-2.2436884) q[1];
sx q[1];
rz(-2.0170246) q[1];
rz(-pi) q[2];
rz(-1.2372381) q[3];
sx q[3];
rz(-1.3106489) q[3];
sx q[3];
rz(3.1317657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0984829) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(-1.8919224) q[2];
rz(-1.194687) q[3];
sx q[3];
rz(-1.0281111) q[3];
sx q[3];
rz(-0.67888129) q[3];
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
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2018305) q[0];
sx q[0];
rz(-1.285346) q[0];
sx q[0];
rz(3.1125267) q[0];
rz(-1.4020231) q[1];
sx q[1];
rz(-2.7840835) q[1];
sx q[1];
rz(-2.8129541) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6304566) q[0];
sx q[0];
rz(-1.6726613) q[0];
sx q[0];
rz(2.0589552) q[0];
x q[1];
rz(0.26942307) q[2];
sx q[2];
rz(-1.4002443) q[2];
sx q[2];
rz(1.558687) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.29341104) q[1];
sx q[1];
rz(-2.003258) q[1];
sx q[1];
rz(1.5007988) q[1];
rz(0.5248431) q[3];
sx q[3];
rz(-0.81916029) q[3];
sx q[3];
rz(-0.51845779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.11792004) q[2];
sx q[2];
rz(-2.6866044) q[2];
sx q[2];
rz(2.9809791) q[2];
rz(1.9953856) q[3];
sx q[3];
rz(-1.8278154) q[3];
sx q[3];
rz(0.62079352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49204957) q[0];
sx q[0];
rz(-2.2821125) q[0];
sx q[0];
rz(-0.78654003) q[0];
rz(-2.7644764) q[1];
sx q[1];
rz(-0.84723324) q[1];
sx q[1];
rz(-0.4424817) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839358) q[0];
sx q[0];
rz(-1.2120976) q[0];
sx q[0];
rz(-1.962612) q[0];
rz(1.2636678) q[2];
sx q[2];
rz(-1.0182292) q[2];
sx q[2];
rz(-2.8657258) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9554555) q[1];
sx q[1];
rz(-1.8147787) q[1];
sx q[1];
rz(1.9013491) q[1];
x q[2];
rz(0.093591452) q[3];
sx q[3];
rz(-1.2564661) q[3];
sx q[3];
rz(2.4623354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4601712) q[2];
sx q[2];
rz(-2.5632016) q[2];
sx q[2];
rz(-0.9712514) q[2];
rz(-0.69532895) q[3];
sx q[3];
rz(-2.9054502) q[3];
sx q[3];
rz(0.42738459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5044395) q[0];
sx q[0];
rz(-1.210286) q[0];
sx q[0];
rz(-2.1321645) q[0];
rz(2.5908453) q[1];
sx q[1];
rz(-1.8661205) q[1];
sx q[1];
rz(1.1632464) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79272072) q[0];
sx q[0];
rz(-2.7808245) q[0];
sx q[0];
rz(1.1611206) q[0];
rz(0.86261729) q[2];
sx q[2];
rz(-1.7628981) q[2];
sx q[2];
rz(-0.62193279) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0102331) q[1];
sx q[1];
rz(-1.7903923) q[1];
sx q[1];
rz(-1.4036199) q[1];
rz(-pi) q[2];
rz(2.8323152) q[3];
sx q[3];
rz(-2.2429372) q[3];
sx q[3];
rz(-2.6291763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78272468) q[2];
sx q[2];
rz(-1.9903368) q[2];
sx q[2];
rz(2.3434095) q[2];
rz(-2.362137) q[3];
sx q[3];
rz(-0.53648406) q[3];
sx q[3];
rz(1.9688169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1786132) q[0];
sx q[0];
rz(-2.6586752) q[0];
sx q[0];
rz(0.12284199) q[0];
rz(3.0154862) q[1];
sx q[1];
rz(-1.6364731) q[1];
sx q[1];
rz(1.925148) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1232237) q[0];
sx q[0];
rz(-2.1835678) q[0];
sx q[0];
rz(-1.0606674) q[0];
rz(-pi) q[1];
rz(-1.4037651) q[2];
sx q[2];
rz(-2.0207496) q[2];
sx q[2];
rz(3.0917167) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6718037) q[1];
sx q[1];
rz(-2.7046013) q[1];
sx q[1];
rz(2.6595518) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7580732) q[3];
sx q[3];
rz(-2.1547199) q[3];
sx q[3];
rz(-0.54959471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.2609666) q[2];
sx q[2];
rz(-2.5239021) q[2];
sx q[2];
rz(3.0984666) q[2];
rz(-0.17523266) q[3];
sx q[3];
rz(-2.2641116) q[3];
sx q[3];
rz(-1.5568679) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50424987) q[0];
sx q[0];
rz(-0.054333996) q[0];
sx q[0];
rz(-2.9328226) q[0];
rz(-1.6437221) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(-1.1245022) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5303406) q[0];
sx q[0];
rz(-1.4350622) q[0];
sx q[0];
rz(1.7087357) q[0];
x q[1];
rz(2.4268609) q[2];
sx q[2];
rz(-1.3681612) q[2];
sx q[2];
rz(1.5243901) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.81793663) q[1];
sx q[1];
rz(-1.4374104) q[1];
sx q[1];
rz(-0.0011841983) q[1];
rz(-pi) q[2];
rz(-0.48891588) q[3];
sx q[3];
rz(-0.93547869) q[3];
sx q[3];
rz(-2.2152701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0008529) q[2];
sx q[2];
rz(-0.72379392) q[2];
sx q[2];
rz(-2.9635584) q[2];
rz(-1.595165) q[3];
sx q[3];
rz(-1.833257) q[3];
sx q[3];
rz(2.6509638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35266018) q[0];
sx q[0];
rz(-1.1505609) q[0];
sx q[0];
rz(2.1066522) q[0];
rz(2.3433698) q[1];
sx q[1];
rz(-0.98033506) q[1];
sx q[1];
rz(-0.14990526) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7222544) q[0];
sx q[0];
rz(-0.94473487) q[0];
sx q[0];
rz(-1.1115848) q[0];
rz(-pi) q[1];
rz(-0.055585102) q[2];
sx q[2];
rz(-0.28700799) q[2];
sx q[2];
rz(0.74319786) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18969892) q[1];
sx q[1];
rz(-2.4687597) q[1];
sx q[1];
rz(0.876902) q[1];
rz(-pi) q[2];
rz(-0.67930119) q[3];
sx q[3];
rz(-1.7656109) q[3];
sx q[3];
rz(2.4904345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7297111) q[2];
sx q[2];
rz(-1.9103266) q[2];
sx q[2];
rz(-2.7330772) q[2];
rz(0.92489964) q[3];
sx q[3];
rz(-2.3291406) q[3];
sx q[3];
rz(-2.6598721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.9020486) q[0];
sx q[0];
rz(-1.5477381) q[0];
sx q[0];
rz(1.5292194) q[0];
rz(1.7655903) q[1];
sx q[1];
rz(-1.1767495) q[1];
sx q[1];
rz(-1.8935988) q[1];
rz(0.79509905) q[2];
sx q[2];
rz(-1.3394525) q[2];
sx q[2];
rz(0.48223334) q[2];
rz(-2.5216907) q[3];
sx q[3];
rz(-1.7696268) q[3];
sx q[3];
rz(2.2110248) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

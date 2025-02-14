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
rz(-0.37201878) q[0];
sx q[0];
rz(-2.7899185) q[0];
sx q[0];
rz(-0.055211842) q[0];
rz(-1.6959603) q[1];
sx q[1];
rz(-2.2382325) q[1];
sx q[1];
rz(0.13394314) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75823821) q[0];
sx q[0];
rz(-1.3516434) q[0];
sx q[0];
rz(1.6732277) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7705461) q[2];
sx q[2];
rz(-0.99989519) q[2];
sx q[2];
rz(-0.3112682) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6559927) q[1];
sx q[1];
rz(-0.47708407) q[1];
sx q[1];
rz(0.78039767) q[1];
rz(-1.6381959) q[3];
sx q[3];
rz(-2.615228) q[3];
sx q[3];
rz(1.1200861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.17406164) q[2];
sx q[2];
rz(-1.920819) q[2];
sx q[2];
rz(-2.3270712) q[2];
rz(-2.9461765) q[3];
sx q[3];
rz(-2.312909) q[3];
sx q[3];
rz(-2.101208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8622417) q[0];
sx q[0];
rz(-2.8793654) q[0];
sx q[0];
rz(-2.1694515) q[0];
rz(-2.5402918) q[1];
sx q[1];
rz(-2.1763132) q[1];
sx q[1];
rz(-1.345529) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20939553) q[0];
sx q[0];
rz(-0.69898134) q[0];
sx q[0];
rz(0.90468927) q[0];
x q[1];
rz(0.43429476) q[2];
sx q[2];
rz(-2.5611097) q[2];
sx q[2];
rz(-2.8934997) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7899155) q[1];
sx q[1];
rz(-0.43444217) q[1];
sx q[1];
rz(1.9138359) q[1];
rz(0.18109326) q[3];
sx q[3];
rz(-1.5646184) q[3];
sx q[3];
rz(-2.1712077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1874275) q[2];
sx q[2];
rz(-1.2792055) q[2];
sx q[2];
rz(0.97174755) q[2];
rz(2.1238972) q[3];
sx q[3];
rz(-1.965799) q[3];
sx q[3];
rz(0.051232256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.510842) q[0];
sx q[0];
rz(-2.172281) q[0];
sx q[0];
rz(-2.4208659) q[0];
rz(-0.12598704) q[1];
sx q[1];
rz(-2.4469913) q[1];
sx q[1];
rz(2.7350977) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0796156) q[0];
sx q[0];
rz(-0.70545371) q[0];
sx q[0];
rz(2.0464226) q[0];
rz(-pi) q[1];
rz(-2.1265246) q[2];
sx q[2];
rz(-0.90733084) q[2];
sx q[2];
rz(0.86554722) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49660044) q[1];
sx q[1];
rz(-0.82372249) q[1];
sx q[1];
rz(2.3053667) q[1];
rz(-pi) q[2];
rz(-1.2062827) q[3];
sx q[3];
rz(-1.8407093) q[3];
sx q[3];
rz(-2.4382537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5410109) q[2];
sx q[2];
rz(-1.3449679) q[2];
sx q[2];
rz(0.37008944) q[2];
rz(-0.078941405) q[3];
sx q[3];
rz(-2.168096) q[3];
sx q[3];
rz(-0.48648849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9321891) q[0];
sx q[0];
rz(-1.7904733) q[0];
sx q[0];
rz(-0.47602794) q[0];
rz(2.027482) q[1];
sx q[1];
rz(-2.1962491) q[1];
sx q[1];
rz(1.6260737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6081656) q[0];
sx q[0];
rz(-1.2100056) q[0];
sx q[0];
rz(1.7955304) q[0];
x q[1];
rz(-2.5427688) q[2];
sx q[2];
rz(-0.90617563) q[2];
sx q[2];
rz(1.7395626) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5596931) q[1];
sx q[1];
rz(-1.4861614) q[1];
sx q[1];
rz(-0.35178784) q[1];
rz(-0.8533303) q[3];
sx q[3];
rz(-1.322041) q[3];
sx q[3];
rz(-1.0568108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6473306) q[2];
sx q[2];
rz(-1.5212955) q[2];
sx q[2];
rz(-2.1935513) q[2];
rz(-2.7028132) q[3];
sx q[3];
rz(-2.572757) q[3];
sx q[3];
rz(-0.89890629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5206443) q[0];
sx q[0];
rz(-2.181894) q[0];
sx q[0];
rz(0.4278675) q[0];
rz(-2.1828792) q[1];
sx q[1];
rz(-1.2370647) q[1];
sx q[1];
rz(-1.0520891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89969544) q[0];
sx q[0];
rz(-2.1546531) q[0];
sx q[0];
rz(0.59052278) q[0];
rz(1.0361996) q[2];
sx q[2];
rz(-1.9048759) q[2];
sx q[2];
rz(2.6697347) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10994831) q[1];
sx q[1];
rz(-2.5606025) q[1];
sx q[1];
rz(-1.7419001) q[1];
rz(-pi) q[2];
rz(0.56628312) q[3];
sx q[3];
rz(-1.1193491) q[3];
sx q[3];
rz(-2.9873893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9713356) q[2];
sx q[2];
rz(-1.8361788) q[2];
sx q[2];
rz(2.7830284) q[2];
rz(1.1809008) q[3];
sx q[3];
rz(-2.2660393) q[3];
sx q[3];
rz(-0.53926474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3233258) q[0];
sx q[0];
rz(-1.8805255) q[0];
sx q[0];
rz(-0.64875025) q[0];
rz(-2.5541041) q[1];
sx q[1];
rz(-1.3222539) q[1];
sx q[1];
rz(2.9041451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64024583) q[0];
sx q[0];
rz(-3.0398439) q[0];
sx q[0];
rz(-2.8215088) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35570972) q[2];
sx q[2];
rz(-0.28537073) q[2];
sx q[2];
rz(-1.732638) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4202423) q[1];
sx q[1];
rz(-1.5974177) q[1];
sx q[1];
rz(-0.74810352) q[1];
rz(1.009672) q[3];
sx q[3];
rz(-2.3497407) q[3];
sx q[3];
rz(-0.71969024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6610403) q[2];
sx q[2];
rz(-0.66142267) q[2];
sx q[2];
rz(2.1858369) q[2];
rz(-1.7620979) q[3];
sx q[3];
rz(-2.5199514) q[3];
sx q[3];
rz(0.1563589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5697923) q[0];
sx q[0];
rz(-3.1389696) q[0];
sx q[0];
rz(-3.0963335) q[0];
rz(0.7943925) q[1];
sx q[1];
rz(-2.2519799) q[1];
sx q[1];
rz(-2.9387567) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5595374) q[0];
sx q[0];
rz(-1.7235316) q[0];
sx q[0];
rz(2.5688897) q[0];
x q[1];
rz(-2.9309996) q[2];
sx q[2];
rz(-0.74264975) q[2];
sx q[2];
rz(2.3499678) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0675373) q[1];
sx q[1];
rz(-2.592784) q[1];
sx q[1];
rz(-2.9850053) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84879227) q[3];
sx q[3];
rz(-1.5429075) q[3];
sx q[3];
rz(-1.3031808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5016735) q[2];
sx q[2];
rz(-1.4036274) q[2];
sx q[2];
rz(0.61140927) q[2];
rz(1.1008788) q[3];
sx q[3];
rz(-0.63488638) q[3];
sx q[3];
rz(-0.27975217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6708577) q[0];
sx q[0];
rz(-1.9453456) q[0];
sx q[0];
rz(-2.0580976) q[0];
rz(-2.1776543) q[1];
sx q[1];
rz(-1.0716535) q[1];
sx q[1];
rz(2.6458157) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14228798) q[0];
sx q[0];
rz(-1.5919627) q[0];
sx q[0];
rz(-1.643996) q[0];
rz(1.8476358) q[2];
sx q[2];
rz(-1.3949035) q[2];
sx q[2];
rz(-1.7664282) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.878022) q[1];
sx q[1];
rz(-1.9847893) q[1];
sx q[1];
rz(-2.1946226) q[1];
x q[2];
rz(-1.5257902) q[3];
sx q[3];
rz(-0.58820217) q[3];
sx q[3];
rz(2.3308995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.701391) q[2];
sx q[2];
rz(-2.3728366) q[2];
sx q[2];
rz(-2.4526147) q[2];
rz(1.6280599) q[3];
sx q[3];
rz(-0.83008927) q[3];
sx q[3];
rz(-2.1400863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5643519) q[0];
sx q[0];
rz(-1.5654726) q[0];
sx q[0];
rz(-1.0871357) q[0];
rz(2.3785036) q[1];
sx q[1];
rz(-1.3975846) q[1];
sx q[1];
rz(0.22689247) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2076465) q[0];
sx q[0];
rz(-1.3469704) q[0];
sx q[0];
rz(2.2730519) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2490602) q[2];
sx q[2];
rz(-1.2870803) q[2];
sx q[2];
rz(-0.64280451) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5530738) q[1];
sx q[1];
rz(-1.5344193) q[1];
sx q[1];
rz(0.22471551) q[1];
rz(-pi) q[2];
rz(-2.1710728) q[3];
sx q[3];
rz(-1.1423938) q[3];
sx q[3];
rz(2.4267933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.521296) q[2];
sx q[2];
rz(-0.93775788) q[2];
sx q[2];
rz(2.1916981) q[2];
rz(2.1096443) q[3];
sx q[3];
rz(-0.87939206) q[3];
sx q[3];
rz(0.43752813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6501605) q[0];
sx q[0];
rz(-3.037945) q[0];
sx q[0];
rz(-1.1520804) q[0];
rz(-3.0488455) q[1];
sx q[1];
rz(-2.4536965) q[1];
sx q[1];
rz(0.50382096) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4277305) q[0];
sx q[0];
rz(-0.52885884) q[0];
sx q[0];
rz(1.4725757) q[0];
x q[1];
rz(-0.26911084) q[2];
sx q[2];
rz(-1.7935441) q[2];
sx q[2];
rz(-2.3151223) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8246333) q[1];
sx q[1];
rz(-1.0976296) q[1];
sx q[1];
rz(1.5231569) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88155611) q[3];
sx q[3];
rz(-2.6668352) q[3];
sx q[3];
rz(0.17533824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6759701) q[2];
sx q[2];
rz(-2.0268107) q[2];
sx q[2];
rz(-2.5176804) q[2];
rz(-0.55656773) q[3];
sx q[3];
rz(-0.68816319) q[3];
sx q[3];
rz(0.19779675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9482166) q[0];
sx q[0];
rz(-2.2103136) q[0];
sx q[0];
rz(0.92934004) q[0];
rz(0.14840645) q[1];
sx q[1];
rz(-1.2444617) q[1];
sx q[1];
rz(-0.1958227) q[1];
rz(2.6779867) q[2];
sx q[2];
rz(-1.4624034) q[2];
sx q[2];
rz(0.69093888) q[2];
rz(1.8529057) q[3];
sx q[3];
rz(-0.28280453) q[3];
sx q[3];
rz(-0.40550532) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

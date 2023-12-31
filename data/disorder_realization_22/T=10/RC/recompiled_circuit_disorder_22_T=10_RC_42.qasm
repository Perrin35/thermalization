OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4322296) q[0];
sx q[0];
rz(-0.95786434) q[0];
sx q[0];
rz(-2.9971478) q[0];
rz(-2.5748409) q[1];
sx q[1];
rz(-2.6161939) q[1];
sx q[1];
rz(2.1638343) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35740556) q[0];
sx q[0];
rz(-1.6561243) q[0];
sx q[0];
rz(2.8822495) q[0];
rz(-pi) q[1];
rz(2.0830886) q[2];
sx q[2];
rz(-0.57527486) q[2];
sx q[2];
rz(0.011205999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8618776) q[1];
sx q[1];
rz(-1.4289083) q[1];
sx q[1];
rz(-0.14076294) q[1];
rz(-pi) q[2];
rz(-2.3542777) q[3];
sx q[3];
rz(-1.4437321) q[3];
sx q[3];
rz(2.3576759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.67291659) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(0.93227512) q[2];
rz(2.9428234) q[3];
sx q[3];
rz(-1.1105744) q[3];
sx q[3];
rz(0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.7070049) q[0];
sx q[0];
rz(-0.90536896) q[0];
sx q[0];
rz(0.36112753) q[0];
rz(1.4350285) q[1];
sx q[1];
rz(-1.7838493) q[1];
sx q[1];
rz(0.8180058) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084413962) q[0];
sx q[0];
rz(-1.1857496) q[0];
sx q[0];
rz(-1.8684698) q[0];
x q[1];
rz(1.1992707) q[2];
sx q[2];
rz(-1.4563592) q[2];
sx q[2];
rz(-2.3105846) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2745143) q[1];
sx q[1];
rz(-0.78650219) q[1];
sx q[1];
rz(-0.92814501) q[1];
rz(-pi) q[2];
rz(3.0844968) q[3];
sx q[3];
rz(-1.9340056) q[3];
sx q[3];
rz(-2.925194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8890185) q[2];
sx q[2];
rz(-0.44712862) q[2];
sx q[2];
rz(-1.4206295) q[2];
rz(-1.3160926) q[3];
sx q[3];
rz(-0.75794739) q[3];
sx q[3];
rz(0.38823286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4017568) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(2.2706568) q[0];
rz(2.8254106) q[1];
sx q[1];
rz(-0.28156391) q[1];
sx q[1];
rz(2.8443764) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57097829) q[0];
sx q[0];
rz(-1.1665205) q[0];
sx q[0];
rz(1.6398318) q[0];
rz(-pi) q[1];
rz(-2.4386114) q[2];
sx q[2];
rz(-2.2361122) q[2];
sx q[2];
rz(-2.2192628) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7797459) q[1];
sx q[1];
rz(-1.6184813) q[1];
sx q[1];
rz(1.7576799) q[1];
rz(-pi) q[2];
rz(-2.475297) q[3];
sx q[3];
rz(-2.1994281) q[3];
sx q[3];
rz(2.6578238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3729942) q[2];
sx q[2];
rz(-0.82295376) q[2];
sx q[2];
rz(0.42759839) q[2];
rz(-1.9528495) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(2.6141613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7052085) q[0];
sx q[0];
rz(-1.5773062) q[0];
sx q[0];
rz(-2.3676681) q[0];
rz(2.4286843) q[1];
sx q[1];
rz(-1.0293101) q[1];
sx q[1];
rz(-0.68177044) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6503158) q[0];
sx q[0];
rz(-2.804545) q[0];
sx q[0];
rz(-2.2469673) q[0];
rz(-pi) q[1];
rz(-0.96659987) q[2];
sx q[2];
rz(-1.7503498) q[2];
sx q[2];
rz(-2.0008848) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.53182787) q[1];
sx q[1];
rz(-2.1064639) q[1];
sx q[1];
rz(0.51171724) q[1];
rz(-pi) q[2];
rz(-1.2069615) q[3];
sx q[3];
rz(-1.893265) q[3];
sx q[3];
rz(-0.22842562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.008808) q[2];
sx q[2];
rz(-0.71295732) q[2];
sx q[2];
rz(-2.183389) q[2];
rz(2.0751674) q[3];
sx q[3];
rz(-1.8582148) q[3];
sx q[3];
rz(0.3796033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5565857) q[0];
sx q[0];
rz(-1.7946694) q[0];
sx q[0];
rz(-2.5812896) q[0];
rz(2.141748) q[1];
sx q[1];
rz(-0.20345774) q[1];
sx q[1];
rz(-1.6220185) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27846562) q[0];
sx q[0];
rz(-1.0882049) q[0];
sx q[0];
rz(1.0653853) q[0];
rz(-pi) q[1];
rz(-2.5298169) q[2];
sx q[2];
rz(-2.0423186) q[2];
sx q[2];
rz(-2.7437291) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.49415627) q[1];
sx q[1];
rz(-1.372822) q[1];
sx q[1];
rz(-2.3784749) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12445478) q[3];
sx q[3];
rz(-2.180047) q[3];
sx q[3];
rz(2.5973158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4524298) q[2];
sx q[2];
rz(-2.5897557) q[2];
sx q[2];
rz(-0.2229283) q[2];
rz(3.1068504) q[3];
sx q[3];
rz(-1.7545173) q[3];
sx q[3];
rz(0.071578659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8054304) q[0];
sx q[0];
rz(-0.31328377) q[0];
sx q[0];
rz(2.0715332) q[0];
rz(-1.7806212) q[1];
sx q[1];
rz(-0.37934163) q[1];
sx q[1];
rz(1.8575352) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491935) q[0];
sx q[0];
rz(-1.0214897) q[0];
sx q[0];
rz(-0.71676371) q[0];
rz(-1.1805004) q[2];
sx q[2];
rz(-1.3109866) q[2];
sx q[2];
rz(0.63894546) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18404993) q[1];
sx q[1];
rz(-1.8584538) q[1];
sx q[1];
rz(0.48344739) q[1];
rz(1.0809903) q[3];
sx q[3];
rz(-2.5488857) q[3];
sx q[3];
rz(-2.5268775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6665035) q[2];
sx q[2];
rz(-1.44375) q[2];
sx q[2];
rz(-0.63759032) q[2];
rz(-0.83667886) q[3];
sx q[3];
rz(-1.1058608) q[3];
sx q[3];
rz(-1.3440514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5752983) q[0];
sx q[0];
rz(-2.7116382) q[0];
sx q[0];
rz(-2.5740525) q[0];
rz(2.7138846) q[1];
sx q[1];
rz(-1.6141012) q[1];
sx q[1];
rz(2.2033851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0465614) q[0];
sx q[0];
rz(-1.1633658) q[0];
sx q[0];
rz(-0.33336063) q[0];
rz(-2.4110255) q[2];
sx q[2];
rz(-1.1267203) q[2];
sx q[2];
rz(-2.4161352) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5349622) q[1];
sx q[1];
rz(-0.76430799) q[1];
sx q[1];
rz(-1.0431837) q[1];
rz(-0.35258099) q[3];
sx q[3];
rz(-2.4028824) q[3];
sx q[3];
rz(2.0032361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.87970916) q[2];
sx q[2];
rz(-1.9677013) q[2];
sx q[2];
rz(-1.3860469) q[2];
rz(1.8188247) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(3.1125606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26738527) q[0];
sx q[0];
rz(-0.33339849) q[0];
sx q[0];
rz(1.7077131) q[0];
rz(1.2738312) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(0.83126718) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3812105) q[0];
sx q[0];
rz(-1.7179278) q[0];
sx q[0];
rz(0.40572625) q[0];
x q[1];
rz(0.93584658) q[2];
sx q[2];
rz(-1.4668462) q[2];
sx q[2];
rz(-0.39633358) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.89325209) q[1];
sx q[1];
rz(-0.89683956) q[1];
sx q[1];
rz(2.2158951) q[1];
rz(-pi) q[2];
rz(-0.23037489) q[3];
sx q[3];
rz(-1.3532234) q[3];
sx q[3];
rz(2.2090467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5635809) q[2];
sx q[2];
rz(-1.7687904) q[2];
sx q[2];
rz(2.1772299) q[2];
rz(-1.1635121) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(-2.4826629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7678541) q[0];
sx q[0];
rz(-2.4066194) q[0];
sx q[0];
rz(-2.1642165) q[0];
rz(1.3865698) q[1];
sx q[1];
rz(-1.8354548) q[1];
sx q[1];
rz(-1.1057373) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10406636) q[0];
sx q[0];
rz(-2.6797047) q[0];
sx q[0];
rz(2.877263) q[0];
rz(-1.703891) q[2];
sx q[2];
rz(-0.67099748) q[2];
sx q[2];
rz(0.69588307) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8399664) q[1];
sx q[1];
rz(-1.7558388) q[1];
sx q[1];
rz(1.806083) q[1];
rz(-2.7902778) q[3];
sx q[3];
rz(-1.2993386) q[3];
sx q[3];
rz(0.16004496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.608312) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(-0.14979714) q[2];
rz(-1.3730565) q[3];
sx q[3];
rz(-1.405973) q[3];
sx q[3];
rz(-1.8201374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6479284) q[0];
sx q[0];
rz(-2.6239008) q[0];
sx q[0];
rz(-3.0143484) q[0];
rz(1.4808902) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(-0.1677992) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4807178) q[0];
sx q[0];
rz(-1.6010451) q[0];
sx q[0];
rz(1.589993) q[0];
rz(1.7255515) q[2];
sx q[2];
rz(-2.1733279) q[2];
sx q[2];
rz(3.0263911) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.64393109) q[1];
sx q[1];
rz(-1.5639468) q[1];
sx q[1];
rz(-0.58845206) q[1];
rz(-pi) q[2];
rz(3.0182748) q[3];
sx q[3];
rz(-1.4535558) q[3];
sx q[3];
rz(2.5532212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.026713513) q[2];
sx q[2];
rz(-2.202704) q[2];
sx q[2];
rz(-0.76114571) q[2];
rz(3.051565) q[3];
sx q[3];
rz(-2.138425) q[3];
sx q[3];
rz(-0.95054039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54820838) q[0];
sx q[0];
rz(-1.1593288) q[0];
sx q[0];
rz(2.819084) q[0];
rz(-2.7535915) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(-0.8969174) q[2];
sx q[2];
rz(-2.2581836) q[2];
sx q[2];
rz(2.7473292) q[2];
rz(-1.8518944) q[3];
sx q[3];
rz(-0.55681183) q[3];
sx q[3];
rz(-1.1435215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

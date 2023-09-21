OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(-0.12806211) q[0];
sx q[0];
rz(-2.3242216) q[0];
rz(-2.1583537) q[1];
sx q[1];
rz(-2.6020738) q[1];
sx q[1];
rz(-1.9411545) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1737328) q[0];
sx q[0];
rz(-1.1636486) q[0];
sx q[0];
rz(1.834154) q[0];
rz(-1.0317694) q[2];
sx q[2];
rz(-1.6845778) q[2];
sx q[2];
rz(-3.0381418) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5719205) q[1];
sx q[1];
rz(-1.8351646) q[1];
sx q[1];
rz(-0.7260679) q[1];
x q[2];
rz(-2.3463763) q[3];
sx q[3];
rz(-1.8024369) q[3];
sx q[3];
rz(2.8455646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0212705) q[2];
sx q[2];
rz(-2.5734084) q[2];
sx q[2];
rz(-1.583064) q[2];
rz(-0.99672404) q[3];
sx q[3];
rz(-2.6895027) q[3];
sx q[3];
rz(0.42580095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5581756) q[0];
sx q[0];
rz(-0.75220627) q[0];
sx q[0];
rz(-0.054071991) q[0];
rz(-1.9460829) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(2.6057459) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0756404) q[0];
sx q[0];
rz(-1.6951121) q[0];
sx q[0];
rz(-2.1527704) q[0];
rz(2.4321062) q[2];
sx q[2];
rz(-1.3777133) q[2];
sx q[2];
rz(2.513934) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.10837308) q[1];
sx q[1];
rz(-2.6544826) q[1];
sx q[1];
rz(2.5750722) q[1];
rz(1.7137394) q[3];
sx q[3];
rz(-1.6030451) q[3];
sx q[3];
rz(-1.3168207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1198931) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(-2.1832441) q[2];
rz(-0.066453233) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(-0.44979969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7217343) q[0];
sx q[0];
rz(-1.9323213) q[0];
sx q[0];
rz(0.15047519) q[0];
rz(-0.45723215) q[1];
sx q[1];
rz(-2.229264) q[1];
sx q[1];
rz(-3.1157852) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31608554) q[0];
sx q[0];
rz(-1.8206017) q[0];
sx q[0];
rz(3.120963) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.843156) q[2];
sx q[2];
rz(-0.32264999) q[2];
sx q[2];
rz(-1.1622365) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4986213) q[1];
sx q[1];
rz(-1.1585304) q[1];
sx q[1];
rz(0.69795124) q[1];
x q[2];
rz(-1.953655) q[3];
sx q[3];
rz(-1.448505) q[3];
sx q[3];
rz(1.0505291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1228483) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(-2.5562111) q[2];
rz(2.9600926) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(1.5649149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.240775) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(2.9649819) q[0];
rz(2.2606842) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(-0.53612971) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6484084) q[0];
sx q[0];
rz(-0.83183653) q[0];
sx q[0];
rz(1.0113082) q[0];
rz(-pi) q[1];
rz(-1.5136112) q[2];
sx q[2];
rz(-1.3344889) q[2];
sx q[2];
rz(1.1255217) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.22524658) q[1];
sx q[1];
rz(-0.9496453) q[1];
sx q[1];
rz(-2.0398554) q[1];
rz(2.4651277) q[3];
sx q[3];
rz(-1.8026661) q[3];
sx q[3];
rz(-1.7959309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46999103) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(-2.0969351) q[2];
rz(2.4345543) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(-2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2351284) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(1.0513603) q[0];
rz(-1.6479187) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(-0.043118127) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2152104) q[0];
sx q[0];
rz(-2.1676817) q[0];
sx q[0];
rz(-0.2290639) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0518603) q[2];
sx q[2];
rz(-1.8199925) q[2];
sx q[2];
rz(-0.11617004) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.67152126) q[1];
sx q[1];
rz(-0.38362353) q[1];
sx q[1];
rz(-0.77312153) q[1];
rz(-pi) q[2];
rz(2.99302) q[3];
sx q[3];
rz(-0.79509495) q[3];
sx q[3];
rz(0.038392301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.12895) q[2];
sx q[2];
rz(-0.49762112) q[2];
sx q[2];
rz(0.35219231) q[2];
rz(2.5514065) q[3];
sx q[3];
rz(-2.6679109) q[3];
sx q[3];
rz(-2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6234289) q[0];
sx q[0];
rz(-1.2151027) q[0];
sx q[0];
rz(-1.1556926) q[0];
rz(0.75025264) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(2.0828784) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4689894) q[0];
sx q[0];
rz(-2.6575343) q[0];
sx q[0];
rz(0.46033916) q[0];
x q[1];
rz(-2.6201453) q[2];
sx q[2];
rz(-1.8910742) q[2];
sx q[2];
rz(2.5495868) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.42029542) q[1];
sx q[1];
rz(-1.2696206) q[1];
sx q[1];
rz(0.23423127) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9162354) q[3];
sx q[3];
rz(-2.4212824) q[3];
sx q[3];
rz(2.5845137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6713312) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(-1.2188101) q[2];
rz(1.9865215) q[3];
sx q[3];
rz(-2.9030436) q[3];
sx q[3];
rz(-1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041615151) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(-1.51145) q[0];
rz(1.7639683) q[1];
sx q[1];
rz(-0.310985) q[1];
sx q[1];
rz(-0.84164936) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8182897) q[0];
sx q[0];
rz(-1.0536195) q[0];
sx q[0];
rz(1.3099758) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.878703) q[2];
sx q[2];
rz(-2.4979696) q[2];
sx q[2];
rz(0.4292683) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0591653) q[1];
sx q[1];
rz(-1.4945684) q[1];
sx q[1];
rz(0.04870292) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8089201) q[3];
sx q[3];
rz(-2.3861285) q[3];
sx q[3];
rz(-1.8734224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1402309) q[2];
sx q[2];
rz(-1.7732239) q[2];
sx q[2];
rz(0.049953071) q[2];
rz(-0.66155457) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426303) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(-1.7393973) q[0];
rz(0.095480355) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(0.41762525) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0324875) q[0];
sx q[0];
rz(-2.0628477) q[0];
sx q[0];
rz(2.0672654) q[0];
rz(-pi) q[1];
rz(-1.9365963) q[2];
sx q[2];
rz(-1.4636283) q[2];
sx q[2];
rz(-2.0049713) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8466134) q[1];
sx q[1];
rz(-1.1049005) q[1];
sx q[1];
rz(2.2094775) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29406677) q[3];
sx q[3];
rz(-1.9293647) q[3];
sx q[3];
rz(2.8420574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3461356) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(-1.0428838) q[2];
rz(-2.4677094) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(-0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(1.8326571) q[0];
sx q[0];
rz(-1.7333663) q[0];
sx q[0];
rz(0.62966627) q[0];
rz(0.57485238) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(2.1946857) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6553584) q[0];
sx q[0];
rz(-1.8706733) q[0];
sx q[0];
rz(2.8006058) q[0];
x q[1];
rz(-1.6106748) q[2];
sx q[2];
rz(-1.4738184) q[2];
sx q[2];
rz(3.1181042) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2448333) q[1];
sx q[1];
rz(-1.4866801) q[1];
sx q[1];
rz(0.33478488) q[1];
rz(2.9506748) q[3];
sx q[3];
rz(-1.1357765) q[3];
sx q[3];
rz(0.85489475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.56069121) q[2];
sx q[2];
rz(-2.0118482) q[2];
sx q[2];
rz(1.2488731) q[2];
rz(0.71436626) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0749851) q[0];
sx q[0];
rz(-0.62244901) q[0];
sx q[0];
rz(-0.65504909) q[0];
rz(-0.89637268) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(2.4972829) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8048332) q[0];
sx q[0];
rz(-1.6258437) q[0];
sx q[0];
rz(2.9677662) q[0];
x q[1];
rz(-2.9115453) q[2];
sx q[2];
rz(-0.86798475) q[2];
sx q[2];
rz(2.7237797) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.88313738) q[1];
sx q[1];
rz(-1.2346134) q[1];
sx q[1];
rz(-1.8717481) q[1];
rz(-pi) q[2];
rz(-1.0993016) q[3];
sx q[3];
rz(-1.1700556) q[3];
sx q[3];
rz(-0.61939643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7908988) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(-2.541686) q[2];
rz(0.89896262) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(-1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.8469289) q[0];
sx q[0];
rz(-1.3273205) q[0];
sx q[0];
rz(2.474665) q[0];
rz(2.9121493) q[1];
sx q[1];
rz(-2.2506917) q[1];
sx q[1];
rz(-3.0058203) q[1];
rz(2.1951998) q[2];
sx q[2];
rz(-1.6740435) q[2];
sx q[2];
rz(-1.2798535) q[2];
rz(2.7502144) q[3];
sx q[3];
rz(-2.2073675) q[3];
sx q[3];
rz(-1.7195306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

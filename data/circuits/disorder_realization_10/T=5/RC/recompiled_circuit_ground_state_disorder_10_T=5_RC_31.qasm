OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6634231) q[0];
sx q[0];
rz(-2.2278251) q[0];
sx q[0];
rz(-1.0650286) q[0];
rz(-1.5364667) q[1];
sx q[1];
rz(-2.0361587) q[1];
sx q[1];
rz(3.0450568) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8138344) q[0];
sx q[0];
rz(-1.7760491) q[0];
sx q[0];
rz(-1.9987518) q[0];
rz(-1.3209302) q[2];
sx q[2];
rz(-2.0786946) q[2];
sx q[2];
rz(2.8405416) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.3693856) q[1];
sx q[1];
rz(-2.188538) q[1];
sx q[1];
rz(-1.2031612) q[1];
x q[2];
rz(-1.7126585) q[3];
sx q[3];
rz(-0.58462287) q[3];
sx q[3];
rz(1.4277924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1221293) q[2];
sx q[2];
rz(-2.9563642) q[2];
sx q[2];
rz(1.2629925) q[2];
rz(-2.7428135) q[3];
sx q[3];
rz(-2.0006657) q[3];
sx q[3];
rz(-2.1122475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1102092) q[0];
sx q[0];
rz(-2.5094014) q[0];
sx q[0];
rz(-1.6803918) q[0];
rz(-0.070143135) q[1];
sx q[1];
rz(-1.5488011) q[1];
sx q[1];
rz(2.7102914) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9316783) q[0];
sx q[0];
rz(-3.0793858) q[0];
sx q[0];
rz(-1.4181523) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1761364) q[2];
sx q[2];
rz(-0.49412974) q[2];
sx q[2];
rz(-2.663341) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0924393) q[1];
sx q[1];
rz(-0.91040694) q[1];
sx q[1];
rz(2.4699442) q[1];
rz(2.7106695) q[3];
sx q[3];
rz(-0.049073372) q[3];
sx q[3];
rz(1.2653093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9554319) q[2];
sx q[2];
rz(-1.7308851) q[2];
sx q[2];
rz(-0.53039941) q[2];
rz(-1.0839328) q[3];
sx q[3];
rz(-2.0518186) q[3];
sx q[3];
rz(-1.8918022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0386061) q[0];
sx q[0];
rz(-0.57612053) q[0];
sx q[0];
rz(-1.9392133) q[0];
rz(2.1309958) q[1];
sx q[1];
rz(-1.3253515) q[1];
sx q[1];
rz(-0.030166322) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9011616) q[0];
sx q[0];
rz(-1.5053006) q[0];
sx q[0];
rz(-1.5237048) q[0];
rz(-pi) q[1];
rz(-1.6232477) q[2];
sx q[2];
rz(-1.5546335) q[2];
sx q[2];
rz(2.9132622) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67913342) q[1];
sx q[1];
rz(-2.7736385) q[1];
sx q[1];
rz(0.36104135) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10134931) q[3];
sx q[3];
rz(-2.1411485) q[3];
sx q[3];
rz(-1.0265019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5304337) q[2];
sx q[2];
rz(-0.23808372) q[2];
sx q[2];
rz(-2.7670624) q[2];
rz(-2.020906) q[3];
sx q[3];
rz(-1.6569542) q[3];
sx q[3];
rz(-2.4499031) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94367868) q[0];
sx q[0];
rz(-1.2926084) q[0];
sx q[0];
rz(1.9669272) q[0];
rz(-1.3150175) q[1];
sx q[1];
rz(-2.0349793) q[1];
sx q[1];
rz(-0.16474251) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68130601) q[0];
sx q[0];
rz(-0.78382712) q[0];
sx q[0];
rz(-0.71390217) q[0];
rz(-pi) q[1];
rz(1.4061178) q[2];
sx q[2];
rz(-1.0633755) q[2];
sx q[2];
rz(-2.8230482) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8604728) q[1];
sx q[1];
rz(-1.9083284) q[1];
sx q[1];
rz(-2.0915477) q[1];
rz(-pi) q[2];
rz(3.1026671) q[3];
sx q[3];
rz(-1.64216) q[3];
sx q[3];
rz(1.4369278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.69079196) q[2];
sx q[2];
rz(-1.6054634) q[2];
sx q[2];
rz(2.6353321) q[2];
rz(2.6311503) q[3];
sx q[3];
rz(-0.78845316) q[3];
sx q[3];
rz(-2.3300664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6798169) q[0];
sx q[0];
rz(-2.8051069) q[0];
sx q[0];
rz(-2.3315954) q[0];
rz(-3.0432155) q[1];
sx q[1];
rz(-0.38946238) q[1];
sx q[1];
rz(2.817165) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8036486) q[0];
sx q[0];
rz(-1.4483101) q[0];
sx q[0];
rz(-0.0019367812) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4348898) q[2];
sx q[2];
rz(-1.3506442) q[2];
sx q[2];
rz(-2.1957601) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4091332) q[1];
sx q[1];
rz(-0.9995102) q[1];
sx q[1];
rz(-1.7833227) q[1];
rz(-1.3611322) q[3];
sx q[3];
rz(-1.1034605) q[3];
sx q[3];
rz(-2.1356104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7847932) q[2];
sx q[2];
rz(-1.5609799) q[2];
sx q[2];
rz(-2.5830833) q[2];
rz(1.8465346) q[3];
sx q[3];
rz(-1.7241155) q[3];
sx q[3];
rz(0.3524802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2297939) q[0];
sx q[0];
rz(-2.7138382) q[0];
sx q[0];
rz(-1.7708923) q[0];
rz(-1.7478583) q[1];
sx q[1];
rz(-2.1262028) q[1];
sx q[1];
rz(0.14399354) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27731652) q[0];
sx q[0];
rz(-1.4668408) q[0];
sx q[0];
rz(0.070965537) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1156111) q[2];
sx q[2];
rz(-1.3521638) q[2];
sx q[2];
rz(1.8393379) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.6274144) q[1];
sx q[1];
rz(-2.5251221) q[1];
sx q[1];
rz(-1.3759088) q[1];
x q[2];
rz(-0.64231974) q[3];
sx q[3];
rz(-1.910247) q[3];
sx q[3];
rz(2.385115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7250942) q[2];
sx q[2];
rz(-2.7855253) q[2];
sx q[2];
rz(0.98102513) q[2];
rz(0.52870005) q[3];
sx q[3];
rz(-1.3542391) q[3];
sx q[3];
rz(-2.5813812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33787802) q[0];
sx q[0];
rz(-0.80334544) q[0];
sx q[0];
rz(1.2038318) q[0];
rz(-3.1375258) q[1];
sx q[1];
rz(-1.4264359) q[1];
sx q[1];
rz(2.8003661) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1004286) q[0];
sx q[0];
rz(-0.89027864) q[0];
sx q[0];
rz(-1.0170223) q[0];
rz(-pi) q[1];
rz(-2.262142) q[2];
sx q[2];
rz(-1.571889) q[2];
sx q[2];
rz(-3.080883) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3292405) q[1];
sx q[1];
rz(-1.4218378) q[1];
sx q[1];
rz(-1.5801058) q[1];
rz(-pi) q[2];
rz(0.39088313) q[3];
sx q[3];
rz(-0.36827403) q[3];
sx q[3];
rz(-1.0523018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61871201) q[2];
sx q[2];
rz(-1.8992385) q[2];
sx q[2];
rz(0.27113554) q[2];
rz(-1.6297657) q[3];
sx q[3];
rz(-2.1745493) q[3];
sx q[3];
rz(-0.40201521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4790344) q[0];
sx q[0];
rz(-2.5738578) q[0];
sx q[0];
rz(2.973279) q[0];
rz(2.1315101) q[1];
sx q[1];
rz(-1.1770266) q[1];
sx q[1];
rz(0.10028663) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94674695) q[0];
sx q[0];
rz(-1.8528209) q[0];
sx q[0];
rz(-0.32026024) q[0];
rz(0.19817721) q[2];
sx q[2];
rz(-0.74286425) q[2];
sx q[2];
rz(0.6720378) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5096693) q[1];
sx q[1];
rz(-2.9384216) q[1];
sx q[1];
rz(-3.0198899) q[1];
rz(-pi) q[2];
rz(-2.4504205) q[3];
sx q[3];
rz(-0.6414957) q[3];
sx q[3];
rz(-0.79090276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4526796) q[2];
sx q[2];
rz(-2.6426688) q[2];
sx q[2];
rz(-1.4208687) q[2];
rz(1.0639327) q[3];
sx q[3];
rz(-1.7199793) q[3];
sx q[3];
rz(-0.091778278) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7559779) q[0];
sx q[0];
rz(-1.6843963) q[0];
sx q[0];
rz(-3.1047367) q[0];
rz(-1.256975) q[1];
sx q[1];
rz(-1.6391862) q[1];
sx q[1];
rz(1.3704376) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94266701) q[0];
sx q[0];
rz(-1.725993) q[0];
sx q[0];
rz(-1.8465081) q[0];
rz(2.7671934) q[2];
sx q[2];
rz(-2.7485195) q[2];
sx q[2];
rz(-3.0469325) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6681666) q[1];
sx q[1];
rz(-1.8500195) q[1];
sx q[1];
rz(-0.88069005) q[1];
rz(-pi) q[2];
rz(1.1537136) q[3];
sx q[3];
rz(-1.5305102) q[3];
sx q[3];
rz(-1.861544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2799985) q[2];
sx q[2];
rz(-1.3759321) q[2];
sx q[2];
rz(2.9691248) q[2];
rz(0.55781588) q[3];
sx q[3];
rz(-0.54111257) q[3];
sx q[3];
rz(-1.7012168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7114792) q[0];
sx q[0];
rz(-0.63876307) q[0];
sx q[0];
rz(-0.3399671) q[0];
rz(-0.66554794) q[1];
sx q[1];
rz(-1.7219209) q[1];
sx q[1];
rz(1.8284214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68379867) q[0];
sx q[0];
rz(-2.9556985) q[0];
sx q[0];
rz(2.8838975) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4943141) q[2];
sx q[2];
rz(-1.8308365) q[2];
sx q[2];
rz(-1.2225012) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1024627) q[1];
sx q[1];
rz(-0.38929554) q[1];
sx q[1];
rz(-2.2401458) q[1];
rz(0.40295593) q[3];
sx q[3];
rz(-1.7149441) q[3];
sx q[3];
rz(0.47553167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7790935) q[2];
sx q[2];
rz(-0.43467793) q[2];
sx q[2];
rz(2.4707826) q[2];
rz(-0.32495156) q[3];
sx q[3];
rz(-0.81348014) q[3];
sx q[3];
rz(-2.1285098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.308607) q[0];
sx q[0];
rz(-1.2815463) q[0];
sx q[0];
rz(1.4776342) q[0];
rz(-1.0302522) q[1];
sx q[1];
rz(-1.4719084) q[1];
sx q[1];
rz(1.7869064) q[1];
rz(-1.0679097) q[2];
sx q[2];
rz(-1.1008762) q[2];
sx q[2];
rz(-1.8864143) q[2];
rz(0.53976157) q[3];
sx q[3];
rz(-2.4884347) q[3];
sx q[3];
rz(1.6275737) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(0.69005203) q[0];
sx q[0];
rz(8.5741841) q[0];
sx q[0];
rz(9.2246715) q[0];
rz(-0.19835681) q[1];
sx q[1];
rz(-0.52302066) q[1];
sx q[1];
rz(1.3318292) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057077335) q[0];
sx q[0];
rz(-2.6383218) q[0];
sx q[0];
rz(1.0751192) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5426504) q[2];
sx q[2];
rz(-1.9859196) q[2];
sx q[2];
rz(0.51905635) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.15439776) q[1];
sx q[1];
rz(-1.618865) q[1];
sx q[1];
rz(0.67832077) q[1];
rz(-pi) q[2];
rz(-1.7302527) q[3];
sx q[3];
rz(-2.660369) q[3];
sx q[3];
rz(0.10018292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.072210463) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(-1.8270095) q[2];
rz(-2.2383111) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(2.7956853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2717993) q[0];
sx q[0];
rz(-1.6558187) q[0];
sx q[0];
rz(-2.2858009) q[0];
rz(2.3675512) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(-2.3537297) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44875328) q[0];
sx q[0];
rz(-1.8215382) q[0];
sx q[0];
rz(3.1186799) q[0];
rz(-pi) q[1];
rz(-1.6131975) q[2];
sx q[2];
rz(-1.2268492) q[2];
sx q[2];
rz(2.6478772) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1545002) q[1];
sx q[1];
rz(-2.5804511) q[1];
sx q[1];
rz(-0.89548703) q[1];
rz(-pi) q[2];
rz(1.3608906) q[3];
sx q[3];
rz(-2.2271101) q[3];
sx q[3];
rz(-0.24649749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0756125) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(-0.39609972) q[2];
rz(1.5710477) q[3];
sx q[3];
rz(-1.6377662) q[3];
sx q[3];
rz(-1.3554696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1876672) q[0];
sx q[0];
rz(-1.6463771) q[0];
sx q[0];
rz(-3.0809825) q[0];
rz(-2.4568779) q[1];
sx q[1];
rz(-1.741332) q[1];
sx q[1];
rz(0.33904591) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44797541) q[0];
sx q[0];
rz(-2.1537499) q[0];
sx q[0];
rz(-1.5956123) q[0];
rz(-1.4475559) q[2];
sx q[2];
rz(-1.3831105) q[2];
sx q[2];
rz(-1.0836243) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3587957) q[1];
sx q[1];
rz(-0.41886815) q[1];
sx q[1];
rz(-2.2936506) q[1];
rz(0.77683461) q[3];
sx q[3];
rz(-1.562444) q[3];
sx q[3];
rz(1.2119105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.48245779) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(-2.6711312) q[2];
rz(0.86236924) q[3];
sx q[3];
rz(-0.44895288) q[3];
sx q[3];
rz(1.580015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793133) q[0];
sx q[0];
rz(-1.1184432) q[0];
sx q[0];
rz(-0.40801868) q[0];
rz(-1.7904003) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(-1.8203576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21463284) q[0];
sx q[0];
rz(-2.0713191) q[0];
sx q[0];
rz(-1.2248125) q[0];
x q[1];
rz(0.91561134) q[2];
sx q[2];
rz(-1.639099) q[2];
sx q[2];
rz(-3.0223522) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.8157573) q[1];
sx q[1];
rz(-1.6920631) q[1];
sx q[1];
rz(1.497529) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7144954) q[3];
sx q[3];
rz(-1.1833041) q[3];
sx q[3];
rz(0.96086335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2230175) q[2];
sx q[2];
rz(-1.431798) q[2];
sx q[2];
rz(2.4187386) q[2];
rz(0.84960788) q[3];
sx q[3];
rz(-1.1725715) q[3];
sx q[3];
rz(2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1866813) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(-2.5816259) q[0];
rz(-0.2050744) q[1];
sx q[1];
rz(-1.8849025) q[1];
sx q[1];
rz(-0.29108873) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7322892) q[0];
sx q[0];
rz(-1.8390053) q[0];
sx q[0];
rz(1.7464253) q[0];
rz(-pi) q[1];
rz(-0.061441378) q[2];
sx q[2];
rz(-0.80020088) q[2];
sx q[2];
rz(2.7903008) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8644789) q[1];
sx q[1];
rz(-0.8833589) q[1];
sx q[1];
rz(0.25798016) q[1];
x q[2];
rz(1.7311312) q[3];
sx q[3];
rz(-2.3162127) q[3];
sx q[3];
rz(-1.3813409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2228955) q[2];
sx q[2];
rz(-1.6632068) q[2];
sx q[2];
rz(-0.29829868) q[2];
rz(-1.9722624) q[3];
sx q[3];
rz(-2.1746641) q[3];
sx q[3];
rz(1.5112618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8339612) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(0.6024012) q[0];
rz(2.7700453) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(-1.0317624) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4435972) q[0];
sx q[0];
rz(-1.5866318) q[0];
sx q[0];
rz(-1.3965194) q[0];
rz(-pi) q[1];
rz(2.736892) q[2];
sx q[2];
rz(-1.3794624) q[2];
sx q[2];
rz(1.9218685) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.885434) q[1];
sx q[1];
rz(-0.51894655) q[1];
sx q[1];
rz(2.7536254) q[1];
rz(-pi) q[2];
rz(0.86182819) q[3];
sx q[3];
rz(-0.78360451) q[3];
sx q[3];
rz(-0.48101048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3370257) q[2];
sx q[2];
rz(-1.1013384) q[2];
sx q[2];
rz(1.6229013) q[2];
rz(-2.2931781) q[3];
sx q[3];
rz(-2.266326) q[3];
sx q[3];
rz(-3.1019822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(1.1161716) q[0];
sx q[0];
rz(-3.0905368) q[0];
sx q[0];
rz(0.99789944) q[0];
rz(-2.8669224) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(1.7707228) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8557963) q[0];
sx q[0];
rz(-1.5573475) q[0];
sx q[0];
rz(0.6310668) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9267124) q[2];
sx q[2];
rz(-1.4265665) q[2];
sx q[2];
rz(-0.027050935) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.22860195) q[1];
sx q[1];
rz(-0.96768846) q[1];
sx q[1];
rz(1.2772389) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79973914) q[3];
sx q[3];
rz(-0.99501401) q[3];
sx q[3];
rz(0.48279253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.85713282) q[2];
sx q[2];
rz(-0.94228116) q[2];
sx q[2];
rz(-0.029646309) q[2];
rz(2.5263785) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(2.7992547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(3.0520332) q[0];
sx q[0];
rz(-2.679306) q[0];
sx q[0];
rz(-0.7830559) q[0];
rz(-1.1198593) q[1];
sx q[1];
rz(-0.66711396) q[1];
sx q[1];
rz(-1.3979744) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8496765) q[0];
sx q[0];
rz(-1.9789985) q[0];
sx q[0];
rz(-2.0565536) q[0];
x q[1];
rz(2.2538848) q[2];
sx q[2];
rz(-2.7811433) q[2];
sx q[2];
rz(0.53149022) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3009792) q[1];
sx q[1];
rz(-0.66794315) q[1];
sx q[1];
rz(1.5850369) q[1];
rz(1.3574202) q[3];
sx q[3];
rz(-2.3436574) q[3];
sx q[3];
rz(-0.40003451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8457501) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(-1.3207377) q[2];
rz(-1.0197506) q[3];
sx q[3];
rz(-1.7176065) q[3];
sx q[3];
rz(-1.7216871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9418697) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(-0.95440188) q[0];
rz(-1.9873387) q[1];
sx q[1];
rz(-2.4955165) q[1];
sx q[1];
rz(-1.0844213) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0823776) q[0];
sx q[0];
rz(-2.431916) q[0];
sx q[0];
rz(1.5060471) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83553548) q[2];
sx q[2];
rz(-0.37097574) q[2];
sx q[2];
rz(2.6717466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0965914) q[1];
sx q[1];
rz(-0.47811478) q[1];
sx q[1];
rz(1.4131312) q[1];
rz(-pi) q[2];
rz(-1.7550657) q[3];
sx q[3];
rz(-0.83067719) q[3];
sx q[3];
rz(-2.4746462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.083100975) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(-2.7044738) q[2];
rz(-1.3433749) q[3];
sx q[3];
rz(-1.790204) q[3];
sx q[3];
rz(-3.0962211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042353543) q[0];
sx q[0];
rz(-2.1411242) q[0];
sx q[0];
rz(-1.4105256) q[0];
rz(0.22077416) q[1];
sx q[1];
rz(-2.6103554) q[1];
sx q[1];
rz(-0.9332307) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2483053) q[0];
sx q[0];
rz(-2.3961005) q[0];
sx q[0];
rz(-3.1214691) q[0];
x q[1];
rz(-0.48671203) q[2];
sx q[2];
rz(-2.463474) q[2];
sx q[2];
rz(0.93914062) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9737967) q[1];
sx q[1];
rz(-2.3575511) q[1];
sx q[1];
rz(-1.7586437) q[1];
x q[2];
rz(0.1118999) q[3];
sx q[3];
rz(-2.2004232) q[3];
sx q[3];
rz(0.96523413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.71020469) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(-2.5101275) q[2];
rz(-2.4225875) q[3];
sx q[3];
rz(-1.3471666) q[3];
sx q[3];
rz(2.9001111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63856335) q[0];
sx q[0];
rz(-1.4894435) q[0];
sx q[0];
rz(1.3110934) q[0];
rz(-0.44166625) q[1];
sx q[1];
rz(-1.3413981) q[1];
sx q[1];
rz(-1.9955019) q[1];
rz(1.491811) q[2];
sx q[2];
rz(-2.7597703) q[2];
sx q[2];
rz(-2.4673354) q[2];
rz(2.9819103) q[3];
sx q[3];
rz(-2.5661712) q[3];
sx q[3];
rz(-2.4721091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-0.85059387) q[0];
sx q[0];
rz(0.2001065) q[0];
rz(2.9432358) q[1];
sx q[1];
rz(-2.618572) q[1];
sx q[1];
rz(1.8097635) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057077335) q[0];
sx q[0];
rz(-2.6383218) q[0];
sx q[0];
rz(-2.0664735) q[0];
x q[1];
rz(2.0609786) q[2];
sx q[2];
rz(-1.0287544) q[2];
sx q[2];
rz(-1.8212183) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7847474) q[1];
sx q[1];
rz(-2.4618399) q[1];
sx q[1];
rz(-3.0650782) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0588714) q[3];
sx q[3];
rz(-1.0961856) q[3];
sx q[3];
rz(-0.079291346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0693822) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(1.3145831) q[2];
rz(-2.2383111) q[3];
sx q[3];
rz(-1.3638146) q[3];
sx q[3];
rz(0.34590736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8697934) q[0];
sx q[0];
rz(-1.6558187) q[0];
sx q[0];
rz(2.2858009) q[0];
rz(2.3675512) q[1];
sx q[1];
rz(-1.9776521) q[1];
sx q[1];
rz(-0.78786293) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5408527) q[0];
sx q[0];
rz(-2.8898281) q[0];
sx q[0];
rz(-1.4815848) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.79736) q[2];
sx q[2];
rz(-1.5308799) q[2];
sx q[2];
rz(1.0913864) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3968421) q[1];
sx q[1];
rz(-1.1424658) q[1];
sx q[1];
rz(0.37439268) q[1];
rz(-0.66707261) q[3];
sx q[3];
rz(-1.7366341) q[3];
sx q[3];
rz(-1.1950243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.065980109) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(0.39609972) q[2];
rz(-1.5710477) q[3];
sx q[3];
rz(-1.5038265) q[3];
sx q[3];
rz(1.7861231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-1.4002607) q[1];
sx q[1];
rz(-0.33904591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.402917) q[0];
sx q[0];
rz(-2.5581723) q[0];
sx q[0];
rz(-0.03761272) q[0];
rz(-1.4475559) q[2];
sx q[2];
rz(-1.3831105) q[2];
sx q[2];
rz(-1.0836243) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0319034) q[1];
sx q[1];
rz(-1.2983783) q[1];
sx q[1];
rz(1.8930356) q[1];
rz(-pi) q[2];
rz(2.364758) q[3];
sx q[3];
rz(-1.5791487) q[3];
sx q[3];
rz(1.2119105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6591349) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(0.47046146) q[2];
rz(-2.2792234) q[3];
sx q[3];
rz(-2.6926398) q[3];
sx q[3];
rz(1.5615777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793133) q[0];
sx q[0];
rz(-2.0231495) q[0];
sx q[0];
rz(0.40801868) q[0];
rz(-1.7904003) q[1];
sx q[1];
rz(-1.8441169) q[1];
sx q[1];
rz(-1.321235) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42962881) q[0];
sx q[0];
rz(-0.59998673) q[0];
sx q[0];
rz(-2.5866051) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2259813) q[2];
sx q[2];
rz(-1.639099) q[2];
sx q[2];
rz(-3.0223522) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.76391782) q[1];
sx q[1];
rz(-1.498068) q[1];
sx q[1];
rz(0.12158981) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9923237) q[3];
sx q[3];
rz(-1.9644004) q[3];
sx q[3];
rz(-0.43962653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2230175) q[2];
sx q[2];
rz(-1.431798) q[2];
sx q[2];
rz(-0.72285405) q[2];
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
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9549114) q[0];
sx q[0];
rz(-2.6246922) q[0];
sx q[0];
rz(0.55996672) q[0];
rz(2.9365183) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(0.29108873) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9331074) q[0];
sx q[0];
rz(-1.4015084) q[0];
sx q[0];
rz(2.8693958) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3423357) q[2];
sx q[2];
rz(-1.5267258) q[2];
sx q[2];
rz(-1.26233) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0138058) q[1];
sx q[1];
rz(-1.3723137) q[1];
sx q[1];
rz(2.2747893) q[1];
rz(-1.7311312) q[3];
sx q[3];
rz(-2.3162127) q[3];
sx q[3];
rz(1.3813409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9186972) q[2];
sx q[2];
rz(-1.6632068) q[2];
sx q[2];
rz(-2.843294) q[2];
rz(1.9722624) q[3];
sx q[3];
rz(-2.1746641) q[3];
sx q[3];
rz(1.6303308) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30763141) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(-2.5391915) q[0];
rz(-0.3715474) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(-1.0317624) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2715797) q[0];
sx q[0];
rz(-1.7450512) q[0];
sx q[0];
rz(-0.016079024) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40470064) q[2];
sx q[2];
rz(-1.3794624) q[2];
sx q[2];
rz(1.9218685) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18373016) q[1];
sx q[1];
rz(-1.0938083) q[1];
sx q[1];
rz(-1.7836003) q[1];
rz(2.2183275) q[3];
sx q[3];
rz(-1.0933236) q[3];
sx q[3];
rz(1.5058194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80456698) q[2];
sx q[2];
rz(-1.1013384) q[2];
sx q[2];
rz(-1.6229013) q[2];
rz(-0.8484146) q[3];
sx q[3];
rz(-2.266326) q[3];
sx q[3];
rz(-0.039610473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.0254211) q[0];
sx q[0];
rz(-0.051055901) q[0];
sx q[0];
rz(-2.1436932) q[0];
rz(2.8669224) q[1];
sx q[1];
rz(-2.2614567) q[1];
sx q[1];
rz(1.7707228) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3034015) q[0];
sx q[0];
rz(-2.5104021) q[0];
sx q[0];
rz(3.1188008) q[0];
x q[1];
rz(1.2148803) q[2];
sx q[2];
rz(-1.7150262) q[2];
sx q[2];
rz(0.027050935) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.26105598) q[1];
sx q[1];
rz(-0.66268259) q[1];
sx q[1];
rz(2.7438394) q[1];
rz(-pi) q[2];
rz(2.3418535) q[3];
sx q[3];
rz(-2.1465786) q[3];
sx q[3];
rz(-0.48279253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85713282) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(-0.029646309) q[2];
rz(0.61521411) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(0.342338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0520332) q[0];
sx q[0];
rz(-0.46228662) q[0];
sx q[0];
rz(2.3585368) q[0];
rz(1.1198593) q[1];
sx q[1];
rz(-0.66711396) q[1];
sx q[1];
rz(1.3979744) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92358855) q[0];
sx q[0];
rz(-2.5177885) q[0];
sx q[0];
rz(-0.82360928) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2538848) q[2];
sx q[2];
rz(-2.7811433) q[2];
sx q[2];
rz(-2.6101024) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.84061342) q[1];
sx q[1];
rz(-0.66794315) q[1];
sx q[1];
rz(1.5565558) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9277738) q[3];
sx q[3];
rz(-0.79584661) q[3];
sx q[3];
rz(2.4406432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8457501) q[2];
sx q[2];
rz(-2.5551899) q[2];
sx q[2];
rz(1.8208549) q[2];
rz(2.1218421) q[3];
sx q[3];
rz(-1.7176065) q[3];
sx q[3];
rz(-1.7216871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9418697) q[0];
sx q[0];
rz(-1.0933192) q[0];
sx q[0];
rz(-2.1871908) q[0];
rz(1.9873387) q[1];
sx q[1];
rz(-0.64607611) q[1];
sx q[1];
rz(-1.0844213) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53756489) q[0];
sx q[0];
rz(-1.6129692) q[0];
sx q[0];
rz(-0.86215638) q[0];
x q[1];
rz(0.25523369) q[2];
sx q[2];
rz(-1.298578) q[2];
sx q[2];
rz(-1.9013426) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3855648) q[1];
sx q[1];
rz(-1.4984908) q[1];
sx q[1];
rz(-1.0977618) q[1];
rz(-1.7550657) q[3];
sx q[3];
rz(-2.3109155) q[3];
sx q[3];
rz(2.4746462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0584917) q[2];
sx q[2];
rz(-0.19187555) q[2];
sx q[2];
rz(2.7044738) q[2];
rz(1.3433749) q[3];
sx q[3];
rz(-1.3513887) q[3];
sx q[3];
rz(-3.0962211) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0992391) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(-1.7310671) q[0];
rz(2.9208185) q[1];
sx q[1];
rz(-2.6103554) q[1];
sx q[1];
rz(-2.208362) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2483053) q[0];
sx q[0];
rz(-0.74549216) q[0];
sx q[0];
rz(3.1214691) q[0];
rz(1.9311229) q[2];
sx q[2];
rz(-0.98305741) q[2];
sx q[2];
rz(2.7992942) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9737967) q[1];
sx q[1];
rz(-2.3575511) q[1];
sx q[1];
rz(1.7586437) q[1];
x q[2];
rz(-1.7228863) q[3];
sx q[3];
rz(-2.503431) q[3];
sx q[3];
rz(1.1537976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.431388) q[2];
sx q[2];
rz(-2.7993671) q[2];
sx q[2];
rz(-2.5101275) q[2];
rz(0.7190052) q[3];
sx q[3];
rz(-1.7944261) q[3];
sx q[3];
rz(-2.9001111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5030293) q[0];
sx q[0];
rz(-1.6521492) q[0];
sx q[0];
rz(-1.8304992) q[0];
rz(-2.6999264) q[1];
sx q[1];
rz(-1.8001945) q[1];
sx q[1];
rz(1.1460907) q[1];
rz(1.9515402) q[2];
sx q[2];
rz(-1.5413918) q[2];
sx q[2];
rz(-0.96985758) q[2];
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

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5496991) q[0];
sx q[0];
rz(-1.1130207) q[0];
sx q[0];
rz(0.23166238) q[0];
rz(0.35056937) q[1];
sx q[1];
rz(-0.21023336) q[1];
sx q[1];
rz(-0.76229131) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0526844) q[0];
sx q[0];
rz(-1.9885953) q[0];
sx q[0];
rz(-0.68939986) q[0];
rz(-pi) q[1];
rz(0.48299787) q[2];
sx q[2];
rz(-1.9873426) q[2];
sx q[2];
rz(2.2838468) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8232441) q[1];
sx q[1];
rz(-1.5155412) q[1];
sx q[1];
rz(-0.18624185) q[1];
x q[2];
rz(0.84299327) q[3];
sx q[3];
rz(-0.61474934) q[3];
sx q[3];
rz(-2.691137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1594557) q[2];
sx q[2];
rz(-2.5348713) q[2];
sx q[2];
rz(-2.7098126) q[2];
rz(2.2744956) q[3];
sx q[3];
rz(-2.1853787) q[3];
sx q[3];
rz(-1.506327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.1113488) q[0];
sx q[0];
rz(-2.3115277) q[0];
sx q[0];
rz(1.9507971) q[0];
rz(1.2442773) q[1];
sx q[1];
rz(-1.1294653) q[1];
sx q[1];
rz(-0.58070374) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55739072) q[0];
sx q[0];
rz(-1.4317436) q[0];
sx q[0];
rz(-1.1625421) q[0];
x q[1];
rz(0.4197311) q[2];
sx q[2];
rz(-2.0629305) q[2];
sx q[2];
rz(2.555814) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3327707) q[1];
sx q[1];
rz(-2.1968107) q[1];
sx q[1];
rz(1.1762192) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45240963) q[3];
sx q[3];
rz(-0.93660855) q[3];
sx q[3];
rz(-0.89434552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1397436) q[2];
sx q[2];
rz(-1.4051583) q[2];
sx q[2];
rz(-1.5025567) q[2];
rz(-0.20770811) q[3];
sx q[3];
rz(-1.8307056) q[3];
sx q[3];
rz(2.1343855) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.945632) q[0];
sx q[0];
rz(-1.0769083) q[0];
sx q[0];
rz(2.9929602) q[0];
rz(-2.0678988) q[1];
sx q[1];
rz(-0.37207347) q[1];
sx q[1];
rz(0.78280848) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1590875) q[0];
sx q[0];
rz(-1.7409865) q[0];
sx q[0];
rz(-3.085932) q[0];
rz(-1.6131079) q[2];
sx q[2];
rz(-2.1033005) q[2];
sx q[2];
rz(-1.29262) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.69082123) q[1];
sx q[1];
rz(-1.3714002) q[1];
sx q[1];
rz(-0.52187199) q[1];
rz(-pi) q[2];
rz(-2.8896096) q[3];
sx q[3];
rz(-2.3445233) q[3];
sx q[3];
rz(2.1261724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17915501) q[2];
sx q[2];
rz(-2.16733) q[2];
sx q[2];
rz(-2.0172737) q[2];
rz(-0.026668523) q[3];
sx q[3];
rz(-2.4494438) q[3];
sx q[3];
rz(2.8587225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.35793316) q[0];
sx q[0];
rz(-1.1555576) q[0];
sx q[0];
rz(-1.5643157) q[0];
rz(2.0951001) q[1];
sx q[1];
rz(-2.107403) q[1];
sx q[1];
rz(1.1276721) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30450102) q[0];
sx q[0];
rz(-2.1620885) q[0];
sx q[0];
rz(1.5832157) q[0];
rz(-pi) q[1];
rz(-1.7972184) q[2];
sx q[2];
rz(-0.34219301) q[2];
sx q[2];
rz(2.7777517) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8279523) q[1];
sx q[1];
rz(-1.8871739) q[1];
sx q[1];
rz(-0.72319855) q[1];
rz(-0.87714449) q[3];
sx q[3];
rz(-2.9877285) q[3];
sx q[3];
rz(-2.2779704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.54101888) q[2];
sx q[2];
rz(-1.270673) q[2];
sx q[2];
rz(2.9533271) q[2];
rz(0.49011928) q[3];
sx q[3];
rz(-1.6907588) q[3];
sx q[3];
rz(-1.274444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1566496) q[0];
sx q[0];
rz(-1.5771414) q[0];
sx q[0];
rz(-1.6749325) q[0];
rz(-2.6690392) q[1];
sx q[1];
rz(-2.2652389) q[1];
sx q[1];
rz(-2.1579425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2547823) q[0];
sx q[0];
rz(-1.3894102) q[0];
sx q[0];
rz(-2.7104293) q[0];
rz(0.46323295) q[2];
sx q[2];
rz(-1.7432799) q[2];
sx q[2];
rz(-1.1197937) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90685788) q[1];
sx q[1];
rz(-0.51381293) q[1];
sx q[1];
rz(-1.4545026) q[1];
x q[2];
rz(-1.9463169) q[3];
sx q[3];
rz(-0.28464475) q[3];
sx q[3];
rz(1.9575021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5094362) q[2];
sx q[2];
rz(-1.7551273) q[2];
sx q[2];
rz(0.40718386) q[2];
rz(1.2716278) q[3];
sx q[3];
rz(-0.41912246) q[3];
sx q[3];
rz(2.4841888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.7206409) q[0];
sx q[0];
rz(-1.9176418) q[0];
sx q[0];
rz(-1.9049013) q[0];
rz(2.0115133) q[1];
sx q[1];
rz(-1.7944444) q[1];
sx q[1];
rz(1.1190266) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41997278) q[0];
sx q[0];
rz(-0.16789745) q[0];
sx q[0];
rz(0.41746087) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79622795) q[2];
sx q[2];
rz(-2.0105462) q[2];
sx q[2];
rz(-1.7481182) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5807075) q[1];
sx q[1];
rz(-1.7153772) q[1];
sx q[1];
rz(2.7156262) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72596999) q[3];
sx q[3];
rz(-0.56258241) q[3];
sx q[3];
rz(1.2364622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59469026) q[2];
sx q[2];
rz(-1.12744) q[2];
sx q[2];
rz(2.8315869) q[2];
rz(1.4340495) q[3];
sx q[3];
rz(-1.6238345) q[3];
sx q[3];
rz(1.2448509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.2731648) q[0];
sx q[0];
rz(-2.7719066) q[0];
sx q[0];
rz(1.0721068) q[0];
rz(1.3605236) q[1];
sx q[1];
rz(-0.37630263) q[1];
sx q[1];
rz(-1.2881813) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5736223) q[0];
sx q[0];
rz(-1.0601655) q[0];
sx q[0];
rz(-0.88167067) q[0];
x q[1];
rz(-0.99432722) q[2];
sx q[2];
rz(-1.948945) q[2];
sx q[2];
rz(-2.7989716) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5406193) q[1];
sx q[1];
rz(-0.29635591) q[1];
sx q[1];
rz(-2.5286416) q[1];
rz(-pi) q[2];
rz(-2.6374972) q[3];
sx q[3];
rz(-2.305784) q[3];
sx q[3];
rz(-3.0127762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0518904) q[2];
sx q[2];
rz(-1.5081729) q[2];
sx q[2];
rz(1.4349597) q[2];
rz(0.47576225) q[3];
sx q[3];
rz(-0.69403726) q[3];
sx q[3];
rz(-0.3835052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39110228) q[0];
sx q[0];
rz(-2.7848211) q[0];
sx q[0];
rz(-1.179689) q[0];
rz(-2.7791595) q[1];
sx q[1];
rz(-2.0795627) q[1];
sx q[1];
rz(1.5550044) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3488551) q[0];
sx q[0];
rz(-0.35733435) q[0];
sx q[0];
rz(-2.6941195) q[0];
rz(-pi) q[1];
rz(2.5767127) q[2];
sx q[2];
rz(-2.3776588) q[2];
sx q[2];
rz(2.8500037) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5240979) q[1];
sx q[1];
rz(-2.3008808) q[1];
sx q[1];
rz(0.9690271) q[1];
rz(2.982627) q[3];
sx q[3];
rz(-1.0086802) q[3];
sx q[3];
rz(-1.6127487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4156551) q[2];
sx q[2];
rz(-1.1154117) q[2];
sx q[2];
rz(1.9647145) q[2];
rz(-2.7365909) q[3];
sx q[3];
rz(-0.98086762) q[3];
sx q[3];
rz(-0.39308959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554351) q[0];
sx q[0];
rz(-0.09849184) q[0];
sx q[0];
rz(1.806102) q[0];
rz(0.91241854) q[1];
sx q[1];
rz(-1.4211979) q[1];
sx q[1];
rz(0.040249912) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9837058) q[0];
sx q[0];
rz(-1.608838) q[0];
sx q[0];
rz(1.4020919) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59800541) q[2];
sx q[2];
rz(-2.0705418) q[2];
sx q[2];
rz(2.6464484) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9736824) q[1];
sx q[1];
rz(-1.0492322) q[1];
sx q[1];
rz(2.6534464) q[1];
x q[2];
rz(2.9222708) q[3];
sx q[3];
rz(-0.78028934) q[3];
sx q[3];
rz(-1.2885476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5361629) q[2];
sx q[2];
rz(-1.9998113) q[2];
sx q[2];
rz(-2.6119168) q[2];
rz(-1.4179432) q[3];
sx q[3];
rz(-2.6726186) q[3];
sx q[3];
rz(-2.6403707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1390851) q[0];
sx q[0];
rz(-1.7739828) q[0];
sx q[0];
rz(3.1097155) q[0];
rz(1.2079976) q[1];
sx q[1];
rz(-0.54787689) q[1];
sx q[1];
rz(-2.8338103) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4144205) q[0];
sx q[0];
rz(-0.57283995) q[0];
sx q[0];
rz(-1.0635873) q[0];
rz(-1.8372907) q[2];
sx q[2];
rz(-1.6937758) q[2];
sx q[2];
rz(2.3044555) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5468502) q[1];
sx q[1];
rz(-2.0786656) q[1];
sx q[1];
rz(-0.87981973) q[1];
rz(0.73460893) q[3];
sx q[3];
rz(-2.9501868) q[3];
sx q[3];
rz(0.020343971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53177437) q[2];
sx q[2];
rz(-2.4419624) q[2];
sx q[2];
rz(3.0979274) q[2];
rz(1.6627436) q[3];
sx q[3];
rz(-2.6378529) q[3];
sx q[3];
rz(-1.7536722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8762348) q[0];
sx q[0];
rz(-1.4690514) q[0];
sx q[0];
rz(1.5734191) q[0];
rz(0.028989446) q[1];
sx q[1];
rz(-2.034076) q[1];
sx q[1];
rz(-2.1709002) q[1];
rz(1.3503648) q[2];
sx q[2];
rz(-1.4511118) q[2];
sx q[2];
rz(1.8076937) q[2];
rz(-0.76282145) q[3];
sx q[3];
rz(-1.9424494) q[3];
sx q[3];
rz(1.8123209) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

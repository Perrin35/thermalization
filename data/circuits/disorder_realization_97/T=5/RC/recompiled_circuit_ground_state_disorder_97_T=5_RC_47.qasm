OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.360541) q[0];
sx q[0];
rz(-0.6147576) q[0];
sx q[0];
rz(-1.0714666) q[0];
rz(2.7331424) q[1];
sx q[1];
rz(-0.93745679) q[1];
sx q[1];
rz(1.4484922) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14872257) q[0];
sx q[0];
rz(-2.5186484) q[0];
sx q[0];
rz(-2.2882266) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57401128) q[2];
sx q[2];
rz(-1.8603627) q[2];
sx q[2];
rz(-2.1776108) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.17098022) q[1];
sx q[1];
rz(-0.17380789) q[1];
sx q[1];
rz(2.5998678) q[1];
rz(0.43464147) q[3];
sx q[3];
rz(-0.9081525) q[3];
sx q[3];
rz(1.545411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.94413269) q[2];
sx q[2];
rz(-2.2977273) q[2];
sx q[2];
rz(2.3485363) q[2];
rz(-0.26886764) q[3];
sx q[3];
rz(-1.3258508) q[3];
sx q[3];
rz(2.1813006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.8307513) q[0];
sx q[0];
rz(-0.38250592) q[0];
sx q[0];
rz(2.261396) q[0];
rz(-0.63086069) q[1];
sx q[1];
rz(-2.3289101) q[1];
sx q[1];
rz(-1.0989443) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55148516) q[0];
sx q[0];
rz(-1.5657358) q[0];
sx q[0];
rz(2.9637815) q[0];
rz(2.9756594) q[2];
sx q[2];
rz(-0.95913061) q[2];
sx q[2];
rz(-0.93390229) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.45881715) q[1];
sx q[1];
rz(-1.4084245) q[1];
sx q[1];
rz(2.008325) q[1];
x q[2];
rz(-0.46862015) q[3];
sx q[3];
rz(-2.1543401) q[3];
sx q[3];
rz(1.1380205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.38689303) q[2];
sx q[2];
rz(-1.6563481) q[2];
sx q[2];
rz(2.5999542) q[2];
rz(0.70332876) q[3];
sx q[3];
rz(-0.43916217) q[3];
sx q[3];
rz(-0.37731236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90848732) q[0];
sx q[0];
rz(-2.2623514) q[0];
sx q[0];
rz(-3.044686) q[0];
rz(-0.58755177) q[1];
sx q[1];
rz(-0.89043003) q[1];
sx q[1];
rz(0.60120916) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5542463) q[0];
sx q[0];
rz(-2.4317435) q[0];
sx q[0];
rz(-0.27940936) q[0];
rz(-pi) q[1];
x q[1];
rz(1.634609) q[2];
sx q[2];
rz(-1.7327961) q[2];
sx q[2];
rz(-0.29159233) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2612348) q[1];
sx q[1];
rz(-0.46030342) q[1];
sx q[1];
rz(0.017437915) q[1];
x q[2];
rz(2.7895687) q[3];
sx q[3];
rz(-0.94402891) q[3];
sx q[3];
rz(-0.75923336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2931557) q[2];
sx q[2];
rz(-1.6907254) q[2];
sx q[2];
rz(2.4270774) q[2];
rz(1.6589288) q[3];
sx q[3];
rz(-2.8667993) q[3];
sx q[3];
rz(2.7627435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.42245427) q[0];
sx q[0];
rz(-1.9062573) q[0];
sx q[0];
rz(-1.5041014) q[0];
rz(1.0927041) q[1];
sx q[1];
rz(-1.2810818) q[1];
sx q[1];
rz(-2.5561996) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8849759) q[0];
sx q[0];
rz(-2.4451849) q[0];
sx q[0];
rz(0.43088669) q[0];
rz(-pi) q[1];
rz(-1.4665514) q[2];
sx q[2];
rz(-1.5603258) q[2];
sx q[2];
rz(0.99828966) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5820739) q[1];
sx q[1];
rz(-0.93739707) q[1];
sx q[1];
rz(-0.60057171) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62111698) q[3];
sx q[3];
rz(-0.74291544) q[3];
sx q[3];
rz(2.7662504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4744711) q[2];
sx q[2];
rz(-0.67255628) q[2];
sx q[2];
rz(-2.7453864) q[2];
rz(2.951156) q[3];
sx q[3];
rz(-0.93947828) q[3];
sx q[3];
rz(-0.52085352) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025295479) q[0];
sx q[0];
rz(-2.681356) q[0];
sx q[0];
rz(0.39475557) q[0];
rz(2.1341628) q[1];
sx q[1];
rz(-1.1578683) q[1];
sx q[1];
rz(-1.6068858) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7028977) q[0];
sx q[0];
rz(-0.83897018) q[0];
sx q[0];
rz(-1.1870866) q[0];
rz(2.7274969) q[2];
sx q[2];
rz(-2.3813435) q[2];
sx q[2];
rz(-0.24157005) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4940961) q[1];
sx q[1];
rz(-0.3682043) q[1];
sx q[1];
rz(-2.1406663) q[1];
x q[2];
rz(2.6768489) q[3];
sx q[3];
rz(-1.9806974) q[3];
sx q[3];
rz(-1.003554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.69636238) q[2];
sx q[2];
rz(-0.69914377) q[2];
sx q[2];
rz(-2.9635079) q[2];
rz(0.94545025) q[3];
sx q[3];
rz(-2.8036696) q[3];
sx q[3];
rz(-0.43615714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6984542) q[0];
sx q[0];
rz(-0.66265023) q[0];
sx q[0];
rz(-0.15289256) q[0];
rz(1.0180417) q[1];
sx q[1];
rz(-0.94296229) q[1];
sx q[1];
rz(-0.61000383) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1069053) q[0];
sx q[0];
rz(-1.6993049) q[0];
sx q[0];
rz(-2.3771108) q[0];
rz(-pi) q[1];
rz(1.583408) q[2];
sx q[2];
rz(-2.5676845) q[2];
sx q[2];
rz(-1.877026) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9023561) q[1];
sx q[1];
rz(-0.78241759) q[1];
sx q[1];
rz(-0.33496015) q[1];
rz(-0.69645564) q[3];
sx q[3];
rz(-3.0142733) q[3];
sx q[3];
rz(0.49139532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6792004) q[2];
sx q[2];
rz(-1.8588763) q[2];
sx q[2];
rz(-2.882615) q[2];
rz(-2.9955043) q[3];
sx q[3];
rz(-0.77272213) q[3];
sx q[3];
rz(-0.63905382) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4884969) q[0];
sx q[0];
rz(-0.86576068) q[0];
sx q[0];
rz(2.3897032) q[0];
rz(-2.409626) q[1];
sx q[1];
rz(-1.7611793) q[1];
sx q[1];
rz(-2.6123349) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89443356) q[0];
sx q[0];
rz(-2.7934847) q[0];
sx q[0];
rz(2.0624119) q[0];
x q[1];
rz(0.13649444) q[2];
sx q[2];
rz(-1.1195099) q[2];
sx q[2];
rz(-0.1642483) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.25594357) q[1];
sx q[1];
rz(-0.082271345) q[1];
sx q[1];
rz(1.72797) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.115836) q[3];
sx q[3];
rz(-1.2833724) q[3];
sx q[3];
rz(-1.5446203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8818714) q[2];
sx q[2];
rz(-1.0292116) q[2];
sx q[2];
rz(2.3465346) q[2];
rz(-1.9184387) q[3];
sx q[3];
rz(-1.3431679) q[3];
sx q[3];
rz(2.3651626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20599468) q[0];
sx q[0];
rz(-1.5933651) q[0];
sx q[0];
rz(-0.11496168) q[0];
rz(-0.089275442) q[1];
sx q[1];
rz(-2.142579) q[1];
sx q[1];
rz(2.9866536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0099338) q[0];
sx q[0];
rz(-1.7007052) q[0];
sx q[0];
rz(-3.0402501) q[0];
rz(2.5903715) q[2];
sx q[2];
rz(-0.9829501) q[2];
sx q[2];
rz(-1.2767222) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3111043) q[1];
sx q[1];
rz(-1.957867) q[1];
sx q[1];
rz(-1.8482659) q[1];
rz(1.6525169) q[3];
sx q[3];
rz(-0.72676672) q[3];
sx q[3];
rz(0.86870199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0838919) q[2];
sx q[2];
rz(-1.8996779) q[2];
sx q[2];
rz(-2.4532301) q[2];
rz(-2.5734731) q[3];
sx q[3];
rz(-2.1391684) q[3];
sx q[3];
rz(0.68305558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4532918) q[0];
sx q[0];
rz(-2.330307) q[0];
sx q[0];
rz(1.2055093) q[0];
rz(-1.0909117) q[1];
sx q[1];
rz(-0.58735192) q[1];
sx q[1];
rz(-1.5376512) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4383008) q[0];
sx q[0];
rz(-2.5592761) q[0];
sx q[0];
rz(-0.87076373) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1411016) q[2];
sx q[2];
rz(-1.0028936) q[2];
sx q[2];
rz(2.9141324) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9133105) q[1];
sx q[1];
rz(-2.3582715) q[1];
sx q[1];
rz(-1.4128774) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0991715) q[3];
sx q[3];
rz(-2.1159008) q[3];
sx q[3];
rz(-0.68600149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9566112) q[2];
sx q[2];
rz(-2.0272171) q[2];
sx q[2];
rz(1.1448917) q[2];
rz(-1.4780809) q[3];
sx q[3];
rz(-0.76668113) q[3];
sx q[3];
rz(1.112282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0612653) q[0];
sx q[0];
rz(-0.63617951) q[0];
sx q[0];
rz(1.7058477) q[0];
rz(-2.0044633) q[1];
sx q[1];
rz(-0.69157332) q[1];
sx q[1];
rz(-0.93708509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81737751) q[0];
sx q[0];
rz(-0.31459168) q[0];
sx q[0];
rz(0.28633519) q[0];
rz(-pi) q[1];
rz(2.6330399) q[2];
sx q[2];
rz(-1.9483231) q[2];
sx q[2];
rz(-3.0005531) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4623066) q[1];
sx q[1];
rz(-2.3704154) q[1];
sx q[1];
rz(0.036840082) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9090462) q[3];
sx q[3];
rz(-2.1595621) q[3];
sx q[3];
rz(-1.3061861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.81637853) q[2];
sx q[2];
rz(-0.64322317) q[2];
sx q[2];
rz(0.31036672) q[2];
rz(-0.54272932) q[3];
sx q[3];
rz(-2.2118745) q[3];
sx q[3];
rz(0.31699666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9515789) q[0];
sx q[0];
rz(-1.0931451) q[0];
sx q[0];
rz(-2.6059294) q[0];
rz(0.71470064) q[1];
sx q[1];
rz(-1.2477881) q[1];
sx q[1];
rz(1.4664149) q[1];
rz(-2.428029) q[2];
sx q[2];
rz(-2.5893671) q[2];
sx q[2];
rz(1.6895369) q[2];
rz(-1.5101931) q[3];
sx q[3];
rz(-2.2987859) q[3];
sx q[3];
rz(1.4684341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

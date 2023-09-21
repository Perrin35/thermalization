OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0918026) q[0];
sx q[0];
rz(-3.0135305) q[0];
sx q[0];
rz(-0.81737104) q[0];
rz(-2.1583537) q[1];
sx q[1];
rz(-2.6020738) q[1];
sx q[1];
rz(-1.9411545) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96785986) q[0];
sx q[0];
rz(-1.9779441) q[0];
sx q[0];
rz(-1.3074387) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3517411) q[2];
sx q[2];
rz(-0.54974216) q[2];
sx q[2];
rz(-1.4866536) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9125036) q[1];
sx q[1];
rz(-2.2664245) q[1];
sx q[1];
rz(1.9181262) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2458385) q[3];
sx q[3];
rz(-0.80245362) q[3];
sx q[3];
rz(-1.6368395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0212705) q[2];
sx q[2];
rz(-2.5734084) q[2];
sx q[2];
rz(-1.5585287) q[2];
rz(2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(-0.42580095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581756) q[0];
sx q[0];
rz(-0.75220627) q[0];
sx q[0];
rz(-3.0875207) q[0];
rz(1.9460829) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(0.53584677) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5862522) q[0];
sx q[0];
rz(-0.99389168) q[0];
sx q[0];
rz(-0.14848407) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70948647) q[2];
sx q[2];
rz(-1.7638793) q[2];
sx q[2];
rz(0.62765861) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0332196) q[1];
sx q[1];
rz(-0.48711005) q[1];
sx q[1];
rz(-0.56652041) q[1];
rz(-1.348098) q[3];
sx q[3];
rz(-2.9950812) q[3];
sx q[3];
rz(-3.1080064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1198931) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(-2.1832441) q[2];
rz(-0.066453233) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7217343) q[0];
sx q[0];
rz(-1.9323213) q[0];
sx q[0];
rz(2.9911175) q[0];
rz(0.45723215) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(-3.1157852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8255071) q[0];
sx q[0];
rz(-1.8206017) q[0];
sx q[0];
rz(0.020629701) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.089695887) q[2];
sx q[2];
rz(-1.8811474) q[2];
sx q[2];
rz(-0.87583625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.60359309) q[1];
sx q[1];
rz(-2.2003761) q[1];
sx q[1];
rz(2.0894719) q[1];
x q[2];
rz(-3.0098626) q[3];
sx q[3];
rz(-1.1909435) q[3];
sx q[3];
rz(2.6704138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1228483) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(-0.5853816) q[2];
rz(2.9600926) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90081763) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(-2.9649819) q[0];
rz(2.2606842) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(-2.6054629) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6648383) q[0];
sx q[0];
rz(-1.9739445) q[0];
sx q[0];
rz(2.320015) q[0];
x q[1];
rz(-1.6279814) q[2];
sx q[2];
rz(-1.8071037) q[2];
sx q[2];
rz(1.1255217) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.94169468) q[1];
sx q[1];
rz(-0.75921339) q[1];
sx q[1];
rz(0.56337507) q[1];
rz(1.2767775) q[3];
sx q[3];
rz(-2.2259568) q[3];
sx q[3];
rz(0.40757195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46999103) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(1.0446576) q[2];
rz(0.70703834) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2351284) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(-2.0902324) q[0];
rz(1.4936739) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(-0.043118127) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5331677) q[0];
sx q[0];
rz(-2.5072917) q[0];
sx q[0];
rz(1.8932635) q[0];
x q[1];
rz(1.0518603) q[2];
sx q[2];
rz(-1.8199925) q[2];
sx q[2];
rz(3.0254226) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9777898) q[1];
sx q[1];
rz(-1.8352574) q[1];
sx q[1];
rz(0.28122854) q[1];
rz(-pi) q[2];
rz(-0.14857265) q[3];
sx q[3];
rz(-0.79509495) q[3];
sx q[3];
rz(0.038392301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.12895) q[2];
sx q[2];
rz(-0.49762112) q[2];
sx q[2];
rz(0.35219231) q[2];
rz(2.5514065) q[3];
sx q[3];
rz(-2.6679109) q[3];
sx q[3];
rz(0.56110704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5181638) q[0];
sx q[0];
rz(-1.2151027) q[0];
sx q[0];
rz(-1.9859001) q[0];
rz(-0.75025264) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(1.0587143) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1832702) q[0];
sx q[0];
rz(-2.000862) q[0];
sx q[0];
rz(1.3413315) q[0];
rz(-pi) q[1];
rz(-0.52144737) q[2];
sx q[2];
rz(-1.8910742) q[2];
sx q[2];
rz(-2.5495868) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2211654) q[1];
sx q[1];
rz(-1.3472918) q[1];
sx q[1];
rz(-1.2616874) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7644464) q[3];
sx q[3];
rz(-0.87246694) q[3];
sx q[3];
rz(-2.8805672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47026149) q[2];
sx q[2];
rz(-1.417421) q[2];
sx q[2];
rz(1.2188101) q[2];
rz(1.1550711) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(1.4412122) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0999775) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(1.6301427) q[0];
rz(-1.3776243) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(0.84164936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11628843) q[0];
sx q[0];
rz(-1.3447273) q[0];
sx q[0];
rz(-0.53209214) q[0];
rz(2.878703) q[2];
sx q[2];
rz(-2.4979696) q[2];
sx q[2];
rz(-0.4292683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0824273) q[1];
sx q[1];
rz(-1.4945684) q[1];
sx q[1];
rz(3.0928897) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3326725) q[3];
sx q[3];
rz(-0.75546414) q[3];
sx q[3];
rz(1.8734224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1402309) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(3.0916396) q[2];
rz(0.66155457) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(0.18930999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426303) q[0];
sx q[0];
rz(-2.1147418) q[0];
sx q[0];
rz(-1.4021953) q[0];
rz(-0.095480355) q[1];
sx q[1];
rz(-1.1663368) q[1];
sx q[1];
rz(-2.7239674) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1091052) q[0];
sx q[0];
rz(-1.0787449) q[0];
sx q[0];
rz(2.0672654) q[0];
rz(-pi) q[1];
rz(-3.0268961) q[2];
sx q[2];
rz(-1.2071929) q[2];
sx q[2];
rz(2.7483658) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8466134) q[1];
sx q[1];
rz(-2.0366922) q[1];
sx q[1];
rz(-2.2094775) q[1];
rz(-pi) q[2];
rz(-2.2291064) q[3];
sx q[3];
rz(-0.45966002) q[3];
sx q[3];
rz(-1.0115136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79545704) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(-2.0987089) q[2];
rz(-2.4677094) q[3];
sx q[3];
rz(-1.6582812) q[3];
sx q[3];
rz(0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8326571) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(-2.5119264) q[0];
rz(-2.5667403) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(-0.94690698) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4862343) q[0];
sx q[0];
rz(-1.2709193) q[0];
sx q[0];
rz(-2.8006058) q[0];
rz(1.6106748) q[2];
sx q[2];
rz(-1.6677742) q[2];
sx q[2];
rz(-0.023488451) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.35518256) q[1];
sx q[1];
rz(-1.2372412) q[1];
sx q[1];
rz(-1.4817609) q[1];
rz(2.9506748) q[3];
sx q[3];
rz(-1.1357765) q[3];
sx q[3];
rz(0.85489475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.56069121) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(-1.2488731) q[2];
rz(-2.4272264) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0666075) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(-2.4865436) q[0];
rz(2.24522) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(0.64430976) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33675942) q[0];
sx q[0];
rz(-1.6258437) q[0];
sx q[0];
rz(2.9677662) q[0];
rz(-0.23004736) q[2];
sx q[2];
rz(-2.2736079) q[2];
sx q[2];
rz(2.7237797) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2584553) q[1];
sx q[1];
rz(-1.9069792) q[1];
sx q[1];
rz(1.8717481) q[1];
rz(1.0993016) q[3];
sx q[3];
rz(-1.1700556) q[3];
sx q[3];
rz(0.61939643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3506938) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(2.541686) q[2];
rz(-0.89896262) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8469289) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(-0.22944336) q[1];
sx q[1];
rz(-2.2506917) q[1];
sx q[1];
rz(-3.0058203) q[1];
rz(-1.7462126) q[2];
sx q[2];
rz(-0.63175628) q[2];
sx q[2];
rz(0.14887688) q[2];
rz(0.89623981) q[3];
sx q[3];
rz(-1.2590209) q[3];
sx q[3];
rz(-0.38929064) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

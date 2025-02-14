OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8371589) q[0];
sx q[0];
rz(-1.2959915) q[0];
sx q[0];
rz(2.0438099) q[0];
rz(-0.9552362) q[1];
sx q[1];
rz(7.0248338) q[1];
sx q[1];
rz(6.4344814) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77028217) q[0];
sx q[0];
rz(-1.6960959) q[0];
sx q[0];
rz(0.025733982) q[0];
rz(0.1163255) q[2];
sx q[2];
rz(-2.1172197) q[2];
sx q[2];
rz(2.8265068) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9469493) q[1];
sx q[1];
rz(-2.2990755) q[1];
sx q[1];
rz(2.3792653) q[1];
x q[2];
rz(1.2829078) q[3];
sx q[3];
rz(-2.3514053) q[3];
sx q[3];
rz(1.1367281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1260881) q[2];
sx q[2];
rz(-1.8512923) q[2];
sx q[2];
rz(2.8224831) q[2];
rz(2.0305521) q[3];
sx q[3];
rz(-2.5824661) q[3];
sx q[3];
rz(-1.3148974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1183209) q[0];
sx q[0];
rz(-2.6225704) q[0];
sx q[0];
rz(-0.58854377) q[0];
rz(2.5449246) q[1];
sx q[1];
rz(-1.8110954) q[1];
sx q[1];
rz(-0.24135022) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5198181) q[0];
sx q[0];
rz(-2.1960917) q[0];
sx q[0];
rz(-0.53859512) q[0];
rz(-pi) q[1];
rz(-1.4407519) q[2];
sx q[2];
rz(-0.50419129) q[2];
sx q[2];
rz(-2.2397704) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6728404) q[1];
sx q[1];
rz(-2.6939055) q[1];
sx q[1];
rz(1.3393558) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7123132) q[3];
sx q[3];
rz(-2.5114473) q[3];
sx q[3];
rz(-2.7272448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1179463) q[2];
sx q[2];
rz(-1.4586552) q[2];
sx q[2];
rz(-3.0493951) q[2];
rz(0.89208952) q[3];
sx q[3];
rz(-2.2690319) q[3];
sx q[3];
rz(-2.0766855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7104765) q[0];
sx q[0];
rz(-1.6350063) q[0];
sx q[0];
rz(0.098467501) q[0];
rz(2.5473728) q[1];
sx q[1];
rz(-2.0829945) q[1];
sx q[1];
rz(-0.99064151) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4020549) q[0];
sx q[0];
rz(-1.4841813) q[0];
sx q[0];
rz(1.9461856) q[0];
rz(-pi) q[1];
rz(2.5712396) q[2];
sx q[2];
rz(-2.2945171) q[2];
sx q[2];
rz(3.0511659) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2222683) q[1];
sx q[1];
rz(-1.3263371) q[1];
sx q[1];
rz(0.9217086) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1414755) q[3];
sx q[3];
rz(-0.50715441) q[3];
sx q[3];
rz(-2.0961268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6777665) q[2];
sx q[2];
rz(-0.78654424) q[2];
sx q[2];
rz(0.80219913) q[2];
rz(1.5944611) q[3];
sx q[3];
rz(-2.0764669) q[3];
sx q[3];
rz(-0.86435634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.7441854) q[0];
sx q[0];
rz(-0.35905251) q[0];
sx q[0];
rz(1.8205951) q[0];
rz(-1.5253223) q[1];
sx q[1];
rz(-2.1883712) q[1];
sx q[1];
rz(2.9811409) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65112075) q[0];
sx q[0];
rz(-1.8118389) q[0];
sx q[0];
rz(-0.81220497) q[0];
rz(-pi) q[1];
rz(1.4549667) q[2];
sx q[2];
rz(-1.7421075) q[2];
sx q[2];
rz(2.0320323) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.077549283) q[1];
sx q[1];
rz(-0.66759118) q[1];
sx q[1];
rz(1.8197219) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64063756) q[3];
sx q[3];
rz(-0.63797073) q[3];
sx q[3];
rz(0.20221329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3204331) q[2];
sx q[2];
rz(-1.4904138) q[2];
sx q[2];
rz(2.5725345) q[2];
rz(-1.9197561) q[3];
sx q[3];
rz(-2.8025083) q[3];
sx q[3];
rz(-0.1327742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4959167) q[0];
sx q[0];
rz(-0.78302947) q[0];
sx q[0];
rz(-2.3107279) q[0];
rz(-1.2105385) q[1];
sx q[1];
rz(-1.4930864) q[1];
sx q[1];
rz(0.99162203) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9308335) q[0];
sx q[0];
rz(-0.88085876) q[0];
sx q[0];
rz(2.7372975) q[0];
rz(0.28056552) q[2];
sx q[2];
rz(-0.90350752) q[2];
sx q[2];
rz(1.9633121) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8960921) q[1];
sx q[1];
rz(-0.84907167) q[1];
sx q[1];
rz(-0.20770276) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6001106) q[3];
sx q[3];
rz(-0.94621822) q[3];
sx q[3];
rz(-1.9196212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7603989) q[2];
sx q[2];
rz(-0.92279592) q[2];
sx q[2];
rz(-2.7823616) q[2];
rz(2.4257816) q[3];
sx q[3];
rz(-0.78407136) q[3];
sx q[3];
rz(-1.1904967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0642218) q[0];
sx q[0];
rz(-1.8697898) q[0];
sx q[0];
rz(0.52128681) q[0];
rz(-2.298666) q[1];
sx q[1];
rz(-1.9631674) q[1];
sx q[1];
rz(-3.0311323) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0960707) q[0];
sx q[0];
rz(-1.318249) q[0];
sx q[0];
rz(-2.4341694) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92824061) q[2];
sx q[2];
rz(-1.5470328) q[2];
sx q[2];
rz(1.0500963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46347324) q[1];
sx q[1];
rz(-1.3729401) q[1];
sx q[1];
rz(-1.0956647) q[1];
x q[2];
rz(-1.0486097) q[3];
sx q[3];
rz(-1.7966101) q[3];
sx q[3];
rz(-1.1831533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.20737401) q[2];
sx q[2];
rz(-1.5997581) q[2];
sx q[2];
rz(0.99384394) q[2];
rz(-0.55073109) q[3];
sx q[3];
rz(-2.6271074) q[3];
sx q[3];
rz(-1.5229092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2391424) q[0];
sx q[0];
rz(-2.7681594) q[0];
sx q[0];
rz(-0.19110876) q[0];
rz(-2.7725819) q[1];
sx q[1];
rz(-1.7419107) q[1];
sx q[1];
rz(-2.5083127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.314069) q[0];
sx q[0];
rz(-1.6681328) q[0];
sx q[0];
rz(0.23989664) q[0];
rz(-pi) q[1];
rz(2.0379637) q[2];
sx q[2];
rz(-0.90641253) q[2];
sx q[2];
rz(0.54069041) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7821473) q[1];
sx q[1];
rz(-1.3656865) q[1];
sx q[1];
rz(-2.9444749) q[1];
x q[2];
rz(2.9464989) q[3];
sx q[3];
rz(-0.70835241) q[3];
sx q[3];
rz(-0.1699902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9621027) q[2];
sx q[2];
rz(-1.7682163) q[2];
sx q[2];
rz(-0.38086677) q[2];
rz(-1.3880091) q[3];
sx q[3];
rz(-1.0264779) q[3];
sx q[3];
rz(1.4306205) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66985828) q[0];
sx q[0];
rz(-0.40100455) q[0];
sx q[0];
rz(1.5420472) q[0];
rz(1.3767287) q[1];
sx q[1];
rz(-1.7574666) q[1];
sx q[1];
rz(-0.9309887) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.72351) q[0];
sx q[0];
rz(-2.2049954) q[0];
sx q[0];
rz(1.4618327) q[0];
x q[1];
rz(0.48540326) q[2];
sx q[2];
rz(-2.6377262) q[2];
sx q[2];
rz(-0.47270838) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7610187) q[1];
sx q[1];
rz(-0.82824681) q[1];
sx q[1];
rz(3.1395802) q[1];
x q[2];
rz(-0.92967195) q[3];
sx q[3];
rz(-2.3768209) q[3];
sx q[3];
rz(1.9509893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5652183) q[2];
sx q[2];
rz(-0.56052506) q[2];
sx q[2];
rz(3.0726688) q[2];
rz(-1.2636412) q[3];
sx q[3];
rz(-2.7694323) q[3];
sx q[3];
rz(-2.4431958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7428335) q[0];
sx q[0];
rz(-2.1837406) q[0];
sx q[0];
rz(3.0897019) q[0];
rz(-0.17008153) q[1];
sx q[1];
rz(-0.64337987) q[1];
sx q[1];
rz(-0.35194078) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55052763) q[0];
sx q[0];
rz(-1.9366486) q[0];
sx q[0];
rz(-1.6195253) q[0];
rz(-2.1743618) q[2];
sx q[2];
rz(-2.1082507) q[2];
sx q[2];
rz(-0.010802566) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.62427795) q[1];
sx q[1];
rz(-2.0950965) q[1];
sx q[1];
rz(2.6844527) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0681814) q[3];
sx q[3];
rz(-1.2651099) q[3];
sx q[3];
rz(2.2524407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6841782) q[2];
sx q[2];
rz(-0.64932051) q[2];
sx q[2];
rz(-2.7049098) q[2];
rz(2.7501578) q[3];
sx q[3];
rz(-1.4623564) q[3];
sx q[3];
rz(-2.2552538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52142414) q[0];
sx q[0];
rz(-2.4346209) q[0];
sx q[0];
rz(0.11908764) q[0];
rz(1.2991615) q[1];
sx q[1];
rz(-1.7487339) q[1];
sx q[1];
rz(1.3577168) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8004476) q[0];
sx q[0];
rz(-0.94609944) q[0];
sx q[0];
rz(0.50826061) q[0];
rz(-pi) q[1];
rz(1.18612) q[2];
sx q[2];
rz(-0.78868491) q[2];
sx q[2];
rz(0.91463156) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.57798856) q[1];
sx q[1];
rz(-1.9050026) q[1];
sx q[1];
rz(-2.5086705) q[1];
rz(0.53472729) q[3];
sx q[3];
rz(-2.2486454) q[3];
sx q[3];
rz(-2.0224366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9327717) q[2];
sx q[2];
rz(-1.5893156) q[2];
sx q[2];
rz(0.61441747) q[2];
rz(-1.6299853) q[3];
sx q[3];
rz(-0.67174086) q[3];
sx q[3];
rz(0.083812788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83196249) q[0];
sx q[0];
rz(-0.93292581) q[0];
sx q[0];
rz(0.25136872) q[0];
rz(2.108719) q[1];
sx q[1];
rz(-1.947247) q[1];
sx q[1];
rz(-1.4364545) q[1];
rz(-1.4658374) q[2];
sx q[2];
rz(-1.5009673) q[2];
sx q[2];
rz(-3.0143723) q[2];
rz(-0.82501199) q[3];
sx q[3];
rz(-0.91559436) q[3];
sx q[3];
rz(0.098500266) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

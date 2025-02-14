OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6239983) q[0];
sx q[0];
rz(-0.52596337) q[0];
sx q[0];
rz(0.21831231) q[0];
rz(1.4766308) q[1];
sx q[1];
rz(-2.6638439) q[1];
sx q[1];
rz(-0.55396095) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.565392) q[0];
sx q[0];
rz(-1.0668653) q[0];
sx q[0];
rz(2.6070057) q[0];
rz(0.10138114) q[2];
sx q[2];
rz(-0.53115618) q[2];
sx q[2];
rz(0.44975933) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2077929) q[1];
sx q[1];
rz(-1.7805532) q[1];
sx q[1];
rz(2.2906474) q[1];
rz(-0.64726909) q[3];
sx q[3];
rz(-2.4993901) q[3];
sx q[3];
rz(-2.6489779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.37609425) q[2];
sx q[2];
rz(-0.29355294) q[2];
sx q[2];
rz(-1.8739088) q[2];
rz(-2.3119161) q[3];
sx q[3];
rz(-1.5209578) q[3];
sx q[3];
rz(0.42384306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053452881) q[0];
sx q[0];
rz(-1.8064878) q[0];
sx q[0];
rz(-1.7470737) q[0];
rz(0.8949737) q[1];
sx q[1];
rz(-1.0286237) q[1];
sx q[1];
rz(-1.5825533) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5415618) q[0];
sx q[0];
rz(-0.91225925) q[0];
sx q[0];
rz(0.70777871) q[0];
x q[1];
rz(-0.87662794) q[2];
sx q[2];
rz(-2.6149984) q[2];
sx q[2];
rz(-1.2132298) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29702052) q[1];
sx q[1];
rz(-0.75648844) q[1];
sx q[1];
rz(-2.6761961) q[1];
x q[2];
rz(-1.1432173) q[3];
sx q[3];
rz(-1.2914011) q[3];
sx q[3];
rz(-2.2277149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.47722694) q[2];
sx q[2];
rz(-1.5218647) q[2];
sx q[2];
rz(3.1108372) q[2];
rz(-0.45298806) q[3];
sx q[3];
rz(-2.9121297) q[3];
sx q[3];
rz(2.9636813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7396616) q[0];
sx q[0];
rz(-0.10098305) q[0];
sx q[0];
rz(0.82823753) q[0];
rz(3.0939057) q[1];
sx q[1];
rz(-0.86367718) q[1];
sx q[1];
rz(-1.9140859) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9102893) q[0];
sx q[0];
rz(-0.72839117) q[0];
sx q[0];
rz(0.69407082) q[0];
x q[1];
rz(-1.1647878) q[2];
sx q[2];
rz(-1.7696524) q[2];
sx q[2];
rz(-2.5420497) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8611698) q[1];
sx q[1];
rz(-2.6211328) q[1];
sx q[1];
rz(-2.1153482) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0892598) q[3];
sx q[3];
rz(-0.84283462) q[3];
sx q[3];
rz(2.6731101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0487655) q[2];
sx q[2];
rz(-2.2266677) q[2];
sx q[2];
rz(-2.0406593) q[2];
rz(-0.01384211) q[3];
sx q[3];
rz(-1.775454) q[3];
sx q[3];
rz(0.83465105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5562627) q[0];
sx q[0];
rz(-1.4034554) q[0];
sx q[0];
rz(2.8072667) q[0];
rz(-0.73257929) q[1];
sx q[1];
rz(-0.94894797) q[1];
sx q[1];
rz(0.11071959) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3667612) q[0];
sx q[0];
rz(-0.80935055) q[0];
sx q[0];
rz(1.5498209) q[0];
rz(-pi) q[1];
rz(-0.51038536) q[2];
sx q[2];
rz(-0.18922999) q[2];
sx q[2];
rz(2.5523579) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1724989) q[1];
sx q[1];
rz(-2.2374472) q[1];
sx q[1];
rz(0.73760017) q[1];
rz(-pi) q[2];
rz(-2.3663051) q[3];
sx q[3];
rz(-1.2167769) q[3];
sx q[3];
rz(-2.275685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.4577786) q[2];
sx q[2];
rz(-2.8595698) q[2];
sx q[2];
rz(1.7337743) q[2];
rz(1.9862566) q[3];
sx q[3];
rz(-1.9852394) q[3];
sx q[3];
rz(-0.19481625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69283501) q[0];
sx q[0];
rz(-2.223707) q[0];
sx q[0];
rz(0.99739972) q[0];
rz(-1.6150486) q[1];
sx q[1];
rz(-0.63957447) q[1];
sx q[1];
rz(-1.4195199) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90516312) q[0];
sx q[0];
rz(-1.4078724) q[0];
sx q[0];
rz(-1.7927367) q[0];
x q[1];
rz(-1.4155343) q[2];
sx q[2];
rz(-2.3387032) q[2];
sx q[2];
rz(-0.93198317) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.43806009) q[1];
sx q[1];
rz(-2.4171481) q[1];
sx q[1];
rz(-2.0634335) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2577293) q[3];
sx q[3];
rz(-3.013844) q[3];
sx q[3];
rz(1.7614438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2898966) q[2];
sx q[2];
rz(-3.0471314) q[2];
sx q[2];
rz(-0.32290253) q[2];
rz(1.1139392) q[3];
sx q[3];
rz(-2.0506004) q[3];
sx q[3];
rz(0.41306257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36393976) q[0];
sx q[0];
rz(-2.8033065) q[0];
sx q[0];
rz(2.3349578) q[0];
rz(0.58397645) q[1];
sx q[1];
rz(-1.1176502) q[1];
sx q[1];
rz(1.1899828) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34790137) q[0];
sx q[0];
rz(-0.37037235) q[0];
sx q[0];
rz(-1.2248951) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3596542) q[2];
sx q[2];
rz(-1.0462772) q[2];
sx q[2];
rz(-3.1280096) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8998225) q[1];
sx q[1];
rz(-0.17647753) q[1];
sx q[1];
rz(1.8679669) q[1];
rz(-pi) q[2];
rz(0.1889008) q[3];
sx q[3];
rz(-2.402225) q[3];
sx q[3];
rz(2.1738571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40818647) q[2];
sx q[2];
rz(-1.5032282) q[2];
sx q[2];
rz(0.1850941) q[2];
rz(-1.5271651) q[3];
sx q[3];
rz(-1.7388758) q[3];
sx q[3];
rz(-2.4534498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1899034) q[0];
sx q[0];
rz(-0.95174319) q[0];
sx q[0];
rz(2.9290747) q[0];
rz(-0.7849794) q[1];
sx q[1];
rz(-1.4947944) q[1];
sx q[1];
rz(0.77883887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3991645) q[0];
sx q[0];
rz(-1.0013097) q[0];
sx q[0];
rz(-2.4692901) q[0];
rz(-pi) q[1];
rz(-0.53063993) q[2];
sx q[2];
rz(-2.4839249) q[2];
sx q[2];
rz(0.69537698) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6952371) q[1];
sx q[1];
rz(-1.8255705) q[1];
sx q[1];
rz(-0.22344113) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1206585) q[3];
sx q[3];
rz(-2.3519197) q[3];
sx q[3];
rz(1.1268953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51398858) q[2];
sx q[2];
rz(-0.49391654) q[2];
sx q[2];
rz(-0.40840515) q[2];
rz(2.4397395) q[3];
sx q[3];
rz(-2.2097094) q[3];
sx q[3];
rz(-1.2592038) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34847611) q[0];
sx q[0];
rz(-1.2154673) q[0];
sx q[0];
rz(-3.075573) q[0];
rz(-1.6325715) q[1];
sx q[1];
rz(-2.0965818) q[1];
sx q[1];
rz(-0.95796934) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97236605) q[0];
sx q[0];
rz(-2.027085) q[0];
sx q[0];
rz(1.2537987) q[0];
rz(-1.0220549) q[2];
sx q[2];
rz(-2.3081452) q[2];
sx q[2];
rz(1.9934143) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9878475) q[1];
sx q[1];
rz(-1.8829131) q[1];
sx q[1];
rz(-2.4448443) q[1];
rz(-pi) q[2];
rz(1.7244206) q[3];
sx q[3];
rz(-1.6699176) q[3];
sx q[3];
rz(-2.2528668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8336739) q[2];
sx q[2];
rz(-1.6180399) q[2];
sx q[2];
rz(-1.7306805) q[2];
rz(2.446567) q[3];
sx q[3];
rz(-1.4796939) q[3];
sx q[3];
rz(-3.1237349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0628292) q[0];
sx q[0];
rz(-1.48209) q[0];
sx q[0];
rz(-0.98989809) q[0];
rz(0.46317378) q[1];
sx q[1];
rz(-1.4521234) q[1];
sx q[1];
rz(-1.90082) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2808229) q[0];
sx q[0];
rz(-1.5131356) q[0];
sx q[0];
rz(-1.8660383) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9389589) q[2];
sx q[2];
rz(-2.0053734) q[2];
sx q[2];
rz(-0.045534924) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3966916) q[1];
sx q[1];
rz(-2.7518175) q[1];
sx q[1];
rz(-2.394963) q[1];
x q[2];
rz(-2.1677446) q[3];
sx q[3];
rz(-0.71248369) q[3];
sx q[3];
rz(-1.2013916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8076294) q[2];
sx q[2];
rz(-1.7619851) q[2];
sx q[2];
rz(0.71869746) q[2];
rz(-1.9888318) q[3];
sx q[3];
rz(-1.4216239) q[3];
sx q[3];
rz(2.1271472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54866791) q[0];
sx q[0];
rz(-0.34807006) q[0];
sx q[0];
rz(0.80192178) q[0];
rz(-1.0879263) q[1];
sx q[1];
rz(-1.3366924) q[1];
sx q[1];
rz(-2.4235639) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1205667) q[0];
sx q[0];
rz(-2.2533198) q[0];
sx q[0];
rz(1.989407) q[0];
x q[1];
rz(2.0640578) q[2];
sx q[2];
rz(-2.0414957) q[2];
sx q[2];
rz(-2.5483709) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9982609) q[1];
sx q[1];
rz(-1.5185396) q[1];
sx q[1];
rz(-1.2388171) q[1];
x q[2];
rz(-0.84661412) q[3];
sx q[3];
rz(-0.9808971) q[3];
sx q[3];
rz(-0.55800948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3512909) q[2];
sx q[2];
rz(-2.9238034) q[2];
sx q[2];
rz(2.6289319) q[2];
rz(-2.3971108) q[3];
sx q[3];
rz(-2.1148041) q[3];
sx q[3];
rz(-2.3775533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-2.9164593) q[0];
sx q[0];
rz(-1.2703348) q[0];
sx q[0];
rz(-1.2402007) q[0];
rz(-0.71612877) q[1];
sx q[1];
rz(-0.54812535) q[1];
sx q[1];
rz(-2.5352238) q[1];
rz(-3.0307583) q[2];
sx q[2];
rz(-1.7946984) q[2];
sx q[2];
rz(2.6738965) q[2];
rz(1.8567139) q[3];
sx q[3];
rz(-2.7919522) q[3];
sx q[3];
rz(-1.5408564) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

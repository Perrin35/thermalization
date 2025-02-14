OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.06935057) q[0];
sx q[0];
rz(-2.0877593) q[0];
sx q[0];
rz(1.8928438) q[0];
rz(1.9034003) q[1];
sx q[1];
rz(-0.44013953) q[1];
sx q[1];
rz(-1.8990489) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0322425) q[0];
sx q[0];
rz(-2.1324131) q[0];
sx q[0];
rz(1.1121145) q[0];
rz(1.6679916) q[2];
sx q[2];
rz(-0.77507918) q[2];
sx q[2];
rz(-0.72778406) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.93852167) q[1];
sx q[1];
rz(-1.2483828) q[1];
sx q[1];
rz(-1.5702269) q[1];
rz(-1.9162769) q[3];
sx q[3];
rz(-0.98967797) q[3];
sx q[3];
rz(1.2863359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0480463) q[2];
sx q[2];
rz(-2.8600433) q[2];
sx q[2];
rz(1.9264889) q[2];
rz(-1.9563227) q[3];
sx q[3];
rz(-0.90648854) q[3];
sx q[3];
rz(-1.5426481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1854061) q[0];
sx q[0];
rz(-2.7423999) q[0];
sx q[0];
rz(1.4584374) q[0];
rz(0.1419119) q[1];
sx q[1];
rz(-1.2071573) q[1];
sx q[1];
rz(1.9292319) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1233009) q[0];
sx q[0];
rz(-0.93138501) q[0];
sx q[0];
rz(-0.42808867) q[0];
rz(-pi) q[1];
rz(-0.0040063695) q[2];
sx q[2];
rz(-0.59213439) q[2];
sx q[2];
rz(2.5483745) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91530529) q[1];
sx q[1];
rz(-2.4608646) q[1];
sx q[1];
rz(-1.4207301) q[1];
x q[2];
rz(-0.47745269) q[3];
sx q[3];
rz(-0.67130723) q[3];
sx q[3];
rz(-0.48760133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.096574664) q[2];
sx q[2];
rz(-0.18288945) q[2];
sx q[2];
rz(0.96192399) q[2];
rz(2.9436881) q[3];
sx q[3];
rz(-1.3062545) q[3];
sx q[3];
rz(-1.643868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(2.4417931) q[0];
sx q[0];
rz(-0.72023359) q[0];
sx q[0];
rz(-2.0752456) q[0];
rz(2.9329246) q[1];
sx q[1];
rz(-2.6228948) q[1];
sx q[1];
rz(0.64812237) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4979734) q[0];
sx q[0];
rz(-2.5622246) q[0];
sx q[0];
rz(1.6672883) q[0];
rz(-3.0096613) q[2];
sx q[2];
rz(-2.6474617) q[2];
sx q[2];
rz(0.021683387) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0433139) q[1];
sx q[1];
rz(-2.2873291) q[1];
sx q[1];
rz(-1.0745506) q[1];
x q[2];
rz(-1.4848292) q[3];
sx q[3];
rz(-0.86557612) q[3];
sx q[3];
rz(-1.7253226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2842747) q[2];
sx q[2];
rz(-1.4633353) q[2];
sx q[2];
rz(-2.5939482) q[2];
rz(-1.0037496) q[3];
sx q[3];
rz(-0.32310969) q[3];
sx q[3];
rz(-2.9799262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.52466398) q[0];
sx q[0];
rz(-1.1268317) q[0];
sx q[0];
rz(-0.3666077) q[0];
rz(1.8560575) q[1];
sx q[1];
rz(-1.6211685) q[1];
sx q[1];
rz(2.3191648) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3503338) q[0];
sx q[0];
rz(-0.21203498) q[0];
sx q[0];
rz(-0.37104443) q[0];
x q[1];
rz(-1.8404615) q[2];
sx q[2];
rz(-1.7326151) q[2];
sx q[2];
rz(-1.9105034) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.3318204) q[1];
sx q[1];
rz(-2.0533516) q[1];
sx q[1];
rz(-2.3259374) q[1];
x q[2];
rz(-3.1347012) q[3];
sx q[3];
rz(-2.8662196) q[3];
sx q[3];
rz(-1.7100348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.6377247) q[2];
sx q[2];
rz(-0.34054264) q[2];
sx q[2];
rz(-1.5857504) q[2];
rz(1.2268892) q[3];
sx q[3];
rz(-2.5033958) q[3];
sx q[3];
rz(1.816412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1894839) q[0];
sx q[0];
rz(-1.1186849) q[0];
sx q[0];
rz(2.191191) q[0];
rz(0.19418007) q[1];
sx q[1];
rz(-1.2545398) q[1];
sx q[1];
rz(-2.8607184) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5820582) q[0];
sx q[0];
rz(-0.96243561) q[0];
sx q[0];
rz(-0.31567659) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5834278) q[2];
sx q[2];
rz(-1.4213741) q[2];
sx q[2];
rz(-2.0858425) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5944526) q[1];
sx q[1];
rz(-3.0277589) q[1];
sx q[1];
rz(0.74233858) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0858795) q[3];
sx q[3];
rz(-0.56929811) q[3];
sx q[3];
rz(1.6804939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6516271) q[2];
sx q[2];
rz(-1.6941035) q[2];
sx q[2];
rz(2.4408686) q[2];
rz(1.0939595) q[3];
sx q[3];
rz(-0.9934727) q[3];
sx q[3];
rz(2.3195364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99264282) q[0];
sx q[0];
rz(-1.0473017) q[0];
sx q[0];
rz(-0.5067504) q[0];
rz(1.1243593) q[1];
sx q[1];
rz(-1.9729112) q[1];
sx q[1];
rz(2.7728424) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4418954) q[0];
sx q[0];
rz(-1.7584193) q[0];
sx q[0];
rz(-2.2806067) q[0];
rz(-pi) q[1];
rz(1.9677591) q[2];
sx q[2];
rz(-2.3572174) q[2];
sx q[2];
rz(-2.337817) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3430887) q[1];
sx q[1];
rz(-1.868416) q[1];
sx q[1];
rz(-0.095620015) q[1];
rz(0.42751023) q[3];
sx q[3];
rz(-1.9906133) q[3];
sx q[3];
rz(1.1838208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2345978) q[2];
sx q[2];
rz(-0.58595053) q[2];
sx q[2];
rz(-2.9242945) q[2];
rz(1.4761188) q[3];
sx q[3];
rz(-0.81513351) q[3];
sx q[3];
rz(-0.34172094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9076964) q[0];
sx q[0];
rz(-0.32231575) q[0];
sx q[0];
rz(-0.79793683) q[0];
rz(-1.4403053) q[1];
sx q[1];
rz(-1.0402352) q[1];
sx q[1];
rz(-2.5247916) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2002522) q[0];
sx q[0];
rz(-1.6464707) q[0];
sx q[0];
rz(-2.5067634) q[0];
rz(0.56767733) q[2];
sx q[2];
rz(-1.857215) q[2];
sx q[2];
rz(-2.84252) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0500461) q[1];
sx q[1];
rz(-1.3669479) q[1];
sx q[1];
rz(2.9943104) q[1];
rz(-pi) q[2];
rz(-2.3921778) q[3];
sx q[3];
rz(-0.85725571) q[3];
sx q[3];
rz(-0.38227043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0133609) q[2];
sx q[2];
rz(-1.1527454) q[2];
sx q[2];
rz(2.4398003) q[2];
rz(0.93891406) q[3];
sx q[3];
rz(-0.89812583) q[3];
sx q[3];
rz(-1.8963337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.668648) q[0];
sx q[0];
rz(-2.9845147) q[0];
sx q[0];
rz(-3.0182086) q[0];
rz(2.3911047) q[1];
sx q[1];
rz(-1.4742943) q[1];
sx q[1];
rz(1.0184681) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21105062) q[0];
sx q[0];
rz(-1.3500542) q[0];
sx q[0];
rz(-0.10737355) q[0];
x q[1];
rz(2.0872714) q[2];
sx q[2];
rz(-1.9811842) q[2];
sx q[2];
rz(-0.058319969) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2606517) q[1];
sx q[1];
rz(-0.48704942) q[1];
sx q[1];
rz(-0.080663514) q[1];
rz(-pi) q[2];
rz(0.94208053) q[3];
sx q[3];
rz(-1.6438369) q[3];
sx q[3];
rz(-1.7353417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5156775) q[2];
sx q[2];
rz(-1.4695784) q[2];
sx q[2];
rz(0.53543004) q[2];
rz(-1.911602) q[3];
sx q[3];
rz(-1.001469) q[3];
sx q[3];
rz(-3.1297704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54315058) q[0];
sx q[0];
rz(-2.4813528) q[0];
sx q[0];
rz(1.2793596) q[0];
rz(1.2641501) q[1];
sx q[1];
rz(-1.8708355) q[1];
sx q[1];
rz(0.15170161) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21148602) q[0];
sx q[0];
rz(-1.1729585) q[0];
sx q[0];
rz(-1.6328638) q[0];
x q[1];
rz(1.5909285) q[2];
sx q[2];
rz(-1.7046456) q[2];
sx q[2];
rz(-1.0878022) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6231767) q[1];
sx q[1];
rz(-2.4881426) q[1];
sx q[1];
rz(-2.099224) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1131198) q[3];
sx q[3];
rz(-2.2185225) q[3];
sx q[3];
rz(-2.5437989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.58574) q[2];
sx q[2];
rz(-2.2057605) q[2];
sx q[2];
rz(1.2305416) q[2];
rz(1.3452283) q[3];
sx q[3];
rz(-1.7788922) q[3];
sx q[3];
rz(1.8432553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1718488) q[0];
sx q[0];
rz(-1.8010704) q[0];
sx q[0];
rz(2.6961683) q[0];
rz(2.716966) q[1];
sx q[1];
rz(-1.5716962) q[1];
sx q[1];
rz(2.5158023) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3836437) q[0];
sx q[0];
rz(-2.1326846) q[0];
sx q[0];
rz(-0.60737078) q[0];
x q[1];
rz(2.9731304) q[2];
sx q[2];
rz(-2.0474153) q[2];
sx q[2];
rz(2.1229975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6546331) q[1];
sx q[1];
rz(-2.5506297) q[1];
sx q[1];
rz(2.0396359) q[1];
rz(-pi) q[2];
rz(-2.8938953) q[3];
sx q[3];
rz(-2.9862635) q[3];
sx q[3];
rz(0.65117902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9618535) q[2];
sx q[2];
rz(-1.0280131) q[2];
sx q[2];
rz(1.6801768) q[2];
rz(0.05750582) q[3];
sx q[3];
rz(-1.4534566) q[3];
sx q[3];
rz(2.1632975) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4089324) q[0];
sx q[0];
rz(-1.4807777) q[0];
sx q[0];
rz(-2.5430191) q[0];
rz(-2.9231425) q[1];
sx q[1];
rz(-1.598806) q[1];
sx q[1];
rz(1.7576408) q[1];
rz(-3.0043428) q[2];
sx q[2];
rz(-0.55036442) q[2];
sx q[2];
rz(-1.3868104) q[2];
rz(-1.5254088) q[3];
sx q[3];
rz(-2.3542891) q[3];
sx q[3];
rz(1.8279984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(1.9873729) q[0];
sx q[0];
rz(4.6648751) q[0];
sx q[0];
rz(9.2821791) q[0];
rz(-2.5481186) q[1];
sx q[1];
rz(-0.65294099) q[1];
sx q[1];
rz(-2.0963734) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3616989) q[0];
sx q[0];
rz(-2.4068659) q[0];
sx q[0];
rz(0.72025062) q[0];
x q[1];
rz(3.1047761) q[2];
sx q[2];
rz(-0.93339257) q[2];
sx q[2];
rz(0.76257818) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3821311) q[1];
sx q[1];
rz(-0.18759094) q[1];
sx q[1];
rz(-2.1407024) q[1];
rz(1.7552283) q[3];
sx q[3];
rz(-0.57653713) q[3];
sx q[3];
rz(0.54631305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.97592252) q[2];
sx q[2];
rz(-1.0079577) q[2];
sx q[2];
rz(-0.21318501) q[2];
rz(0.91607696) q[3];
sx q[3];
rz(-2.7824184) q[3];
sx q[3];
rz(-2.7540414) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85584545) q[0];
sx q[0];
rz(-0.88749945) q[0];
sx q[0];
rz(2.6708653) q[0];
rz(-2.0384906) q[1];
sx q[1];
rz(-2.6328937) q[1];
sx q[1];
rz(-0.57977605) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.267201) q[0];
sx q[0];
rz(-1.5109332) q[0];
sx q[0];
rz(-0.52459985) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9026372) q[2];
sx q[2];
rz(-1.0523426) q[2];
sx q[2];
rz(1.7331074) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5845239) q[1];
sx q[1];
rz(-1.955172) q[1];
sx q[1];
rz(-0.95110431) q[1];
x q[2];
rz(-2.1768119) q[3];
sx q[3];
rz(-0.44108118) q[3];
sx q[3];
rz(1.6562605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9819928) q[2];
sx q[2];
rz(-2.2856568) q[2];
sx q[2];
rz(3.0465872) q[2];
rz(0.57560086) q[3];
sx q[3];
rz(-0.78229457) q[3];
sx q[3];
rz(-1.4466205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59548241) q[0];
sx q[0];
rz(-0.70850104) q[0];
sx q[0];
rz(-2.8693759) q[0];
rz(-3.0369924) q[1];
sx q[1];
rz(-2.101254) q[1];
sx q[1];
rz(-1.4629755) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7143216) q[0];
sx q[0];
rz(-2.2939689) q[0];
sx q[0];
rz(-0.066825213) q[0];
rz(-0.78532312) q[2];
sx q[2];
rz(-1.6162655) q[2];
sx q[2];
rz(-0.77982219) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8045878) q[1];
sx q[1];
rz(-0.92403257) q[1];
sx q[1];
rz(-1.9259364) q[1];
rz(-pi) q[2];
rz(-1.6312863) q[3];
sx q[3];
rz(-1.20245) q[3];
sx q[3];
rz(2.9591536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0649123) q[2];
sx q[2];
rz(-1.400759) q[2];
sx q[2];
rz(2.7317375) q[2];
rz(0.05154933) q[3];
sx q[3];
rz(-0.99286538) q[3];
sx q[3];
rz(-0.73369098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16875295) q[0];
sx q[0];
rz(-2.3641455) q[0];
sx q[0];
rz(0.69946104) q[0];
rz(0.93060023) q[1];
sx q[1];
rz(-1.365463) q[1];
sx q[1];
rz(-1.0994937) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.918952) q[0];
sx q[0];
rz(-1.0066766) q[0];
sx q[0];
rz(0.32911862) q[0];
rz(0.24153562) q[2];
sx q[2];
rz(-2.1061828) q[2];
sx q[2];
rz(-2.8098742) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9605254) q[1];
sx q[1];
rz(-1.2511504) q[1];
sx q[1];
rz(1.245371) q[1];
rz(-pi) q[2];
rz(2.4995125) q[3];
sx q[3];
rz(-1.2386981) q[3];
sx q[3];
rz(-1.737271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3804271) q[2];
sx q[2];
rz(-2.0695504) q[2];
sx q[2];
rz(2.2011444) q[2];
rz(0.56194168) q[3];
sx q[3];
rz(-2.1886531) q[3];
sx q[3];
rz(0.57653069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024427323) q[0];
sx q[0];
rz(-0.086240135) q[0];
sx q[0];
rz(-1.0166136) q[0];
rz(2.2158465) q[1];
sx q[1];
rz(-2.3912906) q[1];
sx q[1];
rz(-3.064916) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5163239) q[0];
sx q[0];
rz(-1.7364572) q[0];
sx q[0];
rz(-1.9693841) q[0];
rz(-pi) q[1];
rz(1.8163278) q[2];
sx q[2];
rz(-1.3708198) q[2];
sx q[2];
rz(2.1918346) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58248427) q[1];
sx q[1];
rz(-0.6000207) q[1];
sx q[1];
rz(-1.1901072) q[1];
rz(-pi) q[2];
rz(-0.29239817) q[3];
sx q[3];
rz(-0.83857036) q[3];
sx q[3];
rz(2.3866619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76806796) q[2];
sx q[2];
rz(-2.3536451) q[2];
sx q[2];
rz(0.12316556) q[2];
rz(-0.67433107) q[3];
sx q[3];
rz(-0.47334039) q[3];
sx q[3];
rz(-0.82790747) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3097565) q[0];
sx q[0];
rz(-2.9530544) q[0];
sx q[0];
rz(-0.48102608) q[0];
rz(-2.9006529) q[1];
sx q[1];
rz(-0.53673458) q[1];
sx q[1];
rz(-0.13490881) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8277934) q[0];
sx q[0];
rz(-2.302357) q[0];
sx q[0];
rz(-0.93904943) q[0];
rz(1.4675702) q[2];
sx q[2];
rz(-2.1035668) q[2];
sx q[2];
rz(-0.68160439) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9489158) q[1];
sx q[1];
rz(-2.2971616) q[1];
sx q[1];
rz(2.7706258) q[1];
rz(-pi) q[2];
rz(1.6429796) q[3];
sx q[3];
rz(-2.0197387) q[3];
sx q[3];
rz(-0.23060966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9350819) q[2];
sx q[2];
rz(-2.2057081) q[2];
sx q[2];
rz(1.8180465) q[2];
rz(0.27213085) q[3];
sx q[3];
rz(-0.85796732) q[3];
sx q[3];
rz(-2.9959397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018933522) q[0];
sx q[0];
rz(-0.7974565) q[0];
sx q[0];
rz(2.4830699) q[0];
rz(-1.5918484) q[1];
sx q[1];
rz(-2.1287983) q[1];
sx q[1];
rz(-0.50643593) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69909912) q[0];
sx q[0];
rz(-0.37537071) q[0];
sx q[0];
rz(-0.18113776) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.030722458) q[2];
sx q[2];
rz(-2.0523768) q[2];
sx q[2];
rz(2.5831163) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5750908) q[1];
sx q[1];
rz(-3.021214) q[1];
sx q[1];
rz(-2.4574404) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6122088) q[3];
sx q[3];
rz(-0.99185252) q[3];
sx q[3];
rz(1.4659297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8062313) q[2];
sx q[2];
rz(-2.3956617) q[2];
sx q[2];
rz(-1.8719505) q[2];
rz(-1.3658203) q[3];
sx q[3];
rz(-3.1277872) q[3];
sx q[3];
rz(-2.5891916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.34614554) q[0];
sx q[0];
rz(-0.2427635) q[0];
sx q[0];
rz(2.7224097) q[0];
rz(-0.090713352) q[1];
sx q[1];
rz(-0.9181298) q[1];
sx q[1];
rz(2.1144287) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8801422) q[0];
sx q[0];
rz(-2.1322726) q[0];
sx q[0];
rz(2.1057157) q[0];
x q[1];
rz(0.33454169) q[2];
sx q[2];
rz(-3.1174264) q[2];
sx q[2];
rz(2.6466359) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.60909437) q[1];
sx q[1];
rz(-0.23819345) q[1];
sx q[1];
rz(2.0373809) q[1];
rz(-pi) q[2];
rz(0.054612463) q[3];
sx q[3];
rz(-2.255618) q[3];
sx q[3];
rz(-0.38338654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1586228) q[2];
sx q[2];
rz(-1.0047487) q[2];
sx q[2];
rz(-0.99009222) q[2];
rz(-0.31791911) q[3];
sx q[3];
rz(-0.81219321) q[3];
sx q[3];
rz(-0.038343553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38744277) q[0];
sx q[0];
rz(-2.4328572) q[0];
sx q[0];
rz(2.5842174) q[0];
rz(0.36224657) q[1];
sx q[1];
rz(-1.7049512) q[1];
sx q[1];
rz(-2.974496) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1633269) q[0];
sx q[0];
rz(-1.3045132) q[0];
sx q[0];
rz(-0.31727088) q[0];
x q[1];
rz(0.2163103) q[2];
sx q[2];
rz(-2.4286793) q[2];
sx q[2];
rz(2.0305433) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7726248) q[1];
sx q[1];
rz(-0.73853044) q[1];
sx q[1];
rz(2.6767297) q[1];
rz(-2.2542893) q[3];
sx q[3];
rz(-1.6287107) q[3];
sx q[3];
rz(1.7660727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1780213) q[2];
sx q[2];
rz(-0.19266291) q[2];
sx q[2];
rz(2.7131405) q[2];
rz(0.7306478) q[3];
sx q[3];
rz(-1.080039) q[3];
sx q[3];
rz(2.54125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927602) q[0];
sx q[0];
rz(-1.4880143) q[0];
sx q[0];
rz(-1.1195419) q[0];
rz(-0.66850942) q[1];
sx q[1];
rz(-2.5709277) q[1];
sx q[1];
rz(2.8817435) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.665312) q[0];
sx q[0];
rz(-2.5710433) q[0];
sx q[0];
rz(2.0670939) q[0];
rz(-pi) q[1];
rz(0.38358263) q[2];
sx q[2];
rz(-1.1337987) q[2];
sx q[2];
rz(2.008977) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47626859) q[1];
sx q[1];
rz(-2.5905053) q[1];
sx q[1];
rz(-2.4399806) q[1];
rz(-pi) q[2];
rz(-1.3803452) q[3];
sx q[3];
rz(-0.69882876) q[3];
sx q[3];
rz(2.8052398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6707637) q[2];
sx q[2];
rz(-1.908611) q[2];
sx q[2];
rz(-1.0864351) q[2];
rz(0.49232617) q[3];
sx q[3];
rz(-0.84572518) q[3];
sx q[3];
rz(-2.4194748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5928741) q[0];
sx q[0];
rz(-1.508779) q[0];
sx q[0];
rz(-1.4887703) q[0];
rz(1.8863574) q[1];
sx q[1];
rz(-1.3175169) q[1];
sx q[1];
rz(-3.0138737) q[1];
rz(-1.9699388) q[2];
sx q[2];
rz(-1.658434) q[2];
sx q[2];
rz(0.077736248) q[2];
rz(0.063379824) q[3];
sx q[3];
rz(-0.90919237) q[3];
sx q[3];
rz(-1.6407871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

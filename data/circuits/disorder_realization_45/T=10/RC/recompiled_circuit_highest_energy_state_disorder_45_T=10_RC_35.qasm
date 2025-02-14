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
rz(1.9396012) q[0];
sx q[0];
rz(-0.48299462) q[0];
sx q[0];
rz(-1.510409) q[0];
rz(3.0022439) q[1];
sx q[1];
rz(-2.5820093) q[1];
sx q[1];
rz(-2.411627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2886292) q[0];
sx q[0];
rz(-1.4242111) q[0];
sx q[0];
rz(0.96340553) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5467186) q[2];
sx q[2];
rz(-1.7584137) q[2];
sx q[2];
rz(1.3485731) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0180359) q[1];
sx q[1];
rz(-2.5606321) q[1];
sx q[1];
rz(0.058079795) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9128051) q[3];
sx q[3];
rz(-1.8513894) q[3];
sx q[3];
rz(-2.5737263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0240747) q[2];
sx q[2];
rz(-2.8705609) q[2];
sx q[2];
rz(0.81895858) q[2];
rz(0.0062746127) q[3];
sx q[3];
rz(-1.2333906) q[3];
sx q[3];
rz(2.1591469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821424) q[0];
sx q[0];
rz(-1.5903951) q[0];
sx q[0];
rz(0.5994125) q[0];
rz(0.86743152) q[1];
sx q[1];
rz(-2.1116833) q[1];
sx q[1];
rz(-0.74554602) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4670938) q[0];
sx q[0];
rz(-1.8601396) q[0];
sx q[0];
rz(-1.4379005) q[0];
x q[1];
rz(0.7600766) q[2];
sx q[2];
rz(-1.3323297) q[2];
sx q[2];
rz(2.5395405) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6709969) q[1];
sx q[1];
rz(-1.7299011) q[1];
sx q[1];
rz(-0.63765031) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50314869) q[3];
sx q[3];
rz(-0.9434349) q[3];
sx q[3];
rz(-2.7706551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.45592371) q[2];
sx q[2];
rz(-2.8440639) q[2];
sx q[2];
rz(-1.4073184) q[2];
rz(-2.2972441) q[3];
sx q[3];
rz(-2.307939) q[3];
sx q[3];
rz(2.1048529) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5582964) q[0];
sx q[0];
rz(-1.946227) q[0];
sx q[0];
rz(-2.1287647) q[0];
rz(0.60802513) q[1];
sx q[1];
rz(-1.5826179) q[1];
sx q[1];
rz(-1.2695405) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8931173) q[0];
sx q[0];
rz(-0.60113827) q[0];
sx q[0];
rz(-2.2918743) q[0];
x q[1];
rz(-1.6744587) q[2];
sx q[2];
rz(-2.4680063) q[2];
sx q[2];
rz(-1.4530593) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7439869) q[1];
sx q[1];
rz(-2.3922046) q[1];
sx q[1];
rz(-0.29008643) q[1];
rz(-pi) q[2];
rz(-1.9023015) q[3];
sx q[3];
rz(-1.1024144) q[3];
sx q[3];
rz(-1.6553594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.515392) q[2];
sx q[2];
rz(-1.0582558) q[2];
sx q[2];
rz(-1.8499648) q[2];
rz(-3.1332704) q[3];
sx q[3];
rz(-2.1885927) q[3];
sx q[3];
rz(-0.21361175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0030647) q[0];
sx q[0];
rz(-0.55532885) q[0];
sx q[0];
rz(1.3767161) q[0];
rz(-1.5273013) q[1];
sx q[1];
rz(-1.9381899) q[1];
sx q[1];
rz(-2.6208904) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6277058) q[0];
sx q[0];
rz(-1.2180274) q[0];
sx q[0];
rz(2.0065432) q[0];
x q[1];
rz(-2.1221913) q[2];
sx q[2];
rz(-1.8810836) q[2];
sx q[2];
rz(0.27675584) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5508073) q[1];
sx q[1];
rz(-2.0001162) q[1];
sx q[1];
rz(2.5270259) q[1];
rz(2.2305829) q[3];
sx q[3];
rz(-1.1413594) q[3];
sx q[3];
rz(2.8991606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.5248096) q[2];
sx q[2];
rz(-2.3526134) q[2];
sx q[2];
rz(1.7899803) q[2];
rz(-2.1221519) q[3];
sx q[3];
rz(-2.5717058) q[3];
sx q[3];
rz(1.1714237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.549642) q[0];
sx q[0];
rz(-1.6787981) q[0];
sx q[0];
rz(0.098966448) q[0];
rz(2.4348266) q[1];
sx q[1];
rz(-2.2711429) q[1];
sx q[1];
rz(0.4695355) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5851368) q[0];
sx q[0];
rz(-0.13988189) q[0];
sx q[0];
rz(-1.3993457) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0355613) q[2];
sx q[2];
rz(-1.240584) q[2];
sx q[2];
rz(0.92452985) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3217056) q[1];
sx q[1];
rz(-0.72607909) q[1];
sx q[1];
rz(-2.5426082) q[1];
rz(-2.585586) q[3];
sx q[3];
rz(-1.4463498) q[3];
sx q[3];
rz(2.5693302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88064757) q[2];
sx q[2];
rz(-1.3346883) q[2];
sx q[2];
rz(-1.3060695) q[2];
rz(2.7169054) q[3];
sx q[3];
rz(-1.7644019) q[3];
sx q[3];
rz(2.2556321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11923085) q[0];
sx q[0];
rz(-2.5194118) q[0];
sx q[0];
rz(-1.8192044) q[0];
rz(-1.127683) q[1];
sx q[1];
rz(-1.6219982) q[1];
sx q[1];
rz(1.7599531) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36395006) q[0];
sx q[0];
rz(-2.3780883) q[0];
sx q[0];
rz(-1.6958773) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3076129) q[2];
sx q[2];
rz(-2.4206941) q[2];
sx q[2];
rz(-0.79254675) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7378714) q[1];
sx q[1];
rz(-0.779169) q[1];
sx q[1];
rz(2.8761151) q[1];
rz(-2.840191) q[3];
sx q[3];
rz(-1.5001138) q[3];
sx q[3];
rz(-1.4247198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3809001) q[2];
sx q[2];
rz(-1.6221294) q[2];
sx q[2];
rz(0.48119989) q[2];
rz(2.6324658) q[3];
sx q[3];
rz(-0.23734084) q[3];
sx q[3];
rz(2.7595162) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5612438) q[0];
sx q[0];
rz(-2.6600397) q[0];
sx q[0];
rz(0.089381889) q[0];
rz(-2.0629758) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(2.1481029) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2334796) q[0];
sx q[0];
rz(-2.3589239) q[0];
sx q[0];
rz(3.0872869) q[0];
rz(2.2682796) q[2];
sx q[2];
rz(-2.2255579) q[2];
sx q[2];
rz(-0.26604929) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1878933) q[1];
sx q[1];
rz(-0.39361289) q[1];
sx q[1];
rz(1.962349) q[1];
rz(-pi) q[2];
rz(2.5105623) q[3];
sx q[3];
rz(-1.1388121) q[3];
sx q[3];
rz(0.61842266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79954687) q[2];
sx q[2];
rz(-0.61508238) q[2];
sx q[2];
rz(-2.7395524) q[2];
rz(-2.0461931) q[3];
sx q[3];
rz(-1.2145372) q[3];
sx q[3];
rz(0.57797617) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47356975) q[0];
sx q[0];
rz(-2.3325925) q[0];
sx q[0];
rz(-1.3336257) q[0];
rz(2.2857621) q[1];
sx q[1];
rz(-1.7355093) q[1];
sx q[1];
rz(-0.91748253) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79865341) q[0];
sx q[0];
rz(-0.22428939) q[0];
sx q[0];
rz(2.0464315) q[0];
rz(-pi) q[1];
rz(-1.7352261) q[2];
sx q[2];
rz(-2.5114369) q[2];
sx q[2];
rz(3.117331) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2189797) q[1];
sx q[1];
rz(-0.41480468) q[1];
sx q[1];
rz(-2.1475683) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4948984) q[3];
sx q[3];
rz(-1.9894674) q[3];
sx q[3];
rz(0.44548098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.63460073) q[2];
sx q[2];
rz(-1.736182) q[2];
sx q[2];
rz(1.290192) q[2];
rz(-2.2318132) q[3];
sx q[3];
rz(-1.6981643) q[3];
sx q[3];
rz(-3.1006052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23583394) q[0];
sx q[0];
rz(-1.9941149) q[0];
sx q[0];
rz(-0.27715096) q[0];
rz(-1.9174891) q[1];
sx q[1];
rz(-1.5382907) q[1];
sx q[1];
rz(-1.7701497) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7200206) q[0];
sx q[0];
rz(-2.4457481) q[0];
sx q[0];
rz(-2.8369342) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0061699) q[2];
sx q[2];
rz(-1.0658385) q[2];
sx q[2];
rz(-1.6944885) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0321694) q[1];
sx q[1];
rz(-0.9309097) q[1];
sx q[1];
rz(-2.8830322) q[1];
rz(-1.7684494) q[3];
sx q[3];
rz(-0.73171333) q[3];
sx q[3];
rz(0.55303516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.203043) q[2];
sx q[2];
rz(-0.86157346) q[2];
sx q[2];
rz(-0.68823632) q[2];
rz(0.31442434) q[3];
sx q[3];
rz(-2.7721072) q[3];
sx q[3];
rz(-1.2615874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37453434) q[0];
sx q[0];
rz(-0.86807591) q[0];
sx q[0];
rz(2.5700997) q[0];
rz(-2.4608965) q[1];
sx q[1];
rz(-1.1704159) q[1];
sx q[1];
rz(2.4933955) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8777953) q[0];
sx q[0];
rz(-1.5476942) q[0];
sx q[0];
rz(3.1277411) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66301753) q[2];
sx q[2];
rz(-0.94919357) q[2];
sx q[2];
rz(0.66689516) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1163016) q[1];
sx q[1];
rz(-2.2696583) q[1];
sx q[1];
rz(-1.541733) q[1];
x q[2];
rz(-2.6246214) q[3];
sx q[3];
rz(-1.7423769) q[3];
sx q[3];
rz(-1.7652546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8615243) q[2];
sx q[2];
rz(-2.0678949) q[2];
sx q[2];
rz(1.6746707) q[2];
rz(0.17659771) q[3];
sx q[3];
rz(-0.64591518) q[3];
sx q[3];
rz(-1.2944029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8661154) q[0];
sx q[0];
rz(-1.7452411) q[0];
sx q[0];
rz(-1.2757975) q[0];
rz(0.38446174) q[1];
sx q[1];
rz(-1.5059595) q[1];
sx q[1];
rz(-0.62677871) q[1];
rz(2.7914417) q[2];
sx q[2];
rz(-1.8684917) q[2];
sx q[2];
rz(-2.6960052) q[2];
rz(0.53388673) q[3];
sx q[3];
rz(-2.0697099) q[3];
sx q[3];
rz(2.8768215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0722421) q[0];
sx q[0];
rz(-1.0538333) q[0];
sx q[0];
rz(-1.8928438) q[0];
rz(-1.2381923) q[1];
sx q[1];
rz(3.5817322) q[1];
sx q[1];
rz(11.323827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1093501) q[0];
sx q[0];
rz(-1.0091796) q[0];
sx q[0];
rz(-1.1121145) q[0];
rz(-0.79808198) q[2];
sx q[2];
rz(-1.6387562) q[2];
sx q[2];
rz(-0.77347212) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.93852167) q[1];
sx q[1];
rz(-1.8932098) q[1];
sx q[1];
rz(1.5702269) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5322024) q[3];
sx q[3];
rz(-1.8577788) q[3];
sx q[3];
rz(0.089394102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0480463) q[2];
sx q[2];
rz(-2.8600433) q[2];
sx q[2];
rz(1.2151037) q[2];
rz(-1.18527) q[3];
sx q[3];
rz(-0.90648854) q[3];
sx q[3];
rz(1.5426481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95618653) q[0];
sx q[0];
rz(-0.39919272) q[0];
sx q[0];
rz(-1.6831552) q[0];
rz(-2.9996808) q[1];
sx q[1];
rz(-1.2071573) q[1];
sx q[1];
rz(1.9292319) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6120554) q[0];
sx q[0];
rz(-0.75242311) q[0];
sx q[0];
rz(2.0799532) q[0];
x q[1];
rz(-0.0040063695) q[2];
sx q[2];
rz(-0.59213439) q[2];
sx q[2];
rz(-0.59321813) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77245678) q[1];
sx q[1];
rz(-1.6650272) q[1];
sx q[1];
rz(-2.2460031) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9208012) q[3];
sx q[3];
rz(-0.98549609) q[3];
sx q[3];
rz(-0.096297527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.096574664) q[2];
sx q[2];
rz(-0.18288945) q[2];
sx q[2];
rz(-0.96192399) q[2];
rz(-2.9436881) q[3];
sx q[3];
rz(-1.3062545) q[3];
sx q[3];
rz(1.643868) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69979954) q[0];
sx q[0];
rz(-0.72023359) q[0];
sx q[0];
rz(-2.0752456) q[0];
rz(-2.9329246) q[1];
sx q[1];
rz(-0.51869789) q[1];
sx q[1];
rz(-2.4934703) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84635669) q[0];
sx q[0];
rz(-1.6235678) q[0];
sx q[0];
rz(2.1480302) q[0];
x q[1];
rz(-0.49049536) q[2];
sx q[2];
rz(-1.5083665) q[2];
sx q[2];
rz(-1.4328116) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0982788) q[1];
sx q[1];
rz(-2.2873291) q[1];
sx q[1];
rz(-2.067042) q[1];
x q[2];
rz(-2.4345458) q[3];
sx q[3];
rz(-1.505369) q[3];
sx q[3];
rz(0.21033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.85731792) q[2];
sx q[2];
rz(-1.4633353) q[2];
sx q[2];
rz(-0.54764444) q[2];
rz(1.0037496) q[3];
sx q[3];
rz(-2.818483) q[3];
sx q[3];
rz(0.16166648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6169287) q[0];
sx q[0];
rz(-1.1268317) q[0];
sx q[0];
rz(-2.774985) q[0];
rz(-1.2855351) q[1];
sx q[1];
rz(-1.5204241) q[1];
sx q[1];
rz(-2.3191648) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4124968) q[0];
sx q[0];
rz(-1.3733882) q[0];
sx q[0];
rz(-1.4928994) q[0];
rz(-pi) q[1];
rz(2.1205495) q[2];
sx q[2];
rz(-2.8281191) q[2];
sx q[2];
rz(-0.18804729) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82717035) q[1];
sx q[1];
rz(-0.91835512) q[1];
sx q[1];
rz(-2.5179522) q[1];
x q[2];
rz(-3.1347012) q[3];
sx q[3];
rz(-2.8662196) q[3];
sx q[3];
rz(1.4315578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.503868) q[2];
sx q[2];
rz(-0.34054264) q[2];
sx q[2];
rz(-1.5558422) q[2];
rz(1.2268892) q[3];
sx q[3];
rz(-0.63819686) q[3];
sx q[3];
rz(-1.816412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9521088) q[0];
sx q[0];
rz(-2.0229078) q[0];
sx q[0];
rz(2.191191) q[0];
rz(2.9474126) q[1];
sx q[1];
rz(-1.8870528) q[1];
sx q[1];
rz(0.28087428) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062894363) q[0];
sx q[0];
rz(-0.67606976) q[0];
sx q[0];
rz(-1.1514787) q[0];
rz(-pi) q[1];
rz(-2.9921586) q[2];
sx q[2];
rz(-1.5583056) q[2];
sx q[2];
rz(0.51692671) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7627613) q[1];
sx q[1];
rz(-1.647659) q[1];
sx q[1];
rz(-3.057544) q[1];
rz(-pi) q[2];
rz(-2.0789541) q[3];
sx q[3];
rz(-1.302037) q[3];
sx q[3];
rz(0.33519676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48996553) q[2];
sx q[2];
rz(-1.4474892) q[2];
sx q[2];
rz(-0.70072407) q[2];
rz(-1.0939595) q[3];
sx q[3];
rz(-2.14812) q[3];
sx q[3];
rz(-0.82205621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1489498) q[0];
sx q[0];
rz(-1.0473017) q[0];
sx q[0];
rz(0.5067504) q[0];
rz(1.1243593) q[1];
sx q[1];
rz(-1.9729112) q[1];
sx q[1];
rz(2.7728424) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4294037) q[0];
sx q[0];
rz(-0.87596873) q[0];
sx q[0];
rz(-2.8963228) q[0];
x q[1];
rz(1.9677591) q[2];
sx q[2];
rz(-0.78437524) q[2];
sx q[2];
rz(2.337817) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0269834) q[1];
sx q[1];
rz(-0.31216808) q[1];
sx q[1];
rz(-1.8725558) q[1];
rz(-pi) q[2];
rz(-1.1147797) q[3];
sx q[3];
rz(-1.9590833) q[3];
sx q[3];
rz(0.20336313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.90699482) q[2];
sx q[2];
rz(-0.58595053) q[2];
sx q[2];
rz(-2.9242945) q[2];
rz(-1.4761188) q[3];
sx q[3];
rz(-0.81513351) q[3];
sx q[3];
rz(0.34172094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9076964) q[0];
sx q[0];
rz(-0.32231575) q[0];
sx q[0];
rz(-2.3436558) q[0];
rz(-1.4403053) q[1];
sx q[1];
rz(-2.1013575) q[1];
sx q[1];
rz(-0.61680102) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9413404) q[0];
sx q[0];
rz(-1.495122) q[0];
sx q[0];
rz(2.5067634) q[0];
x q[1];
rz(-1.2347414) q[2];
sx q[2];
rz(-1.0288887) q[2];
sx q[2];
rz(2.0481179) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7239769) q[1];
sx q[1];
rz(-2.8907052) q[1];
sx q[1];
rz(0.9534568) q[1];
x q[2];
rz(2.2374898) q[3];
sx q[3];
rz(-2.1573503) q[3];
sx q[3];
rz(-1.3400638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1282318) q[2];
sx q[2];
rz(-1.1527454) q[2];
sx q[2];
rz(0.70179233) q[2];
rz(0.93891406) q[3];
sx q[3];
rz(-0.89812583) q[3];
sx q[3];
rz(1.2452589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4729446) q[0];
sx q[0];
rz(-0.15707792) q[0];
sx q[0];
rz(-3.0182086) q[0];
rz(2.3911047) q[1];
sx q[1];
rz(-1.4742943) q[1];
sx q[1];
rz(-2.1231245) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8951851) q[0];
sx q[0];
rz(-2.8965046) q[0];
sx q[0];
rz(-2.0163439) q[0];
rz(-pi) q[1];
rz(-2.6776601) q[2];
sx q[2];
rz(-1.1008769) q[2];
sx q[2];
rz(1.7352833) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5230549) q[1];
sx q[1];
rz(-1.6085165) q[1];
sx q[1];
rz(2.655889) q[1];
x q[2];
rz(-0.94208053) q[3];
sx q[3];
rz(-1.4977557) q[3];
sx q[3];
rz(1.4062509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6259152) q[2];
sx q[2];
rz(-1.6720142) q[2];
sx q[2];
rz(0.53543004) q[2];
rz(1.2299906) q[3];
sx q[3];
rz(-1.001469) q[3];
sx q[3];
rz(-3.1297704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5984421) q[0];
sx q[0];
rz(-0.66023985) q[0];
sx q[0];
rz(1.2793596) q[0];
rz(-1.2641501) q[1];
sx q[1];
rz(-1.2707571) q[1];
sx q[1];
rz(0.15170161) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3705419) q[0];
sx q[0];
rz(-0.40239516) q[0];
sx q[0];
rz(-2.9950525) q[0];
rz(-pi) q[1];
rz(2.9931917) q[2];
sx q[2];
rz(-3.0062468) q[2];
sx q[2];
rz(-0.93805185) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1524304) q[1];
sx q[1];
rz(-2.1235211) q[1];
sx q[1];
rz(2.7731882) q[1];
rz(0.72355481) q[3];
sx q[3];
rz(-1.9949942) q[3];
sx q[3];
rz(1.8198609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.58574) q[2];
sx q[2];
rz(-2.2057605) q[2];
sx q[2];
rz(-1.2305416) q[2];
rz(-1.7963643) q[3];
sx q[3];
rz(-1.7788922) q[3];
sx q[3];
rz(-1.2983373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1718488) q[0];
sx q[0];
rz(-1.8010704) q[0];
sx q[0];
rz(-0.44542435) q[0];
rz(-0.42462665) q[1];
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
rz(1.9741083) q[0];
sx q[0];
rz(-1.0667598) q[0];
sx q[0];
rz(0.91668769) q[0];
x q[1];
rz(1.0883415) q[2];
sx q[2];
rz(-1.4212593) q[2];
sx q[2];
rz(-2.66726) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6596131) q[1];
sx q[1];
rz(-1.3163042) q[1];
sx q[1];
rz(-2.1101497) q[1];
x q[2];
rz(0.24769737) q[3];
sx q[3];
rz(-2.9862635) q[3];
sx q[3];
rz(0.65117902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9618535) q[2];
sx q[2];
rz(-2.1135795) q[2];
sx q[2];
rz(-1.4614159) q[2];
rz(0.05750582) q[3];
sx q[3];
rz(-1.4534566) q[3];
sx q[3];
rz(2.1632975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4089324) q[0];
sx q[0];
rz(-1.6608149) q[0];
sx q[0];
rz(0.59857359) q[0];
rz(2.9231425) q[1];
sx q[1];
rz(-1.5427867) q[1];
sx q[1];
rz(-1.3839518) q[1];
rz(-1.6545532) q[2];
sx q[2];
rz(-1.0261921) q[2];
sx q[2];
rz(1.5941317) q[2];
rz(0.045513734) q[3];
sx q[3];
rz(-0.78452605) q[3];
sx q[3];
rz(1.8922643) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

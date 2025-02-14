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
rz(-1.2381923) q[1];
sx q[1];
rz(3.5817322) q[1];
sx q[1];
rz(11.323827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36168594) q[0];
sx q[0];
rz(-0.70915993) q[0];
sx q[0];
rz(0.61320029) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79808198) q[2];
sx q[2];
rz(-1.5028364) q[2];
sx q[2];
rz(-0.77347212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.203071) q[1];
sx q[1];
rz(-1.2483828) q[1];
sx q[1];
rz(1.5702269) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6655198) q[3];
sx q[3];
rz(-2.4758548) q[3];
sx q[3];
rz(-1.8666603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0480463) q[2];
sx q[2];
rz(-2.8600433) q[2];
sx q[2];
rz(-1.9264889) q[2];
rz(-1.18527) q[3];
sx q[3];
rz(-0.90648854) q[3];
sx q[3];
rz(1.5426481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1854061) q[0];
sx q[0];
rz(-2.7423999) q[0];
sx q[0];
rz(1.6831552) q[0];
rz(-0.1419119) q[1];
sx q[1];
rz(-1.2071573) q[1];
sx q[1];
rz(1.2123607) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6120554) q[0];
sx q[0];
rz(-0.75242311) q[0];
sx q[0];
rz(-1.0616395) q[0];
rz(-pi) q[1];
rz(-1.5734912) q[2];
sx q[2];
rz(-0.97866733) q[2];
sx q[2];
rz(0.58838974) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.77245678) q[1];
sx q[1];
rz(-1.6650272) q[1];
sx q[1];
rz(2.2460031) q[1];
x q[2];
rz(-1.9208012) q[3];
sx q[3];
rz(-0.98549609) q[3];
sx q[3];
rz(-0.096297527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.096574664) q[2];
sx q[2];
rz(-2.9587032) q[2];
sx q[2];
rz(-0.96192399) q[2];
rz(0.19790459) q[3];
sx q[3];
rz(-1.8353381) q[3];
sx q[3];
rz(1.4977247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69979954) q[0];
sx q[0];
rz(-0.72023359) q[0];
sx q[0];
rz(1.066347) q[0];
rz(-0.20866808) q[1];
sx q[1];
rz(-2.6228948) q[1];
sx q[1];
rz(0.64812237) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84635669) q[0];
sx q[0];
rz(-1.5180249) q[0];
sx q[0];
rz(0.99356243) q[0];
rz(-pi) q[1];
rz(-1.5000484) q[2];
sx q[2];
rz(-1.081341) q[2];
sx q[2];
rz(-2.9703028) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3274182) q[1];
sx q[1];
rz(-1.9380373) q[1];
sx q[1];
rz(-0.7805853) q[1];
rz(1.6567635) q[3];
sx q[3];
rz(-0.86557612) q[3];
sx q[3];
rz(1.41627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.85731792) q[2];
sx q[2];
rz(-1.6782574) q[2];
sx q[2];
rz(0.54764444) q[2];
rz(2.137843) q[3];
sx q[3];
rz(-2.818483) q[3];
sx q[3];
rz(-0.16166648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6169287) q[0];
sx q[0];
rz(-1.1268317) q[0];
sx q[0];
rz(-0.3666077) q[0];
rz(-1.8560575) q[1];
sx q[1];
rz(-1.5204241) q[1];
sx q[1];
rz(-0.82242781) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7290958) q[0];
sx q[0];
rz(-1.7682045) q[0];
sx q[0];
rz(1.4928994) q[0];
rz(-pi) q[1];
rz(2.9738178) q[2];
sx q[2];
rz(-1.8368524) q[2];
sx q[2];
rz(-2.7573836) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3144223) q[1];
sx q[1];
rz(-0.91835512) q[1];
sx q[1];
rz(-0.62364044) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1347012) q[3];
sx q[3];
rz(-0.27537307) q[3];
sx q[3];
rz(1.4315578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.503868) q[2];
sx q[2];
rz(-2.80105) q[2];
sx q[2];
rz(-1.5857504) q[2];
rz(1.2268892) q[3];
sx q[3];
rz(-0.63819686) q[3];
sx q[3];
rz(-1.816412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1894839) q[0];
sx q[0];
rz(-2.0229078) q[0];
sx q[0];
rz(-2.191191) q[0];
rz(0.19418007) q[1];
sx q[1];
rz(-1.8870528) q[1];
sx q[1];
rz(-0.28087428) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5595345) q[0];
sx q[0];
rz(-2.179157) q[0];
sx q[0];
rz(0.31567659) q[0];
x q[1];
rz(-3.0578857) q[2];
sx q[2];
rz(-2.9916414) q[2];
sx q[2];
rz(0.97109767) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.54714006) q[1];
sx q[1];
rz(-3.0277589) q[1];
sx q[1];
rz(2.3992541) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0789541) q[3];
sx q[3];
rz(-1.8395556) q[3];
sx q[3];
rz(-2.8063959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.48996553) q[2];
sx q[2];
rz(-1.4474892) q[2];
sx q[2];
rz(-0.70072407) q[2];
rz(2.0476332) q[3];
sx q[3];
rz(-0.9934727) q[3];
sx q[3];
rz(-2.3195364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1489498) q[0];
sx q[0];
rz(-2.094291) q[0];
sx q[0];
rz(2.6348422) q[0];
rz(-1.1243593) q[1];
sx q[1];
rz(-1.9729112) q[1];
sx q[1];
rz(-2.7728424) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6996973) q[0];
sx q[0];
rz(-1.7584193) q[0];
sx q[0];
rz(-0.86098598) q[0];
x q[1];
rz(-2.7733621) q[2];
sx q[2];
rz(-0.86129649) q[2];
sx q[2];
rz(1.3384829) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3974187) q[1];
sx q[1];
rz(-1.4793921) q[1];
sx q[1];
rz(1.2718906) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1147797) q[3];
sx q[3];
rz(-1.9590833) q[3];
sx q[3];
rz(0.20336313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2345978) q[2];
sx q[2];
rz(-0.58595053) q[2];
sx q[2];
rz(-2.9242945) q[2];
rz(-1.4761188) q[3];
sx q[3];
rz(-2.3264591) q[3];
sx q[3];
rz(2.7998717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23389626) q[0];
sx q[0];
rz(-2.8192769) q[0];
sx q[0];
rz(-2.3436558) q[0];
rz(-1.4403053) q[1];
sx q[1];
rz(-2.1013575) q[1];
sx q[1];
rz(-0.61680102) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2002522) q[0];
sx q[0];
rz(-1.495122) q[0];
sx q[0];
rz(2.5067634) q[0];
x q[1];
rz(1.2347414) q[2];
sx q[2];
rz(-1.0288887) q[2];
sx q[2];
rz(1.0934747) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.7239769) q[1];
sx q[1];
rz(-2.8907052) q[1];
sx q[1];
rz(-0.9534568) q[1];
x q[2];
rz(0.7019667) q[3];
sx q[3];
rz(-2.1118374) q[3];
sx q[3];
rz(-0.64149414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1282318) q[2];
sx q[2];
rz(-1.9888473) q[2];
sx q[2];
rz(-0.70179233) q[2];
rz(2.2026786) q[3];
sx q[3];
rz(-0.89812583) q[3];
sx q[3];
rz(-1.2452589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4729446) q[0];
sx q[0];
rz(-2.9845147) q[0];
sx q[0];
rz(0.12338403) q[0];
rz(2.3911047) q[1];
sx q[1];
rz(-1.6672983) q[1];
sx q[1];
rz(2.1231245) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8951851) q[0];
sx q[0];
rz(-0.24508805) q[0];
sx q[0];
rz(1.1252488) q[0];
rz(2.0872714) q[2];
sx q[2];
rz(-1.1604084) q[2];
sx q[2];
rz(0.058319969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88094096) q[1];
sx q[1];
rz(-0.48704942) q[1];
sx q[1];
rz(0.080663514) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4470149) q[3];
sx q[3];
rz(-2.5092193) q[3];
sx q[3];
rz(0.06452175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.6259152) q[2];
sx q[2];
rz(-1.4695784) q[2];
sx q[2];
rz(-0.53543004) q[2];
rz(-1.2299906) q[3];
sx q[3];
rz(-2.1401236) q[3];
sx q[3];
rz(0.011822239) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54315058) q[0];
sx q[0];
rz(-0.66023985) q[0];
sx q[0];
rz(-1.2793596) q[0];
rz(-1.8774425) q[1];
sx q[1];
rz(-1.2707571) q[1];
sx q[1];
rz(2.989891) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21148602) q[0];
sx q[0];
rz(-1.1729585) q[0];
sx q[0];
rz(-1.6328638) q[0];
rz(-1.5506641) q[2];
sx q[2];
rz(-1.7046456) q[2];
sx q[2];
rz(-1.0878022) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.38167414) q[1];
sx q[1];
rz(-1.2592788) q[1];
sx q[1];
rz(-0.98656922) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4180378) q[3];
sx q[3];
rz(-1.1465985) q[3];
sx q[3];
rz(1.3217317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.58574) q[2];
sx q[2];
rz(-2.2057605) q[2];
sx q[2];
rz(-1.2305416) q[2];
rz(-1.3452283) q[3];
sx q[3];
rz(-1.3627005) q[3];
sx q[3];
rz(1.8432553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96974385) q[0];
sx q[0];
rz(-1.8010704) q[0];
sx q[0];
rz(-0.44542435) q[0];
rz(-2.716966) q[1];
sx q[1];
rz(-1.5698965) q[1];
sx q[1];
rz(2.5158023) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9741083) q[0];
sx q[0];
rz(-1.0667598) q[0];
sx q[0];
rz(0.91668769) q[0];
rz(-pi) q[1];
rz(1.0883415) q[2];
sx q[2];
rz(-1.7203334) q[2];
sx q[2];
rz(-0.47433269) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48695951) q[1];
sx q[1];
rz(-2.5506297) q[1];
sx q[1];
rz(-1.1019568) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15066093) q[3];
sx q[3];
rz(-1.5328578) q[3];
sx q[3];
rz(-2.466809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9618535) q[2];
sx q[2];
rz(-2.1135795) q[2];
sx q[2];
rz(1.6801768) q[2];
rz(-0.05750582) q[3];
sx q[3];
rz(-1.688136) q[3];
sx q[3];
rz(2.1632975) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4089324) q[0];
sx q[0];
rz(-1.4807777) q[0];
sx q[0];
rz(-2.5430191) q[0];
rz(2.9231425) q[1];
sx q[1];
rz(-1.5427867) q[1];
sx q[1];
rz(-1.3839518) q[1];
rz(1.6545532) q[2];
sx q[2];
rz(-2.1154006) q[2];
sx q[2];
rz(-1.5474609) q[2];
rz(2.3575847) q[3];
sx q[3];
rz(-1.6029458) q[3];
sx q[3];
rz(-2.8523469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

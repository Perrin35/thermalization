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
rz(2.8528557) q[0];
sx q[0];
rz(-0.69536916) q[0];
sx q[0];
rz(-2.8737972) q[0];
rz(-2.7195622) q[1];
sx q[1];
rz(-0.92075092) q[1];
sx q[1];
rz(1.86778) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0604297) q[0];
sx q[0];
rz(-1.8120684) q[0];
sx q[0];
rz(-3.0953498) q[0];
x q[1];
rz(2.535475) q[2];
sx q[2];
rz(-0.96783468) q[2];
sx q[2];
rz(-0.22113344) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5112524) q[1];
sx q[1];
rz(-2.4773438) q[1];
sx q[1];
rz(-3.0306007) q[1];
rz(-1.4273604) q[3];
sx q[3];
rz(-1.6063476) q[3];
sx q[3];
rz(-1.8604904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5219118) q[2];
sx q[2];
rz(-0.64138594) q[2];
sx q[2];
rz(3.0989975) q[2];
rz(-0.28111449) q[3];
sx q[3];
rz(-1.5646076) q[3];
sx q[3];
rz(-0.59578305) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6302781) q[0];
sx q[0];
rz(-1.9673286) q[0];
sx q[0];
rz(-1.0761155) q[0];
rz(-2.1108744) q[1];
sx q[1];
rz(-1.1117671) q[1];
sx q[1];
rz(1.8310742) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1640747) q[0];
sx q[0];
rz(-2.4830472) q[0];
sx q[0];
rz(-1.5367299) q[0];
x q[1];
rz(0.4035455) q[2];
sx q[2];
rz(-2.2650044) q[2];
sx q[2];
rz(-2.780404) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.11721048) q[1];
sx q[1];
rz(-1.6892994) q[1];
sx q[1];
rz(-0.16853965) q[1];
rz(-pi) q[2];
rz(0.53383975) q[3];
sx q[3];
rz(-1.9263785) q[3];
sx q[3];
rz(-1.3188286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.217546) q[2];
sx q[2];
rz(-2.3895538) q[2];
sx q[2];
rz(-2.1853866) q[2];
rz(-1.8337967) q[3];
sx q[3];
rz(-2.3266413) q[3];
sx q[3];
rz(-3.0202878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2995375) q[0];
sx q[0];
rz(-2.6326023) q[0];
sx q[0];
rz(0.036238413) q[0];
rz(2.3402975) q[1];
sx q[1];
rz(-1.5943297) q[1];
sx q[1];
rz(-1.8410199) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0220003) q[0];
sx q[0];
rz(-1.4109857) q[0];
sx q[0];
rz(0.10800604) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6341798) q[2];
sx q[2];
rz(-0.59813872) q[2];
sx q[2];
rz(1.6063521) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50967783) q[1];
sx q[1];
rz(-1.011112) q[1];
sx q[1];
rz(2.6259929) q[1];
rz(2.4518029) q[3];
sx q[3];
rz(-1.2827323) q[3];
sx q[3];
rz(-0.40460872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7654045) q[2];
sx q[2];
rz(-1.7690965) q[2];
sx q[2];
rz(-2.9236531) q[2];
rz(-0.91207063) q[3];
sx q[3];
rz(-1.6488766) q[3];
sx q[3];
rz(0.85165858) q[3];
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
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7243778) q[0];
sx q[0];
rz(-2.0700924) q[0];
sx q[0];
rz(0.81556129) q[0];
rz(1.0527481) q[1];
sx q[1];
rz(-2.2595854) q[1];
sx q[1];
rz(1.7900593) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6725051) q[0];
sx q[0];
rz(-2.6658305) q[0];
sx q[0];
rz(-2.9722054) q[0];
rz(-pi) q[1];
rz(-0.10064023) q[2];
sx q[2];
rz(-2.422524) q[2];
sx q[2];
rz(1.4532064) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7564769) q[1];
sx q[1];
rz(-0.47620344) q[1];
sx q[1];
rz(-1.508273) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4249486) q[3];
sx q[3];
rz(-1.9633246) q[3];
sx q[3];
rz(1.6300622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1286596) q[2];
sx q[2];
rz(-1.9482875) q[2];
sx q[2];
rz(2.7093757) q[2];
rz(2.2991119) q[3];
sx q[3];
rz(-2.1079) q[3];
sx q[3];
rz(-0.99561083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1603482) q[0];
sx q[0];
rz(-0.46537414) q[0];
sx q[0];
rz(-2.5904742) q[0];
rz(0.61800686) q[1];
sx q[1];
rz(-1.8564686) q[1];
sx q[1];
rz(2.0726223) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38902125) q[0];
sx q[0];
rz(-2.4271936) q[0];
sx q[0];
rz(-0.32984372) q[0];
x q[1];
rz(-2.6591797) q[2];
sx q[2];
rz(-1.97792) q[2];
sx q[2];
rz(0.77821748) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.080090962) q[1];
sx q[1];
rz(-1.6964165) q[1];
sx q[1];
rz(-0.83275034) q[1];
rz(-pi) q[2];
rz(0.25311796) q[3];
sx q[3];
rz(-2.6147644) q[3];
sx q[3];
rz(1.4344858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3507639) q[2];
sx q[2];
rz(-0.5841693) q[2];
sx q[2];
rz(-2.186415) q[2];
rz(-2.5480934) q[3];
sx q[3];
rz(-2.1953526) q[3];
sx q[3];
rz(1.3957297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7668358) q[0];
sx q[0];
rz(-0.042348472) q[0];
sx q[0];
rz(-1.391885) q[0];
rz(-1.1426686) q[1];
sx q[1];
rz(-1.3418158) q[1];
sx q[1];
rz(-1.3124189) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56579933) q[0];
sx q[0];
rz(-0.84747073) q[0];
sx q[0];
rz(-2.2295843) q[0];
rz(-pi) q[1];
rz(2.4644971) q[2];
sx q[2];
rz(-2.1506718) q[2];
sx q[2];
rz(3.0486272) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9095972) q[1];
sx q[1];
rz(-2.084752) q[1];
sx q[1];
rz(1.5088874) q[1];
x q[2];
rz(1.7465318) q[3];
sx q[3];
rz(-2.8214957) q[3];
sx q[3];
rz(-1.2552346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.352508) q[2];
sx q[2];
rz(-2.3859873) q[2];
sx q[2];
rz(-1.7279203) q[2];
rz(0.46755725) q[3];
sx q[3];
rz(-0.65842015) q[3];
sx q[3];
rz(2.5889034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6308052) q[0];
sx q[0];
rz(-3.0785705) q[0];
sx q[0];
rz(-0.68156534) q[0];
rz(1.956578) q[1];
sx q[1];
rz(-2.3763035) q[1];
sx q[1];
rz(-0.78651816) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.701292) q[0];
sx q[0];
rz(-1.1806187) q[0];
sx q[0];
rz(-1.2083645) q[0];
rz(-0.925073) q[2];
sx q[2];
rz(-0.73410119) q[2];
sx q[2];
rz(2.2112598) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9481755) q[1];
sx q[1];
rz(-1.6778291) q[1];
sx q[1];
rz(0.89148895) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.041140838) q[3];
sx q[3];
rz(-2.3015071) q[3];
sx q[3];
rz(-1.5276437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6098392) q[2];
sx q[2];
rz(-0.24007758) q[2];
sx q[2];
rz(1.5698203) q[2];
rz(1.8543367) q[3];
sx q[3];
rz(-1.6183034) q[3];
sx q[3];
rz(2.1985445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59931961) q[0];
sx q[0];
rz(-2.4241408) q[0];
sx q[0];
rz(1.188311) q[0];
rz(0.020847281) q[1];
sx q[1];
rz(-2.3884845) q[1];
sx q[1];
rz(-2.5426224) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1513903) q[0];
sx q[0];
rz(-1.6693475) q[0];
sx q[0];
rz(-1.0644142) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.125419) q[2];
sx q[2];
rz(-2.3756304) q[2];
sx q[2];
rz(0.19419032) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7341344) q[1];
sx q[1];
rz(-2.2770513) q[1];
sx q[1];
rz(0.64532537) q[1];
rz(-pi) q[2];
rz(2.0557559) q[3];
sx q[3];
rz(-0.79130486) q[3];
sx q[3];
rz(2.6216061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.282436) q[2];
sx q[2];
rz(-1.891529) q[2];
sx q[2];
rz(2.6263728) q[2];
rz(-2.6089) q[3];
sx q[3];
rz(-2.0566514) q[3];
sx q[3];
rz(2.2044619) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027997967) q[0];
sx q[0];
rz(-2.4585215) q[0];
sx q[0];
rz(2.7711476) q[0];
rz(-1.0796374) q[1];
sx q[1];
rz(-2.7307983) q[1];
sx q[1];
rz(0.058578514) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48029362) q[0];
sx q[0];
rz(-0.96503557) q[0];
sx q[0];
rz(-1.3783216) q[0];
rz(-pi) q[1];
rz(-1.7981766) q[2];
sx q[2];
rz(-1.6140964) q[2];
sx q[2];
rz(-0.82466489) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6296719) q[1];
sx q[1];
rz(-1.6663934) q[1];
sx q[1];
rz(-0.78533502) q[1];
x q[2];
rz(2.305916) q[3];
sx q[3];
rz(-2.0234851) q[3];
sx q[3];
rz(-1.9122461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8687245) q[2];
sx q[2];
rz(-2.3385907) q[2];
sx q[2];
rz(2.8105984) q[2];
rz(0.013414772) q[3];
sx q[3];
rz(-0.9534854) q[3];
sx q[3];
rz(2.0851871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5633504) q[0];
sx q[0];
rz(-0.42911068) q[0];
sx q[0];
rz(-0.44813928) q[0];
rz(3.0478802) q[1];
sx q[1];
rz(-0.30148503) q[1];
sx q[1];
rz(1.9261446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96998065) q[0];
sx q[0];
rz(-1.6002065) q[0];
sx q[0];
rz(-2.0038811) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9860259) q[2];
sx q[2];
rz(-2.3390963) q[2];
sx q[2];
rz(-1.8181096) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3642973) q[1];
sx q[1];
rz(-2.5357995) q[1];
sx q[1];
rz(-2.0531027) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13588174) q[3];
sx q[3];
rz(-1.1757506) q[3];
sx q[3];
rz(-0.53077936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1756246) q[2];
sx q[2];
rz(-1.5529239) q[2];
sx q[2];
rz(2.442181) q[2];
rz(-1.2753963) q[3];
sx q[3];
rz(-1.9857152) q[3];
sx q[3];
rz(1.8705961) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70915837) q[0];
sx q[0];
rz(-2.2812738) q[0];
sx q[0];
rz(-1.9914837) q[0];
rz(1.3954096) q[1];
sx q[1];
rz(-1.9231053) q[1];
sx q[1];
rz(-1.3833192) q[1];
rz(-1.4906314) q[2];
sx q[2];
rz(-0.90654984) q[2];
sx q[2];
rz(-0.64150099) q[2];
rz(-1.4518573) q[3];
sx q[3];
rz(-2.0343067) q[3];
sx q[3];
rz(0.35342356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

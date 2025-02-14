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
rz(0.68501002) q[0];
sx q[0];
rz(-2.4497439) q[0];
sx q[0];
rz(0.43842167) q[0];
rz(3.0216079) q[1];
sx q[1];
rz(-2.9000403) q[1];
sx q[1];
rz(-0.99689364) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6357898) q[0];
sx q[0];
rz(-1.6576901) q[0];
sx q[0];
rz(-1.4181697) q[0];
x q[1];
rz(0.27049944) q[2];
sx q[2];
rz(-1.5322313) q[2];
sx q[2];
rz(-2.6730516) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4516075) q[1];
sx q[1];
rz(-1.3295191) q[1];
sx q[1];
rz(1.1651768) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71019594) q[3];
sx q[3];
rz(-2.8330148) q[3];
sx q[3];
rz(0.86581826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3991656) q[2];
sx q[2];
rz(-0.70968598) q[2];
sx q[2];
rz(-0.7134552) q[2];
rz(-2.3136638) q[3];
sx q[3];
rz(-1.6498339) q[3];
sx q[3];
rz(2.2152065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4543318) q[0];
sx q[0];
rz(-2.5255272) q[0];
sx q[0];
rz(-1.2607505) q[0];
rz(-2.8997391) q[1];
sx q[1];
rz(-1.8459903) q[1];
sx q[1];
rz(2.5557925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59585786) q[0];
sx q[0];
rz(-2.7867466) q[0];
sx q[0];
rz(-2.0442346) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5795779) q[2];
sx q[2];
rz(-1.9987371) q[2];
sx q[2];
rz(1.2195171) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1928382) q[1];
sx q[1];
rz(-0.84872972) q[1];
sx q[1];
rz(-2.7428736) q[1];
x q[2];
rz(1.582566) q[3];
sx q[3];
rz(-1.571072) q[3];
sx q[3];
rz(-2.2451412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5213617) q[2];
sx q[2];
rz(-0.62180454) q[2];
sx q[2];
rz(-2.8561031) q[2];
rz(-0.96389687) q[3];
sx q[3];
rz(-0.6670835) q[3];
sx q[3];
rz(0.17531659) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6730839) q[0];
sx q[0];
rz(-2.5991169) q[0];
sx q[0];
rz(0.28955224) q[0];
rz(0.30329224) q[1];
sx q[1];
rz(-1.2723609) q[1];
sx q[1];
rz(-1.2264651) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2305812) q[0];
sx q[0];
rz(-2.2104461) q[0];
sx q[0];
rz(0.53627642) q[0];
rz(1.8200335) q[2];
sx q[2];
rz(-2.1366883) q[2];
sx q[2];
rz(1.6911473) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.94680601) q[1];
sx q[1];
rz(-0.99954075) q[1];
sx q[1];
rz(-2.8033546) q[1];
x q[2];
rz(-1.4530334) q[3];
sx q[3];
rz(-2.3533245) q[3];
sx q[3];
rz(0.39358172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.44846416) q[2];
sx q[2];
rz(-2.144564) q[2];
sx q[2];
rz(2.2779951) q[2];
rz(-3.038285) q[3];
sx q[3];
rz(-1.3477252) q[3];
sx q[3];
rz(3.0062655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0057664) q[0];
sx q[0];
rz(-0.62060452) q[0];
sx q[0];
rz(-0.045106877) q[0];
rz(0.82334423) q[1];
sx q[1];
rz(-2.7594559) q[1];
sx q[1];
rz(2.8270922) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.088011) q[0];
sx q[0];
rz(-1.4865666) q[0];
sx q[0];
rz(0.65902577) q[0];
x q[1];
rz(-0.41750022) q[2];
sx q[2];
rz(-2.0682149) q[2];
sx q[2];
rz(-2.6365947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.38344774) q[1];
sx q[1];
rz(-0.52844344) q[1];
sx q[1];
rz(-0.37208884) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32415402) q[3];
sx q[3];
rz(-1.7739033) q[3];
sx q[3];
rz(2.5517011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0663674) q[2];
sx q[2];
rz(-1.6547357) q[2];
sx q[2];
rz(2.5328947) q[2];
rz(1.4501976) q[3];
sx q[3];
rz(-0.83289731) q[3];
sx q[3];
rz(-0.47877413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64959127) q[0];
sx q[0];
rz(-2.2719125) q[0];
sx q[0];
rz(-0.087652303) q[0];
rz(-1.2610669) q[1];
sx q[1];
rz(-1.4050211) q[1];
sx q[1];
rz(2.2671949) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1028624) q[0];
sx q[0];
rz(-2.4695354) q[0];
sx q[0];
rz(-1.5182537) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.034406729) q[2];
sx q[2];
rz(-0.98034562) q[2];
sx q[2];
rz(-1.7337944) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72681422) q[1];
sx q[1];
rz(-0.86514478) q[1];
sx q[1];
rz(2.3470641) q[1];
x q[2];
rz(-1.9256853) q[3];
sx q[3];
rz(-2.5657133) q[3];
sx q[3];
rz(-1.3368541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7586691) q[2];
sx q[2];
rz(-2.1047968) q[2];
sx q[2];
rz(-0.60302889) q[2];
rz(2.1488819) q[3];
sx q[3];
rz(-0.38645667) q[3];
sx q[3];
rz(0.64811289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8822766) q[0];
sx q[0];
rz(-0.74463212) q[0];
sx q[0];
rz(2.4477006) q[0];
rz(2.4248185) q[1];
sx q[1];
rz(-1.0697155) q[1];
sx q[1];
rz(-0.96376354) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1070654) q[0];
sx q[0];
rz(-1.1363159) q[0];
sx q[0];
rz(-1.700842) q[0];
x q[1];
rz(1.3974254) q[2];
sx q[2];
rz(-1.7503305) q[2];
sx q[2];
rz(2.4712359) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.239526) q[1];
sx q[1];
rz(-1.0527624) q[1];
sx q[1];
rz(-0.43398989) q[1];
rz(-pi) q[2];
rz(-0.57235897) q[3];
sx q[3];
rz(-1.5941752) q[3];
sx q[3];
rz(3.1413195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0264549) q[2];
sx q[2];
rz(-2.852735) q[2];
sx q[2];
rz(3.0333983) q[2];
rz(0.9555971) q[3];
sx q[3];
rz(-3.1028265) q[3];
sx q[3];
rz(2.5819216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3947802) q[0];
sx q[0];
rz(-2.96947) q[0];
sx q[0];
rz(-3.0346003) q[0];
rz(-0.50210285) q[1];
sx q[1];
rz(-1.5109477) q[1];
sx q[1];
rz(-2.9923901) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0923) q[0];
sx q[0];
rz(-0.23173103) q[0];
sx q[0];
rz(-0.85706237) q[0];
rz(2.4154941) q[2];
sx q[2];
rz(-2.4319785) q[2];
sx q[2];
rz(0.47898705) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2857692) q[1];
sx q[1];
rz(-2.2659784) q[1];
sx q[1];
rz(2.7421363) q[1];
rz(-2.8362422) q[3];
sx q[3];
rz(-2.6943061) q[3];
sx q[3];
rz(-1.1326912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2627829) q[2];
sx q[2];
rz(-0.97027367) q[2];
sx q[2];
rz(0.61857569) q[2];
rz(2.4570019) q[3];
sx q[3];
rz(-2.9134143) q[3];
sx q[3];
rz(-0.98606199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743643) q[0];
sx q[0];
rz(-1.3716797) q[0];
sx q[0];
rz(-0.57269639) q[0];
rz(1.0193846) q[1];
sx q[1];
rz(-0.23713325) q[1];
sx q[1];
rz(3.0651029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.50216) q[0];
sx q[0];
rz(-1.0442088) q[0];
sx q[0];
rz(-2.9798085) q[0];
rz(-pi) q[1];
rz(-0.29757737) q[2];
sx q[2];
rz(-1.7817162) q[2];
sx q[2];
rz(-0.95266137) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2092549) q[1];
sx q[1];
rz(-0.9088974) q[1];
sx q[1];
rz(-0.13929892) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42726699) q[3];
sx q[3];
rz(-0.92822853) q[3];
sx q[3];
rz(1.515732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0387592) q[2];
sx q[2];
rz(-2.2566654) q[2];
sx q[2];
rz(-2.6494675) q[2];
rz(-0.37880185) q[3];
sx q[3];
rz(-0.4327966) q[3];
sx q[3];
rz(0.88703275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73961306) q[0];
sx q[0];
rz(-2.6519863) q[0];
sx q[0];
rz(2.711645) q[0];
rz(-2.6509189) q[1];
sx q[1];
rz(-2.6538167) q[1];
sx q[1];
rz(-1.0027764) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1802964) q[0];
sx q[0];
rz(-1.8334098) q[0];
sx q[0];
rz(1.0680593) q[0];
x q[1];
rz(1.6677708) q[2];
sx q[2];
rz(-0.67296689) q[2];
sx q[2];
rz(0.99373465) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9168612) q[1];
sx q[1];
rz(-1.5229578) q[1];
sx q[1];
rz(-0.46617723) q[1];
x q[2];
rz(0.70062317) q[3];
sx q[3];
rz(-1.2081606) q[3];
sx q[3];
rz(0.36421916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3873202) q[2];
sx q[2];
rz(-0.60101271) q[2];
sx q[2];
rz(0.69166541) q[2];
rz(-3.0129041) q[3];
sx q[3];
rz(-1.5920937) q[3];
sx q[3];
rz(-0.58840978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16828951) q[0];
sx q[0];
rz(-3.0423218) q[0];
sx q[0];
rz(2.5183992) q[0];
rz(2.9469931) q[1];
sx q[1];
rz(-1.1293026) q[1];
sx q[1];
rz(2.2653939) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.122418) q[0];
sx q[0];
rz(-0.12309531) q[0];
sx q[0];
rz(1.207463) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7224465) q[2];
sx q[2];
rz(-2.306621) q[2];
sx q[2];
rz(-2.4089872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.33223486) q[1];
sx q[1];
rz(-1.7064662) q[1];
sx q[1];
rz(2.7103781) q[1];
x q[2];
rz(2.629452) q[3];
sx q[3];
rz(-2.509553) q[3];
sx q[3];
rz(1.8919945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1759922) q[2];
sx q[2];
rz(-2.8922562) q[2];
sx q[2];
rz(0.636379) q[2];
rz(1.5198358) q[3];
sx q[3];
rz(-2.2607925) q[3];
sx q[3];
rz(2.9145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31430055) q[0];
sx q[0];
rz(-1.5999595) q[0];
sx q[0];
rz(-0.40524361) q[0];
rz(1.4018519) q[1];
sx q[1];
rz(-1.4973462) q[1];
sx q[1];
rz(-3.0037465) q[1];
rz(-0.62217106) q[2];
sx q[2];
rz(-1.470574) q[2];
sx q[2];
rz(-2.1604411) q[2];
rz(0.53608175) q[3];
sx q[3];
rz(-1.6012708) q[3];
sx q[3];
rz(-0.8731245) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-1.7412269) q[0];
sx q[0];
rz(-2.9018612) q[0];
sx q[0];
rz(2.8170407) q[0];
rz(5.3311081) q[1];
sx q[1];
rz(6.0573112) q[1];
sx q[1];
rz(4.4499302) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65303923) q[0];
sx q[0];
rz(-1.4346315) q[0];
sx q[0];
rz(0.67739886) q[0];
x q[1];
rz(1.4412854) q[2];
sx q[2];
rz(-1.0900094) q[2];
sx q[2];
rz(-0.12809424) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1349012) q[1];
sx q[1];
rz(-0.54062245) q[1];
sx q[1];
rz(2.0903793) q[1];
rz(-pi) q[2];
rz(-0.69921391) q[3];
sx q[3];
rz(-1.4627856) q[3];
sx q[3];
rz(2.8349769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41210458) q[2];
sx q[2];
rz(-1.8821913) q[2];
sx q[2];
rz(-0.35120249) q[2];
rz(-1.8515733) q[3];
sx q[3];
rz(-1.131564) q[3];
sx q[3];
rz(1.9180408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7555162) q[0];
sx q[0];
rz(-1.2048683) q[0];
sx q[0];
rz(-0.93038857) q[0];
rz(-1.3049841) q[1];
sx q[1];
rz(-2.3001859) q[1];
sx q[1];
rz(-2.5321541) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064548858) q[0];
sx q[0];
rz(-1.6739558) q[0];
sx q[0];
rz(-1.1431689) q[0];
rz(-pi) q[1];
rz(1.8687042) q[2];
sx q[2];
rz(-1.4845843) q[2];
sx q[2];
rz(2.6484368) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3404559) q[1];
sx q[1];
rz(-1.1994455) q[1];
sx q[1];
rz(-2.7494829) q[1];
x q[2];
rz(1.8279161) q[3];
sx q[3];
rz(-2.5559396) q[3];
sx q[3];
rz(-0.7440834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0448407) q[2];
sx q[2];
rz(-0.98016206) q[2];
sx q[2];
rz(-0.074782221) q[2];
rz(0.52037248) q[3];
sx q[3];
rz(-2.4986391) q[3];
sx q[3];
rz(0.61409942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5388913) q[0];
sx q[0];
rz(-0.7170054) q[0];
sx q[0];
rz(-2.8394748) q[0];
rz(2.865454) q[1];
sx q[1];
rz(-1.8243676) q[1];
sx q[1];
rz(1.709323) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30423588) q[0];
sx q[0];
rz(-2.4563027) q[0];
sx q[0];
rz(0.94288148) q[0];
rz(-1.309607) q[2];
sx q[2];
rz(-2.1495594) q[2];
sx q[2];
rz(2.8566993) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36784962) q[1];
sx q[1];
rz(-2.5605695) q[1];
sx q[1];
rz(-1.8598721) q[1];
x q[2];
rz(0.071050771) q[3];
sx q[3];
rz(-1.1059575) q[3];
sx q[3];
rz(-1.1230866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8351195) q[2];
sx q[2];
rz(-0.46572954) q[2];
sx q[2];
rz(-1.2960557) q[2];
rz(2.6895788) q[3];
sx q[3];
rz(-2.0343503) q[3];
sx q[3];
rz(-2.7590416) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6119659) q[0];
sx q[0];
rz(-1.1671678) q[0];
sx q[0];
rz(-1.8091328) q[0];
rz(-2.0924856) q[1];
sx q[1];
rz(-1.0915979) q[1];
sx q[1];
rz(-1.5886935) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73686872) q[0];
sx q[0];
rz(-1.8717614) q[0];
sx q[0];
rz(0.35643453) q[0];
rz(-pi) q[1];
x q[1];
rz(0.013601549) q[2];
sx q[2];
rz(-2.0715393) q[2];
sx q[2];
rz(-2.406943) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3192352) q[1];
sx q[1];
rz(-0.93531448) q[1];
sx q[1];
rz(2.8947919) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96943198) q[3];
sx q[3];
rz(-0.80925377) q[3];
sx q[3];
rz(-1.4331897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66854746) q[2];
sx q[2];
rz(-2.1210402) q[2];
sx q[2];
rz(-2.891053) q[2];
rz(0.9969095) q[3];
sx q[3];
rz(-1.1397811) q[3];
sx q[3];
rz(0.38058773) q[3];
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
rz(pi/2) q[3];
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
rz(-0.32960358) q[0];
sx q[0];
rz(-0.488509) q[0];
sx q[0];
rz(-1.7937775) q[0];
rz(0.40428058) q[1];
sx q[1];
rz(-2.3826022) q[1];
sx q[1];
rz(1.1291198) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2166259) q[0];
sx q[0];
rz(-2.2499086) q[0];
sx q[0];
rz(-2.62519) q[0];
x q[1];
rz(-0.90520827) q[2];
sx q[2];
rz(-1.7844177) q[2];
sx q[2];
rz(-0.40494949) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89896527) q[1];
sx q[1];
rz(-2.4614442) q[1];
sx q[1];
rz(-1.4037474) q[1];
x q[2];
rz(2.8267838) q[3];
sx q[3];
rz(-1.4658576) q[3];
sx q[3];
rz(0.85280692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3471442) q[2];
sx q[2];
rz(-0.930154) q[2];
sx q[2];
rz(1.8974737) q[2];
rz(-1.0602903) q[3];
sx q[3];
rz(-2.5131707) q[3];
sx q[3];
rz(-0.30317831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65619549) q[0];
sx q[0];
rz(-0.10529101) q[0];
sx q[0];
rz(-2.7144077) q[0];
rz(3.1071013) q[1];
sx q[1];
rz(-1.5723615) q[1];
sx q[1];
rz(3.131386) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0611864) q[0];
sx q[0];
rz(-1.5690593) q[0];
sx q[0];
rz(-3.1392359) q[0];
rz(-pi) q[1];
rz(-2.9811007) q[2];
sx q[2];
rz(-2.1969323) q[2];
sx q[2];
rz(-2.6662835) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8275244) q[1];
sx q[1];
rz(-1.6351103) q[1];
sx q[1];
rz(2.2576648) q[1];
x q[2];
rz(0.56513803) q[3];
sx q[3];
rz(-2.0083154) q[3];
sx q[3];
rz(1.138162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5222142) q[2];
sx q[2];
rz(-2.7284315) q[2];
sx q[2];
rz(0.79356066) q[2];
rz(-2.2788952) q[3];
sx q[3];
rz(-1.5669275) q[3];
sx q[3];
rz(2.4376552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4392387) q[0];
sx q[0];
rz(-0.84504253) q[0];
sx q[0];
rz(-1.1335565) q[0];
rz(1.0776862) q[1];
sx q[1];
rz(-2.6262941) q[1];
sx q[1];
rz(2.8394707) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2847452) q[0];
sx q[0];
rz(-2.6805566) q[0];
sx q[0];
rz(0.7326983) q[0];
rz(-pi) q[1];
rz(-0.63150117) q[2];
sx q[2];
rz(-0.69984791) q[2];
sx q[2];
rz(1.8818784) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.97000098) q[1];
sx q[1];
rz(-2.0920894) q[1];
sx q[1];
rz(1.3212417) q[1];
rz(0.68141209) q[3];
sx q[3];
rz(-2.3294994) q[3];
sx q[3];
rz(2.8594869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2332396) q[2];
sx q[2];
rz(-2.2519604) q[2];
sx q[2];
rz(-1.3471777) q[2];
rz(-1.2498648) q[3];
sx q[3];
rz(-1.9439387) q[3];
sx q[3];
rz(-0.10716042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5477448) q[0];
sx q[0];
rz(-0.43314728) q[0];
sx q[0];
rz(-3.0331842) q[0];
rz(2.0938865) q[1];
sx q[1];
rz(-2.3240418) q[1];
sx q[1];
rz(1.0188867) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2813346) q[0];
sx q[0];
rz(-2.0906013) q[0];
sx q[0];
rz(0.61183527) q[0];
rz(-3.1347796) q[2];
sx q[2];
rz(-2.3802813) q[2];
sx q[2];
rz(-0.50466621) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8727303) q[1];
sx q[1];
rz(-0.62727562) q[1];
sx q[1];
rz(2.4226047) q[1];
rz(-pi) q[2];
rz(-2.9874727) q[3];
sx q[3];
rz(-0.74609038) q[3];
sx q[3];
rz(-0.26827677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5041647) q[2];
sx q[2];
rz(-0.97680682) q[2];
sx q[2];
rz(-2.7039779) q[2];
rz(2.6484683) q[3];
sx q[3];
rz(-0.92032856) q[3];
sx q[3];
rz(1.5639308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86579943) q[0];
sx q[0];
rz(-1.9075305) q[0];
sx q[0];
rz(2.8637874) q[0];
rz(-1.5996784) q[1];
sx q[1];
rz(-2.1733687) q[1];
sx q[1];
rz(-2.7154198) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23754643) q[0];
sx q[0];
rz(-1.5874366) q[0];
sx q[0];
rz(-0.017701935) q[0];
x q[1];
rz(-2.4075899) q[2];
sx q[2];
rz(-2.2224769) q[2];
sx q[2];
rz(-0.82254788) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8746169) q[1];
sx q[1];
rz(-0.24493453) q[1];
sx q[1];
rz(2.0744978) q[1];
x q[2];
rz(-2.3506463) q[3];
sx q[3];
rz(-0.90715796) q[3];
sx q[3];
rz(-2.5928465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.49491945) q[2];
sx q[2];
rz(-0.16725954) q[2];
sx q[2];
rz(-2.0885928) q[2];
rz(0.35887512) q[3];
sx q[3];
rz(-1.6152629) q[3];
sx q[3];
rz(1.0927965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7786355) q[0];
sx q[0];
rz(-0.21720049) q[0];
sx q[0];
rz(0.92078513) q[0];
rz(-2.6329363) q[1];
sx q[1];
rz(-0.76728907) q[1];
sx q[1];
rz(-2.7210534) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8201308) q[0];
sx q[0];
rz(-2.252914) q[0];
sx q[0];
rz(2.2338623) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2013046) q[2];
sx q[2];
rz(-1.3740525) q[2];
sx q[2];
rz(1.8041174) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1433761) q[1];
sx q[1];
rz(-0.80165473) q[1];
sx q[1];
rz(-0.90513913) q[1];
rz(-3.0402571) q[3];
sx q[3];
rz(-2.2592415) q[3];
sx q[3];
rz(2.7197414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.73021182) q[2];
sx q[2];
rz(-0.79019848) q[2];
sx q[2];
rz(-0.31663695) q[2];
rz(-2.5824879) q[3];
sx q[3];
rz(-0.63566256) q[3];
sx q[3];
rz(0.81812286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.6560787) q[0];
sx q[0];
rz(-1.0777127) q[0];
sx q[0];
rz(1.0832473) q[0];
rz(0.47646933) q[1];
sx q[1];
rz(-1.0404027) q[1];
sx q[1];
rz(-1.7146005) q[1];
rz(-3.0546679) q[2];
sx q[2];
rz(-2.3200421) q[2];
sx q[2];
rz(2.5210597) q[2];
rz(1.9540167) q[3];
sx q[3];
rz(-1.6578309) q[3];
sx q[3];
rz(2.1014392) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

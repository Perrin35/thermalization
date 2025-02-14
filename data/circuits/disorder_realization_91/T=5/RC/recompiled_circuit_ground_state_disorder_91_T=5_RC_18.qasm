OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.087091669) q[0];
sx q[0];
rz(0.85252419) q[0];
sx q[0];
rz(9.6871992) q[0];
rz(-3.4497058) q[1];
sx q[1];
rz(0.83642712) q[1];
sx q[1];
rz(15.717419) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.165929) q[0];
sx q[0];
rz(-1.0895415) q[0];
sx q[0];
rz(0.4939618) q[0];
x q[1];
rz(-2.1310102) q[2];
sx q[2];
rz(-0.081801266) q[2];
sx q[2];
rz(-2.8103925) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7345924) q[1];
sx q[1];
rz(-2.0577891) q[1];
sx q[1];
rz(-0.89501801) q[1];
x q[2];
rz(2.3289324) q[3];
sx q[3];
rz(-1.201212) q[3];
sx q[3];
rz(-0.11358914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0414163) q[2];
sx q[2];
rz(-0.88465038) q[2];
sx q[2];
rz(-2.6803988) q[2];
rz(-2.6395116) q[3];
sx q[3];
rz(-0.82379782) q[3];
sx q[3];
rz(-1.4095149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0111888) q[0];
sx q[0];
rz(-1.4606425) q[0];
sx q[0];
rz(-2.9050264) q[0];
rz(-2.323281) q[1];
sx q[1];
rz(-1.9383483) q[1];
sx q[1];
rz(-0.94258211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6494076) q[0];
sx q[0];
rz(-1.5680917) q[0];
sx q[0];
rz(1.6016141) q[0];
rz(-pi) q[1];
rz(3.0741092) q[2];
sx q[2];
rz(-1.4500256) q[2];
sx q[2];
rz(-1.6492594) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89211226) q[1];
sx q[1];
rz(-1.5693519) q[1];
sx q[1];
rz(-0.73489983) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0855186) q[3];
sx q[3];
rz(-1.266978) q[3];
sx q[3];
rz(1.3930381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.73611034) q[2];
sx q[2];
rz(-0.63850275) q[2];
sx q[2];
rz(2.45641) q[2];
rz(-0.80125609) q[3];
sx q[3];
rz(-1.5056491) q[3];
sx q[3];
rz(-1.2509468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.41246978) q[0];
sx q[0];
rz(-2.6710489) q[0];
sx q[0];
rz(2.1614918) q[0];
rz(3.0209172) q[1];
sx q[1];
rz(-1.1680892) q[1];
sx q[1];
rz(1.6623704) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9227994) q[0];
sx q[0];
rz(-0.81577071) q[0];
sx q[0];
rz(1.7936617) q[0];
rz(2.0329887) q[2];
sx q[2];
rz(-2.0304907) q[2];
sx q[2];
rz(2.5407956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.66202098) q[1];
sx q[1];
rz(-1.5618854) q[1];
sx q[1];
rz(-1.575995) q[1];
x q[2];
rz(-1.2347264) q[3];
sx q[3];
rz(-1.0072636) q[3];
sx q[3];
rz(-2.270171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4270619) q[2];
sx q[2];
rz(-1.6468628) q[2];
sx q[2];
rz(0.18043268) q[2];
rz(-0.42029542) q[3];
sx q[3];
rz(-0.97729483) q[3];
sx q[3];
rz(0.54496566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20035289) q[0];
sx q[0];
rz(-1.7914597) q[0];
sx q[0];
rz(-3.0923162) q[0];
rz(-2.409528) q[1];
sx q[1];
rz(-0.8492291) q[1];
sx q[1];
rz(-1.3759618) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7017562) q[0];
sx q[0];
rz(-1.6876093) q[0];
sx q[0];
rz(0.18420903) q[0];
rz(-pi) q[1];
rz(-1.7824729) q[2];
sx q[2];
rz(-2.3251901) q[2];
sx q[2];
rz(2.8704081) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1092618) q[1];
sx q[1];
rz(-1.9606661) q[1];
sx q[1];
rz(1.6435739) q[1];
rz(-pi) q[2];
rz(0.84305544) q[3];
sx q[3];
rz(-0.84765654) q[3];
sx q[3];
rz(1.8384733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.10776082) q[2];
sx q[2];
rz(-1.3904479) q[2];
sx q[2];
rz(0.66229171) q[2];
rz(1.3307339) q[3];
sx q[3];
rz(-2.3190506) q[3];
sx q[3];
rz(-1.2983769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0996284) q[0];
sx q[0];
rz(-0.30655107) q[0];
sx q[0];
rz(2.1339259) q[0];
rz(-0.10534605) q[1];
sx q[1];
rz(-0.65709972) q[1];
sx q[1];
rz(-2.9109921) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4661575) q[0];
sx q[0];
rz(-0.49336067) q[0];
sx q[0];
rz(0.81146474) q[0];
rz(-2.9020374) q[2];
sx q[2];
rz(-1.9223681) q[2];
sx q[2];
rz(0.4229751) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.28921738) q[1];
sx q[1];
rz(-1.7354449) q[1];
sx q[1];
rz(-1.4842121) q[1];
rz(-2.496594) q[3];
sx q[3];
rz(-1.4497744) q[3];
sx q[3];
rz(-1.5078441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0838123) q[2];
sx q[2];
rz(-2.1746706) q[2];
sx q[2];
rz(-0.18873611) q[2];
rz(1.2356637) q[3];
sx q[3];
rz(-2.6344447) q[3];
sx q[3];
rz(1.2602826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1756303) q[0];
sx q[0];
rz(-1.2164793) q[0];
sx q[0];
rz(0.29773444) q[0];
rz(-3.0689902) q[1];
sx q[1];
rz(-0.68521348) q[1];
sx q[1];
rz(-2.4593478) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52132359) q[0];
sx q[0];
rz(-0.96858874) q[0];
sx q[0];
rz(2.7839155) q[0];
rz(-pi) q[1];
rz(2.9932026) q[2];
sx q[2];
rz(-0.94252693) q[2];
sx q[2];
rz(-0.65199696) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1156716) q[1];
sx q[1];
rz(-0.37927765) q[1];
sx q[1];
rz(-1.161452) q[1];
rz(-pi) q[2];
rz(-0.97203927) q[3];
sx q[3];
rz(-1.8818731) q[3];
sx q[3];
rz(2.4233203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5785152) q[2];
sx q[2];
rz(-0.35334057) q[2];
sx q[2];
rz(-0.087470857) q[2];
rz(2.2443917) q[3];
sx q[3];
rz(-1.8950491) q[3];
sx q[3];
rz(0.38952601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0422269) q[0];
sx q[0];
rz(-2.0392188) q[0];
sx q[0];
rz(-1.9788096) q[0];
rz(-2.9258264) q[1];
sx q[1];
rz(-2.3240604) q[1];
sx q[1];
rz(-0.93528265) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1638086) q[0];
sx q[0];
rz(-2.3323614) q[0];
sx q[0];
rz(-0.33235312) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14752702) q[2];
sx q[2];
rz(-1.3152244) q[2];
sx q[2];
rz(0.35596656) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9218848) q[1];
sx q[1];
rz(-2.3517736) q[1];
sx q[1];
rz(-3.0303257) q[1];
rz(1.8963845) q[3];
sx q[3];
rz(-1.2712443) q[3];
sx q[3];
rz(1.8166145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0420578) q[2];
sx q[2];
rz(-0.65379405) q[2];
sx q[2];
rz(2.330244) q[2];
rz(0.30284303) q[3];
sx q[3];
rz(-1.1648014) q[3];
sx q[3];
rz(0.6306878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52043668) q[0];
sx q[0];
rz(-0.32670894) q[0];
sx q[0];
rz(0.69123554) q[0];
rz(2.4825545) q[1];
sx q[1];
rz(-1.6233416) q[1];
sx q[1];
rz(2.5504327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6931019) q[0];
sx q[0];
rz(-1.8536942) q[0];
sx q[0];
rz(-0.012040292) q[0];
rz(-pi) q[1];
rz(-2.065584) q[2];
sx q[2];
rz(-2.587834) q[2];
sx q[2];
rz(-2.8477235) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4860079) q[1];
sx q[1];
rz(-2.297245) q[1];
sx q[1];
rz(-2.6996711) q[1];
rz(-pi) q[2];
rz(-3.0847315) q[3];
sx q[3];
rz(-1.8627852) q[3];
sx q[3];
rz(-0.1869456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.14472321) q[2];
sx q[2];
rz(-0.99088061) q[2];
sx q[2];
rz(1.4608176) q[2];
rz(-2.3878494) q[3];
sx q[3];
rz(-1.4494579) q[3];
sx q[3];
rz(-0.47372216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4562456) q[0];
sx q[0];
rz(-1.101475) q[0];
sx q[0];
rz(2.9154678) q[0];
rz(-1.6926758) q[1];
sx q[1];
rz(-0.96335226) q[1];
sx q[1];
rz(2.1400145) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9739363) q[0];
sx q[0];
rz(-1.4261803) q[0];
sx q[0];
rz(2.5517716) q[0];
rz(-pi) q[1];
rz(1.9492769) q[2];
sx q[2];
rz(-0.90736249) q[2];
sx q[2];
rz(1.1089693) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4521617) q[1];
sx q[1];
rz(-2.9254483) q[1];
sx q[1];
rz(0.5139022) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36205451) q[3];
sx q[3];
rz(-0.71740642) q[3];
sx q[3];
rz(2.2256454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.61233026) q[2];
sx q[2];
rz(-0.77028209) q[2];
sx q[2];
rz(-0.15753499) q[2];
rz(-2.9869288) q[3];
sx q[3];
rz(-1.1636584) q[3];
sx q[3];
rz(2.7111588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4100819) q[0];
sx q[0];
rz(-1.2393351) q[0];
sx q[0];
rz(-2.4391158) q[0];
rz(0.53600535) q[1];
sx q[1];
rz(-0.82967007) q[1];
sx q[1];
rz(-1.3943256) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6781021) q[0];
sx q[0];
rz(-2.9759088) q[0];
sx q[0];
rz(-1.413762) q[0];
x q[1];
rz(-0.86081204) q[2];
sx q[2];
rz(-2.99798) q[2];
sx q[2];
rz(1.8710305) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6977384) q[1];
sx q[1];
rz(-1.2875673) q[1];
sx q[1];
rz(-2.1069788) q[1];
rz(-pi) q[2];
rz(-1.570618) q[3];
sx q[3];
rz(-1.9863673) q[3];
sx q[3];
rz(-3.0316169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.006762) q[2];
sx q[2];
rz(-0.68209058) q[2];
sx q[2];
rz(1.4466059) q[2];
rz(0.26099482) q[3];
sx q[3];
rz(-2.1311396) q[3];
sx q[3];
rz(1.9058156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9217054) q[0];
sx q[0];
rz(-2.4867262) q[0];
sx q[0];
rz(2.126271) q[0];
rz(-1.075853) q[1];
sx q[1];
rz(-1.2747819) q[1];
sx q[1];
rz(-1.4631974) q[1];
rz(-1.4599316) q[2];
sx q[2];
rz(-1.9456429) q[2];
sx q[2];
rz(2.6425998) q[2];
rz(0.45810926) q[3];
sx q[3];
rz(-0.81732133) q[3];
sx q[3];
rz(2.2504239) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

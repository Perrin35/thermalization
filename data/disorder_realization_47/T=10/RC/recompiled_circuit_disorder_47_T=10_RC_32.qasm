OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3818504) q[0];
sx q[0];
rz(-0.83431017) q[0];
sx q[0];
rz(-0.068339737) q[0];
rz(2.4812658) q[1];
sx q[1];
rz(-2.2934409) q[1];
sx q[1];
rz(3.1037722) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8605182) q[0];
sx q[0];
rz(-1.7639419) q[0];
sx q[0];
rz(0.079695745) q[0];
x q[1];
rz(-2.2322466) q[2];
sx q[2];
rz(-1.4771909) q[2];
sx q[2];
rz(-2.0441165) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0953857) q[1];
sx q[1];
rz(-2.3933105) q[1];
sx q[1];
rz(-0.074949646) q[1];
rz(-pi) q[2];
rz(-0.1511729) q[3];
sx q[3];
rz(-1.2281979) q[3];
sx q[3];
rz(1.3959988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7906856) q[2];
sx q[2];
rz(-2.6280792) q[2];
sx q[2];
rz(1.8581871) q[2];
rz(-0.12617271) q[3];
sx q[3];
rz(-1.4161371) q[3];
sx q[3];
rz(-0.090601966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44777563) q[0];
sx q[0];
rz(-0.87647804) q[0];
sx q[0];
rz(-2.5449261) q[0];
rz(1.5860575) q[1];
sx q[1];
rz(-1.3417473) q[1];
sx q[1];
rz(-1.3635051) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039149337) q[0];
sx q[0];
rz(-1.1457232) q[0];
sx q[0];
rz(0.58702472) q[0];
rz(-pi) q[1];
rz(-3.0599669) q[2];
sx q[2];
rz(-2.584169) q[2];
sx q[2];
rz(-2.4068085) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.23681444) q[1];
sx q[1];
rz(-1.292255) q[1];
sx q[1];
rz(-0.65794557) q[1];
x q[2];
rz(-0.73804654) q[3];
sx q[3];
rz(-2.1481206) q[3];
sx q[3];
rz(2.9374124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3229225) q[2];
sx q[2];
rz(-1.0007891) q[2];
sx q[2];
rz(0.73454109) q[2];
rz(0.62720403) q[3];
sx q[3];
rz(-1.5137129) q[3];
sx q[3];
rz(-0.74497765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4194141) q[0];
sx q[0];
rz(-0.83291554) q[0];
sx q[0];
rz(0.27221361) q[0];
rz(2.294337) q[1];
sx q[1];
rz(-1.3191185) q[1];
sx q[1];
rz(-2.8289657) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74041849) q[0];
sx q[0];
rz(-2.9156988) q[0];
sx q[0];
rz(-2.9001539) q[0];
rz(-2.0414511) q[2];
sx q[2];
rz(-2.0690284) q[2];
sx q[2];
rz(2.2217896) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9885013) q[1];
sx q[1];
rz(-1.3107436) q[1];
sx q[1];
rz(-3.0696763) q[1];
x q[2];
rz(2.3787141) q[3];
sx q[3];
rz(-1.2925914) q[3];
sx q[3];
rz(2.0335846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.12038885) q[2];
sx q[2];
rz(-2.618232) q[2];
sx q[2];
rz(0.81494251) q[2];
rz(-0.96873823) q[3];
sx q[3];
rz(-1.0453984) q[3];
sx q[3];
rz(-0.43280861) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75354904) q[0];
sx q[0];
rz(-2.9409445) q[0];
sx q[0];
rz(2.5182305) q[0];
rz(-2.3240044) q[1];
sx q[1];
rz(-2.8249884) q[1];
sx q[1];
rz(-3.0923016) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3885956) q[0];
sx q[0];
rz(-1.5224644) q[0];
sx q[0];
rz(-0.81141332) q[0];
x q[1];
rz(-1.9360147) q[2];
sx q[2];
rz(-0.74539241) q[2];
sx q[2];
rz(-0.12860194) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.99187461) q[1];
sx q[1];
rz(-1.3279993) q[1];
sx q[1];
rz(-1.2319831) q[1];
rz(-pi) q[2];
rz(1.0443346) q[3];
sx q[3];
rz(-1.1757848) q[3];
sx q[3];
rz(-2.928424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7307044) q[2];
sx q[2];
rz(-1.5894019) q[2];
sx q[2];
rz(0.98497406) q[2];
rz(-2.2385521) q[3];
sx q[3];
rz(-1.1629546) q[3];
sx q[3];
rz(-2.8682017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7966998) q[0];
sx q[0];
rz(-1.5364237) q[0];
sx q[0];
rz(-3.0539736) q[0];
rz(2.9852729) q[1];
sx q[1];
rz(-2.5904398) q[1];
sx q[1];
rz(-0.87096754) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317469) q[0];
sx q[0];
rz(-1.5994659) q[0];
sx q[0];
rz(0.28733758) q[0];
rz(-pi) q[1];
rz(1.1768612) q[2];
sx q[2];
rz(-1.7125868) q[2];
sx q[2];
rz(-1.0172435) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.761913) q[1];
sx q[1];
rz(-1.8508136) q[1];
sx q[1];
rz(0.67358394) q[1];
rz(-pi) q[2];
rz(-0.80645873) q[3];
sx q[3];
rz(-1.6975003) q[3];
sx q[3];
rz(2.5106631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0481723) q[2];
sx q[2];
rz(-2.293736) q[2];
sx q[2];
rz(-2.7887153) q[2];
rz(-2.2757163) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(-0.29233366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6784994) q[0];
sx q[0];
rz(-1.1266288) q[0];
sx q[0];
rz(3.034806) q[0];
rz(1.1865901) q[1];
sx q[1];
rz(-1.0266961) q[1];
sx q[1];
rz(-2.1170763) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5750676) q[0];
sx q[0];
rz(-1.1199513) q[0];
sx q[0];
rz(-1.2498115) q[0];
rz(-pi) q[1];
rz(0.60501955) q[2];
sx q[2];
rz(-1.1875249) q[2];
sx q[2];
rz(0.59125102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4752794) q[1];
sx q[1];
rz(-0.76994714) q[1];
sx q[1];
rz(2.2837) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2023737) q[3];
sx q[3];
rz(-1.107186) q[3];
sx q[3];
rz(-1.1045477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7496877) q[2];
sx q[2];
rz(-1.2140032) q[2];
sx q[2];
rz(0.8708896) q[2];
rz(1.8093367) q[3];
sx q[3];
rz(-2.199316) q[3];
sx q[3];
rz(-1.4107305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9695327) q[0];
sx q[0];
rz(-1.7344069) q[0];
sx q[0];
rz(2.5906738) q[0];
rz(-0.46539601) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(2.904772) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2315002) q[0];
sx q[0];
rz(-1.6084558) q[0];
sx q[0];
rz(-1.4715657) q[0];
rz(2.9914855) q[2];
sx q[2];
rz(-1.7068212) q[2];
sx q[2];
rz(2.0748236) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.19141087) q[1];
sx q[1];
rz(-1.6457874) q[1];
sx q[1];
rz(-1.8412776) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83644609) q[3];
sx q[3];
rz(-1.340938) q[3];
sx q[3];
rz(-3.0695855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.44234309) q[2];
sx q[2];
rz(-2.4599059) q[2];
sx q[2];
rz(2.701475) q[2];
rz(1.0951428) q[3];
sx q[3];
rz(-1.4705855) q[3];
sx q[3];
rz(1.7949036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8286164) q[0];
sx q[0];
rz(-1.3422817) q[0];
sx q[0];
rz(3.112088) q[0];
rz(2.1942031) q[1];
sx q[1];
rz(-2.3209929) q[1];
sx q[1];
rz(-0.11880076) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9751562) q[0];
sx q[0];
rz(-1.6140037) q[0];
sx q[0];
rz(-0.31063609) q[0];
rz(-pi) q[1];
x q[1];
rz(2.253184) q[2];
sx q[2];
rz(-2.4783122) q[2];
sx q[2];
rz(2.1036489) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6606632) q[1];
sx q[1];
rz(-2.5813563) q[1];
sx q[1];
rz(0.57628298) q[1];
rz(-1.8700065) q[3];
sx q[3];
rz(-2.6282103) q[3];
sx q[3];
rz(1.3083003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7416731) q[2];
sx q[2];
rz(-0.46367773) q[2];
sx q[2];
rz(1.5734394) q[2];
rz(1.167477) q[3];
sx q[3];
rz(-1.4195331) q[3];
sx q[3];
rz(-1.8673816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90299273) q[0];
sx q[0];
rz(-2.3165343) q[0];
sx q[0];
rz(-0.15326823) q[0];
rz(1.0614456) q[1];
sx q[1];
rz(-2.5669284) q[1];
sx q[1];
rz(-1.1538039) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.823258) q[0];
sx q[0];
rz(-1.9360006) q[0];
sx q[0];
rz(-1.2814786) q[0];
rz(-pi) q[1];
rz(-0.54406464) q[2];
sx q[2];
rz(-2.1779034) q[2];
sx q[2];
rz(2.0008759) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8037387) q[1];
sx q[1];
rz(-1.5152463) q[1];
sx q[1];
rz(-1.7800203) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88345248) q[3];
sx q[3];
rz(-1.6646107) q[3];
sx q[3];
rz(-1.0132307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8447421) q[2];
sx q[2];
rz(-1.1721609) q[2];
sx q[2];
rz(1.5273757) q[2];
rz(-0.99689233) q[3];
sx q[3];
rz(-1.6485873) q[3];
sx q[3];
rz(0.87866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8940354) q[0];
sx q[0];
rz(-1.0972801) q[0];
sx q[0];
rz(-0.19009185) q[0];
rz(-0.67063531) q[1];
sx q[1];
rz(-1.1963444) q[1];
sx q[1];
rz(1.6533096) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0310693) q[0];
sx q[0];
rz(-3.0123683) q[0];
sx q[0];
rz(-2.9002225) q[0];
x q[1];
rz(0.44234862) q[2];
sx q[2];
rz(-1.3253691) q[2];
sx q[2];
rz(2.6128212) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.30257672) q[1];
sx q[1];
rz(-0.95514983) q[1];
sx q[1];
rz(2.2018196) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90115746) q[3];
sx q[3];
rz(-0.81448758) q[3];
sx q[3];
rz(2.8709473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1378479) q[2];
sx q[2];
rz(-0.90128428) q[2];
sx q[2];
rz(-2.2424662) q[2];
rz(-0.51504618) q[3];
sx q[3];
rz(-2.6656272) q[3];
sx q[3];
rz(2.9639444) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0777733) q[0];
sx q[0];
rz(-1.780092) q[0];
sx q[0];
rz(-0.5583981) q[0];
rz(0.28941119) q[1];
sx q[1];
rz(-0.90129539) q[1];
sx q[1];
rz(1.7064066) q[1];
rz(-1.2895073) q[2];
sx q[2];
rz(-2.2608902) q[2];
sx q[2];
rz(2.8534941) q[2];
rz(-0.12712052) q[3];
sx q[3];
rz(-1.8416497) q[3];
sx q[3];
rz(1.9491378) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

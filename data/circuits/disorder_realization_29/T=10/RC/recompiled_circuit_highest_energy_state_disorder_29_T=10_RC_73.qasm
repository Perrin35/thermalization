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
rz(0.88275498) q[0];
sx q[0];
rz(3.2864154) q[0];
sx q[0];
rz(10.176131) q[0];
rz(1.524628) q[1];
sx q[1];
rz(-2.1722062) q[1];
sx q[1];
rz(0.17925395) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67231945) q[0];
sx q[0];
rz(-2.0600962) q[0];
sx q[0];
rz(2.8239522) q[0];
rz(-pi) q[1];
rz(-0.037161552) q[2];
sx q[2];
rz(-2.4183309) q[2];
sx q[2];
rz(3.0918372) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4262095) q[1];
sx q[1];
rz(-2.8236964) q[1];
sx q[1];
rz(-0.23178394) q[1];
rz(-0.93710812) q[3];
sx q[3];
rz(-1.5230012) q[3];
sx q[3];
rz(2.5339047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90193191) q[2];
sx q[2];
rz(-1.8472981) q[2];
sx q[2];
rz(-2.530976) q[2];
rz(0.88879746) q[3];
sx q[3];
rz(-0.68032467) q[3];
sx q[3];
rz(-2.4077267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69773847) q[0];
sx q[0];
rz(-2.839851) q[0];
sx q[0];
rz(1.3296211) q[0];
rz(-1.656172) q[1];
sx q[1];
rz(-1.4621567) q[1];
sx q[1];
rz(0.86404538) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5764424) q[0];
sx q[0];
rz(-0.2731384) q[0];
sx q[0];
rz(-1.9594203) q[0];
rz(2.8301622) q[2];
sx q[2];
rz(-1.127178) q[2];
sx q[2];
rz(-1.6964427) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.73693528) q[1];
sx q[1];
rz(-0.22316775) q[1];
sx q[1];
rz(-2.2523802) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1883499) q[3];
sx q[3];
rz(-0.88897486) q[3];
sx q[3];
rz(0.98495959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0580505) q[2];
sx q[2];
rz(-2.4594049) q[2];
sx q[2];
rz(-0.19134276) q[2];
rz(-0.30609104) q[3];
sx q[3];
rz(-1.6813262) q[3];
sx q[3];
rz(-2.2050841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3092344) q[0];
sx q[0];
rz(-1.2795871) q[0];
sx q[0];
rz(2.7171296) q[0];
rz(-1.8008495) q[1];
sx q[1];
rz(-2.2258591) q[1];
sx q[1];
rz(0.65139604) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0754335) q[0];
sx q[0];
rz(-1.4071583) q[0];
sx q[0];
rz(0.0072102116) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7119262) q[2];
sx q[2];
rz(-0.41791195) q[2];
sx q[2];
rz(-1.9442476) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.60247682) q[1];
sx q[1];
rz(-2.0664363) q[1];
sx q[1];
rz(-1.3631166) q[1];
rz(-pi) q[2];
rz(-0.96730729) q[3];
sx q[3];
rz(-2.8684232) q[3];
sx q[3];
rz(0.73707132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.035630781) q[2];
sx q[2];
rz(-1.1392081) q[2];
sx q[2];
rz(2.4448709) q[2];
rz(-2.4524955) q[3];
sx q[3];
rz(-1.1545811) q[3];
sx q[3];
rz(2.885163) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43831393) q[0];
sx q[0];
rz(-0.2916446) q[0];
sx q[0];
rz(1.0003723) q[0];
rz(1.4601624) q[1];
sx q[1];
rz(-1.6078452) q[1];
sx q[1];
rz(1.9283074) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54736894) q[0];
sx q[0];
rz(-1.4495736) q[0];
sx q[0];
rz(1.8531043) q[0];
rz(-pi) q[1];
rz(2.6268509) q[2];
sx q[2];
rz(-2.3433583) q[2];
sx q[2];
rz(-1.2003743) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.91649969) q[1];
sx q[1];
rz(-1.2429534) q[1];
sx q[1];
rz(1.6986548) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4198076) q[3];
sx q[3];
rz(-0.54364341) q[3];
sx q[3];
rz(-1.0796987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1022243) q[2];
sx q[2];
rz(-0.82304707) q[2];
sx q[2];
rz(-0.064112045) q[2];
rz(0.72470775) q[3];
sx q[3];
rz(-1.9331845) q[3];
sx q[3];
rz(0.39976111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0475912) q[0];
sx q[0];
rz(-1.2971017) q[0];
sx q[0];
rz(-0.79175788) q[0];
rz(-2.0296312) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(-1.8066822) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34599538) q[0];
sx q[0];
rz(-0.79774081) q[0];
sx q[0];
rz(2.255385) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15891517) q[2];
sx q[2];
rz(-0.58664413) q[2];
sx q[2];
rz(-1.212709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6688456) q[1];
sx q[1];
rz(-0.18660523) q[1];
sx q[1];
rz(-1.1531468) q[1];
rz(2.9369257) q[3];
sx q[3];
rz(-1.0042666) q[3];
sx q[3];
rz(3.1086904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2935334) q[2];
sx q[2];
rz(-0.89911014) q[2];
sx q[2];
rz(-2.000957) q[2];
rz(0.10661495) q[3];
sx q[3];
rz(-1.5299608) q[3];
sx q[3];
rz(-2.831736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3265729) q[0];
sx q[0];
rz(-2.6928379) q[0];
sx q[0];
rz(0.003224592) q[0];
rz(3.0942753) q[1];
sx q[1];
rz(-1.4553757) q[1];
sx q[1];
rz(-1.3501732) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91302204) q[0];
sx q[0];
rz(-1.6294894) q[0];
sx q[0];
rz(-1.8236158) q[0];
rz(-pi) q[1];
rz(-0.90099868) q[2];
sx q[2];
rz(-2.0428847) q[2];
sx q[2];
rz(0.53874082) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1803169) q[1];
sx q[1];
rz(-2.2928228) q[1];
sx q[1];
rz(0.78603334) q[1];
rz(1.7058701) q[3];
sx q[3];
rz(-0.22264847) q[3];
sx q[3];
rz(-2.6378353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1515767) q[2];
sx q[2];
rz(-0.185597) q[2];
sx q[2];
rz(-3.0604176) q[2];
rz(-2.2422527) q[3];
sx q[3];
rz(-2.3266561) q[3];
sx q[3];
rz(1.8544633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7154295) q[0];
sx q[0];
rz(-1.8303215) q[0];
sx q[0];
rz(-2.6370866) q[0];
rz(-2.1465178) q[1];
sx q[1];
rz(-1.0202531) q[1];
sx q[1];
rz(3.1212433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5581455) q[0];
sx q[0];
rz(-0.37675315) q[0];
sx q[0];
rz(-1.8280297) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3624914) q[2];
sx q[2];
rz(-1.4464738) q[2];
sx q[2];
rz(-0.97036874) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87096244) q[1];
sx q[1];
rz(-0.32901627) q[1];
sx q[1];
rz(0.94493072) q[1];
rz(-1.4844839) q[3];
sx q[3];
rz(-1.619297) q[3];
sx q[3];
rz(-2.8750471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4703579) q[2];
sx q[2];
rz(-1.8378259) q[2];
sx q[2];
rz(-1.7669558) q[2];
rz(1.5291519) q[3];
sx q[3];
rz(-1.0498472) q[3];
sx q[3];
rz(-3.0488739) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2328211) q[0];
sx q[0];
rz(-0.11255539) q[0];
sx q[0];
rz(1.0821279) q[0];
rz(-0.96083653) q[1];
sx q[1];
rz(-1.4832393) q[1];
sx q[1];
rz(0.94295162) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43555194) q[0];
sx q[0];
rz(-2.0461975) q[0];
sx q[0];
rz(-0.33669223) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8846747) q[2];
sx q[2];
rz(-2.4853443) q[2];
sx q[2];
rz(-0.73168025) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8170191) q[1];
sx q[1];
rz(-0.96265154) q[1];
sx q[1];
rz(-2.2961246) q[1];
rz(-3.037463) q[3];
sx q[3];
rz(-1.8767329) q[3];
sx q[3];
rz(0.96398523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1711787) q[2];
sx q[2];
rz(-1.3033988) q[2];
sx q[2];
rz(1.9937493) q[2];
rz(-0.36192274) q[3];
sx q[3];
rz(-1.7636834) q[3];
sx q[3];
rz(-1.743478) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4107133) q[0];
sx q[0];
rz(-0.35922265) q[0];
sx q[0];
rz(-0.10243375) q[0];
rz(2.5275285) q[1];
sx q[1];
rz(-2.1153317) q[1];
sx q[1];
rz(2.3238497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21120223) q[0];
sx q[0];
rz(-1.0387392) q[0];
sx q[0];
rz(2.3389057) q[0];
x q[1];
rz(-3.0258437) q[2];
sx q[2];
rz(-0.90311804) q[2];
sx q[2];
rz(-2.4493461) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9407258) q[1];
sx q[1];
rz(-1.5518909) q[1];
sx q[1];
rz(-2.3000642) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6800167) q[3];
sx q[3];
rz(-2.1510923) q[3];
sx q[3];
rz(0.020315276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7927336) q[2];
sx q[2];
rz(-1.0772971) q[2];
sx q[2];
rz(-0.31361541) q[2];
rz(2.6775728) q[3];
sx q[3];
rz(-2.2950164) q[3];
sx q[3];
rz(2.2217506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87026507) q[0];
sx q[0];
rz(-0.79065228) q[0];
sx q[0];
rz(2.9050997) q[0];
rz(-0.21367167) q[1];
sx q[1];
rz(-0.90286076) q[1];
sx q[1];
rz(0.22458354) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29447039) q[0];
sx q[0];
rz(-2.0531539) q[0];
sx q[0];
rz(2.3806974) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.557569) q[2];
sx q[2];
rz(-1.8796325) q[2];
sx q[2];
rz(-0.14039224) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.869369) q[1];
sx q[1];
rz(-1.0450796) q[1];
sx q[1];
rz(1.438526) q[1];
x q[2];
rz(2.2891232) q[3];
sx q[3];
rz(-1.895088) q[3];
sx q[3];
rz(1.7773903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.97682041) q[2];
sx q[2];
rz(-0.87146622) q[2];
sx q[2];
rz(0.15178794) q[2];
rz(0.031489059) q[3];
sx q[3];
rz(-1.7446691) q[3];
sx q[3];
rz(0.12846863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0781773) q[0];
sx q[0];
rz(-1.5237533) q[0];
sx q[0];
rz(-0.40405003) q[0];
rz(-2.0582485) q[1];
sx q[1];
rz(-1.5768408) q[1];
sx q[1];
rz(1.5595938) q[1];
rz(-2.940098) q[2];
sx q[2];
rz(-1.6215743) q[2];
sx q[2];
rz(-1.8334186) q[2];
rz(-1.6733698) q[3];
sx q[3];
rz(-1.8986618) q[3];
sx q[3];
rz(-0.31019216) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

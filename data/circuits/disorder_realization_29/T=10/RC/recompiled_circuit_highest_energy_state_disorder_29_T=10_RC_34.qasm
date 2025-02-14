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
rz(-1.6169647) q[1];
sx q[1];
rz(-0.96938649) q[1];
sx q[1];
rz(-0.17925395) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0796666) q[0];
sx q[0];
rz(-0.57625895) q[0];
sx q[0];
rz(2.1012596) q[0];
rz(-pi) q[1];
rz(1.6035853) q[2];
sx q[2];
rz(-0.84814397) q[2];
sx q[2];
rz(0.099309534) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1826659) q[1];
sx q[1];
rz(-1.2616871) q[1];
sx q[1];
rz(-1.6462417) q[1];
x q[2];
rz(-3.0823067) q[3];
sx q[3];
rz(-2.2036457) q[3];
sx q[3];
rz(-2.2135753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.90193191) q[2];
sx q[2];
rz(-1.8472981) q[2];
sx q[2];
rz(0.61061668) q[2];
rz(0.88879746) q[3];
sx q[3];
rz(-0.68032467) q[3];
sx q[3];
rz(-2.4077267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69773847) q[0];
sx q[0];
rz(-2.839851) q[0];
sx q[0];
rz(1.8119716) q[0];
rz(1.4854206) q[1];
sx q[1];
rz(-1.4621567) q[1];
sx q[1];
rz(0.86404538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7603908) q[0];
sx q[0];
rz(-1.6731894) q[0];
sx q[0];
rz(-1.824462) q[0];
rz(2.8301622) q[2];
sx q[2];
rz(-1.127178) q[2];
sx q[2];
rz(1.4451499) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4046574) q[1];
sx q[1];
rz(-2.9184249) q[1];
sx q[1];
rz(2.2523802) q[1];
rz(-pi) q[2];
rz(-0.78317376) q[3];
sx q[3];
rz(-1.1044972) q[3];
sx q[3];
rz(0.1649905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0580505) q[2];
sx q[2];
rz(-0.6821878) q[2];
sx q[2];
rz(0.19134276) q[2];
rz(0.30609104) q[3];
sx q[3];
rz(-1.6813262) q[3];
sx q[3];
rz(2.2050841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.8323583) q[0];
sx q[0];
rz(-1.8620055) q[0];
sx q[0];
rz(0.424463) q[0];
rz(-1.8008495) q[1];
sx q[1];
rz(-0.91573358) q[1];
sx q[1];
rz(2.4901966) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1196647) q[0];
sx q[0];
rz(-0.16379539) q[0];
sx q[0];
rz(-1.6144362) q[0];
rz(-pi) q[1];
rz(1.3878787) q[2];
sx q[2];
rz(-1.1929034) q[2];
sx q[2];
rz(1.479508) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9561852) q[1];
sx q[1];
rz(-0.53403097) q[1];
sx q[1];
rz(2.7772481) q[1];
x q[2];
rz(0.15768361) q[3];
sx q[3];
rz(-1.3467977) q[3];
sx q[3];
rz(-1.7832613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1059619) q[2];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7032787) q[0];
sx q[0];
rz(-0.2916446) q[0];
sx q[0];
rz(-1.0003723) q[0];
rz(-1.6814303) q[1];
sx q[1];
rz(-1.6078452) q[1];
sx q[1];
rz(-1.2132852) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1532261) q[0];
sx q[0];
rz(-1.2906162) q[0];
sx q[0];
rz(0.12616726) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0385267) q[2];
sx q[2];
rz(-2.2437895) q[2];
sx q[2];
rz(-1.2601992) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6049543) q[1];
sx q[1];
rz(-0.35105303) q[1];
sx q[1];
rz(2.7829079) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0322047) q[3];
sx q[3];
rz(-1.6486787) q[3];
sx q[3];
rz(0.36161446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0393684) q[2];
sx q[2];
rz(-2.3185456) q[2];
sx q[2];
rz(3.0774806) q[2];
rz(-2.4168849) q[3];
sx q[3];
rz(-1.9331845) q[3];
sx q[3];
rz(-2.7418315) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0940014) q[0];
sx q[0];
rz(-1.8444909) q[0];
sx q[0];
rz(-0.79175788) q[0];
rz(-1.1119615) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(1.8066822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2091141) q[0];
sx q[0];
rz(-2.1585585) q[0];
sx q[0];
rz(2.5665119) q[0];
x q[1];
rz(2.9826775) q[2];
sx q[2];
rz(-0.58664413) q[2];
sx q[2];
rz(1.9288837) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.47274703) q[1];
sx q[1];
rz(-2.9549874) q[1];
sx q[1];
rz(-1.1531468) q[1];
rz(-pi) q[2];
rz(1.2615292) q[3];
sx q[3];
rz(-0.59854186) q[3];
sx q[3];
rz(2.8054597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2935334) q[2];
sx q[2];
rz(-0.89911014) q[2];
sx q[2];
rz(-1.1406356) q[2];
rz(3.0349777) q[3];
sx q[3];
rz(-1.6116319) q[3];
sx q[3];
rz(0.30985668) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3265729) q[0];
sx q[0];
rz(-0.44875479) q[0];
sx q[0];
rz(-3.1383681) q[0];
rz(-3.0942753) q[1];
sx q[1];
rz(-1.686217) q[1];
sx q[1];
rz(-1.3501732) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91302204) q[0];
sx q[0];
rz(-1.6294894) q[0];
sx q[0];
rz(1.3179768) q[0];
rz(-2.240594) q[2];
sx q[2];
rz(-2.0428847) q[2];
sx q[2];
rz(-0.53874082) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9474831) q[1];
sx q[1];
rz(-1.0110185) q[1];
sx q[1];
rz(0.67621381) q[1];
x q[2];
rz(1.7058701) q[3];
sx q[3];
rz(-0.22264847) q[3];
sx q[3];
rz(-2.6378353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1515767) q[2];
sx q[2];
rz(-2.9559957) q[2];
sx q[2];
rz(3.0604176) q[2];
rz(0.89933991) q[3];
sx q[3];
rz(-2.3266561) q[3];
sx q[3];
rz(-1.2871294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4261632) q[0];
sx q[0];
rz(-1.3112712) q[0];
sx q[0];
rz(-0.50450605) q[0];
rz(2.1465178) q[1];
sx q[1];
rz(-2.1213396) q[1];
sx q[1];
rz(-0.02034932) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.394144) q[0];
sx q[0];
rz(-1.6645303) q[0];
sx q[0];
rz(-1.9362437) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7791012) q[2];
sx q[2];
rz(-1.4464738) q[2];
sx q[2];
rz(-0.97036874) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.87096244) q[1];
sx q[1];
rz(-2.8125764) q[1];
sx q[1];
rz(-0.94493072) q[1];
x q[2];
rz(-3.0929111) q[3];
sx q[3];
rz(-1.657007) q[3];
sx q[3];
rz(1.3000559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4703579) q[2];
sx q[2];
rz(-1.3037668) q[2];
sx q[2];
rz(1.7669558) q[2];
rz(1.6124407) q[3];
sx q[3];
rz(-1.0498472) q[3];
sx q[3];
rz(3.0488739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2328211) q[0];
sx q[0];
rz(-0.11255539) q[0];
sx q[0];
rz(-1.0821279) q[0];
rz(0.96083653) q[1];
sx q[1];
rz(-1.6583534) q[1];
sx q[1];
rz(-2.198641) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0531702) q[0];
sx q[0];
rz(-2.5665847) q[0];
sx q[0];
rz(1.0002329) q[0];
x q[1];
rz(0.23343691) q[2];
sx q[2];
rz(-0.95165157) q[2];
sx q[2];
rz(-0.34282986) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.3245736) q[1];
sx q[1];
rz(-2.1789411) q[1];
sx q[1];
rz(-0.84546802) q[1];
rz(-pi) q[2];
rz(0.10412962) q[3];
sx q[3];
rz(-1.2648598) q[3];
sx q[3];
rz(-0.96398523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1711787) q[2];
sx q[2];
rz(-1.3033988) q[2];
sx q[2];
rz(-1.9937493) q[2];
rz(-2.7796699) q[3];
sx q[3];
rz(-1.3779093) q[3];
sx q[3];
rz(1.3981147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4107133) q[0];
sx q[0];
rz(-0.35922265) q[0];
sx q[0];
rz(3.0391589) q[0];
rz(0.61406413) q[1];
sx q[1];
rz(-2.1153317) q[1];
sx q[1];
rz(0.817743) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8431658) q[0];
sx q[0];
rz(-0.90230391) q[0];
sx q[0];
rz(-2.2737204) q[0];
rz(-pi) q[1];
rz(-1.7162343) q[2];
sx q[2];
rz(-0.67611968) q[2];
sx q[2];
rz(0.506625) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.35303822) q[1];
sx q[1];
rz(-0.8416881) q[1];
sx q[1];
rz(-3.1162412) q[1];
rz(-pi) q[2];
rz(-2.1675112) q[3];
sx q[3];
rz(-2.4170205) q[3];
sx q[3];
rz(0.71648471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3488591) q[2];
sx q[2];
rz(-2.0642955) q[2];
sx q[2];
rz(-2.8279772) q[2];
rz(0.46401986) q[3];
sx q[3];
rz(-2.2950164) q[3];
sx q[3];
rz(0.919842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2713276) q[0];
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
rz(2.8471223) q[0];
sx q[0];
rz(-1.0884388) q[0];
sx q[0];
rz(0.76089528) q[0];
x q[1];
rz(0.5245536) q[2];
sx q[2];
rz(-2.4894425) q[2];
sx q[2];
rz(0.99936501) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3652412) q[1];
sx q[1];
rz(-1.4564716) q[1];
sx q[1];
rz(-2.6120595) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8524695) q[3];
sx q[3];
rz(-1.895088) q[3];
sx q[3];
rz(1.7773903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1647722) q[2];
sx q[2];
rz(-0.87146622) q[2];
sx q[2];
rz(0.15178794) q[2];
rz(-3.1101036) q[3];
sx q[3];
rz(-1.7446691) q[3];
sx q[3];
rz(-3.013124) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0634154) q[0];
sx q[0];
rz(-1.5237533) q[0];
sx q[0];
rz(-0.40405003) q[0];
rz(1.0833441) q[1];
sx q[1];
rz(-1.5768408) q[1];
sx q[1];
rz(1.5595938) q[1];
rz(-1.5189717) q[2];
sx q[2];
rz(-1.3695649) q[2];
sx q[2];
rz(2.8893378) q[2];
rz(2.8121171) q[3];
sx q[3];
rz(-1.4737045) q[3];
sx q[3];
rz(1.2274689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.4438348) q[0];
sx q[0];
rz(4.4153089) q[0];
sx q[0];
rz(8.4755228) q[0];
rz(-1.0957837) q[1];
sx q[1];
rz(2.1646808) q[1];
sx q[1];
rz(10.77471) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5744517) q[0];
sx q[0];
rz(-1.0691088) q[0];
sx q[0];
rz(-1.6654963) q[0];
x q[1];
rz(-0.33265661) q[2];
sx q[2];
rz(-1.4321308) q[2];
sx q[2];
rz(2.0227416) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3826582) q[1];
sx q[1];
rz(-1.9449502) q[1];
sx q[1];
rz(1.6439245) q[1];
rz(-2.8676991) q[3];
sx q[3];
rz(-1.9889758) q[3];
sx q[3];
rz(-2.1194715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.533941) q[2];
sx q[2];
rz(-1.093163) q[2];
sx q[2];
rz(2.9489813) q[2];
rz(1.2827778) q[3];
sx q[3];
rz(-2.0359813) q[3];
sx q[3];
rz(-0.58853308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31805661) q[0];
sx q[0];
rz(-0.41724351) q[0];
sx q[0];
rz(0.70660025) q[0];
rz(0.63287863) q[1];
sx q[1];
rz(-1.0989847) q[1];
sx q[1];
rz(2.852827) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59290409) q[0];
sx q[0];
rz(-0.0026772896) q[0];
sx q[0];
rz(0.38627465) q[0];
rz(-pi) q[1];
rz(2.9339004) q[2];
sx q[2];
rz(-1.7294267) q[2];
sx q[2];
rz(-0.52244782) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0763467) q[1];
sx q[1];
rz(-1.7829121) q[1];
sx q[1];
rz(-0.29233039) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4411753) q[3];
sx q[3];
rz(-2.3312093) q[3];
sx q[3];
rz(-1.885351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67549813) q[2];
sx q[2];
rz(-2.0530632) q[2];
sx q[2];
rz(-1.4906073) q[2];
rz(-2.364482) q[3];
sx q[3];
rz(-1.4584352) q[3];
sx q[3];
rz(-1.3597663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8260088) q[0];
sx q[0];
rz(-1.4588139) q[0];
sx q[0];
rz(0.43651954) q[0];
rz(-2.3468158) q[1];
sx q[1];
rz(-1.3705285) q[1];
sx q[1];
rz(-2.2025542) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90715599) q[0];
sx q[0];
rz(-2.5946989) q[0];
sx q[0];
rz(2.7476385) q[0];
rz(-pi) q[1];
rz(-2.5592778) q[2];
sx q[2];
rz(-2.2852201) q[2];
sx q[2];
rz(-2.0272875) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.446167) q[1];
sx q[1];
rz(-1.0040196) q[1];
sx q[1];
rz(1.9480463) q[1];
x q[2];
rz(1.9179929) q[3];
sx q[3];
rz(-2.5592862) q[3];
sx q[3];
rz(0.83707419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7006435) q[2];
sx q[2];
rz(-1.0627397) q[2];
sx q[2];
rz(-2.5062594) q[2];
rz(0.82713953) q[3];
sx q[3];
rz(-0.45378903) q[3];
sx q[3];
rz(1.2686096) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8968935) q[0];
sx q[0];
rz(-0.06299717) q[0];
sx q[0];
rz(0.52198207) q[0];
rz(-2.7105647) q[1];
sx q[1];
rz(-1.6072175) q[1];
sx q[1];
rz(-2.4897051) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7303607) q[0];
sx q[0];
rz(-1.8069977) q[0];
sx q[0];
rz(-2.622261) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4944369) q[2];
sx q[2];
rz(-2.2324413) q[2];
sx q[2];
rz(2.5376157) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.401571) q[1];
sx q[1];
rz(-1.5844939) q[1];
sx q[1];
rz(-1.1915951) q[1];
rz(-pi) q[2];
rz(2.1736657) q[3];
sx q[3];
rz(-2.180392) q[3];
sx q[3];
rz(2.2948007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7151457) q[2];
sx q[2];
rz(-1.8603674) q[2];
sx q[2];
rz(-2.35671) q[2];
rz(0.52810413) q[3];
sx q[3];
rz(-1.2930861) q[3];
sx q[3];
rz(1.2877134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2454979) q[0];
sx q[0];
rz(-0.69480768) q[0];
sx q[0];
rz(2.3520663) q[0];
rz(2.6441669) q[1];
sx q[1];
rz(-1.0405633) q[1];
sx q[1];
rz(-1.8400037) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3606982) q[0];
sx q[0];
rz(-0.62735617) q[0];
sx q[0];
rz(-2.832546) q[0];
rz(0.17080266) q[2];
sx q[2];
rz(-0.34366648) q[2];
sx q[2];
rz(-0.91533585) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6005206) q[1];
sx q[1];
rz(-0.75319911) q[1];
sx q[1];
rz(2.888117) q[1];
x q[2];
rz(1.2034802) q[3];
sx q[3];
rz(-1.5919893) q[3];
sx q[3];
rz(-0.68875865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5975981) q[2];
sx q[2];
rz(-1.602729) q[2];
sx q[2];
rz(-0.27628118) q[2];
rz(2.1918519) q[3];
sx q[3];
rz(-2.8730928) q[3];
sx q[3];
rz(-2.5883519) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62040579) q[0];
sx q[0];
rz(-3.1053472) q[0];
sx q[0];
rz(-2.1976443) q[0];
rz(1.8432519) q[1];
sx q[1];
rz(-1.0669758) q[1];
sx q[1];
rz(2.3640769) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8572814) q[0];
sx q[0];
rz(-2.4696283) q[0];
sx q[0];
rz(2.551033) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8328415) q[2];
sx q[2];
rz(-2.2876559) q[2];
sx q[2];
rz(1.0345296) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7878163) q[1];
sx q[1];
rz(-2.750038) q[1];
sx q[1];
rz(-2.8038674) q[1];
rz(-0.30413587) q[3];
sx q[3];
rz(-1.602293) q[3];
sx q[3];
rz(-0.57284063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.62766084) q[2];
sx q[2];
rz(-1.3769423) q[2];
sx q[2];
rz(-2.7505006) q[2];
rz(2.5202461) q[3];
sx q[3];
rz(-0.73453271) q[3];
sx q[3];
rz(-2.3656316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4710627) q[0];
sx q[0];
rz(-0.10469086) q[0];
sx q[0];
rz(1.4290357) q[0];
rz(-0.96587005) q[1];
sx q[1];
rz(-1.254225) q[1];
sx q[1];
rz(-0.78580725) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87596506) q[0];
sx q[0];
rz(-2.4959491) q[0];
sx q[0];
rz(0.085787817) q[0];
x q[1];
rz(-0.73432095) q[2];
sx q[2];
rz(-1.87748) q[2];
sx q[2];
rz(1.1994565) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.25868402) q[1];
sx q[1];
rz(-1.6526387) q[1];
sx q[1];
rz(2.6805582) q[1];
rz(-2.7034055) q[3];
sx q[3];
rz(-1.8343958) q[3];
sx q[3];
rz(0.98942843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.69661951) q[2];
sx q[2];
rz(-2.5225621) q[2];
sx q[2];
rz(2.6444198) q[2];
rz(-2.2677126) q[3];
sx q[3];
rz(-2.0495575) q[3];
sx q[3];
rz(-1.2018275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1931605) q[0];
sx q[0];
rz(-0.30696294) q[0];
sx q[0];
rz(0.69751414) q[0];
rz(0.3802309) q[1];
sx q[1];
rz(-1.1153328) q[1];
sx q[1];
rz(3.0013705) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8070458) q[0];
sx q[0];
rz(-1.3893034) q[0];
sx q[0];
rz(-1.6491778) q[0];
rz(0.44084844) q[2];
sx q[2];
rz(-1.5356283) q[2];
sx q[2];
rz(2.9174385) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3933691) q[1];
sx q[1];
rz(-2.2294982) q[1];
sx q[1];
rz(-1.0568807) q[1];
rz(-1.3789375) q[3];
sx q[3];
rz(-1.1664806) q[3];
sx q[3];
rz(-1.9003001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2908638) q[2];
sx q[2];
rz(-1.2173434) q[2];
sx q[2];
rz(-0.49120894) q[2];
rz(0.20719191) q[3];
sx q[3];
rz(-2.5184293) q[3];
sx q[3];
rz(1.7306958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9687013) q[0];
sx q[0];
rz(-1.4145114) q[0];
sx q[0];
rz(0.15039314) q[0];
rz(-2.4405759) q[1];
sx q[1];
rz(-2.0471408) q[1];
sx q[1];
rz(1.4775803) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1883066) q[0];
sx q[0];
rz(-2.4948848) q[0];
sx q[0];
rz(-2.9394399) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6871181) q[2];
sx q[2];
rz(-1.4106299) q[2];
sx q[2];
rz(-0.016906658) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3353053) q[1];
sx q[1];
rz(-1.6769969) q[1];
sx q[1];
rz(-0.152529) q[1];
rz(2.8692352) q[3];
sx q[3];
rz(-2.2873169) q[3];
sx q[3];
rz(2.2509433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.56665862) q[2];
sx q[2];
rz(-0.40859544) q[2];
sx q[2];
rz(-1.6179786) q[2];
rz(0.59897113) q[3];
sx q[3];
rz(-0.50061575) q[3];
sx q[3];
rz(2.9536501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69342518) q[0];
sx q[0];
rz(-0.42891112) q[0];
sx q[0];
rz(1.6397788) q[0];
rz(2.2046454) q[1];
sx q[1];
rz(-2.1138771) q[1];
sx q[1];
rz(-1.017259) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.992998) q[0];
sx q[0];
rz(-2.5879597) q[0];
sx q[0];
rz(-0.77405907) q[0];
rz(1.6566562) q[2];
sx q[2];
rz(-1.4659662) q[2];
sx q[2];
rz(-1.3121999) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76495095) q[1];
sx q[1];
rz(-1.3089453) q[1];
sx q[1];
rz(0.034305926) q[1];
rz(0.42933257) q[3];
sx q[3];
rz(-1.7547914) q[3];
sx q[3];
rz(-1.7443378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1169869) q[2];
sx q[2];
rz(-2.459343) q[2];
sx q[2];
rz(1.5218706) q[2];
rz(1.4874602) q[3];
sx q[3];
rz(-1.6493075) q[3];
sx q[3];
rz(-1.7447932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37987729) q[0];
sx q[0];
rz(-1.7424255) q[0];
sx q[0];
rz(-2.4095834) q[0];
rz(-1.4749745) q[1];
sx q[1];
rz(-1.2154308) q[1];
sx q[1];
rz(-2.348127) q[1];
rz(-0.39948612) q[2];
sx q[2];
rz(-2.1151092) q[2];
sx q[2];
rz(-1.5110037) q[2];
rz(1.2670682) q[3];
sx q[3];
rz(-1.8770185) q[3];
sx q[3];
rz(0.087531623) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

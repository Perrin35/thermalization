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
rz(-2.0877617) q[0];
sx q[0];
rz(-0.55019903) q[0];
sx q[0];
rz(1.7734111) q[0];
rz(-0.47222459) q[1];
sx q[1];
rz(-3.0250186) q[1];
sx q[1];
rz(-0.20224686) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7769788) q[0];
sx q[0];
rz(-1.783051) q[0];
sx q[0];
rz(-1.3670761) q[0];
x q[1];
rz(0.80533452) q[2];
sx q[2];
rz(-0.67395681) q[2];
sx q[2];
rz(-1.0588624) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3319407) q[1];
sx q[1];
rz(-2.1686848) q[1];
sx q[1];
rz(-2.7974068) q[1];
rz(-1.6458166) q[3];
sx q[3];
rz(-0.46722355) q[3];
sx q[3];
rz(0.4969767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5998485) q[2];
sx q[2];
rz(-2.2969963) q[2];
sx q[2];
rz(2.2005626) q[2];
rz(-2.8372676) q[3];
sx q[3];
rz(-2.2797238) q[3];
sx q[3];
rz(-2.8743675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.32352725) q[0];
sx q[0];
rz(-0.37779385) q[0];
sx q[0];
rz(-0.7793119) q[0];
rz(-1.3620954) q[1];
sx q[1];
rz(-2.7423488) q[1];
sx q[1];
rz(2.0077226) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88863605) q[0];
sx q[0];
rz(-2.0042814) q[0];
sx q[0];
rz(-0.10826464) q[0];
x q[1];
rz(-0.46102779) q[2];
sx q[2];
rz(-3.010058) q[2];
sx q[2];
rz(1.6361089) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.081959978) q[1];
sx q[1];
rz(-1.6784918) q[1];
sx q[1];
rz(1.3038941) q[1];
rz(1.5690345) q[3];
sx q[3];
rz(-0.50927466) q[3];
sx q[3];
rz(-2.5363896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69089943) q[2];
sx q[2];
rz(-2.1096114) q[2];
sx q[2];
rz(2.9449985) q[2];
rz(-2.172566) q[3];
sx q[3];
rz(-1.5195547) q[3];
sx q[3];
rz(1.903418) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3900688) q[0];
sx q[0];
rz(-0.5539493) q[0];
sx q[0];
rz(2.4098136) q[0];
rz(-2.7732908) q[1];
sx q[1];
rz(-2.383547) q[1];
sx q[1];
rz(-2.7751353) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0822507) q[0];
sx q[0];
rz(-1.8346795) q[0];
sx q[0];
rz(1.7543955) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3442893) q[2];
sx q[2];
rz(-1.8069805) q[2];
sx q[2];
rz(3.0677441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7190105) q[1];
sx q[1];
rz(-1.9937464) q[1];
sx q[1];
rz(-2.1450536) q[1];
x q[2];
rz(0.59753527) q[3];
sx q[3];
rz(-0.96216494) q[3];
sx q[3];
rz(0.017692117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1825819) q[2];
sx q[2];
rz(-2.3704447) q[2];
sx q[2];
rz(0.89243531) q[2];
rz(0.45448947) q[3];
sx q[3];
rz(-0.82388866) q[3];
sx q[3];
rz(-0.23196001) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7682122) q[0];
sx q[0];
rz(-0.1425655) q[0];
sx q[0];
rz(-2.7349803) q[0];
rz(-2.3136102) q[1];
sx q[1];
rz(-1.4031289) q[1];
sx q[1];
rz(0.83555317) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8260561) q[0];
sx q[0];
rz(-2.9565507) q[0];
sx q[0];
rz(1.7868143) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0876483) q[2];
sx q[2];
rz(-2.7479541) q[2];
sx q[2];
rz(-2.7700499) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4728188) q[1];
sx q[1];
rz(-1.4240648) q[1];
sx q[1];
rz(0.43763486) q[1];
x q[2];
rz(1.9560247) q[3];
sx q[3];
rz(-1.480417) q[3];
sx q[3];
rz(-2.4848695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39997175) q[2];
sx q[2];
rz(-2.8394832) q[2];
sx q[2];
rz(2.2201904) q[2];
rz(1.7737927) q[3];
sx q[3];
rz(-0.96145815) q[3];
sx q[3];
rz(-3.0000946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2530186) q[0];
sx q[0];
rz(-1.0971917) q[0];
sx q[0];
rz(-0.15884037) q[0];
rz(-1.6244434) q[1];
sx q[1];
rz(-0.30888638) q[1];
sx q[1];
rz(-1.812017) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3811262) q[0];
sx q[0];
rz(-1.5601349) q[0];
sx q[0];
rz(0.018254781) q[0];
rz(-pi) q[1];
rz(0.23933584) q[2];
sx q[2];
rz(-2.1194601) q[2];
sx q[2];
rz(-1.885883) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98011866) q[1];
sx q[1];
rz(-1.9494281) q[1];
sx q[1];
rz(-2.9254854) q[1];
rz(2.9333889) q[3];
sx q[3];
rz(-2.1284687) q[3];
sx q[3];
rz(-0.69277908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.57881957) q[2];
sx q[2];
rz(-0.69391888) q[2];
sx q[2];
rz(0.55207437) q[2];
rz(-0.79365802) q[3];
sx q[3];
rz(-0.59760439) q[3];
sx q[3];
rz(-2.1551267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3119222) q[0];
sx q[0];
rz(-2.2770918) q[0];
sx q[0];
rz(2.435834) q[0];
rz(-0.13748473) q[1];
sx q[1];
rz(-0.64422137) q[1];
sx q[1];
rz(0.71298832) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5891083) q[0];
sx q[0];
rz(-2.6875067) q[0];
sx q[0];
rz(0.54291351) q[0];
x q[1];
rz(-2.627264) q[2];
sx q[2];
rz(-1.1900969) q[2];
sx q[2];
rz(1.9388388) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2299774) q[1];
sx q[1];
rz(-1.4739081) q[1];
sx q[1];
rz(-1.1712892) q[1];
x q[2];
rz(3.0253382) q[3];
sx q[3];
rz(-1.2796806) q[3];
sx q[3];
rz(1.8965885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2116427) q[2];
sx q[2];
rz(-0.86923081) q[2];
sx q[2];
rz(-0.28406528) q[2];
rz(0.12187135) q[3];
sx q[3];
rz(-0.28212306) q[3];
sx q[3];
rz(1.4819283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.876494) q[0];
sx q[0];
rz(-0.6186741) q[0];
sx q[0];
rz(0.70736831) q[0];
rz(0.44037285) q[1];
sx q[1];
rz(-2.423954) q[1];
sx q[1];
rz(-1.7792938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.52584) q[0];
sx q[0];
rz(-2.7670015) q[0];
sx q[0];
rz(-2.2994141) q[0];
rz(-pi) q[1];
rz(-2.0411891) q[2];
sx q[2];
rz(-2.3572323) q[2];
sx q[2];
rz(0.91299865) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3063115) q[1];
sx q[1];
rz(-2.0871442) q[1];
sx q[1];
rz(2.5912214) q[1];
rz(1.4555638) q[3];
sx q[3];
rz(-1.1172022) q[3];
sx q[3];
rz(0.99059243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97544396) q[2];
sx q[2];
rz(-2.7232309) q[2];
sx q[2];
rz(0.56900209) q[2];
rz(1.0373908) q[3];
sx q[3];
rz(-2.2834957) q[3];
sx q[3];
rz(2.67498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5910832) q[0];
sx q[0];
rz(-0.35174462) q[0];
sx q[0];
rz(1.2048703) q[0];
rz(2.6288746) q[1];
sx q[1];
rz(-0.7690438) q[1];
sx q[1];
rz(-2.9606294) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34535392) q[0];
sx q[0];
rz(-1.1974044) q[0];
sx q[0];
rz(-0.42723421) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7788073) q[2];
sx q[2];
rz(-2.7698667) q[2];
sx q[2];
rz(-0.065601018) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0389591) q[1];
sx q[1];
rz(-0.60188434) q[1];
sx q[1];
rz(-1.9485056) q[1];
x q[2];
rz(0.11060235) q[3];
sx q[3];
rz(-1.6440653) q[3];
sx q[3];
rz(1.6694809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5391431) q[2];
sx q[2];
rz(-0.64843833) q[2];
sx q[2];
rz(0.096972018) q[2];
rz(-0.55475956) q[3];
sx q[3];
rz(-1.8192889) q[3];
sx q[3];
rz(0.35840148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.48034126) q[0];
sx q[0];
rz(-2.3619409) q[0];
sx q[0];
rz(3.0338147) q[0];
rz(0.55094552) q[1];
sx q[1];
rz(-0.44083732) q[1];
sx q[1];
rz(-1.617618) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0227528) q[0];
sx q[0];
rz(-1.700239) q[0];
sx q[0];
rz(1.8047099) q[0];
x q[1];
rz(1.0650159) q[2];
sx q[2];
rz(-0.98978251) q[2];
sx q[2];
rz(2.6083824) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0131993) q[1];
sx q[1];
rz(-2.6775762) q[1];
sx q[1];
rz(1.6676519) q[1];
rz(3.0436936) q[3];
sx q[3];
rz(-0.49360858) q[3];
sx q[3];
rz(0.61033397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23065755) q[2];
sx q[2];
rz(-0.50369889) q[2];
sx q[2];
rz(2.0111734) q[2];
rz(0.23760992) q[3];
sx q[3];
rz(-2.7311324) q[3];
sx q[3];
rz(-2.3835278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1199353) q[0];
sx q[0];
rz(-2.9627242) q[0];
sx q[0];
rz(-2.2005431) q[0];
rz(2.7811116) q[1];
sx q[1];
rz(-1.4070114) q[1];
sx q[1];
rz(1.2233268) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4775765) q[0];
sx q[0];
rz(-1.7682791) q[0];
sx q[0];
rz(-2.5504179) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98127301) q[2];
sx q[2];
rz(-2.2726577) q[2];
sx q[2];
rz(-1.5379932) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2625433) q[1];
sx q[1];
rz(-1.9886262) q[1];
sx q[1];
rz(1.6476739) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1072708) q[3];
sx q[3];
rz(-1.4015028) q[3];
sx q[3];
rz(-0.062549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.90784043) q[2];
sx q[2];
rz(-1.3774104) q[2];
sx q[2];
rz(3.0506548) q[2];
rz(-2.9371373) q[3];
sx q[3];
rz(-2.3369868) q[3];
sx q[3];
rz(2.0596152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7627056) q[0];
sx q[0];
rz(-1.5305516) q[0];
sx q[0];
rz(-1.332921) q[0];
rz(-1.4270205) q[1];
sx q[1];
rz(-2.6418229) q[1];
sx q[1];
rz(-1.5508834) q[1];
rz(-1.0574404) q[2];
sx q[2];
rz(-2.3427137) q[2];
sx q[2];
rz(-2.8165934) q[2];
rz(-0.25973917) q[3];
sx q[3];
rz(-2.9208214) q[3];
sx q[3];
rz(-0.95820273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

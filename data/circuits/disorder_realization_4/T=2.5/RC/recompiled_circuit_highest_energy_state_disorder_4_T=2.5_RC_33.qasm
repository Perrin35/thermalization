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
rz(1.053831) q[0];
sx q[0];
rz(-2.5913936) q[0];
sx q[0];
rz(1.3681816) q[0];
rz(-0.47222459) q[1];
sx q[1];
rz(-3.0250186) q[1];
sx q[1];
rz(-0.20224686) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36461386) q[0];
sx q[0];
rz(-1.3585416) q[0];
sx q[0];
rz(1.7745166) q[0];
rz(1.0482685) q[2];
sx q[2];
rz(-2.0179581) q[2];
sx q[2];
rz(-1.1559522) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.579548) q[1];
sx q[1];
rz(-1.8534396) q[1];
sx q[1];
rz(-2.1971027) q[1];
x q[2];
rz(1.6458166) q[3];
sx q[3];
rz(-2.6743691) q[3];
sx q[3];
rz(-2.644616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5998485) q[2];
sx q[2];
rz(-0.84459633) q[2];
sx q[2];
rz(-2.2005626) q[2];
rz(0.30432501) q[3];
sx q[3];
rz(-2.2797238) q[3];
sx q[3];
rz(-2.8743675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8180654) q[0];
sx q[0];
rz(-2.7637988) q[0];
sx q[0];
rz(-2.3622808) q[0];
rz(-1.3620954) q[1];
sx q[1];
rz(-2.7423488) q[1];
sx q[1];
rz(-1.13387) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2529566) q[0];
sx q[0];
rz(-1.1373113) q[0];
sx q[0];
rz(3.033328) q[0];
rz(-3.0236565) q[2];
sx q[2];
rz(-1.6291766) q[2];
sx q[2];
rz(0.3922677) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2782005) q[1];
sx q[1];
rz(-0.28732936) q[1];
sx q[1];
rz(1.9598239) q[1];
rz(1.0615223) q[3];
sx q[3];
rz(-1.5716553) q[3];
sx q[3];
rz(-2.1775376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4506932) q[2];
sx q[2];
rz(-2.1096114) q[2];
sx q[2];
rz(2.9449985) q[2];
rz(-0.96902668) q[3];
sx q[3];
rz(-1.5195547) q[3];
sx q[3];
rz(-1.903418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(2.3900688) q[0];
sx q[0];
rz(-0.5539493) q[0];
sx q[0];
rz(0.73177904) q[0];
rz(-0.36830184) q[1];
sx q[1];
rz(-0.75804561) q[1];
sx q[1];
rz(-2.7751353) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059341934) q[0];
sx q[0];
rz(-1.8346795) q[0];
sx q[0];
rz(1.3871971) q[0];
x q[1];
rz(1.7973034) q[2];
sx q[2];
rz(-1.3346121) q[2];
sx q[2];
rz(0.073848595) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0302387) q[1];
sx q[1];
rz(-2.0890279) q[1];
sx q[1];
rz(-0.49211647) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59753527) q[3];
sx q[3];
rz(-0.96216494) q[3];
sx q[3];
rz(3.1239005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9590108) q[2];
sx q[2];
rz(-0.77114791) q[2];
sx q[2];
rz(2.2491573) q[2];
rz(-0.45448947) q[3];
sx q[3];
rz(-2.317704) q[3];
sx q[3];
rz(2.9096326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3733805) q[0];
sx q[0];
rz(-2.9990271) q[0];
sx q[0];
rz(2.7349803) q[0];
rz(-2.3136102) q[1];
sx q[1];
rz(-1.4031289) q[1];
sx q[1];
rz(-2.3060395) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8260561) q[0];
sx q[0];
rz(-0.18504194) q[0];
sx q[0];
rz(-1.3547784) q[0];
rz(-pi) q[1];
rz(2.0876483) q[2];
sx q[2];
rz(-2.7479541) q[2];
sx q[2];
rz(0.37154276) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4728188) q[1];
sx q[1];
rz(-1.4240648) q[1];
sx q[1];
rz(-0.43763486) q[1];
rz(-pi) q[2];
rz(1.1855679) q[3];
sx q[3];
rz(-1.480417) q[3];
sx q[3];
rz(2.4848695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7416209) q[2];
sx q[2];
rz(-2.8394832) q[2];
sx q[2];
rz(-0.92140222) q[2];
rz(1.3678) q[3];
sx q[3];
rz(-0.96145815) q[3];
sx q[3];
rz(-0.14149806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2530186) q[0];
sx q[0];
rz(-2.044401) q[0];
sx q[0];
rz(-0.15884037) q[0];
rz(1.6244434) q[1];
sx q[1];
rz(-0.30888638) q[1];
sx q[1];
rz(1.812017) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3388728) q[0];
sx q[0];
rz(-3.1204528) q[0];
sx q[0];
rz(0.52864023) q[0];
rz(-pi) q[1];
rz(0.23933584) q[2];
sx q[2];
rz(-2.1194601) q[2];
sx q[2];
rz(1.2557097) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6255505) q[1];
sx q[1];
rz(-0.4333638) q[1];
sx q[1];
rz(-2.0651555) q[1];
rz(1.8908126) q[3];
sx q[3];
rz(-2.5502) q[3];
sx q[3];
rz(-1.0726014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5627731) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3119222) q[0];
sx q[0];
rz(-2.2770918) q[0];
sx q[0];
rz(0.70575869) q[0];
rz(3.0041079) q[1];
sx q[1];
rz(-0.64422137) q[1];
sx q[1];
rz(-2.4286043) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1437837) q[0];
sx q[0];
rz(-1.1857872) q[0];
sx q[0];
rz(1.3237756) q[0];
rz(0.6829459) q[2];
sx q[2];
rz(-2.5120398) q[2];
sx q[2];
rz(2.1918639) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29999816) q[1];
sx q[1];
rz(-1.1732686) q[1];
sx q[1];
rz(-0.10511153) q[1];
rz(-pi) q[2];
rz(-1.2014404) q[3];
sx q[3];
rz(-2.8287362) q[3];
sx q[3];
rz(-2.2830056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9299499) q[2];
sx q[2];
rz(-0.86923081) q[2];
sx q[2];
rz(0.28406528) q[2];
rz(-0.12187135) q[3];
sx q[3];
rz(-0.28212306) q[3];
sx q[3];
rz(-1.4819283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26509869) q[0];
sx q[0];
rz(-0.6186741) q[0];
sx q[0];
rz(2.4342243) q[0];
rz(-0.44037285) q[1];
sx q[1];
rz(-0.71763867) q[1];
sx q[1];
rz(-1.7792938) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.52584) q[0];
sx q[0];
rz(-0.37459117) q[0];
sx q[0];
rz(-2.2994141) q[0];
rz(-pi) q[1];
rz(-0.42476219) q[2];
sx q[2];
rz(-0.889689) q[2];
sx q[2];
rz(-2.8515138) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55864298) q[1];
sx q[1];
rz(-1.0986277) q[1];
sx q[1];
rz(0.98319816) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23162095) q[3];
sx q[3];
rz(-0.4670139) q[3];
sx q[3];
rz(1.8927595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1661487) q[2];
sx q[2];
rz(-2.7232309) q[2];
sx q[2];
rz(0.56900209) q[2];
rz(-1.0373908) q[3];
sx q[3];
rz(-0.85809696) q[3];
sx q[3];
rz(-0.4666127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5910832) q[0];
sx q[0];
rz(-0.35174462) q[0];
sx q[0];
rz(-1.2048703) q[0];
rz(0.51271802) q[1];
sx q[1];
rz(-0.7690438) q[1];
sx q[1];
rz(2.9606294) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9009243) q[0];
sx q[0];
rz(-2.5818733) q[0];
sx q[0];
rz(0.75738917) q[0];
rz(1.4333192) q[2];
sx q[2];
rz(-1.2243234) q[2];
sx q[2];
rz(-0.32127831) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1026336) q[1];
sx q[1];
rz(-2.5397083) q[1];
sx q[1];
rz(-1.9485056) q[1];
x q[2];
rz(-1.6445141) q[3];
sx q[3];
rz(-1.6811007) q[3];
sx q[3];
rz(-3.0347787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5391431) q[2];
sx q[2];
rz(-2.4931543) q[2];
sx q[2];
rz(0.096972018) q[2];
rz(-0.55475956) q[3];
sx q[3];
rz(-1.8192889) q[3];
sx q[3];
rz(-2.7831912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48034126) q[0];
sx q[0];
rz(-0.77965176) q[0];
sx q[0];
rz(0.10777792) q[0];
rz(0.55094552) q[1];
sx q[1];
rz(-2.7007553) q[1];
sx q[1];
rz(-1.5239747) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.658889) q[0];
sx q[0];
rz(-1.3388757) q[0];
sx q[0];
rz(-3.0085682) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63568398) q[2];
sx q[2];
rz(-0.7506461) q[2];
sx q[2];
rz(0.25668609) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0131993) q[1];
sx q[1];
rz(-2.6775762) q[1];
sx q[1];
rz(1.4739407) q[1];
rz(-2.6499841) q[3];
sx q[3];
rz(-1.5244686) q[3];
sx q[3];
rz(0.87418782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23065755) q[2];
sx q[2];
rz(-0.50369889) q[2];
sx q[2];
rz(1.1304193) q[2];
rz(-0.23760992) q[3];
sx q[3];
rz(-2.7311324) q[3];
sx q[3];
rz(2.3835278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.0216574) q[0];
sx q[0];
rz(-2.9627242) q[0];
sx q[0];
rz(2.2005431) q[0];
rz(-0.36048105) q[1];
sx q[1];
rz(-1.7345813) q[1];
sx q[1];
rz(1.9182659) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4775765) q[0];
sx q[0];
rz(-1.3733136) q[0];
sx q[0];
rz(2.5504179) q[0];
x q[1];
rz(-2.5599096) q[2];
sx q[2];
rz(-0.8831199) q[2];
sx q[2];
rz(-0.73492848) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0666484) q[1];
sx q[1];
rz(-2.7171591) q[1];
sx q[1];
rz(0.17129274) q[1];
rz(2.9452865) q[3];
sx q[3];
rz(-2.098791) q[3];
sx q[3];
rz(1.7332026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2337522) q[2];
sx q[2];
rz(-1.7641822) q[2];
sx q[2];
rz(-3.0506548) q[2];
rz(2.9371373) q[3];
sx q[3];
rz(-0.8046059) q[3];
sx q[3];
rz(-1.0819775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37888708) q[0];
sx q[0];
rz(-1.611041) q[0];
sx q[0];
rz(1.8086717) q[0];
rz(1.7145722) q[1];
sx q[1];
rz(-2.6418229) q[1];
sx q[1];
rz(-1.5508834) q[1];
rz(2.3007921) q[2];
sx q[2];
rz(-1.9304095) q[2];
sx q[2];
rz(-0.87113397) q[2];
rz(1.5132202) q[3];
sx q[3];
rz(-1.3575469) q[3];
sx q[3];
rz(-1.224106) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

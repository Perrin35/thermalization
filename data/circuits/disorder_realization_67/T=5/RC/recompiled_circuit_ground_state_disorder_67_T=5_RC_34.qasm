OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23842731) q[0];
sx q[0];
rz(-2.981346) q[0];
sx q[0];
rz(1.2877553) q[0];
rz(2.2924478) q[1];
sx q[1];
rz(1.1916817) q[1];
sx q[1];
rz(10.183505) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7641974) q[0];
sx q[0];
rz(-1.281651) q[0];
sx q[0];
rz(-1.3290861) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0082194) q[2];
sx q[2];
rz(-2.6133399) q[2];
sx q[2];
rz(-2.4932441) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1114486) q[1];
sx q[1];
rz(-1.2093102) q[1];
sx q[1];
rz(0.8059146) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2686602) q[3];
sx q[3];
rz(-2.2651787) q[3];
sx q[3];
rz(-2.5257655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8410926) q[2];
sx q[2];
rz(-0.79885834) q[2];
sx q[2];
rz(-0.50660261) q[2];
rz(-1.5233585) q[3];
sx q[3];
rz(-2.5169499) q[3];
sx q[3];
rz(0.055559572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27187207) q[0];
sx q[0];
rz(-0.46590081) q[0];
sx q[0];
rz(-3.077935) q[0];
rz(1.1977389) q[1];
sx q[1];
rz(-1.5143737) q[1];
sx q[1];
rz(-1.0187842) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5658582) q[0];
sx q[0];
rz(-1.3745148) q[0];
sx q[0];
rz(-1.422276) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5449824) q[2];
sx q[2];
rz(-2.7368174) q[2];
sx q[2];
rz(0.34133729) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6117759) q[1];
sx q[1];
rz(-0.88068286) q[1];
sx q[1];
rz(-1.2089096) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7131931) q[3];
sx q[3];
rz(-2.3835858) q[3];
sx q[3];
rz(0.49331323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7055052) q[2];
sx q[2];
rz(-0.59385308) q[2];
sx q[2];
rz(0.84211055) q[2];
rz(2.0841133) q[3];
sx q[3];
rz(-0.92856854) q[3];
sx q[3];
rz(-1.1697945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.24404003) q[0];
sx q[0];
rz(-2.121448) q[0];
sx q[0];
rz(1.9496339) q[0];
rz(-0.92487088) q[1];
sx q[1];
rz(-2.4332739) q[1];
sx q[1];
rz(1.8398197) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98604964) q[0];
sx q[0];
rz(-1.0240004) q[0];
sx q[0];
rz(-0.56122551) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5315751) q[2];
sx q[2];
rz(-1.7698145) q[2];
sx q[2];
rz(1.7165556) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5435346) q[1];
sx q[1];
rz(-2.638276) q[1];
sx q[1];
rz(2.7642168) q[1];
rz(-pi) q[2];
rz(-2.202192) q[3];
sx q[3];
rz(-0.41920788) q[3];
sx q[3];
rz(-1.7260176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3785582) q[2];
sx q[2];
rz(-0.083746567) q[2];
sx q[2];
rz(0.087510022) q[2];
rz(1.3148426) q[3];
sx q[3];
rz(-2.0851236) q[3];
sx q[3];
rz(-0.99353138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73115504) q[0];
sx q[0];
rz(-1.9696099) q[0];
sx q[0];
rz(0.80899578) q[0];
rz(1.0357098) q[1];
sx q[1];
rz(-0.46329841) q[1];
sx q[1];
rz(2.1083924) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22543487) q[0];
sx q[0];
rz(-1.5317129) q[0];
sx q[0];
rz(-0.90644159) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0339917) q[2];
sx q[2];
rz(-1.9516757) q[2];
sx q[2];
rz(-1.8465259) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39604998) q[1];
sx q[1];
rz(-2.0337359) q[1];
sx q[1];
rz(-1.7718503) q[1];
rz(2.4304588) q[3];
sx q[3];
rz(-0.67340771) q[3];
sx q[3];
rz(-2.8720958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0674627) q[2];
sx q[2];
rz(-0.16572696) q[2];
sx q[2];
rz(1.6039675) q[2];
rz(-1.7240723) q[3];
sx q[3];
rz(-2.0310903) q[3];
sx q[3];
rz(-1.0874776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9017225) q[0];
sx q[0];
rz(-2.9481695) q[0];
sx q[0];
rz(1.9523917) q[0];
rz(2.4173648) q[1];
sx q[1];
rz(-1.8502356) q[1];
sx q[1];
rz(-1.6667746) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9928707) q[0];
sx q[0];
rz(-1.9970511) q[0];
sx q[0];
rz(2.5919586) q[0];
rz(1.4101661) q[2];
sx q[2];
rz(-2.361627) q[2];
sx q[2];
rz(-2.6218564) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2853299) q[1];
sx q[1];
rz(-1.8736427) q[1];
sx q[1];
rz(-1.21638) q[1];
rz(-0.46012043) q[3];
sx q[3];
rz(-2.0232206) q[3];
sx q[3];
rz(0.89738536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7597947) q[2];
sx q[2];
rz(-0.90961027) q[2];
sx q[2];
rz(-0.46662113) q[2];
rz(-1.7032547) q[3];
sx q[3];
rz(-1.2983026) q[3];
sx q[3];
rz(-0.94021016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.4284215) q[0];
sx q[0];
rz(-2.3155825) q[0];
sx q[0];
rz(-0.010350479) q[0];
rz(-1.6185282) q[1];
sx q[1];
rz(-1.7489988) q[1];
sx q[1];
rz(-1.3759605) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15555602) q[0];
sx q[0];
rz(-1.2335834) q[0];
sx q[0];
rz(-0.2561148) q[0];
x q[1];
rz(-2.544246) q[2];
sx q[2];
rz(-0.22946363) q[2];
sx q[2];
rz(-0.84537431) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79094564) q[1];
sx q[1];
rz(-0.53574359) q[1];
sx q[1];
rz(-2.78113) q[1];
rz(0.67988396) q[3];
sx q[3];
rz(-2.2671642) q[3];
sx q[3];
rz(-1.0395466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7213664) q[2];
sx q[2];
rz(-2.0927636) q[2];
sx q[2];
rz(2.4143207) q[2];
rz(0.51491245) q[3];
sx q[3];
rz(-1.3313096) q[3];
sx q[3];
rz(1.7236727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0338243) q[0];
sx q[0];
rz(-1.533778) q[0];
sx q[0];
rz(0.4069826) q[0];
rz(2.175323) q[1];
sx q[1];
rz(-0.85710183) q[1];
sx q[1];
rz(-3.0028717) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7329368) q[0];
sx q[0];
rz(-1.4545024) q[0];
sx q[0];
rz(3.1058806) q[0];
x q[1];
rz(-2.474276) q[2];
sx q[2];
rz(-2.7936068) q[2];
sx q[2];
rz(1.1225357) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.87200621) q[1];
sx q[1];
rz(-1.2254189) q[1];
sx q[1];
rz(-2.7589382) q[1];
rz(-pi) q[2];
rz(-1.9646264) q[3];
sx q[3];
rz(-2.7714249) q[3];
sx q[3];
rz(-2.1578876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2762642) q[2];
sx q[2];
rz(-1.7636969) q[2];
sx q[2];
rz(-2.830937) q[2];
rz(-1.3079414) q[3];
sx q[3];
rz(-1.5111978) q[3];
sx q[3];
rz(0.37180296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0528316) q[0];
sx q[0];
rz(-0.3449057) q[0];
sx q[0];
rz(-0.088223591) q[0];
rz(-3.131033) q[1];
sx q[1];
rz(-1.5225531) q[1];
sx q[1];
rz(-0.47058502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39401606) q[0];
sx q[0];
rz(-0.86578275) q[0];
sx q[0];
rz(2.2590619) q[0];
x q[1];
rz(-2.3345474) q[2];
sx q[2];
rz(-1.4331487) q[2];
sx q[2];
rz(-1.1040083) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6424017) q[1];
sx q[1];
rz(-0.27398738) q[1];
sx q[1];
rz(-1.8887331) q[1];
rz(-pi) q[2];
rz(-0.13981522) q[3];
sx q[3];
rz(-1.2212409) q[3];
sx q[3];
rz(-2.4619554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18870246) q[2];
sx q[2];
rz(-2.1850093) q[2];
sx q[2];
rz(1.0799705) q[2];
rz(-0.88932577) q[3];
sx q[3];
rz(-1.5417128) q[3];
sx q[3];
rz(-1.5573474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.476986) q[0];
sx q[0];
rz(-2.7598858) q[0];
sx q[0];
rz(-0.81740022) q[0];
rz(0.74642247) q[1];
sx q[1];
rz(-2.4669929) q[1];
sx q[1];
rz(-0.93961632) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8864) q[0];
sx q[0];
rz(-1.5536683) q[0];
sx q[0];
rz(-0.00047943133) q[0];
rz(-pi) q[1];
rz(-2.7857271) q[2];
sx q[2];
rz(-1.8345222) q[2];
sx q[2];
rz(3.0413994) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47245644) q[1];
sx q[1];
rz(-1.3208773) q[1];
sx q[1];
rz(-3.0210698) q[1];
rz(0.21454105) q[3];
sx q[3];
rz(-1.166178) q[3];
sx q[3];
rz(1.4264533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8115936) q[2];
sx q[2];
rz(-2.1758695) q[2];
sx q[2];
rz(-0.3453556) q[2];
rz(-1.5251478) q[3];
sx q[3];
rz(-1.350178) q[3];
sx q[3];
rz(2.7143872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78373194) q[0];
sx q[0];
rz(-1.9851728) q[0];
sx q[0];
rz(-0.082948908) q[0];
rz(-2.9792765) q[1];
sx q[1];
rz(-1.6971308) q[1];
sx q[1];
rz(2.4457795) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7617924) q[0];
sx q[0];
rz(-1.6125935) q[0];
sx q[0];
rz(-2.8618852) q[0];
rz(-1.4415353) q[2];
sx q[2];
rz(-1.5678355) q[2];
sx q[2];
rz(1.2128613) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79244991) q[1];
sx q[1];
rz(-0.44466838) q[1];
sx q[1];
rz(-0.044388958) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8105177) q[3];
sx q[3];
rz(-0.73270117) q[3];
sx q[3];
rz(-2.1520881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42056981) q[2];
sx q[2];
rz(-0.15779725) q[2];
sx q[2];
rz(-1.2201307) q[2];
rz(0.57341352) q[3];
sx q[3];
rz(-1.3308595) q[3];
sx q[3];
rz(0.24513182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86380105) q[0];
sx q[0];
rz(-1.6391123) q[0];
sx q[0];
rz(3.0085051) q[0];
rz(-0.0028903891) q[1];
sx q[1];
rz(-2.7255701) q[1];
sx q[1];
rz(2.4400673) q[1];
rz(0.9695644) q[2];
sx q[2];
rz(-2.8429016) q[2];
sx q[2];
rz(0.071879172) q[2];
rz(1.6297798) q[3];
sx q[3];
rz(-2.3702757) q[3];
sx q[3];
rz(2.6487917) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

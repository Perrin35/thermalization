OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8393505) q[0];
sx q[0];
rz(-1.8678764) q[0];
sx q[0];
rz(-0.94925517) q[0];
rz(-1.0957837) q[1];
sx q[1];
rz(-0.97691184) q[1];
sx q[1];
rz(1.7916602) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.183584) q[0];
sx q[0];
rz(-1.6537979) q[0];
sx q[0];
rz(0.50358332) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40387965) q[2];
sx q[2];
rz(-0.35940659) q[2];
sx q[2];
rz(-0.071381005) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9567555) q[1];
sx q[1];
rz(-0.3809027) q[1];
sx q[1];
rz(0.18397267) q[1];
rz(-pi) q[2];
x q[2];
rz(1.024035) q[3];
sx q[3];
rz(-2.6461678) q[3];
sx q[3];
rz(1.6273496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6076516) q[2];
sx q[2];
rz(-1.093163) q[2];
sx q[2];
rz(2.9489813) q[2];
rz(1.2827778) q[3];
sx q[3];
rz(-1.1056113) q[3];
sx q[3];
rz(0.58853308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.823536) q[0];
sx q[0];
rz(-0.41724351) q[0];
sx q[0];
rz(2.4349924) q[0];
rz(0.63287863) q[1];
sx q[1];
rz(-1.0989847) q[1];
sx q[1];
rz(-0.28876567) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97917999) q[0];
sx q[0];
rz(-1.5683163) q[0];
sx q[0];
rz(-1.5697877) q[0];
rz(-pi) q[1];
rz(0.20769221) q[2];
sx q[2];
rz(-1.4121659) q[2];
sx q[2];
rz(-0.52244782) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.56882492) q[1];
sx q[1];
rz(-1.8563885) q[1];
sx q[1];
rz(-1.3495803) q[1];
x q[2];
rz(1.1488647) q[3];
sx q[3];
rz(-0.85638085) q[3];
sx q[3];
rz(-1.8569291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67549813) q[2];
sx q[2];
rz(-2.0530632) q[2];
sx q[2];
rz(-1.6509854) q[2];
rz(2.364482) q[3];
sx q[3];
rz(-1.4584352) q[3];
sx q[3];
rz(1.3597663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31558388) q[0];
sx q[0];
rz(-1.6827787) q[0];
sx q[0];
rz(2.7050731) q[0];
rz(-2.3468158) q[1];
sx q[1];
rz(-1.3705285) q[1];
sx q[1];
rz(-2.2025542) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3600814) q[0];
sx q[0];
rz(-1.0699125) q[0];
sx q[0];
rz(1.8003746) q[0];
rz(-2.5592778) q[2];
sx q[2];
rz(-2.2852201) q[2];
sx q[2];
rz(-2.0272875) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.69542563) q[1];
sx q[1];
rz(-2.137573) q[1];
sx q[1];
rz(-1.9480463) q[1];
rz(-pi) q[2];
rz(-1.2235997) q[3];
sx q[3];
rz(-0.58230647) q[3];
sx q[3];
rz(-0.83707419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.44094917) q[2];
sx q[2];
rz(-1.0627397) q[2];
sx q[2];
rz(0.6353333) q[2];
rz(2.3144531) q[3];
sx q[3];
rz(-2.6878036) q[3];
sx q[3];
rz(1.2686096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8968935) q[0];
sx q[0];
rz(-0.06299717) q[0];
sx q[0];
rz(-2.6196106) q[0];
rz(-0.43102795) q[1];
sx q[1];
rz(-1.5343752) q[1];
sx q[1];
rz(0.65188754) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2289705) q[0];
sx q[0];
rz(-0.56601277) q[0];
sx q[0];
rz(-0.45155756) q[0];
rz(-pi) q[1];
rz(0.84649936) q[2];
sx q[2];
rz(-1.9545467) q[2];
sx q[2];
rz(-0.64696128) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.74002162) q[1];
sx q[1];
rz(-1.5844939) q[1];
sx q[1];
rz(-1.1915951) q[1];
rz(2.4383955) q[3];
sx q[3];
rz(-1.0872989) q[3];
sx q[3];
rz(2.0421703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.42644694) q[2];
sx q[2];
rz(-1.8603674) q[2];
sx q[2];
rz(0.78488266) q[2];
rz(-0.52810413) q[3];
sx q[3];
rz(-1.2930861) q[3];
sx q[3];
rz(-1.2877134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2454979) q[0];
sx q[0];
rz(-2.446785) q[0];
sx q[0];
rz(0.78952638) q[0];
rz(-0.49742571) q[1];
sx q[1];
rz(-1.0405633) q[1];
sx q[1];
rz(-1.8400037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042846366) q[0];
sx q[0];
rz(-1.7502971) q[0];
sx q[0];
rz(0.60447201) q[0];
x q[1];
rz(2.80255) q[2];
sx q[2];
rz(-1.6280988) q[2];
sx q[2];
rz(0.49446854) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.54107202) q[1];
sx q[1];
rz(-2.3883935) q[1];
sx q[1];
rz(-0.25347565) q[1];
x q[2];
rz(-0.022707247) q[3];
sx q[3];
rz(-1.2035666) q[3];
sx q[3];
rz(-2.2677088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5975981) q[2];
sx q[2];
rz(-1.5388637) q[2];
sx q[2];
rz(0.27628118) q[2];
rz(2.1918519) q[3];
sx q[3];
rz(-2.8730928) q[3];
sx q[3];
rz(0.55324078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5211869) q[0];
sx q[0];
rz(-3.1053472) q[0];
sx q[0];
rz(0.94394839) q[0];
rz(1.2983407) q[1];
sx q[1];
rz(-2.0746168) q[1];
sx q[1];
rz(2.3640769) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3382638) q[0];
sx q[0];
rz(-1.9247806) q[0];
sx q[0];
rz(-0.5838809) q[0];
x q[1];
rz(2.8526659) q[2];
sx q[2];
rz(-0.75519604) q[2];
sx q[2];
rz(-2.4946405) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7168658) q[1];
sx q[1];
rz(-1.2024554) q[1];
sx q[1];
rz(-1.4348381) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8374568) q[3];
sx q[3];
rz(-1.5392996) q[3];
sx q[3];
rz(-2.568752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62766084) q[2];
sx q[2];
rz(-1.7646503) q[2];
sx q[2];
rz(-0.39109209) q[2];
rz(0.62134653) q[3];
sx q[3];
rz(-0.73453271) q[3];
sx q[3];
rz(-0.77596107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67052996) q[0];
sx q[0];
rz(-3.0369018) q[0];
sx q[0];
rz(-1.712557) q[0];
rz(0.96587005) q[1];
sx q[1];
rz(-1.8873676) q[1];
sx q[1];
rz(-0.78580725) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2656276) q[0];
sx q[0];
rz(-0.64564359) q[0];
sx q[0];
rz(3.0558048) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4072717) q[2];
sx q[2];
rz(-1.87748) q[2];
sx q[2];
rz(1.1994565) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.475226) q[1];
sx q[1];
rz(-2.6738648) q[1];
sx q[1];
rz(-0.18233129) q[1];
x q[2];
rz(2.7034055) q[3];
sx q[3];
rz(-1.3071968) q[3];
sx q[3];
rz(0.98942843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.69661951) q[2];
sx q[2];
rz(-0.61903054) q[2];
sx q[2];
rz(2.6444198) q[2];
rz(-2.2677126) q[3];
sx q[3];
rz(-1.0920352) q[3];
sx q[3];
rz(1.2018275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9484321) q[0];
sx q[0];
rz(-2.8346297) q[0];
sx q[0];
rz(-2.4440785) q[0];
rz(2.7613617) q[1];
sx q[1];
rz(-2.0262599) q[1];
sx q[1];
rz(-0.14022216) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74497667) q[0];
sx q[0];
rz(-0.19752398) q[0];
sx q[0];
rz(-2.7382988) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6096787) q[2];
sx q[2];
rz(-1.1302396) q[2];
sx q[2];
rz(1.3300542) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48988399) q[1];
sx q[1];
rz(-1.9700642) q[1];
sx q[1];
rz(-2.4149766) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41922064) q[3];
sx q[3];
rz(-0.44525075) q[3];
sx q[3];
rz(0.78263301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85072881) q[2];
sx q[2];
rz(-1.2173434) q[2];
sx q[2];
rz(0.49120894) q[2];
rz(-2.9344007) q[3];
sx q[3];
rz(-0.62316337) q[3];
sx q[3];
rz(-1.7306958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
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
rz(2.4405759) q[1];
sx q[1];
rz(-2.0471408) q[1];
sx q[1];
rz(-1.4775803) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.701909) q[0];
sx q[0];
rz(-0.93937568) q[0];
sx q[0];
rz(1.4203493) q[0];
rz(2.9803552) q[2];
sx q[2];
rz(-1.45597) q[2];
sx q[2];
rz(-1.6063362) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3353053) q[1];
sx q[1];
rz(-1.6769969) q[1];
sx q[1];
rz(2.9890636) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83563135) q[3];
sx q[3];
rz(-1.7750677) q[3];
sx q[3];
rz(2.6428619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.574934) q[2];
sx q[2];
rz(-0.40859544) q[2];
sx q[2];
rz(1.5236141) q[2];
rz(2.5426215) q[3];
sx q[3];
rz(-2.6409769) q[3];
sx q[3];
rz(2.9536501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69342518) q[0];
sx q[0];
rz(-0.42891112) q[0];
sx q[0];
rz(1.5018139) q[0];
rz(-0.93694726) q[1];
sx q[1];
rz(-1.0277156) q[1];
sx q[1];
rz(1.017259) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.992998) q[0];
sx q[0];
rz(-2.5879597) q[0];
sx q[0];
rz(2.3675336) q[0];
x q[1];
rz(1.4849365) q[2];
sx q[2];
rz(-1.6756264) q[2];
sx q[2];
rz(1.8293928) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3268633) q[1];
sx q[1];
rz(-1.6039324) q[1];
sx q[1];
rz(-1.3087981) q[1];
rz(-pi) q[2];
rz(-0.42933257) q[3];
sx q[3];
rz(-1.3868012) q[3];
sx q[3];
rz(1.3972549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0246058) q[2];
sx q[2];
rz(-2.459343) q[2];
sx q[2];
rz(1.5218706) q[2];
rz(1.4874602) q[3];
sx q[3];
rz(-1.4922851) q[3];
sx q[3];
rz(1.7447932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7617154) q[0];
sx q[0];
rz(-1.7424255) q[0];
sx q[0];
rz(-2.4095834) q[0];
rz(-1.4749745) q[1];
sx q[1];
rz(-1.2154308) q[1];
sx q[1];
rz(-2.348127) q[1];
rz(-0.99967069) q[2];
sx q[2];
rz(-2.4785505) q[2];
sx q[2];
rz(0.94658755) q[2];
rz(-0.75763221) q[3];
sx q[3];
rz(-0.42790596) q[3];
sx q[3];
rz(-0.71732646) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

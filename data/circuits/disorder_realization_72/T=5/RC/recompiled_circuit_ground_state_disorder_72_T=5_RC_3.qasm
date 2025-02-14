OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6943964) q[0];
sx q[0];
rz(-1.6262654) q[0];
sx q[0];
rz(-0.84994999) q[0];
rz(-1.7605468) q[1];
sx q[1];
rz(-2.2989506) q[1];
sx q[1];
rz(3.091264) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.91431) q[0];
sx q[0];
rz(-1.3903872) q[0];
sx q[0];
rz(1.4854027) q[0];
rz(1.3379132) q[2];
sx q[2];
rz(-0.63265975) q[2];
sx q[2];
rz(-0.92756995) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1427372) q[1];
sx q[1];
rz(-1.5010479) q[1];
sx q[1];
rz(0.19474366) q[1];
rz(-pi) q[2];
rz(1.6292455) q[3];
sx q[3];
rz(-2.2104467) q[3];
sx q[3];
rz(0.15098886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7488659) q[2];
sx q[2];
rz(-1.8957081) q[2];
sx q[2];
rz(-2.5457814) q[2];
rz(-2.3303604) q[3];
sx q[3];
rz(-1.8758643) q[3];
sx q[3];
rz(-2.5442512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8314464) q[0];
sx q[0];
rz(-2.6409288) q[0];
sx q[0];
rz(1.798382) q[0];
rz(1.385618) q[1];
sx q[1];
rz(-1.132553) q[1];
sx q[1];
rz(2.0766506) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6166109) q[0];
sx q[0];
rz(-3.0654902) q[0];
sx q[0];
rz(-2.8625001) q[0];
x q[1];
rz(2.166719) q[2];
sx q[2];
rz(-2.5346014) q[2];
sx q[2];
rz(-2.9526504) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7635325) q[1];
sx q[1];
rz(-3.0010975) q[1];
sx q[1];
rz(-3.1011081) q[1];
rz(-pi) q[2];
rz(1.4707758) q[3];
sx q[3];
rz(-1.6006114) q[3];
sx q[3];
rz(0.884197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.500835) q[2];
sx q[2];
rz(-2.1186327) q[2];
sx q[2];
rz(-2.7309321) q[2];
rz(-1.9985577) q[3];
sx q[3];
rz(-0.7468907) q[3];
sx q[3];
rz(-0.66543287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2086585) q[0];
sx q[0];
rz(-1.7993131) q[0];
sx q[0];
rz(-2.4564273) q[0];
rz(-0.67600983) q[1];
sx q[1];
rz(-1.490255) q[1];
sx q[1];
rz(-0.73838678) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3078416) q[0];
sx q[0];
rz(-2.0594119) q[0];
sx q[0];
rz(1.2060382) q[0];
x q[1];
rz(-1.5398847) q[2];
sx q[2];
rz(-2.3031903) q[2];
sx q[2];
rz(0.67275999) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1967839) q[1];
sx q[1];
rz(-2.9843708) q[1];
sx q[1];
rz(1.9095837) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4285391) q[3];
sx q[3];
rz(-1.0950297) q[3];
sx q[3];
rz(-2.4775164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7076149) q[2];
sx q[2];
rz(-1.6699764) q[2];
sx q[2];
rz(-0.99226704) q[2];
rz(-0.17539242) q[3];
sx q[3];
rz(-0.49687353) q[3];
sx q[3];
rz(-0.3106671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53395143) q[0];
sx q[0];
rz(-1.1247922) q[0];
sx q[0];
rz(2.9006309) q[0];
rz(2.2843649) q[1];
sx q[1];
rz(-2.3576184) q[1];
sx q[1];
rz(1.2711752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64528148) q[0];
sx q[0];
rz(-1.7305614) q[0];
sx q[0];
rz(2.0156572) q[0];
rz(-pi) q[1];
rz(2.056678) q[2];
sx q[2];
rz(-2.0419697) q[2];
sx q[2];
rz(-0.65907329) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.037447151) q[1];
sx q[1];
rz(-1.1318143) q[1];
sx q[1];
rz(-0.46260117) q[1];
rz(-pi) q[2];
rz(1.1587093) q[3];
sx q[3];
rz(-2.2818687) q[3];
sx q[3];
rz(0.76896104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6406389) q[2];
sx q[2];
rz(-0.95879889) q[2];
sx q[2];
rz(2.1232088) q[2];
rz(1.9765249) q[3];
sx q[3];
rz(-0.94912052) q[3];
sx q[3];
rz(-0.44786662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0388221) q[0];
sx q[0];
rz(-0.7907246) q[0];
sx q[0];
rz(-0.08654174) q[0];
rz(-1.4718919) q[1];
sx q[1];
rz(-2.4457928) q[1];
sx q[1];
rz(-2.5873628) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80866586) q[0];
sx q[0];
rz(-0.73360591) q[0];
sx q[0];
rz(-2.8085719) q[0];
x q[1];
rz(1.4325028) q[2];
sx q[2];
rz(-2.2364762) q[2];
sx q[2];
rz(2.4709357) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0010595) q[1];
sx q[1];
rz(-0.78246087) q[1];
sx q[1];
rz(2.4665574) q[1];
rz(0.29227726) q[3];
sx q[3];
rz(-1.9950047) q[3];
sx q[3];
rz(0.19355336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.294813) q[2];
sx q[2];
rz(-0.3598992) q[2];
sx q[2];
rz(-2.2100718) q[2];
rz(-0.31342634) q[3];
sx q[3];
rz(-2.4999908) q[3];
sx q[3];
rz(0.48345598) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1035006) q[0];
sx q[0];
rz(-0.8041389) q[0];
sx q[0];
rz(1.7622129) q[0];
rz(0.74229678) q[1];
sx q[1];
rz(-1.392044) q[1];
sx q[1];
rz(-1.0662063) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3138414) q[0];
sx q[0];
rz(-3.0004639) q[0];
sx q[0];
rz(2.9147097) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0920383) q[2];
sx q[2];
rz(-1.4495499) q[2];
sx q[2];
rz(2.9444314) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51681337) q[1];
sx q[1];
rz(-1.8210026) q[1];
sx q[1];
rz(-2.183319) q[1];
x q[2];
rz(1.6558075) q[3];
sx q[3];
rz(-1.3004616) q[3];
sx q[3];
rz(-0.47209376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1833056) q[2];
sx q[2];
rz(-2.1669407) q[2];
sx q[2];
rz(-0.24678123) q[2];
rz(-2.0731481) q[3];
sx q[3];
rz(-1.9717792) q[3];
sx q[3];
rz(-1.0698414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5088155) q[0];
sx q[0];
rz(-2.1822073) q[0];
sx q[0];
rz(-0.36163133) q[0];
rz(1.5879141) q[1];
sx q[1];
rz(-1.5767153) q[1];
sx q[1];
rz(-1.2036926) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7772352) q[0];
sx q[0];
rz(-2.2493752) q[0];
sx q[0];
rz(-2.4410072) q[0];
x q[1];
rz(1.3839528) q[2];
sx q[2];
rz(-2.0645541) q[2];
sx q[2];
rz(0.46390033) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4154012) q[1];
sx q[1];
rz(-0.63242373) q[1];
sx q[1];
rz(1.9167621) q[1];
rz(-pi) q[2];
rz(2.702868) q[3];
sx q[3];
rz(-1.8774021) q[3];
sx q[3];
rz(-1.0170938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2162073) q[2];
sx q[2];
rz(-1.46571) q[2];
sx q[2];
rz(-1.0125796) q[2];
rz(-2.8818164) q[3];
sx q[3];
rz(-0.24661073) q[3];
sx q[3];
rz(-2.7369505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.986056) q[0];
sx q[0];
rz(-1.2825092) q[0];
sx q[0];
rz(0.78722659) q[0];
rz(2.2975445) q[1];
sx q[1];
rz(-0.35875741) q[1];
sx q[1];
rz(-1.2740096) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0975794) q[0];
sx q[0];
rz(-2.79956) q[0];
sx q[0];
rz(-2.9652389) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4018658) q[2];
sx q[2];
rz(-1.8056895) q[2];
sx q[2];
rz(-1.7559831) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2383409) q[1];
sx q[1];
rz(-2.1834186) q[1];
sx q[1];
rz(2.5765221) q[1];
x q[2];
rz(-2.9926821) q[3];
sx q[3];
rz(-0.40904466) q[3];
sx q[3];
rz(-0.27924505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0446223) q[2];
sx q[2];
rz(-1.5771834) q[2];
sx q[2];
rz(-0.38193646) q[2];
rz(1.1218128) q[3];
sx q[3];
rz(-0.33676454) q[3];
sx q[3];
rz(0.063966123) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065351039) q[0];
sx q[0];
rz(-0.92472804) q[0];
sx q[0];
rz(1.4724154) q[0];
rz(-0.30638254) q[1];
sx q[1];
rz(-2.6173321) q[1];
sx q[1];
rz(-2.870097) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58518325) q[0];
sx q[0];
rz(-2.5228254) q[0];
sx q[0];
rz(0.84689253) q[0];
rz(-1.1708529) q[2];
sx q[2];
rz(-0.78151199) q[2];
sx q[2];
rz(-2.0013546) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2047386) q[1];
sx q[1];
rz(-0.7619155) q[1];
sx q[1];
rz(0.73128766) q[1];
rz(-pi) q[2];
rz(-0.90320194) q[3];
sx q[3];
rz(-2.4861504) q[3];
sx q[3];
rz(0.02070657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86758119) q[2];
sx q[2];
rz(-2.0143955) q[2];
sx q[2];
rz(0.66068823) q[2];
rz(0.88179669) q[3];
sx q[3];
rz(-0.56395689) q[3];
sx q[3];
rz(0.86625117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26447403) q[0];
sx q[0];
rz(-1.5333804) q[0];
sx q[0];
rz(1.8033173) q[0];
rz(1.0461461) q[1];
sx q[1];
rz(-1.0271065) q[1];
sx q[1];
rz(0.17790067) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38707458) q[0];
sx q[0];
rz(-2.1505668) q[0];
sx q[0];
rz(0.80106082) q[0];
x q[1];
rz(2.9650745) q[2];
sx q[2];
rz(-2.1452139) q[2];
sx q[2];
rz(0.18744577) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85965751) q[1];
sx q[1];
rz(-1.9023696) q[1];
sx q[1];
rz(-0.34032199) q[1];
rz(-1.4334861) q[3];
sx q[3];
rz(-1.2009283) q[3];
sx q[3];
rz(0.22320492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.88343128) q[2];
sx q[2];
rz(-2.0475755) q[2];
sx q[2];
rz(0.9359614) q[2];
rz(0.59518138) q[3];
sx q[3];
rz(-1.7116356) q[3];
sx q[3];
rz(2.2619251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4277988) q[0];
sx q[0];
rz(-1.6811163) q[0];
sx q[0];
rz(1.2291193) q[0];
rz(-1.8145369) q[1];
sx q[1];
rz(-1.6358903) q[1];
sx q[1];
rz(-1.9930175) q[1];
rz(3.0279866) q[2];
sx q[2];
rz(-1.7470215) q[2];
sx q[2];
rz(2.5241008) q[2];
rz(1.8775599) q[3];
sx q[3];
rz(-1.8246834) q[3];
sx q[3];
rz(1.9273458) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

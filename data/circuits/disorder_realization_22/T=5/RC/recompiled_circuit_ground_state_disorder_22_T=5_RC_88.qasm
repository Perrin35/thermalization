OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.30224213) q[0];
sx q[0];
rz(-1.2737162) q[0];
sx q[0];
rz(-2.1923375) q[0];
rz(-1.0957837) q[1];
sx q[1];
rz(-0.97691184) q[1];
sx q[1];
rz(-1.3499324) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5744517) q[0];
sx q[0];
rz(-1.0691088) q[0];
sx q[0];
rz(1.4760963) q[0];
rz(-pi) q[1];
x q[1];
rz(2.808936) q[2];
sx q[2];
rz(-1.4321308) q[2];
sx q[2];
rz(-1.1188511) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1848371) q[1];
sx q[1];
rz(-0.3809027) q[1];
sx q[1];
rz(2.95762) q[1];
x q[2];
rz(-0.2738936) q[3];
sx q[3];
rz(-1.1526169) q[3];
sx q[3];
rz(-2.1194715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6076516) q[2];
sx q[2];
rz(-1.093163) q[2];
sx q[2];
rz(-2.9489813) q[2];
rz(-1.8588148) q[3];
sx q[3];
rz(-2.0359813) q[3];
sx q[3];
rz(2.5530596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.823536) q[0];
sx q[0];
rz(-0.41724351) q[0];
sx q[0];
rz(-2.4349924) q[0];
rz(-2.508714) q[1];
sx q[1];
rz(-1.0989847) q[1];
sx q[1];
rz(2.852827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5486886) q[0];
sx q[0];
rz(-0.0026772896) q[0];
sx q[0];
rz(-2.755318) q[0];
rz(-1.7328506) q[2];
sx q[2];
rz(-1.3657492) q[2];
sx q[2];
rz(-2.1265202) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.065246) q[1];
sx q[1];
rz(-1.3586805) q[1];
sx q[1];
rz(-2.8492623) q[1];
rz(1.992728) q[3];
sx q[3];
rz(-2.2852118) q[3];
sx q[3];
rz(1.2846636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67549813) q[2];
sx q[2];
rz(-2.0530632) q[2];
sx q[2];
rz(-1.4906073) q[2];
rz(2.364482) q[3];
sx q[3];
rz(-1.6831574) q[3];
sx q[3];
rz(-1.3597663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8260088) q[0];
sx q[0];
rz(-1.6827787) q[0];
sx q[0];
rz(-2.7050731) q[0];
rz(-0.79477683) q[1];
sx q[1];
rz(-1.3705285) q[1];
sx q[1];
rz(-0.93903843) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3600814) q[0];
sx q[0];
rz(-2.0716801) q[0];
sx q[0];
rz(1.341218) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76656966) q[2];
sx q[2];
rz(-1.1423282) q[2];
sx q[2];
rz(-0.049190532) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4758376) q[1];
sx q[1];
rz(-1.8867954) q[1];
sx q[1];
rz(0.60029392) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22041286) q[3];
sx q[3];
rz(-1.0272946) q[3];
sx q[3];
rz(1.2459178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44094917) q[2];
sx q[2];
rz(-2.078853) q[2];
sx q[2];
rz(2.5062594) q[2];
rz(-2.3144531) q[3];
sx q[3];
rz(-2.6878036) q[3];
sx q[3];
rz(1.872983) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24469911) q[0];
sx q[0];
rz(-3.0785955) q[0];
sx q[0];
rz(2.6196106) q[0];
rz(2.7105647) q[1];
sx q[1];
rz(-1.5343752) q[1];
sx q[1];
rz(0.65188754) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4341457) q[0];
sx q[0];
rz(-2.0743255) q[0];
sx q[0];
rz(1.8412526) q[0];
x q[1];
rz(2.6471557) q[2];
sx q[2];
rz(-0.90915138) q[2];
sx q[2];
rz(-2.5376157) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3451772) q[1];
sx q[1];
rz(-2.762156) q[1];
sx q[1];
rz(-1.6077843) q[1];
rz(-2.4383955) q[3];
sx q[3];
rz(-2.0542938) q[3];
sx q[3];
rz(-1.0994224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42644694) q[2];
sx q[2];
rz(-1.8603674) q[2];
sx q[2];
rz(-2.35671) q[2];
rz(2.6134885) q[3];
sx q[3];
rz(-1.8485066) q[3];
sx q[3];
rz(1.2877134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89609471) q[0];
sx q[0];
rz(-2.446785) q[0];
sx q[0];
rz(0.78952638) q[0];
rz(0.49742571) q[1];
sx q[1];
rz(-1.0405633) q[1];
sx q[1];
rz(1.8400037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3606982) q[0];
sx q[0];
rz(-2.5142365) q[0];
sx q[0];
rz(2.832546) q[0];
rz(-1.5100432) q[2];
sx q[2];
rz(-1.9092602) q[2];
sx q[2];
rz(2.0450704) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19986049) q[1];
sx q[1];
rz(-2.2944415) q[1];
sx q[1];
rz(-1.3398712) q[1];
x q[2];
rz(-1.5118408) q[3];
sx q[3];
rz(-0.36789933) q[3];
sx q[3];
rz(0.93705642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54399458) q[2];
sx q[2];
rz(-1.602729) q[2];
sx q[2];
rz(2.8653115) q[2];
rz(2.1918519) q[3];
sx q[3];
rz(-2.8730928) q[3];
sx q[3];
rz(-2.5883519) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5211869) q[0];
sx q[0];
rz(-0.036245417) q[0];
sx q[0];
rz(2.1976443) q[0];
rz(1.2983407) q[1];
sx q[1];
rz(-2.0746168) q[1];
sx q[1];
rz(2.3640769) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80332887) q[0];
sx q[0];
rz(-1.216812) q[0];
sx q[0];
rz(2.5577118) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3087511) q[2];
sx q[2];
rz(-0.85393674) q[2];
sx q[2];
rz(-1.0345296) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7878163) q[1];
sx q[1];
rz(-2.750038) q[1];
sx q[1];
rz(-0.33772525) q[1];
rz(-pi) q[2];
rz(-3.0367682) q[3];
sx q[3];
rz(-2.8358805) q[3];
sx q[3];
rz(2.0436398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62766084) q[2];
sx q[2];
rz(-1.3769423) q[2];
sx q[2];
rz(-2.7505006) q[2];
rz(0.62134653) q[3];
sx q[3];
rz(-0.73453271) q[3];
sx q[3];
rz(-0.77596107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4710627) q[0];
sx q[0];
rz(-0.10469086) q[0];
sx q[0];
rz(1.4290357) q[0];
rz(0.96587005) q[1];
sx q[1];
rz(-1.8873676) q[1];
sx q[1];
rz(2.3557854) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98322372) q[0];
sx q[0];
rz(-0.92792032) q[0];
sx q[0];
rz(-1.5063365) q[0];
rz(-2.4072717) q[2];
sx q[2];
rz(-1.2641126) q[2];
sx q[2];
rz(1.1994565) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.25868402) q[1];
sx q[1];
rz(-1.488954) q[1];
sx q[1];
rz(2.6805582) q[1];
rz(2.7034055) q[3];
sx q[3];
rz(-1.8343958) q[3];
sx q[3];
rz(-0.98942843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4449731) q[2];
sx q[2];
rz(-0.61903054) q[2];
sx q[2];
rz(-0.49717286) q[2];
rz(-0.87388006) q[3];
sx q[3];
rz(-2.0495575) q[3];
sx q[3];
rz(1.2018275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9484321) q[0];
sx q[0];
rz(-0.30696294) q[0];
sx q[0];
rz(0.69751414) q[0];
rz(-0.3802309) q[1];
sx q[1];
rz(-2.0262599) q[1];
sx q[1];
rz(-0.14022216) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8070458) q[0];
sx q[0];
rz(-1.3893034) q[0];
sx q[0];
rz(-1.6491778) q[0];
rz(1.6096787) q[2];
sx q[2];
rz(-1.1302396) q[2];
sx q[2];
rz(-1.3300542) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74822357) q[1];
sx q[1];
rz(-2.2294982) q[1];
sx q[1];
rz(-2.084712) q[1];
rz(-0.41922064) q[3];
sx q[3];
rz(-2.6963419) q[3];
sx q[3];
rz(2.3589596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.85072881) q[2];
sx q[2];
rz(-1.2173434) q[2];
sx q[2];
rz(-0.49120894) q[2];
rz(2.9344007) q[3];
sx q[3];
rz(-2.5184293) q[3];
sx q[3];
rz(-1.7306958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17289138) q[0];
sx q[0];
rz(-1.7270813) q[0];
sx q[0];
rz(-2.9911995) q[0];
rz(2.4405759) q[1];
sx q[1];
rz(-1.0944518) q[1];
sx q[1];
rz(1.4775803) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1883066) q[0];
sx q[0];
rz(-2.4948848) q[0];
sx q[0];
rz(0.20215277) q[0];
x q[1];
rz(1.4544746) q[2];
sx q[2];
rz(-1.7309628) q[2];
sx q[2];
rz(3.124686) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3353053) q[1];
sx q[1];
rz(-1.6769969) q[1];
sx q[1];
rz(-0.152529) q[1];
rz(2.3059613) q[3];
sx q[3];
rz(-1.366525) q[3];
sx q[3];
rz(-0.49873078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56665862) q[2];
sx q[2];
rz(-2.7329972) q[2];
sx q[2];
rz(-1.6179786) q[2];
rz(-0.59897113) q[3];
sx q[3];
rz(-2.6409769) q[3];
sx q[3];
rz(-0.18794255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69342518) q[0];
sx q[0];
rz(-2.7126815) q[0];
sx q[0];
rz(1.6397788) q[0];
rz(-2.2046454) q[1];
sx q[1];
rz(-2.1138771) q[1];
sx q[1];
rz(-2.1243336) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27148025) q[0];
sx q[0];
rz(-1.9471629) q[0];
sx q[0];
rz(0.41618213) q[0];
rz(-pi) q[1];
rz(-1.4849365) q[2];
sx q[2];
rz(-1.4659662) q[2];
sx q[2];
rz(1.8293928) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76495095) q[1];
sx q[1];
rz(-1.8326474) q[1];
sx q[1];
rz(-3.1072867) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42041619) q[3];
sx q[3];
rz(-2.6767593) q[3];
sx q[3];
rz(0.55373389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1169869) q[2];
sx q[2];
rz(-2.459343) q[2];
sx q[2];
rz(-1.619722) q[2];
rz(-1.4874602) q[3];
sx q[3];
rz(-1.4922851) q[3];
sx q[3];
rz(1.3967995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7617154) q[0];
sx q[0];
rz(-1.7424255) q[0];
sx q[0];
rz(-2.4095834) q[0];
rz(1.4749745) q[1];
sx q[1];
rz(-1.9261618) q[1];
sx q[1];
rz(0.7934657) q[1];
rz(0.99967069) q[2];
sx q[2];
rz(-0.66304211) q[2];
sx q[2];
rz(-2.1950051) q[2];
rz(-0.31994896) q[3];
sx q[3];
rz(-1.2816164) q[3];
sx q[3];
rz(1.5641227) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

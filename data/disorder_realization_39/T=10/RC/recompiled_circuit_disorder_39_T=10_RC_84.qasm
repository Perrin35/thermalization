OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0948148) q[0];
sx q[0];
rz(-2.0733641) q[0];
sx q[0];
rz(0.46407035) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(-1.2844515) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099079236) q[0];
sx q[0];
rz(-2.6455542) q[0];
sx q[0];
rz(0.92143671) q[0];
rz(2.651865) q[2];
sx q[2];
rz(-1.5999319) q[2];
sx q[2];
rz(0.65087989) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1844993) q[1];
sx q[1];
rz(-0.59809369) q[1];
sx q[1];
rz(-0.68006541) q[1];
rz(-1.554261) q[3];
sx q[3];
rz(-1.4201122) q[3];
sx q[3];
rz(-2.2076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7444732) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(-0.31952566) q[2];
rz(0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(-2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7330866) q[0];
sx q[0];
rz(-1.5698313) q[0];
sx q[0];
rz(2.615036) q[0];
rz(0.56354848) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(0.79663509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8684727) q[0];
sx q[0];
rz(-2.1935049) q[0];
sx q[0];
rz(2.7362583) q[0];
rz(-pi) q[1];
rz(2.8393306) q[2];
sx q[2];
rz(-0.59603359) q[2];
sx q[2];
rz(2.9994534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4371725) q[1];
sx q[1];
rz(-0.96357268) q[1];
sx q[1];
rz(1.5864658) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1809686) q[3];
sx q[3];
rz(-1.989813) q[3];
sx q[3];
rz(0.39715365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18156302) q[2];
sx q[2];
rz(-1.3788362) q[2];
sx q[2];
rz(0.78655085) q[2];
rz(-0.49318796) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(-2.6942159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.2717473) q[0];
sx q[0];
rz(-0.87690502) q[0];
sx q[0];
rz(1.440381) q[0];
rz(-0.72021833) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(2.4386141) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36613208) q[0];
sx q[0];
rz(-1.4538987) q[0];
sx q[0];
rz(1.2890105) q[0];
rz(-pi) q[1];
rz(-2.5419652) q[2];
sx q[2];
rz(-2.5587974) q[2];
sx q[2];
rz(-1.3404913) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2726047) q[1];
sx q[1];
rz(-1.4613978) q[1];
sx q[1];
rz(1.617856) q[1];
rz(-pi) q[2];
rz(-0.94743518) q[3];
sx q[3];
rz(-1.4294335) q[3];
sx q[3];
rz(-0.75368222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.26677033) q[2];
sx q[2];
rz(-1.4845994) q[2];
sx q[2];
rz(-0.91153574) q[2];
rz(-2.2276145) q[3];
sx q[3];
rz(-1.4303047) q[3];
sx q[3];
rz(-2.3142557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5754958) q[0];
sx q[0];
rz(-2.5056582) q[0];
sx q[0];
rz(2.3655868) q[0];
rz(-1.874118) q[1];
sx q[1];
rz(-2.0327366) q[1];
sx q[1];
rz(0.56328303) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8877836) q[0];
sx q[0];
rz(-2.1953708) q[0];
sx q[0];
rz(0.33872351) q[0];
x q[1];
rz(-0.69022501) q[2];
sx q[2];
rz(-1.8340655) q[2];
sx q[2];
rz(0.0036247591) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.54484425) q[1];
sx q[1];
rz(-1.8425643) q[1];
sx q[1];
rz(1.5640869) q[1];
rz(1.2761649) q[3];
sx q[3];
rz(-0.8751103) q[3];
sx q[3];
rz(0.61616117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48878601) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(-0.164786) q[2];
rz(-0.22848836) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(2.1951108) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27424681) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(-2.2891323) q[0];
rz(-2.7903941) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(-1.16211) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0162109) q[0];
sx q[0];
rz(-2.1124766) q[0];
sx q[0];
rz(2.0340232) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7439793) q[2];
sx q[2];
rz(-2.0587454) q[2];
sx q[2];
rz(0.60148009) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4967242) q[1];
sx q[1];
rz(-0.86360303) q[1];
sx q[1];
rz(1.0134539) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.01708548) q[3];
sx q[3];
rz(-2.4199329) q[3];
sx q[3];
rz(0.34860308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6953485) q[2];
sx q[2];
rz(-0.39351311) q[2];
sx q[2];
rz(-1.8943141) q[2];
rz(3.128483) q[3];
sx q[3];
rz(-1.40991) q[3];
sx q[3];
rz(2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2728249) q[0];
sx q[0];
rz(-1.8936963) q[0];
sx q[0];
rz(0.75063467) q[0];
rz(-1.1095095) q[1];
sx q[1];
rz(-2.7710319) q[1];
sx q[1];
rz(1.3060588) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40697843) q[0];
sx q[0];
rz(-0.66440551) q[0];
sx q[0];
rz(-2.2218496) q[0];
x q[1];
rz(-2.7517031) q[2];
sx q[2];
rz(-0.49553686) q[2];
sx q[2];
rz(1.3930266) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27481025) q[1];
sx q[1];
rz(-0.85046235) q[1];
sx q[1];
rz(1.5606828) q[1];
rz(-pi) q[2];
rz(-1.7197051) q[3];
sx q[3];
rz(-0.9871261) q[3];
sx q[3];
rz(-2.7910809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3874454) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(2.5409017) q[2];
rz(2.1389652) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(0.35282648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3951185) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(-2.0948998) q[0];
rz(-1.5294317) q[1];
sx q[1];
rz(-1.4271095) q[1];
sx q[1];
rz(-2.7244862) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1244303) q[0];
sx q[0];
rz(-0.22148795) q[0];
sx q[0];
rz(-1.4593967) q[0];
x q[1];
rz(1.8646556) q[2];
sx q[2];
rz(-1.8655348) q[2];
sx q[2];
rz(-0.35640946) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4579791) q[1];
sx q[1];
rz(-0.90404592) q[1];
sx q[1];
rz(0.97138202) q[1];
x q[2];
rz(-1.8700637) q[3];
sx q[3];
rz(-2.2021658) q[3];
sx q[3];
rz(-0.33559468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1356915) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(1.7741514) q[2];
rz(-0.55316365) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.816514) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(1.9518071) q[0];
rz(-1.7143543) q[1];
sx q[1];
rz(-2.863527) q[1];
sx q[1];
rz(0.11238012) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.534879) q[0];
sx q[0];
rz(-1.0977854) q[0];
sx q[0];
rz(-2.9472449) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24765315) q[2];
sx q[2];
rz(-0.34341771) q[2];
sx q[2];
rz(2.0695956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1689414) q[1];
sx q[1];
rz(-2.2691138) q[1];
sx q[1];
rz(2.6886743) q[1];
rz(-1.8295248) q[3];
sx q[3];
rz(-1.6885763) q[3];
sx q[3];
rz(3.0974914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0743951) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(-1.7929662) q[2];
rz(-1.9366692) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7444721) q[0];
sx q[0];
rz(-1.4646126) q[0];
sx q[0];
rz(-0.034974139) q[0];
rz(-2.2947625) q[1];
sx q[1];
rz(-1.433082) q[1];
sx q[1];
rz(-0.91167489) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.864894) q[0];
sx q[0];
rz(-1.8407341) q[0];
sx q[0];
rz(0.2322659) q[0];
rz(-pi) q[1];
rz(1.6976835) q[2];
sx q[2];
rz(-2.1107026) q[2];
sx q[2];
rz(-1.1586231) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4166491) q[1];
sx q[1];
rz(-2.0282201) q[1];
sx q[1];
rz(1.7979421) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4089912) q[3];
sx q[3];
rz(-2.4676975) q[3];
sx q[3];
rz(3.1042838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.70790616) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(2.8923477) q[2];
rz(2.3748659) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(-0.39961091) q[3];
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
rz(-pi/2) q[3];
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
rz(-1.7837759) q[0];
sx q[0];
rz(-1.0483085) q[0];
sx q[0];
rz(-0.37208474) q[0];
rz(-2.5601939) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(1.7262329) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2863662) q[0];
sx q[0];
rz(-2.8088673) q[0];
sx q[0];
rz(2.4871728) q[0];
rz(0.62016983) q[2];
sx q[2];
rz(-1.1405186) q[2];
sx q[2];
rz(-2.534453) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.49299875) q[1];
sx q[1];
rz(-0.10222888) q[1];
sx q[1];
rz(-0.74536721) q[1];
rz(1.0269208) q[3];
sx q[3];
rz(-1.9034991) q[3];
sx q[3];
rz(2.016891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8156478) q[2];
sx q[2];
rz(-3.0013913) q[2];
sx q[2];
rz(-0.20467219) q[2];
rz(-1.4137911) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(1.0958825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50080147) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(1.5851371) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(1.1164222) q[2];
sx q[2];
rz(-1.6797671) q[2];
sx q[2];
rz(-1.4327008) q[2];
rz(2.1429569) q[3];
sx q[3];
rz(-1.5012267) q[3];
sx q[3];
rz(-2.5517626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
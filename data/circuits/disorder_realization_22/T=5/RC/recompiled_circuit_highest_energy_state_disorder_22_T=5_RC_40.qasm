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
rz(2.4616315) q[0];
sx q[0];
rz(7.1890561) q[0];
sx q[0];
rz(11.472975) q[0];
rz(-1.6361341) q[1];
sx q[1];
rz(-1.1727762) q[1];
sx q[1];
rz(-0.42458951) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.070804) q[0];
sx q[0];
rz(-1.53264) q[0];
sx q[0];
rz(-1.7724994) q[0];
rz(1.8867887) q[2];
sx q[2];
rz(-2.9796763) q[2];
sx q[2];
rz(-0.94053167) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.54271419) q[1];
sx q[1];
rz(-1.6495566) q[1];
sx q[1];
rz(1.8171189) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2810117) q[3];
sx q[3];
rz(-1.672555) q[3];
sx q[3];
rz(1.4836131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37497416) q[2];
sx q[2];
rz(-1.582229) q[2];
sx q[2];
rz(1.1589104) q[2];
rz(2.9186987) q[3];
sx q[3];
rz(-1.4054207) q[3];
sx q[3];
rz(2.8515653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3463773) q[0];
sx q[0];
rz(-1.5257436) q[0];
sx q[0];
rz(1.079153) q[0];
rz(-1.4214628) q[1];
sx q[1];
rz(-2.454897) q[1];
sx q[1];
rz(3.0770643) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33693211) q[0];
sx q[0];
rz(-1.6580026) q[0];
sx q[0];
rz(0.052636458) q[0];
x q[1];
rz(3.0874443) q[2];
sx q[2];
rz(-0.47685078) q[2];
sx q[2];
rz(2.3091174) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77212438) q[1];
sx q[1];
rz(-1.3288979) q[1];
sx q[1];
rz(0.78472225) q[1];
rz(-pi) q[2];
rz(1.7713624) q[3];
sx q[3];
rz(-0.75955694) q[3];
sx q[3];
rz(2.650039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1530389) q[2];
sx q[2];
rz(-1.3554074) q[2];
sx q[2];
rz(2.0241731) q[2];
rz(0.85876632) q[3];
sx q[3];
rz(-2.358181) q[3];
sx q[3];
rz(-2.7045238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3808463) q[0];
sx q[0];
rz(-2.8570211) q[0];
sx q[0];
rz(-2.9685156) q[0];
rz(-2.1848047) q[1];
sx q[1];
rz(-2.1134977) q[1];
sx q[1];
rz(-1.0622567) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5231402) q[0];
sx q[0];
rz(-1.5644531) q[0];
sx q[0];
rz(0.84756084) q[0];
rz(-pi) q[1];
rz(-0.99750869) q[2];
sx q[2];
rz(-1.0110223) q[2];
sx q[2];
rz(-2.5791175) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.838321) q[1];
sx q[1];
rz(-1.0365937) q[1];
sx q[1];
rz(-0.81373416) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2532387) q[3];
sx q[3];
rz(-2.0815947) q[3];
sx q[3];
rz(2.3470283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.7303077) q[2];
sx q[2];
rz(-1.9399119) q[2];
sx q[2];
rz(2.3968706) q[2];
rz(-2.5005285) q[3];
sx q[3];
rz(-1.3499667) q[3];
sx q[3];
rz(-0.49066576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6915801) q[0];
sx q[0];
rz(-0.55679655) q[0];
sx q[0];
rz(2.2963754) q[0];
rz(-0.40329626) q[1];
sx q[1];
rz(-1.6019628) q[1];
sx q[1];
rz(2.8135615) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.384399) q[0];
sx q[0];
rz(-1.4310657) q[0];
sx q[0];
rz(1.7534885) q[0];
rz(-pi) q[1];
rz(0.26706605) q[2];
sx q[2];
rz(-2.895148) q[2];
sx q[2];
rz(1.3738969) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1650585) q[1];
sx q[1];
rz(-1.6886025) q[1];
sx q[1];
rz(0.56748135) q[1];
x q[2];
rz(0.27356903) q[3];
sx q[3];
rz(-2.0988587) q[3];
sx q[3];
rz(-2.7298965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8959117) q[2];
sx q[2];
rz(-0.74378219) q[2];
sx q[2];
rz(-0.15920676) q[2];
rz(-3.1253452) q[3];
sx q[3];
rz(-1.0280321) q[3];
sx q[3];
rz(-2.4102559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74723393) q[0];
sx q[0];
rz(-1.9488229) q[0];
sx q[0];
rz(-1.0900981) q[0];
rz(-0.12995003) q[1];
sx q[1];
rz(-2.3268685) q[1];
sx q[1];
rz(1.5257588) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1795883) q[0];
sx q[0];
rz(-1.1561014) q[0];
sx q[0];
rz(0.59831043) q[0];
rz(-1.4065625) q[2];
sx q[2];
rz(-1.4964074) q[2];
sx q[2];
rz(1.0094355) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8860717) q[1];
sx q[1];
rz(-0.93734765) q[1];
sx q[1];
rz(1.5367763) q[1];
rz(1.3899743) q[3];
sx q[3];
rz(-1.5360334) q[3];
sx q[3];
rz(-1.8145479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3913564) q[2];
sx q[2];
rz(-0.82166925) q[2];
sx q[2];
rz(0.4772805) q[2];
rz(-0.20398772) q[3];
sx q[3];
rz(-2.9450649) q[3];
sx q[3];
rz(-2.5175214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1962965) q[0];
sx q[0];
rz(-2.3260703) q[0];
sx q[0];
rz(0.77769172) q[0];
rz(2.5241191) q[1];
sx q[1];
rz(-1.4910881) q[1];
sx q[1];
rz(2.2106574) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5054748) q[0];
sx q[0];
rz(-1.9808928) q[0];
sx q[0];
rz(-1.0763041) q[0];
x q[1];
rz(2.7226376) q[2];
sx q[2];
rz(-2.0134996) q[2];
sx q[2];
rz(2.1328164) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9677862) q[1];
sx q[1];
rz(-1.5266982) q[1];
sx q[1];
rz(1.2323061) q[1];
x q[2];
rz(2.10465) q[3];
sx q[3];
rz(-2.4483238) q[3];
sx q[3];
rz(1.1806837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.65333873) q[2];
sx q[2];
rz(-1.808017) q[2];
sx q[2];
rz(-0.2674357) q[2];
rz(-2.2087162) q[3];
sx q[3];
rz(-1.2837912) q[3];
sx q[3];
rz(-0.4944087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.47144181) q[0];
sx q[0];
rz(-1.143456) q[0];
sx q[0];
rz(2.4757521) q[0];
rz(-0.2039856) q[1];
sx q[1];
rz(-2.1603656) q[1];
sx q[1];
rz(1.1873672) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46118051) q[0];
sx q[0];
rz(-2.0392766) q[0];
sx q[0];
rz(1.4673244) q[0];
rz(-pi) q[1];
rz(0.86249528) q[2];
sx q[2];
rz(-1.1581026) q[2];
sx q[2];
rz(-2.9484381) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.59661181) q[1];
sx q[1];
rz(-0.89681936) q[1];
sx q[1];
rz(-1.2294212) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3755619) q[3];
sx q[3];
rz(-0.25849202) q[3];
sx q[3];
rz(3.0743161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54619706) q[2];
sx q[2];
rz(-0.51837817) q[2];
sx q[2];
rz(-0.67016822) q[2];
rz(-2.8403122) q[3];
sx q[3];
rz(-1.8471085) q[3];
sx q[3];
rz(-0.41845751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8381074) q[0];
sx q[0];
rz(-1.0024339) q[0];
sx q[0];
rz(2.0759034) q[0];
rz(2.4844555) q[1];
sx q[1];
rz(-0.70825759) q[1];
sx q[1];
rz(1.7117975) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2524714) q[0];
sx q[0];
rz(-1.9074797) q[0];
sx q[0];
rz(-0.97674979) q[0];
x q[1];
rz(1.7523319) q[2];
sx q[2];
rz(-1.2444656) q[2];
sx q[2];
rz(-1.2943314) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.011466) q[1];
sx q[1];
rz(-1.3523755) q[1];
sx q[1];
rz(2.7232724) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4265238) q[3];
sx q[3];
rz(-1.383009) q[3];
sx q[3];
rz(-0.71545631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.26057217) q[2];
sx q[2];
rz(-1.2994956) q[2];
sx q[2];
rz(-1.0183498) q[2];
rz(1.5923422) q[3];
sx q[3];
rz(-1.5038303) q[3];
sx q[3];
rz(0.75756592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26579478) q[0];
sx q[0];
rz(-2.6261411) q[0];
sx q[0];
rz(2.1739668) q[0];
rz(-0.34010092) q[1];
sx q[1];
rz(-2.3022771) q[1];
sx q[1];
rz(-1.0338773) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63003016) q[0];
sx q[0];
rz(-2.1242737) q[0];
sx q[0];
rz(-2.9787977) q[0];
x q[1];
rz(-2.5754355) q[2];
sx q[2];
rz(-0.70933178) q[2];
sx q[2];
rz(2.4475803) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.441923) q[1];
sx q[1];
rz(-1.1811387) q[1];
sx q[1];
rz(2.3226435) q[1];
rz(-pi) q[2];
rz(-0.57698099) q[3];
sx q[3];
rz(-2.8762102) q[3];
sx q[3];
rz(-0.25091618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.10244441) q[2];
sx q[2];
rz(-1.4463964) q[2];
sx q[2];
rz(-3.1026802) q[2];
rz(-0.14032042) q[3];
sx q[3];
rz(-3.0571627) q[3];
sx q[3];
rz(-1.3154359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4843531) q[0];
sx q[0];
rz(-2.1359279) q[0];
sx q[0];
rz(-1.8835618) q[0];
rz(-2.925442) q[1];
sx q[1];
rz(-2.436147) q[1];
sx q[1];
rz(0.36453882) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70750694) q[0];
sx q[0];
rz(-2.0179406) q[0];
sx q[0];
rz(2.8672877) q[0];
rz(-pi) q[1];
rz(-2.0272354) q[2];
sx q[2];
rz(-2.3203712) q[2];
sx q[2];
rz(-0.26320339) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3033402) q[1];
sx q[1];
rz(-1.5152182) q[1];
sx q[1];
rz(-1.5838632) q[1];
rz(-pi) q[2];
x q[2];
rz(2.632439) q[3];
sx q[3];
rz(-3.0753071) q[3];
sx q[3];
rz(1.2627511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8517427) q[2];
sx q[2];
rz(-1.936329) q[2];
sx q[2];
rz(-2.5002948) q[2];
rz(-1.6614206) q[3];
sx q[3];
rz(-1.2484173) q[3];
sx q[3];
rz(-0.67354584) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73364532) q[0];
sx q[0];
rz(-0.49325627) q[0];
sx q[0];
rz(2.7706964) q[0];
rz(-0.98450487) q[1];
sx q[1];
rz(-1.4823109) q[1];
sx q[1];
rz(-1.0479814) q[1];
rz(-3.0251236) q[2];
sx q[2];
rz(-1.7670198) q[2];
sx q[2];
rz(2.7311538) q[2];
rz(1.0897286) q[3];
sx q[3];
rz(-1.0952598) q[3];
sx q[3];
rz(2.918603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

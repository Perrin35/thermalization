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
rz(2.4438357) q[0];
sx q[0];
rz(-1.9480167) q[0];
sx q[0];
rz(0.20456631) q[0];
rz(2.3808631) q[1];
sx q[1];
rz(-1.7525571) q[1];
sx q[1];
rz(-1.3995481) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22628173) q[0];
sx q[0];
rz(-1.6618414) q[0];
sx q[0];
rz(-1.7928726) q[0];
rz(1.1209773) q[2];
sx q[2];
rz(-1.9391141) q[2];
sx q[2];
rz(-2.8415547) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7508937) q[1];
sx q[1];
rz(-1.1502153) q[1];
sx q[1];
rz(-1.6665568) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4321626) q[3];
sx q[3];
rz(-1.3660079) q[3];
sx q[3];
rz(-0.11573175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7064887) q[2];
sx q[2];
rz(-3.1092643) q[2];
sx q[2];
rz(-0.48933634) q[2];
rz(-2.1183744) q[3];
sx q[3];
rz(-0.018298572) q[3];
sx q[3];
rz(-1.1255012) q[3];
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
rz(3.0959051) q[0];
sx q[0];
rz(-2.4850595) q[0];
sx q[0];
rz(2.3387961) q[0];
rz(-0.071391694) q[1];
sx q[1];
rz(-0.26669058) q[1];
sx q[1];
rz(3.0838222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18081576) q[0];
sx q[0];
rz(-2.1102) q[0];
sx q[0];
rz(2.8404854) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4647015) q[2];
sx q[2];
rz(-1.8558143) q[2];
sx q[2];
rz(0.34044701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.80177414) q[1];
sx q[1];
rz(-1.3388102) q[1];
sx q[1];
rz(1.07527) q[1];
x q[2];
rz(0.85912786) q[3];
sx q[3];
rz(-1.6831846) q[3];
sx q[3];
rz(-2.4939362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1959261) q[2];
sx q[2];
rz(-1.095093) q[2];
sx q[2];
rz(1.8460974) q[2];
rz(2.1748491) q[3];
sx q[3];
rz(-0.77015489) q[3];
sx q[3];
rz(-0.75743341) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1118689) q[0];
sx q[0];
rz(-1.7787378) q[0];
sx q[0];
rz(-1.7080074) q[0];
rz(3.0729821) q[1];
sx q[1];
rz(-1.5674633) q[1];
sx q[1];
rz(2.5624018) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9703909) q[0];
sx q[0];
rz(-0.10848898) q[0];
sx q[0];
rz(-2.1347339) q[0];
x q[1];
rz(-1.5689108) q[2];
sx q[2];
rz(-1.0950739) q[2];
sx q[2];
rz(1.9370417) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80651307) q[1];
sx q[1];
rz(-1.6066222) q[1];
sx q[1];
rz(2.9153009) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81285254) q[3];
sx q[3];
rz(-0.62004706) q[3];
sx q[3];
rz(0.64921415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.27089831) q[2];
sx q[2];
rz(-1.2535932) q[2];
sx q[2];
rz(-2.9581621) q[2];
rz(-2.3373248) q[3];
sx q[3];
rz(-0.99884123) q[3];
sx q[3];
rz(-0.53406322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.9428228) q[0];
sx q[0];
rz(-0.11480055) q[0];
sx q[0];
rz(2.542069) q[0];
rz(2.7291258) q[1];
sx q[1];
rz(-3.1209374) q[1];
sx q[1];
rz(-1.0106769) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.090792716) q[0];
sx q[0];
rz(-0.37380344) q[0];
sx q[0];
rz(2.7888377) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0907107) q[2];
sx q[2];
rz(-2.2023099) q[2];
sx q[2];
rz(0.33755195) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0914826) q[1];
sx q[1];
rz(-0.98467876) q[1];
sx q[1];
rz(0.32663235) q[1];
rz(-pi) q[2];
rz(1.4323727) q[3];
sx q[3];
rz(-0.81973053) q[3];
sx q[3];
rz(2.8734796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7708873) q[2];
sx q[2];
rz(-2.7949896) q[2];
sx q[2];
rz(2.8240805) q[2];
rz(-0.62234771) q[3];
sx q[3];
rz(-2.205866) q[3];
sx q[3];
rz(-0.58421016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7790826) q[0];
sx q[0];
rz(-2.2267987) q[0];
sx q[0];
rz(-1.0269748) q[0];
rz(0.55038553) q[1];
sx q[1];
rz(-3.0774979) q[1];
sx q[1];
rz(-1.9245573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37792045) q[0];
sx q[0];
rz(-0.91912133) q[0];
sx q[0];
rz(2.4125189) q[0];
x q[1];
rz(-1.6486859) q[2];
sx q[2];
rz(-1.0337892) q[2];
sx q[2];
rz(-2.3324049) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0266307) q[1];
sx q[1];
rz(-1.3991881) q[1];
sx q[1];
rz(-2.2725355) q[1];
rz(-pi) q[2];
rz(-1.9868136) q[3];
sx q[3];
rz(-2.3868594) q[3];
sx q[3];
rz(-0.55547914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43651849) q[2];
sx q[2];
rz(-1.7784092) q[2];
sx q[2];
rz(2.3684033) q[2];
rz(-0.1117205) q[3];
sx q[3];
rz(-1.3216647) q[3];
sx q[3];
rz(2.2022061) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6620827) q[0];
sx q[0];
rz(-0.22308068) q[0];
sx q[0];
rz(-0.42698419) q[0];
rz(-0.93049479) q[1];
sx q[1];
rz(-3.1246779) q[1];
sx q[1];
rz(-0.46447909) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1982959) q[0];
sx q[0];
rz(-1.7664599) q[0];
sx q[0];
rz(0.29384675) q[0];
rz(-2.1709178) q[2];
sx q[2];
rz(-1.8541186) q[2];
sx q[2];
rz(-1.4221869) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.61349866) q[1];
sx q[1];
rz(-1.534675) q[1];
sx q[1];
rz(0.36904676) q[1];
rz(-0.031207009) q[3];
sx q[3];
rz(-1.1466807) q[3];
sx q[3];
rz(-1.1348789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53133196) q[2];
sx q[2];
rz(-1.8208296) q[2];
sx q[2];
rz(-0.28826928) q[2];
rz(-1.0432976) q[3];
sx q[3];
rz(-2.5259924) q[3];
sx q[3];
rz(2.402795) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6601324) q[0];
sx q[0];
rz(-1.8539424) q[0];
sx q[0];
rz(2.3840391) q[0];
rz(3.0689012) q[1];
sx q[1];
rz(-3.1158267) q[1];
sx q[1];
rz(-0.048197897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2295727) q[0];
sx q[0];
rz(-2.1549112) q[0];
sx q[0];
rz(-3.106227) q[0];
rz(2.9075895) q[2];
sx q[2];
rz(-2.1562139) q[2];
sx q[2];
rz(2.8319179) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0695677) q[1];
sx q[1];
rz(-1.1434907) q[1];
sx q[1];
rz(0.91520379) q[1];
x q[2];
rz(-1.4675619) q[3];
sx q[3];
rz(-2.6695731) q[3];
sx q[3];
rz(-1.9341521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.66182071) q[2];
sx q[2];
rz(-1.4996108) q[2];
sx q[2];
rz(-3.0721967) q[2];
rz(1.5540468) q[3];
sx q[3];
rz(-0.78726751) q[3];
sx q[3];
rz(-0.21849304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7540392) q[0];
sx q[0];
rz(-2.0856922) q[0];
sx q[0];
rz(-1.7283424) q[0];
rz(-0.82855254) q[1];
sx q[1];
rz(-3.0998402) q[1];
sx q[1];
rz(-2.5989596) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12801816) q[0];
sx q[0];
rz(-1.7153704) q[0];
sx q[0];
rz(-1.6909063) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2130402) q[2];
sx q[2];
rz(-1.925549) q[2];
sx q[2];
rz(0.32217978) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5396351) q[1];
sx q[1];
rz(-1.4550721) q[1];
sx q[1];
rz(-1.5166111) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7067634) q[3];
sx q[3];
rz(-1.2246338) q[3];
sx q[3];
rz(-1.0258254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1594306) q[2];
sx q[2];
rz(-2.7320778) q[2];
sx q[2];
rz(-0.25992599) q[2];
rz(1.0204756) q[3];
sx q[3];
rz(-2.8838938) q[3];
sx q[3];
rz(0.79403383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4644153) q[0];
sx q[0];
rz(-2.9554415) q[0];
sx q[0];
rz(1.6625241) q[0];
rz(-1.5326477) q[1];
sx q[1];
rz(-2.1023991) q[1];
sx q[1];
rz(-2.398568) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31194132) q[0];
sx q[0];
rz(-0.52642979) q[0];
sx q[0];
rz(-1.3408324) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7588927) q[2];
sx q[2];
rz(-1.0596794) q[2];
sx q[2];
rz(1.9267043) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2496532) q[1];
sx q[1];
rz(-1.5158476) q[1];
sx q[1];
rz(1.5049388) q[1];
rz(-pi) q[2];
rz(-0.28609944) q[3];
sx q[3];
rz(-0.10933441) q[3];
sx q[3];
rz(2.5100893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6741901) q[2];
sx q[2];
rz(-0.82618606) q[2];
sx q[2];
rz(0.80545938) q[2];
rz(-1.449466) q[3];
sx q[3];
rz(-1.9129246) q[3];
sx q[3];
rz(-2.3384371) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2529124) q[0];
sx q[0];
rz(-0.54179931) q[0];
sx q[0];
rz(2.3895277) q[0];
rz(1.1532785) q[1];
sx q[1];
rz(-2.2572932) q[1];
sx q[1];
rz(2.8582252) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.008246) q[0];
sx q[0];
rz(-2.1079685) q[0];
sx q[0];
rz(1.6744012) q[0];
rz(-pi) q[1];
rz(2.9937075) q[2];
sx q[2];
rz(-0.82442946) q[2];
sx q[2];
rz(1.3651207) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4979889) q[1];
sx q[1];
rz(-2.170553) q[1];
sx q[1];
rz(1.1685755) q[1];
rz(2.649077) q[3];
sx q[3];
rz(-2.4253143) q[3];
sx q[3];
rz(-0.034958358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.89455426) q[2];
sx q[2];
rz(-0.082823195) q[2];
sx q[2];
rz(1.7029597) q[2];
rz(-2.8476207) q[3];
sx q[3];
rz(-3.1271264) q[3];
sx q[3];
rz(1.0283874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033584874) q[0];
sx q[0];
rz(-1.4172194) q[0];
sx q[0];
rz(-1.5237756) q[0];
rz(-0.53957466) q[1];
sx q[1];
rz(-2.3537666) q[1];
sx q[1];
rz(-2.981577) q[1];
rz(-2.9931184) q[2];
sx q[2];
rz(-1.9273026) q[2];
sx q[2];
rz(-2.9782563) q[2];
rz(2.7593437) q[3];
sx q[3];
rz(-1.8235689) q[3];
sx q[3];
rz(0.27128661) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

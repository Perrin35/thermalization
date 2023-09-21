OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(4.5594112) q[0];
sx q[0];
rz(8.863908) q[0];
rz(-2.0286735) q[1];
sx q[1];
rz(-1.3781883) q[1];
sx q[1];
rz(1.9265494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6998866) q[0];
sx q[0];
rz(-2.9268648) q[0];
sx q[0];
rz(-1.0617274) q[0];
rz(-1.2826074) q[2];
sx q[2];
rz(-0.93027861) q[2];
sx q[2];
rz(-3.1079907) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3455968) q[1];
sx q[1];
rz(-2.0797634) q[1];
sx q[1];
rz(-0.3791581) q[1];
rz(-2.9925572) q[3];
sx q[3];
rz(-1.3081074) q[3];
sx q[3];
rz(1.5460154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.78757301) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-0.18307486) q[2];
rz(2.7637774) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(-0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29782444) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(0.077117292) q[0];
rz(-0.33879694) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(1.5391301) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2167643) q[0];
sx q[0];
rz(-2.1568858) q[0];
sx q[0];
rz(3.0490962) q[0];
rz(-2.402926) q[2];
sx q[2];
rz(-1.7324442) q[2];
sx q[2];
rz(-2.6612298) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.96670818) q[1];
sx q[1];
rz(-1.0918573) q[1];
sx q[1];
rz(2.944988) q[1];
x q[2];
rz(-1.8061403) q[3];
sx q[3];
rz(-1.9841521) q[3];
sx q[3];
rz(-1.8542765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.845528) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(-0.65845931) q[2];
rz(2.9902839) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(-0.69491274) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.113134) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(2.7084896) q[0];
rz(1.9494879) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(-2.5862397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1028324) q[0];
sx q[0];
rz(-2.0534678) q[0];
sx q[0];
rz(2.32248) q[0];
rz(-pi) q[1];
rz(2.136134) q[2];
sx q[2];
rz(-1.8381422) q[2];
sx q[2];
rz(0.29758673) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8086116) q[1];
sx q[1];
rz(-2.4858027) q[1];
sx q[1];
rz(-1.8379184) q[1];
rz(-pi) q[2];
rz(1.2461353) q[3];
sx q[3];
rz(-2.5723296) q[3];
sx q[3];
rz(3.0666921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4042523) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(1.8910485) q[2];
rz(0.2441497) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(-1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26043949) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(2.326791) q[0];
rz(-1.3793777) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(-0.25517685) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4677306) q[0];
sx q[0];
rz(-0.47463372) q[0];
sx q[0];
rz(-2.4509096) q[0];
rz(1.9011263) q[2];
sx q[2];
rz(-0.9579881) q[2];
sx q[2];
rz(-1.4532879) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7123588) q[1];
sx q[1];
rz(-0.36839596) q[1];
sx q[1];
rz(2.1894987) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0905686) q[3];
sx q[3];
rz(-2.0912366) q[3];
sx q[3];
rz(-0.40363064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2531551) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(0.17318428) q[2];
rz(-0.52982461) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(-0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2816876) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(1.7657071) q[0];
rz(1.8638523) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(-3.0854991) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1306886) q[0];
sx q[0];
rz(-2.4405257) q[0];
sx q[0];
rz(2.5581193) q[0];
rz(-pi) q[1];
rz(0.9390097) q[2];
sx q[2];
rz(-0.54982215) q[2];
sx q[2];
rz(-2.3914571) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5406815) q[1];
sx q[1];
rz(-1.8735421) q[1];
sx q[1];
rz(1.5860228) q[1];
rz(-pi) q[2];
rz(-1.4291184) q[3];
sx q[3];
rz(-0.56483993) q[3];
sx q[3];
rz(0.39839881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(-2.7094005) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(1.4661219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0284001) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(-2.4940441) q[0];
rz(-1.8796857) q[1];
sx q[1];
rz(-1.4636661) q[1];
sx q[1];
rz(0.9544968) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19676767) q[0];
sx q[0];
rz(-1.7260572) q[0];
sx q[0];
rz(1.0374271) q[0];
rz(-1.5231832) q[2];
sx q[2];
rz(-0.63112586) q[2];
sx q[2];
rz(-0.39436755) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1986188) q[1];
sx q[1];
rz(-2.0795515) q[1];
sx q[1];
rz(-1.0191304) q[1];
rz(-pi) q[2];
rz(-2.0733842) q[3];
sx q[3];
rz(-0.32940255) q[3];
sx q[3];
rz(0.1474895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.548617) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(2.0992289) q[2];
rz(2.7029165) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(-1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-0.2838659) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(-0.74321157) q[0];
rz(-1.5076393) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(0.61002237) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9305206) q[0];
sx q[0];
rz(-1.0809582) q[0];
sx q[0];
rz(-1.2852438) q[0];
rz(1.7213983) q[2];
sx q[2];
rz(-2.4820231) q[2];
sx q[2];
rz(1.2130376) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1998132) q[1];
sx q[1];
rz(-2.4738414) q[1];
sx q[1];
rz(1.5303395) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9355572) q[3];
sx q[3];
rz(-1.1987975) q[3];
sx q[3];
rz(-0.68148617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8053749) q[2];
sx q[2];
rz(-1.4423794) q[2];
sx q[2];
rz(2.2231893) q[2];
rz(1.5911128) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(-0.38890719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.780705) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(-1.5135182) q[0];
rz(-0.52945119) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(-0.73658529) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5442473) q[0];
sx q[0];
rz(-1.5816507) q[0];
sx q[0];
rz(1.5407451) q[0];
rz(-0.10891624) q[2];
sx q[2];
rz(-1.250259) q[2];
sx q[2];
rz(-0.045217302) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8383933) q[1];
sx q[1];
rz(-2.0337078) q[1];
sx q[1];
rz(1.4569015) q[1];
x q[2];
rz(2.8183476) q[3];
sx q[3];
rz(-1.9976227) q[3];
sx q[3];
rz(-0.95110287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0344051) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(-0.68391189) q[2];
rz(-1.2290139) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1972315) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(0.18572447) q[0];
rz(0.99705237) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(0.7448147) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4208922) q[0];
sx q[0];
rz(-1.8928796) q[0];
sx q[0];
rz(-2.3548404) q[0];
x q[1];
rz(0.3955598) q[2];
sx q[2];
rz(-0.8330847) q[2];
sx q[2];
rz(-1.8849444) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.77764952) q[1];
sx q[1];
rz(-1.3476106) q[1];
sx q[1];
rz(1.9294444) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9635779) q[3];
sx q[3];
rz(-2.3142356) q[3];
sx q[3];
rz(0.82069293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0361438) q[2];
sx q[2];
rz(-2.3858586) q[2];
sx q[2];
rz(-1.194681) q[2];
rz(0.99669325) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(-0.99635807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982518) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(-0.40400305) q[0];
rz(-0.031127302) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(1.9706479) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02689657) q[0];
sx q[0];
rz(-2.7141502) q[0];
sx q[0];
rz(0.052274152) q[0];
x q[1];
rz(1.3002214) q[2];
sx q[2];
rz(-0.8136533) q[2];
sx q[2];
rz(-0.68683456) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78632894) q[1];
sx q[1];
rz(-1.5621645) q[1];
sx q[1];
rz(-1.5934056) q[1];
rz(-0.58376273) q[3];
sx q[3];
rz(-0.49920344) q[3];
sx q[3];
rz(-2.0016058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4460454) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(0.5919624) q[2];
rz(-2.5752318) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(-1.5238354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3175209) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(0.099427632) q[1];
sx q[1];
rz(-1.8933404) q[1];
sx q[1];
rz(1.0642687) q[1];
rz(-2.4026985) q[2];
sx q[2];
rz(-1.0514435) q[2];
sx q[2];
rz(-2.4098868) q[2];
rz(1.9042653) q[3];
sx q[3];
rz(-1.5555391) q[3];
sx q[3];
rz(-2.4939668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

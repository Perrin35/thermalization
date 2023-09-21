OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47473946) q[0];
sx q[0];
rz(-0.82959509) q[0];
sx q[0];
rz(0.15396804) q[0];
rz(-2.3078168) q[1];
sx q[1];
rz(-0.99234617) q[1];
sx q[1];
rz(-2.8032803) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51988039) q[0];
sx q[0];
rz(-1.8147239) q[0];
sx q[0];
rz(1.2432616) q[0];
rz(-2.7726735) q[2];
sx q[2];
rz(-0.92637617) q[2];
sx q[2];
rz(1.3117787) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.073804341) q[1];
sx q[1];
rz(-1.6238302) q[1];
sx q[1];
rz(-2.8552613) q[1];
rz(1.5428513) q[3];
sx q[3];
rz(-2.6290647) q[3];
sx q[3];
rz(-0.41748369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14264318) q[2];
sx q[2];
rz(-2.8012186) q[2];
sx q[2];
rz(1.1738698) q[2];
rz(3.0657892) q[3];
sx q[3];
rz(-1.9971763) q[3];
sx q[3];
rz(-3.048786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3409815) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(-3.0766292) q[0];
rz(2.5669572) q[1];
sx q[1];
rz(-0.42962933) q[1];
sx q[1];
rz(1.8992791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7466465) q[0];
sx q[0];
rz(-1.6992237) q[0];
sx q[0];
rz(-2.8917679) q[0];
rz(-pi) q[1];
rz(-0.1588891) q[2];
sx q[2];
rz(-0.94413589) q[2];
sx q[2];
rz(-1.6275258) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4878792) q[1];
sx q[1];
rz(-0.55760819) q[1];
sx q[1];
rz(0.30456581) q[1];
rz(0.70988016) q[3];
sx q[3];
rz(-1.9544365) q[3];
sx q[3];
rz(-0.82304728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3339281) q[2];
sx q[2];
rz(-2.0662722) q[2];
sx q[2];
rz(-0.58369613) q[2];
rz(0.57404533) q[3];
sx q[3];
rz(-2.0160926) q[3];
sx q[3];
rz(-3.0103502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7211001) q[0];
sx q[0];
rz(-0.89389602) q[0];
sx q[0];
rz(2.4131391) q[0];
rz(1.6473673) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(-2.1247991) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66660488) q[0];
sx q[0];
rz(-3.0525065) q[0];
sx q[0];
rz(-0.41390093) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74480199) q[2];
sx q[2];
rz(-2.475127) q[2];
sx q[2];
rz(1.1643861) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9079202) q[1];
sx q[1];
rz(-2.3111812) q[1];
sx q[1];
rz(0.18632142) q[1];
x q[2];
rz(-1.4025027) q[3];
sx q[3];
rz(-2.8561391) q[3];
sx q[3];
rz(1.9608378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3399405) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(-2.0920848) q[2];
rz(2.5028051) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3291572) q[0];
sx q[0];
rz(-1.8742467) q[0];
sx q[0];
rz(-1.6695492) q[0];
rz(2.4064348) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(-2.8947815) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8762159) q[0];
sx q[0];
rz(-2.2582158) q[0];
sx q[0];
rz(2.9486297) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3601801) q[2];
sx q[2];
rz(-2.0401376) q[2];
sx q[2];
rz(0.59553713) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4775131) q[1];
sx q[1];
rz(-1.941136) q[1];
sx q[1];
rz(1.5143637) q[1];
rz(0.95440063) q[3];
sx q[3];
rz(-1.3351001) q[3];
sx q[3];
rz(2.1637722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4776769) q[2];
sx q[2];
rz(-1.1185948) q[2];
sx q[2];
rz(-1.5412615) q[2];
rz(-0.70704308) q[3];
sx q[3];
rz(-2.0440846) q[3];
sx q[3];
rz(-0.88821205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33070579) q[0];
sx q[0];
rz(-2.4208477) q[0];
sx q[0];
rz(-1.3274308) q[0];
rz(-1.5785626) q[1];
sx q[1];
rz(-2.6674318) q[1];
sx q[1];
rz(2.8932103) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7800956) q[0];
sx q[0];
rz(-2.1108315) q[0];
sx q[0];
rz(-1.5310775) q[0];
rz(-pi) q[1];
rz(-2.4400473) q[2];
sx q[2];
rz(-0.24214673) q[2];
sx q[2];
rz(1.0813576) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76368139) q[1];
sx q[1];
rz(-1.7421725) q[1];
sx q[1];
rz(-0.59314368) q[1];
rz(-pi) q[2];
rz(-0.31110839) q[3];
sx q[3];
rz(-0.6725544) q[3];
sx q[3];
rz(0.5141408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43626943) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(0.79745897) q[2];
rz(-2.752839) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(-2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1335063) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(-1.7957934) q[0];
rz(-2.0603518) q[1];
sx q[1];
rz(-1.9045647) q[1];
sx q[1];
rz(-0.1246917) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42449441) q[0];
sx q[0];
rz(-1.6569123) q[0];
sx q[0];
rz(-2.9647102) q[0];
rz(-2.7331946) q[2];
sx q[2];
rz(-0.64986594) q[2];
sx q[2];
rz(-2.1085395) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.60285073) q[1];
sx q[1];
rz(-2.3067143) q[1];
sx q[1];
rz(1.49453) q[1];
rz(-pi) q[2];
rz(1.7894621) q[3];
sx q[3];
rz(-1.8018186) q[3];
sx q[3];
rz(0.71393379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.51320118) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(-0.61895269) q[2];
rz(-2.0882873) q[3];
sx q[3];
rz(-0.17799938) q[3];
sx q[3];
rz(-0.9296023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-0.58182794) q[0];
sx q[0];
rz(-1.7511837) q[0];
sx q[0];
rz(-2.0986309) q[0];
rz(2.6783121) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(2.0708864) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8541504) q[0];
sx q[0];
rz(-1.5713912) q[0];
sx q[0];
rz(2.6575412) q[0];
rz(1.7037017) q[2];
sx q[2];
rz(-1.226236) q[2];
sx q[2];
rz(0.29495707) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.13094014) q[1];
sx q[1];
rz(-2.3619235) q[1];
sx q[1];
rz(2.9995549) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8093852) q[3];
sx q[3];
rz(-1.5899961) q[3];
sx q[3];
rz(-1.7393877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7730007) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(-2.8179742) q[2];
rz(-2.1598024) q[3];
sx q[3];
rz(-0.86172813) q[3];
sx q[3];
rz(2.0543082) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48802808) q[0];
sx q[0];
rz(-1.4166778) q[0];
sx q[0];
rz(1.9352242) q[0];
rz(-1.2127097) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(-2.1941197) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0757785) q[0];
sx q[0];
rz(-1.6718719) q[0];
sx q[0];
rz(1.0321192) q[0];
rz(0.73080365) q[2];
sx q[2];
rz(-2.9060504) q[2];
sx q[2];
rz(-1.3256324) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.265043) q[1];
sx q[1];
rz(-1.3339086) q[1];
sx q[1];
rz(0.9898647) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2097589) q[3];
sx q[3];
rz(-0.79663888) q[3];
sx q[3];
rz(0.79209298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5802713) q[2];
sx q[2];
rz(-1.7611971) q[2];
sx q[2];
rz(0.46978152) q[2];
rz(-1.8404768) q[3];
sx q[3];
rz(-1.4353292) q[3];
sx q[3];
rz(0.27967134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8895421) q[0];
sx q[0];
rz(-0.38953504) q[0];
sx q[0];
rz(1.8126194) q[0];
rz(-2.3503616) q[1];
sx q[1];
rz(-2.8104517) q[1];
sx q[1];
rz(0.20283094) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3664353) q[0];
sx q[0];
rz(-3.0946819) q[0];
sx q[0];
rz(2.5527918) q[0];
x q[1];
rz(-2.9183396) q[2];
sx q[2];
rz(-1.7749783) q[2];
sx q[2];
rz(0.60165652) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.054246) q[1];
sx q[1];
rz(-0.11212238) q[1];
sx q[1];
rz(-2.0263158) q[1];
x q[2];
rz(-0.30236249) q[3];
sx q[3];
rz(-0.48067579) q[3];
sx q[3];
rz(-0.16929786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7245076) q[2];
sx q[2];
rz(-0.26792002) q[2];
sx q[2];
rz(-0.99651304) q[2];
rz(2.7881682) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(-0.80741185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0614232) q[0];
sx q[0];
rz(-2.3286979) q[0];
sx q[0];
rz(2.9598575) q[0];
rz(-0.043047992) q[1];
sx q[1];
rz(-2.4964066) q[1];
sx q[1];
rz(-0.28082401) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1998394) q[0];
sx q[0];
rz(-1.117525) q[0];
sx q[0];
rz(-1.0928632) q[0];
rz(-pi) q[1];
rz(1.7634723) q[2];
sx q[2];
rz(-1.780605) q[2];
sx q[2];
rz(3.0378621) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22051375) q[1];
sx q[1];
rz(-1.5346569) q[1];
sx q[1];
rz(1.5131348) q[1];
x q[2];
rz(1.0481846) q[3];
sx q[3];
rz(-2.7624353) q[3];
sx q[3];
rz(-0.89695938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8250371) q[2];
sx q[2];
rz(-1.8871769) q[2];
sx q[2];
rz(-0.62310702) q[2];
rz(1.0021707) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(2.5785057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-1.4939209) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(-0.87396809) q[1];
sx q[1];
rz(-2.0762434) q[1];
sx q[1];
rz(0.11003065) q[1];
rz(-0.80550823) q[2];
sx q[2];
rz(-1.1136354) q[2];
sx q[2];
rz(-1.8352933) q[2];
rz(-0.86250967) q[3];
sx q[3];
rz(-0.52048341) q[3];
sx q[3];
rz(-2.4181548) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
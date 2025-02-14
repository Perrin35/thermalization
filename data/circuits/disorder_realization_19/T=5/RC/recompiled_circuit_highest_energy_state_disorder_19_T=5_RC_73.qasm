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
rz(0.49993604) q[0];
sx q[0];
rz(1.9498107) q[0];
sx q[0];
rz(10.796588) q[0];
rz(1.6504047) q[1];
sx q[1];
rz(5.4159309) q[1];
sx q[1];
rz(9.8210788) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9306878) q[0];
sx q[0];
rz(-0.22432029) q[0];
sx q[0];
rz(0.36422642) q[0];
rz(-pi) q[1];
rz(-1.1566341) q[2];
sx q[2];
rz(-0.4067308) q[2];
sx q[2];
rz(-2.1106281) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.010155023) q[1];
sx q[1];
rz(-1.0442631) q[1];
sx q[1];
rz(0.10897691) q[1];
rz(-pi) q[2];
rz(2.7462148) q[3];
sx q[3];
rz(-1.3002743) q[3];
sx q[3];
rz(-2.6879355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0013915) q[2];
sx q[2];
rz(-0.80696693) q[2];
sx q[2];
rz(1.3414471) q[2];
rz(-2.0724824) q[3];
sx q[3];
rz(-0.1461229) q[3];
sx q[3];
rz(1.3816381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57217252) q[0];
sx q[0];
rz(-2.2241346) q[0];
sx q[0];
rz(0.65895748) q[0];
rz(-3.068889) q[1];
sx q[1];
rz(-2.0854988) q[1];
sx q[1];
rz(-1.8812077) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1457659) q[0];
sx q[0];
rz(-2.6814309) q[0];
sx q[0];
rz(-1.7547248) q[0];
x q[1];
rz(2.6953893) q[2];
sx q[2];
rz(-2.0409191) q[2];
sx q[2];
rz(-1.1371431) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.61140984) q[1];
sx q[1];
rz(-1.5059222) q[1];
sx q[1];
rz(-0.26574175) q[1];
x q[2];
rz(0.96003344) q[3];
sx q[3];
rz(-1.6690147) q[3];
sx q[3];
rz(1.4032422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8947123) q[2];
sx q[2];
rz(-2.3983045) q[2];
sx q[2];
rz(1.3512208) q[2];
rz(-2.5249935) q[3];
sx q[3];
rz(-1.0298046) q[3];
sx q[3];
rz(0.29115796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5093444) q[0];
sx q[0];
rz(-2.6743439) q[0];
sx q[0];
rz(-2.9314991) q[0];
rz(-0.083077438) q[1];
sx q[1];
rz(-1.2385085) q[1];
sx q[1];
rz(0.067213623) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2224565) q[0];
sx q[0];
rz(-1.5828776) q[0];
sx q[0];
rz(0.013508114) q[0];
rz(-pi) q[1];
rz(-2.1234197) q[2];
sx q[2];
rz(-0.45028307) q[2];
sx q[2];
rz(2.1076815) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1695425) q[1];
sx q[1];
rz(-1.9055371) q[1];
sx q[1];
rz(0.85595815) q[1];
rz(-pi) q[2];
rz(-3.0359269) q[3];
sx q[3];
rz(-0.70836954) q[3];
sx q[3];
rz(1.8745223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.82789603) q[2];
sx q[2];
rz(-1.6500429) q[2];
sx q[2];
rz(1.776604) q[2];
rz(2.7374173) q[3];
sx q[3];
rz(-1.3173236) q[3];
sx q[3];
rz(-1.1990168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2760524) q[0];
sx q[0];
rz(-1.4006389) q[0];
sx q[0];
rz(-0.47819594) q[0];
rz(2.1887691) q[1];
sx q[1];
rz(-2.9915504) q[1];
sx q[1];
rz(-0.40963867) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38423019) q[0];
sx q[0];
rz(-0.91529796) q[0];
sx q[0];
rz(-1.537498) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9443789) q[2];
sx q[2];
rz(-2.6773415) q[2];
sx q[2];
rz(1.0997538) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.19189056) q[1];
sx q[1];
rz(-0.97964761) q[1];
sx q[1];
rz(-2.3670235) q[1];
rz(-pi) q[2];
rz(1.9799819) q[3];
sx q[3];
rz(-2.7936862) q[3];
sx q[3];
rz(-2.6285841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.59565583) q[2];
sx q[2];
rz(-1.7031534) q[2];
sx q[2];
rz(-1.1379918) q[2];
rz(2.1465837) q[3];
sx q[3];
rz(-2.5840839) q[3];
sx q[3];
rz(3.0549684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.397641) q[0];
sx q[0];
rz(-1.8331563) q[0];
sx q[0];
rz(-0.27807903) q[0];
rz(1.8771578) q[1];
sx q[1];
rz(-1.8587298) q[1];
sx q[1];
rz(-1.5637195) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20533097) q[0];
sx q[0];
rz(-2.7561128) q[0];
sx q[0];
rz(1.9140713) q[0];
rz(1.9327963) q[2];
sx q[2];
rz(-0.81564192) q[2];
sx q[2];
rz(1.0175616) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5016026) q[1];
sx q[1];
rz(-1.6899365) q[1];
sx q[1];
rz(2.5646025) q[1];
x q[2];
rz(3.0487537) q[3];
sx q[3];
rz(-1.412409) q[3];
sx q[3];
rz(3.1166589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9331253) q[2];
sx q[2];
rz(-1.5310023) q[2];
sx q[2];
rz(-0.36766407) q[2];
rz(-2.0465046) q[3];
sx q[3];
rz(-1.0235267) q[3];
sx q[3];
rz(2.9663405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4472189) q[0];
sx q[0];
rz(-2.8419438) q[0];
sx q[0];
rz(1.3391986) q[0];
rz(3.0195492) q[1];
sx q[1];
rz(-1.782878) q[1];
sx q[1];
rz(2.7809714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9741518) q[0];
sx q[0];
rz(-1.1888778) q[0];
sx q[0];
rz(0.72446211) q[0];
x q[1];
rz(-1.6595938) q[2];
sx q[2];
rz(-0.42150233) q[2];
sx q[2];
rz(-0.56383609) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.569446) q[1];
sx q[1];
rz(-2.0098369) q[1];
sx q[1];
rz(-2.8296058) q[1];
rz(-1.5654605) q[3];
sx q[3];
rz(-1.3179165) q[3];
sx q[3];
rz(1.7897814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3660761) q[2];
sx q[2];
rz(-1.6785494) q[2];
sx q[2];
rz(-1.1386846) q[2];
rz(2.8876997) q[3];
sx q[3];
rz(-2.4963899) q[3];
sx q[3];
rz(2.3937288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426386) q[0];
sx q[0];
rz(-2.2518318) q[0];
sx q[0];
rz(2.1658072) q[0];
rz(2.0121393) q[1];
sx q[1];
rz(-2.6934846) q[1];
sx q[1];
rz(-1.3536369) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0746552) q[0];
sx q[0];
rz(-2.2569509) q[0];
sx q[0];
rz(1.1891436) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1466402) q[2];
sx q[2];
rz(-1.643702) q[2];
sx q[2];
rz(-2.6404501) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6401939) q[1];
sx q[1];
rz(-0.36194776) q[1];
sx q[1];
rz(-2.7622403) q[1];
rz(-2.6725572) q[3];
sx q[3];
rz(-1.9538662) q[3];
sx q[3];
rz(1.9214326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5740616) q[2];
sx q[2];
rz(-2.3109544) q[2];
sx q[2];
rz(-2.2243824) q[2];
rz(1.5926825) q[3];
sx q[3];
rz(-2.8062688) q[3];
sx q[3];
rz(-0.77939051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2061283) q[0];
sx q[0];
rz(-1.520282) q[0];
sx q[0];
rz(-0.024854831) q[0];
rz(-1.286233) q[1];
sx q[1];
rz(-1.7846466) q[1];
sx q[1];
rz(1.53055) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0760374) q[0];
sx q[0];
rz(-0.60638705) q[0];
sx q[0];
rz(0.72640149) q[0];
rz(0.31168657) q[2];
sx q[2];
rz(-0.71060813) q[2];
sx q[2];
rz(0.59618261) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87378824) q[1];
sx q[1];
rz(-0.28364983) q[1];
sx q[1];
rz(0.12445085) q[1];
rz(-pi) q[2];
rz(3.0026818) q[3];
sx q[3];
rz(-1.362097) q[3];
sx q[3];
rz(0.025996093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.020869104) q[2];
sx q[2];
rz(-1.0039696) q[2];
sx q[2];
rz(-1.8570159) q[2];
rz(2.8203216) q[3];
sx q[3];
rz(-2.4925241) q[3];
sx q[3];
rz(0.2215213) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016841737) q[0];
sx q[0];
rz(-1.9945972) q[0];
sx q[0];
rz(2.2220213) q[0];
rz(2.549767) q[1];
sx q[1];
rz(-2.070919) q[1];
sx q[1];
rz(-0.33669534) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82430912) q[0];
sx q[0];
rz(-1.705184) q[0];
sx q[0];
rz(-1.8407673) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36328237) q[2];
sx q[2];
rz(-1.4575301) q[2];
sx q[2];
rz(2.670778) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8996657) q[1];
sx q[1];
rz(-1.4414235) q[1];
sx q[1];
rz(2.632377) q[1];
rz(-pi) q[2];
rz(2.982364) q[3];
sx q[3];
rz(-1.319229) q[3];
sx q[3];
rz(-0.23345527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1153974) q[2];
sx q[2];
rz(-0.43033174) q[2];
sx q[2];
rz(2.3587312) q[2];
rz(0.10910263) q[3];
sx q[3];
rz(-2.1249873) q[3];
sx q[3];
rz(-1.8946064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1054909) q[0];
sx q[0];
rz(-1.4694659) q[0];
sx q[0];
rz(2.2176657) q[0];
rz(2.6149514) q[1];
sx q[1];
rz(-1.275332) q[1];
sx q[1];
rz(0.99871666) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0629421) q[0];
sx q[0];
rz(-1.5178158) q[0];
sx q[0];
rz(2.9880301) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0586658) q[2];
sx q[2];
rz(-1.8605663) q[2];
sx q[2];
rz(-1.7343327) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0044788) q[1];
sx q[1];
rz(-0.58885899) q[1];
sx q[1];
rz(1.6331255) q[1];
rz(-2.4630354) q[3];
sx q[3];
rz(-1.863409) q[3];
sx q[3];
rz(1.7998526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0052789) q[2];
sx q[2];
rz(-0.88999358) q[2];
sx q[2];
rz(-2.8078553) q[2];
rz(-2.2809095) q[3];
sx q[3];
rz(-1.656683) q[3];
sx q[3];
rz(3.1349643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46398791) q[0];
sx q[0];
rz(-1.9515568) q[0];
sx q[0];
rz(-2.4757181) q[0];
rz(0.55466501) q[1];
sx q[1];
rz(-1.8634836) q[1];
sx q[1];
rz(-2.9992933) q[1];
rz(-2.7016976) q[2];
sx q[2];
rz(-1.7976201) q[2];
sx q[2];
rz(1.3646361) q[2];
rz(-2.4246115) q[3];
sx q[3];
rz(-2.3752799) q[3];
sx q[3];
rz(0.3175288) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

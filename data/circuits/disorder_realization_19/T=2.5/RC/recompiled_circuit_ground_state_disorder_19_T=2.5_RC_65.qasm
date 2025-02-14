OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.98915339) q[0];
sx q[0];
rz(1.5751155) q[0];
sx q[0];
rz(10.549904) q[0];
rz(-5.2611051) q[1];
sx q[1];
rz(2.4740969) q[1];
sx q[1];
rz(11.752887) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73972244) q[0];
sx q[0];
rz(-2.7474294) q[0];
sx q[0];
rz(1.9456844) q[0];
rz(1.9314693) q[2];
sx q[2];
rz(-1.9934335) q[2];
sx q[2];
rz(0.3061184) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5578545) q[1];
sx q[1];
rz(-1.1314794) q[1];
sx q[1];
rz(1.735461) q[1];
rz(-pi) q[2];
x q[2];
rz(1.068068) q[3];
sx q[3];
rz(-0.18618551) q[3];
sx q[3];
rz(-1.26621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4890613) q[2];
sx q[2];
rz(-1.4514613) q[2];
sx q[2];
rz(2.2129464) q[2];
rz(-1.5993902) q[3];
sx q[3];
rz(-1.3336811) q[3];
sx q[3];
rz(-1.1716051) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9376675) q[0];
sx q[0];
rz(-1.7610022) q[0];
sx q[0];
rz(0.12705886) q[0];
rz(2.1584885) q[1];
sx q[1];
rz(-1.7763205) q[1];
sx q[1];
rz(-2.3703221) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5708904) q[0];
sx q[0];
rz(-1.9237907) q[0];
sx q[0];
rz(-2.3285026) q[0];
rz(-3.0808582) q[2];
sx q[2];
rz(-1.2337451) q[2];
sx q[2];
rz(-0.13523808) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.018056695) q[1];
sx q[1];
rz(-0.50790826) q[1];
sx q[1];
rz(-2.3695709) q[1];
rz(-pi) q[2];
rz(2.3985474) q[3];
sx q[3];
rz(-0.72297572) q[3];
sx q[3];
rz(1.6188177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1473006) q[2];
sx q[2];
rz(-1.2039801) q[2];
sx q[2];
rz(3.1331983) q[2];
rz(0.66347915) q[3];
sx q[3];
rz(-1.2365664) q[3];
sx q[3];
rz(-0.26507637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0405149) q[0];
sx q[0];
rz(-2.3150257) q[0];
sx q[0];
rz(0.44152942) q[0];
rz(-0.98835522) q[1];
sx q[1];
rz(-2.007273) q[1];
sx q[1];
rz(0.13557869) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39865935) q[0];
sx q[0];
rz(-0.017983111) q[0];
sx q[0];
rz(-2.6989486) q[0];
x q[1];
rz(-0.71696059) q[2];
sx q[2];
rz(-1.8646984) q[2];
sx q[2];
rz(0.64168054) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24568711) q[1];
sx q[1];
rz(-1.0789597) q[1];
sx q[1];
rz(-2.7373284) q[1];
rz(-pi) q[2];
rz(1.3531923) q[3];
sx q[3];
rz(-2.3953468) q[3];
sx q[3];
rz(-1.7640132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.93379891) q[2];
sx q[2];
rz(-1.4321045) q[2];
sx q[2];
rz(-0.03820339) q[2];
rz(2.6162052) q[3];
sx q[3];
rz(-0.62756687) q[3];
sx q[3];
rz(2.4562522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66048375) q[0];
sx q[0];
rz(-0.79367343) q[0];
sx q[0];
rz(0.61087459) q[0];
rz(1.8065709) q[1];
sx q[1];
rz(-1.7673312) q[1];
sx q[1];
rz(-2.9023721) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4668381) q[0];
sx q[0];
rz(-1.309899) q[0];
sx q[0];
rz(0.99715085) q[0];
x q[1];
rz(-2.4609354) q[2];
sx q[2];
rz(-0.58154642) q[2];
sx q[2];
rz(-2.8342063) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1924421) q[1];
sx q[1];
rz(-0.56862967) q[1];
sx q[1];
rz(0.62921625) q[1];
rz(-1.3705105) q[3];
sx q[3];
rz(-2.0242175) q[3];
sx q[3];
rz(3.0016921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4227582) q[2];
sx q[2];
rz(-1.4449291) q[2];
sx q[2];
rz(-0.0017496721) q[2];
rz(2.8478029) q[3];
sx q[3];
rz(-1.1894476) q[3];
sx q[3];
rz(2.8747115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2413498) q[0];
sx q[0];
rz(-1.3768063) q[0];
sx q[0];
rz(2.2985261) q[0];
rz(-2.8040366) q[1];
sx q[1];
rz(-1.866021) q[1];
sx q[1];
rz(1.5379803) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085097236) q[0];
sx q[0];
rz(-0.81696327) q[0];
sx q[0];
rz(-0.2635862) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29422167) q[2];
sx q[2];
rz(-2.9595642) q[2];
sx q[2];
rz(-0.95532571) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.91572054) q[1];
sx q[1];
rz(-0.32494007) q[1];
sx q[1];
rz(-2.48808) q[1];
rz(-1.6513651) q[3];
sx q[3];
rz(-1.3617097) q[3];
sx q[3];
rz(2.5029456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7097077) q[2];
sx q[2];
rz(-2.0543435) q[2];
sx q[2];
rz(-2.0951927) q[2];
rz(0.31878582) q[3];
sx q[3];
rz(-0.71109486) q[3];
sx q[3];
rz(-0.54767245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043561291) q[0];
sx q[0];
rz(-1.2579608) q[0];
sx q[0];
rz(0.64796722) q[0];
rz(1.7421534) q[1];
sx q[1];
rz(-1.49767) q[1];
sx q[1];
rz(1.8125777) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32008313) q[0];
sx q[0];
rz(-1.3504986) q[0];
sx q[0];
rz(-3.0295323) q[0];
rz(-pi) q[1];
rz(1.7850661) q[2];
sx q[2];
rz(-1.8545215) q[2];
sx q[2];
rz(1.2905215) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.700579) q[1];
sx q[1];
rz(-0.77475905) q[1];
sx q[1];
rz(-2.0225594) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5446114) q[3];
sx q[3];
rz(-0.59904811) q[3];
sx q[3];
rz(2.0839276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7602188) q[2];
sx q[2];
rz(-1.9332644) q[2];
sx q[2];
rz(1.5488497) q[2];
rz(-0.069843944) q[3];
sx q[3];
rz(-1.2589688) q[3];
sx q[3];
rz(2.2729592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1182564) q[0];
sx q[0];
rz(-1.2160439) q[0];
sx q[0];
rz(-0.94183952) q[0];
rz(2.9815004) q[1];
sx q[1];
rz(-1.4804877) q[1];
sx q[1];
rz(-2.9494185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20142444) q[0];
sx q[0];
rz(-1.3469101) q[0];
sx q[0];
rz(3.1046449) q[0];
x q[1];
rz(-1.0490369) q[2];
sx q[2];
rz(-1.6055487) q[2];
sx q[2];
rz(3.1086189) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3698764) q[1];
sx q[1];
rz(-1.9100921) q[1];
sx q[1];
rz(0.75668613) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0683769) q[3];
sx q[3];
rz(-1.5901774) q[3];
sx q[3];
rz(-1.3598702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3840702) q[2];
sx q[2];
rz(-1.9033868) q[2];
sx q[2];
rz(2.7080217) q[2];
rz(-1.5571669) q[3];
sx q[3];
rz(-1.6470563) q[3];
sx q[3];
rz(0.43237329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6381391) q[0];
sx q[0];
rz(-1.764955) q[0];
sx q[0];
rz(1.7973416) q[0];
rz(1.9380219) q[1];
sx q[1];
rz(-1.9529587) q[1];
sx q[1];
rz(1.3628091) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2575131) q[0];
sx q[0];
rz(-1.5410893) q[0];
sx q[0];
rz(-0.021935181) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3154738) q[2];
sx q[2];
rz(-1.4419793) q[2];
sx q[2];
rz(1.4908189) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.92719936) q[1];
sx q[1];
rz(-1.7991369) q[1];
sx q[1];
rz(-1.4424999) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1892494) q[3];
sx q[3];
rz(-2.4182712) q[3];
sx q[3];
rz(1.2398694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5726996) q[2];
sx q[2];
rz(-0.77304825) q[2];
sx q[2];
rz(-2.1577238) q[2];
rz(2.900506) q[3];
sx q[3];
rz(-2.2086996) q[3];
sx q[3];
rz(-0.69055313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47138658) q[0];
sx q[0];
rz(-2.3455878) q[0];
sx q[0];
rz(0.27467003) q[0];
rz(-1.7999016) q[1];
sx q[1];
rz(-0.62961737) q[1];
sx q[1];
rz(1.0460269) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5513902) q[0];
sx q[0];
rz(-2.1596163) q[0];
sx q[0];
rz(1.1776275) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3730896) q[2];
sx q[2];
rz(-1.4914163) q[2];
sx q[2];
rz(-2.6153836) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.53262344) q[1];
sx q[1];
rz(-1.5057505) q[1];
sx q[1];
rz(0.58397997) q[1];
rz(-3.0996383) q[3];
sx q[3];
rz(-1.0681515) q[3];
sx q[3];
rz(1.5273389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.86299738) q[2];
sx q[2];
rz(-1.7509165) q[2];
sx q[2];
rz(0.39946237) q[2];
rz(-2.5755889) q[3];
sx q[3];
rz(-2.1750906) q[3];
sx q[3];
rz(2.32617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852916) q[0];
sx q[0];
rz(-1.7221907) q[0];
sx q[0];
rz(-2.5592819) q[0];
rz(-0.18320006) q[1];
sx q[1];
rz(-0.87545005) q[1];
sx q[1];
rz(2.3351672) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69003478) q[0];
sx q[0];
rz(-1.4831044) q[0];
sx q[0];
rz(1.466485) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66566531) q[2];
sx q[2];
rz(-2.7241926) q[2];
sx q[2];
rz(1.1102499) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.4848856) q[1];
sx q[1];
rz(-0.84204145) q[1];
sx q[1];
rz(1.2970112) q[1];
x q[2];
rz(-2.7943198) q[3];
sx q[3];
rz(-1.4662305) q[3];
sx q[3];
rz(-2.4109651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9091984) q[2];
sx q[2];
rz(-0.71467233) q[2];
sx q[2];
rz(-0.8052899) q[2];
rz(-0.14690873) q[3];
sx q[3];
rz(-2.3555136) q[3];
sx q[3];
rz(0.73827353) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.88937) q[0];
sx q[0];
rz(-1.7563553) q[0];
sx q[0];
rz(1.8808543) q[0];
rz(1.5421142) q[1];
sx q[1];
rz(-1.5442994) q[1];
sx q[1];
rz(-1.50179) q[1];
rz(1.7003822) q[2];
sx q[2];
rz(-0.72740462) q[2];
sx q[2];
rz(-1.8047756) q[2];
rz(2.6977758) q[3];
sx q[3];
rz(-2.4469821) q[3];
sx q[3];
rz(-1.2367005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(2.3693585) q[0];
sx q[0];
rz(-1.9263664) q[0];
sx q[0];
rz(-1.9908494) q[0];
rz(-0.28252217) q[1];
sx q[1];
rz(4.2673586) q[1];
sx q[1];
rz(10.537108) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.756506) q[0];
sx q[0];
rz(-1.4539363) q[0];
sx q[0];
rz(-2.8782842) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4006279) q[2];
sx q[2];
rz(-2.552816) q[2];
sx q[2];
rz(-0.44786225) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9291275) q[1];
sx q[1];
rz(-2.0501137) q[1];
sx q[1];
rz(2.1838837) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3521117) q[3];
sx q[3];
rz(-2.1095358) q[3];
sx q[3];
rz(-2.1110502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7925966) q[2];
sx q[2];
rz(-2.0670321) q[2];
sx q[2];
rz(-1.7462771) q[2];
rz(-0.058517728) q[3];
sx q[3];
rz(-1.7250215) q[3];
sx q[3];
rz(-1.831656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.42404744) q[0];
sx q[0];
rz(-2.738364) q[0];
sx q[0];
rz(2.6918217) q[0];
rz(-2.6768118) q[1];
sx q[1];
rz(-1.3100781) q[1];
sx q[1];
rz(-2.8890077) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16134444) q[0];
sx q[0];
rz(-0.35156116) q[0];
sx q[0];
rz(-1.5241966) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0130326) q[2];
sx q[2];
rz(-0.95880836) q[2];
sx q[2];
rz(2.649518) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3623533) q[1];
sx q[1];
rz(-1.5679411) q[1];
sx q[1];
rz(-1.8523907) q[1];
rz(-pi) q[2];
rz(-0.59201805) q[3];
sx q[3];
rz(-2.238174) q[3];
sx q[3];
rz(-0.84718466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2437336) q[2];
sx q[2];
rz(-0.89577883) q[2];
sx q[2];
rz(-0.016156999) q[2];
rz(-2.2195418) q[3];
sx q[3];
rz(-0.99370304) q[3];
sx q[3];
rz(-3.0097358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.09318) q[0];
sx q[0];
rz(-1.484363) q[0];
sx q[0];
rz(-2.4600001) q[0];
rz(-2.5438578) q[1];
sx q[1];
rz(-1.455541) q[1];
sx q[1];
rz(-2.8173503) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7549122) q[0];
sx q[0];
rz(-2.5752299) q[0];
sx q[0];
rz(2.4709111) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6347972) q[2];
sx q[2];
rz(-1.629771) q[2];
sx q[2];
rz(2.5304194) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5983008) q[1];
sx q[1];
rz(-0.17379119) q[1];
sx q[1];
rz(1.8164746) q[1];
rz(-pi) q[2];
rz(-2.1799654) q[3];
sx q[3];
rz(-1.4318083) q[3];
sx q[3];
rz(2.9442996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9165667) q[2];
sx q[2];
rz(-2.7498507) q[2];
sx q[2];
rz(2.0908835) q[2];
rz(2.4010036) q[3];
sx q[3];
rz(-1.2935484) q[3];
sx q[3];
rz(0.061323969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6447613) q[0];
sx q[0];
rz(-2.3065541) q[0];
sx q[0];
rz(-2.187425) q[0];
rz(1.1369811) q[1];
sx q[1];
rz(-1.6258207) q[1];
sx q[1];
rz(2.383393) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9228464) q[0];
sx q[0];
rz(-0.2740261) q[0];
sx q[0];
rz(-1.7195527) q[0];
x q[1];
rz(-0.21110837) q[2];
sx q[2];
rz(-2.3866215) q[2];
sx q[2];
rz(-1.8051219) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.41685795) q[1];
sx q[1];
rz(-1.6685889) q[1];
sx q[1];
rz(2.2657402) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73885804) q[3];
sx q[3];
rz(-1.5631873) q[3];
sx q[3];
rz(1.7926907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.053146426) q[2];
sx q[2];
rz(-2.5909178) q[2];
sx q[2];
rz(-1.3147563) q[2];
rz(-0.30086073) q[3];
sx q[3];
rz(-1.8010537) q[3];
sx q[3];
rz(-0.68527591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(1.143173) q[0];
sx q[0];
rz(-1.5564065) q[0];
sx q[0];
rz(1.6402624) q[0];
rz(0.38652626) q[1];
sx q[1];
rz(-2.0730348) q[1];
sx q[1];
rz(-1.5949257) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7435369) q[0];
sx q[0];
rz(-0.72303444) q[0];
sx q[0];
rz(-2.8788752) q[0];
rz(-pi) q[1];
rz(2.8168786) q[2];
sx q[2];
rz(-0.98718671) q[2];
sx q[2];
rz(1.8167855) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9846696) q[1];
sx q[1];
rz(-2.700173) q[1];
sx q[1];
rz(-0.19766165) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5677979) q[3];
sx q[3];
rz(-1.6125896) q[3];
sx q[3];
rz(-2.6190663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.74653643) q[2];
sx q[2];
rz(-2.655513) q[2];
sx q[2];
rz(-2.0716095) q[2];
rz(2.2561031) q[3];
sx q[3];
rz(-2.077379) q[3];
sx q[3];
rz(1.1522256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.92538658) q[0];
sx q[0];
rz(-2.8031271) q[0];
sx q[0];
rz(-2.0881407) q[0];
rz(2.9523051) q[1];
sx q[1];
rz(-0.82866755) q[1];
sx q[1];
rz(-1.6493571) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5976335) q[0];
sx q[0];
rz(-1.9372371) q[0];
sx q[0];
rz(2.9935915) q[0];
x q[1];
rz(-1.8771421) q[2];
sx q[2];
rz(-1.7736846) q[2];
sx q[2];
rz(-2.5032147) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.63000667) q[1];
sx q[1];
rz(-1.4714071) q[1];
sx q[1];
rz(0.50390999) q[1];
x q[2];
rz(-2.3334741) q[3];
sx q[3];
rz(-1.5424154) q[3];
sx q[3];
rz(-2.89464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.80387992) q[2];
sx q[2];
rz(-2.1264075) q[2];
sx q[2];
rz(-1.052915) q[2];
rz(-3.0247011) q[3];
sx q[3];
rz(-1.0630853) q[3];
sx q[3];
rz(1.4611999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9601124) q[0];
sx q[0];
rz(-0.57992613) q[0];
sx q[0];
rz(0.54452288) q[0];
rz(-2.8879884) q[1];
sx q[1];
rz(-0.2519775) q[1];
sx q[1];
rz(2.8869218) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38886759) q[0];
sx q[0];
rz(-1.6427186) q[0];
sx q[0];
rz(0.0062463721) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4843226) q[2];
sx q[2];
rz(-1.9625798) q[2];
sx q[2];
rz(-2.0287152) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0510301) q[1];
sx q[1];
rz(-1.6335618) q[1];
sx q[1];
rz(0.81845567) q[1];
rz(-pi) q[2];
rz(-2.3319844) q[3];
sx q[3];
rz(-2.0008464) q[3];
sx q[3];
rz(0.82287091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9165245) q[2];
sx q[2];
rz(-2.2038286) q[2];
sx q[2];
rz(2.7226105) q[2];
rz(0.032912832) q[3];
sx q[3];
rz(-1.2337647) q[3];
sx q[3];
rz(-2.029443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58881775) q[0];
sx q[0];
rz(-2.3989615) q[0];
sx q[0];
rz(-1.8137929) q[0];
rz(2.1090419) q[1];
sx q[1];
rz(-1.3464709) q[1];
sx q[1];
rz(2.5520777) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27145984) q[0];
sx q[0];
rz(-2.1147595) q[0];
sx q[0];
rz(-2.5556295) q[0];
x q[1];
rz(1.8399182) q[2];
sx q[2];
rz(-1.8523714) q[2];
sx q[2];
rz(0.21965227) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2736328) q[1];
sx q[1];
rz(-0.74791779) q[1];
sx q[1];
rz(-0.91435097) q[1];
rz(2.74624) q[3];
sx q[3];
rz(-1.1046003) q[3];
sx q[3];
rz(0.485093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0716268) q[2];
sx q[2];
rz(-2.4667141) q[2];
sx q[2];
rz(0.91342941) q[2];
rz(-2.6895798) q[3];
sx q[3];
rz(-1.858523) q[3];
sx q[3];
rz(-1.7970596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2695059) q[0];
sx q[0];
rz(-1.0259314) q[0];
sx q[0];
rz(1.6312697) q[0];
rz(-1.3144846) q[1];
sx q[1];
rz(-2.1035106) q[1];
sx q[1];
rz(0.41183919) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.340419) q[0];
sx q[0];
rz(-2.1864751) q[0];
sx q[0];
rz(-2.8007617) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8226003) q[2];
sx q[2];
rz(-1.5249494) q[2];
sx q[2];
rz(-0.09825966) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.85535586) q[1];
sx q[1];
rz(-2.0986631) q[1];
sx q[1];
rz(-1.1638454) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7564323) q[3];
sx q[3];
rz(-1.876014) q[3];
sx q[3];
rz(0.45476828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0044378) q[2];
sx q[2];
rz(-1.6496481) q[2];
sx q[2];
rz(-2.906666) q[2];
rz(2.2233985) q[3];
sx q[3];
rz(-0.72489679) q[3];
sx q[3];
rz(1.0880067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.039577) q[0];
sx q[0];
rz(-0.35440847) q[0];
sx q[0];
rz(1.2603731) q[0];
rz(-1.4541939) q[1];
sx q[1];
rz(-1.4581542) q[1];
sx q[1];
rz(-1.4448602) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27118998) q[0];
sx q[0];
rz(-1.9810441) q[0];
sx q[0];
rz(-0.48727401) q[0];
rz(-0.95506041) q[2];
sx q[2];
rz(-0.79074016) q[2];
sx q[2];
rz(-0.093971595) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9493172) q[1];
sx q[1];
rz(-0.85964627) q[1];
sx q[1];
rz(1.7356731) q[1];
rz(-2.7637611) q[3];
sx q[3];
rz(-1.2595385) q[3];
sx q[3];
rz(-1.6938083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5789648) q[2];
sx q[2];
rz(-2.8426888) q[2];
sx q[2];
rz(0.12358269) q[2];
rz(0.28892162) q[3];
sx q[3];
rz(-1.8653423) q[3];
sx q[3];
rz(0.45946521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4041483) q[0];
sx q[0];
rz(-1.7973719) q[0];
sx q[0];
rz(1.6112882) q[0];
rz(-2.1743446) q[1];
sx q[1];
rz(-2.1687242) q[1];
sx q[1];
rz(0.8868934) q[1];
rz(1.6254454) q[2];
sx q[2];
rz(-0.9722295) q[2];
sx q[2];
rz(0.35884919) q[2];
rz(2.2933949) q[3];
sx q[3];
rz(-2.316939) q[3];
sx q[3];
rz(1.4311781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

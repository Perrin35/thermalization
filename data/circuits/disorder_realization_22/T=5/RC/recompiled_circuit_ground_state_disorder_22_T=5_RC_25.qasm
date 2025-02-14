OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8393505) q[0];
sx q[0];
rz(-1.8678764) q[0];
sx q[0];
rz(-0.94925517) q[0];
rz(-1.0957837) q[1];
sx q[1];
rz(2.1646808) q[1];
sx q[1];
rz(10.77471) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.183584) q[0];
sx q[0];
rz(-1.4877948) q[0];
sx q[0];
rz(-2.6380093) q[0];
rz(2.737713) q[2];
sx q[2];
rz(-2.7821861) q[2];
sx q[2];
rz(3.0702116) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9266859) q[1];
sx q[1];
rz(-1.6388571) q[1];
sx q[1];
rz(2.7665274) q[1];
x q[2];
rz(1.024035) q[3];
sx q[3];
rz(-2.6461678) q[3];
sx q[3];
rz(1.6273496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6076516) q[2];
sx q[2];
rz(-2.0484296) q[2];
sx q[2];
rz(-2.9489813) q[2];
rz(1.2827778) q[3];
sx q[3];
rz(-2.0359813) q[3];
sx q[3];
rz(-0.58853308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31805661) q[0];
sx q[0];
rz(-0.41724351) q[0];
sx q[0];
rz(0.70660025) q[0];
rz(0.63287863) q[1];
sx q[1];
rz(-2.0426079) q[1];
sx q[1];
rz(0.28876567) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59290409) q[0];
sx q[0];
rz(-3.1389154) q[0];
sx q[0];
rz(0.38627465) q[0];
x q[1];
rz(-0.20769221) q[2];
sx q[2];
rz(-1.7294267) q[2];
sx q[2];
rz(2.6191448) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0763467) q[1];
sx q[1];
rz(-1.7829121) q[1];
sx q[1];
rz(2.8492623) q[1];
x q[2];
rz(-0.76008009) q[3];
sx q[3];
rz(-1.2562498) q[3];
sx q[3];
rz(-3.1414978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4660945) q[2];
sx q[2];
rz(-1.0885295) q[2];
sx q[2];
rz(-1.6509854) q[2];
rz(-2.364482) q[3];
sx q[3];
rz(-1.4584352) q[3];
sx q[3];
rz(-1.3597663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31558388) q[0];
sx q[0];
rz(-1.4588139) q[0];
sx q[0];
rz(-0.43651954) q[0];
rz(-0.79477683) q[1];
sx q[1];
rz(-1.3705285) q[1];
sx q[1];
rz(-0.93903843) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2344367) q[0];
sx q[0];
rz(-2.5946989) q[0];
sx q[0];
rz(2.7476385) q[0];
x q[1];
rz(0.58231488) q[2];
sx q[2];
rz(-2.2852201) q[2];
sx q[2];
rz(1.1143052) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3311829) q[1];
sx q[1];
rz(-0.66920921) q[1];
sx q[1];
rz(-0.52468467) q[1];
rz(-1.0163937) q[3];
sx q[3];
rz(-1.7590343) q[3];
sx q[3];
rz(2.7013626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.44094917) q[2];
sx q[2];
rz(-2.078853) q[2];
sx q[2];
rz(0.6353333) q[2];
rz(-2.3144531) q[3];
sx q[3];
rz(-0.45378903) q[3];
sx q[3];
rz(1.2686096) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24469911) q[0];
sx q[0];
rz(-0.06299717) q[0];
sx q[0];
rz(2.6196106) q[0];
rz(0.43102795) q[1];
sx q[1];
rz(-1.6072175) q[1];
sx q[1];
rz(0.65188754) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4341457) q[0];
sx q[0];
rz(-1.0672671) q[0];
sx q[0];
rz(-1.8412526) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0235223) q[2];
sx q[2];
rz(-0.80308417) q[2];
sx q[2];
rz(-1.3241763) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3053599) q[1];
sx q[1];
rz(-1.1916324) q[1];
sx q[1];
rz(-3.1268478) q[1];
rz(0.68199358) q[3];
sx q[3];
rz(-2.3122283) q[3];
sx q[3];
rz(-0.030205848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7151457) q[2];
sx q[2];
rz(-1.8603674) q[2];
sx q[2];
rz(-0.78488266) q[2];
rz(-0.52810413) q[3];
sx q[3];
rz(-1.8485066) q[3];
sx q[3];
rz(1.2877134) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2454979) q[0];
sx q[0];
rz(-2.446785) q[0];
sx q[0];
rz(0.78952638) q[0];
rz(2.6441669) q[1];
sx q[1];
rz(-1.0405633) q[1];
sx q[1];
rz(-1.8400037) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0987463) q[0];
sx q[0];
rz(-1.3912956) q[0];
sx q[0];
rz(-0.60447201) q[0];
x q[1];
rz(-1.6315494) q[2];
sx q[2];
rz(-1.9092602) q[2];
sx q[2];
rz(-2.0450704) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9250946) q[1];
sx q[1];
rz(-1.7431694) q[1];
sx q[1];
rz(-2.4045776) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1188854) q[3];
sx q[3];
rz(-1.9380261) q[3];
sx q[3];
rz(-2.2677088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.54399458) q[2];
sx q[2];
rz(-1.602729) q[2];
sx q[2];
rz(0.27628118) q[2];
rz(-0.9497408) q[3];
sx q[3];
rz(-0.26849982) q[3];
sx q[3];
rz(2.5883519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-2.5211869) q[0];
sx q[0];
rz(-0.036245417) q[0];
sx q[0];
rz(2.1976443) q[0];
rz(-1.2983407) q[1];
sx q[1];
rz(-1.0669758) q[1];
sx q[1];
rz(2.3640769) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28431129) q[0];
sx q[0];
rz(-2.4696283) q[0];
sx q[0];
rz(0.59055968) q[0];
rz(-2.8526659) q[2];
sx q[2];
rz(-0.75519604) q[2];
sx q[2];
rz(2.4946405) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4247269) q[1];
sx q[1];
rz(-1.9391372) q[1];
sx q[1];
rz(-1.7067545) q[1];
x q[2];
rz(0.30413587) q[3];
sx q[3];
rz(-1.5392996) q[3];
sx q[3];
rz(2.568752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.62766084) q[2];
sx q[2];
rz(-1.7646503) q[2];
sx q[2];
rz(2.7505006) q[2];
rz(2.5202461) q[3];
sx q[3];
rz(-2.4070599) q[3];
sx q[3];
rz(2.3656316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67052996) q[0];
sx q[0];
rz(-0.10469086) q[0];
sx q[0];
rz(1.4290357) q[0];
rz(-2.1757226) q[1];
sx q[1];
rz(-1.8873676) q[1];
sx q[1];
rz(-0.78580725) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2656276) q[0];
sx q[0];
rz(-0.64564359) q[0];
sx q[0];
rz(3.0558048) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4072717) q[2];
sx q[2];
rz(-1.87748) q[2];
sx q[2];
rz(-1.9421362) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8829086) q[1];
sx q[1];
rz(-1.488954) q[1];
sx q[1];
rz(-2.6805582) q[1];
rz(1.8604516) q[3];
sx q[3];
rz(-1.148734) q[3];
sx q[3];
rz(-0.45988032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.69661951) q[2];
sx q[2];
rz(-0.61903054) q[2];
sx q[2];
rz(-2.6444198) q[2];
rz(0.87388006) q[3];
sx q[3];
rz(-2.0495575) q[3];
sx q[3];
rz(1.9397651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9484321) q[0];
sx q[0];
rz(-2.8346297) q[0];
sx q[0];
rz(2.4440785) q[0];
rz(0.3802309) q[1];
sx q[1];
rz(-1.1153328) q[1];
sx q[1];
rz(3.0013705) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74497667) q[0];
sx q[0];
rz(-0.19752398) q[0];
sx q[0];
rz(2.7382988) q[0];
x q[1];
rz(-1.6096787) q[2];
sx q[2];
rz(-1.1302396) q[2];
sx q[2];
rz(1.3300542) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4932053) q[1];
sx q[1];
rz(-2.3304061) q[1];
sx q[1];
rz(0.56583515) q[1];
rz(-pi) q[2];
rz(0.41922064) q[3];
sx q[3];
rz(-0.44525075) q[3];
sx q[3];
rz(2.3589596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.85072881) q[2];
sx q[2];
rz(-1.2173434) q[2];
sx q[2];
rz(-2.6503837) q[2];
rz(2.9344007) q[3];
sx q[3];
rz(-2.5184293) q[3];
sx q[3];
rz(-1.7306958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17289138) q[0];
sx q[0];
rz(-1.7270813) q[0];
sx q[0];
rz(0.15039314) q[0];
rz(-0.70101678) q[1];
sx q[1];
rz(-1.0944518) q[1];
sx q[1];
rz(-1.6640123) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9532861) q[0];
sx q[0];
rz(-2.4948848) q[0];
sx q[0];
rz(-2.9394399) q[0];
rz(-pi) q[1];
rz(0.62297627) q[2];
sx q[2];
rz(-2.9439363) q[2];
sx q[2];
rz(2.4923639) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7733998) q[1];
sx q[1];
rz(-2.9559694) q[1];
sx q[1];
rz(0.61180656) q[1];
rz(-pi) q[2];
rz(0.27235741) q[3];
sx q[3];
rz(-2.2873169) q[3];
sx q[3];
rz(0.89064939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.56665862) q[2];
sx q[2];
rz(-0.40859544) q[2];
sx q[2];
rz(1.6179786) q[2];
rz(-0.59897113) q[3];
sx q[3];
rz(-0.50061575) q[3];
sx q[3];
rz(0.18794255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69342518) q[0];
sx q[0];
rz(-0.42891112) q[0];
sx q[0];
rz(-1.5018139) q[0];
rz(-2.2046454) q[1];
sx q[1];
rz(-2.1138771) q[1];
sx q[1];
rz(-2.1243336) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.992998) q[0];
sx q[0];
rz(-2.5879597) q[0];
sx q[0];
rz(-2.3675336) q[0];
rz(-2.4577519) q[2];
sx q[2];
rz(-3.0061887) q[2];
sx q[2];
rz(1.1410448) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.63314519) q[1];
sx q[1];
rz(-2.8775551) q[1];
sx q[1];
rz(1.4435084) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7211765) q[3];
sx q[3];
rz(-2.6767593) q[3];
sx q[3];
rz(-2.5878588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0246058) q[2];
sx q[2];
rz(-2.459343) q[2];
sx q[2];
rz(-1.619722) q[2];
rz(1.6541325) q[3];
sx q[3];
rz(-1.6493075) q[3];
sx q[3];
rz(-1.3967995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7617154) q[0];
sx q[0];
rz(-1.7424255) q[0];
sx q[0];
rz(-2.4095834) q[0];
rz(-1.4749745) q[1];
sx q[1];
rz(-1.2154308) q[1];
sx q[1];
rz(-2.348127) q[1];
rz(-0.39948612) q[2];
sx q[2];
rz(-2.1151092) q[2];
sx q[2];
rz(-1.5110037) q[2];
rz(1.2670682) q[3];
sx q[3];
rz(-1.8770185) q[3];
sx q[3];
rz(0.087531623) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

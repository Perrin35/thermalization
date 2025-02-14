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
rz(0.762142) q[0];
sx q[0];
rz(-2.6317959) q[0];
sx q[0];
rz(2.9708746) q[0];
rz(-0.40387965) q[2];
sx q[2];
rz(-0.35940659) q[2];
sx q[2];
rz(0.071381005) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7589345) q[1];
sx q[1];
rz(-1.1966424) q[1];
sx q[1];
rz(-1.4976682) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1175577) q[3];
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
rz(-1.6076516) q[2];
sx q[2];
rz(-1.093163) q[2];
sx q[2];
rz(0.1926113) q[2];
rz(1.2827778) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31805661) q[0];
sx q[0];
rz(-2.7243491) q[0];
sx q[0];
rz(0.70660025) q[0];
rz(2.508714) q[1];
sx q[1];
rz(-2.0426079) q[1];
sx q[1];
rz(-0.28876567) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59290409) q[0];
sx q[0];
rz(-3.1389154) q[0];
sx q[0];
rz(2.755318) q[0];
x q[1];
rz(2.9339004) q[2];
sx q[2];
rz(-1.7294267) q[2];
sx q[2];
rz(-0.52244782) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0367248) q[1];
sx q[1];
rz(-2.782208) q[1];
sx q[1];
rz(0.64175989) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1488647) q[3];
sx q[3];
rz(-2.2852118) q[3];
sx q[3];
rz(-1.2846636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.67549813) q[2];
sx q[2];
rz(-2.0530632) q[2];
sx q[2];
rz(-1.4906073) q[2];
rz(-2.364482) q[3];
sx q[3];
rz(-1.6831574) q[3];
sx q[3];
rz(-1.7818264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8260088) q[0];
sx q[0];
rz(-1.4588139) q[0];
sx q[0];
rz(-0.43651954) q[0];
rz(-0.79477683) q[1];
sx q[1];
rz(-1.3705285) q[1];
sx q[1];
rz(2.2025542) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7815112) q[0];
sx q[0];
rz(-1.0699125) q[0];
sx q[0];
rz(1.341218) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.375023) q[2];
sx q[2];
rz(-1.9992644) q[2];
sx q[2];
rz(-0.049190532) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.66575501) q[1];
sx q[1];
rz(-1.2547973) q[1];
sx q[1];
rz(2.5412987) q[1];
x q[2];
rz(1.9179929) q[3];
sx q[3];
rz(-2.5592862) q[3];
sx q[3];
rz(0.83707419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7006435) q[2];
sx q[2];
rz(-1.0627397) q[2];
sx q[2];
rz(-2.5062594) q[2];
rz(-2.3144531) q[3];
sx q[3];
rz(-2.6878036) q[3];
sx q[3];
rz(-1.2686096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24469911) q[0];
sx q[0];
rz(-0.06299717) q[0];
sx q[0];
rz(-2.6196106) q[0];
rz(2.7105647) q[1];
sx q[1];
rz(-1.6072175) q[1];
sx q[1];
rz(2.4897051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4341457) q[0];
sx q[0];
rz(-2.0743255) q[0];
sx q[0];
rz(1.8412526) q[0];
rz(-0.4944369) q[2];
sx q[2];
rz(-2.2324413) q[2];
sx q[2];
rz(-0.60397691) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3451772) q[1];
sx q[1];
rz(-2.762156) q[1];
sx q[1];
rz(1.6077843) q[1];
rz(-pi) q[2];
rz(2.4595991) q[3];
sx q[3];
rz(-0.82936433) q[3];
sx q[3];
rz(-0.030205848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.42644694) q[2];
sx q[2];
rz(-1.8603674) q[2];
sx q[2];
rz(2.35671) q[2];
rz(2.6134885) q[3];
sx q[3];
rz(-1.2930861) q[3];
sx q[3];
rz(-1.2877134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2454979) q[0];
sx q[0];
rz(-0.69480768) q[0];
sx q[0];
rz(2.3520663) q[0];
rz(-2.6441669) q[1];
sx q[1];
rz(-2.1010294) q[1];
sx q[1];
rz(1.3015889) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0987463) q[0];
sx q[0];
rz(-1.3912956) q[0];
sx q[0];
rz(0.60447201) q[0];
rz(-2.80255) q[2];
sx q[2];
rz(-1.5134939) q[2];
sx q[2];
rz(0.49446854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6005206) q[1];
sx q[1];
rz(-0.75319911) q[1];
sx q[1];
rz(-0.25347565) q[1];
x q[2];
rz(3.1188854) q[3];
sx q[3];
rz(-1.9380261) q[3];
sx q[3];
rz(-0.87388384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54399458) q[2];
sx q[2];
rz(-1.5388637) q[2];
sx q[2];
rz(2.8653115) q[2];
rz(-2.1918519) q[3];
sx q[3];
rz(-2.8730928) q[3];
sx q[3];
rz(-0.55324078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.62040579) q[0];
sx q[0];
rz(-0.036245417) q[0];
sx q[0];
rz(-0.94394839) q[0];
rz(1.8432519) q[1];
sx q[1];
rz(-2.0746168) q[1];
sx q[1];
rz(-2.3640769) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.992618) q[0];
sx q[0];
rz(-1.0273522) q[0];
sx q[0];
rz(1.1538366) q[0];
rz(0.73410122) q[2];
sx q[2];
rz(-1.7673552) q[2];
sx q[2];
rz(-0.71069709) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7168658) q[1];
sx q[1];
rz(-1.9391372) q[1];
sx q[1];
rz(1.4348381) q[1];
x q[2];
rz(2.8374568) q[3];
sx q[3];
rz(-1.602293) q[3];
sx q[3];
rz(-0.57284063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62766084) q[2];
sx q[2];
rz(-1.7646503) q[2];
sx q[2];
rz(0.39109209) q[2];
rz(-0.62134653) q[3];
sx q[3];
rz(-2.4070599) q[3];
sx q[3];
rz(2.3656316) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67052996) q[0];
sx q[0];
rz(-3.0369018) q[0];
sx q[0];
rz(1.4290357) q[0];
rz(-2.1757226) q[1];
sx q[1];
rz(-1.254225) q[1];
sx q[1];
rz(-2.3557854) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87596506) q[0];
sx q[0];
rz(-2.4959491) q[0];
sx q[0];
rz(3.0558048) q[0];
rz(-pi) q[1];
rz(-2.4072717) q[2];
sx q[2];
rz(-1.87748) q[2];
sx q[2];
rz(-1.1994565) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8700669) q[1];
sx q[1];
rz(-1.1114239) q[1];
sx q[1];
rz(-1.4794631) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8604516) q[3];
sx q[3];
rz(-1.9928586) q[3];
sx q[3];
rz(-0.45988032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4449731) q[2];
sx q[2];
rz(-0.61903054) q[2];
sx q[2];
rz(-2.6444198) q[2];
rz(-2.2677126) q[3];
sx q[3];
rz(-1.0920352) q[3];
sx q[3];
rz(-1.9397651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9484321) q[0];
sx q[0];
rz(-0.30696294) q[0];
sx q[0];
rz(0.69751414) q[0];
rz(-2.7613617) q[1];
sx q[1];
rz(-2.0262599) q[1];
sx q[1];
rz(-3.0013705) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.396616) q[0];
sx q[0];
rz(-0.19752398) q[0];
sx q[0];
rz(-0.4032938) q[0];
x q[1];
rz(-0.082265286) q[2];
sx q[2];
rz(-2.6994355) q[2];
sx q[2];
rz(-1.7205659) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6517087) q[1];
sx q[1];
rz(-1.1715285) q[1];
sx q[1];
rz(0.7266161) q[1];
rz(0.41922064) q[3];
sx q[3];
rz(-0.44525075) q[3];
sx q[3];
rz(-0.78263301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.85072881) q[2];
sx q[2];
rz(-1.9242492) q[2];
sx q[2];
rz(-0.49120894) q[2];
rz(-0.20719191) q[3];
sx q[3];
rz(-0.62316337) q[3];
sx q[3];
rz(-1.4108968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9687013) q[0];
sx q[0];
rz(-1.7270813) q[0];
sx q[0];
rz(0.15039314) q[0];
rz(0.70101678) q[1];
sx q[1];
rz(-1.0944518) q[1];
sx q[1];
rz(1.6640123) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9212338) q[0];
sx q[0];
rz(-1.6920751) q[0];
sx q[0];
rz(-0.63684271) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16123743) q[2];
sx q[2];
rz(-1.6856226) q[2];
sx q[2];
rz(1.6063362) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36819283) q[1];
sx q[1];
rz(-0.18562323) q[1];
sx q[1];
rz(-2.5297861) q[1];
rz(-pi) q[2];
rz(1.2712237) q[3];
sx q[3];
rz(-2.3837187) q[3];
sx q[3];
rz(-1.8488499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.574934) q[2];
sx q[2];
rz(-0.40859544) q[2];
sx q[2];
rz(-1.6179786) q[2];
rz(-0.59897113) q[3];
sx q[3];
rz(-0.50061575) q[3];
sx q[3];
rz(-2.9536501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4481675) q[0];
sx q[0];
rz(-2.7126815) q[0];
sx q[0];
rz(1.6397788) q[0];
rz(-2.2046454) q[1];
sx q[1];
rz(-2.1138771) q[1];
sx q[1];
rz(1.017259) q[1];
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
x q[1];
rz(-1.6566562) q[2];
sx q[2];
rz(-1.4659662) q[2];
sx q[2];
rz(1.3121999) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.81472936) q[1];
sx q[1];
rz(-1.5376602) q[1];
sx q[1];
rz(-1.3087981) q[1];
rz(-pi) q[2];
rz(-1.7726835) q[3];
sx q[3];
rz(-1.9924148) q[3];
sx q[3];
rz(3.0516171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0246058) q[2];
sx q[2];
rz(-0.6822497) q[2];
sx q[2];
rz(-1.619722) q[2];
rz(1.4874602) q[3];
sx q[3];
rz(-1.6493075) q[3];
sx q[3];
rz(1.3967995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37987729) q[0];
sx q[0];
rz(-1.7424255) q[0];
sx q[0];
rz(-2.4095834) q[0];
rz(-1.4749745) q[1];
sx q[1];
rz(-1.2154308) q[1];
sx q[1];
rz(-2.348127) q[1];
rz(0.99967069) q[2];
sx q[2];
rz(-0.66304211) q[2];
sx q[2];
rz(-2.1950051) q[2];
rz(-2.8216437) q[3];
sx q[3];
rz(-1.8599763) q[3];
sx q[3];
rz(-1.57747) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.30224213) q[0];
sx q[0];
rz(-1.2737162) q[0];
sx q[0];
rz(-2.1923375) q[0];
rz(-1.0957837) q[1];
sx q[1];
rz(-0.97691184) q[1];
sx q[1];
rz(-1.3499324) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95800864) q[0];
sx q[0];
rz(-1.4877948) q[0];
sx q[0];
rz(0.50358332) q[0];
rz(0.40387965) q[2];
sx q[2];
rz(-2.7821861) q[2];
sx q[2];
rz(-3.0702116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3826582) q[1];
sx q[1];
rz(-1.9449502) q[1];
sx q[1];
rz(-1.4976682) q[1];
x q[2];
rz(-1.1383406) q[3];
sx q[3];
rz(-1.3210332) q[3];
sx q[3];
rz(-0.43507155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.533941) q[2];
sx q[2];
rz(-2.0484296) q[2];
sx q[2];
rz(2.9489813) q[2];
rz(1.2827778) q[3];
sx q[3];
rz(-1.1056113) q[3];
sx q[3];
rz(-2.5530596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31805661) q[0];
sx q[0];
rz(-0.41724351) q[0];
sx q[0];
rz(2.4349924) q[0];
rz(-2.508714) q[1];
sx q[1];
rz(-1.0989847) q[1];
sx q[1];
rz(-0.28876567) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97917999) q[0];
sx q[0];
rz(-1.5683163) q[0];
sx q[0];
rz(-1.5697877) q[0];
rz(-1.7328506) q[2];
sx q[2];
rz(-1.7758435) q[2];
sx q[2];
rz(-1.0150725) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.56882492) q[1];
sx q[1];
rz(-1.8563885) q[1];
sx q[1];
rz(1.3495803) q[1];
rz(-2.3815126) q[3];
sx q[3];
rz(-1.8853429) q[3];
sx q[3];
rz(9.4903367e-05) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67549813) q[2];
sx q[2];
rz(-1.0885295) q[2];
sx q[2];
rz(1.4906073) q[2];
rz(0.7771107) q[3];
sx q[3];
rz(-1.4584352) q[3];
sx q[3];
rz(1.7818264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31558388) q[0];
sx q[0];
rz(-1.4588139) q[0];
sx q[0];
rz(2.7050731) q[0];
rz(-2.3468158) q[1];
sx q[1];
rz(-1.7710641) q[1];
sx q[1];
rz(-0.93903843) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7815112) q[0];
sx q[0];
rz(-2.0716801) q[0];
sx q[0];
rz(-1.341218) q[0];
rz(0.76656966) q[2];
sx q[2];
rz(-1.9992644) q[2];
sx q[2];
rz(-0.049190532) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.66575501) q[1];
sx q[1];
rz(-1.8867954) q[1];
sx q[1];
rz(2.5412987) q[1];
rz(-pi) q[2];
rz(2.125199) q[3];
sx q[3];
rz(-1.3825584) q[3];
sx q[3];
rz(-2.7013626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.44094917) q[2];
sx q[2];
rz(-1.0627397) q[2];
sx q[2];
rz(-2.5062594) q[2];
rz(2.3144531) q[3];
sx q[3];
rz(-2.6878036) q[3];
sx q[3];
rz(-1.872983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24469911) q[0];
sx q[0];
rz(-3.0785955) q[0];
sx q[0];
rz(2.6196106) q[0];
rz(-2.7105647) q[1];
sx q[1];
rz(-1.5343752) q[1];
sx q[1];
rz(-0.65188754) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2289705) q[0];
sx q[0];
rz(-0.56601277) q[0];
sx q[0];
rz(2.6900351) q[0];
rz(2.2950933) q[2];
sx q[2];
rz(-1.1870459) q[2];
sx q[2];
rz(2.4946314) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.83623278) q[1];
sx q[1];
rz(-1.1916324) q[1];
sx q[1];
rz(3.1268478) q[1];
x q[2];
rz(2.4595991) q[3];
sx q[3];
rz(-2.3122283) q[3];
sx q[3];
rz(-3.1113868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.42644694) q[2];
sx q[2];
rz(-1.8603674) q[2];
sx q[2];
rz(0.78488266) q[2];
rz(0.52810413) q[3];
sx q[3];
rz(-1.2930861) q[3];
sx q[3];
rz(-1.8538792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2454979) q[0];
sx q[0];
rz(-2.446785) q[0];
sx q[0];
rz(0.78952638) q[0];
rz(-0.49742571) q[1];
sx q[1];
rz(-2.1010294) q[1];
sx q[1];
rz(-1.3015889) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042846366) q[0];
sx q[0];
rz(-1.3912956) q[0];
sx q[0];
rz(0.60447201) q[0];
rz(0.33904262) q[2];
sx q[2];
rz(-1.6280988) q[2];
sx q[2];
rz(2.6471241) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6005206) q[1];
sx q[1];
rz(-0.75319911) q[1];
sx q[1];
rz(2.888117) q[1];
rz(-pi) q[2];
rz(-1.6297518) q[3];
sx q[3];
rz(-2.7736933) q[3];
sx q[3];
rz(-2.2045362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5975981) q[2];
sx q[2];
rz(-1.5388637) q[2];
sx q[2];
rz(0.27628118) q[2];
rz(-0.9497408) q[3];
sx q[3];
rz(-0.26849982) q[3];
sx q[3];
rz(-0.55324078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62040579) q[0];
sx q[0];
rz(-0.036245417) q[0];
sx q[0];
rz(2.1976443) q[0];
rz(1.2983407) q[1];
sx q[1];
rz(-1.0669758) q[1];
sx q[1];
rz(-2.3640769) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1489747) q[0];
sx q[0];
rz(-1.0273522) q[0];
sx q[0];
rz(-1.987756) q[0];
rz(-pi) q[1];
rz(0.73410122) q[2];
sx q[2];
rz(-1.3742374) q[2];
sx q[2];
rz(0.71069709) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7168658) q[1];
sx q[1];
rz(-1.9391372) q[1];
sx q[1];
rz(-1.4348381) q[1];
rz(-2.8374568) q[3];
sx q[3];
rz(-1.602293) q[3];
sx q[3];
rz(-2.568752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5139318) q[2];
sx q[2];
rz(-1.3769423) q[2];
sx q[2];
rz(0.39109209) q[2];
rz(-2.5202461) q[3];
sx q[3];
rz(-0.73453271) q[3];
sx q[3];
rz(-0.77596107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.67052996) q[0];
sx q[0];
rz(-3.0369018) q[0];
sx q[0];
rz(1.712557) q[0];
rz(2.1757226) q[1];
sx q[1];
rz(-1.254225) q[1];
sx q[1];
rz(-0.78580725) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5153421) q[0];
sx q[0];
rz(-1.6223755) q[0];
sx q[0];
rz(-0.64387384) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73432095) q[2];
sx q[2];
rz(-1.2641126) q[2];
sx q[2];
rz(1.1994565) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.475226) q[1];
sx q[1];
rz(-2.6738648) q[1];
sx q[1];
rz(-0.18233129) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7034055) q[3];
sx q[3];
rz(-1.3071968) q[3];
sx q[3];
rz(2.1521642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69661951) q[2];
sx q[2];
rz(-2.5225621) q[2];
sx q[2];
rz(2.6444198) q[2];
rz(-0.87388006) q[3];
sx q[3];
rz(-1.0920352) q[3];
sx q[3];
rz(1.9397651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9484321) q[0];
sx q[0];
rz(-2.8346297) q[0];
sx q[0];
rz(-2.4440785) q[0];
rz(-2.7613617) q[1];
sx q[1];
rz(-2.0262599) q[1];
sx q[1];
rz(-3.0013705) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33454681) q[0];
sx q[0];
rz(-1.3893034) q[0];
sx q[0];
rz(-1.4924148) q[0];
rz(1.6096787) q[2];
sx q[2];
rz(-2.011353) q[2];
sx q[2];
rz(-1.8115385) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6483874) q[1];
sx q[1];
rz(-2.3304061) q[1];
sx q[1];
rz(0.56583515) q[1];
x q[2];
rz(-1.7626552) q[3];
sx q[3];
rz(-1.1664806) q[3];
sx q[3];
rz(1.9003001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2908638) q[2];
sx q[2];
rz(-1.2173434) q[2];
sx q[2];
rz(-0.49120894) q[2];
rz(2.9344007) q[3];
sx q[3];
rz(-0.62316337) q[3];
sx q[3];
rz(-1.4108968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-2.0471408) q[1];
sx q[1];
rz(-1.6640123) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22035881) q[0];
sx q[0];
rz(-1.6920751) q[0];
sx q[0];
rz(-2.5047499) q[0];
x q[1];
rz(2.9803552) q[2];
sx q[2];
rz(-1.6856226) q[2];
sx q[2];
rz(1.6063362) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7733998) q[1];
sx q[1];
rz(-2.9559694) q[1];
sx q[1];
rz(0.61180656) q[1];
x q[2];
rz(1.870369) q[3];
sx q[3];
rz(-0.75787395) q[3];
sx q[3];
rz(-1.8488499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.56665862) q[2];
sx q[2];
rz(-0.40859544) q[2];
sx q[2];
rz(1.6179786) q[2];
rz(2.5426215) q[3];
sx q[3];
rz(-0.50061575) q[3];
sx q[3];
rz(0.18794255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4481675) q[0];
sx q[0];
rz(-2.7126815) q[0];
sx q[0];
rz(-1.5018139) q[0];
rz(-2.2046454) q[1];
sx q[1];
rz(-2.1138771) q[1];
sx q[1];
rz(1.017259) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.992998) q[0];
sx q[0];
rz(-0.55363292) q[0];
sx q[0];
rz(-0.77405907) q[0];
rz(-pi) q[1];
rz(1.4849365) q[2];
sx q[2];
rz(-1.6756264) q[2];
sx q[2];
rz(-1.3121999) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3268633) q[1];
sx q[1];
rz(-1.5376602) q[1];
sx q[1];
rz(-1.3087981) q[1];
x q[2];
rz(0.42041619) q[3];
sx q[3];
rz(-0.46483332) q[3];
sx q[3];
rz(-2.5878588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1169869) q[2];
sx q[2];
rz(-0.6822497) q[2];
sx q[2];
rz(1.619722) q[2];
rz(-1.4874602) q[3];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37987729) q[0];
sx q[0];
rz(-1.3991671) q[0];
sx q[0];
rz(0.73200926) q[0];
rz(1.6666182) q[1];
sx q[1];
rz(-1.2154308) q[1];
sx q[1];
rz(-2.348127) q[1];
rz(2.141922) q[2];
sx q[2];
rz(-2.4785505) q[2];
sx q[2];
rz(0.94658755) q[2];
rz(-0.75763221) q[3];
sx q[3];
rz(-0.42790596) q[3];
sx q[3];
rz(-0.71732646) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

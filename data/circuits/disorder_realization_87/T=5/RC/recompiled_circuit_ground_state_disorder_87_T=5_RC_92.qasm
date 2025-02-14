OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4053722) q[0];
sx q[0];
rz(-0.020981941) q[0];
sx q[0];
rz(-1.1809281) q[0];
rz(2.6612072) q[1];
sx q[1];
rz(-1.1552224) q[1];
sx q[1];
rz(0.46673271) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6581527) q[0];
sx q[0];
rz(-2.3077093) q[0];
sx q[0];
rz(1.9574653) q[0];
rz(-2.4619815) q[2];
sx q[2];
rz(-0.70290297) q[2];
sx q[2];
rz(1.694569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8267104) q[1];
sx q[1];
rz(-0.2383543) q[1];
sx q[1];
rz(-1.5895859) q[1];
rz(-pi) q[2];
rz(-1.0084413) q[3];
sx q[3];
rz(-2.1690024) q[3];
sx q[3];
rz(-1.7266718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0537609) q[2];
sx q[2];
rz(-1.6873282) q[2];
sx q[2];
rz(-1.2254747) q[2];
rz(-0.37442225) q[3];
sx q[3];
rz(-1.1332847) q[3];
sx q[3];
rz(0.55330127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-0.18749172) q[0];
sx q[0];
rz(-0.88045374) q[0];
sx q[0];
rz(-2.6075897) q[0];
rz(-0.39132896) q[1];
sx q[1];
rz(-2.6983039) q[1];
sx q[1];
rz(0.88723976) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92838642) q[0];
sx q[0];
rz(-1.22661) q[0];
sx q[0];
rz(-1.8977988) q[0];
rz(2.5975304) q[2];
sx q[2];
rz(-0.83457046) q[2];
sx q[2];
rz(2.6092333) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56480592) q[1];
sx q[1];
rz(-2.8999355) q[1];
sx q[1];
rz(2.5741626) q[1];
rz(-1.0530472) q[3];
sx q[3];
rz(-1.2710921) q[3];
sx q[3];
rz(-0.24417711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5947764) q[2];
sx q[2];
rz(-1.8211326) q[2];
sx q[2];
rz(-0.36396626) q[2];
rz(1.724203) q[3];
sx q[3];
rz(-2.851749) q[3];
sx q[3];
rz(2.3072306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2367547) q[0];
sx q[0];
rz(-0.098243864) q[0];
sx q[0];
rz(2.8114317) q[0];
rz(2.5579021) q[1];
sx q[1];
rz(-0.81283641) q[1];
sx q[1];
rz(-2.6469753) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0000293) q[0];
sx q[0];
rz(-2.7589679) q[0];
sx q[0];
rz(-2.3429246) q[0];
x q[1];
rz(-2.885899) q[2];
sx q[2];
rz(-2.0478559) q[2];
sx q[2];
rz(-3.0050803) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6223294) q[1];
sx q[1];
rz(-2.0262782) q[1];
sx q[1];
rz(-2.141327) q[1];
x q[2];
rz(-0.98683896) q[3];
sx q[3];
rz(-1.47746) q[3];
sx q[3];
rz(0.10658857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3459449) q[2];
sx q[2];
rz(-0.66850418) q[2];
sx q[2];
rz(0.12269679) q[2];
rz(3.0986541) q[3];
sx q[3];
rz(-1.0791082) q[3];
sx q[3];
rz(0.19449657) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73404679) q[0];
sx q[0];
rz(-0.50739822) q[0];
sx q[0];
rz(-2.2518482) q[0];
rz(2.6906158) q[1];
sx q[1];
rz(-2.273592) q[1];
sx q[1];
rz(0.45423347) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.736859) q[0];
sx q[0];
rz(-1.6772207) q[0];
sx q[0];
rz(2.7989796) q[0];
rz(-pi) q[1];
rz(0.003963917) q[2];
sx q[2];
rz(-2.7142453) q[2];
sx q[2];
rz(-1.3766118) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3754546) q[1];
sx q[1];
rz(-1.8667867) q[1];
sx q[1];
rz(2.3691404) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7866335) q[3];
sx q[3];
rz(-1.4982274) q[3];
sx q[3];
rz(-1.5458969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1526327) q[2];
sx q[2];
rz(-0.6344513) q[2];
sx q[2];
rz(-2.3839942) q[2];
rz(-2.9518413) q[3];
sx q[3];
rz(-1.8228143) q[3];
sx q[3];
rz(-0.54722133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9309288) q[0];
sx q[0];
rz(-2.3543816) q[0];
sx q[0];
rz(-1.8481365) q[0];
rz(2.5594607) q[1];
sx q[1];
rz(-1.7030092) q[1];
sx q[1];
rz(-0.17939803) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.005527786) q[0];
sx q[0];
rz(-1.3910157) q[0];
sx q[0];
rz(-1.3639569) q[0];
rz(-1.8151692) q[2];
sx q[2];
rz(-0.91460157) q[2];
sx q[2];
rz(0.49733053) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7475954) q[1];
sx q[1];
rz(-2.4863805) q[1];
sx q[1];
rz(-1.2607695) q[1];
rz(-pi) q[2];
rz(-2.6396181) q[3];
sx q[3];
rz(-2.2036512) q[3];
sx q[3];
rz(2.8090734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.77357117) q[2];
sx q[2];
rz(-2.7715235) q[2];
sx q[2];
rz(0.5698815) q[2];
rz(-0.44024769) q[3];
sx q[3];
rz(-1.8972242) q[3];
sx q[3];
rz(-0.35278916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2920947) q[0];
sx q[0];
rz(-2.8773913) q[0];
sx q[0];
rz(2.0000892) q[0];
rz(1.796272) q[1];
sx q[1];
rz(-2.1514386) q[1];
sx q[1];
rz(-0.17123953) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2889338) q[0];
sx q[0];
rz(-1.4698514) q[0];
sx q[0];
rz(0.48152083) q[0];
rz(-1.0467347) q[2];
sx q[2];
rz(-1.9643133) q[2];
sx q[2];
rz(1.6167906) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6936612) q[1];
sx q[1];
rz(-1.0180001) q[1];
sx q[1];
rz(1.3752244) q[1];
rz(-2.1275131) q[3];
sx q[3];
rz(-1.5814335) q[3];
sx q[3];
rz(2.6680846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77219999) q[2];
sx q[2];
rz(-2.2056396) q[2];
sx q[2];
rz(3.0041223) q[2];
rz(1.8933659) q[3];
sx q[3];
rz(-0.56709254) q[3];
sx q[3];
rz(1.9198157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0565979) q[0];
sx q[0];
rz(-2.5039112) q[0];
sx q[0];
rz(-1.6531264) q[0];
rz(0.78530637) q[1];
sx q[1];
rz(-1.7410802) q[1];
sx q[1];
rz(-1.8280425) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97264719) q[0];
sx q[0];
rz(-0.15014507) q[0];
sx q[0];
rz(-0.96359588) q[0];
rz(-pi) q[1];
rz(-1.684217) q[2];
sx q[2];
rz(-1.3807502) q[2];
sx q[2];
rz(-0.52670497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.32282) q[1];
sx q[1];
rz(-1.5427151) q[1];
sx q[1];
rz(3.0525467) q[1];
rz(0.87487674) q[3];
sx q[3];
rz(-1.531344) q[3];
sx q[3];
rz(-2.6764398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.58275783) q[2];
sx q[2];
rz(-1.02966) q[2];
sx q[2];
rz(-2.6666857) q[2];
rz(2.9389935) q[3];
sx q[3];
rz(-0.83297268) q[3];
sx q[3];
rz(-1.9995662) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1982034) q[0];
sx q[0];
rz(-3.0084963) q[0];
sx q[0];
rz(0.30580795) q[0];
rz(-0.43931475) q[1];
sx q[1];
rz(-0.82249928) q[1];
sx q[1];
rz(1.3154715) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00061247608) q[0];
sx q[0];
rz(-0.68353739) q[0];
sx q[0];
rz(-0.017869259) q[0];
x q[1];
rz(2.010263) q[2];
sx q[2];
rz(-0.60053289) q[2];
sx q[2];
rz(3.0366922) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5852192) q[1];
sx q[1];
rz(-1.4498158) q[1];
sx q[1];
rz(2.35384) q[1];
rz(-pi) q[2];
rz(-2.1602116) q[3];
sx q[3];
rz(-1.9187315) q[3];
sx q[3];
rz(-0.98435005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7266015) q[2];
sx q[2];
rz(-1.5433658) q[2];
sx q[2];
rz(-2.7169363) q[2];
rz(-3.1096544) q[3];
sx q[3];
rz(-1.5141124) q[3];
sx q[3];
rz(2.4922075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53751078) q[0];
sx q[0];
rz(-0.6518971) q[0];
sx q[0];
rz(-0.55484581) q[0];
rz(-0.26508731) q[1];
sx q[1];
rz(-1.4684497) q[1];
sx q[1];
rz(1.2010942) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6938615) q[0];
sx q[0];
rz(-0.70708067) q[0];
sx q[0];
rz(-2.5140544) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5910406) q[2];
sx q[2];
rz(-0.32062396) q[2];
sx q[2];
rz(0.52832802) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8256067) q[1];
sx q[1];
rz(-2.791643) q[1];
sx q[1];
rz(-3.0341329) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3917434) q[3];
sx q[3];
rz(-1.5886652) q[3];
sx q[3];
rz(-0.32040473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.889664) q[2];
sx q[2];
rz(-2.1519075) q[2];
sx q[2];
rz(-0.84452334) q[2];
rz(1.1318413) q[3];
sx q[3];
rz(-0.90513888) q[3];
sx q[3];
rz(-1.3593146) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97245401) q[0];
sx q[0];
rz(-2.7511981) q[0];
sx q[0];
rz(1.9688695) q[0];
rz(-2.4523465) q[1];
sx q[1];
rz(-1.0968364) q[1];
sx q[1];
rz(2.8776339) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.496752) q[0];
sx q[0];
rz(-0.94366108) q[0];
sx q[0];
rz(-2.1261313) q[0];
rz(0.51610871) q[2];
sx q[2];
rz(-0.91876578) q[2];
sx q[2];
rz(2.3626346) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.76494886) q[1];
sx q[1];
rz(-1.7564492) q[1];
sx q[1];
rz(-1.8913881) q[1];
rz(-1.7329151) q[3];
sx q[3];
rz(-0.65026426) q[3];
sx q[3];
rz(-1.4635603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3802203) q[2];
sx q[2];
rz(-1.994588) q[2];
sx q[2];
rz(-1.1205565) q[2];
rz(1.5348966) q[3];
sx q[3];
rz(-2.1963019) q[3];
sx q[3];
rz(-1.5103316) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0705538) q[0];
sx q[0];
rz(-0.65503913) q[0];
sx q[0];
rz(-0.10792637) q[0];
rz(0.72036605) q[1];
sx q[1];
rz(-1.9279059) q[1];
sx q[1];
rz(2.7005213) q[1];
rz(-2.5369562) q[2];
sx q[2];
rz(-2.2199134) q[2];
sx q[2];
rz(1.8668957) q[2];
rz(-1.646044) q[3];
sx q[3];
rz(-0.9552707) q[3];
sx q[3];
rz(0.09494119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

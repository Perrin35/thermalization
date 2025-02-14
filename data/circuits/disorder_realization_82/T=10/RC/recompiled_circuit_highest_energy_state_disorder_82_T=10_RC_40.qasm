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
rz(0.0042301099) q[0];
sx q[0];
rz(-0.63679305) q[0];
sx q[0];
rz(-1.1871583) q[0];
rz(0.98183331) q[1];
sx q[1];
rz(-0.47458664) q[1];
sx q[1];
rz(2.2003953) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1008074) q[0];
sx q[0];
rz(-1.5724764) q[0];
sx q[0];
rz(2.6462862) q[0];
rz(-pi) q[1];
x q[1];
rz(1.773052) q[2];
sx q[2];
rz(-1.6967684) q[2];
sx q[2];
rz(1.6304315) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.57454771) q[1];
sx q[1];
rz(-0.68463782) q[1];
sx q[1];
rz(2.5570832) q[1];
rz(-pi) q[2];
rz(-1.7205129) q[3];
sx q[3];
rz(-2.6805373) q[3];
sx q[3];
rz(1.0284916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5473951) q[2];
sx q[2];
rz(-1.5028468) q[2];
sx q[2];
rz(2.1217864) q[2];
rz(-2.8504573) q[3];
sx q[3];
rz(-0.19164339) q[3];
sx q[3];
rz(-1.8118793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-2.6137467) q[0];
sx q[0];
rz(-2.3044523) q[0];
sx q[0];
rz(-1.5035195) q[0];
rz(-1.9095406) q[1];
sx q[1];
rz(-2.3161395) q[1];
sx q[1];
rz(-1.253461) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4811752) q[0];
sx q[0];
rz(-1.345013) q[0];
sx q[0];
rz(1.9893454) q[0];
rz(-0.49053662) q[2];
sx q[2];
rz(-0.56714688) q[2];
sx q[2];
rz(0.74037742) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8294222) q[1];
sx q[1];
rz(-0.24917425) q[1];
sx q[1];
rz(-2.5214358) q[1];
rz(-pi) q[2];
rz(-2.3971629) q[3];
sx q[3];
rz(-1.9081313) q[3];
sx q[3];
rz(2.6312089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.85679179) q[2];
sx q[2];
rz(-0.54232001) q[2];
sx q[2];
rz(-0.33033237) q[2];
rz(3.0387943) q[3];
sx q[3];
rz(-1.8967352) q[3];
sx q[3];
rz(0.61906329) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5264346) q[0];
sx q[0];
rz(-1.1410843) q[0];
sx q[0];
rz(-2.0779628) q[0];
rz(-3.0377153) q[1];
sx q[1];
rz(-2.3492298) q[1];
sx q[1];
rz(1.1054543) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5828667) q[0];
sx q[0];
rz(-1.5387282) q[0];
sx q[0];
rz(1.523287) q[0];
x q[1];
rz(-2.8145829) q[2];
sx q[2];
rz(-0.78377027) q[2];
sx q[2];
rz(1.4732337) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3893095) q[1];
sx q[1];
rz(-1.0573309) q[1];
sx q[1];
rz(3.0016243) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1298769) q[3];
sx q[3];
rz(-1.5289617) q[3];
sx q[3];
rz(-1.5277629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0217454) q[2];
sx q[2];
rz(-0.80588078) q[2];
sx q[2];
rz(-0.22551192) q[2];
rz(1.3269904) q[3];
sx q[3];
rz(-0.83695379) q[3];
sx q[3];
rz(-2.5007611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91318146) q[0];
sx q[0];
rz(-1.5786194) q[0];
sx q[0];
rz(-0.6382424) q[0];
rz(2.0300716) q[1];
sx q[1];
rz(-0.94739729) q[1];
sx q[1];
rz(2.9543455) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53786862) q[0];
sx q[0];
rz(-2.4477673) q[0];
sx q[0];
rz(1.1129473) q[0];
x q[1];
rz(0.063129877) q[2];
sx q[2];
rz(-2.3597096) q[2];
sx q[2];
rz(-0.031459122) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3447002) q[1];
sx q[1];
rz(-0.53965599) q[1];
sx q[1];
rz(2.9873288) q[1];
x q[2];
rz(-1.4832014) q[3];
sx q[3];
rz(-2.2822652) q[3];
sx q[3];
rz(-0.30184823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.036666544) q[2];
sx q[2];
rz(-1.6380402) q[2];
sx q[2];
rz(0.46910134) q[2];
rz(-0.30096287) q[3];
sx q[3];
rz(-2.8807237) q[3];
sx q[3];
rz(-1.4891967) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24333532) q[0];
sx q[0];
rz(-1.5168334) q[0];
sx q[0];
rz(-0.54018706) q[0];
rz(0.46145269) q[1];
sx q[1];
rz(-1.5216454) q[1];
sx q[1];
rz(0.25925055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.777939) q[0];
sx q[0];
rz(-0.41746062) q[0];
sx q[0];
rz(2.4221942) q[0];
rz(-pi) q[1];
rz(2.0820145) q[2];
sx q[2];
rz(-1.9661742) q[2];
sx q[2];
rz(-1.2653654) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1446602) q[1];
sx q[1];
rz(-1.0029666) q[1];
sx q[1];
rz(2.9564942) q[1];
rz(0.26210014) q[3];
sx q[3];
rz(-1.8912695) q[3];
sx q[3];
rz(2.4985034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3429012) q[2];
sx q[2];
rz(-2.4553802) q[2];
sx q[2];
rz(2.2314824) q[2];
rz(2.5765007) q[3];
sx q[3];
rz(-1.3684973) q[3];
sx q[3];
rz(-2.4250987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27123505) q[0];
sx q[0];
rz(-0.39327455) q[0];
sx q[0];
rz(-1.2275335) q[0];
rz(0.50499376) q[1];
sx q[1];
rz(-1.0161437) q[1];
sx q[1];
rz(-1.203677) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6617463) q[0];
sx q[0];
rz(-1.8176314) q[0];
sx q[0];
rz(-1.6936452) q[0];
rz(-pi) q[1];
rz(2.2199322) q[2];
sx q[2];
rz(-1.9699012) q[2];
sx q[2];
rz(-1.468935) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.54654361) q[1];
sx q[1];
rz(-0.96009925) q[1];
sx q[1];
rz(2.2799019) q[1];
rz(-pi) q[2];
rz(-2.876128) q[3];
sx q[3];
rz(-2.0164312) q[3];
sx q[3];
rz(-1.3973941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0248854) q[2];
sx q[2];
rz(-2.1163157) q[2];
sx q[2];
rz(2.1825979) q[2];
rz(-0.2995019) q[3];
sx q[3];
rz(-3.0193269) q[3];
sx q[3];
rz(1.449077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2285948) q[0];
sx q[0];
rz(-1.2499502) q[0];
sx q[0];
rz(2.6477497) q[0];
rz(-0.85810581) q[1];
sx q[1];
rz(-1.162642) q[1];
sx q[1];
rz(-1.3546622) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6436967) q[0];
sx q[0];
rz(-0.61031155) q[0];
sx q[0];
rz(-1.1417139) q[0];
rz(2.0628711) q[2];
sx q[2];
rz(-1.4495696) q[2];
sx q[2];
rz(1.6693142) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1272116) q[1];
sx q[1];
rz(-0.76103044) q[1];
sx q[1];
rz(-0.10199031) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94791651) q[3];
sx q[3];
rz(-0.77548945) q[3];
sx q[3];
rz(2.754107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.4570423) q[2];
sx q[2];
rz(-2.5164618) q[2];
sx q[2];
rz(0.54614145) q[2];
rz(-0.1977194) q[3];
sx q[3];
rz(-1.0309912) q[3];
sx q[3];
rz(-2.8945967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
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
rz(-1.2208743) q[0];
sx q[0];
rz(-1.9863167) q[0];
sx q[0];
rz(1.1065296) q[0];
rz(2.5571892) q[1];
sx q[1];
rz(-1.9677275) q[1];
sx q[1];
rz(-1.0027142) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1258216) q[0];
sx q[0];
rz(-0.46575156) q[0];
sx q[0];
rz(1.2076103) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8991088) q[2];
sx q[2];
rz(-1.8576745) q[2];
sx q[2];
rz(0.74082021) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7103191) q[1];
sx q[1];
rz(-0.51352507) q[1];
sx q[1];
rz(3.0679614) q[1];
rz(-pi) q[2];
rz(0.61912304) q[3];
sx q[3];
rz(-1.2803708) q[3];
sx q[3];
rz(2.8324236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.57548412) q[2];
sx q[2];
rz(-0.9245975) q[2];
sx q[2];
rz(1.3893348) q[2];
rz(-0.53650457) q[3];
sx q[3];
rz(-1.5080695) q[3];
sx q[3];
rz(-0.84684053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.093512) q[0];
sx q[0];
rz(-0.41963136) q[0];
sx q[0];
rz(2.804948) q[0];
rz(-0.078314217) q[1];
sx q[1];
rz(-1.2093733) q[1];
sx q[1];
rz(-1.1557109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1378095) q[0];
sx q[0];
rz(-2.3562618) q[0];
sx q[0];
rz(-2.0565874) q[0];
x q[1];
rz(1.7215426) q[2];
sx q[2];
rz(-0.41778247) q[2];
sx q[2];
rz(2.5007375) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3382321) q[1];
sx q[1];
rz(-2.6322717) q[1];
sx q[1];
rz(1.2912911) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31678699) q[3];
sx q[3];
rz(-1.5996577) q[3];
sx q[3];
rz(-0.76735088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.081850514) q[2];
sx q[2];
rz(-1.0826702) q[2];
sx q[2];
rz(0.82320631) q[2];
rz(-2.1246223) q[3];
sx q[3];
rz(-2.7401994) q[3];
sx q[3];
rz(3.1034191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80242771) q[0];
sx q[0];
rz(-0.22295727) q[0];
sx q[0];
rz(1.4270576) q[0];
rz(-1.4786221) q[1];
sx q[1];
rz(-2.069811) q[1];
sx q[1];
rz(-2.2907385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9768391) q[0];
sx q[0];
rz(-2.0434336) q[0];
sx q[0];
rz(0.12906277) q[0];
rz(-2.1673604) q[2];
sx q[2];
rz(-2.9880004) q[2];
sx q[2];
rz(0.099791137) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5580299) q[1];
sx q[1];
rz(-2.0940015) q[1];
sx q[1];
rz(-1.5091404) q[1];
rz(1.6353278) q[3];
sx q[3];
rz(-1.7108265) q[3];
sx q[3];
rz(-0.16209312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6829546) q[2];
sx q[2];
rz(-1.3502324) q[2];
sx q[2];
rz(0.82009912) q[2];
rz(-1.594515) q[3];
sx q[3];
rz(-1.7087414) q[3];
sx q[3];
rz(1.9471751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.9783258) q[0];
sx q[0];
rz(-1.6303202) q[0];
sx q[0];
rz(-1.6250961) q[0];
rz(-0.77240472) q[1];
sx q[1];
rz(-2.518387) q[1];
sx q[1];
rz(1.0060681) q[1];
rz(-1.3428691) q[2];
sx q[2];
rz(-2.2493636) q[2];
sx q[2];
rz(2.0668088) q[2];
rz(2.6670868) q[3];
sx q[3];
rz(-1.9941327) q[3];
sx q[3];
rz(1.7607104) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

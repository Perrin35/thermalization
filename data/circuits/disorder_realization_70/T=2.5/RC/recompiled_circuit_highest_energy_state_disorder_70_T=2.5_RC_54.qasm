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
rz(-0.47082585) q[0];
sx q[0];
rz(3.5916632) q[0];
sx q[0];
rz(9.8150742) q[0];
rz(0.14416873) q[1];
sx q[1];
rz(-1.498797) q[1];
sx q[1];
rz(2.1013451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9594355) q[0];
sx q[0];
rz(-2.6819026) q[0];
sx q[0];
rz(-2.2612786) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7261271) q[2];
sx q[2];
rz(-1.9305781) q[2];
sx q[2];
rz(1.8632165) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92683243) q[1];
sx q[1];
rz(-2.4364987) q[1];
sx q[1];
rz(1.2486723) q[1];
rz(-pi) q[2];
rz(-0.066180996) q[3];
sx q[3];
rz(-0.77941637) q[3];
sx q[3];
rz(-2.4090637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1990004) q[2];
sx q[2];
rz(-0.60574836) q[2];
sx q[2];
rz(1.5462297) q[2];
rz(0.47510251) q[3];
sx q[3];
rz(-2.4964156) q[3];
sx q[3];
rz(1.9861541) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5763181) q[0];
sx q[0];
rz(-0.98868889) q[0];
sx q[0];
rz(-2.7069672) q[0];
rz(-2.077153) q[1];
sx q[1];
rz(-0.39355215) q[1];
sx q[1];
rz(-2.2501066) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1284244) q[0];
sx q[0];
rz(-1.1718318) q[0];
sx q[0];
rz(1.1815039) q[0];
rz(-pi) q[1];
rz(0.72618809) q[2];
sx q[2];
rz(-1.9627769) q[2];
sx q[2];
rz(1.7418944) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.52917186) q[1];
sx q[1];
rz(-1.1193706) q[1];
sx q[1];
rz(2.5542817) q[1];
rz(-pi) q[2];
rz(3.1117571) q[3];
sx q[3];
rz(-2.2298457) q[3];
sx q[3];
rz(-1.4277136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.50856227) q[2];
sx q[2];
rz(-0.56266251) q[2];
sx q[2];
rz(-2.0821345) q[2];
rz(1.9153473) q[3];
sx q[3];
rz(-1.6112593) q[3];
sx q[3];
rz(-0.1196158) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0934963) q[0];
sx q[0];
rz(-2.7432848) q[0];
sx q[0];
rz(2.4523729) q[0];
rz(-0.41220176) q[1];
sx q[1];
rz(-1.1176132) q[1];
sx q[1];
rz(-1.9788007) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0242694) q[0];
sx q[0];
rz(-1.4179967) q[0];
sx q[0];
rz(1.8788179) q[0];
rz(-0.96647977) q[2];
sx q[2];
rz(-0.89077026) q[2];
sx q[2];
rz(2.4228408) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7603569) q[1];
sx q[1];
rz(-0.79264005) q[1];
sx q[1];
rz(2.4355678) q[1];
x q[2];
rz(2.338004) q[3];
sx q[3];
rz(-2.7492964) q[3];
sx q[3];
rz(-2.4367849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73843655) q[2];
sx q[2];
rz(-1.7313892) q[2];
sx q[2];
rz(-0.23846826) q[2];
rz(-0.23369914) q[3];
sx q[3];
rz(-0.77478474) q[3];
sx q[3];
rz(0.15271798) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9018263) q[0];
sx q[0];
rz(-0.16401839) q[0];
sx q[0];
rz(2.8357586) q[0];
rz(0.26890525) q[1];
sx q[1];
rz(-0.79335672) q[1];
sx q[1];
rz(2.7893524) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21111317) q[0];
sx q[0];
rz(-1.6105299) q[0];
sx q[0];
rz(1.7046403) q[0];
rz(0.9446804) q[2];
sx q[2];
rz(-1.7292103) q[2];
sx q[2];
rz(0.84412727) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4388386) q[1];
sx q[1];
rz(-1.1795591) q[1];
sx q[1];
rz(-0.66968285) q[1];
rz(-pi) q[2];
rz(0.24573487) q[3];
sx q[3];
rz(-0.80348368) q[3];
sx q[3];
rz(1.5218228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8504146) q[2];
sx q[2];
rz(-1.6213657) q[2];
sx q[2];
rz(2.9746383) q[2];
rz(0.11635612) q[3];
sx q[3];
rz(-0.51893026) q[3];
sx q[3];
rz(1.8714582) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6312234) q[0];
sx q[0];
rz(-0.33741697) q[0];
sx q[0];
rz(1.1153197) q[0];
rz(1.910031) q[1];
sx q[1];
rz(-1.7111338) q[1];
sx q[1];
rz(0.11428741) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5530722) q[0];
sx q[0];
rz(-0.77295303) q[0];
sx q[0];
rz(0.78307481) q[0];
rz(2.8631526) q[2];
sx q[2];
rz(-1.8634708) q[2];
sx q[2];
rz(-3.1029683) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.85897) q[1];
sx q[1];
rz(-2.3630574) q[1];
sx q[1];
rz(-2.5384739) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.435569) q[3];
sx q[3];
rz(-1.5750275) q[3];
sx q[3];
rz(-2.8805681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.085122434) q[2];
sx q[2];
rz(-1.2135999) q[2];
sx q[2];
rz(-1.4465793) q[2];
rz(-0.96212402) q[3];
sx q[3];
rz(-0.4947997) q[3];
sx q[3];
rz(-1.3303293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
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
rz(2.1591448) q[0];
sx q[0];
rz(-1.5626937) q[0];
sx q[0];
rz(-0.53318095) q[0];
rz(1.5348596) q[1];
sx q[1];
rz(-0.77203647) q[1];
sx q[1];
rz(0.74347043) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3822298) q[0];
sx q[0];
rz(-2.8385525) q[0];
sx q[0];
rz(-2.8735301) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1986046) q[2];
sx q[2];
rz(-2.5649338) q[2];
sx q[2];
rz(-1.8221591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.20821321) q[1];
sx q[1];
rz(-0.68720308) q[1];
sx q[1];
rz(2.739868) q[1];
x q[2];
rz(-2.0251158) q[3];
sx q[3];
rz(-2.4719596) q[3];
sx q[3];
rz(2.7463934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47034904) q[2];
sx q[2];
rz(-0.84731421) q[2];
sx q[2];
rz(2.642855) q[2];
rz(2.4399452) q[3];
sx q[3];
rz(-0.68877733) q[3];
sx q[3];
rz(2.4832895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43211234) q[0];
sx q[0];
rz(-1.226959) q[0];
sx q[0];
rz(0.099932583) q[0];
rz(-1.8607032) q[1];
sx q[1];
rz(-0.51191267) q[1];
sx q[1];
rz(0.73922149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59935887) q[0];
sx q[0];
rz(-1.1809491) q[0];
sx q[0];
rz(2.705036) q[0];
rz(-pi) q[1];
rz(1.141369) q[2];
sx q[2];
rz(-1.543056) q[2];
sx q[2];
rz(0.77125473) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9770551) q[1];
sx q[1];
rz(-2.4083369) q[1];
sx q[1];
rz(-0.96340839) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80436665) q[3];
sx q[3];
rz(-1.5782246) q[3];
sx q[3];
rz(2.4480889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7382536) q[2];
sx q[2];
rz(-0.6311987) q[2];
sx q[2];
rz(-1.5367907) q[2];
rz(-1.712045) q[3];
sx q[3];
rz(-1.6481383) q[3];
sx q[3];
rz(-2.3454989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0159863) q[0];
sx q[0];
rz(-2.7637389) q[0];
sx q[0];
rz(1.1146389) q[0];
rz(0.58427018) q[1];
sx q[1];
rz(-0.92002267) q[1];
sx q[1];
rz(0.11411962) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18543359) q[0];
sx q[0];
rz(-1.8562278) q[0];
sx q[0];
rz(2.5117842) q[0];
rz(-pi) q[1];
rz(2.3405206) q[2];
sx q[2];
rz(-1.1566887) q[2];
sx q[2];
rz(1.0528477) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89592592) q[1];
sx q[1];
rz(-1.1628431) q[1];
sx q[1];
rz(2.308564) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8382254) q[3];
sx q[3];
rz(-1.4576313) q[3];
sx q[3];
rz(2.6306977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4383661) q[2];
sx q[2];
rz(-2.4195636) q[2];
sx q[2];
rz(-2.2528193) q[2];
rz(1.599954) q[3];
sx q[3];
rz(-1.9346574) q[3];
sx q[3];
rz(2.3204939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0515161) q[0];
sx q[0];
rz(-0.30802825) q[0];
sx q[0];
rz(-0.28710452) q[0];
rz(-1.9376532) q[1];
sx q[1];
rz(-0.86910373) q[1];
sx q[1];
rz(-1.513011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6490895) q[0];
sx q[0];
rz(-1.2375298) q[0];
sx q[0];
rz(1.4901584) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92219917) q[2];
sx q[2];
rz(-1.9364127) q[2];
sx q[2];
rz(-2.7067685) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2681458) q[1];
sx q[1];
rz(-1.526274) q[1];
sx q[1];
rz(0.73486199) q[1];
x q[2];
rz(1.318526) q[3];
sx q[3];
rz(-3.1384094) q[3];
sx q[3];
rz(-2.0306272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5055351) q[2];
sx q[2];
rz(-1.1670185) q[2];
sx q[2];
rz(0.30910811) q[2];
rz(1.2067893) q[3];
sx q[3];
rz(-2.759582) q[3];
sx q[3];
rz(2.8822854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1768271) q[0];
sx q[0];
rz(-2.4137156) q[0];
sx q[0];
rz(0.65560174) q[0];
rz(-2.5128095) q[1];
sx q[1];
rz(-1.7318334) q[1];
sx q[1];
rz(-1.7579196) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9533328) q[0];
sx q[0];
rz(-1.6142705) q[0];
sx q[0];
rz(-0.047634634) q[0];
x q[1];
rz(0.45510095) q[2];
sx q[2];
rz(-2.2544207) q[2];
sx q[2];
rz(-0.065082642) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68798472) q[1];
sx q[1];
rz(-1.5888643) q[1];
sx q[1];
rz(-0.14931909) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8877533) q[3];
sx q[3];
rz(-0.51341265) q[3];
sx q[3];
rz(1.1644582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3722966) q[2];
sx q[2];
rz(-1.5089704) q[2];
sx q[2];
rz(2.8995635) q[2];
rz(0.58487839) q[3];
sx q[3];
rz(-0.67348981) q[3];
sx q[3];
rz(-2.6175595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181716) q[0];
sx q[0];
rz(-0.62234288) q[0];
sx q[0];
rz(-0.55707669) q[0];
rz(0.88049018) q[1];
sx q[1];
rz(-1.7387895) q[1];
sx q[1];
rz(-2.1008076) q[1];
rz(-1.3646094) q[2];
sx q[2];
rz(-1.577081) q[2];
sx q[2];
rz(-1.1834363) q[2];
rz(0.09481341) q[3];
sx q[3];
rz(-2.1349813) q[3];
sx q[3];
rz(2.3053942) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

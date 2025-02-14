OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.40795657) q[0];
sx q[0];
rz(3.3269296) q[0];
sx q[0];
rz(11.067631) q[0];
rz(2.9474131) q[1];
sx q[1];
rz(-2.7979538) q[1];
sx q[1];
rz(-2.763881) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35750264) q[0];
sx q[0];
rz(-1.329485) q[0];
sx q[0];
rz(1.3426128) q[0];
rz(-1.1324563) q[2];
sx q[2];
rz(-1.9563663) q[2];
sx q[2];
rz(1.7913246) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6785665) q[1];
sx q[1];
rz(-0.65788473) q[1];
sx q[1];
rz(-3.0856864) q[1];
rz(1.6781647) q[3];
sx q[3];
rz(-2.8753548) q[3];
sx q[3];
rz(-3.0913946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7517884) q[2];
sx q[2];
rz(-1.7916388) q[2];
sx q[2];
rz(-2.8765615) q[2];
rz(2.7122279) q[3];
sx q[3];
rz(-1.2484442) q[3];
sx q[3];
rz(-0.00092367729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4877743) q[0];
sx q[0];
rz(-0.14350292) q[0];
sx q[0];
rz(-1.300746) q[0];
rz(-3.0744413) q[1];
sx q[1];
rz(-2.2323699) q[1];
sx q[1];
rz(-0.86004177) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9165827) q[0];
sx q[0];
rz(-0.89677484) q[0];
sx q[0];
rz(-0.81153481) q[0];
rz(0.76171909) q[2];
sx q[2];
rz(-2.2713619) q[2];
sx q[2];
rz(-0.127244) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7262267) q[1];
sx q[1];
rz(-0.72390717) q[1];
sx q[1];
rz(2.2103146) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8710526) q[3];
sx q[3];
rz(-1.3949864) q[3];
sx q[3];
rz(0.45775698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5106875) q[2];
sx q[2];
rz(-0.61669934) q[2];
sx q[2];
rz(2.5751233) q[2];
rz(-0.72930068) q[3];
sx q[3];
rz(-2.2972378) q[3];
sx q[3];
rz(2.9384889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1931964) q[0];
sx q[0];
rz(-2.0553135) q[0];
sx q[0];
rz(-0.36619827) q[0];
rz(1.5858448) q[1];
sx q[1];
rz(-1.0428753) q[1];
sx q[1];
rz(0.050447024) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4719066) q[0];
sx q[0];
rz(-1.6499053) q[0];
sx q[0];
rz(-3.0997653) q[0];
rz(-pi) q[1];
rz(-0.57951219) q[2];
sx q[2];
rz(-0.77634927) q[2];
sx q[2];
rz(0.001359847) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0512426) q[1];
sx q[1];
rz(-1.7900677) q[1];
sx q[1];
rz(1.3996814) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8837711) q[3];
sx q[3];
rz(-2.69776) q[3];
sx q[3];
rz(2.7601506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.067165) q[2];
sx q[2];
rz(-1.0600435) q[2];
sx q[2];
rz(1.3107497) q[2];
rz(1.7835167) q[3];
sx q[3];
rz(-0.57099968) q[3];
sx q[3];
rz(1.8324435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8686789) q[0];
sx q[0];
rz(-0.22628117) q[0];
sx q[0];
rz(-1.3847466) q[0];
rz(1.9028496) q[1];
sx q[1];
rz(-1.2307931) q[1];
sx q[1];
rz(1.967427) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0525111) q[0];
sx q[0];
rz(-1.3818372) q[0];
sx q[0];
rz(-2.004839) q[0];
rz(-0.27710813) q[2];
sx q[2];
rz(-1.3302251) q[2];
sx q[2];
rz(1.6024557) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3191022) q[1];
sx q[1];
rz(-0.61885364) q[1];
sx q[1];
rz(2.1560505) q[1];
rz(-2.2635095) q[3];
sx q[3];
rz(-1.9290079) q[3];
sx q[3];
rz(2.3830551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43526402) q[2];
sx q[2];
rz(-1.3202983) q[2];
sx q[2];
rz(0.038979385) q[2];
rz(2.4987761) q[3];
sx q[3];
rz(-1.0182074) q[3];
sx q[3];
rz(-0.62079287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83127999) q[0];
sx q[0];
rz(-2.1273002) q[0];
sx q[0];
rz(0.020462791) q[0];
rz(2.9488355) q[1];
sx q[1];
rz(-2.5242476) q[1];
sx q[1];
rz(1.9761168) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11949018) q[0];
sx q[0];
rz(-2.179232) q[0];
sx q[0];
rz(-1.0761989) q[0];
rz(-pi) q[1];
rz(1.1212249) q[2];
sx q[2];
rz(-0.94749852) q[2];
sx q[2];
rz(2.4433608) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7811657) q[1];
sx q[1];
rz(-1.5505704) q[1];
sx q[1];
rz(1.0670877) q[1];
x q[2];
rz(-1.7609414) q[3];
sx q[3];
rz(-1.645429) q[3];
sx q[3];
rz(0.41342218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8274902) q[2];
sx q[2];
rz(-0.9919439) q[2];
sx q[2];
rz(-2.6143383) q[2];
rz(-0.54316795) q[3];
sx q[3];
rz(-0.68734622) q[3];
sx q[3];
rz(-0.34272042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-1.0602033) q[0];
sx q[0];
rz(-2.9855766) q[0];
sx q[0];
rz(-0.40670893) q[0];
rz(1.1609062) q[1];
sx q[1];
rz(-0.91861594) q[1];
sx q[1];
rz(-2.1312174) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.062101) q[0];
sx q[0];
rz(-1.6915503) q[0];
sx q[0];
rz(-2.1202205) q[0];
x q[1];
rz(-0.40596227) q[2];
sx q[2];
rz(-2.7926707) q[2];
sx q[2];
rz(-0.59970784) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95376172) q[1];
sx q[1];
rz(-1.7825025) q[1];
sx q[1];
rz(-0.37734887) q[1];
rz(-1.3090735) q[3];
sx q[3];
rz(-0.78262586) q[3];
sx q[3];
rz(-1.4306376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0747718) q[2];
sx q[2];
rz(-1.2819042) q[2];
sx q[2];
rz(1.034896) q[2];
rz(1.3567989) q[3];
sx q[3];
rz(-0.59485888) q[3];
sx q[3];
rz(0.6234197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9676301) q[0];
sx q[0];
rz(-1.5246464) q[0];
sx q[0];
rz(0.3279283) q[0];
rz(-3.0294042) q[1];
sx q[1];
rz(-1.2022377) q[1];
sx q[1];
rz(0.97253886) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2043861) q[0];
sx q[0];
rz(-2.4582259) q[0];
sx q[0];
rz(2.1478189) q[0];
rz(2.1489086) q[2];
sx q[2];
rz(-0.51380605) q[2];
sx q[2];
rz(-0.93573278) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.211208) q[1];
sx q[1];
rz(-0.58961419) q[1];
sx q[1];
rz(-3.1077887) q[1];
x q[2];
rz(2.6554606) q[3];
sx q[3];
rz(-2.4555169) q[3];
sx q[3];
rz(-3.0239575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5358676) q[2];
sx q[2];
rz(-2.2650227) q[2];
sx q[2];
rz(2.4884339) q[2];
rz(-2.7086835) q[3];
sx q[3];
rz(-2.631729) q[3];
sx q[3];
rz(-2.9292817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9369649) q[0];
sx q[0];
rz(-0.84119868) q[0];
sx q[0];
rz(0.2247819) q[0];
rz(2.3460491) q[1];
sx q[1];
rz(-1.8276151) q[1];
sx q[1];
rz(2.0733817) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7601677) q[0];
sx q[0];
rz(-2.0695218) q[0];
sx q[0];
rz(-2.4295761) q[0];
rz(2.5184987) q[2];
sx q[2];
rz(-0.87298191) q[2];
sx q[2];
rz(-2.0562003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8074016) q[1];
sx q[1];
rz(-2.4535123) q[1];
sx q[1];
rz(3.1013558) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7925451) q[3];
sx q[3];
rz(-1.5132679) q[3];
sx q[3];
rz(-1.3105931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3855359) q[2];
sx q[2];
rz(-1.7588561) q[2];
sx q[2];
rz(-2.5059911) q[2];
rz(-0.33291891) q[3];
sx q[3];
rz(-1.0115441) q[3];
sx q[3];
rz(2.6788768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8024837) q[0];
sx q[0];
rz(-0.026263069) q[0];
sx q[0];
rz(3.0185757) q[0];
rz(-1.7440375) q[1];
sx q[1];
rz(-1.7536283) q[1];
sx q[1];
rz(2.3618598) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1542669) q[0];
sx q[0];
rz(-0.92181081) q[0];
sx q[0];
rz(2.799688) q[0];
rz(1.4238724) q[2];
sx q[2];
rz(-2.2511852) q[2];
sx q[2];
rz(0.18934393) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2998979) q[1];
sx q[1];
rz(-1.3720241) q[1];
sx q[1];
rz(-1.190462) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10752435) q[3];
sx q[3];
rz(-1.4443099) q[3];
sx q[3];
rz(1.3832987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.66703779) q[2];
sx q[2];
rz(-1.8755308) q[2];
sx q[2];
rz(3.0511268) q[2];
rz(-0.41680923) q[3];
sx q[3];
rz(-2.3495245) q[3];
sx q[3];
rz(-2.1478103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6592634) q[0];
sx q[0];
rz(-1.3220795) q[0];
sx q[0];
rz(-0.089476712) q[0];
rz(2.3327475) q[1];
sx q[1];
rz(-0.50061148) q[1];
sx q[1];
rz(-2.083875) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1093921) q[0];
sx q[0];
rz(-1.1940178) q[0];
sx q[0];
rz(-0.44372875) q[0];
rz(3.1359768) q[2];
sx q[2];
rz(-1.0188017) q[2];
sx q[2];
rz(2.698026) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8174929) q[1];
sx q[1];
rz(-1.1760532) q[1];
sx q[1];
rz(2.7604719) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1154409) q[3];
sx q[3];
rz(-2.3847178) q[3];
sx q[3];
rz(0.43318403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7384501) q[2];
sx q[2];
rz(-2.6046627) q[2];
sx q[2];
rz(0.37330791) q[2];
rz(0.25660723) q[3];
sx q[3];
rz(-1.720287) q[3];
sx q[3];
rz(0.074450113) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8858717) q[0];
sx q[0];
rz(-1.5789565) q[0];
sx q[0];
rz(-1.4174905) q[0];
rz(-1.4295084) q[1];
sx q[1];
rz(-2.7919339) q[1];
sx q[1];
rz(-2.3333593) q[1];
rz(-1.5612372) q[2];
sx q[2];
rz(-2.0295967) q[2];
sx q[2];
rz(1.5461736) q[2];
rz(0.96961602) q[3];
sx q[3];
rz(-1.6027228) q[3];
sx q[3];
rz(3.0016196) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

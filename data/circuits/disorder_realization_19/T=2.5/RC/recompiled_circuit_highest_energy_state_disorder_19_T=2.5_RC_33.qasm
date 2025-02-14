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
rz(2.1228696) q[0];
sx q[0];
rz(-2.2824204) q[0];
sx q[0];
rz(0.81508842) q[0];
rz(-1.1852784) q[1];
sx q[1];
rz(-1.4108682) q[1];
sx q[1];
rz(1.0676395) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0370705) q[0];
sx q[0];
rz(-2.7094954) q[0];
sx q[0];
rz(1.5140526) q[0];
x q[1];
rz(0.52402988) q[2];
sx q[2];
rz(-1.2175778) q[2];
sx q[2];
rz(-2.2729682) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31413919) q[1];
sx q[1];
rz(-1.2541391) q[1];
sx q[1];
rz(3.0400671) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2566363) q[3];
sx q[3];
rz(-1.3643294) q[3];
sx q[3];
rz(-1.5790758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.69433576) q[2];
sx q[2];
rz(-0.7434291) q[2];
sx q[2];
rz(1.8667963) q[2];
rz(2.6796807) q[3];
sx q[3];
rz(-0.67449823) q[3];
sx q[3];
rz(2.0089202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3502515) q[0];
sx q[0];
rz(-2.8332062) q[0];
sx q[0];
rz(-1.2530918) q[0];
rz(0.14532267) q[1];
sx q[1];
rz(-1.7456313) q[1];
sx q[1];
rz(-1.0911509) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017830124) q[0];
sx q[0];
rz(-1.995867) q[0];
sx q[0];
rz(0.2353038) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2107466) q[2];
sx q[2];
rz(-1.1384247) q[2];
sx q[2];
rz(-2.6170309) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1205463) q[1];
sx q[1];
rz(-1.7851549) q[1];
sx q[1];
rz(1.7477186) q[1];
rz(-pi) q[2];
rz(-1.1717779) q[3];
sx q[3];
rz(-1.1091091) q[3];
sx q[3];
rz(-0.53129133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.48065177) q[2];
sx q[2];
rz(-1.3903684) q[2];
sx q[2];
rz(-2.6118028) q[2];
rz(-0.79331136) q[3];
sx q[3];
rz(-1.5373693) q[3];
sx q[3];
rz(-2.5180499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6563501) q[0];
sx q[0];
rz(-2.1915477) q[0];
sx q[0];
rz(2.6128838) q[0];
rz(2.5953925) q[1];
sx q[1];
rz(-0.9587973) q[1];
sx q[1];
rz(0.34034696) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3211283) q[0];
sx q[0];
rz(-1.2437151) q[0];
sx q[0];
rz(-1.2820679) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7163926) q[2];
sx q[2];
rz(-1.3383838) q[2];
sx q[2];
rz(-0.53047859) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9593411) q[1];
sx q[1];
rz(-1.4762344) q[1];
sx q[1];
rz(1.574081) q[1];
rz(-pi) q[2];
rz(2.3119218) q[3];
sx q[3];
rz(-1.2874914) q[3];
sx q[3];
rz(2.0633351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1722374) q[2];
sx q[2];
rz(-2.7293971) q[2];
sx q[2];
rz(0.42984143) q[2];
rz(-0.75628453) q[3];
sx q[3];
rz(-0.1736621) q[3];
sx q[3];
rz(-0.82591301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38465685) q[0];
sx q[0];
rz(-1.5058368) q[0];
sx q[0];
rz(1.6480308) q[0];
rz(-2.3736296) q[1];
sx q[1];
rz(-0.47052828) q[1];
sx q[1];
rz(-0.54642645) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7554005) q[0];
sx q[0];
rz(-1.4996756) q[0];
sx q[0];
rz(3.0767308) q[0];
rz(-pi) q[1];
rz(0.77026639) q[2];
sx q[2];
rz(-2.4663743) q[2];
sx q[2];
rz(0.16083052) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5855693) q[1];
sx q[1];
rz(-1.6680191) q[1];
sx q[1];
rz(-0.14859622) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91937842) q[3];
sx q[3];
rz(-2.1438144) q[3];
sx q[3];
rz(-2.1185377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7118608) q[2];
sx q[2];
rz(-1.6542566) q[2];
sx q[2];
rz(-0.30409733) q[2];
rz(-1.7763304) q[3];
sx q[3];
rz(-1.920776) q[3];
sx q[3];
rz(2.064866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51233184) q[0];
sx q[0];
rz(-0.4442232) q[0];
sx q[0];
rz(-3.1410134) q[0];
rz(-2.0594788) q[1];
sx q[1];
rz(-2.6172456) q[1];
sx q[1];
rz(2.2023315) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1056553) q[0];
sx q[0];
rz(-2.2024184) q[0];
sx q[0];
rz(1.0211591) q[0];
rz(-2.1823857) q[2];
sx q[2];
rz(-1.2805174) q[2];
sx q[2];
rz(-0.91735754) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1402005) q[1];
sx q[1];
rz(-1.0085229) q[1];
sx q[1];
rz(-2.8881025) q[1];
rz(-pi) q[2];
rz(1.2573832) q[3];
sx q[3];
rz(-1.7872918) q[3];
sx q[3];
rz(-2.8683942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.038736343) q[2];
sx q[2];
rz(-1.2612217) q[2];
sx q[2];
rz(2.8803414) q[2];
rz(-2.2767565) q[3];
sx q[3];
rz(-1.4120925) q[3];
sx q[3];
rz(-0.29204667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84457266) q[0];
sx q[0];
rz(-1.427587) q[0];
sx q[0];
rz(-2.936506) q[0];
rz(1.0467485) q[1];
sx q[1];
rz(-1.3958684) q[1];
sx q[1];
rz(-2.7632025) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9721824) q[0];
sx q[0];
rz(-0.43146389) q[0];
sx q[0];
rz(1.1110825) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7360052) q[2];
sx q[2];
rz(-1.1045611) q[2];
sx q[2];
rz(-1.8777443) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.69120317) q[1];
sx q[1];
rz(-0.65237633) q[1];
sx q[1];
rz(-0.97263402) q[1];
rz(-pi) q[2];
rz(-1.1388402) q[3];
sx q[3];
rz(-2.216385) q[3];
sx q[3];
rz(1.0243675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64688524) q[2];
sx q[2];
rz(-1.9834221) q[2];
sx q[2];
rz(-1.0542487) q[2];
rz(-2.4833637) q[3];
sx q[3];
rz(-0.64526486) q[3];
sx q[3];
rz(-2.6575507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9996416) q[0];
sx q[0];
rz(-1.208409) q[0];
sx q[0];
rz(-2.5010338) q[0];
rz(-2.7236252) q[1];
sx q[1];
rz(-1.9211946) q[1];
sx q[1];
rz(2.3366065) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60123309) q[0];
sx q[0];
rz(-0.86358294) q[0];
sx q[0];
rz(1.5401031) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80193582) q[2];
sx q[2];
rz(-1.871043) q[2];
sx q[2];
rz(1.173347) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1272443) q[1];
sx q[1];
rz(-1.5509203) q[1];
sx q[1];
rz(-1.5424506) q[1];
rz(1.8667614) q[3];
sx q[3];
rz(-2.1419542) q[3];
sx q[3];
rz(-0.35863245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54887041) q[2];
sx q[2];
rz(-2.1465116) q[2];
sx q[2];
rz(-2.3804046) q[2];
rz(-2.393764) q[3];
sx q[3];
rz(-1.3176094) q[3];
sx q[3];
rz(2.4650011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6250896) q[0];
sx q[0];
rz(-2.857132) q[0];
sx q[0];
rz(1.4768584) q[0];
rz(0.45627108) q[1];
sx q[1];
rz(-1.4197333) q[1];
sx q[1];
rz(-2.2705073) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3251299) q[0];
sx q[0];
rz(-1.5583546) q[0];
sx q[0];
rz(2.5274656) q[0];
rz(-pi) q[1];
rz(-0.41042491) q[2];
sx q[2];
rz(-2.7685809) q[2];
sx q[2];
rz(-0.57508627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.58486119) q[1];
sx q[1];
rz(-1.24354) q[1];
sx q[1];
rz(1.3393988) q[1];
rz(-pi) q[2];
rz(-1.1444451) q[3];
sx q[3];
rz(-1.6255857) q[3];
sx q[3];
rz(0.74758119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2386834) q[2];
sx q[2];
rz(-0.52672714) q[2];
sx q[2];
rz(0.61863679) q[2];
rz(-0.020261852) q[3];
sx q[3];
rz(-2.198115) q[3];
sx q[3];
rz(2.3616135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063754931) q[0];
sx q[0];
rz(-1.5564593) q[0];
sx q[0];
rz(2.7177287) q[0];
rz(-0.18691143) q[1];
sx q[1];
rz(-0.82004768) q[1];
sx q[1];
rz(-1.4580457) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33015051) q[0];
sx q[0];
rz(-1.4528265) q[0];
sx q[0];
rz(-1.5850204) q[0];
rz(-pi) q[1];
rz(2.0844314) q[2];
sx q[2];
rz(-1.3328526) q[2];
sx q[2];
rz(0.4713716) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9002237) q[1];
sx q[1];
rz(-1.4321064) q[1];
sx q[1];
rz(1.6226416) q[1];
rz(1.6601172) q[3];
sx q[3];
rz(-1.2161939) q[3];
sx q[3];
rz(-1.8338628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3905048) q[2];
sx q[2];
rz(-1.0672528) q[2];
sx q[2];
rz(1.0941774) q[2];
rz(1.68082) q[3];
sx q[3];
rz(-1.7655617) q[3];
sx q[3];
rz(1.8028397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3222892) q[0];
sx q[0];
rz(-1.7604473) q[0];
sx q[0];
rz(1.639701) q[0];
rz(0.27885258) q[1];
sx q[1];
rz(-0.96062213) q[1];
sx q[1];
rz(-1.0241114) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49660027) q[0];
sx q[0];
rz(-1.4787714) q[0];
sx q[0];
rz(-0.49019469) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7488519) q[2];
sx q[2];
rz(-0.2789871) q[2];
sx q[2];
rz(1.6131608) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38899079) q[1];
sx q[1];
rz(-1.3704469) q[1];
sx q[1];
rz(0.2255535) q[1];
rz(0.34124438) q[3];
sx q[3];
rz(-2.2856718) q[3];
sx q[3];
rz(1.9174066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9378822) q[2];
sx q[2];
rz(-1.2869765) q[2];
sx q[2];
rz(2.6046275) q[2];
rz(-1.9133866) q[3];
sx q[3];
rz(-2.1824586) q[3];
sx q[3];
rz(2.8502407) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0252329) q[0];
sx q[0];
rz(-1.7740842) q[0];
sx q[0];
rz(-2.0728003) q[0];
rz(0.7069201) q[1];
sx q[1];
rz(-2.0126577) q[1];
sx q[1];
rz(-0.001002034) q[1];
rz(0.42264414) q[2];
sx q[2];
rz(-2.7012237) q[2];
sx q[2];
rz(1.2958432) q[2];
rz(-1.7892006) q[3];
sx q[3];
rz(-1.7159749) q[3];
sx q[3];
rz(0.80930474) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

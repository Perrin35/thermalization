OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34713137) q[0];
sx q[0];
rz(5.2678582) q[0];
sx q[0];
rz(9.8922748) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(1.7953035) q[1];
sx q[1];
rz(10.105159) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43359783) q[0];
sx q[0];
rz(-2.2888219) q[0];
sx q[0];
rz(-2.9684767) q[0];
rz(-pi) q[1];
x q[1];
rz(1.275185) q[2];
sx q[2];
rz(-1.7586853) q[2];
sx q[2];
rz(2.0762028) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2441794) q[1];
sx q[1];
rz(-1.3379828) q[1];
sx q[1];
rz(-2.0884872) q[1];
rz(2.666113) q[3];
sx q[3];
rz(-0.37215044) q[3];
sx q[3];
rz(-1.1572157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9900069) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(3.0541259) q[2];
rz(2.4123689) q[3];
sx q[3];
rz(-2.7681523) q[3];
sx q[3];
rz(0.27174404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7049578) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(-2.1402284) q[0];
rz(-0.17240605) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(0.52406812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8952626) q[0];
sx q[0];
rz(-1.2753914) q[0];
sx q[0];
rz(0.88460478) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0753724) q[2];
sx q[2];
rz(-2.3790857) q[2];
sx q[2];
rz(-2.5410595) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14428917) q[1];
sx q[1];
rz(-1.3376437) q[1];
sx q[1];
rz(-0.38645978) q[1];
x q[2];
rz(2.8855521) q[3];
sx q[3];
rz(-0.55169332) q[3];
sx q[3];
rz(-1.3224758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9841763) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(1.8015507) q[2];
rz(0.79483461) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(-3.1047344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79725093) q[0];
sx q[0];
rz(-0.69673711) q[0];
sx q[0];
rz(2.5168193) q[0];
rz(-2.1773188) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(2.952081) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7150869) q[0];
sx q[0];
rz(-2.3622892) q[0];
sx q[0];
rz(2.7515609) q[0];
rz(-pi) q[1];
rz(2.7663642) q[2];
sx q[2];
rz(-1.669075) q[2];
sx q[2];
rz(-2.3972942) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.0057919766) q[1];
sx q[1];
rz(-1.6264919) q[1];
sx q[1];
rz(0.23975753) q[1];
rz(0.96954815) q[3];
sx q[3];
rz(-1.4253972) q[3];
sx q[3];
rz(-1.5090514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0345962) q[2];
sx q[2];
rz(-0.84356374) q[2];
sx q[2];
rz(-1.3872046) q[2];
rz(-2.710279) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(-1.0881902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4819734) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(-1.5455998) q[0];
rz(-1.1384456) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(1.0901573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91256234) q[0];
sx q[0];
rz(-1.1704485) q[0];
sx q[0];
rz(1.1926665) q[0];
rz(2.1074739) q[2];
sx q[2];
rz(-0.82674971) q[2];
sx q[2];
rz(-2.2817051) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.66149368) q[1];
sx q[1];
rz(-1.1575932) q[1];
sx q[1];
rz(-0.6306298) q[1];
x q[2];
rz(2.7100536) q[3];
sx q[3];
rz(-2.6772237) q[3];
sx q[3];
rz(-0.73016703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.99889341) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(-2.2857655) q[2];
rz(1.9479729) q[3];
sx q[3];
rz(-1.5193628) q[3];
sx q[3];
rz(-2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6511433) q[0];
sx q[0];
rz(-2.1676846) q[0];
sx q[0];
rz(0.14973101) q[0];
rz(2.1504452) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(-1.9715462) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2802551) q[0];
sx q[0];
rz(-0.72685188) q[0];
sx q[0];
rz(-1.0759541) q[0];
rz(-pi) q[1];
rz(2.8180426) q[2];
sx q[2];
rz(-2.2252482) q[2];
sx q[2];
rz(-2.8138585) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.40844803) q[1];
sx q[1];
rz(-1.1322347) q[1];
sx q[1];
rz(-2.0088197) q[1];
x q[2];
rz(1.7647469) q[3];
sx q[3];
rz(-2.100482) q[3];
sx q[3];
rz(2.291631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.871792) q[2];
sx q[2];
rz(-1.5503927) q[2];
sx q[2];
rz(-3.0043547) q[2];
rz(-1.3373226) q[3];
sx q[3];
rz(-0.58115712) q[3];
sx q[3];
rz(1.8491245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542434) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(1.7911918) q[0];
rz(-0.84287914) q[1];
sx q[1];
rz(-2.402669) q[1];
sx q[1];
rz(-0.67289105) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56129365) q[0];
sx q[0];
rz(-1.6742047) q[0];
sx q[0];
rz(0.99594492) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.599309) q[2];
sx q[2];
rz(-0.71454222) q[2];
sx q[2];
rz(1.6376094) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3402965) q[1];
sx q[1];
rz(-0.9749037) q[1];
sx q[1];
rz(1.8227541) q[1];
rz(-2.9665885) q[3];
sx q[3];
rz(-2.724218) q[3];
sx q[3];
rz(1.1328896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.3488397) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(2.7065281) q[2];
rz(-1.3600291) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(2.8939261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9687987) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(1.0621747) q[0];
rz(2.0299714) q[1];
sx q[1];
rz(-1.9042791) q[1];
sx q[1];
rz(-1.7395082) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7715493) q[0];
sx q[0];
rz(-1.84329) q[0];
sx q[0];
rz(2.7605961) q[0];
x q[1];
rz(1.3349322) q[2];
sx q[2];
rz(-0.81996041) q[2];
sx q[2];
rz(1.6279398) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3558823) q[1];
sx q[1];
rz(-1.4175698) q[1];
sx q[1];
rz(-2.7760837) q[1];
rz(-0.13774638) q[3];
sx q[3];
rz(-2.1546954) q[3];
sx q[3];
rz(3.0018842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7193675) q[2];
sx q[2];
rz(-0.56754595) q[2];
sx q[2];
rz(-0.68022234) q[2];
rz(-0.42823544) q[3];
sx q[3];
rz(-1.2546344) q[3];
sx q[3];
rz(1.7817106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.632804) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(1.5493786) q[0];
rz(-0.24958615) q[1];
sx q[1];
rz(-1.1523749) q[1];
sx q[1];
rz(0.54135281) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68529785) q[0];
sx q[0];
rz(-2.1349499) q[0];
sx q[0];
rz(-1.5575404) q[0];
rz(2.6893822) q[2];
sx q[2];
rz(-1.8794249) q[2];
sx q[2];
rz(1.897915) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7975446) q[1];
sx q[1];
rz(-0.86595264) q[1];
sx q[1];
rz(1.0182347) q[1];
rz(3.1369282) q[3];
sx q[3];
rz(-0.38032535) q[3];
sx q[3];
rz(2.5120146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0970739) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(-2.3642335) q[2];
rz(0.85123953) q[3];
sx q[3];
rz(-2.080353) q[3];
sx q[3];
rz(1.3210993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38354307) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(0.0099649075) q[0];
rz(-1.0154356) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(-1.6962956) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1709258) q[0];
sx q[0];
rz(-2.2609841) q[0];
sx q[0];
rz(-0.70113457) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2465835) q[2];
sx q[2];
rz(-1.0969321) q[2];
sx q[2];
rz(-2.5788139) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0547486) q[1];
sx q[1];
rz(-2.5111685) q[1];
sx q[1];
rz(0.026442095) q[1];
rz(-pi) q[2];
rz(-0.61049283) q[3];
sx q[3];
rz(-2.1497512) q[3];
sx q[3];
rz(-1.0510774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4426667) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(0.60738579) q[2];
rz(-1.4390885) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8059175) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(0.6138531) q[0];
rz(1.0461668) q[1];
sx q[1];
rz(-2.8764953) q[1];
sx q[1];
rz(0.39224958) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8426725) q[0];
sx q[0];
rz(-0.93085104) q[0];
sx q[0];
rz(-2.0755033) q[0];
rz(-2.5160772) q[2];
sx q[2];
rz(-2.054347) q[2];
sx q[2];
rz(0.12550437) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0543921) q[1];
sx q[1];
rz(-1.1322081) q[1];
sx q[1];
rz(-0.88263504) q[1];
x q[2];
rz(-1.5997821) q[3];
sx q[3];
rz(-1.3658804) q[3];
sx q[3];
rz(-0.25587413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.067651) q[2];
sx q[2];
rz(-1.757688) q[2];
sx q[2];
rz(-0.60047853) q[2];
rz(2.0742119) q[3];
sx q[3];
rz(-1.821527) q[3];
sx q[3];
rz(1.0935121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7077211) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(2.1451163) q[1];
sx q[1];
rz(-1.6468208) q[1];
sx q[1];
rz(-1.5368808) q[1];
rz(-1.9267351) q[2];
sx q[2];
rz(-0.83649737) q[2];
sx q[2];
rz(2.5964824) q[2];
rz(2.2511803) q[3];
sx q[3];
rz(-2.012841) q[3];
sx q[3];
rz(-2.8750318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

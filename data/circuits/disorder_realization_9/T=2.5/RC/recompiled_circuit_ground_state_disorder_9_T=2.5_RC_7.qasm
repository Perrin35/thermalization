OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0661434) q[0];
sx q[0];
rz(-2.0976522) q[0];
sx q[0];
rz(3.1312842) q[0];
rz(2.161624) q[1];
sx q[1];
rz(-1.4449395) q[1];
sx q[1];
rz(2.565032) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4701472) q[0];
sx q[0];
rz(-1.2630442) q[0];
sx q[0];
rz(1.4438011) q[0];
rz(-0.61277436) q[2];
sx q[2];
rz(-0.98721993) q[2];
sx q[2];
rz(0.60223168) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.7652055) q[1];
sx q[1];
rz(-2.58316) q[1];
sx q[1];
rz(-2.8064578) q[1];
rz(-2.303894) q[3];
sx q[3];
rz(-1.695096) q[3];
sx q[3];
rz(1.4364786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.52446857) q[2];
sx q[2];
rz(-1.6813797) q[2];
sx q[2];
rz(1.2855533) q[2];
rz(-1.589795) q[3];
sx q[3];
rz(-1.0714622) q[3];
sx q[3];
rz(1.0514528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0153506) q[0];
sx q[0];
rz(-1.9168357) q[0];
sx q[0];
rz(2.8958564) q[0];
rz(-2.083678) q[1];
sx q[1];
rz(-1.825288) q[1];
sx q[1];
rz(0.33448514) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5366556) q[0];
sx q[0];
rz(-0.70022445) q[0];
sx q[0];
rz(2.3930153) q[0];
rz(-1.1673141) q[2];
sx q[2];
rz(-2.2902787) q[2];
sx q[2];
rz(-0.011842273) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1002754) q[1];
sx q[1];
rz(-1.9875437) q[1];
sx q[1];
rz(0.4112501) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7354419) q[3];
sx q[3];
rz(-0.63780071) q[3];
sx q[3];
rz(2.8075308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1841396) q[2];
sx q[2];
rz(-0.89749557) q[2];
sx q[2];
rz(-1.3857566) q[2];
rz(-0.98006788) q[3];
sx q[3];
rz(-2.6762784) q[3];
sx q[3];
rz(-3.0906299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4362713) q[0];
sx q[0];
rz(-0.956981) q[0];
sx q[0];
rz(1.3901688) q[0];
rz(0.35274371) q[1];
sx q[1];
rz(-2.2309062) q[1];
sx q[1];
rz(0.62044755) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1245956) q[0];
sx q[0];
rz(-0.11717883) q[0];
sx q[0];
rz(-2.4524053) q[0];
rz(-pi) q[1];
rz(-2.1587055) q[2];
sx q[2];
rz(-1.6574142) q[2];
sx q[2];
rz(-0.69332214) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.31252334) q[1];
sx q[1];
rz(-1.8644973) q[1];
sx q[1];
rz(-1.2306661) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1716072) q[3];
sx q[3];
rz(-0.78819345) q[3];
sx q[3];
rz(-1.5662686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1227526) q[2];
sx q[2];
rz(-0.41308013) q[2];
sx q[2];
rz(1.2472461) q[2];
rz(-2.9774169) q[3];
sx q[3];
rz(-1.7069867) q[3];
sx q[3];
rz(2.5119787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4267047) q[0];
sx q[0];
rz(-0.96788228) q[0];
sx q[0];
rz(1.2870652) q[0];
rz(-2.5413068) q[1];
sx q[1];
rz(-1.7900107) q[1];
sx q[1];
rz(-2.2379025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.569904) q[0];
sx q[0];
rz(-1.4635509) q[0];
sx q[0];
rz(1.5041758) q[0];
rz(-1.6389334) q[2];
sx q[2];
rz(-1.1916416) q[2];
sx q[2];
rz(3.1284297) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4111705) q[1];
sx q[1];
rz(-0.31283411) q[1];
sx q[1];
rz(0.61928796) q[1];
rz(-0.41208668) q[3];
sx q[3];
rz(-1.2743605) q[3];
sx q[3];
rz(-0.81338289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6220182) q[2];
sx q[2];
rz(-2.3081686) q[2];
sx q[2];
rz(-0.77345094) q[2];
rz(1.8526239) q[3];
sx q[3];
rz(-0.52923146) q[3];
sx q[3];
rz(-0.73498631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89161038) q[0];
sx q[0];
rz(-1.0050499) q[0];
sx q[0];
rz(-0.49215677) q[0];
rz(2.6026169) q[1];
sx q[1];
rz(-1.4424126) q[1];
sx q[1];
rz(-2.0416562) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5163706) q[0];
sx q[0];
rz(-1.3110135) q[0];
sx q[0];
rz(-0.35210877) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95910889) q[2];
sx q[2];
rz(-1.5796184) q[2];
sx q[2];
rz(2.9438643) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.089842794) q[1];
sx q[1];
rz(-2.3266092) q[1];
sx q[1];
rz(-3.0137193) q[1];
rz(1.8995011) q[3];
sx q[3];
rz(-1.1100169) q[3];
sx q[3];
rz(0.24390175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.34053549) q[2];
sx q[2];
rz(-2.8057782) q[2];
sx q[2];
rz(2.578793) q[2];
rz(-2.4417012) q[3];
sx q[3];
rz(-1.2908582) q[3];
sx q[3];
rz(0.25115299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3444779) q[0];
sx q[0];
rz(-0.7239224) q[0];
sx q[0];
rz(2.7864454) q[0];
rz(-1.2190602) q[1];
sx q[1];
rz(-1.2390169) q[1];
sx q[1];
rz(2.6893137) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4773524) q[0];
sx q[0];
rz(-1.1731804) q[0];
sx q[0];
rz(2.2009322) q[0];
rz(-pi) q[1];
rz(1.9777279) q[2];
sx q[2];
rz(-1.2135047) q[2];
sx q[2];
rz(-0.99321625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.86231316) q[1];
sx q[1];
rz(-1.24839) q[1];
sx q[1];
rz(-3.1161948) q[1];
rz(-2.5127453) q[3];
sx q[3];
rz(-1.4553242) q[3];
sx q[3];
rz(-0.94146282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5111115) q[2];
sx q[2];
rz(-0.52383542) q[2];
sx q[2];
rz(-0.13709489) q[2];
rz(-1.5591722) q[3];
sx q[3];
rz(-1.1538006) q[3];
sx q[3];
rz(2.3649575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1906076) q[0];
sx q[0];
rz(-2.3793716) q[0];
sx q[0];
rz(-1.5451587) q[0];
rz(2.6243788) q[1];
sx q[1];
rz(-1.1511753) q[1];
sx q[1];
rz(2.1545765) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4518829) q[0];
sx q[0];
rz(-1.5919884) q[0];
sx q[0];
rz(0.13084335) q[0];
rz(-pi) q[1];
rz(0.9632684) q[2];
sx q[2];
rz(-1.844256) q[2];
sx q[2];
rz(0.69404049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.034711866) q[1];
sx q[1];
rz(-1.3415404) q[1];
sx q[1];
rz(-2.8696612) q[1];
x q[2];
rz(0.084264755) q[3];
sx q[3];
rz(-1.8883421) q[3];
sx q[3];
rz(1.0061044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78642693) q[2];
sx q[2];
rz(-2.1106909) q[2];
sx q[2];
rz(-3.1332341) q[2];
rz(-0.23056325) q[3];
sx q[3];
rz(-1.2059261) q[3];
sx q[3];
rz(-1.4768538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.1281328) q[0];
sx q[0];
rz(-3.1210493) q[0];
sx q[0];
rz(-1.531456) q[0];
rz(1.6186591) q[1];
sx q[1];
rz(-1.5956655) q[1];
sx q[1];
rz(-1.5787554) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42199907) q[0];
sx q[0];
rz(-2.5223456) q[0];
sx q[0];
rz(2.8887755) q[0];
rz(1.7152856) q[2];
sx q[2];
rz(-1.943271) q[2];
sx q[2];
rz(2.5691751) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2819689) q[1];
sx q[1];
rz(-1.109237) q[1];
sx q[1];
rz(-2.0412316) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2501026) q[3];
sx q[3];
rz(-2.2134292) q[3];
sx q[3];
rz(1.2144517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83850399) q[2];
sx q[2];
rz(-1.9946626) q[2];
sx q[2];
rz(-0.037503555) q[2];
rz(2.904902) q[3];
sx q[3];
rz(-0.37875566) q[3];
sx q[3];
rz(-1.2934575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26838747) q[0];
sx q[0];
rz(-2.3516042) q[0];
sx q[0];
rz(-2.068212) q[0];
rz(2.7997596) q[1];
sx q[1];
rz(-2.5624202) q[1];
sx q[1];
rz(-0.62172186) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9699696) q[0];
sx q[0];
rz(-1.8299654) q[0];
sx q[0];
rz(0.39879946) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.694181) q[2];
sx q[2];
rz(-2.2272791) q[2];
sx q[2];
rz(2.3036164) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4704807) q[1];
sx q[1];
rz(-2.5501304) q[1];
sx q[1];
rz(-0.5670814) q[1];
x q[2];
rz(1.4831495) q[3];
sx q[3];
rz(-1.7109799) q[3];
sx q[3];
rz(0.46938458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5198034) q[2];
sx q[2];
rz(-0.18814627) q[2];
sx q[2];
rz(-0.56358799) q[2];
rz(-1.8300736) q[3];
sx q[3];
rz(-1.2389641) q[3];
sx q[3];
rz(-0.81609503) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9594864) q[0];
sx q[0];
rz(-2.2725548) q[0];
sx q[0];
rz(2.2667789) q[0];
rz(-1.733571) q[1];
sx q[1];
rz(-0.66649109) q[1];
sx q[1];
rz(-0.63546884) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68671526) q[0];
sx q[0];
rz(-1.8679163) q[0];
sx q[0];
rz(2.3700506) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1383912) q[2];
sx q[2];
rz(-0.48241189) q[2];
sx q[2];
rz(-2.8495827) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8567896) q[1];
sx q[1];
rz(-1.9724047) q[1];
sx q[1];
rz(1.0076341) q[1];
x q[2];
rz(0.54896574) q[3];
sx q[3];
rz(-2.0655051) q[3];
sx q[3];
rz(3.0567808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94329876) q[2];
sx q[2];
rz(-1.0903) q[2];
sx q[2];
rz(0.80149209) q[2];
rz(-1.9258026) q[3];
sx q[3];
rz(-0.91167584) q[3];
sx q[3];
rz(-3.099814) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53032482) q[0];
sx q[0];
rz(-1.7014736) q[0];
sx q[0];
rz(0.3076719) q[0];
rz(1.2904185) q[1];
sx q[1];
rz(-0.94480521) q[1];
sx q[1];
rz(0.31029846) q[1];
rz(0.63772884) q[2];
sx q[2];
rz(-2.9165033) q[2];
sx q[2];
rz(1.8115911) q[2];
rz(-0.65677662) q[3];
sx q[3];
rz(-1.6776424) q[3];
sx q[3];
rz(1.3959891) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

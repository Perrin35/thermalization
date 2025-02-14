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
rz(-1.0626592) q[0];
sx q[0];
rz(-0.81544977) q[0];
sx q[0];
rz(0.71001473) q[0];
rz(2.8275936) q[1];
sx q[1];
rz(-2.2065838) q[1];
sx q[1];
rz(-1.8097872) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3001271) q[0];
sx q[0];
rz(-1.3464768) q[0];
sx q[0];
rz(-2.8594144) q[0];
x q[1];
rz(-0.62628905) q[2];
sx q[2];
rz(-0.41538996) q[2];
sx q[2];
rz(0.20520575) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.59939042) q[1];
sx q[1];
rz(-2.1676461) q[1];
sx q[1];
rz(0.29920293) q[1];
rz(-pi) q[2];
rz(-2.2043742) q[3];
sx q[3];
rz(-1.9301747) q[3];
sx q[3];
rz(-1.0868974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72267246) q[2];
sx q[2];
rz(-1.924943) q[2];
sx q[2];
rz(-0.81532064) q[2];
rz(2.3319862) q[3];
sx q[3];
rz(-1.6428734) q[3];
sx q[3];
rz(2.0577551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0631436) q[0];
sx q[0];
rz(-2.0922631) q[0];
sx q[0];
rz(-3.0310042) q[0];
rz(2.1965006) q[1];
sx q[1];
rz(-2.7022305) q[1];
sx q[1];
rz(-1.5792712) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7720582) q[0];
sx q[0];
rz(-0.86449776) q[0];
sx q[0];
rz(0.95916551) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7209372) q[2];
sx q[2];
rz(-1.9845982) q[2];
sx q[2];
rz(-2.5977573) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4936714) q[1];
sx q[1];
rz(-2.646138) q[1];
sx q[1];
rz(-2.9205771) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82206313) q[3];
sx q[3];
rz(-1.2044319) q[3];
sx q[3];
rz(-1.0613393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9940146) q[2];
sx q[2];
rz(-2.9958604) q[2];
sx q[2];
rz(-2.6914524) q[2];
rz(-1.133793) q[3];
sx q[3];
rz(-1.4530051) q[3];
sx q[3];
rz(-2.1639737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56871539) q[0];
sx q[0];
rz(-2.1935538) q[0];
sx q[0];
rz(-0.28133389) q[0];
rz(1.7549134) q[1];
sx q[1];
rz(-1.7117056) q[1];
sx q[1];
rz(-2.5414355) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2413797) q[0];
sx q[0];
rz(-0.29415392) q[0];
sx q[0];
rz(2.9808729) q[0];
rz(-pi) q[1];
x q[1];
rz(0.05608989) q[2];
sx q[2];
rz(-2.1768835) q[2];
sx q[2];
rz(-2.1457399) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.197613) q[1];
sx q[1];
rz(-0.67786067) q[1];
sx q[1];
rz(2.6070836) q[1];
rz(2.568432) q[3];
sx q[3];
rz(-1.7148682) q[3];
sx q[3];
rz(-2.746563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0448138) q[2];
sx q[2];
rz(-1.1324977) q[2];
sx q[2];
rz(-2.0951994) q[2];
rz(0.3977631) q[3];
sx q[3];
rz(-0.64266509) q[3];
sx q[3];
rz(-0.1325632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9699049) q[0];
sx q[0];
rz(-3.1145018) q[0];
sx q[0];
rz(1.0525674) q[0];
rz(-1.9453847) q[1];
sx q[1];
rz(-1.2095249) q[1];
sx q[1];
rz(0.41950163) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9297816) q[0];
sx q[0];
rz(-2.503184) q[0];
sx q[0];
rz(0.45244658) q[0];
rz(-pi) q[1];
rz(2.9151462) q[2];
sx q[2];
rz(-2.1927877) q[2];
sx q[2];
rz(0.43606191) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7345404) q[1];
sx q[1];
rz(-2.190935) q[1];
sx q[1];
rz(2.6728515) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.64873) q[3];
sx q[3];
rz(-0.56418428) q[3];
sx q[3];
rz(-1.197937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9516051) q[2];
sx q[2];
rz(-1.1185458) q[2];
sx q[2];
rz(2.0513963) q[2];
rz(-0.53127855) q[3];
sx q[3];
rz(-0.71610206) q[3];
sx q[3];
rz(1.5268911) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7200318) q[0];
sx q[0];
rz(-1.5694542) q[0];
sx q[0];
rz(-1.3979727) q[0];
rz(-3.0086503) q[1];
sx q[1];
rz(-1.839919) q[1];
sx q[1];
rz(3.1403819) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.557363) q[0];
sx q[0];
rz(-1.5869291) q[0];
sx q[0];
rz(-2.7610012) q[0];
rz(-pi) q[1];
rz(-2.9820061) q[2];
sx q[2];
rz(-1.288646) q[2];
sx q[2];
rz(2.7196376) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19211719) q[1];
sx q[1];
rz(-1.8743974) q[1];
sx q[1];
rz(2.4785701) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3962383) q[3];
sx q[3];
rz(-2.9896185) q[3];
sx q[3];
rz(-2.0106955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91263897) q[2];
sx q[2];
rz(-1.8043844) q[2];
sx q[2];
rz(0.21052989) q[2];
rz(-1.0362961) q[3];
sx q[3];
rz(-2.3573037) q[3];
sx q[3];
rz(1.9265296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-1.9482816) q[0];
sx q[0];
rz(-1.223215) q[0];
sx q[0];
rz(-0.40570983) q[0];
rz(1.3544719) q[1];
sx q[1];
rz(-0.82258075) q[1];
sx q[1];
rz(2.2192661) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8278481) q[0];
sx q[0];
rz(-1.8593615) q[0];
sx q[0];
rz(0.11103481) q[0];
rz(-pi) q[1];
rz(-0.90409235) q[2];
sx q[2];
rz(-1.6518001) q[2];
sx q[2];
rz(-1.9503649) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4757918) q[1];
sx q[1];
rz(-2.5195055) q[1];
sx q[1];
rz(1.3740963) q[1];
rz(1.9590553) q[3];
sx q[3];
rz(-2.0960596) q[3];
sx q[3];
rz(-1.0356667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7102082) q[2];
sx q[2];
rz(-1.2146543) q[2];
sx q[2];
rz(0.40194884) q[2];
rz(1.2881783) q[3];
sx q[3];
rz(-0.86789075) q[3];
sx q[3];
rz(-1.8624381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61889082) q[0];
sx q[0];
rz(-2.0382477) q[0];
sx q[0];
rz(3.0772305) q[0];
rz(-2.3035658) q[1];
sx q[1];
rz(-2.3152654) q[1];
sx q[1];
rz(-1.3305957) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9143599) q[0];
sx q[0];
rz(-1.6961251) q[0];
sx q[0];
rz(-2.5663239) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66489545) q[2];
sx q[2];
rz(-0.722363) q[2];
sx q[2];
rz(1.6160053) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45695282) q[1];
sx q[1];
rz(-1.350953) q[1];
sx q[1];
rz(2.6067642) q[1];
rz(-1.8136386) q[3];
sx q[3];
rz(-0.62627072) q[3];
sx q[3];
rz(-2.0988718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4055206) q[2];
sx q[2];
rz(-1.5084927) q[2];
sx q[2];
rz(-0.68354496) q[2];
rz(-0.73379597) q[3];
sx q[3];
rz(-1.4132696) q[3];
sx q[3];
rz(2.2939513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.7965294) q[0];
sx q[0];
rz(-1.0836443) q[0];
sx q[0];
rz(-1.1731359) q[0];
rz(2.7660811) q[1];
sx q[1];
rz(-2.651732) q[1];
sx q[1];
rz(2.1048022) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0309854) q[0];
sx q[0];
rz(-1.5944423) q[0];
sx q[0];
rz(1.551669) q[0];
x q[1];
rz(-2.1027742) q[2];
sx q[2];
rz(-2.8237282) q[2];
sx q[2];
rz(1.3587855) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1634961) q[1];
sx q[1];
rz(-1.4304203) q[1];
sx q[1];
rz(0.76444334) q[1];
rz(2.6508207) q[3];
sx q[3];
rz(-1.0939071) q[3];
sx q[3];
rz(-2.8618438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5658687) q[2];
sx q[2];
rz(-0.91369358) q[2];
sx q[2];
rz(-2.4110528) q[2];
rz(0.20914397) q[3];
sx q[3];
rz(-2.6181965) q[3];
sx q[3];
rz(1.8714347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65649477) q[0];
sx q[0];
rz(-0.42335835) q[0];
sx q[0];
rz(-0.04743162) q[0];
rz(2.9196396) q[1];
sx q[1];
rz(-1.0739948) q[1];
sx q[1];
rz(-2.2304631) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.52235) q[0];
sx q[0];
rz(-1.5372835) q[0];
sx q[0];
rz(1.5438118) q[0];
rz(-pi) q[1];
x q[1];
rz(1.59052) q[2];
sx q[2];
rz(-2.4998186) q[2];
sx q[2];
rz(-2.3373147) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4998656) q[1];
sx q[1];
rz(-1.2033101) q[1];
sx q[1];
rz(1.5843452) q[1];
x q[2];
rz(-2.9373475) q[3];
sx q[3];
rz(-1.0139272) q[3];
sx q[3];
rz(-2.000252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2719416) q[2];
sx q[2];
rz(-0.38001529) q[2];
sx q[2];
rz(2.9730049) q[2];
rz(-0.21909675) q[3];
sx q[3];
rz(-0.85182652) q[3];
sx q[3];
rz(-1.8144089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.99437) q[0];
sx q[0];
rz(-2.8012025) q[0];
sx q[0];
rz(-0.45743531) q[0];
rz(1.2262729) q[1];
sx q[1];
rz(-2.8044082) q[1];
sx q[1];
rz(1.7785243) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5209893) q[0];
sx q[0];
rz(-0.61043569) q[0];
sx q[0];
rz(3.0362398) q[0];
rz(1.326637) q[2];
sx q[2];
rz(-2.1863822) q[2];
sx q[2];
rz(0.89060874) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.41238989) q[1];
sx q[1];
rz(-2.7250184) q[1];
sx q[1];
rz(-2.7135486) q[1];
rz(-1.4665589) q[3];
sx q[3];
rz(-0.6851894) q[3];
sx q[3];
rz(0.70885056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7676131) q[2];
sx q[2];
rz(-2.5869936) q[2];
sx q[2];
rz(-0.45200959) q[2];
rz(-0.40397817) q[3];
sx q[3];
rz(-1.2510866) q[3];
sx q[3];
rz(-2.0154791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44611888) q[0];
sx q[0];
rz(-1.6071381) q[0];
sx q[0];
rz(0.18679609) q[0];
rz(2.6192464) q[1];
sx q[1];
rz(-0.53032395) q[1];
sx q[1];
rz(1.2293336) q[1];
rz(1.8277373) q[2];
sx q[2];
rz(-1.8443454) q[2];
sx q[2];
rz(0.92583427) q[2];
rz(-1.7200851) q[3];
sx q[3];
rz(-1.310077) q[3];
sx q[3];
rz(2.6826912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1801017) q[0];
sx q[0];
rz(-0.56668007) q[0];
sx q[0];
rz(2.3508747) q[0];
rz(2.2136731) q[1];
sx q[1];
rz(-3.1065431) q[1];
sx q[1];
rz(-0.30326453) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68214455) q[0];
sx q[0];
rz(-1.7984227) q[0];
sx q[0];
rz(-2.4683003) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35776414) q[2];
sx q[2];
rz(-1.9310037) q[2];
sx q[2];
rz(1.2586762) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11285148) q[1];
sx q[1];
rz(-1.0489221) q[1];
sx q[1];
rz(2.584409) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4879151) q[3];
sx q[3];
rz(-0.83717665) q[3];
sx q[3];
rz(2.8745911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.37320331) q[2];
sx q[2];
rz(-0.36036569) q[2];
sx q[2];
rz(-2.4599794) q[2];
rz(2.5810589) q[3];
sx q[3];
rz(-0.78442854) q[3];
sx q[3];
rz(-2.4973629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19175567) q[0];
sx q[0];
rz(-2.2332709) q[0];
sx q[0];
rz(2.7857696) q[0];
rz(0.25335723) q[1];
sx q[1];
rz(-2.2919877) q[1];
sx q[1];
rz(1.8960948) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6424774) q[0];
sx q[0];
rz(-1.0986975) q[0];
sx q[0];
rz(-1.319122) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0959082) q[2];
sx q[2];
rz(-2.137299) q[2];
sx q[2];
rz(-2.1691466) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8207606) q[1];
sx q[1];
rz(-2.09873) q[1];
sx q[1];
rz(1.072701) q[1];
x q[2];
rz(-0.50470501) q[3];
sx q[3];
rz(-2.0619806) q[3];
sx q[3];
rz(-2.1757656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7443098) q[2];
sx q[2];
rz(-0.77703589) q[2];
sx q[2];
rz(-3.0139319) q[2];
rz(2.4658261) q[3];
sx q[3];
rz(-2.6830169) q[3];
sx q[3];
rz(-2.7202386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0965939) q[0];
sx q[0];
rz(-2.8969722) q[0];
sx q[0];
rz(2.007572) q[0];
rz(0.89316142) q[1];
sx q[1];
rz(-2.2375219) q[1];
sx q[1];
rz(1.2718511) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24090996) q[0];
sx q[0];
rz(-2.5782452) q[0];
sx q[0];
rz(1.6202089) q[0];
rz(-1.8142305) q[2];
sx q[2];
rz(-2.8231502) q[2];
sx q[2];
rz(1.3322347) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33128502) q[1];
sx q[1];
rz(-2.2687817) q[1];
sx q[1];
rz(-0.26878243) q[1];
x q[2];
rz(1.351139) q[3];
sx q[3];
rz(-1.7336415) q[3];
sx q[3];
rz(1.6514678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21358718) q[2];
sx q[2];
rz(-2.2760133) q[2];
sx q[2];
rz(-1.340284) q[2];
rz(-1.853893) q[3];
sx q[3];
rz(-1.6940593) q[3];
sx q[3];
rz(-0.48937669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5260148) q[0];
sx q[0];
rz(-0.57681334) q[0];
sx q[0];
rz(2.4082129) q[0];
rz(-2.3461657) q[1];
sx q[1];
rz(-2.9454102) q[1];
sx q[1];
rz(-1.1682074) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35287898) q[0];
sx q[0];
rz(-1.9813992) q[0];
sx q[0];
rz(-0.11524185) q[0];
rz(-pi) q[1];
rz(-2.2890786) q[2];
sx q[2];
rz(-2.3688683) q[2];
sx q[2];
rz(-1.0393522) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.296098) q[1];
sx q[1];
rz(-1.6137527) q[1];
sx q[1];
rz(1.6011493) q[1];
rz(0.73428728) q[3];
sx q[3];
rz(-0.81273116) q[3];
sx q[3];
rz(1.6664291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1802133) q[2];
sx q[2];
rz(-2.3258379) q[2];
sx q[2];
rz(0.95480314) q[2];
rz(-2.09156) q[3];
sx q[3];
rz(-2.3539216) q[3];
sx q[3];
rz(-2.5900456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0410205) q[0];
sx q[0];
rz(-2.6028778) q[0];
sx q[0];
rz(-2.4464497) q[0];
rz(1.2160542) q[1];
sx q[1];
rz(-1.0168409) q[1];
sx q[1];
rz(3.1207808) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6489844) q[0];
sx q[0];
rz(-2.2243315) q[0];
sx q[0];
rz(0.18006353) q[0];
rz(-pi) q[1];
rz(1.928442) q[2];
sx q[2];
rz(-0.76947533) q[2];
sx q[2];
rz(0.77126743) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0328524) q[1];
sx q[1];
rz(-1.5592062) q[1];
sx q[1];
rz(-1.3316915) q[1];
rz(-pi) q[2];
rz(0.5848123) q[3];
sx q[3];
rz(-0.11128128) q[3];
sx q[3];
rz(2.5113784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9713328) q[2];
sx q[2];
rz(-0.33997619) q[2];
sx q[2];
rz(-2.5374832) q[2];
rz(-0.64770925) q[3];
sx q[3];
rz(-0.52102399) q[3];
sx q[3];
rz(0.46085301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0172417) q[0];
sx q[0];
rz(-2.0483973) q[0];
sx q[0];
rz(-0.35853115) q[0];
rz(-2.5784946) q[1];
sx q[1];
rz(-2.2393354) q[1];
sx q[1];
rz(0.30555746) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2222851) q[0];
sx q[0];
rz(-0.74757517) q[0];
sx q[0];
rz(-2.5966552) q[0];
x q[1];
rz(-0.75336918) q[2];
sx q[2];
rz(-2.8940331) q[2];
sx q[2];
rz(-2.2581429) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.36781947) q[1];
sx q[1];
rz(-2.1065284) q[1];
sx q[1];
rz(2.1526835) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0642387) q[3];
sx q[3];
rz(-1.6408444) q[3];
sx q[3];
rz(-1.182568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87865889) q[2];
sx q[2];
rz(-2.1808251) q[2];
sx q[2];
rz(3.089454) q[2];
rz(2.7122998) q[3];
sx q[3];
rz(-1.7310127) q[3];
sx q[3];
rz(-2.8894292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.6843863) q[0];
sx q[0];
rz(-2.6173213) q[0];
sx q[0];
rz(-2.9065409) q[0];
rz(-2.5382407) q[1];
sx q[1];
rz(-0.82913202) q[1];
sx q[1];
rz(2.5650909) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.910945) q[0];
sx q[0];
rz(-1.0663138) q[0];
sx q[0];
rz(-0.77433681) q[0];
rz(-pi) q[1];
rz(-2.3660955) q[2];
sx q[2];
rz(-1.3243054) q[2];
sx q[2];
rz(1.6131439) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.41771463) q[1];
sx q[1];
rz(-0.52199328) q[1];
sx q[1];
rz(2.2531409) q[1];
rz(-pi) q[2];
rz(-1.0888388) q[3];
sx q[3];
rz(-1.4432943) q[3];
sx q[3];
rz(1.8021999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0779886) q[2];
sx q[2];
rz(-2.3006907) q[2];
sx q[2];
rz(-1.7951175) q[2];
rz(0.50955647) q[3];
sx q[3];
rz(-1.2600803) q[3];
sx q[3];
rz(2.8349561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6154196) q[0];
sx q[0];
rz(-1.0831447) q[0];
sx q[0];
rz(-0.35588595) q[0];
rz(0.7016167) q[1];
sx q[1];
rz(-2.9569148) q[1];
sx q[1];
rz(0.96431771) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.334737) q[0];
sx q[0];
rz(-0.47602326) q[0];
sx q[0];
rz(-1.8658338) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.020267296) q[2];
sx q[2];
rz(-1.8635443) q[2];
sx q[2];
rz(1.0525289) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6881975) q[1];
sx q[1];
rz(-1.5904237) q[1];
sx q[1];
rz(-0.49491306) q[1];
x q[2];
rz(-2.8300222) q[3];
sx q[3];
rz(-1.4928679) q[3];
sx q[3];
rz(3.1337217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4926766) q[2];
sx q[2];
rz(-0.70768386) q[2];
sx q[2];
rz(2.3259582) q[2];
rz(-2.2260769) q[3];
sx q[3];
rz(-0.86796498) q[3];
sx q[3];
rz(-0.22667949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30064073) q[0];
sx q[0];
rz(-0.17423593) q[0];
sx q[0];
rz(1.9006282) q[0];
rz(-0.10494431) q[1];
sx q[1];
rz(-0.34067708) q[1];
sx q[1];
rz(2.6822283) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37840415) q[0];
sx q[0];
rz(-1.8109461) q[0];
sx q[0];
rz(1.9394919) q[0];
x q[1];
rz(2.8723889) q[2];
sx q[2];
rz(-1.2503617) q[2];
sx q[2];
rz(-1.073369) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8018559) q[1];
sx q[1];
rz(-1.1366399) q[1];
sx q[1];
rz(-1.3375907) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9271803) q[3];
sx q[3];
rz(-2.5567856) q[3];
sx q[3];
rz(-2.6406276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0505872) q[2];
sx q[2];
rz(-0.688474) q[2];
sx q[2];
rz(3.0509994) q[2];
rz(0.76720864) q[3];
sx q[3];
rz(-1.6588666) q[3];
sx q[3];
rz(-1.4513133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7200658) q[0];
sx q[0];
rz(-0.92362112) q[0];
sx q[0];
rz(-1.6044755) q[0];
rz(-0.9556669) q[1];
sx q[1];
rz(-1.6803398) q[1];
sx q[1];
rz(2.59424) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5794649) q[0];
sx q[0];
rz(-1.5766931) q[0];
sx q[0];
rz(1.6505169) q[0];
rz(-pi) q[1];
rz(1.1905306) q[2];
sx q[2];
rz(-2.0693738) q[2];
sx q[2];
rz(-0.42344365) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9319084) q[1];
sx q[1];
rz(-2.3386777) q[1];
sx q[1];
rz(-1.241339) q[1];
rz(2.5741816) q[3];
sx q[3];
rz(-2.3443522) q[3];
sx q[3];
rz(-0.85464961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0305816) q[2];
sx q[2];
rz(-2.9647398) q[2];
sx q[2];
rz(2.492823) q[2];
rz(-2.9099921) q[3];
sx q[3];
rz(-2.2571371) q[3];
sx q[3];
rz(-3.054936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34407525) q[0];
sx q[0];
rz(-1.2114914) q[0];
sx q[0];
rz(-1.399566) q[0];
rz(0.65921417) q[1];
sx q[1];
rz(-1.2792239) q[1];
sx q[1];
rz(-1.4469133) q[1];
rz(-1.7991015) q[2];
sx q[2];
rz(-1.5226428) q[2];
sx q[2];
rz(-1.5626283) q[2];
rz(-3.0869879) q[3];
sx q[3];
rz(-1.7163897) q[3];
sx q[3];
rz(-0.39672273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

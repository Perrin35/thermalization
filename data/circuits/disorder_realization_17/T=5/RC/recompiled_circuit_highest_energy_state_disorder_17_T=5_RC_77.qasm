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
rz(-2.4315779) q[0];
rz(2.8275936) q[1];
sx q[1];
rz(-2.2065838) q[1];
sx q[1];
rz(1.3318055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0667257) q[0];
sx q[0];
rz(-0.35860379) q[0];
sx q[0];
rz(0.68645607) q[0];
rz(-pi) q[1];
rz(-0.34320404) q[2];
sx q[2];
rz(-1.3319974) q[2];
sx q[2];
rz(-0.78087872) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0441452) q[1];
sx q[1];
rz(-2.4821979) q[1];
sx q[1];
rz(1.9800746) q[1];
rz(-2.136257) q[3];
sx q[3];
rz(-0.71601235) q[3];
sx q[3];
rz(2.211192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4189202) q[2];
sx q[2];
rz(-1.2166497) q[2];
sx q[2];
rz(2.326272) q[2];
rz(0.80960649) q[3];
sx q[3];
rz(-1.6428734) q[3];
sx q[3];
rz(1.0838375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0631436) q[0];
sx q[0];
rz(-2.0922631) q[0];
sx q[0];
rz(0.11058841) q[0];
rz(-0.94509205) q[1];
sx q[1];
rz(-2.7022305) q[1];
sx q[1];
rz(-1.5792712) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54540578) q[0];
sx q[0];
rz(-2.2431787) q[0];
sx q[0];
rz(2.5491712) q[0];
rz(-pi) q[1];
rz(-1.1223824) q[2];
sx q[2];
rz(-1.9539991) q[2];
sx q[2];
rz(0.84898432) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4936714) q[1];
sx q[1];
rz(-0.4954547) q[1];
sx q[1];
rz(-0.22101553) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.14757806) q[2];
sx q[2];
rz(-2.9958604) q[2];
sx q[2];
rz(-0.4501403) q[2];
rz(-1.133793) q[3];
sx q[3];
rz(-1.4530051) q[3];
sx q[3];
rz(-2.1639737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.5728773) q[0];
sx q[0];
rz(-0.94803888) q[0];
sx q[0];
rz(0.28133389) q[0];
rz(1.7549134) q[1];
sx q[1];
rz(-1.7117056) q[1];
sx q[1];
rz(0.6001572) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0735787) q[0];
sx q[0];
rz(-1.2805443) q[0];
sx q[0];
rz(1.5223548) q[0];
rz(-0.96397206) q[2];
sx q[2];
rz(-1.5247048) q[2];
sx q[2];
rz(2.5986236) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0827209) q[1];
sx q[1];
rz(-1.8959672) q[1];
sx q[1];
rz(-0.60589686) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3998109) q[3];
sx q[3];
rz(-2.1372843) q[3];
sx q[3];
rz(-1.2681707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0967789) q[2];
sx q[2];
rz(-2.009095) q[2];
sx q[2];
rz(1.0463932) q[2];
rz(0.3977631) q[3];
sx q[3];
rz(-0.64266509) q[3];
sx q[3];
rz(-0.1325632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9699049) q[0];
sx q[0];
rz(-3.1145018) q[0];
sx q[0];
rz(-1.0525674) q[0];
rz(1.196208) q[1];
sx q[1];
rz(-1.2095249) q[1];
sx q[1];
rz(0.41950163) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1284258) q[0];
sx q[0];
rz(-1.3072398) q[0];
sx q[0];
rz(-2.5530784) q[0];
rz(2.9151462) q[2];
sx q[2];
rz(-0.94880494) q[2];
sx q[2];
rz(-0.43606191) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7345404) q[1];
sx q[1];
rz(-0.95065763) q[1];
sx q[1];
rz(2.6728515) q[1];
rz(1.2798843) q[3];
sx q[3];
rz(-2.0613163) q[3];
sx q[3];
rz(-1.377493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18998751) q[2];
sx q[2];
rz(-2.0230468) q[2];
sx q[2];
rz(2.0513963) q[2];
rz(-2.6103141) q[3];
sx q[3];
rz(-2.4254906) q[3];
sx q[3];
rz(1.5268911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4215609) q[0];
sx q[0];
rz(-1.5694542) q[0];
sx q[0];
rz(-1.74362) q[0];
rz(-3.0086503) q[1];
sx q[1];
rz(-1.3016737) q[1];
sx q[1];
rz(0.001210777) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1346136) q[0];
sx q[0];
rz(-1.9513357) q[0];
sx q[0];
rz(-1.5881722) q[0];
x q[1];
rz(-2.0722449) q[2];
sx q[2];
rz(-0.32310383) q[2];
sx q[2];
rz(2.1954775) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1493133) q[1];
sx q[1];
rz(-0.94301254) q[1];
sx q[1];
rz(1.1924442) q[1];
rz(0.026592606) q[3];
sx q[3];
rz(-1.4211492) q[3];
sx q[3];
rz(-0.95434556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2289537) q[2];
sx q[2];
rz(-1.3372083) q[2];
sx q[2];
rz(2.9310628) q[2];
rz(-2.1052965) q[3];
sx q[3];
rz(-2.3573037) q[3];
sx q[3];
rz(1.2150631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1933111) q[0];
sx q[0];
rz(-1.223215) q[0];
sx q[0];
rz(2.7358828) q[0];
rz(1.7871208) q[1];
sx q[1];
rz(-2.3190119) q[1];
sx q[1];
rz(2.2192661) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2012537) q[0];
sx q[0];
rz(-2.8329599) q[0];
sx q[0];
rz(-1.9280532) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7013266) q[2];
sx q[2];
rz(-0.67085999) q[2];
sx q[2];
rz(-2.6595569) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.90627024) q[1];
sx q[1];
rz(-2.1791237) q[1];
sx q[1];
rz(-3.0023605) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5629696) q[3];
sx q[3];
rz(-2.4994183) q[3];
sx q[3];
rz(-0.35143055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7102082) q[2];
sx q[2];
rz(-1.9269383) q[2];
sx q[2];
rz(2.7396438) q[2];
rz(-1.8534144) q[3];
sx q[3];
rz(-2.2737019) q[3];
sx q[3];
rz(1.8624381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5227018) q[0];
sx q[0];
rz(-2.0382477) q[0];
sx q[0];
rz(-3.0772305) q[0];
rz(2.3035658) q[1];
sx q[1];
rz(-0.82632724) q[1];
sx q[1];
rz(1.8109969) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7171519) q[0];
sx q[0];
rz(-2.1409875) q[0];
sx q[0];
rz(1.7198404) q[0];
x q[1];
rz(2.4766972) q[2];
sx q[2];
rz(-0.722363) q[2];
sx q[2];
rz(1.6160053) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6846398) q[1];
sx q[1];
rz(-1.7906396) q[1];
sx q[1];
rz(-2.6067642) q[1];
x q[2];
rz(1.8136386) q[3];
sx q[3];
rz(-2.5153219) q[3];
sx q[3];
rz(-2.0988718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4055206) q[2];
sx q[2];
rz(-1.6330999) q[2];
sx q[2];
rz(2.4580477) q[2];
rz(-0.73379597) q[3];
sx q[3];
rz(-1.7283231) q[3];
sx q[3];
rz(0.84764135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.34506327) q[0];
sx q[0];
rz(-2.0579484) q[0];
sx q[0];
rz(-1.9684568) q[0];
rz(-2.7660811) q[1];
sx q[1];
rz(-0.48986062) q[1];
sx q[1];
rz(2.1048022) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6818559) q[0];
sx q[0];
rz(-1.5516743) q[0];
sx q[0];
rz(-3.1179423) q[0];
rz(-1.2944968) q[2];
sx q[2];
rz(-1.7299998) q[2];
sx q[2];
rz(-2.4197848) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1634961) q[1];
sx q[1];
rz(-1.4304203) q[1];
sx q[1];
rz(-2.3771493) q[1];
rz(-pi) q[2];
x q[2];
rz(2.100714) q[3];
sx q[3];
rz(-1.138759) q[3];
sx q[3];
rz(2.0910859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5658687) q[2];
sx q[2];
rz(-2.2278991) q[2];
sx q[2];
rz(2.4110528) q[2];
rz(0.20914397) q[3];
sx q[3];
rz(-2.6181965) q[3];
sx q[3];
rz(1.8714347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65649477) q[0];
sx q[0];
rz(-0.42335835) q[0];
sx q[0];
rz(0.04743162) q[0];
rz(0.221953) q[1];
sx q[1];
rz(-2.0675979) q[1];
sx q[1];
rz(0.91112959) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.297356) q[0];
sx q[0];
rz(-0.043023303) q[0];
sx q[0];
rz(2.4639315) q[0];
rz(-2.2124771) q[2];
sx q[2];
rz(-1.5826028) q[2];
sx q[2];
rz(-2.3908743) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.67942372) q[1];
sx q[1];
rz(-0.36772455) q[1];
sx q[1];
rz(-0.035178758) q[1];
rz(0.20424517) q[3];
sx q[3];
rz(-1.0139272) q[3];
sx q[3];
rz(-2.000252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.86965108) q[2];
sx q[2];
rz(-2.7615774) q[2];
sx q[2];
rz(0.16858777) q[2];
rz(-0.21909675) q[3];
sx q[3];
rz(-0.85182652) q[3];
sx q[3];
rz(-1.8144089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1472226) q[0];
sx q[0];
rz(-0.34039012) q[0];
sx q[0];
rz(0.45743531) q[0];
rz(1.9153197) q[1];
sx q[1];
rz(-2.8044082) q[1];
sx q[1];
rz(1.3630684) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03662388) q[0];
sx q[0];
rz(-1.6311121) q[0];
sx q[0];
rz(0.6078267) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62984939) q[2];
sx q[2];
rz(-1.3721264) q[2];
sx q[2];
rz(-0.82306403) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0911965) q[1];
sx q[1];
rz(-1.1938057) q[1];
sx q[1];
rz(-1.3891549) q[1];
rz(-pi) q[2];
rz(1.6750338) q[3];
sx q[3];
rz(-0.6851894) q[3];
sx q[3];
rz(-2.4327421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7676131) q[2];
sx q[2];
rz(-0.55459905) q[2];
sx q[2];
rz(-2.6895831) q[2];
rz(2.7376145) q[3];
sx q[3];
rz(-1.890506) q[3];
sx q[3];
rz(-1.1261136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6954738) q[0];
sx q[0];
rz(-1.5344545) q[0];
sx q[0];
rz(-2.9547966) q[0];
rz(0.52234621) q[1];
sx q[1];
rz(-2.6112687) q[1];
sx q[1];
rz(-1.912259) q[1];
rz(2.4056388) q[2];
sx q[2];
rz(-2.7685168) q[2];
sx q[2];
rz(-1.4442486) q[2];
rz(1.4215076) q[3];
sx q[3];
rz(-1.310077) q[3];
sx q[3];
rz(2.6826912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

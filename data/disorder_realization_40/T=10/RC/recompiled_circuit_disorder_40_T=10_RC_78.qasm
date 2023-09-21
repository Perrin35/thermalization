OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6089132) q[0];
sx q[0];
rz(-0.37663868) q[0];
sx q[0];
rz(-3.0298046) q[0];
rz(1.6821661) q[1];
sx q[1];
rz(-1.4844126) q[1];
sx q[1];
rz(2.9878374) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4610447) q[0];
sx q[0];
rz(-1.8616315) q[0];
sx q[0];
rz(0.44797795) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2674238) q[2];
sx q[2];
rz(-2.8890925) q[2];
sx q[2];
rz(1.4361824) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7003324) q[1];
sx q[1];
rz(-1.2520737) q[1];
sx q[1];
rz(3.1096427) q[1];
x q[2];
rz(0.33711707) q[3];
sx q[3];
rz(-1.1532591) q[3];
sx q[3];
rz(1.9330213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6979606) q[2];
sx q[2];
rz(-1.4322832) q[2];
sx q[2];
rz(1.704818) q[2];
rz(0.73389655) q[3];
sx q[3];
rz(-1.5489483) q[3];
sx q[3];
rz(-2.6255887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2519418) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(-0.92457986) q[0];
rz(0.997116) q[1];
sx q[1];
rz(-0.50874248) q[1];
sx q[1];
rz(1.3234214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2572718) q[0];
sx q[0];
rz(-1.1182251) q[0];
sx q[0];
rz(-1.0898468) q[0];
x q[1];
rz(2.7919865) q[2];
sx q[2];
rz(-1.6125624) q[2];
sx q[2];
rz(-1.2431527) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4393649) q[1];
sx q[1];
rz(-2.339748) q[1];
sx q[1];
rz(2.761809) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0074535) q[3];
sx q[3];
rz(-0.65348071) q[3];
sx q[3];
rz(2.7622278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9374342) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(2.3834174) q[2];
rz(-0.6289064) q[3];
sx q[3];
rz(-2.7401676) q[3];
sx q[3];
rz(1.988407) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73308289) q[0];
sx q[0];
rz(-1.874431) q[0];
sx q[0];
rz(-2.4531903) q[0];
rz(3.0738661) q[1];
sx q[1];
rz(-1.3893145) q[1];
sx q[1];
rz(-2.6115131) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.006223) q[0];
sx q[0];
rz(-0.3100937) q[0];
sx q[0];
rz(-1.2139411) q[0];
rz(-0.42874239) q[2];
sx q[2];
rz(-1.0849761) q[2];
sx q[2];
rz(-0.57810099) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0878151) q[1];
sx q[1];
rz(-0.74790819) q[1];
sx q[1];
rz(-1.0653711) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3313053) q[3];
sx q[3];
rz(-1.665984) q[3];
sx q[3];
rz(-2.5483607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80785859) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(-0.90144908) q[2];
rz(0.83550134) q[3];
sx q[3];
rz(-1.6136026) q[3];
sx q[3];
rz(1.3114595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.19105844) q[0];
sx q[0];
rz(-1.5427417) q[0];
sx q[0];
rz(-2.3572671) q[0];
rz(-3.0803608) q[1];
sx q[1];
rz(-0.71413723) q[1];
sx q[1];
rz(0.13664666) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9204606) q[0];
sx q[0];
rz(-1.7014628) q[0];
sx q[0];
rz(1.0871825) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8881477) q[2];
sx q[2];
rz(-2.0145406) q[2];
sx q[2];
rz(-1.2197942) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9347958) q[1];
sx q[1];
rz(-0.072428457) q[1];
sx q[1];
rz(-1.2345242) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0487719) q[3];
sx q[3];
rz(-2.1448359) q[3];
sx q[3];
rz(-2.1294347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0121997) q[2];
sx q[2];
rz(-2.2012074) q[2];
sx q[2];
rz(2.5811035) q[2];
rz(-3.1292606) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(2.0509317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7085768) q[0];
sx q[0];
rz(-0.60537678) q[0];
sx q[0];
rz(-2.3204455) q[0];
rz(0.87617809) q[1];
sx q[1];
rz(-0.89996243) q[1];
sx q[1];
rz(1.7339773) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8713184) q[0];
sx q[0];
rz(-2.6042013) q[0];
sx q[0];
rz(0.73211615) q[0];
x q[1];
rz(-2.1528835) q[2];
sx q[2];
rz(-1.8277797) q[2];
sx q[2];
rz(-0.35494057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.059543) q[1];
sx q[1];
rz(-1.7051011) q[1];
sx q[1];
rz(2.1422269) q[1];
x q[2];
rz(-2.9191454) q[3];
sx q[3];
rz(-1.2192093) q[3];
sx q[3];
rz(-2.7001911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.451482) q[2];
sx q[2];
rz(-1.9268945) q[2];
sx q[2];
rz(-3.0991128) q[2];
rz(-2.5111607) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(-0.49155864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522488) q[0];
sx q[0];
rz(-1.2279953) q[0];
sx q[0];
rz(2.65843) q[0];
rz(-1.0522316) q[1];
sx q[1];
rz(-1.1455043) q[1];
sx q[1];
rz(2.5767456) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4120867) q[0];
sx q[0];
rz(-0.54486638) q[0];
sx q[0];
rz(1.8338404) q[0];
rz(-pi) q[1];
rz(1.738756) q[2];
sx q[2];
rz(-0.42617455) q[2];
sx q[2];
rz(-0.053152823) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5833252) q[1];
sx q[1];
rz(-2.7100483) q[1];
sx q[1];
rz(-1.4285018) q[1];
rz(-pi) q[2];
rz(-0.58638339) q[3];
sx q[3];
rz(-2.8660503) q[3];
sx q[3];
rz(1.1875718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.57006449) q[2];
sx q[2];
rz(-2.0740985) q[2];
sx q[2];
rz(-1.8035536) q[2];
rz(1.8367052) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(2.4664972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17334443) q[0];
sx q[0];
rz(-1.4302379) q[0];
sx q[0];
rz(2.5937953) q[0];
rz(0.785218) q[1];
sx q[1];
rz(-1.806587) q[1];
sx q[1];
rz(-2.8731667) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1584394) q[0];
sx q[0];
rz(-1.8162677) q[0];
sx q[0];
rz(1.4863187) q[0];
x q[1];
rz(-0.7873017) q[2];
sx q[2];
rz(-0.95677081) q[2];
sx q[2];
rz(-1.3348483) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35560289) q[1];
sx q[1];
rz(-2.5463748) q[1];
sx q[1];
rz(1.7463513) q[1];
rz(-pi) q[2];
rz(-0.3018474) q[3];
sx q[3];
rz(-2.0904623) q[3];
sx q[3];
rz(-0.93320751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.61775529) q[2];
sx q[2];
rz(-2.338151) q[2];
sx q[2];
rz(-0.8141554) q[2];
rz(-0.37627775) q[3];
sx q[3];
rz(-1.1637996) q[3];
sx q[3];
rz(3.0686839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.5557264) q[0];
sx q[0];
rz(-1.7365475) q[0];
sx q[0];
rz(-1.0193753) q[0];
rz(-2.2881919) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(2.6928435) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42496029) q[0];
sx q[0];
rz(-0.53629959) q[0];
sx q[0];
rz(1.1015571) q[0];
x q[1];
rz(2.2737695) q[2];
sx q[2];
rz(-2.7542369) q[2];
sx q[2];
rz(-0.58434904) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.031520695) q[1];
sx q[1];
rz(-1.8247461) q[1];
sx q[1];
rz(2.5961848) q[1];
rz(2.2711146) q[3];
sx q[3];
rz(-0.92910367) q[3];
sx q[3];
rz(0.93922797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86924187) q[2];
sx q[2];
rz(-1.3808455) q[2];
sx q[2];
rz(-0.31420079) q[2];
rz(-0.82434404) q[3];
sx q[3];
rz(-2.6894675) q[3];
sx q[3];
rz(2.3468988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.0712414) q[0];
sx q[0];
rz(-3.0817139) q[0];
sx q[0];
rz(-1.2605793) q[0];
rz(2.4977327) q[1];
sx q[1];
rz(-1.2327797) q[1];
sx q[1];
rz(3.1226645) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8277055) q[0];
sx q[0];
rz(-0.2642309) q[0];
sx q[0];
rz(-0.068473579) q[0];
rz(-pi) q[1];
x q[1];
rz(1.47255) q[2];
sx q[2];
rz(-1.4155404) q[2];
sx q[2];
rz(1.9948024) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0156118) q[1];
sx q[1];
rz(-0.19128448) q[1];
sx q[1];
rz(0.26856883) q[1];
rz(2.4631259) q[3];
sx q[3];
rz(-0.36558357) q[3];
sx q[3];
rz(1.6775223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54946047) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(2.4712759) q[2];
rz(-0.51236764) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5877514) q[0];
sx q[0];
rz(-1.2441664) q[0];
sx q[0];
rz(-2.642139) q[0];
rz(-1.5669426) q[1];
sx q[1];
rz(-0.27856871) q[1];
sx q[1];
rz(2.0589028) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1854541) q[0];
sx q[0];
rz(-0.98422613) q[0];
sx q[0];
rz(1.5120718) q[0];
rz(-2.322299) q[2];
sx q[2];
rz(-1.594992) q[2];
sx q[2];
rz(-2.1404612) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5233366) q[1];
sx q[1];
rz(-2.5976631) q[1];
sx q[1];
rz(0.036332794) q[1];
rz(-pi) q[2];
rz(2.4389078) q[3];
sx q[3];
rz(-0.48326884) q[3];
sx q[3];
rz(-1.0815222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3988951) q[2];
sx q[2];
rz(-1.9766786) q[2];
sx q[2];
rz(-0.51188525) q[2];
rz(0.39294696) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(-1.0569364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95505161) q[0];
sx q[0];
rz(-1.825009) q[0];
sx q[0];
rz(0.64074989) q[0];
rz(0.74116771) q[1];
sx q[1];
rz(-0.82294958) q[1];
sx q[1];
rz(-0.23946147) q[1];
rz(-2.7411186) q[2];
sx q[2];
rz(-0.99292143) q[2];
sx q[2];
rz(-1.6510486) q[2];
rz(0.18556553) q[3];
sx q[3];
rz(-0.75510988) q[3];
sx q[3];
rz(1.818944) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
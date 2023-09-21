OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5326795) q[0];
sx q[0];
rz(-2.764954) q[0];
sx q[0];
rz(-0.11178804) q[0];
rz(-1.4594266) q[1];
sx q[1];
rz(-1.6571801) q[1];
sx q[1];
rz(-2.9878374) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35225866) q[0];
sx q[0];
rz(-0.52868045) q[0];
sx q[0];
rz(-0.60469158) q[0];
rz(-1.329374) q[2];
sx q[2];
rz(-1.6454988) q[2];
sx q[2];
rz(0.15969294) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.598688) q[1];
sx q[1];
rz(-2.821327) q[1];
sx q[1];
rz(-1.6673052) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2114803) q[3];
sx q[3];
rz(-2.6112587) q[3];
sx q[3];
rz(1.2202642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.443632) q[2];
sx q[2];
rz(-1.7093095) q[2];
sx q[2];
rz(1.4367746) q[2];
rz(-2.4076961) q[3];
sx q[3];
rz(-1.5926444) q[3];
sx q[3];
rz(-0.51600391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88965082) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(-0.92457986) q[0];
rz(-2.1444767) q[1];
sx q[1];
rz(-0.50874248) q[1];
sx q[1];
rz(-1.8181713) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2307777) q[0];
sx q[0];
rz(-1.1417023) q[0];
sx q[0];
rz(2.6399415) q[0];
rz(1.5263444) q[2];
sx q[2];
rz(-1.2215081) q[2];
sx q[2];
rz(0.312422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1393226) q[1];
sx q[1];
rz(-1.3011258) q[1];
sx q[1];
rz(-0.76489277) q[1];
rz(-2.1341392) q[3];
sx q[3];
rz(-0.65348071) q[3];
sx q[3];
rz(2.7622278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.20415846) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(-2.3834174) q[2];
rz(2.5126863) q[3];
sx q[3];
rz(-0.40142504) q[3];
sx q[3];
rz(-1.988407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4085098) q[0];
sx q[0];
rz(-1.2671616) q[0];
sx q[0];
rz(2.4531903) q[0];
rz(3.0738661) q[1];
sx q[1];
rz(-1.7522782) q[1];
sx q[1];
rz(2.6115131) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5085174) q[0];
sx q[0];
rz(-1.2808262) q[0];
sx q[0];
rz(3.0301208) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90383756) q[2];
sx q[2];
rz(-2.505216) q[2];
sx q[2];
rz(0.19665502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5487919) q[1];
sx q[1];
rz(-0.93344102) q[1];
sx q[1];
rz(-0.42216502) q[1];
rz(-pi) q[2];
rz(1.8102874) q[3];
sx q[3];
rz(-1.4756087) q[3];
sx q[3];
rz(-2.5483607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3337341) q[2];
sx q[2];
rz(-3.1299751) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9505342) q[0];
sx q[0];
rz(-1.5427417) q[0];
sx q[0];
rz(0.78432551) q[0];
rz(-3.0803608) q[1];
sx q[1];
rz(-0.71413723) q[1];
sx q[1];
rz(0.13664666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89350677) q[0];
sx q[0];
rz(-0.49960217) q[0];
sx q[0];
rz(1.8462371) q[0];
rz(-pi) q[1];
x q[1];
rz(1.253445) q[2];
sx q[2];
rz(-1.1270521) q[2];
sx q[2];
rz(1.9217984) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54388753) q[1];
sx q[1];
rz(-1.6391616) q[1];
sx q[1];
rz(0.023936546) q[1];
rz(-1.7131545) q[3];
sx q[3];
rz(-2.5609303) q[3];
sx q[3];
rz(2.2992087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0121997) q[2];
sx q[2];
rz(-2.2012074) q[2];
sx q[2];
rz(-0.56048918) q[2];
rz(-0.012332049) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(-2.0509317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43301582) q[0];
sx q[0];
rz(-2.5362159) q[0];
sx q[0];
rz(-0.82114712) q[0];
rz(-2.2654146) q[1];
sx q[1];
rz(-2.2416302) q[1];
sx q[1];
rz(-1.7339773) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8713184) q[0];
sx q[0];
rz(-0.5373913) q[0];
sx q[0];
rz(2.4094765) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1248963) q[2];
sx q[2];
rz(-2.5113528) q[2];
sx q[2];
rz(1.5843887) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.059543) q[1];
sx q[1];
rz(-1.4364916) q[1];
sx q[1];
rz(-0.99936578) q[1];
x q[2];
rz(-0.22244723) q[3];
sx q[3];
rz(-1.9223833) q[3];
sx q[3];
rz(-2.7001911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.451482) q[2];
sx q[2];
rz(-1.9268945) q[2];
sx q[2];
rz(0.042479854) q[2];
rz(0.63043198) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522488) q[0];
sx q[0];
rz(-1.9135973) q[0];
sx q[0];
rz(-2.65843) q[0];
rz(-1.0522316) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(0.56484708) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4120867) q[0];
sx q[0];
rz(-0.54486638) q[0];
sx q[0];
rz(-1.3077523) q[0];
x q[1];
rz(-1.4028366) q[2];
sx q[2];
rz(-0.42617455) q[2];
sx q[2];
rz(3.0884398) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55826742) q[1];
sx q[1];
rz(-0.43154432) q[1];
sx q[1];
rz(1.4285018) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7259898) q[3];
sx q[3];
rz(-1.7994013) q[3];
sx q[3];
rz(-0.58333635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5715282) q[2];
sx q[2];
rz(-2.0740985) q[2];
sx q[2];
rz(-1.8035536) q[2];
rz(-1.8367052) q[3];
sx q[3];
rz(-1.0675425) q[3];
sx q[3];
rz(0.6750955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17334443) q[0];
sx q[0];
rz(-1.7113547) q[0];
sx q[0];
rz(-0.54779732) q[0];
rz(-0.785218) q[1];
sx q[1];
rz(-1.3350057) q[1];
sx q[1];
rz(0.26842591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1584394) q[0];
sx q[0];
rz(-1.8162677) q[0];
sx q[0];
rz(-1.655274) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3556049) q[2];
sx q[2];
rz(-2.18835) q[2];
sx q[2];
rz(2.3812889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.06936) q[1];
sx q[1];
rz(-1.4727122) q[1];
sx q[1];
rz(-0.98274883) q[1];
rz(-1.0309585) q[3];
sx q[3];
rz(-1.8317878) q[3];
sx q[3];
rz(-2.6574082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5238374) q[2];
sx q[2];
rz(-0.80344168) q[2];
sx q[2];
rz(2.3274373) q[2];
rz(2.7653149) q[3];
sx q[3];
rz(-1.1637996) q[3];
sx q[3];
rz(-0.072908727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58586621) q[0];
sx q[0];
rz(-1.4050452) q[0];
sx q[0];
rz(2.1222173) q[0];
rz(2.2881919) q[1];
sx q[1];
rz(-1.9995721) q[1];
sx q[1];
rz(2.6928435) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9578581) q[0];
sx q[0];
rz(-2.0439889) q[0];
sx q[0];
rz(-0.26259043) q[0];
x q[1];
rz(2.2737695) q[2];
sx q[2];
rz(-0.38735577) q[2];
sx q[2];
rz(0.58434904) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.031520695) q[1];
sx q[1];
rz(-1.8247461) q[1];
sx q[1];
rz(-2.5961848) q[1];
rz(-0.7738503) q[3];
sx q[3];
rz(-2.1132831) q[3];
sx q[3];
rz(1.098793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2723508) q[2];
sx q[2];
rz(-1.3808455) q[2];
sx q[2];
rz(2.8273919) q[2];
rz(-2.3172486) q[3];
sx q[3];
rz(-0.45212513) q[3];
sx q[3];
rz(-0.79469386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0712414) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(-1.8810133) q[0];
rz(0.64385995) q[1];
sx q[1];
rz(-1.2327797) q[1];
sx q[1];
rz(-3.1226645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32301329) q[0];
sx q[0];
rz(-1.5529263) q[0];
sx q[0];
rz(-0.26364003) q[0];
x q[1];
rz(2.9855965) q[2];
sx q[2];
rz(-1.4737355) q[2];
sx q[2];
rz(2.7023466) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4328806) q[1];
sx q[1];
rz(-1.6212665) q[1];
sx q[1];
rz(2.9570079) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8065679) q[3];
sx q[3];
rz(-1.2887495) q[3];
sx q[3];
rz(-0.75197938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54946047) q[2];
sx q[2];
rz(-2.7533054) q[2];
sx q[2];
rz(-2.4712759) q[2];
rz(-0.51236764) q[3];
sx q[3];
rz(-1.3918326) q[3];
sx q[3];
rz(-1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5877514) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(-0.49945369) q[0];
rz(-1.5746501) q[1];
sx q[1];
rz(-0.27856871) q[1];
sx q[1];
rz(-2.0589028) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1854541) q[0];
sx q[0];
rz(-2.1573665) q[0];
sx q[0];
rz(-1.5120718) q[0];
rz(-2.322299) q[2];
sx q[2];
rz(-1.5466006) q[2];
sx q[2];
rz(2.1404612) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.618256) q[1];
sx q[1];
rz(-0.54392951) q[1];
sx q[1];
rz(-0.036332794) q[1];
rz(-pi) q[2];
rz(-0.38090221) q[3];
sx q[3];
rz(-1.2657832) q[3];
sx q[3];
rz(-2.0088793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7426976) q[2];
sx q[2];
rz(-1.9766786) q[2];
sx q[2];
rz(-2.6297074) q[2];
rz(-2.7486457) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(2.0846562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95505161) q[0];
sx q[0];
rz(-1.825009) q[0];
sx q[0];
rz(0.64074989) q[0];
rz(-2.4004249) q[1];
sx q[1];
rz(-0.82294958) q[1];
sx q[1];
rz(-0.23946147) q[1];
rz(2.1096061) q[2];
sx q[2];
rz(-2.451755) q[2];
sx q[2];
rz(-2.3103726) q[2];
rz(-2.3951204) q[3];
sx q[3];
rz(-1.4440047) q[3];
sx q[3];
rz(0.11228893) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
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
rz(1.7684608) q[0];
sx q[0];
rz(-1.7113577) q[0];
sx q[0];
rz(-2.2447684) q[0];
rz(1.3431909) q[1];
sx q[1];
rz(-0.78044432) q[1];
sx q[1];
rz(1.1358383) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61124046) q[0];
sx q[0];
rz(-1.4119524) q[0];
sx q[0];
rz(0.071119951) q[0];
rz(-pi) q[1];
rz(-0.99857371) q[2];
sx q[2];
rz(-2.3731447) q[2];
sx q[2];
rz(-0.060800663) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9998127) q[1];
sx q[1];
rz(-2.2408812) q[1];
sx q[1];
rz(-0.72061909) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4189931) q[3];
sx q[3];
rz(-0.54908592) q[3];
sx q[3];
rz(0.3749519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4844126) q[2];
sx q[2];
rz(-1.5766532) q[2];
sx q[2];
rz(2.8489992) q[2];
rz(1.6491133) q[3];
sx q[3];
rz(-1.7238659) q[3];
sx q[3];
rz(-2.9610146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0507616) q[0];
sx q[0];
rz(-1.2240336) q[0];
sx q[0];
rz(-2.5378788) q[0];
rz(2.1049888) q[1];
sx q[1];
rz(-2.7991468) q[1];
sx q[1];
rz(2.2110914) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0261352) q[0];
sx q[0];
rz(-1.8107426) q[0];
sx q[0];
rz(-1.351993) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23463191) q[2];
sx q[2];
rz(-0.57387239) q[2];
sx q[2];
rz(-3.0374276) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89255561) q[1];
sx q[1];
rz(-2.3904789) q[1];
sx q[1];
rz(-0.19601258) q[1];
rz(-pi) q[2];
rz(-0.33035461) q[3];
sx q[3];
rz(-2.1393015) q[3];
sx q[3];
rz(-1.1363883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36640627) q[2];
sx q[2];
rz(-1.2440871) q[2];
sx q[2];
rz(1.5001971) q[2];
rz(-0.90534798) q[3];
sx q[3];
rz(-1.2572181) q[3];
sx q[3];
rz(-0.94178981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6833078) q[0];
sx q[0];
rz(-1.6907121) q[0];
sx q[0];
rz(2.2990551) q[0];
rz(1.204035) q[1];
sx q[1];
rz(-1.8572109) q[1];
sx q[1];
rz(-1.1010928) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8592671) q[0];
sx q[0];
rz(-2.8045033) q[0];
sx q[0];
rz(1.8021148) q[0];
x q[1];
rz(-0.21761009) q[2];
sx q[2];
rz(-2.3869385) q[2];
sx q[2];
rz(-0.77384743) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2223375) q[1];
sx q[1];
rz(-0.65523096) q[1];
sx q[1];
rz(-2.3690556) q[1];
rz(-2.7277953) q[3];
sx q[3];
rz(-0.99552508) q[3];
sx q[3];
rz(0.12518342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.6847685) q[2];
sx q[2];
rz(-2.5825239) q[2];
sx q[2];
rz(3.0341201) q[2];
rz(-2.0939854) q[3];
sx q[3];
rz(-1.8045629) q[3];
sx q[3];
rz(1.5260772) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90077129) q[0];
sx q[0];
rz(-2.695684) q[0];
sx q[0];
rz(1.6915503) q[0];
rz(-2.5598473) q[1];
sx q[1];
rz(-0.83256871) q[1];
sx q[1];
rz(1.1004826) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2028401) q[0];
sx q[0];
rz(-2.0838408) q[0];
sx q[0];
rz(1.130453) q[0];
rz(-pi) q[1];
rz(-2.5610177) q[2];
sx q[2];
rz(-2.5126804) q[2];
sx q[2];
rz(3.1286503) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2128513) q[1];
sx q[1];
rz(-0.87777661) q[1];
sx q[1];
rz(-1.5726967) q[1];
x q[2];
rz(2.2868312) q[3];
sx q[3];
rz(-1.3472036) q[3];
sx q[3];
rz(1.8098698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.45740286) q[2];
sx q[2];
rz(-2.0997014) q[2];
sx q[2];
rz(-0.602496) q[2];
rz(1.4495133) q[3];
sx q[3];
rz(-1.7159729) q[3];
sx q[3];
rz(-0.95363936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7736037) q[0];
sx q[0];
rz(-1.8743176) q[0];
sx q[0];
rz(-1.8220655) q[0];
rz(-0.52917448) q[1];
sx q[1];
rz(-1.2472943) q[1];
sx q[1];
rz(0.41437638) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0841211) q[0];
sx q[0];
rz(-1.9638279) q[0];
sx q[0];
rz(-1.074206) q[0];
rz(-0.90535449) q[2];
sx q[2];
rz(-2.2756409) q[2];
sx q[2];
rz(-0.030817835) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.16769174) q[1];
sx q[1];
rz(-1.138559) q[1];
sx q[1];
rz(0.28829379) q[1];
x q[2];
rz(0.46196725) q[3];
sx q[3];
rz(-0.35502343) q[3];
sx q[3];
rz(-2.9649761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6759912) q[2];
sx q[2];
rz(-2.4310515) q[2];
sx q[2];
rz(0.65593925) q[2];
rz(1.9115619) q[3];
sx q[3];
rz(-2.0345119) q[3];
sx q[3];
rz(0.51903498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32992724) q[0];
sx q[0];
rz(-2.1274607) q[0];
sx q[0];
rz(-2.2299679) q[0];
rz(1.7650334) q[1];
sx q[1];
rz(-2.7477317) q[1];
sx q[1];
rz(0.54042655) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9003284) q[0];
sx q[0];
rz(-2.6200128) q[0];
sx q[0];
rz(-2.5425884) q[0];
rz(-pi) q[1];
rz(1.592268) q[2];
sx q[2];
rz(-1.2228857) q[2];
sx q[2];
rz(-2.3644115) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6669861) q[1];
sx q[1];
rz(-1.9691393) q[1];
sx q[1];
rz(1.6380659) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8761108) q[3];
sx q[3];
rz(-2.1830227) q[3];
sx q[3];
rz(1.4683362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2889169) q[2];
sx q[2];
rz(-1.6923075) q[2];
sx q[2];
rz(1.1360315) q[2];
rz(1.8063258) q[3];
sx q[3];
rz(-1.5832333) q[3];
sx q[3];
rz(-3.0291338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54953185) q[0];
sx q[0];
rz(-2.9887152) q[0];
sx q[0];
rz(0.98912799) q[0];
rz(1.9391183) q[1];
sx q[1];
rz(-1.4710642) q[1];
sx q[1];
rz(0.30754009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3002994) q[0];
sx q[0];
rz(-1.3430149) q[0];
sx q[0];
rz(-2.862597) q[0];
rz(-1.9498242) q[2];
sx q[2];
rz(-2.267365) q[2];
sx q[2];
rz(-1.7988169) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.16368281) q[1];
sx q[1];
rz(-1.3505757) q[1];
sx q[1];
rz(1.5942595) q[1];
x q[2];
rz(-0.6102517) q[3];
sx q[3];
rz(-2.0008127) q[3];
sx q[3];
rz(0.56320923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62693739) q[2];
sx q[2];
rz(-0.81521002) q[2];
sx q[2];
rz(1.6132272) q[2];
rz(-0.79742399) q[3];
sx q[3];
rz(-1.593109) q[3];
sx q[3];
rz(-3.048866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02354694) q[0];
sx q[0];
rz(-2.4674456) q[0];
sx q[0];
rz(2.4349037) q[0];
rz(2.8382909) q[1];
sx q[1];
rz(-1.3677596) q[1];
sx q[1];
rz(1.2695262) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72699914) q[0];
sx q[0];
rz(-1.4395066) q[0];
sx q[0];
rz(-2.8069228) q[0];
rz(-pi) q[1];
rz(-1.3236155) q[2];
sx q[2];
rz(-2.0529785) q[2];
sx q[2];
rz(-3.1076933) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8757585) q[1];
sx q[1];
rz(-2.7175329) q[1];
sx q[1];
rz(1.1316749) q[1];
x q[2];
rz(-1.4702099) q[3];
sx q[3];
rz(-1.0142418) q[3];
sx q[3];
rz(-0.33410698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.068744008) q[2];
sx q[2];
rz(-2.7903283) q[2];
sx q[2];
rz(-2.6212485) q[2];
rz(-1.8955692) q[3];
sx q[3];
rz(-1.7249853) q[3];
sx q[3];
rz(-2.0791159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50616566) q[0];
sx q[0];
rz(-1.4816254) q[0];
sx q[0];
rz(0.51314276) q[0];
rz(-0.48140934) q[1];
sx q[1];
rz(-0.82999271) q[1];
sx q[1];
rz(-2.381348) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82448563) q[0];
sx q[0];
rz(-3.1104507) q[0];
sx q[0];
rz(-1.1481773) q[0];
rz(-pi) q[1];
rz(-1.3584601) q[2];
sx q[2];
rz(-0.80695769) q[2];
sx q[2];
rz(-2.2495861) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1380242) q[1];
sx q[1];
rz(-1.6449182) q[1];
sx q[1];
rz(-2.4155492) q[1];
rz(-0.097858345) q[3];
sx q[3];
rz(-1.9966085) q[3];
sx q[3];
rz(-1.875836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61797872) q[2];
sx q[2];
rz(-0.60680497) q[2];
sx q[2];
rz(-0.65110171) q[2];
rz(-1.2111604) q[3];
sx q[3];
rz(-1.159472) q[3];
sx q[3];
rz(2.7291362) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3245658) q[0];
sx q[0];
rz(-0.71837076) q[0];
sx q[0];
rz(2.4107362) q[0];
rz(-2.6279347) q[1];
sx q[1];
rz(-0.81413236) q[1];
sx q[1];
rz(-2.6252852) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56591153) q[0];
sx q[0];
rz(-1.9464419) q[0];
sx q[0];
rz(-0.89976914) q[0];
rz(-2.7514156) q[2];
sx q[2];
rz(-1.220229) q[2];
sx q[2];
rz(-2.1250909) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.67077209) q[1];
sx q[1];
rz(-1.9665475) q[1];
sx q[1];
rz(-0.20352965) q[1];
rz(-pi) q[2];
rz(-1.794471) q[3];
sx q[3];
rz(-1.8590611) q[3];
sx q[3];
rz(2.1065421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98661304) q[2];
sx q[2];
rz(-1.2519138) q[2];
sx q[2];
rz(2.4775141) q[2];
rz(-2.1443071) q[3];
sx q[3];
rz(-1.7136796) q[3];
sx q[3];
rz(-0.47521457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.3288788) q[0];
sx q[0];
rz(-1.4007778) q[0];
sx q[0];
rz(3.1123871) q[0];
rz(-3.008814) q[1];
sx q[1];
rz(-1.4135502) q[1];
sx q[1];
rz(2.2739364) q[1];
rz(1.221413) q[2];
sx q[2];
rz(-1.0950015) q[2];
sx q[2];
rz(-2.0910043) q[2];
rz(2.7326591) q[3];
sx q[3];
rz(-1.144617) q[3];
sx q[3];
rz(-1.495818) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

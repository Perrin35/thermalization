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
rz(-3.0566293) q[0];
sx q[0];
rz(-0.30240107) q[0];
sx q[0];
rz(0.032057134) q[0];
rz(3.1337466) q[1];
sx q[1];
rz(-2.6377331) q[1];
sx q[1];
rz(-2.388968) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40850484) q[0];
sx q[0];
rz(-1.2512387) q[0];
sx q[0];
rz(-1.7825141) q[0];
x q[1];
rz(1.9448124) q[2];
sx q[2];
rz(-1.2311282) q[2];
sx q[2];
rz(2.938478) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9702985) q[1];
sx q[1];
rz(-0.42810218) q[1];
sx q[1];
rz(2.5004205) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2129477) q[3];
sx q[3];
rz(-1.3697624) q[3];
sx q[3];
rz(1.4154289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96693119) q[2];
sx q[2];
rz(-1.175468) q[2];
sx q[2];
rz(-0.82007972) q[2];
rz(0.44102937) q[3];
sx q[3];
rz(-2.0377772) q[3];
sx q[3];
rz(-0.57344121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1287307) q[0];
sx q[0];
rz(-2.2282889) q[0];
sx q[0];
rz(-0.41020694) q[0];
rz(2.7606616) q[1];
sx q[1];
rz(-1.610264) q[1];
sx q[1];
rz(1.7832696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5429496) q[0];
sx q[0];
rz(-2.1457167) q[0];
sx q[0];
rz(0.47358124) q[0];
rz(-pi) q[1];
rz(1.0274506) q[2];
sx q[2];
rz(-2.0110235) q[2];
sx q[2];
rz(1.5433951) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7100923) q[1];
sx q[1];
rz(-0.18431252) q[1];
sx q[1];
rz(1.9955562) q[1];
rz(-pi) q[2];
rz(-0.15492691) q[3];
sx q[3];
rz(-2.2048031) q[3];
sx q[3];
rz(2.4661494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7831948) q[2];
sx q[2];
rz(-2.7326549) q[2];
sx q[2];
rz(-2.6386293) q[2];
rz(1.8424235) q[3];
sx q[3];
rz(-0.75853577) q[3];
sx q[3];
rz(2.9097596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7314887) q[0];
sx q[0];
rz(-2.6298611) q[0];
sx q[0];
rz(2.1657535) q[0];
rz(-0.93636912) q[1];
sx q[1];
rz(-1.5872637) q[1];
sx q[1];
rz(0.13793129) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5463211) q[0];
sx q[0];
rz(-1.0965276) q[0];
sx q[0];
rz(-2.7300937) q[0];
x q[1];
rz(1.6866272) q[2];
sx q[2];
rz(-1.1755014) q[2];
sx q[2];
rz(-2.9651053) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5974849) q[1];
sx q[1];
rz(-0.710809) q[1];
sx q[1];
rz(1.3686468) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6343752) q[3];
sx q[3];
rz(-2.6196369) q[3];
sx q[3];
rz(2.0824661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5522573) q[2];
sx q[2];
rz(-1.5312803) q[2];
sx q[2];
rz(1.4790347) q[2];
rz(-2.3432483) q[3];
sx q[3];
rz(-0.84404293) q[3];
sx q[3];
rz(3.121283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8524356) q[0];
sx q[0];
rz(-3.030179) q[0];
sx q[0];
rz(1.0668466) q[0];
rz(1.5929219) q[1];
sx q[1];
rz(-1.2499864) q[1];
sx q[1];
rz(-0.41935316) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8521312) q[0];
sx q[0];
rz(-1.7999819) q[0];
sx q[0];
rz(0.16332345) q[0];
rz(-0.52813645) q[2];
sx q[2];
rz(-1.4924876) q[2];
sx q[2];
rz(-2.9700043) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4511299) q[1];
sx q[1];
rz(-1.2138543) q[1];
sx q[1];
rz(-0.068515645) q[1];
rz(-pi) q[2];
rz(-0.97030117) q[3];
sx q[3];
rz(-0.56823778) q[3];
sx q[3];
rz(2.0214391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7833917) q[2];
sx q[2];
rz(-0.50938598) q[2];
sx q[2];
rz(-1.5436714) q[2];
rz(1.915043) q[3];
sx q[3];
rz(-1.2164601) q[3];
sx q[3];
rz(1.290087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.876038) q[0];
sx q[0];
rz(-0.72172481) q[0];
sx q[0];
rz(-2.3853886) q[0];
rz(-1.2528231) q[1];
sx q[1];
rz(-0.93106657) q[1];
sx q[1];
rz(-1.4160215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5042345) q[0];
sx q[0];
rz(-1.7880482) q[0];
sx q[0];
rz(-0.603032) q[0];
x q[1];
rz(1.6966693) q[2];
sx q[2];
rz(-1.2294266) q[2];
sx q[2];
rz(1.9370796) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6452613) q[1];
sx q[1];
rz(-2.1095536) q[1];
sx q[1];
rz(-1.5790198) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46964271) q[3];
sx q[3];
rz(-2.7386463) q[3];
sx q[3];
rz(-0.5136036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.34117928) q[2];
sx q[2];
rz(-2.4004553) q[2];
sx q[2];
rz(2.9608012) q[2];
rz(1.6527269) q[3];
sx q[3];
rz(-1.7427665) q[3];
sx q[3];
rz(1.5159877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7086696) q[0];
sx q[0];
rz(-2.0661418) q[0];
sx q[0];
rz(-2.3413626) q[0];
rz(-0.284614) q[1];
sx q[1];
rz(-1.0601284) q[1];
sx q[1];
rz(-3.0175993) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02888814) q[0];
sx q[0];
rz(-1.7005973) q[0];
sx q[0];
rz(3.1129254) q[0];
rz(-1.335786) q[2];
sx q[2];
rz(-1.628865) q[2];
sx q[2];
rz(-0.11561671) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1734487) q[1];
sx q[1];
rz(-0.3160797) q[1];
sx q[1];
rz(-0.26559243) q[1];
rz(2.769737) q[3];
sx q[3];
rz(-2.1346492) q[3];
sx q[3];
rz(1.2363889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93512145) q[2];
sx q[2];
rz(-1.7266885) q[2];
sx q[2];
rz(-0.07621152) q[2];
rz(1.291409) q[3];
sx q[3];
rz(-1.094787) q[3];
sx q[3];
rz(-0.61856234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9718219) q[0];
sx q[0];
rz(-1.562028) q[0];
sx q[0];
rz(-0.75702697) q[0];
rz(-2.3454759) q[1];
sx q[1];
rz(-1.7346104) q[1];
sx q[1];
rz(2.1597791) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1025006) q[0];
sx q[0];
rz(-1.8538686) q[0];
sx q[0];
rz(0.82067482) q[0];
rz(0.99346353) q[2];
sx q[2];
rz(-3.0138123) q[2];
sx q[2];
rz(-1.3410695) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.031739) q[1];
sx q[1];
rz(-0.64424911) q[1];
sx q[1];
rz(2.2883313) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8514093) q[3];
sx q[3];
rz(-1.7095437) q[3];
sx q[3];
rz(2.8900103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.41165274) q[2];
sx q[2];
rz(-0.81092683) q[2];
sx q[2];
rz(-0.5438424) q[2];
rz(0.020921556) q[3];
sx q[3];
rz(-1.436751) q[3];
sx q[3];
rz(-0.81834832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2127317) q[0];
sx q[0];
rz(-2.2305771) q[0];
sx q[0];
rz(2.5182356) q[0];
rz(2.8202672) q[1];
sx q[1];
rz(-2.5265381) q[1];
sx q[1];
rz(-0.26434937) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0118222) q[0];
sx q[0];
rz(-1.3589824) q[0];
sx q[0];
rz(0.84902936) q[0];
rz(1.4521763) q[2];
sx q[2];
rz(-2.8037694) q[2];
sx q[2];
rz(-1.092697) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.26786727) q[1];
sx q[1];
rz(-1.8494938) q[1];
sx q[1];
rz(0.81588094) q[1];
rz(-1.1846011) q[3];
sx q[3];
rz(-2.8004146) q[3];
sx q[3];
rz(-0.5285078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.136772) q[2];
sx q[2];
rz(-0.68237582) q[2];
sx q[2];
rz(2.6668059) q[2];
rz(0.36561203) q[3];
sx q[3];
rz(-2.6545299) q[3];
sx q[3];
rz(-0.18317187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.8356165) q[0];
sx q[0];
rz(-0.37537471) q[0];
sx q[0];
rz(2.6557652) q[0];
rz(1.0659418) q[1];
sx q[1];
rz(-1.7168047) q[1];
sx q[1];
rz(2.8071383) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1565103) q[0];
sx q[0];
rz(-1.872439) q[0];
sx q[0];
rz(1.3972052) q[0];
x q[1];
rz(2.8489001) q[2];
sx q[2];
rz(-2.6862157) q[2];
sx q[2];
rz(0.51560452) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.124467) q[1];
sx q[1];
rz(-1.4216058) q[1];
sx q[1];
rz(-0.45254947) q[1];
rz(-pi) q[2];
rz(0.45463698) q[3];
sx q[3];
rz(-2.2918923) q[3];
sx q[3];
rz(2.7434512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36755422) q[2];
sx q[2];
rz(-0.84711051) q[2];
sx q[2];
rz(0.96662194) q[2];
rz(-1.6537846) q[3];
sx q[3];
rz(-2.8130468) q[3];
sx q[3];
rz(-0.070913471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8129355) q[0];
sx q[0];
rz(-2.2870977) q[0];
sx q[0];
rz(2.8712414) q[0];
rz(1.6290172) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(2.5152452) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26890818) q[0];
sx q[0];
rz(-1.2546526) q[0];
sx q[0];
rz(2.288398) q[0];
rz(-pi) q[1];
rz(-0.93339351) q[2];
sx q[2];
rz(-1.2167756) q[2];
sx q[2];
rz(-0.4590946) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69520187) q[1];
sx q[1];
rz(-0.099041136) q[1];
sx q[1];
rz(2.7161069) q[1];
rz(1.8966497) q[3];
sx q[3];
rz(-1.2030461) q[3];
sx q[3];
rz(-2.8830584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7865929) q[2];
sx q[2];
rz(-0.91675106) q[2];
sx q[2];
rz(1.5724486) q[2];
rz(-2.3908424) q[3];
sx q[3];
rz(-0.34199491) q[3];
sx q[3];
rz(-0.87740889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.5759721) q[0];
sx q[0];
rz(-1.2979869) q[0];
sx q[0];
rz(-0.57814231) q[0];
rz(-1.7987953) q[1];
sx q[1];
rz(-2.8248351) q[1];
sx q[1];
rz(0.14541365) q[1];
rz(3.1411896) q[2];
sx q[2];
rz(-3.0192791) q[2];
sx q[2];
rz(3.064439) q[2];
rz(-0.20082898) q[3];
sx q[3];
rz(-0.26293892) q[3];
sx q[3];
rz(0.50504167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

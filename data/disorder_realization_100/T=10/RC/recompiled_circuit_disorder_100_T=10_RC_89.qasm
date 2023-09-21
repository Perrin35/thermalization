OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.202012) q[0];
sx q[0];
rz(-2.7913845) q[0];
sx q[0];
rz(0.36663088) q[0];
rz(0.86759138) q[1];
sx q[1];
rz(-2.4974513) q[1];
sx q[1];
rz(1.4555414) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5397352) q[0];
sx q[0];
rz(-0.55668133) q[0];
sx q[0];
rz(0.88168983) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8800814) q[2];
sx q[2];
rz(-0.72578428) q[2];
sx q[2];
rz(2.6860565) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.213233) q[1];
sx q[1];
rz(-1.302049) q[1];
sx q[1];
rz(2.1409722) q[1];
rz(-pi) q[2];
rz(-2.9204908) q[3];
sx q[3];
rz(-0.049115291) q[3];
sx q[3];
rz(2.5887095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.036844) q[2];
sx q[2];
rz(-1.8061183) q[2];
sx q[2];
rz(-0.39164266) q[2];
rz(0.13970217) q[3];
sx q[3];
rz(-0.67290664) q[3];
sx q[3];
rz(0.02123775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7249107) q[0];
sx q[0];
rz(-1.7371438) q[0];
sx q[0];
rz(1.8649944) q[0];
rz(-0.87031594) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(1.2044027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63277584) q[0];
sx q[0];
rz(-2.5460089) q[0];
sx q[0];
rz(2.3217208) q[0];
x q[1];
rz(-1.7877098) q[2];
sx q[2];
rz(-2.4241944) q[2];
sx q[2];
rz(-2.4153828) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8783427) q[1];
sx q[1];
rz(-1.6701356) q[1];
sx q[1];
rz(-2.7753825) q[1];
rz(-0.28537206) q[3];
sx q[3];
rz(-2.223569) q[3];
sx q[3];
rz(-2.5170346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.63699547) q[2];
sx q[2];
rz(-0.4578788) q[2];
sx q[2];
rz(-1.9223928) q[2];
rz(-0.35456625) q[3];
sx q[3];
rz(-0.48186007) q[3];
sx q[3];
rz(-1.569081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5511659) q[0];
sx q[0];
rz(-1.5274436) q[0];
sx q[0];
rz(0.61808008) q[0];
rz(-0.016050054) q[1];
sx q[1];
rz(-2.2647104) q[1];
sx q[1];
rz(1.191167) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55710775) q[0];
sx q[0];
rz(-2.7911721) q[0];
sx q[0];
rz(-2.0273897) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8148003) q[2];
sx q[2];
rz(-1.2952431) q[2];
sx q[2];
rz(-0.22305605) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1515854) q[1];
sx q[1];
rz(-1.3602435) q[1];
sx q[1];
rz(-1.7819808) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61471625) q[3];
sx q[3];
rz(-0.63496548) q[3];
sx q[3];
rz(-2.281651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1742192) q[2];
sx q[2];
rz(-2.1790395) q[2];
sx q[2];
rz(-2.1276786) q[2];
rz(1.7845456) q[3];
sx q[3];
rz(-2.3116528) q[3];
sx q[3];
rz(2.1092265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33092609) q[0];
sx q[0];
rz(-1.7174915) q[0];
sx q[0];
rz(2.9372835) q[0];
rz(-1.3775685) q[1];
sx q[1];
rz(-1.4148477) q[1];
sx q[1];
rz(-0.92299443) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5343691) q[0];
sx q[0];
rz(-1.9019706) q[0];
sx q[0];
rz(1.2481199) q[0];
rz(-pi) q[1];
rz(-1.5825908) q[2];
sx q[2];
rz(-1.6147385) q[2];
sx q[2];
rz(-0.61461385) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35034414) q[1];
sx q[1];
rz(-1.6748168) q[1];
sx q[1];
rz(-1.9700325) q[1];
rz(-3.119167) q[3];
sx q[3];
rz(-0.86655819) q[3];
sx q[3];
rz(1.9989597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6874281) q[2];
sx q[2];
rz(-1.585008) q[2];
sx q[2];
rz(-2.6320809) q[2];
rz(-0.64309684) q[3];
sx q[3];
rz(-2.0431079) q[3];
sx q[3];
rz(-0.96737868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5287857) q[0];
sx q[0];
rz(-1.9354154) q[0];
sx q[0];
rz(0.91941419) q[0];
rz(-0.98584229) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(-1.7808419) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082598535) q[0];
sx q[0];
rz(-1.7599802) q[0];
sx q[0];
rz(0.22342213) q[0];
rz(0.64702101) q[2];
sx q[2];
rz(-2.2685452) q[2];
sx q[2];
rz(2.9014212) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9779224) q[1];
sx q[1];
rz(-1.1174035) q[1];
sx q[1];
rz(-1.3400673) q[1];
x q[2];
rz(-0.97134437) q[3];
sx q[3];
rz(-1.506862) q[3];
sx q[3];
rz(1.5435227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78648606) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(-0.11486593) q[2];
rz(2.6857175) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(-1.280064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546656) q[0];
sx q[0];
rz(-1.2446612) q[0];
sx q[0];
rz(1.2448357) q[0];
rz(2.4339829) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(2.8663666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9305785) q[0];
sx q[0];
rz(-1.0761346) q[0];
sx q[0];
rz(-2.7633694) q[0];
rz(-0.92991021) q[2];
sx q[2];
rz(-1.9800131) q[2];
sx q[2];
rz(-0.90604679) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.59906193) q[1];
sx q[1];
rz(-1.2236569) q[1];
sx q[1];
rz(-0.51862006) q[1];
rz(1.4554548) q[3];
sx q[3];
rz(-0.62386688) q[3];
sx q[3];
rz(-1.3184402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.747921) q[2];
sx q[2];
rz(-1.5966406) q[2];
sx q[2];
rz(2.6829524) q[2];
rz(-0.7115055) q[3];
sx q[3];
rz(-2.4578874) q[3];
sx q[3];
rz(0.59035629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28073072) q[0];
sx q[0];
rz(-1.9545398) q[0];
sx q[0];
rz(-1.3635427) q[0];
rz(-1.7215615) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(0.71969676) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73368209) q[0];
sx q[0];
rz(-0.94607293) q[0];
sx q[0];
rz(1.2172933) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3709488) q[2];
sx q[2];
rz(-2.1241803) q[2];
sx q[2];
rz(1.9945952) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2232659) q[1];
sx q[1];
rz(-1.0263138) q[1];
sx q[1];
rz(1.5189927) q[1];
x q[2];
rz(1.8411438) q[3];
sx q[3];
rz(-2.2336604) q[3];
sx q[3];
rz(-0.28966749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.65934962) q[2];
sx q[2];
rz(-1.6094001) q[2];
sx q[2];
rz(2.4534295) q[2];
rz(-0.041953772) q[3];
sx q[3];
rz(-1.099702) q[3];
sx q[3];
rz(-2.8614614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91350895) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(2.9550609) q[0];
rz(-2.5371011) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(-1.4950745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4743487) q[0];
sx q[0];
rz(-1.7843887) q[0];
sx q[0];
rz(-2.7937907) q[0];
rz(2.9096793) q[2];
sx q[2];
rz(-2.4319318) q[2];
sx q[2];
rz(2.2556925) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.013156803) q[1];
sx q[1];
rz(-1.4669347) q[1];
sx q[1];
rz(-0.26058148) q[1];
rz(-pi) q[2];
rz(-0.64077326) q[3];
sx q[3];
rz(-0.91859222) q[3];
sx q[3];
rz(2.3280502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.24215332) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(3.126826) q[2];
rz(1.8642558) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(1.4130672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9799141) q[0];
sx q[0];
rz(-2.4665687) q[0];
sx q[0];
rz(-2.263608) q[0];
rz(0.19628482) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(1.7810129) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36688766) q[0];
sx q[0];
rz(-2.5686712) q[0];
sx q[0];
rz(-1.5013055) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0948823) q[2];
sx q[2];
rz(-0.94418664) q[2];
sx q[2];
rz(-1.4057297) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.18371836) q[1];
sx q[1];
rz(-2.6596333) q[1];
sx q[1];
rz(-0.79224371) q[1];
rz(-pi) q[2];
rz(-1.2336041) q[3];
sx q[3];
rz(-2.7003151) q[3];
sx q[3];
rz(2.6233167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.62375162) q[2];
sx q[2];
rz(-1.5072284) q[2];
sx q[2];
rz(2.2488135) q[2];
rz(0.039285224) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(2.3915496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2980625) q[0];
sx q[0];
rz(-0.003304464) q[0];
sx q[0];
rz(-3.0950586) q[0];
rz(-2.2325366) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(0.7199026) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7881308) q[0];
sx q[0];
rz(-0.25665584) q[0];
sx q[0];
rz(1.8887595) q[0];
x q[1];
rz(0.017059762) q[2];
sx q[2];
rz(-0.3088726) q[2];
sx q[2];
rz(-2.942254) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9626999) q[1];
sx q[1];
rz(-1.3239412) q[1];
sx q[1];
rz(-0.015951338) q[1];
rz(0.7474483) q[3];
sx q[3];
rz(-0.8851074) q[3];
sx q[3];
rz(1.3626584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4204734) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(0.026467888) q[2];
rz(-1.8376393) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3532886) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(-0.2164671) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(-0.36454501) q[2];
sx q[2];
rz(-2.400763) q[2];
sx q[2];
rz(2.7999067) q[2];
rz(-1.5218232) q[3];
sx q[3];
rz(-0.58717849) q[3];
sx q[3];
rz(-1.0321454) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
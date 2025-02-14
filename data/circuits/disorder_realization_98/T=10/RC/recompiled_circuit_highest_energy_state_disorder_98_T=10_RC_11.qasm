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
rz(0.75818169) q[0];
sx q[0];
rz(-2.7739006) q[0];
sx q[0];
rz(-2.0453069) q[0];
rz(0.81049377) q[1];
sx q[1];
rz(-0.23263045) q[1];
sx q[1];
rz(1.1059603) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152412) q[0];
sx q[0];
rz(-2.5509074) q[0];
sx q[0];
rz(-0.75449852) q[0];
rz(-pi) q[1];
rz(2.8739086) q[2];
sx q[2];
rz(-1.6639198) q[2];
sx q[2];
rz(2.6665319) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3333862) q[1];
sx q[1];
rz(-1.574834) q[1];
sx q[1];
rz(1.5334849) q[1];
rz(-pi) q[2];
rz(-1.5638142) q[3];
sx q[3];
rz(-1.2864224) q[3];
sx q[3];
rz(-2.6740536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37630633) q[2];
sx q[2];
rz(-0.47069612) q[2];
sx q[2];
rz(-0.37962309) q[2];
rz(-1.6342573) q[3];
sx q[3];
rz(-2.0416656) q[3];
sx q[3];
rz(-0.96989337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4955502) q[0];
sx q[0];
rz(-0.33691418) q[0];
sx q[0];
rz(2.2406793) q[0];
rz(-0.13149978) q[1];
sx q[1];
rz(-1.9410746) q[1];
sx q[1];
rz(-2.6995755) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95831224) q[0];
sx q[0];
rz(-1.4889476) q[0];
sx q[0];
rz(-0.32572066) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9346048) q[2];
sx q[2];
rz(-0.0047193165) q[2];
sx q[2];
rz(-2.8568792) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1902679) q[1];
sx q[1];
rz(-2.3735078) q[1];
sx q[1];
rz(2.4332341) q[1];
rz(-pi) q[2];
rz(-2.8825661) q[3];
sx q[3];
rz(-0.4126848) q[3];
sx q[3];
rz(-0.55283062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3758731) q[2];
sx q[2];
rz(-1.4292382) q[2];
sx q[2];
rz(0.6089375) q[2];
rz(-0.40306148) q[3];
sx q[3];
rz(-1.9654704) q[3];
sx q[3];
rz(2.6398931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8159863) q[0];
sx q[0];
rz(-0.41929647) q[0];
sx q[0];
rz(2.0035279) q[0];
rz(1.552938) q[1];
sx q[1];
rz(-0.76858968) q[1];
sx q[1];
rz(2.3345711) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53483665) q[0];
sx q[0];
rz(-1.5436158) q[0];
sx q[0];
rz(-3.0423714) q[0];
rz(-pi) q[1];
rz(-1.4078232) q[2];
sx q[2];
rz(-2.0103243) q[2];
sx q[2];
rz(-1.3541612) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5182636) q[1];
sx q[1];
rz(-0.22364549) q[1];
sx q[1];
rz(2.5751051) q[1];
rz(-0.28636264) q[3];
sx q[3];
rz(-1.935528) q[3];
sx q[3];
rz(1.3160365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.011586172) q[2];
sx q[2];
rz(-0.61379543) q[2];
sx q[2];
rz(1.9278795) q[2];
rz(-0.064149292) q[3];
sx q[3];
rz(-1.4870653) q[3];
sx q[3];
rz(-3.1411662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59378687) q[0];
sx q[0];
rz(-1.2603899) q[0];
sx q[0];
rz(0.88687819) q[0];
rz(-0.15874323) q[1];
sx q[1];
rz(-0.49247772) q[1];
sx q[1];
rz(-1.8633206) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0905545) q[0];
sx q[0];
rz(-2.0060385) q[0];
sx q[0];
rz(2.3440775) q[0];
x q[1];
rz(2.5910225) q[2];
sx q[2];
rz(-1.9801557) q[2];
sx q[2];
rz(-2.0602202) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5090829) q[1];
sx q[1];
rz(-0.69493587) q[1];
sx q[1];
rz(-1.664723) q[1];
rz(2.9904891) q[3];
sx q[3];
rz(-2.1204815) q[3];
sx q[3];
rz(-2.0696039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0025803) q[2];
sx q[2];
rz(-0.85550344) q[2];
sx q[2];
rz(0.95376897) q[2];
rz(-1.194713) q[3];
sx q[3];
rz(-2.298893) q[3];
sx q[3];
rz(0.4755303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82763982) q[0];
sx q[0];
rz(-0.71641818) q[0];
sx q[0];
rz(2.5415976) q[0];
rz(-2.4586239) q[1];
sx q[1];
rz(-0.78790793) q[1];
sx q[1];
rz(0.33040985) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67067388) q[0];
sx q[0];
rz(-2.3944811) q[0];
sx q[0];
rz(-1.8404191) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1938964) q[2];
sx q[2];
rz(-0.27505878) q[2];
sx q[2];
rz(0.81240053) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.85530419) q[1];
sx q[1];
rz(-2.0518612) q[1];
sx q[1];
rz(-0.21800133) q[1];
x q[2];
rz(-0.41466576) q[3];
sx q[3];
rz(-2.2288481) q[3];
sx q[3];
rz(-2.0882704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7052475) q[2];
sx q[2];
rz(-1.0785495) q[2];
sx q[2];
rz(0.62502965) q[2];
rz(-2.446512) q[3];
sx q[3];
rz(-1.2806226) q[3];
sx q[3];
rz(-0.052791031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92796749) q[0];
sx q[0];
rz(-1.9897505) q[0];
sx q[0];
rz(1.7684162) q[0];
rz(1.3881418) q[1];
sx q[1];
rz(-2.2699247) q[1];
sx q[1];
rz(-1.53481) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2508188) q[0];
sx q[0];
rz(-1.5355331) q[0];
sx q[0];
rz(-1.718344) q[0];
rz(-pi) q[1];
rz(-1.0523782) q[2];
sx q[2];
rz(-0.3568584) q[2];
sx q[2];
rz(-1.269358) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36027137) q[1];
sx q[1];
rz(-2.2151673) q[1];
sx q[1];
rz(-2.086198) q[1];
x q[2];
rz(3.0032773) q[3];
sx q[3];
rz(-1.9683016) q[3];
sx q[3];
rz(-2.60633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9390949) q[2];
sx q[2];
rz(-1.9715344) q[2];
sx q[2];
rz(-1.5411752) q[2];
rz(-1.0209068) q[3];
sx q[3];
rz(-1.7424135) q[3];
sx q[3];
rz(2.8405564) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0114667) q[0];
sx q[0];
rz(-0.13700329) q[0];
sx q[0];
rz(-0.89114183) q[0];
rz(2.4961684) q[1];
sx q[1];
rz(-1.3963457) q[1];
sx q[1];
rz(0.83271629) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7695205) q[0];
sx q[0];
rz(-1.4345048) q[0];
sx q[0];
rz(2.6762677) q[0];
x q[1];
rz(1.4786167) q[2];
sx q[2];
rz(-2.5155009) q[2];
sx q[2];
rz(2.7220059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.558185) q[1];
sx q[1];
rz(-1.3761569) q[1];
sx q[1];
rz(0.50048142) q[1];
rz(-pi) q[2];
rz(-0.54033684) q[3];
sx q[3];
rz(-1.2286006) q[3];
sx q[3];
rz(-2.1006753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.58925313) q[2];
sx q[2];
rz(-2.4746042) q[2];
sx q[2];
rz(1.5480631) q[2];
rz(2.7175236) q[3];
sx q[3];
rz(-1.721902) q[3];
sx q[3];
rz(-0.29485318) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9645204) q[0];
sx q[0];
rz(-1.9449214) q[0];
sx q[0];
rz(-2.3181584) q[0];
rz(-1.9442762) q[1];
sx q[1];
rz(-1.6540534) q[1];
sx q[1];
rz(-1.4607325) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4890503) q[0];
sx q[0];
rz(-0.33777896) q[0];
sx q[0];
rz(-1.2342288) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8180188) q[2];
sx q[2];
rz(-2.5342016) q[2];
sx q[2];
rz(-3.0957019) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3590815) q[1];
sx q[1];
rz(-2.4788279) q[1];
sx q[1];
rz(3.0726391) q[1];
rz(-pi) q[2];
rz(0.72328679) q[3];
sx q[3];
rz(-0.46542612) q[3];
sx q[3];
rz(-1.775072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0921649) q[2];
sx q[2];
rz(-0.66114134) q[2];
sx q[2];
rz(-3.0180422) q[2];
rz(0.83856797) q[3];
sx q[3];
rz(-2.0288012) q[3];
sx q[3];
rz(2.6113094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0274444) q[0];
sx q[0];
rz(-2.1338978) q[0];
sx q[0];
rz(1.4134407) q[0];
rz(-2.24276) q[1];
sx q[1];
rz(-2.2347968) q[1];
sx q[1];
rz(-0.59284219) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57164416) q[0];
sx q[0];
rz(-1.0019394) q[0];
sx q[0];
rz(2.4242899) q[0];
x q[1];
rz(0.18211629) q[2];
sx q[2];
rz(-1.6556634) q[2];
sx q[2];
rz(-3.0013468) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6839016) q[1];
sx q[1];
rz(-1.2448893) q[1];
sx q[1];
rz(-1.1363455) q[1];
rz(-pi) q[2];
rz(-2.4355008) q[3];
sx q[3];
rz(-0.81940813) q[3];
sx q[3];
rz(-3.135526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2841407) q[2];
sx q[2];
rz(-2.393674) q[2];
sx q[2];
rz(-0.22135529) q[2];
rz(1.0434693) q[3];
sx q[3];
rz(-0.9404434) q[3];
sx q[3];
rz(-0.38844696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8429883) q[0];
sx q[0];
rz(-0.28268155) q[0];
sx q[0];
rz(-0.84841949) q[0];
rz(0.09659718) q[1];
sx q[1];
rz(-2.0674457) q[1];
sx q[1];
rz(-2.3883147) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83148513) q[0];
sx q[0];
rz(-1.7979171) q[0];
sx q[0];
rz(3.0390396) q[0];
rz(0.82920544) q[2];
sx q[2];
rz(-1.5463788) q[2];
sx q[2];
rz(-1.3361734) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9946361) q[1];
sx q[1];
rz(-2.5038233) q[1];
sx q[1];
rz(0.90177782) q[1];
rz(0.71452272) q[3];
sx q[3];
rz(-1.2163289) q[3];
sx q[3];
rz(-2.9936341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0775371) q[2];
sx q[2];
rz(-0.83751837) q[2];
sx q[2];
rz(2.688664) q[2];
rz(2.5405267) q[3];
sx q[3];
rz(-1.6744924) q[3];
sx q[3];
rz(2.1452904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.8375028) q[0];
sx q[0];
rz(-1.4488198) q[0];
sx q[0];
rz(-2.6813843) q[0];
rz(0.17644633) q[1];
sx q[1];
rz(-0.23527589) q[1];
sx q[1];
rz(-1.489524) q[1];
rz(-2.7929581) q[2];
sx q[2];
rz(-1.5713816) q[2];
sx q[2];
rz(-1.1272507) q[2];
rz(2.413977) q[3];
sx q[3];
rz(-1.1066827) q[3];
sx q[3];
rz(-0.91409693) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6018574) q[0];
sx q[0];
rz(-2.5849113) q[0];
sx q[0];
rz(-0.88168983) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2615112) q[2];
sx q[2];
rz(-2.4158084) q[2];
sx q[2];
rz(-0.45553614) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.18892203) q[1];
sx q[1];
rz(-2.1181207) q[1];
sx q[1];
rz(0.31618936) q[1];
x q[2];
rz(1.5815758) q[3];
sx q[3];
rz(-1.5228776) q[3];
sx q[3];
rz(0.77424327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.1047487) q[2];
sx q[2];
rz(-1.8061183) q[2];
sx q[2];
rz(2.74995) q[2];
rz(-0.13970217) q[3];
sx q[3];
rz(-2.468686) q[3];
sx q[3];
rz(0.02123775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4166819) q[0];
sx q[0];
rz(-1.7371438) q[0];
sx q[0];
rz(-1.2765983) q[0];
rz(0.87031594) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(-1.2044027) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5088168) q[0];
sx q[0];
rz(-0.59558376) q[0];
sx q[0];
rz(2.3217208) q[0];
rz(-1.7877098) q[2];
sx q[2];
rz(-2.4241944) q[2];
sx q[2];
rz(0.72620981) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.054489) q[1];
sx q[1];
rz(-2.7627355) q[1];
sx q[1];
rz(-2.8701251) q[1];
x q[2];
rz(-2.2435917) q[3];
sx q[3];
rz(-1.7963396) q[3];
sx q[3];
rz(-2.0190092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5045972) q[2];
sx q[2];
rz(-0.4578788) q[2];
sx q[2];
rz(1.2191999) q[2];
rz(-2.7870264) q[3];
sx q[3];
rz(-0.48186007) q[3];
sx q[3];
rz(-1.5725117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59042674) q[0];
sx q[0];
rz(-1.6141491) q[0];
sx q[0];
rz(-2.5235126) q[0];
rz(-3.1255426) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(-1.9504257) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5844849) q[0];
sx q[0];
rz(-0.35042052) q[0];
sx q[0];
rz(2.0273897) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3267924) q[2];
sx q[2];
rz(-1.8463496) q[2];
sx q[2];
rz(-0.22305605) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.788113) q[1];
sx q[1];
rz(-0.2971) q[1];
sx q[1];
rz(0.77570813) q[1];
x q[2];
rz(-2.5997945) q[3];
sx q[3];
rz(-1.9199315) q[3];
sx q[3];
rz(2.9475714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.96737343) q[2];
sx q[2];
rz(-0.9625532) q[2];
sx q[2];
rz(-1.013914) q[2];
rz(-1.3570471) q[3];
sx q[3];
rz(-0.8299399) q[3];
sx q[3];
rz(1.0323662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8106666) q[0];
sx q[0];
rz(-1.4241011) q[0];
sx q[0];
rz(0.20430918) q[0];
rz(1.7640242) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(-2.2185982) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6072236) q[0];
sx q[0];
rz(-1.9019706) q[0];
sx q[0];
rz(-1.8934728) q[0];
rz(0.043945233) q[2];
sx q[2];
rz(-1.5590132) q[2];
sx q[2];
rz(2.1859283) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7912485) q[1];
sx q[1];
rz(-1.6748168) q[1];
sx q[1];
rz(1.1715602) q[1];
rz(-pi) q[2];
x q[2];
rz(3.119167) q[3];
sx q[3];
rz(-0.86655819) q[3];
sx q[3];
rz(1.1426329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4541645) q[2];
sx q[2];
rz(-1.585008) q[2];
sx q[2];
rz(-0.50951177) q[2];
rz(-2.4984958) q[3];
sx q[3];
rz(-2.0431079) q[3];
sx q[3];
rz(0.96737868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.61280695) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(0.91941419) q[0];
rz(-2.1557504) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(1.7808419) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0589941) q[0];
sx q[0];
rz(-1.3816125) q[0];
sx q[0];
rz(0.22342213) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1941357) q[2];
sx q[2];
rz(-0.91295487) q[2];
sx q[2];
rz(2.5156977) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1636703) q[1];
sx q[1];
rz(-2.0241892) q[1];
sx q[1];
rz(1.8015253) q[1];
rz(-pi) q[2];
rz(-2.1702483) q[3];
sx q[3];
rz(-1.6347307) q[3];
sx q[3];
rz(1.5435227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3551066) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(3.0267267) q[2];
rz(2.6857175) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(-1.280064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24546656) q[0];
sx q[0];
rz(-1.8969314) q[0];
sx q[0];
rz(1.2448357) q[0];
rz(-2.4339829) q[1];
sx q[1];
rz(-1.6376303) q[1];
sx q[1];
rz(-0.27522603) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21101418) q[0];
sx q[0];
rz(-1.0761346) q[0];
sx q[0];
rz(0.37822322) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92991021) q[2];
sx q[2];
rz(-1.1615796) q[2];
sx q[2];
rz(2.2355459) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5425307) q[1];
sx q[1];
rz(-1.2236569) q[1];
sx q[1];
rz(0.51862006) q[1];
rz(-2.1915057) q[3];
sx q[3];
rz(-1.6380777) q[3];
sx q[3];
rz(-0.34611191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.747921) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(2.6829524) q[2];
rz(0.7115055) q[3];
sx q[3];
rz(-0.68370521) q[3];
sx q[3];
rz(-2.5512364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8608619) q[0];
sx q[0];
rz(-1.9545398) q[0];
sx q[0];
rz(1.3635427) q[0];
rz(1.4200312) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(-2.4218959) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5170383) q[0];
sx q[0];
rz(-1.2861684) q[0];
sx q[0];
rz(2.4863003) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4160956) q[2];
sx q[2];
rz(-2.2273387) q[2];
sx q[2];
rz(3.0692435) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4622725) q[1];
sx q[1];
rz(-1.5264891) q[1];
sx q[1];
rz(0.54507749) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8411438) q[3];
sx q[3];
rz(-2.2336604) q[3];
sx q[3];
rz(0.28966749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.482243) q[2];
sx q[2];
rz(-1.6094001) q[2];
sx q[2];
rz(-2.4534295) q[2];
rz(-3.0996389) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(0.28013128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2280837) q[0];
sx q[0];
rz(-1.1627473) q[0];
sx q[0];
rz(-0.18653175) q[0];
rz(-0.60449156) q[1];
sx q[1];
rz(-2.1291321) q[1];
sx q[1];
rz(1.6465181) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.667244) q[0];
sx q[0];
rz(-1.357204) q[0];
sx q[0];
rz(2.7937907) q[0];
rz(2.9096793) q[2];
sx q[2];
rz(-0.70966087) q[2];
sx q[2];
rz(-2.2556925) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5563158) q[1];
sx q[1];
rz(-1.3116515) q[1];
sx q[1];
rz(1.4633333) q[1];
rz(-pi) q[2];
rz(0.64077326) q[3];
sx q[3];
rz(-2.2230004) q[3];
sx q[3];
rz(-0.81354248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8994393) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(-3.126826) q[2];
rz(1.8642558) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(1.4130672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16167851) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(2.263608) q[0];
rz(0.19628482) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(-1.3605798) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9961062) q[0];
sx q[0];
rz(-1.5331475) q[0];
sx q[0];
rz(-0.99897511) q[0];
rz(-pi) q[1];
rz(-0.94366818) q[2];
sx q[2];
rz(-1.5329648) q[2];
sx q[2];
rz(0.19247069) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.66879241) q[1];
sx q[1];
rz(-1.2392514) q[1];
sx q[1];
rz(-1.214295) q[1];
rz(-0.15501539) q[3];
sx q[3];
rz(-1.1559556) q[3];
sx q[3];
rz(-2.993194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.62375162) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(2.2488135) q[2];
rz(-0.039285224) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(0.75004309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2980625) q[0];
sx q[0];
rz(-0.003304464) q[0];
sx q[0];
rz(0.046534006) q[0];
rz(-0.90905601) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(-0.7199026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4601446) q[0];
sx q[0];
rz(-1.8143192) q[0];
sx q[0];
rz(-0.081865099) q[0];
x q[1];
rz(-0.30883046) q[2];
sx q[2];
rz(-1.575982) q[2];
sx q[2];
rz(-1.3877102) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7457909) q[1];
sx q[1];
rz(-1.5862641) q[1];
sx q[1];
rz(-1.323911) q[1];
x q[2];
rz(2.4107237) q[3];
sx q[3];
rz(-2.1248397) q[3];
sx q[3];
rz(2.4027367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7211192) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(3.1151248) q[2];
rz(1.3039533) q[3];
sx q[3];
rz(-0.81726685) q[3];
sx q[3];
rz(1.8855689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883041) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(-0.2164671) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(2.7770476) q[2];
sx q[2];
rz(-2.400763) q[2];
sx q[2];
rz(2.7999067) q[2];
rz(-0.032565928) q[3];
sx q[3];
rz(-0.98441549) q[3];
sx q[3];
rz(2.050642) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
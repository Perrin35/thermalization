OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(-0.35020819) q[0];
sx q[0];
rz(-0.36663088) q[0];
rz(-2.2740013) q[1];
sx q[1];
rz(-0.64414135) q[1];
sx q[1];
rz(1.6860513) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35858425) q[0];
sx q[0];
rz(-1.2281679) q[0];
sx q[0];
rz(1.1230099) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8777962) q[2];
sx q[2];
rz(-2.2552239) q[2];
sx q[2];
rz(-0.85927187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75016025) q[1];
sx q[1];
rz(-0.62392985) q[1];
sx q[1];
rz(-2.0425914) q[1];
x q[2];
rz(-2.9204908) q[3];
sx q[3];
rz(-0.049115291) q[3];
sx q[3];
rz(-0.55288314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1047487) q[2];
sx q[2];
rz(-1.3354744) q[2];
sx q[2];
rz(2.74995) q[2];
rz(-3.0018905) q[3];
sx q[3];
rz(-2.468686) q[3];
sx q[3];
rz(3.1203549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7249107) q[0];
sx q[0];
rz(-1.4044489) q[0];
sx q[0];
rz(-1.8649944) q[0];
rz(0.87031594) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(-1.2044027) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6635839) q[0];
sx q[0];
rz(-1.9933797) q[0];
sx q[0];
rz(-0.43310662) q[0];
x q[1];
rz(-2.9559829) q[2];
sx q[2];
rz(-0.87366548) q[2];
sx q[2];
rz(-0.44167659) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7960297) q[1];
sx q[1];
rz(-1.2064762) q[1];
sx q[1];
rz(1.6771392) q[1];
x q[2];
rz(-2.2435917) q[3];
sx q[3];
rz(-1.3452531) q[3];
sx q[3];
rz(-1.1225835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63699547) q[2];
sx q[2];
rz(-2.6837139) q[2];
sx q[2];
rz(-1.2191999) q[2];
rz(2.7870264) q[3];
sx q[3];
rz(-0.48186007) q[3];
sx q[3];
rz(-1.569081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59042674) q[0];
sx q[0];
rz(-1.6141491) q[0];
sx q[0];
rz(0.61808008) q[0];
rz(-3.1255426) q[1];
sx q[1];
rz(-2.2647104) q[1];
sx q[1];
rz(-1.191167) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5601657) q[0];
sx q[0];
rz(-1.7227356) q[0];
sx q[0];
rz(1.8877958) q[0];
x q[1];
rz(1.3267924) q[2];
sx q[2];
rz(-1.2952431) q[2];
sx q[2];
rz(2.9185366) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.788113) q[1];
sx q[1];
rz(-0.2971) q[1];
sx q[1];
rz(-0.77570813) q[1];
rz(-1.1690087) q[3];
sx q[3];
rz(-1.0649293) q[3];
sx q[3];
rz(-1.5617621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1742192) q[2];
sx q[2];
rz(-2.1790395) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33092609) q[0];
sx q[0];
rz(-1.4241011) q[0];
sx q[0];
rz(-2.9372835) q[0];
rz(1.3775685) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(-0.92299443) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92813338) q[0];
sx q[0];
rz(-1.8753578) q[0];
sx q[0];
rz(-0.34780985) q[0];
rz(-pi) q[1];
rz(1.5590018) q[2];
sx q[2];
rz(-1.5268541) q[2];
sx q[2];
rz(-2.5269788) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35034414) q[1];
sx q[1];
rz(-1.6748168) q[1];
sx q[1];
rz(-1.1715602) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5971848) q[3];
sx q[3];
rz(-2.4370586) q[3];
sx q[3];
rz(-1.964331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6874281) q[2];
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
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6960984) q[0];
sx q[0];
rz(-1.7901663) q[0];
sx q[0];
rz(1.7646837) q[0];
rz(-pi) q[1];
rz(-2.4945716) q[2];
sx q[2];
rz(-0.87304742) q[2];
sx q[2];
rz(0.24017142) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9779224) q[1];
sx q[1];
rz(-2.0241892) q[1];
sx q[1];
rz(-1.8015253) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0642062) q[3];
sx q[3];
rz(-0.97273982) q[3];
sx q[3];
rz(-0.070904562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78648606) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(-3.0267267) q[2];
rz(-0.45587513) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(-1.280064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961261) q[0];
sx q[0];
rz(-1.8969314) q[0];
sx q[0];
rz(-1.8967569) q[0];
rz(-2.4339829) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(-2.8663666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5462287) q[0];
sx q[0];
rz(-1.2397791) q[0];
sx q[0];
rz(2.0966895) q[0];
rz(-pi) q[1];
rz(-0.49595828) q[2];
sx q[2];
rz(-0.99018103) q[2];
sx q[2];
rz(-0.95326391) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6323159) q[1];
sx q[1];
rz(-2.5264611) q[1];
sx q[1];
rz(-2.5110911) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.082645881) q[3];
sx q[3];
rz(-0.9517037) q[3];
sx q[3];
rz(-1.9649399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.39367166) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(-0.45864027) q[2];
rz(-2.4300872) q[3];
sx q[3];
rz(-2.4578874) q[3];
sx q[3];
rz(2.5512364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28073072) q[0];
sx q[0];
rz(-1.1870528) q[0];
sx q[0];
rz(-1.3635427) q[0];
rz(-1.7215615) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(0.71969676) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4079106) q[0];
sx q[0];
rz(-0.94607293) q[0];
sx q[0];
rz(1.2172933) q[0];
rz(2.2816706) q[2];
sx q[2];
rz(-0.93647525) q[2];
sx q[2];
rz(-2.2459522) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4622725) q[1];
sx q[1];
rz(-1.6151036) q[1];
sx q[1];
rz(0.54507749) q[1];
rz(-pi) q[2];
rz(-2.8119874) q[3];
sx q[3];
rz(-0.70809396) q[3];
sx q[3];
rz(0.13347382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.65934962) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(2.4534295) q[2];
rz(3.0996389) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(2.8614614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.2280837) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(-0.18653175) q[0];
rz(-2.5371011) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(-1.4950745) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.667244) q[0];
sx q[0];
rz(-1.7843887) q[0];
sx q[0];
rz(0.347802) q[0];
rz(-pi) q[1];
rz(1.375884) q[2];
sx q[2];
rz(-0.88390985) q[2];
sx q[2];
rz(-1.1877103) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.013156803) q[1];
sx q[1];
rz(-1.4669347) q[1];
sx q[1];
rz(2.8810112) q[1];
rz(-pi) q[2];
rz(2.2349615) q[3];
sx q[3];
rz(-0.88007054) q[3];
sx q[3];
rz(0.074113473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24215332) q[2];
sx q[2];
rz(-2.184325) q[2];
sx q[2];
rz(0.014766679) q[2];
rz(1.8642558) q[3];
sx q[3];
rz(-1.8321962) q[3];
sx q[3];
rz(1.7285255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16167851) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(2.263608) q[0];
rz(-0.19628482) q[1];
sx q[1];
rz(-1.9613962) q[1];
sx q[1];
rz(-1.3605798) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1454865) q[0];
sx q[0];
rz(-1.6084451) q[0];
sx q[0];
rz(0.99897511) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5063862) q[2];
sx q[2];
rz(-0.62811479) q[2];
sx q[2];
rz(1.8154085) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4728002) q[1];
sx q[1];
rz(-1.9023412) q[1];
sx q[1];
rz(-1.214295) q[1];
rz(1.1514879) q[3];
sx q[3];
rz(-1.4290223) q[3];
sx q[3];
rz(-1.7820953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.62375162) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(0.8927792) q[2];
rz(0.039285224) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(2.3915496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2980625) q[0];
sx q[0];
rz(-3.1382882) q[0];
sx q[0];
rz(0.046534006) q[0];
rz(0.90905601) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(-2.4216901) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68144804) q[0];
sx q[0];
rz(-1.3272734) q[0];
sx q[0];
rz(3.0597276) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1245329) q[2];
sx q[2];
rz(-2.8327201) q[2];
sx q[2];
rz(2.942254) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1788927) q[1];
sx q[1];
rz(-1.3239412) q[1];
sx q[1];
rz(-3.1256413) q[1];
rz(2.2640957) q[3];
sx q[3];
rz(-0.96713669) q[3];
sx q[3];
rz(2.7503777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4204734) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(-3.1151248) q[2];
rz(-1.3039533) q[3];
sx q[3];
rz(-0.81726685) q[3];
sx q[3];
rz(1.2560237) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7883041) q[0];
sx q[0];
rz(-1.371654) q[0];
sx q[0];
rz(-1.3375682) q[0];
rz(-2.9251255) q[1];
sx q[1];
rz(-1.4420061) q[1];
sx q[1];
rz(-1.9056086) q[1];
rz(2.7770476) q[2];
sx q[2];
rz(-2.400763) q[2];
sx q[2];
rz(2.7999067) q[2];
rz(0.032565928) q[3];
sx q[3];
rz(-2.1571772) q[3];
sx q[3];
rz(-1.0909506) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
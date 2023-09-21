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
rz(2.7749618) q[0];
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
rz(-1.9134247) q[0];
sx q[0];
rz(2.0185828) q[0];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9526706) q[1];
sx q[1];
rz(-2.1181207) q[1];
sx q[1];
rz(0.31618936) q[1];
rz(-pi) q[2];
rz(0.22110181) q[3];
sx q[3];
rz(-3.0924774) q[3];
sx q[3];
rz(0.55288314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.036844) q[2];
sx q[2];
rz(-1.3354744) q[2];
sx q[2];
rz(-0.39164266) q[2];
rz(-0.13970217) q[3];
sx q[3];
rz(-0.67290664) q[3];
sx q[3];
rz(3.1203549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4166819) q[0];
sx q[0];
rz(-1.4044489) q[0];
sx q[0];
rz(1.8649944) q[0];
rz(-0.87031594) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(1.2044027) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28019529) q[0];
sx q[0];
rz(-1.9635927) q[0];
sx q[0];
rz(2.0307721) q[0];
rz(-0.86512489) q[2];
sx q[2];
rz(-1.428831) q[2];
sx q[2];
rz(-1.0091458) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.345563) q[1];
sx q[1];
rz(-1.2064762) q[1];
sx q[1];
rz(1.4644535) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28537206) q[3];
sx q[3];
rz(-0.91802363) q[3];
sx q[3];
rz(0.6245581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.63699547) q[2];
sx q[2];
rz(-2.6837139) q[2];
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
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59042674) q[0];
sx q[0];
rz(-1.5274436) q[0];
sx q[0];
rz(2.5235126) q[0];
rz(-3.1255426) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(1.191167) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5601657) q[0];
sx q[0];
rz(-1.418857) q[0];
sx q[0];
rz(1.2537969) q[0];
x q[1];
rz(0.28352719) q[2];
sx q[2];
rz(-1.3361738) q[2];
sx q[2];
rz(-1.7262176) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.788113) q[1];
sx q[1];
rz(-0.2971) q[1];
sx q[1];
rz(0.77570813) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5997945) q[3];
sx q[3];
rz(-1.9199315) q[3];
sx q[3];
rz(2.9475714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96737343) q[2];
sx q[2];
rz(-2.1790395) q[2];
sx q[2];
rz(1.013914) q[2];
rz(1.7845456) q[3];
sx q[3];
rz(-0.8299399) q[3];
sx q[3];
rz(-2.1092265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.4148477) q[1];
sx q[1];
rz(2.2185982) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6072236) q[0];
sx q[0];
rz(-1.2396221) q[0];
sx q[0];
rz(1.2481199) q[0];
x q[1];
rz(-0.043945233) q[2];
sx q[2];
rz(-1.5590132) q[2];
sx q[2];
rz(0.95566434) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.877362) q[1];
sx q[1];
rz(-1.1738395) q[1];
sx q[1];
rz(0.11282679) q[1];
rz(0.8664341) q[3];
sx q[3];
rz(-1.5537062) q[3];
sx q[3];
rz(-2.7279502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6874281) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(2.6320809) q[2];
rz(2.4984958) q[3];
sx q[3];
rz(-1.0984848) q[3];
sx q[3];
rz(0.96737868) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61280695) q[0];
sx q[0];
rz(-1.9354154) q[0];
sx q[0];
rz(-2.2221785) q[0];
rz(2.1557504) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(-1.7808419) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0589941) q[0];
sx q[0];
rz(-1.7599802) q[0];
sx q[0];
rz(0.22342213) q[0];
rz(-pi) q[1];
rz(2.4945716) q[2];
sx q[2];
rz(-2.2685452) q[2];
sx q[2];
rz(0.24017142) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1636703) q[1];
sx q[1];
rz(-1.1174035) q[1];
sx q[1];
rz(-1.8015253) q[1];
x q[2];
rz(-1.4578044) q[3];
sx q[3];
rz(-2.5391573) q[3];
sx q[3];
rz(-0.065954176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3551066) q[2];
sx q[2];
rz(-1.3856709) q[2];
sx q[2];
rz(-3.0267267) q[2];
rz(0.45587513) q[3];
sx q[3];
rz(-1.8084278) q[3];
sx q[3];
rz(-1.280064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24546656) q[0];
sx q[0];
rz(-1.2446612) q[0];
sx q[0];
rz(1.2448357) q[0];
rz(0.70760977) q[1];
sx q[1];
rz(-1.5039624) q[1];
sx q[1];
rz(-2.8663666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9305785) q[0];
sx q[0];
rz(-2.0654581) q[0];
sx q[0];
rz(-2.7633694) q[0];
rz(0.94324865) q[2];
sx q[2];
rz(-2.3970282) q[2];
sx q[2];
rz(-0.17472357) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.59906193) q[1];
sx q[1];
rz(-1.9179357) q[1];
sx q[1];
rz(0.51862006) q[1];
rz(-pi) q[2];
rz(3.0589468) q[3];
sx q[3];
rz(-0.9517037) q[3];
sx q[3];
rz(1.1766528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39367166) q[2];
sx q[2];
rz(-1.5449521) q[2];
sx q[2];
rz(0.45864027) q[2];
rz(0.7115055) q[3];
sx q[3];
rz(-0.68370521) q[3];
sx q[3];
rz(0.59035629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.8608619) q[0];
sx q[0];
rz(-1.9545398) q[0];
sx q[0];
rz(1.3635427) q[0];
rz(-1.4200312) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(2.4218959) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8450711) q[0];
sx q[0];
rz(-2.4356027) q[0];
sx q[0];
rz(-0.4476053) q[0];
x q[1];
rz(-0.77064387) q[2];
sx q[2];
rz(-1.0174123) q[2];
sx q[2];
rz(1.9945952) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.67932018) q[1];
sx q[1];
rz(-1.5264891) q[1];
sx q[1];
rz(2.5965152) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4607055) q[3];
sx q[3];
rz(-1.7828935) q[3];
sx q[3];
rz(-1.6915481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.65934962) q[2];
sx q[2];
rz(-1.6094001) q[2];
sx q[2];
rz(-0.68816319) q[2];
rz(-0.041953772) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(2.8614614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2280837) q[0];
sx q[0];
rz(-1.1627473) q[0];
sx q[0];
rz(2.9550609) q[0];
rz(2.5371011) q[1];
sx q[1];
rz(-2.1291321) q[1];
sx q[1];
rz(-1.4950745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5160617) q[0];
sx q[0];
rz(-2.7357258) q[0];
sx q[0];
rz(-2.5748475) q[0];
rz(1.375884) q[2];
sx q[2];
rz(-0.88390985) q[2];
sx q[2];
rz(1.9538823) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5852768) q[1];
sx q[1];
rz(-1.3116515) q[1];
sx q[1];
rz(-1.4633333) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2349615) q[3];
sx q[3];
rz(-0.88007054) q[3];
sx q[3];
rz(-3.0674792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24215332) q[2];
sx q[2];
rz(-2.184325) q[2];
sx q[2];
rz(-0.014766679) q[2];
rz(-1.2773369) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(1.4130672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16167851) q[0];
sx q[0];
rz(-0.67502397) q[0];
sx q[0];
rz(-0.87798464) q[0];
rz(2.9453078) q[1];
sx q[1];
rz(-1.9613962) q[1];
sx q[1];
rz(1.7810129) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44952794) q[0];
sx q[0];
rz(-0.99943107) q[0];
sx q[0];
rz(3.0968303) q[0];
x q[1];
rz(-1.5063862) q[2];
sx q[2];
rz(-0.62811479) q[2];
sx q[2];
rz(-1.3261842) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.18371836) q[1];
sx q[1];
rz(-0.4819594) q[1];
sx q[1];
rz(-2.3493489) q[1];
rz(2.9865773) q[3];
sx q[3];
rz(-1.985637) q[3];
sx q[3];
rz(-0.14839867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.517841) q[2];
sx q[2];
rz(-1.5072284) q[2];
sx q[2];
rz(-0.8927792) q[2];
rz(-3.1023074) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(2.3915496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2980625) q[0];
sx q[0];
rz(-0.003304464) q[0];
sx q[0];
rz(-0.046534006) q[0];
rz(0.90905601) q[1];
sx q[1];
rz(-0.87288705) q[1];
sx q[1];
rz(-0.7199026) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90912949) q[0];
sx q[0];
rz(-1.4913519) q[0];
sx q[0];
rz(-1.3264873) q[0];
x q[1];
rz(0.30883046) q[2];
sx q[2];
rz(-1.575982) q[2];
sx q[2];
rz(-1.7538824) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.1788927) q[1];
sx q[1];
rz(-1.8176515) q[1];
sx q[1];
rz(-3.1256413) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87749691) q[3];
sx q[3];
rz(-0.96713669) q[3];
sx q[3];
rz(-0.391215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7211192) q[2];
sx q[2];
rz(-1.7420008) q[2];
sx q[2];
rz(-3.1151248) q[2];
rz(1.8376393) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(-1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3532886) q[0];
sx q[0];
rz(-1.7699387) q[0];
sx q[0];
rz(1.8040245) q[0];
rz(-2.9251255) q[1];
sx q[1];
rz(-1.4420061) q[1];
sx q[1];
rz(-1.9056086) q[1];
rz(2.4344865) q[2];
sx q[2];
rz(-1.8137992) q[2];
sx q[2];
rz(0.9546311) q[2];
rz(0.98417102) q[3];
sx q[3];
rz(-1.5979206) q[3];
sx q[3];
rz(0.4978705) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

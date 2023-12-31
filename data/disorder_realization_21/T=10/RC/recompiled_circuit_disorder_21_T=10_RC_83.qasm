OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0857467) q[0];
sx q[0];
rz(-0.081781713) q[0];
sx q[0];
rz(-2.6401289) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(0.3224386) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8034536) q[0];
sx q[0];
rz(-2.8087466) q[0];
sx q[0];
rz(1.3071609) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5311004) q[2];
sx q[2];
rz(-0.98449003) q[2];
sx q[2];
rz(2.5551978) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0281801) q[1];
sx q[1];
rz(-2.1537188) q[1];
sx q[1];
rz(-2.7229573) q[1];
rz(-pi) q[2];
rz(-0.93482165) q[3];
sx q[3];
rz(-1.1879731) q[3];
sx q[3];
rz(-0.32360199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6364608) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(2.5855529) q[2];
rz(-0.83267823) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(0.94579831) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44822025) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(2.9843176) q[0];
rz(2.8804624) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(-3.0325586) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63576525) q[0];
sx q[0];
rz(-0.286239) q[0];
sx q[0];
rz(0.92606996) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5017545) q[2];
sx q[2];
rz(-0.32391641) q[2];
sx q[2];
rz(-1.4878291) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8728767) q[1];
sx q[1];
rz(-2.0375588) q[1];
sx q[1];
rz(2.4750535) q[1];
rz(0.29626366) q[3];
sx q[3];
rz(-0.54013541) q[3];
sx q[3];
rz(-1.4512645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9033501) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(-1.2634574) q[2];
rz(0.3271099) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(1.9272778) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064421244) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(1.3431312) q[0];
rz(-0.24761565) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(2.6599191) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7307229) q[0];
sx q[0];
rz(-1.0856837) q[0];
sx q[0];
rz(-2.7184125) q[0];
rz(-1.1142251) q[2];
sx q[2];
rz(-2.4764428) q[2];
sx q[2];
rz(0.84012023) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8183221) q[1];
sx q[1];
rz(-2.1774877) q[1];
sx q[1];
rz(-2.2491749) q[1];
x q[2];
rz(2.2948202) q[3];
sx q[3];
rz(-2.0623042) q[3];
sx q[3];
rz(-1.5566952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8032916) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(0.49989191) q[2];
rz(2.5806184) q[3];
sx q[3];
rz(-1.2596954) q[3];
sx q[3];
rz(-1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9445779) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(2.5894077) q[0];
rz(1.588297) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(1.8968556) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8053999) q[0];
sx q[0];
rz(-1.5724143) q[0];
sx q[0];
rz(2.8021115) q[0];
x q[1];
rz(1.7719901) q[2];
sx q[2];
rz(-1.6987213) q[2];
sx q[2];
rz(-2.163137) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4501805) q[1];
sx q[1];
rz(-2.1069063) q[1];
sx q[1];
rz(0.24713534) q[1];
x q[2];
rz(-2.7987715) q[3];
sx q[3];
rz(-0.59844136) q[3];
sx q[3];
rz(3.0330021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84919471) q[2];
sx q[2];
rz(-1.8820102) q[2];
sx q[2];
rz(1.9909031) q[2];
rz(1.6644647) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(2.6586444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078159049) q[0];
sx q[0];
rz(-2.3796191) q[0];
sx q[0];
rz(-0.081469014) q[0];
rz(-3.07913) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(1.5030456) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8836356) q[0];
sx q[0];
rz(-0.99146087) q[0];
sx q[0];
rz(2.9888319) q[0];
x q[1];
rz(0.91707768) q[2];
sx q[2];
rz(-0.13609016) q[2];
sx q[2];
rz(1.2166785) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.45080966) q[1];
sx q[1];
rz(-2.8125617) q[1];
sx q[1];
rz(-1.2500709) q[1];
rz(-pi) q[2];
rz(0.50932192) q[3];
sx q[3];
rz(-1.8511726) q[3];
sx q[3];
rz(-2.544739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6333255) q[2];
sx q[2];
rz(-0.94255629) q[2];
sx q[2];
rz(-1.903669) q[2];
rz(1.1226908) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(-2.6200263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6102819) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(0.26671985) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(-0.7985324) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0621588) q[0];
sx q[0];
rz(-1.845813) q[0];
sx q[0];
rz(-0.14607231) q[0];
rz(-2.6449634) q[2];
sx q[2];
rz(-1.2611258) q[2];
sx q[2];
rz(1.3336099) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9763899) q[1];
sx q[1];
rz(-2.5678647) q[1];
sx q[1];
rz(1.4742673) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3409307) q[3];
sx q[3];
rz(-1.1122397) q[3];
sx q[3];
rz(0.82160219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9810527) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(0.36671656) q[2];
rz(1.8803053) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(-0.096207531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40903184) q[0];
sx q[0];
rz(-2.2213187) q[0];
sx q[0];
rz(-2.5352056) q[0];
rz(-2.9442893) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(-0.46404776) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2460829) q[0];
sx q[0];
rz(-1.9946788) q[0];
sx q[0];
rz(-0.8830107) q[0];
rz(-1.4887965) q[2];
sx q[2];
rz(-2.900015) q[2];
sx q[2];
rz(1.8919924) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0276427) q[1];
sx q[1];
rz(-1.0877891) q[1];
sx q[1];
rz(0.34160683) q[1];
rz(-pi) q[2];
rz(1.6528321) q[3];
sx q[3];
rz(-1.8063407) q[3];
sx q[3];
rz(1.9627278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2074034) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(-0.25804538) q[2];
rz(-1.9559654) q[3];
sx q[3];
rz(-1.6262755) q[3];
sx q[3];
rz(-3.1380838) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19514062) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(0.38129693) q[0];
rz(-3.0463468) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(-1.4415178) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9863319) q[0];
sx q[0];
rz(-0.065247029) q[0];
sx q[0];
rz(-0.23373993) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61261119) q[2];
sx q[2];
rz(-1.5393886) q[2];
sx q[2];
rz(0.14234662) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0563593) q[1];
sx q[1];
rz(-1.4409522) q[1];
sx q[1];
rz(-1.4443881) q[1];
rz(-1.5278682) q[3];
sx q[3];
rz(-1.9044442) q[3];
sx q[3];
rz(1.4294525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9986481) q[2];
sx q[2];
rz(-2.7286077) q[2];
sx q[2];
rz(0.22658919) q[2];
rz(-2.6930124) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(-2.3118238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.7609693) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(1.3990078) q[0];
rz(-2.8245068) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(0.98659602) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9771377) q[0];
sx q[0];
rz(-1.4129606) q[0];
sx q[0];
rz(-0.49789238) q[0];
rz(-pi) q[1];
rz(0.95790205) q[2];
sx q[2];
rz(-0.79682486) q[2];
sx q[2];
rz(-0.56607841) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4793195) q[1];
sx q[1];
rz(-2.5179177) q[1];
sx q[1];
rz(-2.8981266) q[1];
rz(-pi) q[2];
rz(-1.1986198) q[3];
sx q[3];
rz(-2.6937006) q[3];
sx q[3];
rz(-0.4263634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9265147) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(-0.40965664) q[2];
rz(0.26327291) q[3];
sx q[3];
rz(-1.8299088) q[3];
sx q[3];
rz(0.51945654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7423994) q[0];
sx q[0];
rz(-0.078646794) q[0];
sx q[0];
rz(1.7364527) q[0];
rz(-2.3204904) q[1];
sx q[1];
rz(-0.91870538) q[1];
sx q[1];
rz(1.7260889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6575359) q[0];
sx q[0];
rz(-2.7298379) q[0];
sx q[0];
rz(-0.62396892) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2070933) q[2];
sx q[2];
rz(-1.3377681) q[2];
sx q[2];
rz(-1.9050913) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0714598) q[1];
sx q[1];
rz(-1.3576641) q[1];
sx q[1];
rz(-0.15503426) q[1];
rz(-pi) q[2];
rz(-1.8652925) q[3];
sx q[3];
rz(-2.6101972) q[3];
sx q[3];
rz(2.6607799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.71904174) q[2];
sx q[2];
rz(-0.30095235) q[2];
sx q[2];
rz(-0.12410513) q[2];
rz(-2.1758046) q[3];
sx q[3];
rz(-1.6671168) q[3];
sx q[3];
rz(0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60349764) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(0.30766906) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(-1.3901426) q[2];
sx q[2];
rz(-1.8271108) q[2];
sx q[2];
rz(1.9432632) q[2];
rz(-0.55660558) q[3];
sx q[3];
rz(-1.442933) q[3];
sx q[3];
rz(1.8419151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

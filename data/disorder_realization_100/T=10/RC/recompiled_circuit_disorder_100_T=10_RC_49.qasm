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
rz(0.86759138) q[1];
sx q[1];
rz(-2.4974513) q[1];
sx q[1];
rz(-1.6860513) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35858425) q[0];
sx q[0];
rz(-1.2281679) q[0];
sx q[0];
rz(-1.1230099) q[0];
rz(-pi) q[1];
rz(0.86906616) q[2];
sx q[2];
rz(-1.7742187) q[2];
sx q[2];
rz(-2.2609401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9283596) q[1];
sx q[1];
rz(-1.302049) q[1];
sx q[1];
rz(1.0006204) q[1];
rz(-pi) q[2];
rz(-2.9204908) q[3];
sx q[3];
rz(-0.049115291) q[3];
sx q[3];
rz(2.5887095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.4166819) q[0];
sx q[0];
rz(-1.4044489) q[0];
sx q[0];
rz(1.2765983) q[0];
rz(-0.87031594) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(-1.93719) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6635839) q[0];
sx q[0];
rz(-1.9933797) q[0];
sx q[0];
rz(-2.708486) q[0];
rz(2.9559829) q[2];
sx q[2];
rz(-2.2679272) q[2];
sx q[2];
rz(-0.44167659) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26324998) q[1];
sx q[1];
rz(-1.471457) q[1];
sx q[1];
rz(-0.36621014) q[1];
rz(-pi) q[2];
rz(2.2435917) q[3];
sx q[3];
rz(-1.7963396) q[3];
sx q[3];
rz(2.0190092) q[3];
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
rz(2.7870264) q[3];
sx q[3];
rz(-2.6597326) q[3];
sx q[3];
rz(-1.5725117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5511659) q[0];
sx q[0];
rz(-1.6141491) q[0];
sx q[0];
rz(0.61808008) q[0];
rz(-3.1255426) q[1];
sx q[1];
rz(-2.2647104) q[1];
sx q[1];
rz(1.9504257) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55710775) q[0];
sx q[0];
rz(-2.7911721) q[0];
sx q[0];
rz(2.0273897) q[0];
rz(1.3267924) q[2];
sx q[2];
rz(-1.8463496) q[2];
sx q[2];
rz(-2.9185366) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9900073) q[1];
sx q[1];
rz(-1.3602435) q[1];
sx q[1];
rz(1.7819808) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54179811) q[3];
sx q[3];
rz(-1.2216611) q[3];
sx q[3];
rz(-0.19402129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1742192) q[2];
sx q[2];
rz(-2.1790395) q[2];
sx q[2];
rz(-1.013914) q[2];
rz(-1.3570471) q[3];
sx q[3];
rz(-2.3116528) q[3];
sx q[3];
rz(2.1092265) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8106666) q[0];
sx q[0];
rz(-1.4241011) q[0];
sx q[0];
rz(-0.20430918) q[0];
rz(1.7640242) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(-2.2185982) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8078634) q[0];
sx q[0];
rz(-0.45818746) q[0];
sx q[0];
rz(0.74497594) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.043945233) q[2];
sx q[2];
rz(-1.5590132) q[2];
sx q[2];
rz(-2.1859283) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.35034414) q[1];
sx q[1];
rz(-1.4667759) q[1];
sx q[1];
rz(1.9700325) q[1];
rz(-0.022425671) q[3];
sx q[3];
rz(-2.2750345) q[3];
sx q[3];
rz(1.9989597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6874281) q[2];
sx q[2];
rz(-1.585008) q[2];
sx q[2];
rz(-2.6320809) q[2];
rz(-0.64309684) q[3];
sx q[3];
rz(-1.0984848) q[3];
sx q[3];
rz(-2.174214) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4454942) q[0];
sx q[0];
rz(-1.3514263) q[0];
sx q[0];
rz(-1.3769089) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64702101) q[2];
sx q[2];
rz(-2.2685452) q[2];
sx q[2];
rz(-2.9014212) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1636703) q[1];
sx q[1];
rz(-1.1174035) q[1];
sx q[1];
rz(-1.3400673) q[1];
rz(-pi) q[2];
rz(-1.6837882) q[3];
sx q[3];
rz(-0.60243536) q[3];
sx q[3];
rz(3.0756385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.78648606) q[2];
sx q[2];
rz(-1.7559218) q[2];
sx q[2];
rz(0.11486593) q[2];
rz(-2.6857175) q[3];
sx q[3];
rz(-1.8084278) q[3];
sx q[3];
rz(-1.280064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(-1.6376303) q[1];
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
x q[1];
rz(2.198344) q[2];
sx q[2];
rz(-0.74456442) q[2];
sx q[2];
rz(-0.17472357) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5425307) q[1];
sx q[1];
rz(-1.2236569) q[1];
sx q[1];
rz(-0.51862006) q[1];
rz(-1.6861378) q[3];
sx q[3];
rz(-2.5177258) q[3];
sx q[3];
rz(-1.8231525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.747921) q[2];
sx q[2];
rz(-1.5966406) q[2];
sx q[2];
rz(2.6829524) q[2];
rz(0.7115055) q[3];
sx q[3];
rz(-2.4578874) q[3];
sx q[3];
rz(2.5512364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.8608619) q[0];
sx q[0];
rz(-1.9545398) q[0];
sx q[0];
rz(-1.77805) q[0];
rz(1.4200312) q[1];
sx q[1];
rz(-0.51912156) q[1];
sx q[1];
rz(-2.4218959) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73368209) q[0];
sx q[0];
rz(-2.1955197) q[0];
sx q[0];
rz(-1.9242994) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2816706) q[2];
sx q[2];
rz(-0.93647525) q[2];
sx q[2];
rz(0.89564043) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.67932018) q[1];
sx q[1];
rz(-1.5264891) q[1];
sx q[1];
rz(-2.5965152) q[1];
rz(-1.8411438) q[3];
sx q[3];
rz(-2.2336604) q[3];
sx q[3];
rz(0.28966749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.65934962) q[2];
sx q[2];
rz(-1.5321926) q[2];
sx q[2];
rz(-2.4534295) q[2];
rz(0.041953772) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(0.28013128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91350895) q[0];
sx q[0];
rz(-1.1627473) q[0];
sx q[0];
rz(-2.9550609) q[0];
rz(0.60449156) q[1];
sx q[1];
rz(-1.0124606) q[1];
sx q[1];
rz(1.6465181) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5160617) q[0];
sx q[0];
rz(-2.7357258) q[0];
sx q[0];
rz(-0.56674515) q[0];
x q[1];
rz(0.23191339) q[2];
sx q[2];
rz(-0.70966087) q[2];
sx q[2];
rz(2.2556925) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1868134) q[1];
sx q[1];
rz(-2.8615132) q[1];
sx q[1];
rz(2.7571452) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5008194) q[3];
sx q[3];
rz(-2.2230004) q[3];
sx q[3];
rz(-0.81354248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8994393) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(0.014766679) q[2];
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
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9799141) q[0];
sx q[0];
rz(-2.4665687) q[0];
sx q[0];
rz(-0.87798464) q[0];
rz(2.9453078) q[1];
sx q[1];
rz(-1.1801964) q[1];
sx q[1];
rz(1.3605798) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44952794) q[0];
sx q[0];
rz(-0.99943107) q[0];
sx q[0];
rz(-3.0968303) q[0];
rz(-pi) q[1];
rz(2.1979245) q[2];
sx q[2];
rz(-1.5329648) q[2];
sx q[2];
rz(0.19247069) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9578743) q[1];
sx q[1];
rz(-0.4819594) q[1];
sx q[1];
rz(0.79224371) q[1];
x q[2];
rz(-0.15501539) q[3];
sx q[3];
rz(-1.1559556) q[3];
sx q[3];
rz(0.14839867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.62375162) q[2];
sx q[2];
rz(-1.6343642) q[2];
sx q[2];
rz(-2.2488135) q[2];
rz(-0.039285224) q[3];
sx q[3];
rz(-1.6623442) q[3];
sx q[3];
rz(-2.3915496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84353012) q[0];
sx q[0];
rz(-0.003304464) q[0];
sx q[0];
rz(3.0950586) q[0];
rz(0.90905601) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(0.7199026) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7881308) q[0];
sx q[0];
rz(-0.25665584) q[0];
sx q[0];
rz(-1.8887595) q[0];
rz(-pi) q[1];
rz(-0.017059762) q[2];
sx q[2];
rz(-0.3088726) q[2];
sx q[2];
rz(-0.19933867) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.1137005) q[1];
sx q[1];
rz(-0.24735951) q[1];
sx q[1];
rz(1.6340096) q[1];
rz(-0.7474483) q[3];
sx q[3];
rz(-2.2564853) q[3];
sx q[3];
rz(-1.7789343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7211192) q[2];
sx q[2];
rz(-1.3995918) q[2];
sx q[2];
rz(-0.026467888) q[2];
rz(-1.8376393) q[3];
sx q[3];
rz(-2.3243258) q[3];
sx q[3];
rz(1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(0.2164671) q[1];
sx q[1];
rz(-1.4420061) q[1];
sx q[1];
rz(-1.9056086) q[1];
rz(0.70710612) q[2];
sx q[2];
rz(-1.3277935) q[2];
sx q[2];
rz(-2.1869616) q[2];
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

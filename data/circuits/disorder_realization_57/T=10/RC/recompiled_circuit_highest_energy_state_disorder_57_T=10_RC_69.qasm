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
rz(1.1433831) q[0];
sx q[0];
rz(1.5513865) q[0];
sx q[0];
rz(9.0976465) q[0];
rz(1.7786572) q[1];
sx q[1];
rz(6.8805334) q[1];
sx q[1];
rz(9.7272275) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10714794) q[0];
sx q[0];
rz(-1.4054789) q[0];
sx q[0];
rz(-1.3470696) q[0];
x q[1];
rz(0.91384943) q[2];
sx q[2];
rz(-0.7759594) q[2];
sx q[2];
rz(1.0010819) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0781344) q[1];
sx q[1];
rz(-1.7094104) q[1];
sx q[1];
rz(1.8891508) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6685295) q[3];
sx q[3];
rz(-1.6636208) q[3];
sx q[3];
rz(2.0833833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7678307) q[2];
sx q[2];
rz(-3.0475898) q[2];
sx q[2];
rz(0.2730228) q[2];
rz(-2.7772969) q[3];
sx q[3];
rz(-1.4231057) q[3];
sx q[3];
rz(-3.1209768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.6881707) q[0];
sx q[0];
rz(-2.0437129) q[0];
sx q[0];
rz(-1.4349487) q[0];
rz(-0.20225987) q[1];
sx q[1];
rz(-2.5113998) q[1];
sx q[1];
rz(2.8369301) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1218095) q[0];
sx q[0];
rz(-0.87723786) q[0];
sx q[0];
rz(0.32525639) q[0];
x q[1];
rz(-1.1528913) q[2];
sx q[2];
rz(-1.8455681) q[2];
sx q[2];
rz(2.1407549) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2167336) q[1];
sx q[1];
rz(-1.9113411) q[1];
sx q[1];
rz(1.3751283) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6913139) q[3];
sx q[3];
rz(-1.6778062) q[3];
sx q[3];
rz(-1.4329239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8649586) q[2];
sx q[2];
rz(-0.43143299) q[2];
sx q[2];
rz(1.5383833) q[2];
rz(-0.81613427) q[3];
sx q[3];
rz(-2.3147801) q[3];
sx q[3];
rz(2.2468467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64033878) q[0];
sx q[0];
rz(-0.29359874) q[0];
sx q[0];
rz(1.5119934) q[0];
rz(1.8005796) q[1];
sx q[1];
rz(-2.2684542) q[1];
sx q[1];
rz(0.27892932) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4607083) q[0];
sx q[0];
rz(-1.4620004) q[0];
sx q[0];
rz(0.1291424) q[0];
rz(0.90143335) q[2];
sx q[2];
rz(-2.3877904) q[2];
sx q[2];
rz(1.9447226) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3027569) q[1];
sx q[1];
rz(-2.347337) q[1];
sx q[1];
rz(1.1357186) q[1];
rz(-pi) q[2];
rz(-3.1380021) q[3];
sx q[3];
rz(-0.85455214) q[3];
sx q[3];
rz(-2.8479171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.638008) q[2];
sx q[2];
rz(-1.6562485) q[2];
sx q[2];
rz(-0.33835641) q[2];
rz(-1.7342957) q[3];
sx q[3];
rz(-2.0140078) q[3];
sx q[3];
rz(0.97051632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.19891837) q[0];
sx q[0];
rz(-3.119454) q[0];
sx q[0];
rz(2.4778147) q[0];
rz(-2.5860419) q[1];
sx q[1];
rz(-1.6352362) q[1];
sx q[1];
rz(-1.3551855) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4595918) q[0];
sx q[0];
rz(-0.92456901) q[0];
sx q[0];
rz(1.610397) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.023472114) q[2];
sx q[2];
rz(-1.9292574) q[2];
sx q[2];
rz(-1.9752928) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7807676) q[1];
sx q[1];
rz(-1.0064565) q[1];
sx q[1];
rz(2.6884262) q[1];
rz(-pi) q[2];
rz(3.074792) q[3];
sx q[3];
rz(-1.1054174) q[3];
sx q[3];
rz(-2.3540879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3163471) q[2];
sx q[2];
rz(-0.97858846) q[2];
sx q[2];
rz(0.52581954) q[2];
rz(2.53287) q[3];
sx q[3];
rz(-2.5059097) q[3];
sx q[3];
rz(-0.70613247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71259251) q[0];
sx q[0];
rz(-1.9569995) q[0];
sx q[0];
rz(2.1767966) q[0];
rz(-2.7305799) q[1];
sx q[1];
rz(-2.1844468) q[1];
sx q[1];
rz(2.029665) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72090688) q[0];
sx q[0];
rz(-1.5941275) q[0];
sx q[0];
rz(0.015714808) q[0];
rz(-pi) q[1];
x q[1];
rz(1.13222) q[2];
sx q[2];
rz(-0.91735943) q[2];
sx q[2];
rz(2.6613622) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3810242) q[1];
sx q[1];
rz(-2.3443017) q[1];
sx q[1];
rz(-2.3098195) q[1];
rz(0.49647496) q[3];
sx q[3];
rz(-2.8395445) q[3];
sx q[3];
rz(-2.1453017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8370221) q[2];
sx q[2];
rz(-1.3304173) q[2];
sx q[2];
rz(-2.0988317) q[2];
rz(-1.1834831) q[3];
sx q[3];
rz(-0.8684929) q[3];
sx q[3];
rz(-0.97064251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9633164) q[0];
sx q[0];
rz(-1.9177508) q[0];
sx q[0];
rz(-0.29603145) q[0];
rz(3.1405247) q[1];
sx q[1];
rz(-1.8031392) q[1];
sx q[1];
rz(-1.0329049) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84802062) q[0];
sx q[0];
rz(-1.5762934) q[0];
sx q[0];
rz(-1.574541) q[0];
x q[1];
rz(-1.270711) q[2];
sx q[2];
rz(-0.3870766) q[2];
sx q[2];
rz(-2.84969) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.34570899) q[1];
sx q[1];
rz(-2.8879037) q[1];
sx q[1];
rz(-2.9182052) q[1];
x q[2];
rz(1.9719129) q[3];
sx q[3];
rz(-2.2925431) q[3];
sx q[3];
rz(-2.7665319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7615776) q[2];
sx q[2];
rz(-1.4662611) q[2];
sx q[2];
rz(1.638691) q[2];
rz(-0.6014398) q[3];
sx q[3];
rz(-0.99965874) q[3];
sx q[3];
rz(-2.7433266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0346506) q[0];
sx q[0];
rz(-1.5926188) q[0];
sx q[0];
rz(-0.73367992) q[0];
rz(-0.94379464) q[1];
sx q[1];
rz(-1.5422041) q[1];
sx q[1];
rz(-2.5234047) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0252903) q[0];
sx q[0];
rz(-2.9553614) q[0];
sx q[0];
rz(-1.2982606) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84039495) q[2];
sx q[2];
rz(-1.5638509) q[2];
sx q[2];
rz(2.9027129) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.427278) q[1];
sx q[1];
rz(-2.4763417) q[1];
sx q[1];
rz(-1.5284431) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9312956) q[3];
sx q[3];
rz(-2.1749718) q[3];
sx q[3];
rz(-1.9856356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1980629) q[2];
sx q[2];
rz(-1.1897831) q[2];
sx q[2];
rz(-2.5243536) q[2];
rz(1.1715568) q[3];
sx q[3];
rz(-2.7660683) q[3];
sx q[3];
rz(2.530781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1307369) q[0];
sx q[0];
rz(-0.78748381) q[0];
sx q[0];
rz(-2.6686344) q[0];
rz(1.6821776) q[1];
sx q[1];
rz(-1.727203) q[1];
sx q[1];
rz(-3.1088631) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9851042) q[0];
sx q[0];
rz(-1.5384983) q[0];
sx q[0];
rz(1.558139) q[0];
rz(-pi) q[1];
rz(2.031139) q[2];
sx q[2];
rz(-2.2810069) q[2];
sx q[2];
rz(1.4905765) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1188726) q[1];
sx q[1];
rz(-2.0141083) q[1];
sx q[1];
rz(2.9037649) q[1];
x q[2];
rz(1.8440014) q[3];
sx q[3];
rz(-3.0298067) q[3];
sx q[3];
rz(2.8120086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.98256436) q[2];
sx q[2];
rz(-0.50243598) q[2];
sx q[2];
rz(-2.0459797) q[2];
rz(2.380373) q[3];
sx q[3];
rz(-1.2844362) q[3];
sx q[3];
rz(2.718954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465076) q[0];
sx q[0];
rz(-3.0618771) q[0];
sx q[0];
rz(-2.4327143) q[0];
rz(-0.31157663) q[1];
sx q[1];
rz(-2.0496924) q[1];
sx q[1];
rz(0.30837217) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4080183) q[0];
sx q[0];
rz(-0.43034205) q[0];
sx q[0];
rz(0.4056613) q[0];
rz(-1.9498655) q[2];
sx q[2];
rz(-1.1806025) q[2];
sx q[2];
rz(-1.6782111) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4235792) q[1];
sx q[1];
rz(-2.1732959) q[1];
sx q[1];
rz(1.275255) q[1];
rz(1.8294556) q[3];
sx q[3];
rz(-0.80052278) q[3];
sx q[3];
rz(-2.0009934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0333905) q[2];
sx q[2];
rz(-1.6609001) q[2];
sx q[2];
rz(0.22076503) q[2];
rz(-2.0681785) q[3];
sx q[3];
rz(-0.99075166) q[3];
sx q[3];
rz(-2.3555135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49755001) q[0];
sx q[0];
rz(-0.86174091) q[0];
sx q[0];
rz(2.1047237) q[0];
rz(1.4700302) q[1];
sx q[1];
rz(-1.0187047) q[1];
sx q[1];
rz(1.1933614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84813335) q[0];
sx q[0];
rz(-0.443123) q[0];
sx q[0];
rz(-0.070290914) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6380874) q[2];
sx q[2];
rz(-0.76447884) q[2];
sx q[2];
rz(1.6703005) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2630849) q[1];
sx q[1];
rz(-2.2495338) q[1];
sx q[1];
rz(-1.3287067) q[1];
rz(2.0641238) q[3];
sx q[3];
rz(-0.26255373) q[3];
sx q[3];
rz(-1.5880104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6679473) q[2];
sx q[2];
rz(-1.8569943) q[2];
sx q[2];
rz(0.13772193) q[2];
rz(1.600945) q[3];
sx q[3];
rz(-2.9100304) q[3];
sx q[3];
rz(0.93129492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66844983) q[0];
sx q[0];
rz(-1.6669597) q[0];
sx q[0];
rz(-1.054635) q[0];
rz(-2.6558381) q[1];
sx q[1];
rz(-1.9172485) q[1];
sx q[1];
rz(1.6419372) q[1];
rz(-0.68436868) q[2];
sx q[2];
rz(-1.3791313) q[2];
sx q[2];
rz(-2.287938) q[2];
rz(0.53268643) q[3];
sx q[3];
rz(-0.83704994) q[3];
sx q[3];
rz(-0.34216135) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

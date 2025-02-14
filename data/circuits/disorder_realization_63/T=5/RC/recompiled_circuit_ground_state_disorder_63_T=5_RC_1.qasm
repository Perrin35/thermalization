OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37501332) q[0];
sx q[0];
rz(-1.902268) q[0];
sx q[0];
rz(-1.0898606) q[0];
rz(-2.0774948) q[1];
sx q[1];
rz(-1.2531333) q[1];
sx q[1];
rz(-0.90322948) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82683691) q[0];
sx q[0];
rz(-0.96393782) q[0];
sx q[0];
rz(2.3640102) q[0];
rz(-pi) q[1];
rz(-1.5888693) q[2];
sx q[2];
rz(-0.68328349) q[2];
sx q[2];
rz(-2.9073213) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1841539) q[1];
sx q[1];
rz(-1.9136962) q[1];
sx q[1];
rz(-2.5162906) q[1];
x q[2];
rz(1.7890268) q[3];
sx q[3];
rz(-2.6317843) q[3];
sx q[3];
rz(0.17579432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.696306) q[2];
sx q[2];
rz(-0.13662766) q[2];
sx q[2];
rz(2.7083) q[2];
rz(0.28996921) q[3];
sx q[3];
rz(-0.86499298) q[3];
sx q[3];
rz(-0.32052952) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1746154) q[0];
sx q[0];
rz(-2.2791635) q[0];
sx q[0];
rz(1.2906661) q[0];
rz(-2.4621452) q[1];
sx q[1];
rz(-2.3652855) q[1];
sx q[1];
rz(-2.2490833) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3474059) q[0];
sx q[0];
rz(-1.5958324) q[0];
sx q[0];
rz(1.3935318) q[0];
rz(-pi) q[1];
rz(-2.638916) q[2];
sx q[2];
rz(-2.423175) q[2];
sx q[2];
rz(-0.65332149) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.034169056) q[1];
sx q[1];
rz(-0.71798199) q[1];
sx q[1];
rz(-2.7173244) q[1];
rz(-0.99231798) q[3];
sx q[3];
rz(-2.7035993) q[3];
sx q[3];
rz(2.7051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.14900011) q[2];
sx q[2];
rz(-2.0626103) q[2];
sx q[2];
rz(-1.1326257) q[2];
rz(2.4752786) q[3];
sx q[3];
rz(-1.3601114) q[3];
sx q[3];
rz(1.9168436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3742974) q[0];
sx q[0];
rz(-0.39300028) q[0];
sx q[0];
rz(2.1597916) q[0];
rz(-2.1994195) q[1];
sx q[1];
rz(-1.3092382) q[1];
sx q[1];
rz(-0.91032496) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86493409) q[0];
sx q[0];
rz(-0.66682839) q[0];
sx q[0];
rz(-0.72526057) q[0];
rz(1.6907755) q[2];
sx q[2];
rz(-1.2960805) q[2];
sx q[2];
rz(-0.071823013) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7187319) q[1];
sx q[1];
rz(-1.2694468) q[1];
sx q[1];
rz(1.2217997) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1317838) q[3];
sx q[3];
rz(-1.3512423) q[3];
sx q[3];
rz(-1.1846402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2008449) q[2];
sx q[2];
rz(-2.7767599) q[2];
sx q[2];
rz(-2.9546837) q[2];
rz(2.9777891) q[3];
sx q[3];
rz(-2.2713594) q[3];
sx q[3];
rz(-3.0189309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1214889) q[0];
sx q[0];
rz(-1.7262456) q[0];
sx q[0];
rz(0.69766587) q[0];
rz(-2.5949219) q[1];
sx q[1];
rz(-0.29269871) q[1];
sx q[1];
rz(1.1950511) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23604017) q[0];
sx q[0];
rz(-0.70193203) q[0];
sx q[0];
rz(2.1208956) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5279453) q[2];
sx q[2];
rz(-2.9093008) q[2];
sx q[2];
rz(2.7181546) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66142144) q[1];
sx q[1];
rz(-1.464559) q[1];
sx q[1];
rz(1.7376971) q[1];
rz(-pi) q[2];
rz(-1.8175587) q[3];
sx q[3];
rz(-2.8487848) q[3];
sx q[3];
rz(-2.4532401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4635072) q[2];
sx q[2];
rz(-1.7184075) q[2];
sx q[2];
rz(2.9441693) q[2];
rz(-0.10489634) q[3];
sx q[3];
rz(-1.1095122) q[3];
sx q[3];
rz(0.088002861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5941493) q[0];
sx q[0];
rz(-2.6056885) q[0];
sx q[0];
rz(-2.1323668) q[0];
rz(-3.1127473) q[1];
sx q[1];
rz(-1.4776769) q[1];
sx q[1];
rz(0.65863329) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.490059) q[0];
sx q[0];
rz(-2.5680827) q[0];
sx q[0];
rz(0.83322462) q[0];
rz(-pi) q[1];
rz(-2.057107) q[2];
sx q[2];
rz(-1.7520112) q[2];
sx q[2];
rz(-1.043037) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7004905) q[1];
sx q[1];
rz(-1.6630739) q[1];
sx q[1];
rz(-2.9725595) q[1];
rz(-0.27401383) q[3];
sx q[3];
rz(-0.7139132) q[3];
sx q[3];
rz(1.4307724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.37904) q[2];
sx q[2];
rz(-2.5950044) q[2];
sx q[2];
rz(-2.9075882) q[2];
rz(1.0903357) q[3];
sx q[3];
rz(-1.5666311) q[3];
sx q[3];
rz(2.0719297) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91319084) q[0];
sx q[0];
rz(-2.985432) q[0];
sx q[0];
rz(-1.8432023) q[0];
rz(0.78650728) q[1];
sx q[1];
rz(-0.85580099) q[1];
sx q[1];
rz(-1.7826805) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3284861) q[0];
sx q[0];
rz(-2.4755602) q[0];
sx q[0];
rz(0.12909992) q[0];
rz(-pi) q[1];
rz(0.69245371) q[2];
sx q[2];
rz(-2.7992749) q[2];
sx q[2];
rz(-2.3794425) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37005597) q[1];
sx q[1];
rz(-0.0020469804) q[1];
sx q[1];
rz(-2.5084346) q[1];
rz(2.5547736) q[3];
sx q[3];
rz(-0.53884655) q[3];
sx q[3];
rz(-0.19308819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7274373) q[2];
sx q[2];
rz(-1.4667908) q[2];
sx q[2];
rz(-2.9690913) q[2];
rz(2.6532459) q[3];
sx q[3];
rz(-1.9988873) q[3];
sx q[3];
rz(2.02777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97335029) q[0];
sx q[0];
rz(-2.2689447) q[0];
sx q[0];
rz(0.33997047) q[0];
rz(-2.8151457) q[1];
sx q[1];
rz(-0.68845922) q[1];
sx q[1];
rz(-1.2350157) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8064855) q[0];
sx q[0];
rz(-2.3705784) q[0];
sx q[0];
rz(2.5946027) q[0];
rz(-2.4900774) q[2];
sx q[2];
rz(-0.55547248) q[2];
sx q[2];
rz(3.0998203) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81081796) q[1];
sx q[1];
rz(-2.0502809) q[1];
sx q[1];
rz(1.5768087) q[1];
x q[2];
rz(2.6458287) q[3];
sx q[3];
rz(-0.44124441) q[3];
sx q[3];
rz(-3.038681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0299224) q[2];
sx q[2];
rz(-2.3387574) q[2];
sx q[2];
rz(2.5862582) q[2];
rz(0.16767821) q[3];
sx q[3];
rz(-1.5372814) q[3];
sx q[3];
rz(-0.47784561) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70917201) q[0];
sx q[0];
rz(-0.92996159) q[0];
sx q[0];
rz(-2.0322556) q[0];
rz(1.1322016) q[1];
sx q[1];
rz(-0.39674509) q[1];
sx q[1];
rz(-0.94165492) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1718194) q[0];
sx q[0];
rz(-2.1226774) q[0];
sx q[0];
rz(0.18751796) q[0];
rz(-pi) q[1];
rz(-0.15502013) q[2];
sx q[2];
rz(-1.2885254) q[2];
sx q[2];
rz(-1.0978497) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4577427) q[1];
sx q[1];
rz(-2.2096429) q[1];
sx q[1];
rz(-0.15775494) q[1];
rz(-pi) q[2];
rz(1.3579943) q[3];
sx q[3];
rz(-1.9486685) q[3];
sx q[3];
rz(1.4546284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38810101) q[2];
sx q[2];
rz(-1.0270303) q[2];
sx q[2];
rz(1.5388185) q[2];
rz(-1.7704891) q[3];
sx q[3];
rz(-1.1857827) q[3];
sx q[3];
rz(-0.051430844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082315363) q[0];
sx q[0];
rz(-0.61510724) q[0];
sx q[0];
rz(-2.1737461) q[0];
rz(-1.8702501) q[1];
sx q[1];
rz(-1.8495193) q[1];
sx q[1];
rz(0.06180067) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74154749) q[0];
sx q[0];
rz(-1.9572659) q[0];
sx q[0];
rz(0.58048141) q[0];
rz(1.3702622) q[2];
sx q[2];
rz(-1.3167183) q[2];
sx q[2];
rz(0.24607436) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3810515) q[1];
sx q[1];
rz(-0.88812056) q[1];
sx q[1];
rz(1.6939916) q[1];
rz(-pi) q[2];
rz(-0.33225191) q[3];
sx q[3];
rz(-1.2945172) q[3];
sx q[3];
rz(2.0656697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9515848) q[2];
sx q[2];
rz(-0.20415674) q[2];
sx q[2];
rz(-0.07587138) q[2];
rz(-1.1672945) q[3];
sx q[3];
rz(-2.0397525) q[3];
sx q[3];
rz(-2.8372138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9645914) q[0];
sx q[0];
rz(-0.27934203) q[0];
sx q[0];
rz(-2.8569073) q[0];
rz(0.074706569) q[1];
sx q[1];
rz(-1.9120522) q[1];
sx q[1];
rz(1.7237192) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95686326) q[0];
sx q[0];
rz(-1.7635582) q[0];
sx q[0];
rz(-0.75988976) q[0];
rz(-pi) q[1];
rz(0.14947628) q[2];
sx q[2];
rz(-1.0188531) q[2];
sx q[2];
rz(-1.9665444) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3813044) q[1];
sx q[1];
rz(-0.78236303) q[1];
sx q[1];
rz(-2.5336877) q[1];
rz(-pi) q[2];
rz(1.6154434) q[3];
sx q[3];
rz(-2.4952125) q[3];
sx q[3];
rz(0.97585362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3384) q[2];
sx q[2];
rz(-1.1407547) q[2];
sx q[2];
rz(1.2791862) q[2];
rz(0.98215669) q[3];
sx q[3];
rz(-1.4582062) q[3];
sx q[3];
rz(2.2568683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2411156) q[0];
sx q[0];
rz(-1.3237088) q[0];
sx q[0];
rz(1.5644912) q[0];
rz(-1.0152394) q[1];
sx q[1];
rz(-1.1493586) q[1];
sx q[1];
rz(-1.1372067) q[1];
rz(-1.3080636) q[2];
sx q[2];
rz(-1.4663525) q[2];
sx q[2];
rz(2.732085) q[2];
rz(-0.7393403) q[3];
sx q[3];
rz(-1.1237595) q[3];
sx q[3];
rz(2.9343284) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

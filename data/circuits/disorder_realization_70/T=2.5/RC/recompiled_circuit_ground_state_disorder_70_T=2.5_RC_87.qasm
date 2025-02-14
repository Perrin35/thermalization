OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.81808972) q[0];
sx q[0];
rz(-0.40677318) q[0];
sx q[0];
rz(0.50954252) q[0];
rz(2.8264363) q[1];
sx q[1];
rz(-0.1935614) q[1];
sx q[1];
rz(-1.0055746) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60765841) q[0];
sx q[0];
rz(-2.9276121) q[0];
sx q[0];
rz(2.7817552) q[0];
x q[1];
rz(2.2174913) q[2];
sx q[2];
rz(-1.2900724) q[2];
sx q[2];
rz(-1.4626056) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7889346) q[1];
sx q[1];
rz(-0.96436687) q[1];
sx q[1];
rz(-0.85373099) q[1];
rz(-pi) q[2];
rz(-0.77862424) q[3];
sx q[3];
rz(-1.1168241) q[3];
sx q[3];
rz(0.046128143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.84844184) q[2];
sx q[2];
rz(-1.8895431) q[2];
sx q[2];
rz(-1.1085917) q[2];
rz(-2.0075924) q[3];
sx q[3];
rz(-1.0568551) q[3];
sx q[3];
rz(-0.081324287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14652458) q[0];
sx q[0];
rz(-2.0159371) q[0];
sx q[0];
rz(2.8296237) q[0];
rz(2.9106855) q[1];
sx q[1];
rz(-2.0377908) q[1];
sx q[1];
rz(-1.1280967) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9474079) q[0];
sx q[0];
rz(-1.1028521) q[0];
sx q[0];
rz(-0.2121966) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6083469) q[2];
sx q[2];
rz(-2.4002078) q[2];
sx q[2];
rz(-0.35517755) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23975325) q[1];
sx q[1];
rz(-1.1247207) q[1];
sx q[1];
rz(-1.6935936) q[1];
x q[2];
rz(-0.76448582) q[3];
sx q[3];
rz(-1.8967046) q[3];
sx q[3];
rz(-2.9284262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4565304) q[2];
sx q[2];
rz(-1.637746) q[2];
sx q[2];
rz(0.064229639) q[2];
rz(-1.8188933) q[3];
sx q[3];
rz(-0.6367681) q[3];
sx q[3];
rz(-0.88671154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5927758) q[0];
sx q[0];
rz(-2.9058457) q[0];
sx q[0];
rz(1.2491666) q[0];
rz(-2.8104172) q[1];
sx q[1];
rz(-1.1391897) q[1];
sx q[1];
rz(1.1963199) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6702658) q[0];
sx q[0];
rz(-1.571079) q[0];
sx q[0];
rz(3.1413881) q[0];
rz(-pi) q[1];
rz(-2.8941089) q[2];
sx q[2];
rz(-1.5759517) q[2];
sx q[2];
rz(2.2725819) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0532562) q[1];
sx q[1];
rz(-0.98829809) q[1];
sx q[1];
rz(-2.1964369) q[1];
rz(0.41091856) q[3];
sx q[3];
rz(-0.79686368) q[3];
sx q[3];
rz(1.7601907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.62200481) q[2];
sx q[2];
rz(-2.376611) q[2];
sx q[2];
rz(0.91274846) q[2];
rz(1.4384455) q[3];
sx q[3];
rz(-1.5104975) q[3];
sx q[3];
rz(1.8057757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66252935) q[0];
sx q[0];
rz(-1.1071858) q[0];
sx q[0];
rz(2.4712439) q[0];
rz(2.4743075) q[1];
sx q[1];
rz(-2.1332462) q[1];
sx q[1];
rz(-2.537312) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8025148) q[0];
sx q[0];
rz(-2.3065563) q[0];
sx q[0];
rz(-1.5123532) q[0];
rz(-pi) q[1];
rz(0.38197869) q[2];
sx q[2];
rz(-1.7151217) q[2];
sx q[2];
rz(1.195418) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.62897077) q[1];
sx q[1];
rz(-1.781946) q[1];
sx q[1];
rz(-1.8308012) q[1];
x q[2];
rz(1.8205582) q[3];
sx q[3];
rz(-0.74640025) q[3];
sx q[3];
rz(-1.5045034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0303354) q[2];
sx q[2];
rz(-1.5735441) q[2];
sx q[2];
rz(-2.7643909) q[2];
rz(-3.0180569) q[3];
sx q[3];
rz(-1.7081407) q[3];
sx q[3];
rz(-1.4089353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1971624) q[0];
sx q[0];
rz(-2.5640709) q[0];
sx q[0];
rz(-0.29931983) q[0];
rz(-2.0206644) q[1];
sx q[1];
rz(-1.9363554) q[1];
sx q[1];
rz(2.1790806) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0960353) q[0];
sx q[0];
rz(-1.0518414) q[0];
sx q[0];
rz(2.9039498) q[0];
x q[1];
rz(-0.99330866) q[2];
sx q[2];
rz(-2.3984512) q[2];
sx q[2];
rz(-0.20692839) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.968041) q[1];
sx q[1];
rz(-2.3690201) q[1];
sx q[1];
rz(-0.79141683) q[1];
rz(0.59904091) q[3];
sx q[3];
rz(-2.0192462) q[3];
sx q[3];
rz(-0.34235172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0472497) q[2];
sx q[2];
rz(-2.3361358) q[2];
sx q[2];
rz(2.1979525) q[2];
rz(-1.0654248) q[3];
sx q[3];
rz(-1.3506972) q[3];
sx q[3];
rz(1.0984727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0672673) q[0];
sx q[0];
rz(-1.4370947) q[0];
sx q[0];
rz(-3.1332916) q[0];
rz(-2.6240194) q[1];
sx q[1];
rz(-2.4868496) q[1];
sx q[1];
rz(-1.5483206) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28866523) q[0];
sx q[0];
rz(-1.2199243) q[0];
sx q[0];
rz(3.1159668) q[0];
rz(-2.6363434) q[2];
sx q[2];
rz(-1.8191686) q[2];
sx q[2];
rz(-2.5587683) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1344117) q[1];
sx q[1];
rz(-2.6406277) q[1];
sx q[1];
rz(2.3208614) q[1];
x q[2];
rz(-2.1861784) q[3];
sx q[3];
rz(-1.9166167) q[3];
sx q[3];
rz(2.7503848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3090618) q[2];
sx q[2];
rz(-1.4893724) q[2];
sx q[2];
rz(-1.3365411) q[2];
rz(-0.055770326) q[3];
sx q[3];
rz(-2.1499108) q[3];
sx q[3];
rz(-0.43558863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.5585153) q[0];
sx q[0];
rz(-2.6641088) q[0];
sx q[0];
rz(0.68786311) q[0];
rz(1.6784003) q[1];
sx q[1];
rz(-2.4122767) q[1];
sx q[1];
rz(-2.8942143) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5106544) q[0];
sx q[0];
rz(-1.3888089) q[0];
sx q[0];
rz(-2.8923558) q[0];
rz(-pi) q[1];
rz(2.192657) q[2];
sx q[2];
rz(-2.1315711) q[2];
sx q[2];
rz(1.6552629) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2823729) q[1];
sx q[1];
rz(-2.3973208) q[1];
sx q[1];
rz(-0.93282916) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.081253) q[3];
sx q[3];
rz(-2.8791028) q[3];
sx q[3];
rz(2.9406527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4443724) q[2];
sx q[2];
rz(-1.3345382) q[2];
sx q[2];
rz(-2.7174301) q[2];
rz(1.0423202) q[3];
sx q[3];
rz(-0.76516953) q[3];
sx q[3];
rz(-0.55646363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1384077) q[0];
sx q[0];
rz(-2.1557032) q[0];
sx q[0];
rz(0.61088046) q[0];
rz(-2.8914087) q[1];
sx q[1];
rz(-1.4004204) q[1];
sx q[1];
rz(-0.47666034) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2395916) q[0];
sx q[0];
rz(-1.6098813) q[0];
sx q[0];
rz(-2.0530967) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.119461) q[2];
sx q[2];
rz(-0.20344606) q[2];
sx q[2];
rz(0.97285336) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6338773) q[1];
sx q[1];
rz(-2.4832279) q[1];
sx q[1];
rz(-2.4509199) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1564756) q[3];
sx q[3];
rz(-1.9352311) q[3];
sx q[3];
rz(-1.5353312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.99159795) q[2];
sx q[2];
rz(-2.1074882) q[2];
sx q[2];
rz(2.4617713) q[2];
rz(-2.6796135) q[3];
sx q[3];
rz(-1.446412) q[3];
sx q[3];
rz(-1.6405039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3221472) q[0];
sx q[0];
rz(-1.3752022) q[0];
sx q[0];
rz(2.9875901) q[0];
rz(-0.18140659) q[1];
sx q[1];
rz(-2.0712974) q[1];
sx q[1];
rz(-0.12012404) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2699566) q[0];
sx q[0];
rz(-1.1678732) q[0];
sx q[0];
rz(3.0771034) q[0];
rz(-2.207802) q[2];
sx q[2];
rz(-0.32378886) q[2];
sx q[2];
rz(-2.0227014) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6473173) q[1];
sx q[1];
rz(-2.9396012) q[1];
sx q[1];
rz(1.0624733) q[1];
rz(0.52773962) q[3];
sx q[3];
rz(-1.4459918) q[3];
sx q[3];
rz(-3.0230056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4274365) q[2];
sx q[2];
rz(-1.7566046) q[2];
sx q[2];
rz(-0.45905054) q[2];
rz(-2.6311724) q[3];
sx q[3];
rz(-0.84158689) q[3];
sx q[3];
rz(-2.4418805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0521312) q[0];
sx q[0];
rz(-0.94877807) q[0];
sx q[0];
rz(2.7556038) q[0];
rz(-1.5929068) q[1];
sx q[1];
rz(-2.6451151) q[1];
sx q[1];
rz(1.5492424) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.693394) q[0];
sx q[0];
rz(-1.5992461) q[0];
sx q[0];
rz(1.7398119) q[0];
rz(0.73076325) q[2];
sx q[2];
rz(-0.77789069) q[2];
sx q[2];
rz(-0.90532263) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1116699) q[1];
sx q[1];
rz(-0.1055461) q[1];
sx q[1];
rz(-1.3131857) q[1];
rz(-pi) q[2];
rz(-1.4172961) q[3];
sx q[3];
rz(-1.5046538) q[3];
sx q[3];
rz(1.9754174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7768895) q[2];
sx q[2];
rz(-0.93901712) q[2];
sx q[2];
rz(-1.4466064) q[2];
rz(-2.1052836) q[3];
sx q[3];
rz(-0.88007897) q[3];
sx q[3];
rz(-3.115263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25630367) q[0];
sx q[0];
rz(-0.97923179) q[0];
sx q[0];
rz(-1.6543065) q[0];
rz(2.9215095) q[1];
sx q[1];
rz(-1.1047803) q[1];
sx q[1];
rz(1.6688375) q[1];
rz(0.34192495) q[2];
sx q[2];
rz(-2.3734063) q[2];
sx q[2];
rz(-2.0362324) q[2];
rz(-0.55752936) q[3];
sx q[3];
rz(-2.5685608) q[3];
sx q[3];
rz(1.3309042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

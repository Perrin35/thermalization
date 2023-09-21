OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2064535) q[0];
sx q[0];
rz(-0.78092617) q[0];
sx q[0];
rz(2.934802) q[0];
rz(-2.7517154) q[1];
sx q[1];
rz(-2.0808527) q[1];
sx q[1];
rz(0.18016711) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6102403) q[0];
sx q[0];
rz(-2.5533479) q[0];
sx q[0];
rz(0.63740762) q[0];
x q[1];
rz(-0.51282672) q[2];
sx q[2];
rz(-1.4564118) q[2];
sx q[2];
rz(-1.0637384) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85045056) q[1];
sx q[1];
rz(-2.2524815) q[1];
sx q[1];
rz(-0.095231685) q[1];
x q[2];
rz(2.4597635) q[3];
sx q[3];
rz(-1.4273943) q[3];
sx q[3];
rz(1.6553866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41123286) q[2];
sx q[2];
rz(-2.2683472) q[2];
sx q[2];
rz(1.8475378) q[2];
rz(2.7358352) q[3];
sx q[3];
rz(-1.5016705) q[3];
sx q[3];
rz(-0.40677235) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2186573) q[0];
sx q[0];
rz(-2.1056392) q[0];
sx q[0];
rz(3.0157715) q[0];
rz(-0.80548349) q[1];
sx q[1];
rz(-2.3352354) q[1];
sx q[1];
rz(-1.3719826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67324154) q[0];
sx q[0];
rz(-1.6325608) q[0];
sx q[0];
rz(-0.67586918) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20184529) q[2];
sx q[2];
rz(-2.5864374) q[2];
sx q[2];
rz(1.3943878) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0731197) q[1];
sx q[1];
rz(-1.8205376) q[1];
sx q[1];
rz(1.6176893) q[1];
rz(-pi) q[2];
rz(-2.6729229) q[3];
sx q[3];
rz(-1.7892924) q[3];
sx q[3];
rz(3.1359429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8877318) q[2];
sx q[2];
rz(-0.68887201) q[2];
sx q[2];
rz(0.99622336) q[2];
rz(-2.0455202) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(0.85038275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.6753733) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(-1.0590142) q[0];
rz(1.1478708) q[1];
sx q[1];
rz(-2.4025326) q[1];
sx q[1];
rz(-2.0770729) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36944593) q[0];
sx q[0];
rz(-1.1755953) q[0];
sx q[0];
rz(-2.4482083) q[0];
rz(0.15757615) q[2];
sx q[2];
rz(-1.9489947) q[2];
sx q[2];
rz(0.45730293) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4050582) q[1];
sx q[1];
rz(-1.5591803) q[1];
sx q[1];
rz(-1.2606773) q[1];
x q[2];
rz(0.84633175) q[3];
sx q[3];
rz(-1.4457821) q[3];
sx q[3];
rz(-0.66305977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.069783) q[2];
sx q[2];
rz(-1.1867384) q[2];
sx q[2];
rz(-2.8783669) q[2];
rz(-2.0227382) q[3];
sx q[3];
rz(-1.8656732) q[3];
sx q[3];
rz(-2.3222205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082821) q[0];
sx q[0];
rz(-3.0968956) q[0];
sx q[0];
rz(-1.3942962) q[0];
rz(1.0385723) q[1];
sx q[1];
rz(-1.4515406) q[1];
sx q[1];
rz(2.0746453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2571714) q[0];
sx q[0];
rz(-1.0413678) q[0];
sx q[0];
rz(3.0966395) q[0];
x q[1];
rz(0.26950403) q[2];
sx q[2];
rz(-1.4925033) q[2];
sx q[2];
rz(1.6646202) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2414788) q[1];
sx q[1];
rz(-0.80726868) q[1];
sx q[1];
rz(0.57358731) q[1];
rz(-1.1569818) q[3];
sx q[3];
rz(-2.6226351) q[3];
sx q[3];
rz(-1.3299691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2944494) q[2];
sx q[2];
rz(-1.3501945) q[2];
sx q[2];
rz(-0.28953826) q[2];
rz(-0.52044049) q[3];
sx q[3];
rz(-1.8013022) q[3];
sx q[3];
rz(2.382544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(1.6510058) q[0];
sx q[0];
rz(-0.0087954272) q[0];
sx q[0];
rz(-0.83475137) q[0];
rz(-1.897215) q[1];
sx q[1];
rz(-1.2507739) q[1];
sx q[1];
rz(1.429819) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5127038) q[0];
sx q[0];
rz(-1.5682724) q[0];
sx q[0];
rz(0.062688962) q[0];
rz(0.69447563) q[2];
sx q[2];
rz(-0.36437964) q[2];
sx q[2];
rz(2.7642872) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6068078) q[1];
sx q[1];
rz(-0.77960289) q[1];
sx q[1];
rz(-1.603655) q[1];
rz(-1.3423213) q[3];
sx q[3];
rz(-0.67111525) q[3];
sx q[3];
rz(-1.445678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66118583) q[2];
sx q[2];
rz(-2.1358261) q[2];
sx q[2];
rz(-1.2832114) q[2];
rz(-3.0094106) q[3];
sx q[3];
rz(-0.32871267) q[3];
sx q[3];
rz(-0.33974084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4955687) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(2.6640889) q[0];
rz(-1.6409138) q[1];
sx q[1];
rz(-1.6093107) q[1];
sx q[1];
rz(0.57055155) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30734277) q[0];
sx q[0];
rz(-1.2120005) q[0];
sx q[0];
rz(2.8310199) q[0];
x q[1];
rz(1.9448026) q[2];
sx q[2];
rz(-1.9960969) q[2];
sx q[2];
rz(-0.60356319) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.693336) q[1];
sx q[1];
rz(-0.92474557) q[1];
sx q[1];
rz(-2.4962884) q[1];
rz(1.5485351) q[3];
sx q[3];
rz(-0.92261693) q[3];
sx q[3];
rz(-2.9472586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.70696124) q[2];
sx q[2];
rz(-1.0501477) q[2];
sx q[2];
rz(-0.79664191) q[2];
rz(0.32026511) q[3];
sx q[3];
rz(-1.0864778) q[3];
sx q[3];
rz(1.2789352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0758078) q[0];
sx q[0];
rz(-1.0362352) q[0];
sx q[0];
rz(2.7835223) q[0];
rz(2.8529196) q[1];
sx q[1];
rz(-0.52983785) q[1];
sx q[1];
rz(1.3105062) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.34328) q[0];
sx q[0];
rz(-2.5710921) q[0];
sx q[0];
rz(-2.7891854) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99408044) q[2];
sx q[2];
rz(-2.1834282) q[2];
sx q[2];
rz(0.47555579) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9999034) q[1];
sx q[1];
rz(-2.018376) q[1];
sx q[1];
rz(-0.56772851) q[1];
rz(-pi) q[2];
rz(-0.88498022) q[3];
sx q[3];
rz(-1.0776057) q[3];
sx q[3];
rz(-0.54515884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.33401176) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(2.1441148) q[2];
rz(-2.5618662) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(1.4355481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8758133) q[0];
sx q[0];
rz(-1.4900603) q[0];
sx q[0];
rz(-2.8009801) q[0];
rz(1.2365201) q[1];
sx q[1];
rz(-0.75729901) q[1];
sx q[1];
rz(-1.9514726) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1163568) q[0];
sx q[0];
rz(-1.802117) q[0];
sx q[0];
rz(0.063409253) q[0];
x q[1];
rz(-0.868452) q[2];
sx q[2];
rz(-0.10570279) q[2];
sx q[2];
rz(-0.58123523) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5684143) q[1];
sx q[1];
rz(-1.8548994) q[1];
sx q[1];
rz(1.2969639) q[1];
rz(-pi) q[2];
rz(1.0134775) q[3];
sx q[3];
rz(-1.1790457) q[3];
sx q[3];
rz(-1.980892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.74165806) q[2];
sx q[2];
rz(-0.88230336) q[2];
sx q[2];
rz(-0.41440543) q[2];
rz(-1.3828145) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25093108) q[0];
sx q[0];
rz(-2.119976) q[0];
sx q[0];
rz(2.9898306) q[0];
rz(1.7157308) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(-0.94917667) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5700891) q[0];
sx q[0];
rz(-2.6766258) q[0];
sx q[0];
rz(-0.31682195) q[0];
rz(0.62990909) q[2];
sx q[2];
rz(-1.3917149) q[2];
sx q[2];
rz(-1.5392898) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6712499) q[1];
sx q[1];
rz(-1.3312695) q[1];
sx q[1];
rz(-2.4688979) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2373958) q[3];
sx q[3];
rz(-0.70611806) q[3];
sx q[3];
rz(1.0977942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7342547) q[2];
sx q[2];
rz(-0.39525017) q[2];
sx q[2];
rz(-2.7091743) q[2];
rz(-1.5405103) q[3];
sx q[3];
rz(-1.6316905) q[3];
sx q[3];
rz(-0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0703053) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(-2.7684257) q[0];
rz(0.35692731) q[1];
sx q[1];
rz(-0.86078763) q[1];
sx q[1];
rz(2.4180791) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05941885) q[0];
sx q[0];
rz(-0.84616236) q[0];
sx q[0];
rz(-1.3979777) q[0];
x q[1];
rz(1.4805484) q[2];
sx q[2];
rz(-2.2482578) q[2];
sx q[2];
rz(1.2573164) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.53612667) q[1];
sx q[1];
rz(-1.8408753) q[1];
sx q[1];
rz(0.77401604) q[1];
rz(-pi) q[2];
rz(-0.27042737) q[3];
sx q[3];
rz(-1.5661097) q[3];
sx q[3];
rz(-3.1386496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.92131203) q[2];
sx q[2];
rz(-1.7582515) q[2];
sx q[2];
rz(-0.55345654) q[2];
rz(-2.4297595) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(-2.0994983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9217459) q[0];
sx q[0];
rz(-0.60315673) q[0];
sx q[0];
rz(-3.0043816) q[0];
rz(-2.8593821) q[1];
sx q[1];
rz(-0.78411513) q[1];
sx q[1];
rz(0.93906739) q[1];
rz(0.7278022) q[2];
sx q[2];
rz(-2.5143378) q[2];
sx q[2];
rz(0.85744748) q[2];
rz(0.59755748) q[3];
sx q[3];
rz(-2.6664824) q[3];
sx q[3];
rz(0.53371724) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

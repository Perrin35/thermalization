OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6157827) q[0];
sx q[0];
rz(-1.4178185) q[0];
sx q[0];
rz(0.56086993) q[0];
rz(-2.0286735) q[1];
sx q[1];
rz(-1.3781883) q[1];
sx q[1];
rz(-1.2150432) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44170609) q[0];
sx q[0];
rz(-0.21472782) q[0];
sx q[0];
rz(2.0798652) q[0];
rz(-pi) q[1];
rz(0.36429976) q[2];
sx q[2];
rz(-0.69395739) q[2];
sx q[2];
rz(-2.7147164) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.558555) q[1];
sx q[1];
rz(-1.2416632) q[1];
sx q[1];
rz(-1.0298883) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0753161) q[3];
sx q[3];
rz(-0.30116943) q[3];
sx q[3];
rz(2.1198213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3540196) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-2.9585178) q[2];
rz(0.37781528) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.29782444) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(3.0644754) q[0];
rz(2.8027957) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(-1.6024626) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8388226) q[0];
sx q[0];
rz(-1.4937703) q[0];
sx q[0];
rz(2.1588615) q[0];
rz(2.402926) q[2];
sx q[2];
rz(-1.7324442) q[2];
sx q[2];
rz(-0.48036286) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7669249) q[1];
sx q[1];
rz(-0.51480773) q[1];
sx q[1];
rz(-1.9306081) q[1];
rz(-pi) q[2];
rz(-1.3354524) q[3];
sx q[3];
rz(-1.1574405) q[3];
sx q[3];
rz(1.2873161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(0.65845931) q[2];
rz(2.9902839) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(-0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028458683) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(-0.43310305) q[0];
rz(-1.9494879) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(0.55535299) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0387602) q[0];
sx q[0];
rz(-2.0534678) q[0];
sx q[0];
rz(2.32248) q[0];
x q[1];
rz(1.0054587) q[2];
sx q[2];
rz(-1.8381422) q[2];
sx q[2];
rz(-0.29758673) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8086116) q[1];
sx q[1];
rz(-2.4858027) q[1];
sx q[1];
rz(-1.8379184) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.025612) q[3];
sx q[3];
rz(-1.7435929) q[3];
sx q[3];
rz(-1.2196513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4042523) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(-1.8910485) q[2];
rz(-0.2441497) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(1.4499433) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26043949) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(2.326791) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(-0.25517685) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2164453) q[0];
sx q[0];
rz(-1.2108004) q[0];
sx q[0];
rz(1.2544592) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5023979) q[2];
sx q[2];
rz(-1.8393469) q[2];
sx q[2];
rz(-3.064379) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.081101) q[1];
sx q[1];
rz(-1.2730518) q[1];
sx q[1];
rz(-0.22025073) q[1];
rz(-1.049794) q[3];
sx q[3];
rz(-1.5265326) q[3];
sx q[3];
rz(-1.1417768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2531551) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(-0.17318428) q[2];
rz(0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(-0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2816876) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(1.7657071) q[0];
rz(-1.2777404) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(0.056093562) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.010904) q[0];
sx q[0];
rz(-2.4405257) q[0];
sx q[0];
rz(2.5581193) q[0];
rz(-pi) q[1];
x q[1];
rz(2.202583) q[2];
sx q[2];
rz(-2.5917705) q[2];
sx q[2];
rz(0.75013559) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.034654831) q[1];
sx q[1];
rz(-1.5562623) q[1];
sx q[1];
rz(2.8388139) q[1];
rz(-pi) q[2];
rz(2.131093) q[3];
sx q[3];
rz(-1.4951402) q[3];
sx q[3];
rz(1.8492941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(0.43219217) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0284001) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(-1.8796857) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(2.1870959) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.944825) q[0];
sx q[0];
rz(-1.4155354) q[0];
sx q[0];
rz(-1.0374271) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2013821) q[2];
sx q[2];
rz(-1.5427089) q[2];
sx q[2];
rz(2.0036151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2973605) q[1];
sx q[1];
rz(-2.4095222) q[1];
sx q[1];
rz(-2.3872603) q[1];
x q[2];
rz(1.2797221) q[3];
sx q[3];
rz(-1.4143412) q[3];
sx q[3];
rz(-0.94369704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.548617) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(2.0992289) q[2];
rz(2.7029165) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(-1.3180805) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8577268) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(-2.3983811) q[0];
rz(-1.6339533) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(2.5315703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3726495) q[0];
sx q[0];
rz(-2.5805051) q[0];
sx q[0];
rz(0.48604301) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1158175) q[2];
sx q[2];
rz(-0.91997416) q[2];
sx q[2];
rz(1.4027558) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4808018) q[1];
sx q[1];
rz(-1.5457488) q[1];
sx q[1];
rz(-0.90344306) q[1];
rz(-1.1915019) q[3];
sx q[3];
rz(-1.3790352) q[3];
sx q[3];
rz(-0.8134884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(0.91840333) q[2];
rz(-1.5911128) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(-0.38890719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.780705) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(-1.5135182) q[0];
rz(0.52945119) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(-2.4050074) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026222762) q[0];
sx q[0];
rz(-1.6008458) q[0];
sx q[0];
rz(-0.010859246) q[0];
rz(-1.2543711) q[2];
sx q[2];
rz(-2.8036615) q[2];
sx q[2];
rz(-0.28883176) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8383933) q[1];
sx q[1];
rz(-2.0337078) q[1];
sx q[1];
rz(1.4569015) q[1];
x q[2];
rz(-2.1804817) q[3];
sx q[3];
rz(-2.612252) q[3];
sx q[3];
rz(1.6314268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0344051) q[2];
sx q[2];
rz(-1.9341058) q[2];
sx q[2];
rz(0.68391189) q[2];
rz(-1.9125787) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1972315) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(-0.18572447) q[0];
rz(0.99705237) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(2.396778) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1574402) q[0];
sx q[0];
rz(-2.3072349) q[0];
sx q[0];
rz(-2.0122583) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3955598) q[2];
sx q[2];
rz(-2.308508) q[2];
sx q[2];
rz(1.2566483) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8150755) q[1];
sx q[1];
rz(-2.7217367) q[1];
sx q[1];
rz(-0.99680568) q[1];
x q[2];
rz(-0.39448491) q[3];
sx q[3];
rz(-2.3186765) q[3];
sx q[3];
rz(-0.2713954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1054489) q[2];
sx q[2];
rz(-2.3858586) q[2];
sx q[2];
rz(1.194681) q[2];
rz(2.1448994) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(-2.1452346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.3982518) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(2.7375896) q[0];
rz(-0.031127302) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(-1.9706479) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084328018) q[0];
sx q[0];
rz(-1.143976) q[0];
sx q[0];
rz(-1.5945934) q[0];
x q[1];
rz(-0.2756341) q[2];
sx q[2];
rz(-2.3466957) q[2];
sx q[2];
rz(2.8385712) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3552637) q[1];
sx q[1];
rz(-1.5794282) q[1];
sx q[1];
rz(1.5481871) q[1];
rz(-pi) q[2];
rz(1.2788494) q[3];
sx q[3];
rz(-1.9817096) q[3];
sx q[3];
rz(1.7850072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6955473) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(0.5919624) q[2];
rz(-2.5752318) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(-1.5238354) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82407172) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(-3.042165) q[1];
sx q[1];
rz(-1.8933404) q[1];
sx q[1];
rz(1.0642687) q[1];
rz(-0.70384937) q[2];
sx q[2];
rz(-2.2675632) q[2];
sx q[2];
rz(2.8013196) q[2];
rz(-1.6173784) q[3];
sx q[3];
rz(-2.8077879) q[3];
sx q[3];
rz(2.2624364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

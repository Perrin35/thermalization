OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47867632) q[0];
sx q[0];
rz(7.8362099) q[0];
sx q[0];
rz(10.683164) q[0];
rz(-5.0212669) q[1];
sx q[1];
rz(6.8016383) q[1];
sx q[1];
rz(6.3432884) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10190554) q[0];
sx q[0];
rz(-1.6015669) q[0];
sx q[0];
rz(-3.0512848) q[0];
x q[1];
rz(-1.41729) q[2];
sx q[2];
rz(-2.058284) q[2];
sx q[2];
rz(1.1660898) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93906885) q[1];
sx q[1];
rz(-1.2092522) q[1];
sx q[1];
rz(-0.099379813) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1677371) q[3];
sx q[3];
rz(-1.1789448) q[3];
sx q[3];
rz(-0.69249707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1897159) q[2];
sx q[2];
rz(-2.3592301) q[2];
sx q[2];
rz(0.18048364) q[2];
rz(-2.8372724) q[3];
sx q[3];
rz(-2.2173827) q[3];
sx q[3];
rz(-0.2271823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3825398) q[0];
sx q[0];
rz(-2.5681684) q[0];
sx q[0];
rz(0.10391129) q[0];
rz(-2.3567764) q[1];
sx q[1];
rz(-1.3094614) q[1];
sx q[1];
rz(0.58194247) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4020445) q[0];
sx q[0];
rz(-2.6091837) q[0];
sx q[0];
rz(0.043921434) q[0];
rz(1.9230827) q[2];
sx q[2];
rz(-1.0880044) q[2];
sx q[2];
rz(-0.60868516) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71680333) q[1];
sx q[1];
rz(-1.5967073) q[1];
sx q[1];
rz(-1.0938563) q[1];
x q[2];
rz(-0.76762565) q[3];
sx q[3];
rz(-2.5968142) q[3];
sx q[3];
rz(1.2380074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71501422) q[2];
sx q[2];
rz(-2.7535186) q[2];
sx q[2];
rz(1.0283872) q[2];
rz(2.5905124) q[3];
sx q[3];
rz(-1.9469399) q[3];
sx q[3];
rz(-2.6330131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5904163) q[0];
sx q[0];
rz(-1.506378) q[0];
sx q[0];
rz(-1.2580385) q[0];
rz(-0.60351562) q[1];
sx q[1];
rz(-1.8305093) q[1];
sx q[1];
rz(0.18361941) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98417898) q[0];
sx q[0];
rz(-1.3065728) q[0];
sx q[0];
rz(-1.9963032) q[0];
rz(0.48086353) q[2];
sx q[2];
rz(-1.4521952) q[2];
sx q[2];
rz(-3.1080217) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.078122) q[1];
sx q[1];
rz(-2.2666711) q[1];
sx q[1];
rz(2.1433926) q[1];
x q[2];
rz(2.9345064) q[3];
sx q[3];
rz(-1.6112176) q[3];
sx q[3];
rz(-2.25349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6151578) q[2];
sx q[2];
rz(-1.7650471) q[2];
sx q[2];
rz(-2.5395425) q[2];
rz(0.43271068) q[3];
sx q[3];
rz(-1.5475464) q[3];
sx q[3];
rz(-1.4759147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3848307) q[0];
sx q[0];
rz(-0.47465208) q[0];
sx q[0];
rz(0.62622825) q[0];
rz(1.2681883) q[1];
sx q[1];
rz(-1.405193) q[1];
sx q[1];
rz(0.38571206) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46661147) q[0];
sx q[0];
rz(-0.85133305) q[0];
sx q[0];
rz(-0.58746679) q[0];
x q[1];
rz(-2.0247632) q[2];
sx q[2];
rz(-1.0359284) q[2];
sx q[2];
rz(-1.7810389) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6059162) q[1];
sx q[1];
rz(-1.7711346) q[1];
sx q[1];
rz(1.6388333) q[1];
rz(-pi) q[2];
rz(-1.3703501) q[3];
sx q[3];
rz(-2.5894937) q[3];
sx q[3];
rz(0.31262661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4422153) q[2];
sx q[2];
rz(-2.0065887) q[2];
sx q[2];
rz(-0.29083148) q[2];
rz(1.1902827) q[3];
sx q[3];
rz(-0.61605993) q[3];
sx q[3];
rz(-2.0975838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46665835) q[0];
sx q[0];
rz(-3.0789154) q[0];
sx q[0];
rz(1.689893) q[0];
rz(0.85707227) q[1];
sx q[1];
rz(-1.8430201) q[1];
sx q[1];
rz(0.42795408) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8030386) q[0];
sx q[0];
rz(-2.0074275) q[0];
sx q[0];
rz(3.0019041) q[0];
x q[1];
rz(-0.54938118) q[2];
sx q[2];
rz(-2.3231988) q[2];
sx q[2];
rz(-0.1521509) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.9792069) q[1];
sx q[1];
rz(-0.89079327) q[1];
sx q[1];
rz(2.0335781) q[1];
x q[2];
rz(0.6537207) q[3];
sx q[3];
rz(-0.80721426) q[3];
sx q[3];
rz(1.2098055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2032623) q[2];
sx q[2];
rz(-2.0945956) q[2];
sx q[2];
rz(-2.9948575) q[2];
rz(2.9444368) q[3];
sx q[3];
rz(-1.4898172) q[3];
sx q[3];
rz(-0.76876918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41106221) q[0];
sx q[0];
rz(-1.0163607) q[0];
sx q[0];
rz(0.16832571) q[0];
rz(0.39237818) q[1];
sx q[1];
rz(-2.7057251) q[1];
sx q[1];
rz(1.4515152) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1010676) q[0];
sx q[0];
rz(-0.4372963) q[0];
sx q[0];
rz(-1.6213202) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2757569) q[2];
sx q[2];
rz(-0.74271281) q[2];
sx q[2];
rz(-2.4969522) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.299356) q[1];
sx q[1];
rz(-2.3122687) q[1];
sx q[1];
rz(-1.7620128) q[1];
rz(0.78877641) q[3];
sx q[3];
rz(-1.2433854) q[3];
sx q[3];
rz(-0.57306266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77439857) q[2];
sx q[2];
rz(-0.44782475) q[2];
sx q[2];
rz(-1.9642584) q[2];
rz(0.92159739) q[3];
sx q[3];
rz(-2.2955743) q[3];
sx q[3];
rz(2.6109076) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0717936) q[0];
sx q[0];
rz(-1.256218) q[0];
sx q[0];
rz(2.0846833) q[0];
rz(2.9131556) q[1];
sx q[1];
rz(-0.7799131) q[1];
sx q[1];
rz(0.25100073) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4514654) q[0];
sx q[0];
rz(-1.2862035) q[0];
sx q[0];
rz(-1.837912) q[0];
x q[1];
rz(-2.5152048) q[2];
sx q[2];
rz(-2.1715626) q[2];
sx q[2];
rz(-2.8383672) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31902203) q[1];
sx q[1];
rz(-2.9290556) q[1];
sx q[1];
rz(-2.4478649) q[1];
x q[2];
rz(-2.6282309) q[3];
sx q[3];
rz(-1.8979478) q[3];
sx q[3];
rz(-2.1029496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2061657) q[2];
sx q[2];
rz(-2.3212104) q[2];
sx q[2];
rz(-2.7247735) q[2];
rz(2.5721512) q[3];
sx q[3];
rz(-2.0408401) q[3];
sx q[3];
rz(-2.4429564) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57132974) q[0];
sx q[0];
rz(-1.2705734) q[0];
sx q[0];
rz(0.50462333) q[0];
rz(-0.2568256) q[1];
sx q[1];
rz(-2.3378614) q[1];
sx q[1];
rz(-1.9907192) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6136632) q[0];
sx q[0];
rz(-1.3656034) q[0];
sx q[0];
rz(-2.0733207) q[0];
rz(-2.2058401) q[2];
sx q[2];
rz(-1.2798736) q[2];
sx q[2];
rz(-2.3036602) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.013864036) q[1];
sx q[1];
rz(-1.999966) q[1];
sx q[1];
rz(-0.60529373) q[1];
rz(-0.41696291) q[3];
sx q[3];
rz(-1.4259669) q[3];
sx q[3];
rz(-0.24709283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.63200942) q[2];
sx q[2];
rz(-2.7183967) q[2];
sx q[2];
rz(0.43935856) q[2];
rz(1.1257233) q[3];
sx q[3];
rz(-1.6042387) q[3];
sx q[3];
rz(-1.2996947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87089649) q[0];
sx q[0];
rz(-0.82042158) q[0];
sx q[0];
rz(0.46723715) q[0];
rz(-2.1029419) q[1];
sx q[1];
rz(-2.6324582) q[1];
sx q[1];
rz(-1.4899303) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56541967) q[0];
sx q[0];
rz(-1.3955981) q[0];
sx q[0];
rz(0.498226) q[0];
rz(-pi) q[1];
rz(1.6655671) q[2];
sx q[2];
rz(-1.7168003) q[2];
sx q[2];
rz(0.80751792) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6257621) q[1];
sx q[1];
rz(-0.7554248) q[1];
sx q[1];
rz(-1.6095407) q[1];
rz(-pi) q[2];
rz(-2.0155679) q[3];
sx q[3];
rz(-2.091134) q[3];
sx q[3];
rz(2.570278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6963639) q[2];
sx q[2];
rz(-2.7859521) q[2];
sx q[2];
rz(-1.5544372) q[2];
rz(2.1571531) q[3];
sx q[3];
rz(-0.90960228) q[3];
sx q[3];
rz(-1.7279153) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5928891) q[0];
sx q[0];
rz(-2.2991572) q[0];
sx q[0];
rz(-0.62826759) q[0];
rz(-0.20333044) q[1];
sx q[1];
rz(-1.1140099) q[1];
sx q[1];
rz(0.59439739) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0258163) q[0];
sx q[0];
rz(-2.5133142) q[0];
sx q[0];
rz(1.216287) q[0];
x q[1];
rz(0.35697414) q[2];
sx q[2];
rz(-2.5737615) q[2];
sx q[2];
rz(0.50895509) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5446721) q[1];
sx q[1];
rz(-1.2697176) q[1];
sx q[1];
rz(-1.2872417) q[1];
x q[2];
rz(0.10108642) q[3];
sx q[3];
rz(-2.5261263) q[3];
sx q[3];
rz(-3.0973928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9336885) q[2];
sx q[2];
rz(-0.7242569) q[2];
sx q[2];
rz(-0.17624632) q[2];
rz(-1.8267953) q[3];
sx q[3];
rz(-2.3254471) q[3];
sx q[3];
rz(-0.76752457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.15277302) q[0];
sx q[0];
rz(-1.7012699) q[0];
sx q[0];
rz(-1.758601) q[0];
rz(3.026961) q[1];
sx q[1];
rz(-0.75420598) q[1];
sx q[1];
rz(2.4892714) q[1];
rz(0.23723142) q[2];
sx q[2];
rz(-1.5330412) q[2];
sx q[2];
rz(1.588149) q[2];
rz(-1.0788315) q[3];
sx q[3];
rz(-1.9117461) q[3];
sx q[3];
rz(-0.030464813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

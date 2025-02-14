OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3235029) q[0];
sx q[0];
rz(-2.7348195) q[0];
sx q[0];
rz(-0.50954252) q[0];
rz(2.8264363) q[1];
sx q[1];
rz(-0.1935614) q[1];
sx q[1];
rz(-1.0055746) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97516528) q[0];
sx q[0];
rz(-1.7708798) q[0];
sx q[0];
rz(1.4944264) q[0];
rz(0.92410134) q[2];
sx q[2];
rz(-1.8515203) q[2];
sx q[2];
rz(1.678987) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3444958) q[1];
sx q[1];
rz(-2.2386947) q[1];
sx q[1];
rz(-0.75842661) q[1];
rz(0.60723234) q[3];
sx q[3];
rz(-2.2650026) q[3];
sx q[3];
rz(1.1986365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2931508) q[2];
sx q[2];
rz(-1.8895431) q[2];
sx q[2];
rz(2.0330009) q[2];
rz(1.1340002) q[3];
sx q[3];
rz(-1.0568551) q[3];
sx q[3];
rz(3.0602684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9950681) q[0];
sx q[0];
rz(-1.1256555) q[0];
sx q[0];
rz(0.31196892) q[0];
rz(0.23090714) q[1];
sx q[1];
rz(-1.1038019) q[1];
sx q[1];
rz(2.013496) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8618516) q[0];
sx q[0];
rz(-1.759888) q[0];
sx q[0];
rz(2.0479338) q[0];
rz(-pi) q[1];
rz(-0.53324576) q[2];
sx q[2];
rz(-0.74138481) q[2];
sx q[2];
rz(-2.7864151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1027185) q[1];
sx q[1];
rz(-2.6800214) q[1];
sx q[1];
rz(-0.25074236) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76448582) q[3];
sx q[3];
rz(-1.8967046) q[3];
sx q[3];
rz(2.9284262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4565304) q[2];
sx q[2];
rz(-1.5038467) q[2];
sx q[2];
rz(3.077363) q[2];
rz(1.8188933) q[3];
sx q[3];
rz(-2.5048246) q[3];
sx q[3];
rz(-0.88671154) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54881683) q[0];
sx q[0];
rz(-2.9058457) q[0];
sx q[0];
rz(-1.2491666) q[0];
rz(-2.8104172) q[1];
sx q[1];
rz(-1.1391897) q[1];
sx q[1];
rz(-1.9452728) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0421233) q[0];
sx q[0];
rz(-1.5705918) q[0];
sx q[0];
rz(-1.5705137) q[0];
rz(-pi) q[1];
rz(-3.1205503) q[2];
sx q[2];
rz(-2.8940563) q[2];
sx q[2];
rz(-2.4602082) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1332449) q[1];
sx q[1];
rz(-2.3143594) q[1];
sx q[1];
rz(2.4148947) q[1];
rz(-pi) q[2];
rz(-0.41091856) q[3];
sx q[3];
rz(-0.79686368) q[3];
sx q[3];
rz(1.381402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5195878) q[2];
sx q[2];
rz(-2.376611) q[2];
sx q[2];
rz(-0.91274846) q[2];
rz(-1.4384455) q[3];
sx q[3];
rz(-1.6310952) q[3];
sx q[3];
rz(1.8057757) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4790633) q[0];
sx q[0];
rz(-2.0344069) q[0];
sx q[0];
rz(0.67034876) q[0];
rz(-0.66728512) q[1];
sx q[1];
rz(-1.0083464) q[1];
sx q[1];
rz(-0.60428062) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7155565) q[0];
sx q[0];
rz(-2.4039488) q[0];
sx q[0];
rz(-0.064427388) q[0];
x q[1];
rz(2.759614) q[2];
sx q[2];
rz(-1.7151217) q[2];
sx q[2];
rz(-1.195418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5126219) q[1];
sx q[1];
rz(-1.3596467) q[1];
sx q[1];
rz(-1.8308012) q[1];
rz(-pi) q[2];
rz(-0.84008645) q[3];
sx q[3];
rz(-1.7394251) q[3];
sx q[3];
rz(0.11883277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0303354) q[2];
sx q[2];
rz(-1.5735441) q[2];
sx q[2];
rz(-0.37720171) q[2];
rz(3.0180569) q[3];
sx q[3];
rz(-1.433452) q[3];
sx q[3];
rz(1.7326573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.1971624) q[0];
sx q[0];
rz(-0.57752174) q[0];
sx q[0];
rz(0.29931983) q[0];
rz(2.0206644) q[1];
sx q[1];
rz(-1.9363554) q[1];
sx q[1];
rz(-2.1790806) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40874562) q[0];
sx q[0];
rz(-0.56622177) q[0];
sx q[0];
rz(-1.179856) q[0];
rz(0.99330866) q[2];
sx q[2];
rz(-2.3984512) q[2];
sx q[2];
rz(-2.9346643) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1735517) q[1];
sx q[1];
rz(-0.7725726) q[1];
sx q[1];
rz(2.3501758) q[1];
rz(-pi) q[2];
rz(-2.0983134) q[3];
sx q[3];
rz(-1.0378222) q[3];
sx q[3];
rz(-1.5162374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0943429) q[2];
sx q[2];
rz(-0.80545682) q[2];
sx q[2];
rz(0.94364014) q[2];
rz(1.0654248) q[3];
sx q[3];
rz(-1.7908955) q[3];
sx q[3];
rz(-2.0431199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074325398) q[0];
sx q[0];
rz(-1.4370947) q[0];
sx q[0];
rz(3.1332916) q[0];
rz(-2.6240194) q[1];
sx q[1];
rz(-0.65474302) q[1];
sx q[1];
rz(1.5483206) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36309886) q[0];
sx q[0];
rz(-0.35176793) q[0];
sx q[0];
rz(-1.5009053) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50524925) q[2];
sx q[2];
rz(-1.3224241) q[2];
sx q[2];
rz(-2.5587683) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0071809) q[1];
sx q[1];
rz(-0.500965) q[1];
sx q[1];
rz(-2.3208614) q[1];
x q[2];
rz(0.95541422) q[3];
sx q[3];
rz(-1.2249759) q[3];
sx q[3];
rz(0.39120787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3090618) q[2];
sx q[2];
rz(-1.4893724) q[2];
sx q[2];
rz(-1.8050516) q[2];
rz(-3.0858223) q[3];
sx q[3];
rz(-2.1499108) q[3];
sx q[3];
rz(0.43558863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5585153) q[0];
sx q[0];
rz(-0.4774839) q[0];
sx q[0];
rz(-2.4537295) q[0];
rz(-1.4631924) q[1];
sx q[1];
rz(-2.4122767) q[1];
sx q[1];
rz(0.24737839) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63093827) q[0];
sx q[0];
rz(-1.7527837) q[0];
sx q[0];
rz(2.8923558) q[0];
x q[1];
rz(0.94893564) q[2];
sx q[2];
rz(-1.0100216) q[2];
sx q[2];
rz(-1.4863297) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2108442) q[1];
sx q[1];
rz(-1.9860876) q[1];
sx q[1];
rz(-2.2077435) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12567802) q[3];
sx q[3];
rz(-1.8018556) q[3];
sx q[3];
rz(-0.30320689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6972202) q[2];
sx q[2];
rz(-1.3345382) q[2];
sx q[2];
rz(2.7174301) q[2];
rz(2.0992725) q[3];
sx q[3];
rz(-2.3764231) q[3];
sx q[3];
rz(2.585129) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0031849) q[0];
sx q[0];
rz(-2.1557032) q[0];
sx q[0];
rz(2.5307122) q[0];
rz(0.25018397) q[1];
sx q[1];
rz(-1.4004204) q[1];
sx q[1];
rz(2.6649323) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7358897) q[0];
sx q[0];
rz(-2.6578356) q[0];
sx q[0];
rz(-1.4866845) q[0];
rz(1.3872434) q[2];
sx q[2];
rz(-1.4825562) q[2];
sx q[2];
rz(-2.1004408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6338773) q[1];
sx q[1];
rz(-2.4832279) q[1];
sx q[1];
rz(0.69067278) q[1];
rz(-pi) q[2];
rz(2.7122981) q[3];
sx q[3];
rz(-2.1135114) q[3];
sx q[3];
rz(-2.8739704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.99159795) q[2];
sx q[2];
rz(-2.1074882) q[2];
sx q[2];
rz(2.4617713) q[2];
rz(0.46197915) q[3];
sx q[3];
rz(-1.446412) q[3];
sx q[3];
rz(-1.6405039) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81944549) q[0];
sx q[0];
rz(-1.3752022) q[0];
sx q[0];
rz(-0.15400259) q[0];
rz(2.9601861) q[1];
sx q[1];
rz(-2.0712974) q[1];
sx q[1];
rz(-0.12012404) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2699566) q[0];
sx q[0];
rz(-1.9737195) q[0];
sx q[0];
rz(-3.0771034) q[0];
rz(-pi) q[1];
rz(1.8343049) q[2];
sx q[2];
rz(-1.7611836) q[2];
sx q[2];
rz(0.15979494) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5654345) q[1];
sx q[1];
rz(-1.6685969) q[1];
sx q[1];
rz(-1.7478155) q[1];
x q[2];
rz(-1.7150015) q[3];
sx q[3];
rz(-2.0940082) q[3];
sx q[3];
rz(1.5246403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7141562) q[2];
sx q[2];
rz(-1.7566046) q[2];
sx q[2];
rz(2.6825421) q[2];
rz(-0.51042026) q[3];
sx q[3];
rz(-0.84158689) q[3];
sx q[3];
rz(2.4418805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0521312) q[0];
sx q[0];
rz(-0.94877807) q[0];
sx q[0];
rz(-2.7556038) q[0];
rz(1.5929068) q[1];
sx q[1];
rz(-0.49647757) q[1];
sx q[1];
rz(-1.5923502) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-1.5992461) q[0];
sx q[0];
rz(1.4017808) q[0];
rz(-pi) q[1];
rz(2.5086918) q[2];
sx q[2];
rz(-1.0833246) q[2];
sx q[2];
rz(-0.097214708) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2889115) q[1];
sx q[1];
rz(-1.6728472) q[1];
sx q[1];
rz(3.1146088) q[1];
rz(-1.7242966) q[3];
sx q[3];
rz(-1.5046538) q[3];
sx q[3];
rz(-1.9754174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7768895) q[2];
sx q[2];
rz(-2.2025755) q[2];
sx q[2];
rz(1.4466064) q[2];
rz(-1.0363091) q[3];
sx q[3];
rz(-0.88007897) q[3];
sx q[3];
rz(3.115263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.885289) q[0];
sx q[0];
rz(-0.97923179) q[0];
sx q[0];
rz(-1.6543065) q[0];
rz(-2.9215095) q[1];
sx q[1];
rz(-2.0368123) q[1];
sx q[1];
rz(-1.4727551) q[1];
rz(2.7996677) q[2];
sx q[2];
rz(-0.76818633) q[2];
sx q[2];
rz(1.1053602) q[2];
rz(1.2418048) q[3];
sx q[3];
rz(-1.0927148) q[3];
sx q[3];
rz(0.69251251) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

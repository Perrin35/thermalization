OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0176528) q[0];
sx q[0];
rz(4.058429) q[0];
sx q[0];
rz(9.8596758) q[0];
rz(-1.544156) q[1];
sx q[1];
rz(-0.51900744) q[1];
sx q[1];
rz(-2.5595698) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8911444) q[0];
sx q[0];
rz(-1.8316557) q[0];
sx q[0];
rz(-2.9820739) q[0];
x q[1];
rz(0.97442128) q[2];
sx q[2];
rz(-2.4108464) q[2];
sx q[2];
rz(-1.3517018) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6710606) q[1];
sx q[1];
rz(-1.8166421) q[1];
sx q[1];
rz(-2.0162986) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9635756) q[3];
sx q[3];
rz(-1.2937163) q[3];
sx q[3];
rz(-1.2704364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.55328289) q[2];
sx q[2];
rz(-1.2574235) q[2];
sx q[2];
rz(-0.29169875) q[2];
rz(2.6907673) q[3];
sx q[3];
rz(-0.25142938) q[3];
sx q[3];
rz(-0.60744557) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4400103) q[0];
sx q[0];
rz(-0.99047438) q[0];
sx q[0];
rz(-2.0462346) q[0];
rz(-1.1913242) q[1];
sx q[1];
rz(-2.2215863) q[1];
sx q[1];
rz(-0.90257588) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0046834) q[0];
sx q[0];
rz(-2.1666514) q[0];
sx q[0];
rz(-0.76120241) q[0];
rz(-pi) q[1];
rz(-1.9368725) q[2];
sx q[2];
rz(-2.4395025) q[2];
sx q[2];
rz(-3.0099208) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2486156) q[1];
sx q[1];
rz(-1.4505523) q[1];
sx q[1];
rz(-0.29386947) q[1];
rz(-2.6379073) q[3];
sx q[3];
rz(-1.0102444) q[3];
sx q[3];
rz(0.47893804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9537182) q[2];
sx q[2];
rz(-0.97430054) q[2];
sx q[2];
rz(0.11889674) q[2];
rz(0.64905727) q[3];
sx q[3];
rz(-0.18638149) q[3];
sx q[3];
rz(1.4547179) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4897937) q[0];
sx q[0];
rz(-1.2540023) q[0];
sx q[0];
rz(2.6258262) q[0];
rz(-1.0924529) q[1];
sx q[1];
rz(-0.40320435) q[1];
sx q[1];
rz(1.8663503) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37947734) q[0];
sx q[0];
rz(-1.5526875) q[0];
sx q[0];
rz(0.041569592) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3408717) q[2];
sx q[2];
rz(-1.2650239) q[2];
sx q[2];
rz(1.1912322) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.470388) q[1];
sx q[1];
rz(-0.71031308) q[1];
sx q[1];
rz(1.2673668) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4687863) q[3];
sx q[3];
rz(-0.77931858) q[3];
sx q[3];
rz(1.9106227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0755997) q[2];
sx q[2];
rz(-1.3174572) q[2];
sx q[2];
rz(1.9817748) q[2];
rz(-2.6521111) q[3];
sx q[3];
rz(-1.5882086) q[3];
sx q[3];
rz(2.421853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2751665) q[0];
sx q[0];
rz(-0.32298276) q[0];
sx q[0];
rz(0.48165709) q[0];
rz(2.2109168) q[1];
sx q[1];
rz(-1.296867) q[1];
sx q[1];
rz(-3.1226588) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9809119) q[0];
sx q[0];
rz(-1.5917516) q[0];
sx q[0];
rz(1.5207855) q[0];
rz(-1.5140947) q[2];
sx q[2];
rz(-0.66326521) q[2];
sx q[2];
rz(-0.5723638) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0731813) q[1];
sx q[1];
rz(-0.8340237) q[1];
sx q[1];
rz(0.72552105) q[1];
x q[2];
rz(2.3036257) q[3];
sx q[3];
rz(-2.1893756) q[3];
sx q[3];
rz(-0.092242084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21861741) q[2];
sx q[2];
rz(-0.3564035) q[2];
sx q[2];
rz(0.74529988) q[2];
rz(0.37799147) q[3];
sx q[3];
rz(-1.9972921) q[3];
sx q[3];
rz(-2.2530344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9208263) q[0];
sx q[0];
rz(-0.31196088) q[0];
sx q[0];
rz(-1.1908603) q[0];
rz(2.6530755) q[1];
sx q[1];
rz(-0.91594511) q[1];
sx q[1];
rz(0.8078422) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4427748) q[0];
sx q[0];
rz(-1.4801868) q[0];
sx q[0];
rz(1.8987937) q[0];
x q[1];
rz(2.9132782) q[2];
sx q[2];
rz(-2.0234152) q[2];
sx q[2];
rz(-1.376898) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1792308) q[1];
sx q[1];
rz(-2.6768502) q[1];
sx q[1];
rz(0.19728139) q[1];
rz(-0.94781117) q[3];
sx q[3];
rz(-1.4075507) q[3];
sx q[3];
rz(-0.80006525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1390344) q[2];
sx q[2];
rz(-2.148874) q[2];
sx q[2];
rz(-2.9774418) q[2];
rz(-0.74603355) q[3];
sx q[3];
rz(-1.371871) q[3];
sx q[3];
rz(3.0677838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98110759) q[0];
sx q[0];
rz(-1.8758513) q[0];
sx q[0];
rz(0.72738457) q[0];
rz(0.47850594) q[1];
sx q[1];
rz(-1.6886657) q[1];
sx q[1];
rz(-2.4047638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19664431) q[0];
sx q[0];
rz(-1.6262615) q[0];
sx q[0];
rz(-1.4329628) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13104266) q[2];
sx q[2];
rz(-1.8497582) q[2];
sx q[2];
rz(0.38273465) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.7718862) q[1];
sx q[1];
rz(-0.73644887) q[1];
sx q[1];
rz(1.0095897) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9584452) q[3];
sx q[3];
rz(-1.5836925) q[3];
sx q[3];
rz(1.8356334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9937146) q[2];
sx q[2];
rz(-2.1174049) q[2];
sx q[2];
rz(-0.68515879) q[2];
rz(-1.0088751) q[3];
sx q[3];
rz(-2.4515371) q[3];
sx q[3];
rz(2.8907997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20550263) q[0];
sx q[0];
rz(-0.093955366) q[0];
sx q[0];
rz(1.093338) q[0];
rz(-3.1178442) q[1];
sx q[1];
rz(-0.7370342) q[1];
sx q[1];
rz(-2.645983) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6709061) q[0];
sx q[0];
rz(-1.5829594) q[0];
sx q[0];
rz(3.0755416) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43735403) q[2];
sx q[2];
rz(-1.3833356) q[2];
sx q[2];
rz(2.1700493) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.62140761) q[1];
sx q[1];
rz(-1.9193135) q[1];
sx q[1];
rz(1.7640616) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.43066671) q[3];
sx q[3];
rz(-2.0206631) q[3];
sx q[3];
rz(1.7783742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2451943) q[2];
sx q[2];
rz(-2.3485025) q[2];
sx q[2];
rz(-0.29655656) q[2];
rz(0.5101997) q[3];
sx q[3];
rz(-1.6161796) q[3];
sx q[3];
rz(-0.27142522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1475913) q[0];
sx q[0];
rz(-0.95518249) q[0];
sx q[0];
rz(1.1543132) q[0];
rz(1.9860024) q[1];
sx q[1];
rz(-2.4842333) q[1];
sx q[1];
rz(-1.4591699) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66508913) q[0];
sx q[0];
rz(-0.71528331) q[0];
sx q[0];
rz(3.1066549) q[0];
rz(-pi) q[1];
rz(0.84648561) q[2];
sx q[2];
rz(-1.5716388) q[2];
sx q[2];
rz(3.0924606) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0399122) q[1];
sx q[1];
rz(-1.3915359) q[1];
sx q[1];
rz(1.0615361) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0351712) q[3];
sx q[3];
rz(-1.0058306) q[3];
sx q[3];
rz(1.21711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3079754) q[2];
sx q[2];
rz(-2.5994382) q[2];
sx q[2];
rz(1.7657492) q[2];
rz(-0.086325072) q[3];
sx q[3];
rz(-0.83266801) q[3];
sx q[3];
rz(-1.8858006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.9432705) q[0];
sx q[0];
rz(-0.69381303) q[0];
sx q[0];
rz(2.2721403) q[0];
rz(2.4122639) q[1];
sx q[1];
rz(-0.84356934) q[1];
sx q[1];
rz(-2.9076911) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37321842) q[0];
sx q[0];
rz(-1.1541751) q[0];
sx q[0];
rz(1.6229902) q[0];
rz(-pi) q[1];
rz(0.42371427) q[2];
sx q[2];
rz(-0.22988453) q[2];
sx q[2];
rz(0.73424852) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.42185509) q[1];
sx q[1];
rz(-2.6736587) q[1];
sx q[1];
rz(-2.809932) q[1];
rz(-pi) q[2];
rz(-1.9332658) q[3];
sx q[3];
rz(-0.7504645) q[3];
sx q[3];
rz(-0.01736162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0539315) q[2];
sx q[2];
rz(-1.4940741) q[2];
sx q[2];
rz(1.9653448) q[2];
rz(-0.67115274) q[3];
sx q[3];
rz(-1.8920369) q[3];
sx q[3];
rz(1.6254856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7072356) q[0];
sx q[0];
rz(-0.86152995) q[0];
sx q[0];
rz(-1.0435411) q[0];
rz(-3.0442944) q[1];
sx q[1];
rz(-1.0568591) q[1];
sx q[1];
rz(1.1204488) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0138884) q[0];
sx q[0];
rz(-2.5426425) q[0];
sx q[0];
rz(-1.4636135) q[0];
rz(-0.78033041) q[2];
sx q[2];
rz(-0.95338168) q[2];
sx q[2];
rz(1.2517901) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.079890117) q[1];
sx q[1];
rz(-1.7482867) q[1];
sx q[1];
rz(-1.6203141) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0886285) q[3];
sx q[3];
rz(-1.4034162) q[3];
sx q[3];
rz(0.24576223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6843159) q[2];
sx q[2];
rz(-2.5280759) q[2];
sx q[2];
rz(1.786001) q[2];
rz(1.5654303) q[3];
sx q[3];
rz(-2.0578945) q[3];
sx q[3];
rz(-1.0160329) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9044357) q[0];
sx q[0];
rz(-0.56023993) q[0];
sx q[0];
rz(2.3518363) q[0];
rz(2.4850028) q[1];
sx q[1];
rz(-2.0273392) q[1];
sx q[1];
rz(-0.56710342) q[1];
rz(-2.1661027) q[2];
sx q[2];
rz(-0.75303034) q[2];
sx q[2];
rz(-1.2122214) q[2];
rz(-1.3764894) q[3];
sx q[3];
rz(-1.8560709) q[3];
sx q[3];
rz(-0.85154497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

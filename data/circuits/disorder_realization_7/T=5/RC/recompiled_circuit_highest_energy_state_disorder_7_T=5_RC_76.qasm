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
rz(0.12662521) q[0];
sx q[0];
rz(-1.5669444) q[0];
sx q[0];
rz(-0.5156762) q[0];
rz(-2.9134143) q[1];
sx q[1];
rz(-2.394634) q[1];
sx q[1];
rz(-0.42626122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3328457) q[0];
sx q[0];
rz(-1.7431743) q[0];
sx q[0];
rz(-0.7953978) q[0];
rz(-pi) q[1];
rz(-0.50349109) q[2];
sx q[2];
rz(-1.3174743) q[2];
sx q[2];
rz(-1.5905141) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8378046) q[1];
sx q[1];
rz(-0.96988867) q[1];
sx q[1];
rz(1.6407042) q[1];
x q[2];
rz(-0.32483806) q[3];
sx q[3];
rz(-1.0110613) q[3];
sx q[3];
rz(-1.7369651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8969741) q[2];
sx q[2];
rz(-1.3099058) q[2];
sx q[2];
rz(-2.8986325) q[2];
rz(2.7729559) q[3];
sx q[3];
rz(-2.536085) q[3];
sx q[3];
rz(-1.2322371) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.826137) q[0];
sx q[0];
rz(-2.9189126) q[0];
sx q[0];
rz(2.7742703) q[0];
rz(-2.3452554) q[1];
sx q[1];
rz(-2.0834736) q[1];
sx q[1];
rz(-2.499089) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135408) q[0];
sx q[0];
rz(-2.0407676) q[0];
sx q[0];
rz(-1.096482) q[0];
x q[1];
rz(-2.5902469) q[2];
sx q[2];
rz(-1.9187821) q[2];
sx q[2];
rz(2.6642193) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5629329) q[1];
sx q[1];
rz(-1.1242492) q[1];
sx q[1];
rz(0.62942735) q[1];
rz(-pi) q[2];
rz(-1.0890555) q[3];
sx q[3];
rz(-1.6828487) q[3];
sx q[3];
rz(-0.58603906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.502304) q[2];
sx q[2];
rz(-2.2059811) q[2];
sx q[2];
rz(1.7712234) q[2];
rz(2.9141407) q[3];
sx q[3];
rz(-1.8919614) q[3];
sx q[3];
rz(1.9500218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74533904) q[0];
sx q[0];
rz(-2.9662913) q[0];
sx q[0];
rz(1.1981717) q[0];
rz(2.0612969) q[1];
sx q[1];
rz(-0.21427576) q[1];
sx q[1];
rz(-3.0467196) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8137275) q[0];
sx q[0];
rz(-1.828232) q[0];
sx q[0];
rz(2.1290921) q[0];
rz(-pi) q[1];
rz(2.7327161) q[2];
sx q[2];
rz(-2.1890702) q[2];
sx q[2];
rz(2.8586819) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2641622) q[1];
sx q[1];
rz(-2.1070711) q[1];
sx q[1];
rz(0.047604851) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1867939) q[3];
sx q[3];
rz(-1.3130762) q[3];
sx q[3];
rz(-2.6741762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0723116) q[2];
sx q[2];
rz(-2.5371964) q[2];
sx q[2];
rz(-2.7351232) q[2];
rz(1.1786002) q[3];
sx q[3];
rz(-1.7013763) q[3];
sx q[3];
rz(1.1627722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.5793295) q[0];
sx q[0];
rz(-3.0070906) q[0];
sx q[0];
rz(0.33600268) q[0];
rz(2.621189) q[1];
sx q[1];
rz(-0.86417472) q[1];
sx q[1];
rz(0.10428183) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9606303) q[0];
sx q[0];
rz(-1.2237566) q[0];
sx q[0];
rz(0.35023679) q[0];
rz(-1.5849131) q[2];
sx q[2];
rz(-1.7272564) q[2];
sx q[2];
rz(0.9597646) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.49890624) q[1];
sx q[1];
rz(-0.79267529) q[1];
sx q[1];
rz(-1.0426056) q[1];
rz(-pi) q[2];
rz(2.9041192) q[3];
sx q[3];
rz(-0.4646315) q[3];
sx q[3];
rz(-0.97685087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5229554) q[2];
sx q[2];
rz(-1.0982265) q[2];
sx q[2];
rz(1.1445507) q[2];
rz(-2.5942904) q[3];
sx q[3];
rz(-1.9132883) q[3];
sx q[3];
rz(-2.0539637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2523786) q[0];
sx q[0];
rz(-0.19110282) q[0];
sx q[0];
rz(-2.5872173) q[0];
rz(-0.91122183) q[1];
sx q[1];
rz(-1.2634042) q[1];
sx q[1];
rz(-2.8401781) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2269991) q[0];
sx q[0];
rz(-2.1560139) q[0];
sx q[0];
rz(0.27497681) q[0];
x q[1];
rz(0.47573702) q[2];
sx q[2];
rz(-2.7103376) q[2];
sx q[2];
rz(-0.83275821) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0423066) q[1];
sx q[1];
rz(-1.9507244) q[1];
sx q[1];
rz(-2.9126588) q[1];
rz(0.88121342) q[3];
sx q[3];
rz(-2.6699319) q[3];
sx q[3];
rz(-0.72615964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.021412795) q[2];
sx q[2];
rz(-0.18099774) q[2];
sx q[2];
rz(-0.0069228355) q[2];
rz(2.210468) q[3];
sx q[3];
rz(-1.1313063) q[3];
sx q[3];
rz(1.1091703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50865006) q[0];
sx q[0];
rz(-0.83704346) q[0];
sx q[0];
rz(1.5337926) q[0];
rz(0.7026698) q[1];
sx q[1];
rz(-0.53121316) q[1];
sx q[1];
rz(0.30219561) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55794278) q[0];
sx q[0];
rz(-1.4690225) q[0];
sx q[0];
rz(2.2324865) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7415761) q[2];
sx q[2];
rz(-0.75655327) q[2];
sx q[2];
rz(0.053002593) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8540878) q[1];
sx q[1];
rz(-2.4937428) q[1];
sx q[1];
rz(-1.8604398) q[1];
x q[2];
rz(1.4190361) q[3];
sx q[3];
rz(-1.5192864) q[3];
sx q[3];
rz(-0.53799483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.89221421) q[2];
sx q[2];
rz(-1.0503146) q[2];
sx q[2];
rz(-2.2155679) q[2];
rz(-1.9269491) q[3];
sx q[3];
rz(-0.49321431) q[3];
sx q[3];
rz(0.27629575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-2.4183384) q[0];
sx q[0];
rz(-2.3746018) q[0];
sx q[0];
rz(-2.0330698) q[0];
rz(-0.33379894) q[1];
sx q[1];
rz(-1.326694) q[1];
sx q[1];
rz(1.0858067) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1340116) q[0];
sx q[0];
rz(-1.3659048) q[0];
sx q[0];
rz(1.6519068) q[0];
rz(1.9721617) q[2];
sx q[2];
rz(-1.737672) q[2];
sx q[2];
rz(-2.3292975) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95744229) q[1];
sx q[1];
rz(-1.0556907) q[1];
sx q[1];
rz(0.95179291) q[1];
rz(0.13768519) q[3];
sx q[3];
rz(-2.2028649) q[3];
sx q[3];
rz(1.1732303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6256025) q[2];
sx q[2];
rz(-2.7770999) q[2];
sx q[2];
rz(-1.1642574) q[2];
rz(0.26816756) q[3];
sx q[3];
rz(-1.2428913) q[3];
sx q[3];
rz(-2.3837762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6044354) q[0];
sx q[0];
rz(-2.6217961) q[0];
sx q[0];
rz(1.9926158) q[0];
rz(-0.2941429) q[1];
sx q[1];
rz(-1.5831169) q[1];
sx q[1];
rz(1.0424967) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9669078) q[0];
sx q[0];
rz(-0.82414851) q[0];
sx q[0];
rz(-2.0281726) q[0];
rz(-pi) q[1];
rz(-2.9127321) q[2];
sx q[2];
rz(-2.3386152) q[2];
sx q[2];
rz(-0.91509089) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5120705) q[1];
sx q[1];
rz(-0.40813706) q[1];
sx q[1];
rz(2.6916608) q[1];
x q[2];
rz(1.6076261) q[3];
sx q[3];
rz(-2.4772518) q[3];
sx q[3];
rz(-0.0078594154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5440172) q[2];
sx q[2];
rz(-2.47051) q[2];
sx q[2];
rz(-1.7768804) q[2];
rz(-2.4890066) q[3];
sx q[3];
rz(-0.79129523) q[3];
sx q[3];
rz(0.64835382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.515601) q[0];
sx q[0];
rz(-0.37726548) q[0];
sx q[0];
rz(-0.01734497) q[0];
rz(3.1248202) q[1];
sx q[1];
rz(-0.78712946) q[1];
sx q[1];
rz(2.9877072) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2543593) q[0];
sx q[0];
rz(-1.6463929) q[0];
sx q[0];
rz(0.31799728) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71435228) q[2];
sx q[2];
rz(-2.0660984) q[2];
sx q[2];
rz(0.53955267) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8887043) q[1];
sx q[1];
rz(-1.9333464) q[1];
sx q[1];
rz(1.2555372) q[1];
rz(-pi) q[2];
x q[2];
rz(0.023970402) q[3];
sx q[3];
rz(-1.4949189) q[3];
sx q[3];
rz(-2.6872203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.953557) q[2];
sx q[2];
rz(-2.3253658) q[2];
sx q[2];
rz(2.6970862) q[2];
rz(0.19017531) q[3];
sx q[3];
rz(-1.5835652) q[3];
sx q[3];
rz(1.3727413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1403777) q[0];
sx q[0];
rz(-1.8920521) q[0];
sx q[0];
rz(1.0706527) q[0];
rz(3.095678) q[1];
sx q[1];
rz(-1.4713947) q[1];
sx q[1];
rz(-2.3505223) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96595461) q[0];
sx q[0];
rz(-1.483206) q[0];
sx q[0];
rz(1.9071294) q[0];
rz(-pi) q[1];
rz(-1.8395109) q[2];
sx q[2];
rz(-1.1860606) q[2];
sx q[2];
rz(0.045298227) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6562499) q[1];
sx q[1];
rz(-2.5963077) q[1];
sx q[1];
rz(-0.7363018) q[1];
rz(-1.3338575) q[3];
sx q[3];
rz(-0.51722368) q[3];
sx q[3];
rz(-2.5631529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4967686) q[2];
sx q[2];
rz(-2.9569929) q[2];
sx q[2];
rz(1.6688639) q[2];
rz(-1.6139) q[3];
sx q[3];
rz(-2.0852641) q[3];
sx q[3];
rz(1.24019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.938217) q[0];
sx q[0];
rz(-1.6435517) q[0];
sx q[0];
rz(2.4284651) q[0];
rz(0.51364246) q[1];
sx q[1];
rz(-1.7084264) q[1];
sx q[1];
rz(1.4485566) q[1];
rz(0.41507872) q[2];
sx q[2];
rz(-2.5998301) q[2];
sx q[2];
rz(2.0143581) q[2];
rz(0.80908262) q[3];
sx q[3];
rz(-1.2590903) q[3];
sx q[3];
rz(2.1303582) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

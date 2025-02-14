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
rz(0.88275498) q[0];
sx q[0];
rz(-2.9967699) q[0];
sx q[0];
rz(0.75135279) q[0];
rz(-1.6169647) q[1];
sx q[1];
rz(-0.96938649) q[1];
sx q[1];
rz(-0.17925395) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0796666) q[0];
sx q[0];
rz(-2.5653337) q[0];
sx q[0];
rz(1.040333) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5380074) q[2];
sx q[2];
rz(-0.84814397) q[2];
sx q[2];
rz(-0.099309534) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4262095) q[1];
sx q[1];
rz(-0.31789624) q[1];
sx q[1];
rz(-2.9098087) q[1];
x q[2];
rz(-3.0823067) q[3];
sx q[3];
rz(-2.2036457) q[3];
sx q[3];
rz(0.92801731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90193191) q[2];
sx q[2];
rz(-1.2942945) q[2];
sx q[2];
rz(0.61061668) q[2];
rz(0.88879746) q[3];
sx q[3];
rz(-2.461268) q[3];
sx q[3];
rz(2.4077267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4438542) q[0];
sx q[0];
rz(-0.30174169) q[0];
sx q[0];
rz(-1.8119716) q[0];
rz(1.4854206) q[1];
sx q[1];
rz(-1.4621567) q[1];
sx q[1];
rz(0.86404538) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1631016) q[0];
sx q[0];
rz(-1.3184883) q[0];
sx q[0];
rz(0.10575328) q[0];
rz(0.99807605) q[2];
sx q[2];
rz(-2.6055899) q[2];
sx q[2];
rz(2.0886476) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.73693528) q[1];
sx q[1];
rz(-0.22316775) q[1];
sx q[1];
rz(0.88921247) q[1];
rz(2.3584189) q[3];
sx q[3];
rz(-1.1044972) q[3];
sx q[3];
rz(-2.9766022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0580505) q[2];
sx q[2];
rz(-2.4594049) q[2];
sx q[2];
rz(-0.19134276) q[2];
rz(-2.8355016) q[3];
sx q[3];
rz(-1.4602665) q[3];
sx q[3];
rz(0.93650854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8323583) q[0];
sx q[0];
rz(-1.2795871) q[0];
sx q[0];
rz(-0.424463) q[0];
rz(1.3407432) q[1];
sx q[1];
rz(-2.2258591) q[1];
sx q[1];
rz(-2.4901966) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.021928) q[0];
sx q[0];
rz(-2.9777973) q[0];
sx q[0];
rz(-1.5271565) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3878787) q[2];
sx q[2];
rz(-1.9486893) q[2];
sx q[2];
rz(1.479508) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5391158) q[1];
sx q[1];
rz(-2.0664363) q[1];
sx q[1];
rz(-1.7784761) q[1];
rz(1.3440787) q[3];
sx q[3];
rz(-1.4170839) q[3];
sx q[3];
rz(2.8938229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1059619) q[2];
sx q[2];
rz(-1.1392081) q[2];
sx q[2];
rz(0.6967217) q[2];
rz(0.68909711) q[3];
sx q[3];
rz(-1.1545811) q[3];
sx q[3];
rz(-0.2564297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43831393) q[0];
sx q[0];
rz(-2.8499481) q[0];
sx q[0];
rz(2.1412204) q[0];
rz(-1.4601624) q[1];
sx q[1];
rz(-1.5337475) q[1];
sx q[1];
rz(1.9283074) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.723169) q[0];
sx q[0];
rz(-2.8349986) q[0];
sx q[0];
rz(-1.1585537) q[0];
rz(-2.0385267) q[2];
sx q[2];
rz(-2.2437895) q[2];
sx q[2];
rz(-1.2601992) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.91649969) q[1];
sx q[1];
rz(-1.2429534) q[1];
sx q[1];
rz(1.4429379) q[1];
rz(-pi) q[2];
rz(1.7217851) q[3];
sx q[3];
rz(-0.54364341) q[3];
sx q[3];
rz(-1.0796987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0393684) q[2];
sx q[2];
rz(-2.3185456) q[2];
sx q[2];
rz(0.064112045) q[2];
rz(-2.4168849) q[3];
sx q[3];
rz(-1.2084081) q[3];
sx q[3];
rz(-0.39976111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0940014) q[0];
sx q[0];
rz(-1.2971017) q[0];
sx q[0];
rz(-2.3498348) q[0];
rz(-2.0296312) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(-1.8066822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4348818) q[0];
sx q[0];
rz(-2.0405053) q[0];
sx q[0];
rz(-0.89969866) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6755988) q[2];
sx q[2];
rz(-2.1490879) q[2];
sx q[2];
rz(2.1190018) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.68685442) q[1];
sx q[1];
rz(-1.4954741) q[1];
sx q[1];
rz(1.7416864) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8800635) q[3];
sx q[3];
rz(-0.59854186) q[3];
sx q[3];
rz(-2.8054597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.84805924) q[2];
sx q[2];
rz(-0.89911014) q[2];
sx q[2];
rz(2.000957) q[2];
rz(-3.0349777) q[3];
sx q[3];
rz(-1.5299608) q[3];
sx q[3];
rz(-2.831736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8150197) q[0];
sx q[0];
rz(-2.6928379) q[0];
sx q[0];
rz(0.003224592) q[0];
rz(0.047317304) q[1];
sx q[1];
rz(-1.686217) q[1];
sx q[1];
rz(1.7914194) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2605476) q[0];
sx q[0];
rz(-0.25940093) q[0];
sx q[0];
rz(1.8015284) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57729849) q[2];
sx q[2];
rz(-2.1566628) q[2];
sx q[2];
rz(2.455204) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1941095) q[1];
sx q[1];
rz(-2.1305741) q[1];
sx q[1];
rz(-0.67621381) q[1];
rz(-1.7058701) q[3];
sx q[3];
rz(-0.22264847) q[3];
sx q[3];
rz(2.6378353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99001592) q[2];
sx q[2];
rz(-0.185597) q[2];
sx q[2];
rz(-0.081175096) q[2];
rz(-0.89933991) q[3];
sx q[3];
rz(-0.81493655) q[3];
sx q[3];
rz(-1.2871294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7154295) q[0];
sx q[0];
rz(-1.3112712) q[0];
sx q[0];
rz(2.6370866) q[0];
rz(2.1465178) q[1];
sx q[1];
rz(-2.1213396) q[1];
sx q[1];
rz(3.1212433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824469) q[0];
sx q[0];
rz(-1.2070281) q[0];
sx q[0];
rz(-0.10031853) q[0];
x q[1];
rz(-2.11436) q[2];
sx q[2];
rz(-2.8994716) q[2];
sx q[2];
rz(1.1309792) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2706302) q[1];
sx q[1];
rz(-2.8125764) q[1];
sx q[1];
rz(-0.94493072) q[1];
rz(-pi) q[2];
rz(2.0836104) q[3];
sx q[3];
rz(-3.0426164) q[3];
sx q[3];
rz(-1.81497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.67123479) q[2];
sx q[2];
rz(-1.3037668) q[2];
sx q[2];
rz(1.7669558) q[2];
rz(-1.6124407) q[3];
sx q[3];
rz(-1.0498472) q[3];
sx q[3];
rz(0.092718743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90877157) q[0];
sx q[0];
rz(-3.0290373) q[0];
sx q[0];
rz(1.0821279) q[0];
rz(2.1807561) q[1];
sx q[1];
rz(-1.6583534) q[1];
sx q[1];
rz(2.198641) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0884224) q[0];
sx q[0];
rz(-0.57500792) q[0];
sx q[0];
rz(-1.0002329) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8846747) q[2];
sx q[2];
rz(-0.65624833) q[2];
sx q[2];
rz(-0.73168025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8170191) q[1];
sx q[1];
rz(-2.1789411) q[1];
sx q[1];
rz(0.84546802) q[1];
rz(-pi) q[2];
rz(-3.037463) q[3];
sx q[3];
rz(-1.2648598) q[3];
sx q[3];
rz(2.1776074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1711787) q[2];
sx q[2];
rz(-1.8381939) q[2];
sx q[2];
rz(1.9937493) q[2];
rz(-2.7796699) q[3];
sx q[3];
rz(-1.3779093) q[3];
sx q[3];
rz(1.3981147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73087937) q[0];
sx q[0];
rz(-0.35922265) q[0];
sx q[0];
rz(0.10243375) q[0];
rz(0.61406413) q[1];
sx q[1];
rz(-1.026261) q[1];
sx q[1];
rz(-0.817743) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21120223) q[0];
sx q[0];
rz(-1.0387392) q[0];
sx q[0];
rz(-2.3389057) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7162343) q[2];
sx q[2];
rz(-2.465473) q[2];
sx q[2];
rz(-0.506625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9407258) q[1];
sx q[1];
rz(-1.5897017) q[1];
sx q[1];
rz(0.84152842) q[1];
rz(-pi) q[2];
rz(2.1675112) q[3];
sx q[3];
rz(-2.4170205) q[3];
sx q[3];
rz(2.4251079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7927336) q[2];
sx q[2];
rz(-1.0772971) q[2];
sx q[2];
rz(-0.31361541) q[2];
rz(2.6775728) q[3];
sx q[3];
rz(-2.2950164) q[3];
sx q[3];
rz(-0.919842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2713276) q[0];
sx q[0];
rz(-2.3509404) q[0];
sx q[0];
rz(2.9050997) q[0];
rz(0.21367167) q[1];
sx q[1];
rz(-0.90286076) q[1];
sx q[1];
rz(2.9170091) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.41193) q[0];
sx q[0];
rz(-0.87422919) q[0];
sx q[0];
rz(-0.64944546) q[0];
x q[1];
rz(-2.6170391) q[2];
sx q[2];
rz(-0.65215014) q[2];
sx q[2];
rz(-0.99936501) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1285291) q[1];
sx q[1];
rz(-2.6010102) q[1];
sx q[1];
rz(-0.22352) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8524695) q[3];
sx q[3];
rz(-1.2465047) q[3];
sx q[3];
rz(-1.7773903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1647722) q[2];
sx q[2];
rz(-2.2701264) q[2];
sx q[2];
rz(2.9898047) q[2];
rz(3.1101036) q[3];
sx q[3];
rz(-1.7446691) q[3];
sx q[3];
rz(-0.12846863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0781773) q[0];
sx q[0];
rz(-1.5237533) q[0];
sx q[0];
rz(-0.40405003) q[0];
rz(-2.0582485) q[1];
sx q[1];
rz(-1.5768408) q[1];
sx q[1];
rz(1.5595938) q[1];
rz(-2.940098) q[2];
sx q[2];
rz(-1.6215743) q[2];
sx q[2];
rz(-1.8334186) q[2];
rz(2.8121171) q[3];
sx q[3];
rz(-1.4737045) q[3];
sx q[3];
rz(1.2274689) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

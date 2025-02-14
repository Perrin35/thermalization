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
rz(-1.3731319) q[0];
sx q[0];
rz(-1.4302349) q[0];
sx q[0];
rz(-0.89682427) q[0];
rz(-1.7984017) q[1];
sx q[1];
rz(3.922037) q[1];
sx q[1];
rz(8.2889397) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97082393) q[0];
sx q[0];
rz(-1.6410195) q[0];
sx q[0];
rz(-1.4115566) q[0];
rz(-pi) q[1];
rz(-0.48223038) q[2];
sx q[2];
rz(-2.1948174) q[2];
sx q[2];
rz(0.66971794) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.070277409) q[1];
sx q[1];
rz(-1.0272553) q[1];
sx q[1];
rz(-0.75884968) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1862983) q[3];
sx q[3];
rz(-1.1685598) q[3];
sx q[3];
rz(-1.1768794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4844126) q[2];
sx q[2];
rz(-1.5766532) q[2];
sx q[2];
rz(-2.8489992) q[2];
rz(-1.6491133) q[3];
sx q[3];
rz(-1.7238659) q[3];
sx q[3];
rz(2.9610146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-3.0507616) q[0];
sx q[0];
rz(-1.917559) q[0];
sx q[0];
rz(-0.60371387) q[0];
rz(-1.0366038) q[1];
sx q[1];
rz(-0.34244582) q[1];
sx q[1];
rz(0.93050122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11545746) q[0];
sx q[0];
rz(-1.8107426) q[0];
sx q[0];
rz(-1.351993) q[0];
x q[1];
rz(2.5803135) q[2];
sx q[2];
rz(-1.6973472) q[2];
sx q[2];
rz(1.2685219) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.249037) q[1];
sx q[1];
rz(-2.3904789) q[1];
sx q[1];
rz(2.9455801) q[1];
x q[2];
rz(-2.811238) q[3];
sx q[3];
rz(-1.0022911) q[3];
sx q[3];
rz(-1.1363883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7751864) q[2];
sx q[2];
rz(-1.8975056) q[2];
sx q[2];
rz(-1.5001971) q[2];
rz(-0.90534798) q[3];
sx q[3];
rz(-1.8843745) q[3];
sx q[3];
rz(-2.1998028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4582849) q[0];
sx q[0];
rz(-1.6907121) q[0];
sx q[0];
rz(-2.2990551) q[0];
rz(-1.9375577) q[1];
sx q[1];
rz(-1.8572109) q[1];
sx q[1];
rz(-1.1010928) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0718449) q[0];
sx q[0];
rz(-1.4948973) q[0];
sx q[0];
rz(-1.2420446) q[0];
x q[1];
rz(1.3705047) q[2];
sx q[2];
rz(-2.3034987) q[2];
sx q[2];
rz(1.0685234) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8070911) q[1];
sx q[1];
rz(-1.1192296) q[1];
sx q[1];
rz(-1.0784763) q[1];
rz(2.1870062) q[3];
sx q[3];
rz(-1.2266739) q[3];
sx q[3];
rz(1.6801429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4568242) q[2];
sx q[2];
rz(-2.5825239) q[2];
sx q[2];
rz(-0.10747257) q[2];
rz(2.0939854) q[3];
sx q[3];
rz(-1.3370297) q[3];
sx q[3];
rz(-1.6155155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2408214) q[0];
sx q[0];
rz(-2.695684) q[0];
sx q[0];
rz(-1.4500424) q[0];
rz(-0.58174539) q[1];
sx q[1];
rz(-0.83256871) q[1];
sx q[1];
rz(-1.1004826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9387526) q[0];
sx q[0];
rz(-1.0577518) q[0];
sx q[0];
rz(-2.0111397) q[0];
rz(-pi) q[1];
rz(1.9504488) q[2];
sx q[2];
rz(-1.0565524) q[2];
sx q[2];
rz(0.66858993) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2128513) q[1];
sx q[1];
rz(-2.263816) q[1];
sx q[1];
rz(-1.5726967) q[1];
x q[2];
rz(-0.29275818) q[3];
sx q[3];
rz(-2.2653711) q[3];
sx q[3];
rz(3.0931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.45740286) q[2];
sx q[2];
rz(-2.0997014) q[2];
sx q[2];
rz(-2.5390967) q[2];
rz(-1.4495133) q[3];
sx q[3];
rz(-1.4256198) q[3];
sx q[3];
rz(2.1879533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36798894) q[0];
sx q[0];
rz(-1.8743176) q[0];
sx q[0];
rz(1.3195272) q[0];
rz(-0.52917448) q[1];
sx q[1];
rz(-1.8942984) q[1];
sx q[1];
rz(-0.41437638) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0396578) q[0];
sx q[0];
rz(-0.62292719) q[0];
sx q[0];
rz(-0.85467546) q[0];
rz(2.2362382) q[2];
sx q[2];
rz(-0.86595172) q[2];
sx q[2];
rz(-3.1107748) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9739009) q[1];
sx q[1];
rz(-1.138559) q[1];
sx q[1];
rz(0.28829379) q[1];
x q[2];
rz(-0.32043362) q[3];
sx q[3];
rz(-1.7263573) q[3];
sx q[3];
rz(-1.8309586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6759912) q[2];
sx q[2];
rz(-0.71054116) q[2];
sx q[2];
rz(0.65593925) q[2];
rz(-1.2300308) q[3];
sx q[3];
rz(-1.1070808) q[3];
sx q[3];
rz(-0.51903498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32992724) q[0];
sx q[0];
rz(-2.1274607) q[0];
sx q[0];
rz(-2.2299679) q[0];
rz(1.3765593) q[1];
sx q[1];
rz(-2.7477317) q[1];
sx q[1];
rz(-0.54042655) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86397643) q[0];
sx q[0];
rz(-1.8555529) q[0];
sx q[0];
rz(-0.44313125) q[0];
rz(-pi) q[1];
rz(-3.0824605) q[2];
sx q[2];
rz(-0.34854564) q[2];
sx q[2];
rz(2.4273171) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0192761) q[1];
sx q[1];
rz(-1.5088006) q[1];
sx q[1];
rz(-2.7424395) q[1];
x q[2];
rz(-1.9283617) q[3];
sx q[3];
rz(-0.66047664) q[3];
sx q[3];
rz(-1.9102526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8526758) q[2];
sx q[2];
rz(-1.4492852) q[2];
sx q[2];
rz(2.0055611) q[2];
rz(1.3352669) q[3];
sx q[3];
rz(-1.5583594) q[3];
sx q[3];
rz(-3.0291338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54953185) q[0];
sx q[0];
rz(-0.15287748) q[0];
sx q[0];
rz(-0.98912799) q[0];
rz(-1.2024744) q[1];
sx q[1];
rz(-1.6705284) q[1];
sx q[1];
rz(2.8340526) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0796442) q[0];
sx q[0];
rz(-2.7833118) q[0];
sx q[0];
rz(-0.69967358) q[0];
x q[1];
rz(0.73300377) q[2];
sx q[2];
rz(-1.8585732) q[2];
sx q[2];
rz(2.6633898) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7293528) q[1];
sx q[1];
rz(-1.5478999) q[1];
sx q[1];
rz(2.9213133) q[1];
x q[2];
rz(0.6102517) q[3];
sx q[3];
rz(-2.0008127) q[3];
sx q[3];
rz(-0.56320923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5146553) q[2];
sx q[2];
rz(-2.3263826) q[2];
sx q[2];
rz(1.6132272) q[2];
rz(-2.3441687) q[3];
sx q[3];
rz(-1.593109) q[3];
sx q[3];
rz(3.048866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02354694) q[0];
sx q[0];
rz(-0.67414701) q[0];
sx q[0];
rz(2.4349037) q[0];
rz(-2.8382909) q[1];
sx q[1];
rz(-1.3677596) q[1];
sx q[1];
rz(-1.2695262) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2523018) q[0];
sx q[0];
rz(-1.902475) q[0];
sx q[0];
rz(-1.7097019) q[0];
rz(-pi) q[1];
rz(0.49497382) q[2];
sx q[2];
rz(-1.7893063) q[2];
sx q[2];
rz(-1.72118) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4321338) q[1];
sx q[1];
rz(-1.7466326) q[1];
sx q[1];
rz(-1.9587112) q[1];
rz(-pi) q[2];
rz(-0.16002051) q[3];
sx q[3];
rz(-0.56462641) q[3];
sx q[3];
rz(-0.52288632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.068744008) q[2];
sx q[2];
rz(-0.35126433) q[2];
sx q[2];
rz(2.6212485) q[2];
rz(1.2460234) q[3];
sx q[3];
rz(-1.4166074) q[3];
sx q[3];
rz(2.0791159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50616566) q[0];
sx q[0];
rz(-1.4816254) q[0];
sx q[0];
rz(2.6284499) q[0];
rz(0.48140934) q[1];
sx q[1];
rz(-0.82999271) q[1];
sx q[1];
rz(-0.76024461) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8943065) q[0];
sx q[0];
rz(-1.5991976) q[0];
sx q[0];
rz(3.1288163) q[0];
x q[1];
rz(-0.77519007) q[2];
sx q[2];
rz(-1.7235867) q[2];
sx q[2];
rz(0.82681954) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7920391) q[1];
sx q[1];
rz(-2.4124618) q[1];
sx q[1];
rz(-0.11138536) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3586274) q[3];
sx q[3];
rz(-0.43624076) q[3];
sx q[3];
rz(1.0324163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5236139) q[2];
sx q[2];
rz(-0.60680497) q[2];
sx q[2];
rz(2.4904909) q[2];
rz(1.9304322) q[3];
sx q[3];
rz(-1.9821207) q[3];
sx q[3];
rz(-2.7291362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8170269) q[0];
sx q[0];
rz(-0.71837076) q[0];
sx q[0];
rz(0.73085648) q[0];
rz(0.51365799) q[1];
sx q[1];
rz(-2.3274603) q[1];
sx q[1];
rz(-0.51630744) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8532787) q[0];
sx q[0];
rz(-0.95398879) q[0];
sx q[0];
rz(0.46648394) q[0];
x q[1];
rz(1.1942719) q[2];
sx q[2];
rz(-1.9360844) q[2];
sx q[2];
rz(0.69460642) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4708206) q[1];
sx q[1];
rz(-1.1750452) q[1];
sx q[1];
rz(2.938063) q[1];
x q[2];
rz(0.29521355) q[3];
sx q[3];
rz(-1.3564979) q[3];
sx q[3];
rz(-0.47116531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.98661304) q[2];
sx q[2];
rz(-1.2519138) q[2];
sx q[2];
rz(0.66407859) q[2];
rz(-0.99728552) q[3];
sx q[3];
rz(-1.4279131) q[3];
sx q[3];
rz(-0.47521457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3288788) q[0];
sx q[0];
rz(-1.7408149) q[0];
sx q[0];
rz(-0.029205532) q[0];
rz(3.008814) q[1];
sx q[1];
rz(-1.7280424) q[1];
sx q[1];
rz(-0.8676563) q[1];
rz(2.6399707) q[2];
sx q[2];
rz(-1.8799964) q[2];
sx q[2];
rz(2.7867405) q[2];
rz(2.7326591) q[3];
sx q[3];
rz(-1.144617) q[3];
sx q[3];
rz(-1.495818) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

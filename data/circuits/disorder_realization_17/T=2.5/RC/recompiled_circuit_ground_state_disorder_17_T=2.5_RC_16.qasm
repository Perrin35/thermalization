OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.75093961) q[0];
sx q[0];
rz(-2.5680464) q[0];
sx q[0];
rz(0.030800495) q[0];
rz(6.6703086) q[1];
sx q[1];
rz(3.3422254) q[1];
sx q[1];
rz(2.2320342) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46900374) q[0];
sx q[0];
rz(-1.0116972) q[0];
sx q[0];
rz(-0.78796537) q[0];
x q[1];
rz(2.4417905) q[2];
sx q[2];
rz(-1.095311) q[2];
sx q[2];
rz(-0.83719992) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2647029) q[1];
sx q[1];
rz(-0.24395058) q[1];
sx q[1];
rz(0.96744858) q[1];
rz(-pi) q[2];
rz(-1.2991139) q[3];
sx q[3];
rz(-1.6713023) q[3];
sx q[3];
rz(-2.1372014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6261141) q[2];
sx q[2];
rz(-1.6823744) q[2];
sx q[2];
rz(0.98575753) q[2];
rz(1.9043026) q[3];
sx q[3];
rz(-0.83043778) q[3];
sx q[3];
rz(-1.5498836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80376959) q[0];
sx q[0];
rz(-1.060312) q[0];
sx q[0];
rz(0.34326237) q[0];
rz(2.901851) q[1];
sx q[1];
rz(-2.7570351) q[1];
sx q[1];
rz(-0.029031001) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0800901) q[0];
sx q[0];
rz(-2.6250589) q[0];
sx q[0];
rz(-0.37688984) q[0];
x q[1];
rz(1.5324872) q[2];
sx q[2];
rz(-0.84332436) q[2];
sx q[2];
rz(-3.0366355) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1085384) q[1];
sx q[1];
rz(-1.5863451) q[1];
sx q[1];
rz(-0.59301807) q[1];
x q[2];
rz(0.47201158) q[3];
sx q[3];
rz(-2.4084512) q[3];
sx q[3];
rz(2.5705119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.79641882) q[2];
sx q[2];
rz(-2.5721305) q[2];
sx q[2];
rz(0.84863895) q[2];
rz(-0.38550115) q[3];
sx q[3];
rz(-1.4717088) q[3];
sx q[3];
rz(2.8482385) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.327453) q[0];
sx q[0];
rz(-1.7610981) q[0];
sx q[0];
rz(3.0019548) q[0];
rz(-0.26593727) q[1];
sx q[1];
rz(-1.9640924) q[1];
sx q[1];
rz(2.1279443) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.06482) q[0];
sx q[0];
rz(-2.1401323) q[0];
sx q[0];
rz(2.6808617) q[0];
rz(-0.43439867) q[2];
sx q[2];
rz(-1.1936099) q[2];
sx q[2];
rz(2.3247256) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7820236) q[1];
sx q[1];
rz(-1.4064565) q[1];
sx q[1];
rz(-0.4811299) q[1];
x q[2];
rz(-3.072282) q[3];
sx q[3];
rz(-1.4678174) q[3];
sx q[3];
rz(-1.6268831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2024112) q[2];
sx q[2];
rz(-0.95042578) q[2];
sx q[2];
rz(-2.6964296) q[2];
rz(-2.4116481) q[3];
sx q[3];
rz(-1.688262) q[3];
sx q[3];
rz(-1.2030407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79778033) q[0];
sx q[0];
rz(-0.36425632) q[0];
sx q[0];
rz(1.9601747) q[0];
rz(-1.5161318) q[1];
sx q[1];
rz(-1.5189891) q[1];
sx q[1];
rz(-0.48798645) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054783527) q[0];
sx q[0];
rz(-1.8300608) q[0];
sx q[0];
rz(-1.7340585) q[0];
rz(0.96089604) q[2];
sx q[2];
rz(-0.95003613) q[2];
sx q[2];
rz(-1.4767978) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1688465) q[1];
sx q[1];
rz(-2.1543145) q[1];
sx q[1];
rz(0.31202392) q[1];
x q[2];
rz(2.7236688) q[3];
sx q[3];
rz(-1.2439787) q[3];
sx q[3];
rz(-0.55909294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2130012) q[2];
sx q[2];
rz(-1.4139621) q[2];
sx q[2];
rz(-2.6334527) q[2];
rz(2.7022434) q[3];
sx q[3];
rz(-2.5033247) q[3];
sx q[3];
rz(-0.58258575) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64684922) q[0];
sx q[0];
rz(-1.2127533) q[0];
sx q[0];
rz(0.61432046) q[0];
rz(2.1454504) q[1];
sx q[1];
rz(-0.83729815) q[1];
sx q[1];
rz(-0.070929758) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9324294) q[0];
sx q[0];
rz(-2.4485561) q[0];
sx q[0];
rz(-2.0636051) q[0];
x q[1];
rz(1.1267012) q[2];
sx q[2];
rz(-1.3108062) q[2];
sx q[2];
rz(0.6636493) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.48551895) q[1];
sx q[1];
rz(-0.89886256) q[1];
sx q[1];
rz(-2.333332) q[1];
rz(-pi) q[2];
rz(-2.4627389) q[3];
sx q[3];
rz(-1.0104138) q[3];
sx q[3];
rz(2.7075775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7190242) q[2];
sx q[2];
rz(-1.7302128) q[2];
sx q[2];
rz(-3.0651429) q[2];
rz(-3.0494173) q[3];
sx q[3];
rz(-1.1276827) q[3];
sx q[3];
rz(2.5341212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3344264) q[0];
sx q[0];
rz(-2.9521515) q[0];
sx q[0];
rz(-0.46519753) q[0];
rz(-3.0472962) q[1];
sx q[1];
rz(-1.0153208) q[1];
sx q[1];
rz(1.2634855) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7955728) q[0];
sx q[0];
rz(-1.4515522) q[0];
sx q[0];
rz(-1.9247965) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0025778254) q[2];
sx q[2];
rz(-0.63226223) q[2];
sx q[2];
rz(-0.31429502) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.81663075) q[1];
sx q[1];
rz(-1.1405319) q[1];
sx q[1];
rz(2.5041786) q[1];
rz(-1.1912548) q[3];
sx q[3];
rz(-2.0058245) q[3];
sx q[3];
rz(-1.7779999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5491526) q[2];
sx q[2];
rz(-2.0024039) q[2];
sx q[2];
rz(0.69063866) q[2];
rz(-1.45951) q[3];
sx q[3];
rz(-0.032460902) q[3];
sx q[3];
rz(-0.261511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.14677793) q[0];
sx q[0];
rz(-1.0228461) q[0];
sx q[0];
rz(2.6652375) q[0];
rz(-1.1789383) q[1];
sx q[1];
rz(-2.4504688) q[1];
sx q[1];
rz(-1.8144511) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9479312) q[0];
sx q[0];
rz(-1.9597988) q[0];
sx q[0];
rz(2.6575628) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.671538) q[2];
sx q[2];
rz(-2.6664554) q[2];
sx q[2];
rz(2.6784999) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9500014) q[1];
sx q[1];
rz(-2.0954544) q[1];
sx q[1];
rz(1.4895951) q[1];
rz(2.5267739) q[3];
sx q[3];
rz(-1.58883) q[3];
sx q[3];
rz(-0.65476894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2771161) q[2];
sx q[2];
rz(-1.4863374) q[2];
sx q[2];
rz(-0.76249301) q[2];
rz(-1.3660733) q[3];
sx q[3];
rz(-1.6906747) q[3];
sx q[3];
rz(-0.33179992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45258006) q[0];
sx q[0];
rz(-0.068050139) q[0];
sx q[0];
rz(-2.3633603) q[0];
rz(-1.4100086) q[1];
sx q[1];
rz(-0.86763132) q[1];
sx q[1];
rz(0.57077879) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3616791) q[0];
sx q[0];
rz(-0.81244367) q[0];
sx q[0];
rz(-2.1606584) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0887509) q[2];
sx q[2];
rz(-1.2882917) q[2];
sx q[2];
rz(-1.5517915) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0062756) q[1];
sx q[1];
rz(-1.0500188) q[1];
sx q[1];
rz(0.09002491) q[1];
x q[2];
rz(-0.024939288) q[3];
sx q[3];
rz(-0.93526269) q[3];
sx q[3];
rz(-0.81944114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1125696) q[2];
sx q[2];
rz(-1.7971635) q[2];
sx q[2];
rz(-2.6018108) q[2];
rz(0.44627407) q[3];
sx q[3];
rz(-0.73791426) q[3];
sx q[3];
rz(2.3257183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0453294) q[0];
sx q[0];
rz(-2.6427866) q[0];
sx q[0];
rz(-2.2742284) q[0];
rz(1.7034886) q[1];
sx q[1];
rz(-0.73355621) q[1];
sx q[1];
rz(0.47384706) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6561013) q[0];
sx q[0];
rz(-0.37603077) q[0];
sx q[0];
rz(-2.4185527) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7452186) q[2];
sx q[2];
rz(-2.1389942) q[2];
sx q[2];
rz(-0.098357226) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.074939097) q[1];
sx q[1];
rz(-0.84006305) q[1];
sx q[1];
rz(-1.380062) q[1];
rz(-0.97029347) q[3];
sx q[3];
rz(-2.3564995) q[3];
sx q[3];
rz(-0.11150211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69502407) q[2];
sx q[2];
rz(-2.5769672) q[2];
sx q[2];
rz(1.1308283) q[2];
rz(1.1787777) q[3];
sx q[3];
rz(-1.7480353) q[3];
sx q[3];
rz(0.47275561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9216565) q[0];
sx q[0];
rz(-0.49880767) q[0];
sx q[0];
rz(2.1529799) q[0];
rz(1.0722718) q[1];
sx q[1];
rz(-1.5360473) q[1];
sx q[1];
rz(-1.4354717) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6349115) q[0];
sx q[0];
rz(-1.5456063) q[0];
sx q[0];
rz(-2.1468625) q[0];
rz(-3.1085787) q[2];
sx q[2];
rz(-0.22840127) q[2];
sx q[2];
rz(1.8637125) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2769673) q[1];
sx q[1];
rz(-1.5706148) q[1];
sx q[1];
rz(-1.6305883) q[1];
rz(-pi) q[2];
rz(1.9063453) q[3];
sx q[3];
rz(-2.23508) q[3];
sx q[3];
rz(1.2227525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2995305) q[2];
sx q[2];
rz(-2.5008423) q[2];
sx q[2];
rz(0.93471849) q[2];
rz(1.0457906) q[3];
sx q[3];
rz(-0.17242923) q[3];
sx q[3];
rz(-0.24833965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1187779) q[0];
sx q[0];
rz(-1.3640484) q[0];
sx q[0];
rz(1.7400297) q[0];
rz(3.0339495) q[1];
sx q[1];
rz(-1.5168774) q[1];
sx q[1];
rz(0.37227896) q[1];
rz(-1.8423745) q[2];
sx q[2];
rz(-0.80119802) q[2];
sx q[2];
rz(-0.3225633) q[2];
rz(-0.94318642) q[3];
sx q[3];
rz(-1.9780157) q[3];
sx q[3];
rz(-2.1383022) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

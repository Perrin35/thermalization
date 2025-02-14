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
rz(1.1531416) q[0];
sx q[0];
rz(-0.81557953) q[0];
sx q[0];
rz(-0.75818914) q[0];
rz(-0.012501333) q[1];
sx q[1];
rz(-1.3036417) q[1];
sx q[1];
rz(1.5703896) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9725295) q[0];
sx q[0];
rz(-0.95989812) q[0];
sx q[0];
rz(-1.1283895) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.067367359) q[2];
sx q[2];
rz(-1.172003) q[2];
sx q[2];
rz(0.68902868) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.058152288) q[1];
sx q[1];
rz(-1.5136763) q[1];
sx q[1];
rz(1.9316462) q[1];
rz(-pi) q[2];
rz(1.8800354) q[3];
sx q[3];
rz(-1.5133563) q[3];
sx q[3];
rz(-0.32306898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5300753) q[2];
sx q[2];
rz(-3.1334183) q[2];
sx q[2];
rz(0.44069904) q[2];
rz(-0.079744451) q[3];
sx q[3];
rz(-0.00010448797) q[3];
sx q[3];
rz(-2.0174446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.9663064) q[0];
sx q[0];
rz(-0.056862406) q[0];
sx q[0];
rz(2.9735907) q[0];
rz(3.1213144) q[1];
sx q[1];
rz(-2.8333277) q[1];
sx q[1];
rz(1.6049989) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56242053) q[0];
sx q[0];
rz(-1.354166) q[0];
sx q[0];
rz(-0.24026339) q[0];
rz(-pi) q[1];
rz(-3.1390983) q[2];
sx q[2];
rz(-1.5809142) q[2];
sx q[2];
rz(-1.5960787) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78129301) q[1];
sx q[1];
rz(-1.5661245) q[1];
sx q[1];
rz(-0.0010990573) q[1];
rz(3.0992989) q[3];
sx q[3];
rz(-1.5797857) q[3];
sx q[3];
rz(-0.90153722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7961879) q[2];
sx q[2];
rz(-2.2296843) q[2];
sx q[2];
rz(-1.7339285) q[2];
rz(1.0450854) q[3];
sx q[3];
rz(-3.0920691) q[3];
sx q[3];
rz(-0.27200562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.8318091) q[0];
sx q[0];
rz(-0.97667664) q[0];
sx q[0];
rz(0.56104863) q[0];
rz(-0.27684119) q[1];
sx q[1];
rz(-0.012877348) q[1];
sx q[1];
rz(-1.8337839) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1501859) q[0];
sx q[0];
rz(-1.6117818) q[0];
sx q[0];
rz(1.2738373) q[0];
rz(-pi) q[1];
rz(-1.5707914) q[2];
sx q[2];
rz(-1.563407) q[2];
sx q[2];
rz(1.0634729) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0289291) q[1];
sx q[1];
rz(-1.6293007) q[1];
sx q[1];
rz(-2.1444291) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4588608) q[3];
sx q[3];
rz(-0.94597497) q[3];
sx q[3];
rz(-0.030908728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7967367) q[2];
sx q[2];
rz(-3.1414746) q[2];
sx q[2];
rz(-2.5522088) q[2];
rz(1.2423337) q[3];
sx q[3];
rz(-0.012367736) q[3];
sx q[3];
rz(-1.7813659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0068483343) q[0];
sx q[0];
rz(-0.51181781) q[0];
sx q[0];
rz(-1.7929329) q[0];
rz(-0.0066561247) q[1];
sx q[1];
rz(-1.321188) q[1];
sx q[1];
rz(-0.033500813) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2804256) q[0];
sx q[0];
rz(-2.4651732) q[0];
sx q[0];
rz(1.8489714) q[0];
x q[1];
rz(3.0219565) q[2];
sx q[2];
rz(-1.575261) q[2];
sx q[2];
rz(-1.9890832) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0211612) q[1];
sx q[1];
rz(-2.8725) q[1];
sx q[1];
rz(-1.5325559) q[1];
rz(-pi) q[2];
rz(1.4064404) q[3];
sx q[3];
rz(-1.6911611) q[3];
sx q[3];
rz(-0.09935483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.032430705) q[2];
sx q[2];
rz(-0.0065294821) q[2];
sx q[2];
rz(2.8429441) q[2];
rz(1.8715035) q[3];
sx q[3];
rz(-3.1254369) q[3];
sx q[3];
rz(0.0531918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.3700767) q[0];
sx q[0];
rz(-1.6271485) q[0];
sx q[0];
rz(0.44684967) q[0];
rz(-2.9554548) q[1];
sx q[1];
rz(-0.061807241) q[1];
sx q[1];
rz(-1.4131379) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22767775) q[0];
sx q[0];
rz(-0.23328885) q[0];
sx q[0];
rz(1.7167709) q[0];
rz(1.9930219) q[2];
sx q[2];
rz(-1.3729501) q[2];
sx q[2];
rz(-2.5289218) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.062033) q[1];
sx q[1];
rz(-1.6158982) q[1];
sx q[1];
rz(1.6203887) q[1];
rz(-pi) q[2];
rz(2.917971) q[3];
sx q[3];
rz(-1.2960172) q[3];
sx q[3];
rz(-1.4449121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3073005) q[2];
sx q[2];
rz(-1.5520381) q[2];
sx q[2];
rz(-0.49868047) q[2];
rz(-0.57791609) q[3];
sx q[3];
rz(-2.6579865) q[3];
sx q[3];
rz(-0.59670603) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1805098) q[0];
sx q[0];
rz(-1.120765) q[0];
sx q[0];
rz(-0.34641308) q[0];
rz(0.60180426) q[1];
sx q[1];
rz(-1.5806942) q[1];
sx q[1];
rz(-0.7535038) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44418535) q[0];
sx q[0];
rz(-1.3988528) q[0];
sx q[0];
rz(-2.9471022) q[0];
rz(-pi) q[1];
rz(2.1624203) q[2];
sx q[2];
rz(-0.1340999) q[2];
sx q[2];
rz(-0.025452415) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.063976377) q[1];
sx q[1];
rz(-1.9850296) q[1];
sx q[1];
rz(-0.86552455) q[1];
x q[2];
rz(-1.0585045) q[3];
sx q[3];
rz(-2.9550094) q[3];
sx q[3];
rz(0.032840289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.57009131) q[2];
sx q[2];
rz(-3.1381021) q[2];
sx q[2];
rz(1.6174779) q[2];
rz(-3.0213455) q[3];
sx q[3];
rz(-3.1382939) q[3];
sx q[3];
rz(0.53774589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26871249) q[0];
sx q[0];
rz(-2.217642) q[0];
sx q[0];
rz(-0.22126108) q[0];
rz(-1.6809173) q[1];
sx q[1];
rz(-0.9333846) q[1];
sx q[1];
rz(-0.077300765) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0635522) q[0];
sx q[0];
rz(-1.5693328) q[0];
sx q[0];
rz(1.5290676) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5796698) q[2];
sx q[2];
rz(-1.5755782) q[2];
sx q[2];
rz(1.3591131) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.12880023) q[1];
sx q[1];
rz(-0.18928738) q[1];
sx q[1];
rz(2.7526593) q[1];
x q[2];
rz(-1.7847925) q[3];
sx q[3];
rz(-1.6955175) q[3];
sx q[3];
rz(0.41984841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7920502) q[2];
sx q[2];
rz(-0.011186102) q[2];
sx q[2];
rz(0.95996094) q[2];
rz(-0.32944426) q[3];
sx q[3];
rz(-0.0080778413) q[3];
sx q[3];
rz(-2.2913057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(2.195381) q[0];
sx q[0];
rz(-0.61310261) q[0];
sx q[0];
rz(0.10928133) q[0];
rz(-2.7682313) q[1];
sx q[1];
rz(-0.80972087) q[1];
sx q[1];
rz(-1.2304617) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62723535) q[0];
sx q[0];
rz(-2.1684596) q[0];
sx q[0];
rz(2.3742832) q[0];
rz(-pi) q[1];
rz(2.3745499) q[2];
sx q[2];
rz(-0.27077507) q[2];
sx q[2];
rz(-2.3363638) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8671678) q[1];
sx q[1];
rz(-1.6376357) q[1];
sx q[1];
rz(-0.011456077) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91992141) q[3];
sx q[3];
rz(-1.3215142) q[3];
sx q[3];
rz(1.1703619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5750778) q[2];
sx q[2];
rz(-1.2357624) q[2];
sx q[2];
rz(-1.3285948) q[2];
rz(1.396842) q[3];
sx q[3];
rz(-0.003740398) q[3];
sx q[3];
rz(-1.0218792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063865572) q[0];
sx q[0];
rz(-1.4430178) q[0];
sx q[0];
rz(-0.57300895) q[0];
rz(-0.30814463) q[1];
sx q[1];
rz(-0.40987086) q[1];
sx q[1];
rz(-2.134197) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44237374) q[0];
sx q[0];
rz(-1.464252) q[0];
sx q[0];
rz(-3.1359948) q[0];
x q[1];
rz(-1.4332268) q[2];
sx q[2];
rz(-1.614288) q[2];
sx q[2];
rz(1.5303591) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7511661) q[1];
sx q[1];
rz(-1.6114283) q[1];
sx q[1];
rz(1.4875814) q[1];
rz(1.9898207) q[3];
sx q[3];
rz(-0.033009987) q[3];
sx q[3];
rz(0.21402436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.31743) q[2];
sx q[2];
rz(-2.514826) q[2];
sx q[2];
rz(0.38995788) q[2];
rz(3.0725078) q[3];
sx q[3];
rz(-0.0091113541) q[3];
sx q[3];
rz(0.33153427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020141715) q[0];
sx q[0];
rz(-2.3905601) q[0];
sx q[0];
rz(2.6556515) q[0];
rz(0.87156975) q[1];
sx q[1];
rz(-1.3078682) q[1];
sx q[1];
rz(1.647324) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93372351) q[0];
sx q[0];
rz(-1.2656801) q[0];
sx q[0];
rz(0.56201571) q[0];
x q[1];
rz(-2.5267692) q[2];
sx q[2];
rz(-1.5346926) q[2];
sx q[2];
rz(1.5272128) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4758265) q[1];
sx q[1];
rz(-1.8987055) q[1];
sx q[1];
rz(1.2535415) q[1];
rz(-pi) q[2];
rz(-0.11313862) q[3];
sx q[3];
rz(-1.5290878) q[3];
sx q[3];
rz(2.0737518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5668874) q[2];
sx q[2];
rz(-0.042782728) q[2];
sx q[2];
rz(3.1097143) q[2];
rz(0.77830642) q[3];
sx q[3];
rz(-0.0068155546) q[3];
sx q[3];
rz(-0.2955029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.42302172) q[0];
sx q[0];
rz(-1.6091249) q[0];
sx q[0];
rz(-1.3269497) q[0];
rz(3.0146535) q[1];
sx q[1];
rz(-2.9025684) q[1];
sx q[1];
rz(-2.92166) q[1];
rz(-1.710464) q[2];
sx q[2];
rz(-1.5652547) q[2];
sx q[2];
rz(1.8143285) q[2];
rz(-0.89613468) q[3];
sx q[3];
rz(-1.5265161) q[3];
sx q[3];
rz(2.0731887) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

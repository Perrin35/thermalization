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
rz(1.6506305) q[0];
sx q[0];
rz(-1.4181674) q[0];
sx q[0];
rz(11.081063) q[0];
rz(1.0506884) q[1];
sx q[1];
rz(-1.7424072) q[1];
sx q[1];
rz(-0.73837003) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90621072) q[0];
sx q[0];
rz(-2.8538508) q[0];
sx q[0];
rz(0.55858992) q[0];
x q[1];
rz(0.15282571) q[2];
sx q[2];
rz(-1.9597561) q[2];
sx q[2];
rz(-2.9940384) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8833137) q[1];
sx q[1];
rz(-1.4313233) q[1];
sx q[1];
rz(0.44513227) q[1];
rz(-pi) q[2];
rz(-0.61951903) q[3];
sx q[3];
rz(-1.8219007) q[3];
sx q[3];
rz(2.8015603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.38306132) q[2];
sx q[2];
rz(-2.8585275) q[2];
sx q[2];
rz(0.4134678) q[2];
rz(-2.5773898) q[3];
sx q[3];
rz(-1.4671624) q[3];
sx q[3];
rz(-1.1845425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7269932) q[0];
sx q[0];
rz(-2.0818384) q[0];
sx q[0];
rz(-0.10398908) q[0];
rz(3.0005786) q[1];
sx q[1];
rz(-0.68599373) q[1];
sx q[1];
rz(-0.41735059) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7234897) q[0];
sx q[0];
rz(-1.8716629) q[0];
sx q[0];
rz(1.0194433) q[0];
rz(-pi) q[1];
rz(-1.9345788) q[2];
sx q[2];
rz(-1.5593537) q[2];
sx q[2];
rz(-1.4909084) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.88637832) q[1];
sx q[1];
rz(-0.30579771) q[1];
sx q[1];
rz(-1.7096108) q[1];
rz(-pi) q[2];
rz(1.0052698) q[3];
sx q[3];
rz(-2.0291174) q[3];
sx q[3];
rz(1.2456196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.523681) q[2];
sx q[2];
rz(-2.1126426) q[2];
sx q[2];
rz(0.20393142) q[2];
rz(1.6265053) q[3];
sx q[3];
rz(-0.47658673) q[3];
sx q[3];
rz(-1.6319857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.5390891) q[0];
sx q[0];
rz(-1.4668377) q[0];
sx q[0];
rz(0.3983101) q[0];
rz(1.0771982) q[1];
sx q[1];
rz(-0.64882433) q[1];
sx q[1];
rz(0.55346742) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7837778) q[0];
sx q[0];
rz(-0.73448616) q[0];
sx q[0];
rz(-3.1260043) q[0];
x q[1];
rz(-2.8936549) q[2];
sx q[2];
rz(-0.72681475) q[2];
sx q[2];
rz(-0.13927973) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.043442) q[1];
sx q[1];
rz(-0.69313184) q[1];
sx q[1];
rz(0.56482238) q[1];
x q[2];
rz(-2.9514246) q[3];
sx q[3];
rz(-1.7432508) q[3];
sx q[3];
rz(1.8784539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.074097721) q[2];
sx q[2];
rz(-2.7390538) q[2];
sx q[2];
rz(1.925776) q[2];
rz(2.1408234) q[3];
sx q[3];
rz(-1.7583022) q[3];
sx q[3];
rz(1.2522987) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85553402) q[0];
sx q[0];
rz(-2.0834041) q[0];
sx q[0];
rz(-1.4372987) q[0];
rz(1.8057436) q[1];
sx q[1];
rz(-0.62594405) q[1];
sx q[1];
rz(-2.6864973) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2360596) q[0];
sx q[0];
rz(-0.62232557) q[0];
sx q[0];
rz(-0.13747352) q[0];
rz(2.0680769) q[2];
sx q[2];
rz(-0.4182564) q[2];
sx q[2];
rz(-0.69415316) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7981209) q[1];
sx q[1];
rz(-1.2373072) q[1];
sx q[1];
rz(1.8300959) q[1];
rz(-pi) q[2];
rz(-0.20858553) q[3];
sx q[3];
rz(-1.4781532) q[3];
sx q[3];
rz(-1.988033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2745634) q[2];
sx q[2];
rz(-0.44963351) q[2];
sx q[2];
rz(2.3095798) q[2];
rz(2.4882107) q[3];
sx q[3];
rz(-0.91975206) q[3];
sx q[3];
rz(0.67460361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5117383) q[0];
sx q[0];
rz(-1.3146223) q[0];
sx q[0];
rz(2.4526556) q[0];
rz(-2.9298933) q[1];
sx q[1];
rz(-1.7023106) q[1];
sx q[1];
rz(-2.6874218) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70063299) q[0];
sx q[0];
rz(-1.129375) q[0];
sx q[0];
rz(-0.9671797) q[0];
rz(2.7709097) q[2];
sx q[2];
rz(-0.96763583) q[2];
sx q[2];
rz(-2.2070845) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.59133021) q[1];
sx q[1];
rz(-2.2095517) q[1];
sx q[1];
rz(0.62229054) q[1];
rz(0.16726475) q[3];
sx q[3];
rz(-1.8297046) q[3];
sx q[3];
rz(2.6791999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6058558) q[2];
sx q[2];
rz(-1.7822632) q[2];
sx q[2];
rz(0.5109171) q[2];
rz(1.2153252) q[3];
sx q[3];
rz(-1.5646224) q[3];
sx q[3];
rz(-0.63466614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11584347) q[0];
sx q[0];
rz(-0.31084335) q[0];
sx q[0];
rz(-2.6281443) q[0];
rz(-0.91649857) q[1];
sx q[1];
rz(-1.1767574) q[1];
sx q[1];
rz(1.7880012) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2368456) q[0];
sx q[0];
rz(-1.488027) q[0];
sx q[0];
rz(-2.6968357) q[0];
rz(-1.9802092) q[2];
sx q[2];
rz(-1.231603) q[2];
sx q[2];
rz(-1.8095176) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60823764) q[1];
sx q[1];
rz(-1.5441193) q[1];
sx q[1];
rz(-0.4877301) q[1];
rz(-pi) q[2];
rz(1.4818125) q[3];
sx q[3];
rz(-1.6355733) q[3];
sx q[3];
rz(1.806206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6650271) q[2];
sx q[2];
rz(-0.55321425) q[2];
sx q[2];
rz(0.033585699) q[2];
rz(2.636886) q[3];
sx q[3];
rz(-1.7028156) q[3];
sx q[3];
rz(-2.3869042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1228834) q[0];
sx q[0];
rz(-1.9174734) q[0];
sx q[0];
rz(-1.9710185) q[0];
rz(2.9508044) q[1];
sx q[1];
rz(-2.6759594) q[1];
sx q[1];
rz(0.64291397) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2633154) q[0];
sx q[0];
rz(-1.8366792) q[0];
sx q[0];
rz(-3.0356867) q[0];
rz(-0.25919886) q[2];
sx q[2];
rz(-1.2507696) q[2];
sx q[2];
rz(-0.3896524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9289405) q[1];
sx q[1];
rz(-1.8096605) q[1];
sx q[1];
rz(0.3122621) q[1];
rz(2.8689485) q[3];
sx q[3];
rz(-2.7231403) q[3];
sx q[3];
rz(-1.2116644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1428895) q[2];
sx q[2];
rz(-0.10991749) q[2];
sx q[2];
rz(-1.1419123) q[2];
rz(-0.029684639) q[3];
sx q[3];
rz(-1.6073062) q[3];
sx q[3];
rz(-2.3619385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2340045) q[0];
sx q[0];
rz(-1.9588082) q[0];
sx q[0];
rz(-1.8335861) q[0];
rz(-2.0246778) q[1];
sx q[1];
rz(-1.6447379) q[1];
sx q[1];
rz(2.9294779) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.708688) q[0];
sx q[0];
rz(-1.2554662) q[0];
sx q[0];
rz(-0.84724119) q[0];
rz(-pi) q[1];
rz(0.64977744) q[2];
sx q[2];
rz(-0.37576518) q[2];
sx q[2];
rz(2.3670769) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6740283) q[1];
sx q[1];
rz(-1.0564359) q[1];
sx q[1];
rz(2.6326724) q[1];
rz(2.4808933) q[3];
sx q[3];
rz(-2.1815119) q[3];
sx q[3];
rz(-0.92044059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9461296) q[2];
sx q[2];
rz(-1.7886432) q[2];
sx q[2];
rz(-1.2809666) q[2];
rz(1.403275) q[3];
sx q[3];
rz(-0.91126982) q[3];
sx q[3];
rz(2.4173071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69287777) q[0];
sx q[0];
rz(-0.71664482) q[0];
sx q[0];
rz(2.2591059) q[0];
rz(-2.8485883) q[1];
sx q[1];
rz(-0.92331433) q[1];
sx q[1];
rz(-2.015347) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6544271) q[0];
sx q[0];
rz(-1.9357271) q[0];
sx q[0];
rz(-2.4449744) q[0];
rz(-0.0056929767) q[2];
sx q[2];
rz(-2.5416592) q[2];
sx q[2];
rz(-1.1545187) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4418151) q[1];
sx q[1];
rz(-2.338577) q[1];
sx q[1];
rz(-3.0670428) q[1];
rz(-1.8873439) q[3];
sx q[3];
rz(-1.1831565) q[3];
sx q[3];
rz(0.46771184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61665159) q[2];
sx q[2];
rz(-2.5280819) q[2];
sx q[2];
rz(-2.9743527) q[2];
rz(-3.1104769) q[3];
sx q[3];
rz(-1.1789221) q[3];
sx q[3];
rz(2.4955366) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6328218) q[0];
sx q[0];
rz(-1.333586) q[0];
sx q[0];
rz(2.3311145) q[0];
rz(1.3088016) q[1];
sx q[1];
rz(-0.93593132) q[1];
sx q[1];
rz(3.0442309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33695147) q[0];
sx q[0];
rz(-1.7367678) q[0];
sx q[0];
rz(2.354855) q[0];
x q[1];
rz(1.6426769) q[2];
sx q[2];
rz(-2.7595363) q[2];
sx q[2];
rz(-2.9660513) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3607935) q[1];
sx q[1];
rz(-0.9876087) q[1];
sx q[1];
rz(1.6287644) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1646268) q[3];
sx q[3];
rz(-2.8517234) q[3];
sx q[3];
rz(-2.4932662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.29791609) q[2];
sx q[2];
rz(-2.127779) q[2];
sx q[2];
rz(1.8664912) q[2];
rz(-2.7013333) q[3];
sx q[3];
rz(-2.2831254) q[3];
sx q[3];
rz(0.27537235) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1210099) q[0];
sx q[0];
rz(-0.6211716) q[0];
sx q[0];
rz(2.4564263) q[0];
rz(0.62200017) q[1];
sx q[1];
rz(-1.8432462) q[1];
sx q[1];
rz(-1.8274399) q[1];
rz(-2.1828281) q[2];
sx q[2];
rz(-3.0059881) q[2];
sx q[2];
rz(-0.81632951) q[2];
rz(-0.22777423) q[3];
sx q[3];
rz(-2.4854599) q[3];
sx q[3];
rz(1.5067185) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

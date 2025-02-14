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
rz(2.3447073) q[0];
sx q[0];
rz(-1.7303884) q[0];
sx q[0];
rz(2.222173) q[0];
rz(1.3031651) q[1];
sx q[1];
rz(-0.04920955) q[1];
sx q[1];
rz(2.4278909) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588736) q[0];
sx q[0];
rz(-1.2964411) q[0];
sx q[0];
rz(-2.9885942) q[0];
x q[1];
rz(-1.4422779) q[2];
sx q[2];
rz(-1.6331731) q[2];
sx q[2];
rz(-0.85278748) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1459058) q[1];
sx q[1];
rz(-0.45098454) q[1];
sx q[1];
rz(0.0031975071) q[1];
rz(-pi) q[2];
rz(1.0994751) q[3];
sx q[3];
rz(-1.5402392) q[3];
sx q[3];
rz(-1.9224642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.087622341) q[2];
sx q[2];
rz(-0.38603187) q[2];
sx q[2];
rz(0.36960441) q[2];
rz(-1.9965648) q[3];
sx q[3];
rz(-1.6686882) q[3];
sx q[3];
rz(-2.7692774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0358148) q[0];
sx q[0];
rz(-1.6749629) q[0];
sx q[0];
rz(-2.5685487) q[0];
rz(-2.0447958) q[1];
sx q[1];
rz(-1.6005452) q[1];
sx q[1];
rz(-0.47164741) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7268611) q[0];
sx q[0];
rz(-1.6808173) q[0];
sx q[0];
rz(-0.96745205) q[0];
x q[1];
rz(-0.69219442) q[2];
sx q[2];
rz(-1.0613777) q[2];
sx q[2];
rz(0.30864) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9976366) q[1];
sx q[1];
rz(-2.156971) q[1];
sx q[1];
rz(0.8950986) q[1];
rz(1.0466659) q[3];
sx q[3];
rz(-1.9283623) q[3];
sx q[3];
rz(-2.3563354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4196709) q[2];
sx q[2];
rz(-0.93908834) q[2];
sx q[2];
rz(-2.0527077) q[2];
rz(2.0770843) q[3];
sx q[3];
rz(-2.0239315) q[3];
sx q[3];
rz(-0.3256807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0786781) q[0];
sx q[0];
rz(-2.9525472) q[0];
sx q[0];
rz(3.1245226) q[0];
rz(0.71006376) q[1];
sx q[1];
rz(-2.3080669) q[1];
sx q[1];
rz(2.5753218) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6951272) q[0];
sx q[0];
rz(-2.4003865) q[0];
sx q[0];
rz(-1.5873853) q[0];
rz(2.429661) q[2];
sx q[2];
rz(-1.997479) q[2];
sx q[2];
rz(-0.79244765) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8108845) q[1];
sx q[1];
rz(-0.20527923) q[1];
sx q[1];
rz(2.7122871) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8999103) q[3];
sx q[3];
rz(-2.8688326) q[3];
sx q[3];
rz(-0.54995104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.337734) q[2];
sx q[2];
rz(-2.2457819) q[2];
sx q[2];
rz(-3.1324978) q[2];
rz(-3.005262) q[3];
sx q[3];
rz(-0.74458849) q[3];
sx q[3];
rz(-1.9280619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.424054) q[0];
sx q[0];
rz(-0.67988765) q[0];
sx q[0];
rz(-2.9754382) q[0];
rz(-2.0280139) q[1];
sx q[1];
rz(-2.6515617) q[1];
sx q[1];
rz(2.9794433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55570468) q[0];
sx q[0];
rz(-1.7506208) q[0];
sx q[0];
rz(0.094688133) q[0];
rz(-3.1012898) q[2];
sx q[2];
rz(-0.54931123) q[2];
sx q[2];
rz(-1.2072762) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2723586) q[1];
sx q[1];
rz(-1.2125535) q[1];
sx q[1];
rz(-0.76613249) q[1];
rz(-pi) q[2];
rz(2.6147551) q[3];
sx q[3];
rz(-0.73606811) q[3];
sx q[3];
rz(1.3730622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1589511) q[2];
sx q[2];
rz(-1.0794159) q[2];
sx q[2];
rz(-0.52784935) q[2];
rz(2.2900901) q[3];
sx q[3];
rz(-0.32367555) q[3];
sx q[3];
rz(-0.94312704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8714137) q[0];
sx q[0];
rz(-1.5988007) q[0];
sx q[0];
rz(0.80108368) q[0];
rz(-2.357645) q[1];
sx q[1];
rz(-0.74257094) q[1];
sx q[1];
rz(2.1606826) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50612569) q[0];
sx q[0];
rz(-2.6775595) q[0];
sx q[0];
rz(-2.6543174) q[0];
rz(-pi) q[1];
rz(-1.0641685) q[2];
sx q[2];
rz(-1.2384103) q[2];
sx q[2];
rz(-1.6695963) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66554932) q[1];
sx q[1];
rz(-0.11075039) q[1];
sx q[1];
rz(1.0057463) q[1];
rz(-pi) q[2];
rz(1.0193056) q[3];
sx q[3];
rz(-1.3634342) q[3];
sx q[3];
rz(-0.99121782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0826147) q[2];
sx q[2];
rz(-1.5951944) q[2];
sx q[2];
rz(-2.2583466) q[2];
rz(-0.20032459) q[3];
sx q[3];
rz(-0.89503461) q[3];
sx q[3];
rz(-2.0506355) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9969295) q[0];
sx q[0];
rz(-1.0522333) q[0];
sx q[0];
rz(-2.1494179) q[0];
rz(-0.80157533) q[1];
sx q[1];
rz(-1.293332) q[1];
sx q[1];
rz(1.4206402) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0376037) q[0];
sx q[0];
rz(-2.1891356) q[0];
sx q[0];
rz(-1.9644587) q[0];
rz(-1.0384485) q[2];
sx q[2];
rz(-1.6611929) q[2];
sx q[2];
rz(0.20648512) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3657377) q[1];
sx q[1];
rz(-0.72682646) q[1];
sx q[1];
rz(2.6353542) q[1];
rz(3.104131) q[3];
sx q[3];
rz(-1.6679296) q[3];
sx q[3];
rz(-2.284632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.27986032) q[2];
sx q[2];
rz(-1.4272775) q[2];
sx q[2];
rz(-2.6045065) q[2];
rz(-3.1308657) q[3];
sx q[3];
rz(-2.3948632) q[3];
sx q[3];
rz(-0.85095325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2456197) q[0];
sx q[0];
rz(-2.8697822) q[0];
sx q[0];
rz(-2.8017092) q[0];
rz(-1.8396359) q[1];
sx q[1];
rz(-2.1959031) q[1];
sx q[1];
rz(1.9702912) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7599277) q[0];
sx q[0];
rz(-1.6070843) q[0];
sx q[0];
rz(2.8445303) q[0];
rz(-pi) q[1];
rz(-1.3572071) q[2];
sx q[2];
rz(-1.3554975) q[2];
sx q[2];
rz(-2.4395296) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.724546) q[1];
sx q[1];
rz(-1.5929211) q[1];
sx q[1];
rz(-1.6032277) q[1];
rz(2.5170691) q[3];
sx q[3];
rz(-1.7134662) q[3];
sx q[3];
rz(2.4017911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.12477144) q[2];
sx q[2];
rz(-1.4171968) q[2];
sx q[2];
rz(-2.4701414) q[2];
rz(-0.5365544) q[3];
sx q[3];
rz(-0.96009976) q[3];
sx q[3];
rz(-1.051739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3143828) q[0];
sx q[0];
rz(-2.0874513) q[0];
sx q[0];
rz(-0.89624727) q[0];
rz(-0.25262901) q[1];
sx q[1];
rz(-0.8600421) q[1];
sx q[1];
rz(2.0910697) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69347113) q[0];
sx q[0];
rz(-0.52883178) q[0];
sx q[0];
rz(3.0419087) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4765313) q[2];
sx q[2];
rz(-1.6740693) q[2];
sx q[2];
rz(-2.2879083) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.76681644) q[1];
sx q[1];
rz(-1.4207867) q[1];
sx q[1];
rz(-2.0376202) q[1];
rz(-pi) q[2];
rz(-3.0862634) q[3];
sx q[3];
rz(-0.75859514) q[3];
sx q[3];
rz(-0.38842312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5184021) q[2];
sx q[2];
rz(-2.9727327) q[2];
sx q[2];
rz(-1.5625578) q[2];
rz(1.7744428) q[3];
sx q[3];
rz(-1.046215) q[3];
sx q[3];
rz(0.76213837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52687454) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(0.38994625) q[0];
rz(0.82603106) q[1];
sx q[1];
rz(-0.69458687) q[1];
sx q[1];
rz(1.0521851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56985939) q[0];
sx q[0];
rz(-0.89825199) q[0];
sx q[0];
rz(-2.74617) q[0];
rz(-pi) q[1];
rz(2.7446943) q[2];
sx q[2];
rz(-1.4087311) q[2];
sx q[2];
rz(1.4805178) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4349367) q[1];
sx q[1];
rz(-1.860992) q[1];
sx q[1];
rz(1.1654084) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1901699) q[3];
sx q[3];
rz(-2.4149899) q[3];
sx q[3];
rz(2.8539336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0612001) q[2];
sx q[2];
rz(-1.281851) q[2];
sx q[2];
rz(-0.52919394) q[2];
rz(2.3906294) q[3];
sx q[3];
rz(-2.1801528) q[3];
sx q[3];
rz(2.9474337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.879409) q[0];
sx q[0];
rz(-0.59469596) q[0];
sx q[0];
rz(1.1645114) q[0];
rz(1.2871845) q[1];
sx q[1];
rz(-0.6856122) q[1];
sx q[1];
rz(0.50416344) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39162779) q[0];
sx q[0];
rz(-2.7996783) q[0];
sx q[0];
rz(-2.2750399) q[0];
rz(-pi) q[1];
rz(-2.7648534) q[2];
sx q[2];
rz(-1.5887056) q[2];
sx q[2];
rz(1.1056545) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2626942) q[1];
sx q[1];
rz(-2.8126908) q[1];
sx q[1];
rz(-0.94543381) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21728094) q[3];
sx q[3];
rz(-2.2437527) q[3];
sx q[3];
rz(-0.5289883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.21504271) q[2];
sx q[2];
rz(-1.5554917) q[2];
sx q[2];
rz(3.1070993) q[2];
rz(1.5283594) q[3];
sx q[3];
rz(-0.72051636) q[3];
sx q[3];
rz(0.31149402) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16019776) q[0];
sx q[0];
rz(-1.5934516) q[0];
sx q[0];
rz(-0.39504575) q[0];
rz(-2.1389217) q[1];
sx q[1];
rz(-1.4613338) q[1];
sx q[1];
rz(2.7636539) q[1];
rz(-2.7591697) q[2];
sx q[2];
rz(-0.48880807) q[2];
sx q[2];
rz(0.53570408) q[2];
rz(-1.5080084) q[3];
sx q[3];
rz(-1.8973402) q[3];
sx q[3];
rz(-0.8473581) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

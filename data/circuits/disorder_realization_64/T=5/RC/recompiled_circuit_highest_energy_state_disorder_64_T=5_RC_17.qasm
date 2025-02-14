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
rz(-0.79688537) q[0];
sx q[0];
rz(4.8719811) q[0];
sx q[0];
rz(13.48579) q[0];
rz(1.3031651) q[1];
sx q[1];
rz(-0.04920955) q[1];
sx q[1];
rz(-0.71370178) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2588736) q[0];
sx q[0];
rz(-1.2964411) q[0];
sx q[0];
rz(-0.15299847) q[0];
x q[1];
rz(-3.0786985) q[2];
sx q[2];
rz(-1.4425292) q[2];
sx q[2];
rz(-2.4155282) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9992396) q[1];
sx q[1];
rz(-1.1198143) q[1];
sx q[1];
rz(-1.5723448) q[1];
x q[2];
rz(3.1072989) q[3];
sx q[3];
rz(-2.0418797) q[3];
sx q[3];
rz(-0.33609875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0539703) q[2];
sx q[2];
rz(-2.7555608) q[2];
sx q[2];
rz(0.36960441) q[2];
rz(1.1450279) q[3];
sx q[3];
rz(-1.4729045) q[3];
sx q[3];
rz(-0.37231529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10577781) q[0];
sx q[0];
rz(-1.4666297) q[0];
sx q[0];
rz(-0.57304397) q[0];
rz(2.0447958) q[1];
sx q[1];
rz(-1.6005452) q[1];
sx q[1];
rz(0.47164741) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9100138) q[0];
sx q[0];
rz(-0.97161228) q[0];
sx q[0];
rz(-3.0082361) q[0];
x q[1];
rz(-0.69219442) q[2];
sx q[2];
rz(-2.080215) q[2];
sx q[2];
rz(2.8329527) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9976366) q[1];
sx q[1];
rz(-2.156971) q[1];
sx q[1];
rz(-2.2464941) q[1];
rz(2.0949267) q[3];
sx q[3];
rz(-1.2132304) q[3];
sx q[3];
rz(0.78525728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4196709) q[2];
sx q[2];
rz(-0.93908834) q[2];
sx q[2];
rz(-1.088885) q[2];
rz(2.0770843) q[3];
sx q[3];
rz(-1.1176611) q[3];
sx q[3];
rz(-2.815912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0629145) q[0];
sx q[0];
rz(-0.18904541) q[0];
sx q[0];
rz(3.1245226) q[0];
rz(2.4315289) q[1];
sx q[1];
rz(-0.83352572) q[1];
sx q[1];
rz(2.5753218) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44646548) q[0];
sx q[0];
rz(-2.4003865) q[0];
sx q[0];
rz(1.5542073) q[0];
rz(1.0300358) q[2];
sx q[2];
rz(-0.93387253) q[2];
sx q[2];
rz(-2.7062397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.480174) q[1];
sx q[1];
rz(-1.4858477) q[1];
sx q[1];
rz(-0.18710356) q[1];
rz(2.8763883) q[3];
sx q[3];
rz(-1.6353161) q[3];
sx q[3];
rz(-1.2539188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.337734) q[2];
sx q[2];
rz(-2.2457819) q[2];
sx q[2];
rz(0.0090948661) q[2];
rz(3.005262) q[3];
sx q[3];
rz(-0.74458849) q[3];
sx q[3];
rz(1.9280619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71753865) q[0];
sx q[0];
rz(-0.67988765) q[0];
sx q[0];
rz(0.16615443) q[0];
rz(-2.0280139) q[1];
sx q[1];
rz(-2.6515617) q[1];
sx q[1];
rz(-0.16214935) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.585888) q[0];
sx q[0];
rz(-1.7506208) q[0];
sx q[0];
rz(-3.0469045) q[0];
rz(-1.5954564) q[2];
sx q[2];
rz(-2.1196105) q[2];
sx q[2];
rz(-1.8870712) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.869234) q[1];
sx q[1];
rz(-1.2125535) q[1];
sx q[1];
rz(0.76613249) q[1];
rz(1.1433853) q[3];
sx q[3];
rz(-2.1899438) q[3];
sx q[3];
rz(-1.1030847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98264155) q[2];
sx q[2];
rz(-2.0621767) q[2];
sx q[2];
rz(2.6137433) q[2];
rz(-0.85150254) q[3];
sx q[3];
rz(-2.8179171) q[3];
sx q[3];
rz(-2.1984656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8714137) q[0];
sx q[0];
rz(-1.542792) q[0];
sx q[0];
rz(0.80108368) q[0];
rz(0.78394765) q[1];
sx q[1];
rz(-2.3990217) q[1];
sx q[1];
rz(0.98091006) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5194395) q[0];
sx q[0];
rz(-1.781917) q[0];
sx q[0];
rz(-0.4163752) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0774242) q[2];
sx q[2];
rz(-1.2384103) q[2];
sx q[2];
rz(1.6695963) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4675216) q[1];
sx q[1];
rz(-1.630012) q[1];
sx q[1];
rz(1.4771512) q[1];
rz(-pi) q[2];
rz(1.1889691) q[3];
sx q[3];
rz(-2.5562048) q[3];
sx q[3];
rz(-0.25661925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.058978) q[2];
sx q[2];
rz(-1.5951944) q[2];
sx q[2];
rz(-2.2583466) q[2];
rz(-2.9412681) q[3];
sx q[3];
rz(-0.89503461) q[3];
sx q[3];
rz(2.0506355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14466318) q[0];
sx q[0];
rz(-1.0522333) q[0];
sx q[0];
rz(2.1494179) q[0];
rz(-0.80157533) q[1];
sx q[1];
rz(-1.293332) q[1];
sx q[1];
rz(-1.7209524) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4822749) q[0];
sx q[0];
rz(-2.4226696) q[0];
sx q[0];
rz(2.6470967) q[0];
rz(-pi) q[1];
rz(3.0367766) q[2];
sx q[2];
rz(-2.1007406) q[2];
sx q[2];
rz(1.8304093) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0040505) q[1];
sx q[1];
rz(-0.95064771) q[1];
sx q[1];
rz(1.9778848) q[1];
x q[2];
rz(1.2038369) q[3];
sx q[3];
rz(-3.0375069) q[3];
sx q[3];
rz(2.6534125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8617323) q[2];
sx q[2];
rz(-1.4272775) q[2];
sx q[2];
rz(-2.6045065) q[2];
rz(3.1308657) q[3];
sx q[3];
rz(-2.3948632) q[3];
sx q[3];
rz(0.85095325) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.895973) q[0];
sx q[0];
rz(-2.8697822) q[0];
sx q[0];
rz(0.33988345) q[0];
rz(-1.3019568) q[1];
sx q[1];
rz(-2.1959031) q[1];
sx q[1];
rz(1.1713015) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7599277) q[0];
sx q[0];
rz(-1.6070843) q[0];
sx q[0];
rz(0.29706232) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9214462) q[2];
sx q[2];
rz(-1.7793806) q[2];
sx q[2];
rz(2.2265547) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.75222941) q[1];
sx q[1];
rz(-3.1023355) q[1];
sx q[1];
rz(-2.1696349) q[1];
rz(2.5170691) q[3];
sx q[3];
rz(-1.4281264) q[3];
sx q[3];
rz(-2.4017911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0168212) q[2];
sx q[2];
rz(-1.4171968) q[2];
sx q[2];
rz(-0.67145124) q[2];
rz(0.5365544) q[3];
sx q[3];
rz(-2.1814929) q[3];
sx q[3];
rz(-1.051739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8272098) q[0];
sx q[0];
rz(-1.0541414) q[0];
sx q[0];
rz(-2.2453454) q[0];
rz(-2.8889636) q[1];
sx q[1];
rz(-2.2815506) q[1];
sx q[1];
rz(2.0910697) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3327961) q[0];
sx q[0];
rz(-1.044863) q[0];
sx q[0];
rz(-1.6288847) q[0];
rz(-pi) q[1];
rz(-0.66506135) q[2];
sx q[2];
rz(-1.6740693) q[2];
sx q[2];
rz(-0.85368431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3747762) q[1];
sx q[1];
rz(-1.720806) q[1];
sx q[1];
rz(2.0376202) q[1];
rz(-3.0862634) q[3];
sx q[3];
rz(-0.75859514) q[3];
sx q[3];
rz(2.7531695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5184021) q[2];
sx q[2];
rz(-2.9727327) q[2];
sx q[2];
rz(-1.5790348) q[2];
rz(1.7744428) q[3];
sx q[3];
rz(-1.046215) q[3];
sx q[3];
rz(0.76213837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6147181) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(2.7516464) q[0];
rz(2.3155616) q[1];
sx q[1];
rz(-0.69458687) q[1];
sx q[1];
rz(2.0894076) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56985939) q[0];
sx q[0];
rz(-0.89825199) q[0];
sx q[0];
rz(-2.74617) q[0];
rz(2.7446943) q[2];
sx q[2];
rz(-1.7328615) q[2];
sx q[2];
rz(-1.4805178) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6893951) q[1];
sx q[1];
rz(-0.493825) q[1];
sx q[1];
rz(-0.92269519) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47635079) q[3];
sx q[3];
rz(-2.1423376) q[3];
sx q[3];
rz(2.0921354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0803926) q[2];
sx q[2];
rz(-1.8597417) q[2];
sx q[2];
rz(2.6123987) q[2];
rz(2.3906294) q[3];
sx q[3];
rz(-2.1801528) q[3];
sx q[3];
rz(2.9474337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26218364) q[0];
sx q[0];
rz(-2.5468967) q[0];
sx q[0];
rz(-1.1645114) q[0];
rz(-1.2871845) q[1];
sx q[1];
rz(-0.6856122) q[1];
sx q[1];
rz(2.6374292) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39162779) q[0];
sx q[0];
rz(-2.7996783) q[0];
sx q[0];
rz(0.86655273) q[0];
x q[1];
rz(-2.7648534) q[2];
sx q[2];
rz(-1.5887056) q[2];
sx q[2];
rz(-2.0359382) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6109687) q[1];
sx q[1];
rz(-1.8357616) q[1];
sx q[1];
rz(-0.19719657) q[1];
rz(-pi) q[2];
rz(0.88621288) q[3];
sx q[3];
rz(-1.7401812) q[3];
sx q[3];
rz(2.2365295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9265499) q[2];
sx q[2];
rz(-1.5861009) q[2];
sx q[2];
rz(0.034493383) q[2];
rz(-1.5283594) q[3];
sx q[3];
rz(-0.72051636) q[3];
sx q[3];
rz(2.8300986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
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
rz(2.7591697) q[2];
sx q[2];
rz(-2.6527846) q[2];
sx q[2];
rz(-2.6058886) q[2];
rz(2.9583951) q[3];
sx q[3];
rz(-2.8092794) q[3];
sx q[3];
rz(2.487779) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

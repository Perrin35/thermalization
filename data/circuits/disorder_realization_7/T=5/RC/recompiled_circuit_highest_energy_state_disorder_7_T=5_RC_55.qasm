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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58881271) q[0];
sx q[0];
rz(-2.3511887) q[0];
sx q[0];
rz(-1.3270204) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6381016) q[2];
sx q[2];
rz(-1.8241183) q[2];
sx q[2];
rz(1.5510786) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8350153) q[1];
sx q[1];
rz(-1.6284429) q[1];
sx q[1];
rz(-0.60204864) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0996885) q[3];
sx q[3];
rz(-2.503241) q[3];
sx q[3];
rz(-1.1717351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2446186) q[2];
sx q[2];
rz(-1.3099058) q[2];
sx q[2];
rz(0.2429602) q[2];
rz(2.7729559) q[3];
sx q[3];
rz(-0.60550767) q[3];
sx q[3];
rz(-1.9093556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.826137) q[0];
sx q[0];
rz(-0.22268) q[0];
sx q[0];
rz(-2.7742703) q[0];
rz(2.3452554) q[1];
sx q[1];
rz(-1.0581191) q[1];
sx q[1];
rz(-2.499089) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8272907) q[0];
sx q[0];
rz(-1.9901942) q[0];
sx q[0];
rz(-2.6227996) q[0];
rz(-2.5359277) q[2];
sx q[2];
rz(-0.64222756) q[2];
sx q[2];
rz(-0.58711827) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6141908) q[1];
sx q[1];
rz(-2.3878015) q[1];
sx q[1];
rz(-0.68282737) q[1];
rz(-3.0152937) q[3];
sx q[3];
rz(-1.0923315) q[3];
sx q[3];
rz(-2.0984405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.6392886) q[2];
sx q[2];
rz(-2.2059811) q[2];
sx q[2];
rz(1.3703692) q[2];
rz(-0.227452) q[3];
sx q[3];
rz(-1.2496313) q[3];
sx q[3];
rz(-1.9500218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74533904) q[0];
sx q[0];
rz(-2.9662913) q[0];
sx q[0];
rz(-1.1981717) q[0];
rz(1.0802957) q[1];
sx q[1];
rz(-0.21427576) q[1];
sx q[1];
rz(3.0467196) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.085233) q[0];
sx q[0];
rz(-1.0329536) q[0];
sx q[0];
rz(0.30098029) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0804786) q[2];
sx q[2];
rz(-0.72619263) q[2];
sx q[2];
rz(-2.2167571) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.87743041) q[1];
sx q[1];
rz(-2.1070711) q[1];
sx q[1];
rz(-0.047604851) q[1];
rz(-pi) q[2];
rz(2.8646144) q[3];
sx q[3];
rz(-1.2001032) q[3];
sx q[3];
rz(1.0007559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0723116) q[2];
sx q[2];
rz(-2.5371964) q[2];
sx q[2];
rz(2.7351232) q[2];
rz(-1.1786002) q[3];
sx q[3];
rz(-1.4402163) q[3];
sx q[3];
rz(-1.9788205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5793295) q[0];
sx q[0];
rz(-0.13450204) q[0];
sx q[0];
rz(-0.33600268) q[0];
rz(-0.52040368) q[1];
sx q[1];
rz(-2.2774179) q[1];
sx q[1];
rz(-0.10428183) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18096237) q[0];
sx q[0];
rz(-1.917836) q[0];
sx q[0];
rz(-0.35023679) q[0];
rz(-3.0523446) q[2];
sx q[2];
rz(-0.1570905) q[2];
sx q[2];
rz(2.2721827) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1946009) q[1];
sx q[1];
rz(-2.2333988) q[1];
sx q[1];
rz(-2.6688982) q[1];
x q[2];
rz(0.4533259) q[3];
sx q[3];
rz(-1.6764055) q[3];
sx q[3];
rz(0.80704122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5229554) q[2];
sx q[2];
rz(-2.0433661) q[2];
sx q[2];
rz(-1.1445507) q[2];
rz(-2.5942904) q[3];
sx q[3];
rz(-1.9132883) q[3];
sx q[3];
rz(1.087629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2523786) q[0];
sx q[0];
rz(-2.9504898) q[0];
sx q[0];
rz(-2.5872173) q[0];
rz(0.91122183) q[1];
sx q[1];
rz(-1.8781885) q[1];
sx q[1];
rz(0.30141452) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49839334) q[0];
sx q[0];
rz(-1.3424771) q[0];
sx q[0];
rz(-0.96781815) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6658556) q[2];
sx q[2];
rz(-0.43125501) q[2];
sx q[2];
rz(-2.3088344) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44228783) q[1];
sx q[1];
rz(-1.7831452) q[1];
sx q[1];
rz(1.9599171) q[1];
rz(2.2603792) q[3];
sx q[3];
rz(-0.47166079) q[3];
sx q[3];
rz(-0.72615964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1201799) q[2];
sx q[2];
rz(-2.9605949) q[2];
sx q[2];
rz(-0.0069228355) q[2];
rz(-0.93112469) q[3];
sx q[3];
rz(-2.0102863) q[3];
sx q[3];
rz(-1.1091703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329426) q[0];
sx q[0];
rz(-2.3045492) q[0];
sx q[0];
rz(-1.5337926) q[0];
rz(-0.7026698) q[1];
sx q[1];
rz(-2.6103795) q[1];
sx q[1];
rz(-2.839397) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.142573) q[0];
sx q[0];
rz(-2.4732865) q[0];
sx q[0];
rz(-1.4060941) q[0];
rz(2.7415761) q[2];
sx q[2];
rz(-2.3850394) q[2];
sx q[2];
rz(-0.053002593) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2875048) q[1];
sx q[1];
rz(-2.4937428) q[1];
sx q[1];
rz(1.2811529) q[1];
rz(0.052107776) q[3];
sx q[3];
rz(-1.7223538) q[3];
sx q[3];
rz(-1.0406756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89221421) q[2];
sx q[2];
rz(-1.0503146) q[2];
sx q[2];
rz(-0.92602473) q[2];
rz(-1.2146436) q[3];
sx q[3];
rz(-2.6483783) q[3];
sx q[3];
rz(0.27629575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72325426) q[0];
sx q[0];
rz(-2.3746018) q[0];
sx q[0];
rz(-1.1085229) q[0];
rz(0.33379894) q[1];
sx q[1];
rz(-1.8148986) q[1];
sx q[1];
rz(1.0858067) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.007581) q[0];
sx q[0];
rz(-1.3659048) q[0];
sx q[0];
rz(1.4896859) q[0];
rz(-pi) q[1];
rz(1.1637229) q[2];
sx q[2];
rz(-0.43292871) q[2];
sx q[2];
rz(2.0100398) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27583308) q[1];
sx q[1];
rz(-2.1001022) q[1];
sx q[1];
rz(2.5342026) q[1];
rz(-1.3855188) q[3];
sx q[3];
rz(-2.4967125) q[3];
sx q[3];
rz(1.737999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6256025) q[2];
sx q[2];
rz(-2.7770999) q[2];
sx q[2];
rz(1.1642574) q[2];
rz(2.8734251) q[3];
sx q[3];
rz(-1.2428913) q[3];
sx q[3];
rz(-0.75781649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5371573) q[0];
sx q[0];
rz(-2.6217961) q[0];
sx q[0];
rz(1.9926158) q[0];
rz(0.2941429) q[1];
sx q[1];
rz(-1.5831169) q[1];
sx q[1];
rz(2.099096) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5940018) q[0];
sx q[0];
rz(-0.85193513) q[0];
sx q[0];
rz(-0.44525624) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78977079) q[2];
sx q[2];
rz(-1.4068479) q[2];
sx q[2];
rz(-2.3254834) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6175872) q[1];
sx q[1];
rz(-1.7442787) q[1];
sx q[1];
rz(0.37131997) q[1];
rz(-pi) q[2];
rz(1.5339666) q[3];
sx q[3];
rz(-2.4772518) q[3];
sx q[3];
rz(0.0078594154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5440172) q[2];
sx q[2];
rz(-0.67108265) q[2];
sx q[2];
rz(-1.3647122) q[2];
rz(2.4890066) q[3];
sx q[3];
rz(-2.3502974) q[3];
sx q[3];
rz(0.64835382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.515601) q[0];
sx q[0];
rz(-0.37726548) q[0];
sx q[0];
rz(3.1242477) q[0];
rz(-0.01677244) q[1];
sx q[1];
rz(-2.3544632) q[1];
sx q[1];
rz(0.15388547) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2543593) q[0];
sx q[0];
rz(-1.6463929) q[0];
sx q[0];
rz(2.8235954) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4520313) q[2];
sx q[2];
rz(-0.84367263) q[2];
sx q[2];
rz(2.6117556) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9388403) q[1];
sx q[1];
rz(-1.8649264) q[1];
sx q[1];
rz(-0.37962706) q[1];
rz(-pi) q[2];
x q[2];
rz(0.023970402) q[3];
sx q[3];
rz(-1.4949189) q[3];
sx q[3];
rz(0.45437231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.953557) q[2];
sx q[2];
rz(-0.81622684) q[2];
sx q[2];
rz(0.44450644) q[2];
rz(-2.9514173) q[3];
sx q[3];
rz(-1.5580274) q[3];
sx q[3];
rz(1.7688513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1403777) q[0];
sx q[0];
rz(-1.8920521) q[0];
sx q[0];
rz(2.0709399) q[0];
rz(-0.045914687) q[1];
sx q[1];
rz(-1.4713947) q[1];
sx q[1];
rz(-2.3505223) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57426737) q[0];
sx q[0];
rz(-1.2358032) q[0];
sx q[0];
rz(0.092760249) q[0];
rz(-2.5612381) q[2];
sx q[2];
rz(-0.46541801) q[2];
sx q[2];
rz(0.67829715) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6562499) q[1];
sx q[1];
rz(-0.54528499) q[1];
sx q[1];
rz(0.7363018) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13274712) q[3];
sx q[3];
rz(-1.069396) q[3];
sx q[3];
rz(0.30745922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4967686) q[2];
sx q[2];
rz(-2.9569929) q[2];
sx q[2];
rz(1.4727288) q[2];
rz(-1.6139) q[3];
sx q[3];
rz(-1.0563285) q[3];
sx q[3];
rz(-1.24019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2033757) q[0];
sx q[0];
rz(-1.4980409) q[0];
sx q[0];
rz(-0.71312755) q[0];
rz(-0.51364246) q[1];
sx q[1];
rz(-1.4331663) q[1];
sx q[1];
rz(-1.6930361) q[1];
rz(-2.7265139) q[2];
sx q[2];
rz(-2.5998301) q[2];
sx q[2];
rz(2.0143581) q[2];
rz(2.722693) q[3];
sx q[3];
rz(-2.2875026) q[3];
sx q[3];
rz(0.84411375) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

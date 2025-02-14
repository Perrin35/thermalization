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
rz(-2.1751997) q[0];
sx q[0];
rz(4.8973358) q[0];
sx q[0];
rz(10.020221) q[0];
rz(-1.8499941) q[1];
sx q[1];
rz(2.5505677) q[1];
sx q[1];
rz(12.443065) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5714471) q[0];
sx q[0];
rz(-2.4295904) q[0];
sx q[0];
rz(-0.90306905) q[0];
x q[1];
rz(2.3594728) q[2];
sx q[2];
rz(-0.95172666) q[2];
sx q[2];
rz(-1.182488) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.0088657503) q[1];
sx q[1];
rz(-1.6648597) q[1];
sx q[1];
rz(1.3558741) q[1];
x q[2];
rz(1.7233347) q[3];
sx q[3];
rz(-1.9328874) q[3];
sx q[3];
rz(-0.83521508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3143602) q[2];
sx q[2];
rz(-1.102697) q[2];
sx q[2];
rz(3.0296791) q[2];
rz(-1.3917475) q[3];
sx q[3];
rz(-1.1575969) q[3];
sx q[3];
rz(0.30645034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7258485) q[0];
sx q[0];
rz(-0.84472504) q[0];
sx q[0];
rz(0.14990212) q[0];
rz(-0.28597486) q[1];
sx q[1];
rz(-1.7466702) q[1];
sx q[1];
rz(-2.012595) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6647346) q[0];
sx q[0];
rz(-1.6010512) q[0];
sx q[0];
rz(1.2824351) q[0];
rz(-1.8104042) q[2];
sx q[2];
rz(-0.061366038) q[2];
sx q[2];
rz(2.1642016) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.45980487) q[1];
sx q[1];
rz(-0.39783172) q[1];
sx q[1];
rz(0.43983207) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82650696) q[3];
sx q[3];
rz(-0.77946589) q[3];
sx q[3];
rz(-2.731088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4631606) q[2];
sx q[2];
rz(-2.2191935) q[2];
sx q[2];
rz(1.1506259) q[2];
rz(-1.9567418) q[3];
sx q[3];
rz(-1.7246282) q[3];
sx q[3];
rz(-3.1209893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-0.48035574) q[0];
sx q[0];
rz(-0.24599563) q[0];
sx q[0];
rz(-0.38159698) q[0];
rz(-1.2941788) q[1];
sx q[1];
rz(-2.0353863) q[1];
sx q[1];
rz(0.19827422) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1564045) q[0];
sx q[0];
rz(-1.745178) q[0];
sx q[0];
rz(-1.8013493) q[0];
x q[1];
rz(-2.2402693) q[2];
sx q[2];
rz(-1.633989) q[2];
sx q[2];
rz(2.2128999) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5570017) q[1];
sx q[1];
rz(-1.7832568) q[1];
sx q[1];
rz(2.9032488) q[1];
x q[2];
rz(-1.395969) q[3];
sx q[3];
rz(-1.5906189) q[3];
sx q[3];
rz(-2.1682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76564378) q[2];
sx q[2];
rz(-0.15268923) q[2];
sx q[2];
rz(0.51885968) q[2];
rz(1.2505924) q[3];
sx q[3];
rz(-1.0873245) q[3];
sx q[3];
rz(-0.58108228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4984703) q[0];
sx q[0];
rz(-2.9415218) q[0];
sx q[0];
rz(-2.6923687) q[0];
rz(1.4462224) q[1];
sx q[1];
rz(-0.64165533) q[1];
sx q[1];
rz(0.62072388) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45263824) q[0];
sx q[0];
rz(-1.7151378) q[0];
sx q[0];
rz(-2.4792433) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16731842) q[2];
sx q[2];
rz(-0.64019055) q[2];
sx q[2];
rz(-0.13969914) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6171332) q[1];
sx q[1];
rz(-1.7726328) q[1];
sx q[1];
rz(-2.733888) q[1];
rz(-pi) q[2];
rz(2.1093217) q[3];
sx q[3];
rz(-1.5251951) q[3];
sx q[3];
rz(-0.63889438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0351403) q[2];
sx q[2];
rz(-1.2713212) q[2];
sx q[2];
rz(-2.980496) q[2];
rz(-2.7437239) q[3];
sx q[3];
rz(-1.4784003) q[3];
sx q[3];
rz(-0.14959344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.38254) q[0];
sx q[0];
rz(-1.362514) q[0];
sx q[0];
rz(0.25704849) q[0];
rz(-2.2611179) q[1];
sx q[1];
rz(-2.5011261) q[1];
sx q[1];
rz(1.2535198) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3220164) q[0];
sx q[0];
rz(-1.5347693) q[0];
sx q[0];
rz(-3.1234804) q[0];
x q[1];
rz(-0.80998151) q[2];
sx q[2];
rz(-1.4000386) q[2];
sx q[2];
rz(-0.25639889) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.67318007) q[1];
sx q[1];
rz(-2.367503) q[1];
sx q[1];
rz(-1.6353398) q[1];
rz(-0.0048794172) q[3];
sx q[3];
rz(-2.7944428) q[3];
sx q[3];
rz(-2.7018869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72838655) q[2];
sx q[2];
rz(-1.9325958) q[2];
sx q[2];
rz(-1.2665292) q[2];
rz(-2.5201216) q[3];
sx q[3];
rz(-0.60756835) q[3];
sx q[3];
rz(-2.6004041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7414339) q[0];
sx q[0];
rz(-2.3024237) q[0];
sx q[0];
rz(-0.025064502) q[0];
rz(-1.9301682) q[1];
sx q[1];
rz(-1.2592659) q[1];
sx q[1];
rz(0.2549583) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0023277442) q[0];
sx q[0];
rz(-3.0920017) q[0];
sx q[0];
rz(3.0399486) q[0];
rz(1.0827052) q[2];
sx q[2];
rz(-1.0399264) q[2];
sx q[2];
rz(1.4029897) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3820262) q[1];
sx q[1];
rz(-0.81795482) q[1];
sx q[1];
rz(-2.3809529) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99825777) q[3];
sx q[3];
rz(-1.3911595) q[3];
sx q[3];
rz(1.6999753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8956464) q[2];
sx q[2];
rz(-2.0993555) q[2];
sx q[2];
rz(-0.45905217) q[2];
rz(-1.4970655) q[3];
sx q[3];
rz(-0.11681695) q[3];
sx q[3];
rz(-3.0204401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6727305) q[0];
sx q[0];
rz(-2.6519096) q[0];
sx q[0];
rz(0.086061867) q[0];
rz(2.22279) q[1];
sx q[1];
rz(-1.1163534) q[1];
sx q[1];
rz(-1.2109717) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0491517) q[0];
sx q[0];
rz(-1.7730129) q[0];
sx q[0];
rz(-2.8114955) q[0];
x q[1];
rz(-0.93386704) q[2];
sx q[2];
rz(-2.4709765) q[2];
sx q[2];
rz(-3.0702555) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1001491) q[1];
sx q[1];
rz(-2.7584834) q[1];
sx q[1];
rz(0.99975296) q[1];
rz(-1.9350697) q[3];
sx q[3];
rz(-1.2281872) q[3];
sx q[3];
rz(2.2628257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.81898895) q[2];
sx q[2];
rz(-2.6456867) q[2];
sx q[2];
rz(0.71713478) q[2];
rz(-1.0273733) q[3];
sx q[3];
rz(-1.7337948) q[3];
sx q[3];
rz(2.7023442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9913919) q[0];
sx q[0];
rz(-0.99656492) q[0];
sx q[0];
rz(-0.014658654) q[0];
rz(2.390059) q[1];
sx q[1];
rz(-1.2208168) q[1];
sx q[1];
rz(-1.6709447) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96318564) q[0];
sx q[0];
rz(-1.7589594) q[0];
sx q[0];
rz(1.7117731) q[0];
rz(-pi) q[1];
rz(1.4906497) q[2];
sx q[2];
rz(-1.8266457) q[2];
sx q[2];
rz(2.3639774) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4336509) q[1];
sx q[1];
rz(-1.6835064) q[1];
sx q[1];
rz(1.1658843) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1378391) q[3];
sx q[3];
rz(-1.9158746) q[3];
sx q[3];
rz(2.0451289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8810001) q[2];
sx q[2];
rz(-1.1129817) q[2];
sx q[2];
rz(2.8680958) q[2];
rz(-1.986844) q[3];
sx q[3];
rz(-2.4879849) q[3];
sx q[3];
rz(-0.65034136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038789373) q[0];
sx q[0];
rz(-2.0162855) q[0];
sx q[0];
rz(-1.0754841) q[0];
rz(-0.12818809) q[1];
sx q[1];
rz(-0.85103858) q[1];
sx q[1];
rz(2.7959965) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3227279) q[0];
sx q[0];
rz(-0.68536192) q[0];
sx q[0];
rz(1.8249874) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1394452) q[2];
sx q[2];
rz(-1.2969742) q[2];
sx q[2];
rz(-1.3783) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4491475) q[1];
sx q[1];
rz(-1.0303823) q[1];
sx q[1];
rz(1.6442597) q[1];
rz(-pi) q[2];
rz(0.59806268) q[3];
sx q[3];
rz(-1.4334213) q[3];
sx q[3];
rz(-0.029702317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0868316) q[2];
sx q[2];
rz(-2.0589477) q[2];
sx q[2];
rz(-1.5209939) q[2];
rz(1.3711551) q[3];
sx q[3];
rz(-2.3833279) q[3];
sx q[3];
rz(-0.41485205) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3138251) q[0];
sx q[0];
rz(-2.1791552) q[0];
sx q[0];
rz(2.9058822) q[0];
rz(-0.57304263) q[1];
sx q[1];
rz(-1.0083116) q[1];
sx q[1];
rz(-1.3689573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1397676) q[0];
sx q[0];
rz(-0.84472504) q[0];
sx q[0];
rz(0.98066179) q[0];
rz(2.4339381) q[2];
sx q[2];
rz(-0.81205149) q[2];
sx q[2];
rz(2.6132513) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2409093) q[1];
sx q[1];
rz(-0.59421173) q[1];
sx q[1];
rz(1.2628984) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7810516) q[3];
sx q[3];
rz(-1.6586132) q[3];
sx q[3];
rz(0.58780625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3706751) q[2];
sx q[2];
rz(-1.3266027) q[2];
sx q[2];
rz(-3.0885922) q[2];
rz(-0.75183374) q[3];
sx q[3];
rz(-1.0337318) q[3];
sx q[3];
rz(2.2161765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62405217) q[0];
sx q[0];
rz(-1.0160099) q[0];
sx q[0];
rz(2.1852063) q[0];
rz(-1.7908295) q[1];
sx q[1];
rz(-0.39719926) q[1];
sx q[1];
rz(-0.15695708) q[1];
rz(1.440329) q[2];
sx q[2];
rz(-2.2638276) q[2];
sx q[2];
rz(2.9439415) q[2];
rz(2.1696321) q[3];
sx q[3];
rz(-0.85799118) q[3];
sx q[3];
rz(0.21256577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

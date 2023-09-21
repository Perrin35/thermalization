OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82436615) q[0];
sx q[0];
rz(5.1685652) q[0];
sx q[0];
rz(9.4246372) q[0];
rz(1.3340985) q[1];
sx q[1];
rz(4.1058022) q[1];
sx q[1];
rz(10.618187) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50085917) q[0];
sx q[0];
rz(-1.2841604) q[0];
sx q[0];
rz(0.1851693) q[0];
x q[1];
rz(-0.5483746) q[2];
sx q[2];
rz(-1.8273241) q[2];
sx q[2];
rz(-2.246849) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6686033) q[1];
sx q[1];
rz(-2.0358634) q[1];
sx q[1];
rz(2.4238062) q[1];
rz(-1.1080997) q[3];
sx q[3];
rz(-2.8994312) q[3];
sx q[3];
rz(2.7778181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6821735) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(-1.2288644) q[2];
rz(1.7284283) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(1.4878954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380149) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(-2.1287825) q[0];
rz(3.1139328) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(1.123463) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7682122) q[0];
sx q[0];
rz(-2.0313615) q[0];
sx q[0];
rz(-0.064231355) q[0];
x q[1];
rz(0.77983071) q[2];
sx q[2];
rz(-1.5260328) q[2];
sx q[2];
rz(1.1080527) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1883205) q[1];
sx q[1];
rz(-2.3768432) q[1];
sx q[1];
rz(-2.3482167) q[1];
x q[2];
rz(-2.4554159) q[3];
sx q[3];
rz(-2.3585329) q[3];
sx q[3];
rz(0.99883294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79364395) q[2];
sx q[2];
rz(-2.0517893) q[2];
sx q[2];
rz(-2.222555) q[2];
rz(-2.4675026) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(-1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27750257) q[0];
sx q[0];
rz(-2.9798177) q[0];
sx q[0];
rz(1.8664237) q[0];
rz(2.4480942) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(-2.0085874) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5223761) q[0];
sx q[0];
rz(-1.4285061) q[0];
sx q[0];
rz(3.1382568) q[0];
rz(-pi) q[1];
rz(-2.0793545) q[2];
sx q[2];
rz(-2.2937751) q[2];
sx q[2];
rz(-1.522097) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9085711) q[1];
sx q[1];
rz(-2.1070707) q[1];
sx q[1];
rz(2.8051393) q[1];
x q[2];
rz(-2.1902309) q[3];
sx q[3];
rz(-1.0361443) q[3];
sx q[3];
rz(-1.149328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(0.039316468) q[3];
sx q[3];
rz(-1.9226363) q[3];
sx q[3];
rz(-1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2599729) q[0];
sx q[0];
rz(-0.078475229) q[0];
sx q[0];
rz(1.1608634) q[0];
rz(2.2456031) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(-3.0060351) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7134705) q[0];
sx q[0];
rz(-0.8849511) q[0];
sx q[0];
rz(1.9786406) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9549471) q[2];
sx q[2];
rz(-1.5015366) q[2];
sx q[2];
rz(2.4406976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1783274) q[1];
sx q[1];
rz(-1.7566924) q[1];
sx q[1];
rz(-0.16190298) q[1];
rz(-0.1578219) q[3];
sx q[3];
rz(-1.8421679) q[3];
sx q[3];
rz(1.8872758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.23665145) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(-2.2616852) q[2];
rz(3.0974292) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(0.28863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0376461) q[0];
sx q[0];
rz(-0.3750616) q[0];
sx q[0];
rz(-2.1283545) q[0];
rz(-0.049731072) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(2.0577046) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3736553) q[0];
sx q[0];
rz(-0.3070139) q[0];
sx q[0];
rz(-2.2178749) q[0];
rz(2.3124144) q[2];
sx q[2];
rz(-0.36311705) q[2];
sx q[2];
rz(1.6104289) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0368082) q[1];
sx q[1];
rz(-2.0953062) q[1];
sx q[1];
rz(0.12496897) q[1];
rz(-pi) q[2];
x q[2];
rz(0.089133457) q[3];
sx q[3];
rz(-2.6260178) q[3];
sx q[3];
rz(2.6775132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23285398) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(2.8971635) q[2];
rz(-0.43236732) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(-0.50306815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2844834) q[0];
sx q[0];
rz(-1.4215707) q[0];
sx q[0];
rz(-0.094141468) q[0];
rz(2.969818) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(0.89541268) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0462414) q[0];
sx q[0];
rz(-1.6086676) q[0];
sx q[0];
rz(0.34237679) q[0];
rz(1.3307829) q[2];
sx q[2];
rz(-2.1905106) q[2];
sx q[2];
rz(1.0735219) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9858866) q[1];
sx q[1];
rz(-2.967318) q[1];
sx q[1];
rz(1.8272912) q[1];
rz(1.6454837) q[3];
sx q[3];
rz(-1.5684621) q[3];
sx q[3];
rz(2.1722349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0078997) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(0.80319476) q[2];
rz(-1.1903654) q[3];
sx q[3];
rz(-1.9093711) q[3];
sx q[3];
rz(-0.41263321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0728834) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(2.6224526) q[0];
rz(2.5601162) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(1.2566459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4706659) q[0];
sx q[0];
rz(-1.5304655) q[0];
sx q[0];
rz(1.8223902) q[0];
rz(-pi) q[1];
x q[1];
rz(1.16876) q[2];
sx q[2];
rz(-0.49905825) q[2];
sx q[2];
rz(-1.7390342) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1756806) q[1];
sx q[1];
rz(-1.6351846) q[1];
sx q[1];
rz(-2.8596911) q[1];
rz(-pi) q[2];
rz(-1.890896) q[3];
sx q[3];
rz(-0.88722908) q[3];
sx q[3];
rz(0.86910955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98823035) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(-1.777565) q[2];
rz(2.2310232) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751188) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(2.8334154) q[0];
rz(-0.072487436) q[1];
sx q[1];
rz(-1.0132353) q[1];
sx q[1];
rz(0.38696188) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68650866) q[0];
sx q[0];
rz(-2.5065656) q[0];
sx q[0];
rz(-1.1839266) q[0];
rz(-pi) q[1];
rz(-0.87256356) q[2];
sx q[2];
rz(-1.1485032) q[2];
sx q[2];
rz(-0.97908212) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8291694) q[1];
sx q[1];
rz(-2.2409391) q[1];
sx q[1];
rz(-1.4751242) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3513226) q[3];
sx q[3];
rz(-2.05728) q[3];
sx q[3];
rz(-0.51318491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5179634) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(0.44000885) q[2];
rz(-0.7157588) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(-2.0429042) q[0];
rz(0.72775841) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(3.0922906) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54287275) q[0];
sx q[0];
rz(-0.96519404) q[0];
sx q[0];
rz(-2.4332895) q[0];
rz(0.98408913) q[2];
sx q[2];
rz(-2.3684635) q[2];
sx q[2];
rz(0.099260515) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94433632) q[1];
sx q[1];
rz(-1.1068871) q[1];
sx q[1];
rz(-2.0891561) q[1];
rz(2.6155869) q[3];
sx q[3];
rz(-0.81323871) q[3];
sx q[3];
rz(-0.40532743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0631642) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(2.4196529) q[2];
rz(2.1980964) q[3];
sx q[3];
rz(-0.75073457) q[3];
sx q[3];
rz(-0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14770517) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(-2.06185) q[0];
rz(-2.0823157) q[1];
sx q[1];
rz(-2.9187027) q[1];
sx q[1];
rz(-1.7396897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7125177) q[0];
sx q[0];
rz(-1.4445462) q[0];
sx q[0];
rz(-1.7849126) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9981668) q[2];
sx q[2];
rz(-0.18866814) q[2];
sx q[2];
rz(-0.27486899) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68114963) q[1];
sx q[1];
rz(-2.3773758) q[1];
sx q[1];
rz(-1.1230254) q[1];
x q[2];
rz(0.63125061) q[3];
sx q[3];
rz(-1.7389986) q[3];
sx q[3];
rz(1.9025308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.4828651) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(-1.6213017) q[2];
rz(-2.5907717) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(2.4441392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14810066) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(0.91611721) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(-2.312071) q[2];
sx q[2];
rz(-2.1567232) q[2];
sx q[2];
rz(0.17270252) q[2];
rz(-1.5394474) q[3];
sx q[3];
rz(-2.6242704) q[3];
sx q[3];
rz(-0.27192413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
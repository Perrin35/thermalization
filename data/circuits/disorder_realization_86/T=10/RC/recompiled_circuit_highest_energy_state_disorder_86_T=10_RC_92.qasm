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
rz(-1.2404233) q[0];
sx q[0];
rz(1.9440396) q[0];
sx q[0];
rz(9.64111) q[0];
rz(4.0989838) q[1];
sx q[1];
rz(5.6070072) q[1];
sx q[1];
rz(11.054872) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0759685) q[0];
sx q[0];
rz(-1.0145717) q[0];
sx q[0];
rz(-1.5594622) q[0];
rz(-0.47343238) q[2];
sx q[2];
rz(-1.8634999) q[2];
sx q[2];
rz(0.44854376) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2701157) q[1];
sx q[1];
rz(-1.9028712) q[1];
sx q[1];
rz(0.52712743) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0262296) q[3];
sx q[3];
rz(-2.3854227) q[3];
sx q[3];
rz(-2.2952473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.034885255) q[2];
sx q[2];
rz(-0.35197508) q[2];
sx q[2];
rz(2.0931639) q[2];
rz(2.9599221) q[3];
sx q[3];
rz(-0.964966) q[3];
sx q[3];
rz(2.488193) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.492391) q[0];
sx q[0];
rz(-0.94434706) q[0];
sx q[0];
rz(-2.6948068) q[0];
rz(-0.88042879) q[1];
sx q[1];
rz(-1.3648405) q[1];
sx q[1];
rz(2.3588038) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44109684) q[0];
sx q[0];
rz(-1.6320328) q[0];
sx q[0];
rz(1.815531) q[0];
rz(-1.9982463) q[2];
sx q[2];
rz(-2.1437217) q[2];
sx q[2];
rz(-1.7597511) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70322733) q[1];
sx q[1];
rz(-1.2709054) q[1];
sx q[1];
rz(-0.48734003) q[1];
x q[2];
rz(-1.0083593) q[3];
sx q[3];
rz(-0.76518067) q[3];
sx q[3];
rz(-0.4972813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.11144) q[2];
sx q[2];
rz(-1.4806662) q[2];
sx q[2];
rz(2.1622369) q[2];
rz(0.3848981) q[3];
sx q[3];
rz(-1.220547) q[3];
sx q[3];
rz(-0.16429193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6556743) q[0];
sx q[0];
rz(-0.05412183) q[0];
sx q[0];
rz(-2.3517877) q[0];
rz(-2.9549331) q[1];
sx q[1];
rz(-1.7402382) q[1];
sx q[1];
rz(2.1479215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45897608) q[0];
sx q[0];
rz(-2.1135215) q[0];
sx q[0];
rz(1.8629462) q[0];
rz(-0.31618677) q[2];
sx q[2];
rz(-0.63717604) q[2];
sx q[2];
rz(0.030046163) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51998617) q[1];
sx q[1];
rz(-2.4944759) q[1];
sx q[1];
rz(-2.3220329) q[1];
rz(-pi) q[2];
x q[2];
rz(2.42071) q[3];
sx q[3];
rz(-1.3647121) q[3];
sx q[3];
rz(2.4997847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8197202) q[2];
sx q[2];
rz(-1.1178144) q[2];
sx q[2];
rz(-2.7703088) q[2];
rz(0.34879455) q[3];
sx q[3];
rz(-2.0472417) q[3];
sx q[3];
rz(-0.5955407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6868941) q[0];
sx q[0];
rz(-2.1214387) q[0];
sx q[0];
rz(-1.4042847) q[0];
rz(0.67717254) q[1];
sx q[1];
rz(-1.155747) q[1];
sx q[1];
rz(1.4926532) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66144511) q[0];
sx q[0];
rz(-0.3269402) q[0];
sx q[0];
rz(-0.72823712) q[0];
rz(-pi) q[1];
rz(-0.29784305) q[2];
sx q[2];
rz(-1.1561511) q[2];
sx q[2];
rz(-0.21098247) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.75912913) q[1];
sx q[1];
rz(-1.8378403) q[1];
sx q[1];
rz(1.8014531) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31217137) q[3];
sx q[3];
rz(-0.87066423) q[3];
sx q[3];
rz(-2.2362102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6877785) q[2];
sx q[2];
rz(-0.9684338) q[2];
sx q[2];
rz(-0.34823927) q[2];
rz(-1.6866775) q[3];
sx q[3];
rz(-1.7125407) q[3];
sx q[3];
rz(-1.3454364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1763879) q[0];
sx q[0];
rz(-2.5768953) q[0];
sx q[0];
rz(0.89214605) q[0];
rz(2.6773894) q[1];
sx q[1];
rz(-1.2449539) q[1];
sx q[1];
rz(1.2947327) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1071467) q[0];
sx q[0];
rz(-2.4305516) q[0];
sx q[0];
rz(-0.76816316) q[0];
rz(-pi) q[1];
rz(-1.6991529) q[2];
sx q[2];
rz(-1.0367298) q[2];
sx q[2];
rz(-2.4160224) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6695822) q[1];
sx q[1];
rz(-1.8947487) q[1];
sx q[1];
rz(-1.1060017) q[1];
rz(-pi) q[2];
rz(-0.89035676) q[3];
sx q[3];
rz(-2.0687752) q[3];
sx q[3];
rz(-2.339956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3225473) q[2];
sx q[2];
rz(-0.61083856) q[2];
sx q[2];
rz(0.56582212) q[2];
rz(0.081341751) q[3];
sx q[3];
rz(-0.9809202) q[3];
sx q[3];
rz(-2.5210023) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.464798) q[0];
sx q[0];
rz(-3.0749574) q[0];
sx q[0];
rz(-1.5555405) q[0];
rz(1.0606891) q[1];
sx q[1];
rz(-1.565275) q[1];
sx q[1];
rz(2.5097844) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6752967) q[0];
sx q[0];
rz(-1.491437) q[0];
sx q[0];
rz(-2.9064889) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2802629) q[2];
sx q[2];
rz(-1.1785186) q[2];
sx q[2];
rz(2.5226468) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7673847) q[1];
sx q[1];
rz(-2.4904618) q[1];
sx q[1];
rz(2.1506449) q[1];
rz(-0.68354179) q[3];
sx q[3];
rz(-1.6850796) q[3];
sx q[3];
rz(1.6268693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.98171988) q[2];
sx q[2];
rz(-2.0231415) q[2];
sx q[2];
rz(2.6427606) q[2];
rz(-1.8286797) q[3];
sx q[3];
rz(-2.4098318) q[3];
sx q[3];
rz(1.7074728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.6071534) q[0];
sx q[0];
rz(-2.7767015) q[0];
sx q[0];
rz(0.69951192) q[0];
rz(0.46547678) q[1];
sx q[1];
rz(-0.87124467) q[1];
sx q[1];
rz(2.2043998) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58317157) q[0];
sx q[0];
rz(-1.7390378) q[0];
sx q[0];
rz(2.9267163) q[0];
rz(-2.7363051) q[2];
sx q[2];
rz(-2.8817085) q[2];
sx q[2];
rz(-2.0854307) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4129137) q[1];
sx q[1];
rz(-1.2866486) q[1];
sx q[1];
rz(-0.22032622) q[1];
rz(-2.8392467) q[3];
sx q[3];
rz(-1.8519326) q[3];
sx q[3];
rz(-0.1880364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.84404868) q[2];
sx q[2];
rz(-1.8397477) q[2];
sx q[2];
rz(0.37332264) q[2];
rz(-1.1897872) q[3];
sx q[3];
rz(-2.6326284) q[3];
sx q[3];
rz(0.51026195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3968286) q[0];
sx q[0];
rz(-1.0823534) q[0];
sx q[0];
rz(-0.74791351) q[0];
rz(-2.3751936) q[1];
sx q[1];
rz(-2.8728569) q[1];
sx q[1];
rz(0.0029729923) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82389861) q[0];
sx q[0];
rz(-2.4269773) q[0];
sx q[0];
rz(1.8665642) q[0];
rz(-1.6417129) q[2];
sx q[2];
rz(-1.5512964) q[2];
sx q[2];
rz(1.8766581) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5103448) q[1];
sx q[1];
rz(-1.5273792) q[1];
sx q[1];
rz(-1.7249291) q[1];
rz(-pi) q[2];
x q[2];
rz(2.946704) q[3];
sx q[3];
rz(-1.7656544) q[3];
sx q[3];
rz(-0.080435924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.11091867) q[2];
sx q[2];
rz(-1.7577533) q[2];
sx q[2];
rz(-2.0745011) q[2];
rz(0.083960697) q[3];
sx q[3];
rz(-0.48560086) q[3];
sx q[3];
rz(2.3626204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79754168) q[0];
sx q[0];
rz(-2.2028956) q[0];
sx q[0];
rz(-3.0352266) q[0];
rz(-0.97995177) q[1];
sx q[1];
rz(-1.4915024) q[1];
sx q[1];
rz(-2.3756557) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.879012) q[0];
sx q[0];
rz(-1.1404788) q[0];
sx q[0];
rz(2.5613214) q[0];
x q[1];
rz(-2.7226549) q[2];
sx q[2];
rz(-2.9235002) q[2];
sx q[2];
rz(0.60148394) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.83726604) q[1];
sx q[1];
rz(-0.89511739) q[1];
sx q[1];
rz(-1.6153687) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8838653) q[3];
sx q[3];
rz(-0.014583909) q[3];
sx q[3];
rz(-2.7590883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.47436675) q[2];
sx q[2];
rz(-2.6211278) q[2];
sx q[2];
rz(0.70029798) q[2];
rz(-2.4750366) q[3];
sx q[3];
rz(-1.8066112) q[3];
sx q[3];
rz(-2.0535645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(2.4661082) q[0];
sx q[0];
rz(-2.1145144) q[0];
sx q[0];
rz(1.3960557) q[0];
rz(-1.667977) q[1];
sx q[1];
rz(-1.8831848) q[1];
sx q[1];
rz(0.6667164) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4748866) q[0];
sx q[0];
rz(-0.96893725) q[0];
sx q[0];
rz(-0.78420297) q[0];
rz(-pi) q[1];
rz(-1.855895) q[2];
sx q[2];
rz(-0.93109967) q[2];
sx q[2];
rz(2.6486047) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0565225) q[1];
sx q[1];
rz(-1.4398972) q[1];
sx q[1];
rz(-0.31899778) q[1];
x q[2];
rz(2.1176841) q[3];
sx q[3];
rz(-1.0677665) q[3];
sx q[3];
rz(-1.5707113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0746158) q[2];
sx q[2];
rz(-2.4846027) q[2];
sx q[2];
rz(2.2508049) q[2];
rz(2.7095419) q[3];
sx q[3];
rz(-1.0703577) q[3];
sx q[3];
rz(0.74705684) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4074832) q[0];
sx q[0];
rz(-1.643184) q[0];
sx q[0];
rz(-1.2947422) q[0];
rz(-1.8941849) q[1];
sx q[1];
rz(-0.71949646) q[1];
sx q[1];
rz(1.234642) q[1];
rz(-2.9802889) q[2];
sx q[2];
rz(-1.9111173) q[2];
sx q[2];
rz(2.5572122) q[2];
rz(3.0790764) q[3];
sx q[3];
rz(-2.6517031) q[3];
sx q[3];
rz(-0.22417886) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

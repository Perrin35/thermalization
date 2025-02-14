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
rz(-0.15859088) q[0];
sx q[0];
rz(2.6245485) q[0];
sx q[0];
rz(11.256097) q[0];
rz(-0.34215555) q[1];
sx q[1];
rz(-1.181239) q[1];
sx q[1];
rz(0.96460834) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2268596) q[0];
sx q[0];
rz(-0.64382416) q[0];
sx q[0];
rz(-0.22636803) q[0];
rz(-pi) q[1];
rz(1.3753424) q[2];
sx q[2];
rz(-1.9103622) q[2];
sx q[2];
rz(1.5528785) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8339616) q[1];
sx q[1];
rz(-1.2738859) q[1];
sx q[1];
rz(2.3351257) q[1];
rz(0.11537376) q[3];
sx q[3];
rz(-2.2386207) q[3];
sx q[3];
rz(1.7911719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1356807) q[2];
sx q[2];
rz(-0.19401208) q[2];
sx q[2];
rz(1.676959) q[2];
rz(1.3095193) q[3];
sx q[3];
rz(-2.194761) q[3];
sx q[3];
rz(2.8215698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2784073) q[0];
sx q[0];
rz(-2.0018556) q[0];
sx q[0];
rz(-1.7742668) q[0];
rz(1.4890081) q[1];
sx q[1];
rz(-2.3431578) q[1];
sx q[1];
rz(2.8489825) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83787888) q[0];
sx q[0];
rz(-0.3488763) q[0];
sx q[0];
rz(0.29033355) q[0];
x q[1];
rz(0.16186773) q[2];
sx q[2];
rz(-1.6965995) q[2];
sx q[2];
rz(2.6243072) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83963385) q[1];
sx q[1];
rz(-1.3657161) q[1];
sx q[1];
rz(2.9006216) q[1];
rz(-pi) q[2];
rz(-1.5047362) q[3];
sx q[3];
rz(-1.8365897) q[3];
sx q[3];
rz(3.0817666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9802398) q[2];
sx q[2];
rz(-0.62675256) q[2];
sx q[2];
rz(0.18388595) q[2];
rz(0.45267496) q[3];
sx q[3];
rz(-1.5225182) q[3];
sx q[3];
rz(-3.0116853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(1.1995131) q[0];
sx q[0];
rz(-0.71490723) q[0];
sx q[0];
rz(-2.3364501) q[0];
rz(1.0792271) q[1];
sx q[1];
rz(-2.3715623) q[1];
sx q[1];
rz(-1.315518) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0502687) q[0];
sx q[0];
rz(-1.5116232) q[0];
sx q[0];
rz(2.4135804) q[0];
x q[1];
rz(-0.82683021) q[2];
sx q[2];
rz(-1.2203802) q[2];
sx q[2];
rz(-2.0598799) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6650369) q[1];
sx q[1];
rz(-1.2893234) q[1];
sx q[1];
rz(-3.046077) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27045111) q[3];
sx q[3];
rz(-2.8282249) q[3];
sx q[3];
rz(-0.53821401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9341854) q[2];
sx q[2];
rz(-1.4855874) q[2];
sx q[2];
rz(2.5918813) q[2];
rz(-3.1304729) q[3];
sx q[3];
rz(-2.34237) q[3];
sx q[3];
rz(1.3219249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2928807) q[0];
sx q[0];
rz(-0.95893812) q[0];
sx q[0];
rz(0.30353656) q[0];
rz(-0.15829463) q[1];
sx q[1];
rz(-1.0946495) q[1];
sx q[1];
rz(-2.1319481) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9859814) q[0];
sx q[0];
rz(-1.2878875) q[0];
sx q[0];
rz(-1.1213746) q[0];
rz(-0.43731205) q[2];
sx q[2];
rz(-0.6491937) q[2];
sx q[2];
rz(-1.5199666) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4588488) q[1];
sx q[1];
rz(-1.7284596) q[1];
sx q[1];
rz(0.47994061) q[1];
rz(-pi) q[2];
x q[2];
rz(0.03831717) q[3];
sx q[3];
rz(-2.6434757) q[3];
sx q[3];
rz(-1.8565053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2046854) q[2];
sx q[2];
rz(-2.3202809) q[2];
sx q[2];
rz(2.8724907) q[2];
rz(-1.4540295) q[3];
sx q[3];
rz(-1.5498091) q[3];
sx q[3];
rz(-1.5788186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.27075416) q[0];
sx q[0];
rz(-2.0814867) q[0];
sx q[0];
rz(0.36218542) q[0];
rz(-2.4902792) q[1];
sx q[1];
rz(-1.6505227) q[1];
sx q[1];
rz(-1.6168894) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7633879) q[0];
sx q[0];
rz(-1.7001061) q[0];
sx q[0];
rz(-0.61152258) q[0];
rz(-pi) q[1];
rz(-2.1936622) q[2];
sx q[2];
rz(-1.6375223) q[2];
sx q[2];
rz(-0.78363505) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0901047) q[1];
sx q[1];
rz(-1.0516775) q[1];
sx q[1];
rz(-2.725802) q[1];
rz(-0.42585856) q[3];
sx q[3];
rz(-2.3506864) q[3];
sx q[3];
rz(0.87825852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0321956) q[2];
sx q[2];
rz(-2.1869662) q[2];
sx q[2];
rz(1.7677914) q[2];
rz(-2.7023756) q[3];
sx q[3];
rz(-0.64911157) q[3];
sx q[3];
rz(-2.5325328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1010901) q[0];
sx q[0];
rz(-2.4093565) q[0];
sx q[0];
rz(2.8771583) q[0];
rz(2.4624372) q[1];
sx q[1];
rz(-0.87398386) q[1];
sx q[1];
rz(-0.98495475) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0593582) q[0];
sx q[0];
rz(-1.2684039) q[0];
sx q[0];
rz(-1.9134269) q[0];
rz(-2.1169871) q[2];
sx q[2];
rz(-2.6272079) q[2];
sx q[2];
rz(2.3895604) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90675844) q[1];
sx q[1];
rz(-2.0528767) q[1];
sx q[1];
rz(-1.158403) q[1];
x q[2];
rz(1.156928) q[3];
sx q[3];
rz(-2.664251) q[3];
sx q[3];
rz(-1.1787924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1017477) q[2];
sx q[2];
rz(-0.16182772) q[2];
sx q[2];
rz(2.7092773) q[2];
rz(2.1281706) q[3];
sx q[3];
rz(-1.6067959) q[3];
sx q[3];
rz(-2.5920946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48563114) q[0];
sx q[0];
rz(-0.39325842) q[0];
sx q[0];
rz(1.1095169) q[0];
rz(2.3484777) q[1];
sx q[1];
rz(-1.3536645) q[1];
sx q[1];
rz(1.5951593) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6731747) q[0];
sx q[0];
rz(-1.0145079) q[0];
sx q[0];
rz(-2.2493717) q[0];
rz(-pi) q[1];
rz(-1.3605453) q[2];
sx q[2];
rz(-2.2814007) q[2];
sx q[2];
rz(2.9637314) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7132414) q[1];
sx q[1];
rz(-0.71231406) q[1];
sx q[1];
rz(1.8800432) q[1];
x q[2];
rz(3.0199354) q[3];
sx q[3];
rz(-2.4215048) q[3];
sx q[3];
rz(-0.64313221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.82211295) q[2];
sx q[2];
rz(-0.18222465) q[2];
sx q[2];
rz(1.3694084) q[2];
rz(-0.93044126) q[3];
sx q[3];
rz(-1.3481827) q[3];
sx q[3];
rz(1.5038917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1996138) q[0];
sx q[0];
rz(-1.9604585) q[0];
sx q[0];
rz(-0.30366316) q[0];
rz(-2.1619201) q[1];
sx q[1];
rz(-0.62398463) q[1];
sx q[1];
rz(2.7365275) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4265012) q[0];
sx q[0];
rz(-0.42761567) q[0];
sx q[0];
rz(-2.6602488) q[0];
rz(-pi) q[1];
rz(-1.1978537) q[2];
sx q[2];
rz(-0.44633807) q[2];
sx q[2];
rz(0.17910236) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0822337) q[1];
sx q[1];
rz(-2.5463366) q[1];
sx q[1];
rz(-2.8221376) q[1];
x q[2];
rz(1.9937421) q[3];
sx q[3];
rz(-2.6788524) q[3];
sx q[3];
rz(0.014960814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8651809) q[2];
sx q[2];
rz(-0.96651912) q[2];
sx q[2];
rz(-0.69918862) q[2];
rz(0.80896038) q[3];
sx q[3];
rz(-1.8374846) q[3];
sx q[3];
rz(-2.8185524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.9998099) q[0];
sx q[0];
rz(-1.9843822) q[0];
sx q[0];
rz(-3.004177) q[0];
rz(0.049086463) q[1];
sx q[1];
rz(-2.0259435) q[1];
sx q[1];
rz(2.633599) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.346283) q[0];
sx q[0];
rz(-0.49581832) q[0];
sx q[0];
rz(3.0853372) q[0];
rz(-pi) q[1];
rz(1.5586583) q[2];
sx q[2];
rz(-1.7430787) q[2];
sx q[2];
rz(0.42086312) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55785364) q[1];
sx q[1];
rz(-2.5370829) q[1];
sx q[1];
rz(0.97911759) q[1];
x q[2];
rz(0.93799641) q[3];
sx q[3];
rz(-2.9909903) q[3];
sx q[3];
rz(-1.6737398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.55234838) q[2];
sx q[2];
rz(-1.9231223) q[2];
sx q[2];
rz(-0.98769665) q[2];
rz(-2.6622631) q[3];
sx q[3];
rz(-1.597581) q[3];
sx q[3];
rz(-2.2492669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4722897) q[0];
sx q[0];
rz(-0.23858128) q[0];
sx q[0];
rz(-0.14078374) q[0];
rz(-2.7760778) q[1];
sx q[1];
rz(-0.82217685) q[1];
sx q[1];
rz(-0.64868322) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1876342) q[0];
sx q[0];
rz(-1.9184858) q[0];
sx q[0];
rz(2.4566922) q[0];
rz(-pi) q[1];
rz(0.76463215) q[2];
sx q[2];
rz(-1.6551541) q[2];
sx q[2];
rz(0.81071883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5488447) q[1];
sx q[1];
rz(-2.5429568) q[1];
sx q[1];
rz(3.0416802) q[1];
rz(-pi) q[2];
rz(-2.5124536) q[3];
sx q[3];
rz(-0.39689343) q[3];
sx q[3];
rz(1.6861738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4623798) q[2];
sx q[2];
rz(-1.0054192) q[2];
sx q[2];
rz(3.0688378) q[2];
rz(-2.4788729) q[3];
sx q[3];
rz(-1.3136256) q[3];
sx q[3];
rz(-2.976118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0872021) q[0];
sx q[0];
rz(-1.5329755) q[0];
sx q[0];
rz(-1.6736915) q[0];
rz(1.5351334) q[1];
sx q[1];
rz(-0.88773334) q[1];
sx q[1];
rz(1.6233374) q[1];
rz(1.4521815) q[2];
sx q[2];
rz(-1.7977503) q[2];
sx q[2];
rz(1.1888421) q[2];
rz(2.2519464) q[3];
sx q[3];
rz(-1.8915265) q[3];
sx q[3];
rz(-1.693207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

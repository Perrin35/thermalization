OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39188448) q[0];
sx q[0];
rz(-0.19667974) q[0];
sx q[0];
rz(-1.952202) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(2.300188) q[1];
sx q[1];
rz(9.0471164) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1188287) q[0];
sx q[0];
rz(-1.6435992) q[0];
sx q[0];
rz(1.4490436) q[0];
x q[1];
rz(-1.6701277) q[2];
sx q[2];
rz(-2.8566395) q[2];
sx q[2];
rz(-0.003665912) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.79757534) q[1];
sx q[1];
rz(-2.0098149) q[1];
sx q[1];
rz(-2.3637191) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7226726) q[3];
sx q[3];
rz(-1.4421316) q[3];
sx q[3];
rz(-3.0959689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.010014023) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(0.67260355) q[2];
rz(-0.16942313) q[3];
sx q[3];
rz(-2.7526581) q[3];
sx q[3];
rz(1.3385564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.910903) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(-2.9130274) q[0];
rz(-2.9810492) q[1];
sx q[1];
rz(-1.4385782) q[1];
sx q[1];
rz(-0.28796089) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1624958) q[0];
sx q[0];
rz(-1.393035) q[0];
sx q[0];
rz(1.7642679) q[0];
rz(-2.1638853) q[2];
sx q[2];
rz(-1.8988673) q[2];
sx q[2];
rz(-1.6239945) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.49111734) q[1];
sx q[1];
rz(-0.75956356) q[1];
sx q[1];
rz(-2.526545) q[1];
x q[2];
rz(-2.7912451) q[3];
sx q[3];
rz(-1.5347893) q[3];
sx q[3];
rz(-2.2909174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5471389) q[2];
sx q[2];
rz(-1.889166) q[2];
sx q[2];
rz(-1.4734369) q[2];
rz(2.2041221) q[3];
sx q[3];
rz(-0.44527403) q[3];
sx q[3];
rz(0.37500769) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32524747) q[0];
sx q[0];
rz(-0.49961093) q[0];
sx q[0];
rz(-0.26741272) q[0];
rz(1.422241) q[1];
sx q[1];
rz(-2.0098675) q[1];
sx q[1];
rz(-0.95169383) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0824453) q[0];
sx q[0];
rz(-1.379231) q[0];
sx q[0];
rz(-3.0491769) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0864262) q[2];
sx q[2];
rz(-1.3736758) q[2];
sx q[2];
rz(-0.65161639) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.051085) q[1];
sx q[1];
rz(-1.8790434) q[1];
sx q[1];
rz(0.80866637) q[1];
rz(-pi) q[2];
rz(0.94727256) q[3];
sx q[3];
rz(-1.1215278) q[3];
sx q[3];
rz(2.505213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0114228) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(-0.34417957) q[2];
rz(2.2551645) q[3];
sx q[3];
rz(-0.27799806) q[3];
sx q[3];
rz(2.5755431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0043871) q[0];
sx q[0];
rz(-1.6442278) q[0];
sx q[0];
rz(-1.8970867) q[0];
rz(-3.0124774) q[1];
sx q[1];
rz(-1.7872417) q[1];
sx q[1];
rz(-2.7688162) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9739649) q[0];
sx q[0];
rz(-2.3229257) q[0];
sx q[0];
rz(1.7276093) q[0];
rz(-pi) q[1];
rz(-0.77809019) q[2];
sx q[2];
rz(-1.8847244) q[2];
sx q[2];
rz(-2.2997466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2671632) q[1];
sx q[1];
rz(-2.0544555) q[1];
sx q[1];
rz(2.8664385) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9863015) q[3];
sx q[3];
rz(-0.55404514) q[3];
sx q[3];
rz(-2.154175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.63561511) q[2];
sx q[2];
rz(-1.4900692) q[2];
sx q[2];
rz(0.90488952) q[2];
rz(-2.7010226) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(0.99159616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6624517) q[0];
sx q[0];
rz(-0.15563706) q[0];
sx q[0];
rz(-1.8537846) q[0];
rz(-1.4783391) q[1];
sx q[1];
rz(-1.9961424) q[1];
sx q[1];
rz(2.6073661) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6144845) q[0];
sx q[0];
rz(-1.8728349) q[0];
sx q[0];
rz(-1.8943647) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0787557) q[2];
sx q[2];
rz(-1.8740219) q[2];
sx q[2];
rz(-2.6360896) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2248762) q[1];
sx q[1];
rz(-1.4955048) q[1];
sx q[1];
rz(0.063898357) q[1];
rz(-0.71877919) q[3];
sx q[3];
rz(-1.4712417) q[3];
sx q[3];
rz(-2.0841141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6405032) q[2];
sx q[2];
rz(-1.2155632) q[2];
sx q[2];
rz(-2.9980998) q[2];
rz(1.7701373) q[3];
sx q[3];
rz(-1.2797132) q[3];
sx q[3];
rz(-2.9218856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4390398) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(-2.904536) q[0];
rz(-1.9006231) q[1];
sx q[1];
rz(-2.3126912) q[1];
sx q[1];
rz(-1.3751078) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28778827) q[0];
sx q[0];
rz(-1.8332991) q[0];
sx q[0];
rz(2.7909423) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3926922) q[2];
sx q[2];
rz(-1.5178711) q[2];
sx q[2];
rz(0.63268328) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4543537) q[1];
sx q[1];
rz(-0.41077405) q[1];
sx q[1];
rz(-0.06128581) q[1];
rz(-pi) q[2];
rz(1.9938019) q[3];
sx q[3];
rz(-2.0241963) q[3];
sx q[3];
rz(-2.6231433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6087626) q[2];
sx q[2];
rz(-2.4217748) q[2];
sx q[2];
rz(2.9928845) q[2];
rz(-3.1249629) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(0.19255157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6454813) q[0];
sx q[0];
rz(-0.2922903) q[0];
sx q[0];
rz(-2.2684229) q[0];
rz(-2.219615) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(0.08392863) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.98691) q[0];
sx q[0];
rz(-1.0202408) q[0];
sx q[0];
rz(-2.8793392) q[0];
rz(-1.7882118) q[2];
sx q[2];
rz(-1.916421) q[2];
sx q[2];
rz(0.84920041) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2578656) q[1];
sx q[1];
rz(-2.1818433) q[1];
sx q[1];
rz(0.74531598) q[1];
rz(-pi) q[2];
x q[2];
rz(2.300802) q[3];
sx q[3];
rz(-1.519324) q[3];
sx q[3];
rz(-1.3703025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7579047) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(2.596358) q[2];
rz(-0.43045726) q[3];
sx q[3];
rz(-1.1377708) q[3];
sx q[3];
rz(0.30495131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6298744) q[0];
sx q[0];
rz(-3.1261303) q[0];
sx q[0];
rz(-1.2114725) q[0];
rz(-0.97310549) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(2.1957695) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4771172) q[0];
sx q[0];
rz(-2.9070589) q[0];
sx q[0];
rz(2.7995336) q[0];
rz(-pi) q[1];
rz(-0.68265712) q[2];
sx q[2];
rz(-2.1652512) q[2];
sx q[2];
rz(-0.97460954) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.91499) q[1];
sx q[1];
rz(-1.9493002) q[1];
sx q[1];
rz(-1.6016017) q[1];
rz(-pi) q[2];
rz(3.1315348) q[3];
sx q[3];
rz(-2.4754482) q[3];
sx q[3];
rz(-2.2705164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8699845) q[2];
sx q[2];
rz(-1.4166069) q[2];
sx q[2];
rz(2.8611709) q[2];
rz(-2.5366606) q[3];
sx q[3];
rz(-2.0623902) q[3];
sx q[3];
rz(2.159507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1542926) q[0];
sx q[0];
rz(-2.2315318) q[0];
sx q[0];
rz(-0.61532414) q[0];
rz(2.212021) q[1];
sx q[1];
rz(-2.2832182) q[1];
sx q[1];
rz(2.5659134) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8262779) q[0];
sx q[0];
rz(-1.1018254) q[0];
sx q[0];
rz(-2.4995575) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0909537) q[2];
sx q[2];
rz(-2.9165604) q[2];
sx q[2];
rz(-2.2573543) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.74734028) q[1];
sx q[1];
rz(-1.5963975) q[1];
sx q[1];
rz(-1.6781374) q[1];
rz(-pi) q[2];
rz(-2.1115033) q[3];
sx q[3];
rz(-1.4095777) q[3];
sx q[3];
rz(-0.89564039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3118887) q[2];
sx q[2];
rz(-2.3788033) q[2];
sx q[2];
rz(0.021961948) q[2];
rz(0.17045505) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(0.26486614) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72572529) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(2.5073994) q[0];
rz(-0.028907396) q[1];
sx q[1];
rz(-2.3529265) q[1];
sx q[1];
rz(-2.172487) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6432188) q[0];
sx q[0];
rz(-1.6341097) q[0];
sx q[0];
rz(-1.4988585) q[0];
x q[1];
rz(2.3949941) q[2];
sx q[2];
rz(-1.714163) q[2];
sx q[2];
rz(-1.5310841) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9173968) q[1];
sx q[1];
rz(-0.90985137) q[1];
sx q[1];
rz(-1.3437273) q[1];
rz(2.5892331) q[3];
sx q[3];
rz(-0.5061572) q[3];
sx q[3];
rz(1.6278704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6897631) q[2];
sx q[2];
rz(-1.4844866) q[2];
sx q[2];
rz(-2.7815212) q[2];
rz(-1.5277956) q[3];
sx q[3];
rz(-0.66643047) q[3];
sx q[3];
rz(-2.578919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2071028) q[0];
sx q[0];
rz(-1.5705382) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(3.0974401) q[1];
sx q[1];
rz(-1.4587198) q[1];
sx q[1];
rz(-1.1062467) q[1];
rz(1.1883696) q[2];
sx q[2];
rz(-1.8845176) q[2];
sx q[2];
rz(0.11558576) q[2];
rz(-2.5526657) q[3];
sx q[3];
rz(-0.52994655) q[3];
sx q[3];
rz(0.51701057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
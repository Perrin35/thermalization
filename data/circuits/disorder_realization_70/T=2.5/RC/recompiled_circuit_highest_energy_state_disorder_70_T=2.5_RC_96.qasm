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
rz(2.6707668) q[0];
sx q[0];
rz(-0.45007053) q[0];
sx q[0];
rz(-0.39029628) q[0];
rz(0.14416873) q[1];
sx q[1];
rz(-1.498797) q[1];
sx q[1];
rz(2.1013451) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24867985) q[0];
sx q[0];
rz(-1.8572766) q[0];
sx q[0];
rz(1.9353959) q[0];
x q[1];
rz(-2.7778012) q[2];
sx q[2];
rz(-1.7161088) q[2];
sx q[2];
rz(-0.34749588) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.39506787) q[1];
sx q[1];
rz(-1.7774425) q[1];
sx q[1];
rz(2.2499491) q[1];
x q[2];
rz(1.5055429) q[3];
sx q[3];
rz(-0.79354078) q[3];
sx q[3];
rz(-0.82543594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1990004) q[2];
sx q[2];
rz(-0.60574836) q[2];
sx q[2];
rz(-1.5462297) q[2];
rz(-0.47510251) q[3];
sx q[3];
rz(-0.64517704) q[3];
sx q[3];
rz(1.9861541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5652745) q[0];
sx q[0];
rz(-0.98868889) q[0];
sx q[0];
rz(-2.7069672) q[0];
rz(1.0644396) q[1];
sx q[1];
rz(-0.39355215) q[1];
sx q[1];
rz(-2.2501066) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28434174) q[0];
sx q[0];
rz(-1.9280757) q[0];
sx q[0];
rz(-0.42755677) q[0];
x q[1];
rz(1.0657677) q[2];
sx q[2];
rz(-2.2314851) q[2];
sx q[2];
rz(-2.6434174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6124208) q[1];
sx q[1];
rz(-2.0222221) q[1];
sx q[1];
rz(-0.58731095) q[1];
rz(-pi) q[2];
rz(-1.5323029) q[3];
sx q[3];
rz(-0.6596237) q[3];
sx q[3];
rz(1.4764113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.50856227) q[2];
sx q[2];
rz(-0.56266251) q[2];
sx q[2];
rz(-1.0594581) q[2];
rz(1.2262454) q[3];
sx q[3];
rz(-1.6112593) q[3];
sx q[3];
rz(-3.0219769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0934963) q[0];
sx q[0];
rz(-2.7432848) q[0];
sx q[0];
rz(-2.4523729) q[0];
rz(-2.7293909) q[1];
sx q[1];
rz(-1.1176132) q[1];
sx q[1];
rz(-1.162792) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6434518) q[0];
sx q[0];
rz(-1.8751133) q[0];
sx q[0];
rz(-2.9813719) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96647977) q[2];
sx q[2];
rz(-0.89077026) q[2];
sx q[2];
rz(0.71875188) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4913669) q[1];
sx q[1];
rz(-1.0904445) q[1];
sx q[1];
rz(-2.4841289) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80358867) q[3];
sx q[3];
rz(-0.39229624) q[3];
sx q[3];
rz(2.4367849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4031561) q[2];
sx q[2];
rz(-1.7313892) q[2];
sx q[2];
rz(2.9031244) q[2];
rz(-2.9078935) q[3];
sx q[3];
rz(-2.3668079) q[3];
sx q[3];
rz(0.15271798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23976633) q[0];
sx q[0];
rz(-0.16401839) q[0];
sx q[0];
rz(0.30583403) q[0];
rz(2.8726874) q[1];
sx q[1];
rz(-0.79335672) q[1];
sx q[1];
rz(-2.7893524) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3650318) q[0];
sx q[0];
rz(-1.7045341) q[0];
sx q[0];
rz(-0.040091776) q[0];
x q[1];
rz(-1.8369434) q[2];
sx q[2];
rz(-2.4983642) q[2];
sx q[2];
rz(-0.51189724) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8246987) q[1];
sx q[1];
rz(-2.3815063) q[1];
sx q[1];
rz(-2.555067) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78823678) q[3];
sx q[3];
rz(-1.394789) q[3];
sx q[3];
rz(0.22135425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29117808) q[2];
sx q[2];
rz(-1.520227) q[2];
sx q[2];
rz(-2.9746383) q[2];
rz(-0.11635612) q[3];
sx q[3];
rz(-0.51893026) q[3];
sx q[3];
rz(1.2701344) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51036924) q[0];
sx q[0];
rz(-0.33741697) q[0];
sx q[0];
rz(-1.1153197) q[0];
rz(-1.2315617) q[1];
sx q[1];
rz(-1.4304588) q[1];
sx q[1];
rz(3.0273052) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63681215) q[0];
sx q[0];
rz(-2.0858602) q[0];
sx q[0];
rz(-0.60487855) q[0];
x q[1];
rz(-0.27844001) q[2];
sx q[2];
rz(-1.8634708) q[2];
sx q[2];
rz(0.038624374) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7440577) q[1];
sx q[1];
rz(-1.1611132) q[1];
sx q[1];
rz(-0.68222218) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5763567) q[3];
sx q[3];
rz(-2.2768124) q[3];
sx q[3];
rz(1.3061641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.085122434) q[2];
sx q[2];
rz(-1.2135999) q[2];
sx q[2];
rz(-1.6950133) q[2];
rz(0.96212402) q[3];
sx q[3];
rz(-2.646793) q[3];
sx q[3];
rz(1.8112633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98244786) q[0];
sx q[0];
rz(-1.5626937) q[0];
sx q[0];
rz(2.6084117) q[0];
rz(1.5348596) q[1];
sx q[1];
rz(-2.3695562) q[1];
sx q[1];
rz(2.3981222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3822298) q[0];
sx q[0];
rz(-2.8385525) q[0];
sx q[0];
rz(2.8735301) q[0];
rz(-pi) q[1];
rz(2.1155223) q[2];
sx q[2];
rz(-1.3711978) q[2];
sx q[2];
rz(-0.064909086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.20821321) q[1];
sx q[1];
rz(-2.4543896) q[1];
sx q[1];
rz(2.739868) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1164769) q[3];
sx q[3];
rz(-2.4719596) q[3];
sx q[3];
rz(0.39519924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47034904) q[2];
sx q[2];
rz(-2.2942784) q[2];
sx q[2];
rz(-2.642855) q[2];
rz(-2.4399452) q[3];
sx q[3];
rz(-0.68877733) q[3];
sx q[3];
rz(-2.4832895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7094803) q[0];
sx q[0];
rz(-1.9146336) q[0];
sx q[0];
rz(0.099932583) q[0];
rz(-1.2808895) q[1];
sx q[1];
rz(-0.51191267) q[1];
sx q[1];
rz(2.4023712) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5422338) q[0];
sx q[0];
rz(-1.9606435) q[0];
sx q[0];
rz(-0.43655661) q[0];
rz(-pi) q[1];
rz(1.6373424) q[2];
sx q[2];
rz(-0.4302667) q[2];
sx q[2];
rz(-2.2815548) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91659613) q[1];
sx q[1];
rz(-0.9889305) q[1];
sx q[1];
rz(0.47486979) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5815063) q[3];
sx q[3];
rz(-2.3751343) q[3];
sx q[3];
rz(0.86957726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40333906) q[2];
sx q[2];
rz(-2.5103939) q[2];
sx q[2];
rz(1.6048019) q[2];
rz(1.4295476) q[3];
sx q[3];
rz(-1.6481383) q[3];
sx q[3];
rz(0.79609377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0159863) q[0];
sx q[0];
rz(-2.7637389) q[0];
sx q[0];
rz(2.0269537) q[0];
rz(0.58427018) q[1];
sx q[1];
rz(-2.22157) q[1];
sx q[1];
rz(-0.11411962) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7540365) q[0];
sx q[0];
rz(-2.4582131) q[0];
sx q[0];
rz(-0.46221931) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.592351) q[2];
sx q[2];
rz(-2.261544) q[2];
sx q[2];
rz(-0.88976414) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0551968) q[1];
sx q[1];
rz(-2.317531) q[1];
sx q[1];
rz(-2.1419129) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9769642) q[3];
sx q[3];
rz(-0.28985786) q[3];
sx q[3];
rz(1.690762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4383661) q[2];
sx q[2];
rz(-2.4195636) q[2];
sx q[2];
rz(-0.88877338) q[2];
rz(1.5416386) q[3];
sx q[3];
rz(-1.9346574) q[3];
sx q[3];
rz(-2.3204939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0515161) q[0];
sx q[0];
rz(-2.8335644) q[0];
sx q[0];
rz(2.8544881) q[0];
rz(1.2039394) q[1];
sx q[1];
rz(-2.2724889) q[1];
sx q[1];
rz(-1.6285816) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0368691) q[0];
sx q[0];
rz(-1.6469886) q[0];
sx q[0];
rz(0.33427396) q[0];
rz(-pi) q[1];
rz(2.2193935) q[2];
sx q[2];
rz(-1.2051799) q[2];
sx q[2];
rz(-0.43482414) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7375499) q[1];
sx q[1];
rz(-2.3047631) q[1];
sx q[1];
rz(1.5108193) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8230667) q[3];
sx q[3];
rz(-3.1384094) q[3];
sx q[3];
rz(2.0306272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5055351) q[2];
sx q[2];
rz(-1.9745741) q[2];
sx q[2];
rz(2.8324845) q[2];
rz(-1.2067893) q[3];
sx q[3];
rz(-0.38201067) q[3];
sx q[3];
rz(-0.25930723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9647656) q[0];
sx q[0];
rz(-0.72787705) q[0];
sx q[0];
rz(0.65560174) q[0];
rz(-0.62878311) q[1];
sx q[1];
rz(-1.7318334) q[1];
sx q[1];
rz(-1.383673) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7611279) q[0];
sx q[0];
rz(-1.5232067) q[0];
sx q[0];
rz(1.6143198) q[0];
rz(-pi) q[1];
rz(-1.0760154) q[2];
sx q[2];
rz(-2.3411334) q[2];
sx q[2];
rz(-0.59409522) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1392581) q[1];
sx q[1];
rz(-2.9911925) q[1];
sx q[1];
rz(-3.0207182) q[1];
rz(-pi) q[2];
rz(0.25383935) q[3];
sx q[3];
rz(-0.51341265) q[3];
sx q[3];
rz(-1.9771345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3722966) q[2];
sx q[2];
rz(-1.5089704) q[2];
sx q[2];
rz(-0.24202913) q[2];
rz(0.58487839) q[3];
sx q[3];
rz(-2.4681028) q[3];
sx q[3];
rz(2.6175595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4181716) q[0];
sx q[0];
rz(-2.5192498) q[0];
sx q[0];
rz(2.584516) q[0];
rz(2.2611025) q[1];
sx q[1];
rz(-1.4028032) q[1];
sx q[1];
rz(1.0407851) q[1];
rz(-1.7769832) q[2];
sx q[2];
rz(-1.5645116) q[2];
sx q[2];
rz(1.9581563) q[2];
rz(-2.1370173) q[3];
sx q[3];
rz(-1.6508816) q[3];
sx q[3];
rz(-2.3561867) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

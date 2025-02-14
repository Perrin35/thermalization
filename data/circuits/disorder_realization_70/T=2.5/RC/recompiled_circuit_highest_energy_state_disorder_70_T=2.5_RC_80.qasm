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
rz(-0.47082585) q[0];
sx q[0];
rz(-2.6915221) q[0];
sx q[0];
rz(0.39029628) q[0];
rz(0.14416873) q[1];
sx q[1];
rz(-1.498797) q[1];
sx q[1];
rz(2.1013451) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1821572) q[0];
sx q[0];
rz(-2.6819026) q[0];
sx q[0];
rz(0.88031405) q[0];
x q[1];
rz(-1.7261271) q[2];
sx q[2];
rz(-1.2110146) q[2];
sx q[2];
rz(1.8632165) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3398709) q[1];
sx q[1];
rz(-0.90870181) q[1];
sx q[1];
rz(-0.26317308) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.066180996) q[3];
sx q[3];
rz(-2.3621763) q[3];
sx q[3];
rz(-0.73252892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1990004) q[2];
sx q[2];
rz(-2.5358443) q[2];
sx q[2];
rz(-1.5462297) q[2];
rz(2.6664901) q[3];
sx q[3];
rz(-2.4964156) q[3];
sx q[3];
rz(-1.9861541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5652745) q[0];
sx q[0];
rz(-2.1529038) q[0];
sx q[0];
rz(-2.7069672) q[0];
rz(-2.077153) q[1];
sx q[1];
rz(-2.7480405) q[1];
sx q[1];
rz(-0.89148608) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9410649) q[0];
sx q[0];
rz(-2.5915798) q[0];
sx q[0];
rz(0.7329697) q[0];
x q[1];
rz(2.4154046) q[2];
sx q[2];
rz(-1.1788158) q[2];
sx q[2];
rz(1.7418944) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6124208) q[1];
sx q[1];
rz(-1.1193706) q[1];
sx q[1];
rz(2.5542817) q[1];
x q[2];
rz(1.6092898) q[3];
sx q[3];
rz(-0.6596237) q[3];
sx q[3];
rz(-1.6651814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6330304) q[2];
sx q[2];
rz(-0.56266251) q[2];
sx q[2];
rz(-2.0821345) q[2];
rz(1.2262454) q[3];
sx q[3];
rz(-1.5303333) q[3];
sx q[3];
rz(3.0219769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0480963) q[0];
sx q[0];
rz(-2.7432848) q[0];
sx q[0];
rz(0.68921971) q[0];
rz(-0.41220176) q[1];
sx q[1];
rz(-1.1176132) q[1];
sx q[1];
rz(1.162792) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0242694) q[0];
sx q[0];
rz(-1.7235959) q[0];
sx q[0];
rz(1.8788179) q[0];
x q[1];
rz(2.5291205) q[2];
sx q[2];
rz(-0.87650132) q[2];
sx q[2];
rz(1.5907254) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.38123576) q[1];
sx q[1];
rz(-2.3489526) q[1];
sx q[1];
rz(2.4355678) q[1];
rz(-pi) q[2];
x q[2];
rz(2.338004) q[3];
sx q[3];
rz(-2.7492964) q[3];
sx q[3];
rz(0.70480777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.73843655) q[2];
sx q[2];
rz(-1.4102035) q[2];
sx q[2];
rz(0.23846826) q[2];
rz(-0.23369914) q[3];
sx q[3];
rz(-0.77478474) q[3];
sx q[3];
rz(-2.9888747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9018263) q[0];
sx q[0];
rz(-0.16401839) q[0];
sx q[0];
rz(2.8357586) q[0];
rz(0.26890525) q[1];
sx q[1];
rz(-0.79335672) q[1];
sx q[1];
rz(2.7893524) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7765609) q[0];
sx q[0];
rz(-1.4370586) q[0];
sx q[0];
rz(0.040091776) q[0];
x q[1];
rz(2.1969123) q[2];
sx q[2];
rz(-1.4123823) q[2];
sx q[2];
rz(0.84412727) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.316894) q[1];
sx q[1];
rz(-0.76008633) q[1];
sx q[1];
rz(-2.555067) q[1];
rz(-0.24573487) q[3];
sx q[3];
rz(-2.338109) q[3];
sx q[3];
rz(1.5218228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8504146) q[2];
sx q[2];
rz(-1.6213657) q[2];
sx q[2];
rz(2.9746383) q[2];
rz(0.11635612) q[3];
sx q[3];
rz(-2.6226624) q[3];
sx q[3];
rz(-1.8714582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6312234) q[0];
sx q[0];
rz(-0.33741697) q[0];
sx q[0];
rz(2.0262729) q[0];
rz(-1.2315617) q[1];
sx q[1];
rz(-1.4304588) q[1];
sx q[1];
rz(-0.11428741) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5530722) q[0];
sx q[0];
rz(-0.77295303) q[0];
sx q[0];
rz(-2.3585178) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2670934) q[2];
sx q[2];
rz(-1.8371008) q[2];
sx q[2];
rz(-1.5271306) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3975349) q[1];
sx q[1];
rz(-1.9804795) q[1];
sx q[1];
rz(-0.68222218) q[1];
rz(-pi) q[2];
rz(-0.0065213477) q[3];
sx q[3];
rz(-2.4355585) q[3];
sx q[3];
rz(-1.8268585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.085122434) q[2];
sx q[2];
rz(-1.9279927) q[2];
sx q[2];
rz(1.6950133) q[2];
rz(0.96212402) q[3];
sx q[3];
rz(-0.4947997) q[3];
sx q[3];
rz(-1.8112633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98244786) q[0];
sx q[0];
rz(-1.5788989) q[0];
sx q[0];
rz(0.53318095) q[0];
rz(1.5348596) q[1];
sx q[1];
rz(-0.77203647) q[1];
sx q[1];
rz(0.74347043) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0395775) q[0];
sx q[0];
rz(-1.2789037) q[0];
sx q[0];
rz(1.4881698) q[0];
rz(-pi) q[1];
rz(-1.1986046) q[2];
sx q[2];
rz(-0.57665885) q[2];
sx q[2];
rz(1.8221591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.20821321) q[1];
sx q[1];
rz(-0.68720308) q[1];
sx q[1];
rz(-0.40172462) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1164769) q[3];
sx q[3];
rz(-0.66963306) q[3];
sx q[3];
rz(-0.39519924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6712436) q[2];
sx q[2];
rz(-2.2942784) q[2];
sx q[2];
rz(-0.49873763) q[2];
rz(2.4399452) q[3];
sx q[3];
rz(-2.4528153) q[3];
sx q[3];
rz(-2.4832895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43211234) q[0];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4866132) q[0];
sx q[0];
rz(-2.5648022) q[0];
sx q[0];
rz(-0.77107112) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1110837) q[2];
sx q[2];
rz(-1.1415452) q[2];
sx q[2];
rz(0.78684083) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2249965) q[1];
sx q[1];
rz(-0.9889305) q[1];
sx q[1];
rz(-0.47486979) q[1];
rz(-pi) q[2];
x q[2];
rz(0.010311237) q[3];
sx q[3];
rz(-2.3371994) q[3];
sx q[3];
rz(-0.88444405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40333906) q[2];
sx q[2];
rz(-0.6311987) q[2];
sx q[2];
rz(-1.5367907) q[2];
rz(1.712045) q[3];
sx q[3];
rz(-1.6481383) q[3];
sx q[3];
rz(-0.79609377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1256063) q[0];
sx q[0];
rz(-2.7637389) q[0];
sx q[0];
rz(-1.1146389) q[0];
rz(-0.58427018) q[1];
sx q[1];
rz(-0.92002267) q[1];
sx q[1];
rz(-0.11411962) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9561591) q[0];
sx q[0];
rz(-1.8562278) q[0];
sx q[0];
rz(2.5117842) q[0];
rz(0.54924168) q[2];
sx q[2];
rz(-0.88004862) q[2];
sx q[2];
rz(0.88976414) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.32875672) q[1];
sx q[1];
rz(-0.90531534) q[1];
sx q[1];
rz(-0.52862851) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8382254) q[3];
sx q[3];
rz(-1.6839613) q[3];
sx q[3];
rz(-0.51089493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4383661) q[2];
sx q[2];
rz(-0.72202903) q[2];
sx q[2];
rz(0.88877338) q[2];
rz(-1.5416386) q[3];
sx q[3];
rz(-1.2069353) q[3];
sx q[3];
rz(0.82109872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0515161) q[0];
sx q[0];
rz(-0.30802825) q[0];
sx q[0];
rz(2.8544881) q[0];
rz(1.9376532) q[1];
sx q[1];
rz(-0.86910373) q[1];
sx q[1];
rz(1.513011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2503165) q[0];
sx q[0];
rz(-0.34252942) q[0];
sx q[0];
rz(2.9129759) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44782717) q[2];
sx q[2];
rz(-0.97140233) q[2];
sx q[2];
rz(-2.2702655) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2681458) q[1];
sx q[1];
rz(-1.526274) q[1];
sx q[1];
rz(-2.4067307) q[1];
rz(1.5677139) q[3];
sx q[3];
rz(-1.5700018) q[3];
sx q[3];
rz(-2.4294927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.63605753) q[2];
sx q[2];
rz(-1.9745741) q[2];
sx q[2];
rz(-2.8324845) q[2];
rz(1.9348034) q[3];
sx q[3];
rz(-0.38201067) q[3];
sx q[3];
rz(2.8822854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1768271) q[0];
sx q[0];
rz(-2.4137156) q[0];
sx q[0];
rz(2.4859909) q[0];
rz(-0.62878311) q[1];
sx q[1];
rz(-1.4097593) q[1];
sx q[1];
rz(-1.7579196) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7611279) q[0];
sx q[0];
rz(-1.6183859) q[0];
sx q[0];
rz(1.5272729) q[0];
x q[1];
rz(-0.45510095) q[2];
sx q[2];
rz(-2.2544207) q[2];
sx q[2];
rz(0.065082642) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88009363) q[1];
sx q[1];
rz(-1.4215018) q[1];
sx q[1];
rz(-1.5525251) q[1];
rz(-pi) q[2];
rz(2.6419956) q[3];
sx q[3];
rz(-1.6944505) q[3];
sx q[3];
rz(2.5130003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76929602) q[2];
sx q[2];
rz(-1.6326222) q[2];
sx q[2];
rz(-2.8995635) q[2];
rz(-2.5567143) q[3];
sx q[3];
rz(-2.4681028) q[3];
sx q[3];
rz(2.6175595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.723421) q[0];
sx q[0];
rz(-0.62234288) q[0];
sx q[0];
rz(-0.55707669) q[0];
rz(2.2611025) q[1];
sx q[1];
rz(-1.4028032) q[1];
sx q[1];
rz(1.0407851) q[1];
rz(1.3646094) q[2];
sx q[2];
rz(-1.5645116) q[2];
sx q[2];
rz(1.9581563) q[2];
rz(-1.7193033) q[3];
sx q[3];
rz(-2.5703493) q[3];
sx q[3];
rz(-0.66019365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

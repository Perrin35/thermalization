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
rz(2.7512964) q[0];
rz(0.14416873) q[1];
sx q[1];
rz(4.7843883) q[1];
sx q[1];
rz(11.526123) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24867985) q[0];
sx q[0];
rz(-1.2843161) q[0];
sx q[0];
rz(1.2061968) q[0];
x q[1];
rz(-0.36379142) q[2];
sx q[2];
rz(-1.4254839) q[2];
sx q[2];
rz(2.7940968) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3398709) q[1];
sx q[1];
rz(-2.2328908) q[1];
sx q[1];
rz(2.8784196) q[1];
rz(-pi) q[2];
rz(-1.6360498) q[3];
sx q[3];
rz(-2.3480519) q[3];
sx q[3];
rz(-2.3161567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1990004) q[2];
sx q[2];
rz(-2.5358443) q[2];
sx q[2];
rz(1.5462297) q[2];
rz(0.47510251) q[3];
sx q[3];
rz(-2.4964156) q[3];
sx q[3];
rz(-1.1554385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.5763181) q[0];
sx q[0];
rz(-2.1529038) q[0];
sx q[0];
rz(-0.43462547) q[0];
rz(1.0644396) q[1];
sx q[1];
rz(-2.7480405) q[1];
sx q[1];
rz(-0.89148608) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9410649) q[0];
sx q[0];
rz(-0.55001282) q[0];
sx q[0];
rz(0.7329697) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0758249) q[2];
sx q[2];
rz(-2.2314851) q[2];
sx q[2];
rz(2.6434174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5198377) q[1];
sx q[1];
rz(-0.72417604) q[1];
sx q[1];
rz(-2.4228079) q[1];
rz(-1.5323029) q[3];
sx q[3];
rz(-2.481969) q[3];
sx q[3];
rz(1.6651814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6330304) q[2];
sx q[2];
rz(-0.56266251) q[2];
sx q[2];
rz(-1.0594581) q[2];
rz(-1.2262454) q[3];
sx q[3];
rz(-1.5303333) q[3];
sx q[3];
rz(0.1196158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0934963) q[0];
sx q[0];
rz(-0.39830783) q[0];
sx q[0];
rz(2.4523729) q[0];
rz(0.41220176) q[1];
sx q[1];
rz(-2.0239794) q[1];
sx q[1];
rz(1.162792) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6434518) q[0];
sx q[0];
rz(-1.8751133) q[0];
sx q[0];
rz(-2.9813719) q[0];
rz(-0.61247214) q[2];
sx q[2];
rz(-0.87650132) q[2];
sx q[2];
rz(-1.5508673) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4913669) q[1];
sx q[1];
rz(-1.0904445) q[1];
sx q[1];
rz(-0.65746376) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80358867) q[3];
sx q[3];
rz(-2.7492964) q[3];
sx q[3];
rz(2.4367849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.73843655) q[2];
sx q[2];
rz(-1.4102035) q[2];
sx q[2];
rz(2.9031244) q[2];
rz(0.23369914) q[3];
sx q[3];
rz(-0.77478474) q[3];
sx q[3];
rz(-0.15271798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9018263) q[0];
sx q[0];
rz(-2.9775743) q[0];
sx q[0];
rz(0.30583403) q[0];
rz(-2.8726874) q[1];
sx q[1];
rz(-0.79335672) q[1];
sx q[1];
rz(-0.3522402) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9304795) q[0];
sx q[0];
rz(-1.5310627) q[0];
sx q[0];
rz(1.7046403) q[0];
x q[1];
rz(-2.9469389) q[2];
sx q[2];
rz(-2.1878864) q[2];
sx q[2];
rz(2.301331) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4388386) q[1];
sx q[1];
rz(-1.9620336) q[1];
sx q[1];
rz(-0.66968285) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8504146) q[2];
sx q[2];
rz(-1.520227) q[2];
sx q[2];
rz(-2.9746383) q[2];
rz(3.0252365) q[3];
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
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6312234) q[0];
sx q[0];
rz(-2.8041757) q[0];
sx q[0];
rz(1.1153197) q[0];
rz(1.910031) q[1];
sx q[1];
rz(-1.4304588) q[1];
sx q[1];
rz(3.0273052) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60576263) q[0];
sx q[0];
rz(-2.0884984) q[0];
sx q[0];
rz(-2.1735031) q[0];
x q[1];
rz(1.2670934) q[2];
sx q[2];
rz(-1.8371008) q[2];
sx q[2];
rz(1.6144621) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3975349) q[1];
sx q[1];
rz(-1.1611132) q[1];
sx q[1];
rz(0.68222218) q[1];
x q[2];
rz(-0.0065213477) q[3];
sx q[3];
rz(-2.4355585) q[3];
sx q[3];
rz(1.3147342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.085122434) q[2];
sx q[2];
rz(-1.9279927) q[2];
sx q[2];
rz(-1.6950133) q[2];
rz(-2.1794686) q[3];
sx q[3];
rz(-2.646793) q[3];
sx q[3];
rz(-1.3303293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1591448) q[0];
sx q[0];
rz(-1.5788989) q[0];
sx q[0];
rz(0.53318095) q[0];
rz(1.606733) q[1];
sx q[1];
rz(-2.3695562) q[1];
sx q[1];
rz(0.74347043) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7593628) q[0];
sx q[0];
rz(-2.8385525) q[0];
sx q[0];
rz(2.8735301) q[0];
rz(-2.9093365) q[2];
sx q[2];
rz(-2.1035367) q[2];
sx q[2];
rz(-1.7552623) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.71072342) q[1];
sx q[1];
rz(-2.1942881) q[1];
sx q[1];
rz(1.2602978) q[1];
rz(-pi) q[2];
rz(-2.1891001) q[3];
sx q[3];
rz(-1.8466766) q[3];
sx q[3];
rz(-2.3316958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.47034904) q[2];
sx q[2];
rz(-0.84731421) q[2];
sx q[2];
rz(2.642855) q[2];
rz(2.4399452) q[3];
sx q[3];
rz(-2.4528153) q[3];
sx q[3];
rz(-2.4832895) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7094803) q[0];
sx q[0];
rz(-1.9146336) q[0];
sx q[0];
rz(3.0416601) q[0];
rz(1.8607032) q[1];
sx q[1];
rz(-0.51191267) q[1];
sx q[1];
rz(2.4023712) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6549795) q[0];
sx q[0];
rz(-2.5648022) q[0];
sx q[0];
rz(-2.3705215) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0002236) q[2];
sx q[2];
rz(-1.5985367) q[2];
sx q[2];
rz(-0.77125473) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9770551) q[1];
sx q[1];
rz(-2.4083369) q[1];
sx q[1];
rz(0.96340839) q[1];
rz(-2.337226) q[3];
sx q[3];
rz(-1.5633681) q[3];
sx q[3];
rz(2.4480889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7382536) q[2];
sx q[2];
rz(-0.6311987) q[2];
sx q[2];
rz(-1.5367907) q[2];
rz(1.4295476) q[3];
sx q[3];
rz(-1.6481383) q[3];
sx q[3];
rz(-2.3454989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0159863) q[0];
sx q[0];
rz(-2.7637389) q[0];
sx q[0];
rz(-1.1146389) q[0];
rz(-2.5573225) q[1];
sx q[1];
rz(-2.22157) q[1];
sx q[1];
rz(3.027473) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9586351) q[0];
sx q[0];
rz(-0.97016969) q[0];
sx q[0];
rz(-1.9191027) q[0];
rz(-pi) q[1];
rz(2.3405206) q[2];
sx q[2];
rz(-1.1566887) q[2];
sx q[2];
rz(1.0528477) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0863959) q[1];
sx q[1];
rz(-2.317531) q[1];
sx q[1];
rz(2.1419129) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11729838) q[3];
sx q[3];
rz(-1.8364732) q[3];
sx q[3];
rz(2.1126216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4383661) q[2];
sx q[2];
rz(-0.72202903) q[2];
sx q[2];
rz(0.88877338) q[2];
rz(-1.599954) q[3];
sx q[3];
rz(-1.9346574) q[3];
sx q[3];
rz(-2.3204939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0515161) q[0];
sx q[0];
rz(-0.30802825) q[0];
sx q[0];
rz(-2.8544881) q[0];
rz(-1.2039394) q[1];
sx q[1];
rz(-0.86910373) q[1];
sx q[1];
rz(1.513011) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10472351) q[0];
sx q[0];
rz(-1.6469886) q[0];
sx q[0];
rz(2.8073187) q[0];
x q[1];
rz(-2.1356629) q[2];
sx q[2];
rz(-0.73137368) q[2];
sx q[2];
rz(-1.5767026) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4934532) q[1];
sx q[1];
rz(-0.73595787) q[1];
sx q[1];
rz(3.0752431) q[1];
rz(-pi) q[2];
rz(1.8230667) q[3];
sx q[3];
rz(-3.1384094) q[3];
sx q[3];
rz(2.0306272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5055351) q[2];
sx q[2];
rz(-1.9745741) q[2];
sx q[2];
rz(-0.30910811) q[2];
rz(1.9348034) q[3];
sx q[3];
rz(-2.759582) q[3];
sx q[3];
rz(0.25930723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.9647656) q[0];
sx q[0];
rz(-2.4137156) q[0];
sx q[0];
rz(-0.65560174) q[0];
rz(0.62878311) q[1];
sx q[1];
rz(-1.7318334) q[1];
sx q[1];
rz(1.383673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38046471) q[0];
sx q[0];
rz(-1.6183859) q[0];
sx q[0];
rz(-1.6143198) q[0];
rz(0.83413275) q[2];
sx q[2];
rz(-1.9185432) q[2];
sx q[2];
rz(1.8054838) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0023345) q[1];
sx q[1];
rz(-0.15040018) q[1];
sx q[1];
rz(-0.12087442) q[1];
rz(-pi) q[2];
rz(-2.8877533) q[3];
sx q[3];
rz(-0.51341265) q[3];
sx q[3];
rz(-1.9771345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3722966) q[2];
sx q[2];
rz(-1.6326222) q[2];
sx q[2];
rz(2.8995635) q[2];
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
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.723421) q[0];
sx q[0];
rz(-2.5192498) q[0];
sx q[0];
rz(2.584516) q[0];
rz(0.88049018) q[1];
sx q[1];
rz(-1.7387895) q[1];
sx q[1];
rz(-2.1008076) q[1];
rz(3.1351719) q[2];
sx q[2];
rz(-1.7769791) q[2];
sx q[2];
rz(0.38867452) q[2];
rz(3.0467792) q[3];
sx q[3];
rz(-1.0066114) q[3];
sx q[3];
rz(-0.8361984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

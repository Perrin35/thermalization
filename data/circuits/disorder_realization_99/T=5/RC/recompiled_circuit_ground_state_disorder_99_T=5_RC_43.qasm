OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91166624) q[0];
sx q[0];
rz(-0.32719964) q[0];
sx q[0];
rz(1.3258452) q[0];
rz(-1.7757379) q[1];
sx q[1];
rz(-2.7471625) q[1];
sx q[1];
rz(2.3529513) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45730725) q[0];
sx q[0];
rz(-2.4198902) q[0];
sx q[0];
rz(-1.0941605) q[0];
rz(-2.3173877) q[2];
sx q[2];
rz(-1.0366777) q[2];
sx q[2];
rz(0.20301486) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44581283) q[1];
sx q[1];
rz(-1.7632293) q[1];
sx q[1];
rz(0.28833994) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12874916) q[3];
sx q[3];
rz(-0.94673079) q[3];
sx q[3];
rz(2.9381813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6476562) q[2];
sx q[2];
rz(-1.6948573) q[2];
sx q[2];
rz(-0.20230618) q[2];
rz(0.10937396) q[3];
sx q[3];
rz(-0.60905639) q[3];
sx q[3];
rz(-0.085453184) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0789455) q[0];
sx q[0];
rz(-2.2181692) q[0];
sx q[0];
rz(-1.9664636) q[0];
rz(1.6773978) q[1];
sx q[1];
rz(-1.0025832) q[1];
sx q[1];
rz(0.35951231) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8749411) q[0];
sx q[0];
rz(-1.4066937) q[0];
sx q[0];
rz(-3.0709188) q[0];
x q[1];
rz(-2.5746042) q[2];
sx q[2];
rz(-2.5051281) q[2];
sx q[2];
rz(-1.8228965) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0766616) q[1];
sx q[1];
rz(-2.1691469) q[1];
sx q[1];
rz(1.9877355) q[1];
x q[2];
rz(0.80163281) q[3];
sx q[3];
rz(-2.3838861) q[3];
sx q[3];
rz(-3.0996292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.29490617) q[2];
sx q[2];
rz(-1.0470231) q[2];
sx q[2];
rz(0.27493757) q[2];
rz(-1.8168195) q[3];
sx q[3];
rz(-1.6144269) q[3];
sx q[3];
rz(1.2980609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0777968) q[0];
sx q[0];
rz(-2.9834788) q[0];
sx q[0];
rz(-2.7445444) q[0];
rz(-0.60931698) q[1];
sx q[1];
rz(-1.3343697) q[1];
sx q[1];
rz(1.5999925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9843053) q[0];
sx q[0];
rz(-1.2799036) q[0];
sx q[0];
rz(0.10910927) q[0];
rz(1.3808525) q[2];
sx q[2];
rz(-3.0751588) q[2];
sx q[2];
rz(-2.9125467) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.065464547) q[1];
sx q[1];
rz(-0.72678002) q[1];
sx q[1];
rz(-1.5453592) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81542374) q[3];
sx q[3];
rz(-2.2924726) q[3];
sx q[3];
rz(-1.7343002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2043173) q[2];
sx q[2];
rz(-1.5949564) q[2];
sx q[2];
rz(-2.9079962) q[2];
rz(1.8448081) q[3];
sx q[3];
rz(-1.9498884) q[3];
sx q[3];
rz(0.72062033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5066756) q[0];
sx q[0];
rz(-1.1577865) q[0];
sx q[0];
rz(-1.7764212) q[0];
rz(-1.8695976) q[1];
sx q[1];
rz(-1.9422266) q[1];
sx q[1];
rz(1.8796399) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.711447) q[0];
sx q[0];
rz(-1.0017348) q[0];
sx q[0];
rz(1.7993159) q[0];
rz(-pi) q[1];
rz(1.2620401) q[2];
sx q[2];
rz(-1.595904) q[2];
sx q[2];
rz(1.0756191) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.024725951) q[1];
sx q[1];
rz(-0.19389175) q[1];
sx q[1];
rz(-2.7996382) q[1];
rz(-0.19660321) q[3];
sx q[3];
rz(-0.8340618) q[3];
sx q[3];
rz(-1.3722668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9103553) q[2];
sx q[2];
rz(-2.1114025) q[2];
sx q[2];
rz(-2.3528698) q[2];
rz(-1.8606868) q[3];
sx q[3];
rz(-1.7773881) q[3];
sx q[3];
rz(-1.0219215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7541517) q[0];
sx q[0];
rz(-1.0197637) q[0];
sx q[0];
rz(0.66106853) q[0];
rz(-2.0263653) q[1];
sx q[1];
rz(-1.8293569) q[1];
sx q[1];
rz(-2.1818395) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57908995) q[0];
sx q[0];
rz(-1.3013562) q[0];
sx q[0];
rz(2.605382) q[0];
x q[1];
rz(-2.9767545) q[2];
sx q[2];
rz(-2.1123114) q[2];
sx q[2];
rz(2.7680226) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5603191) q[1];
sx q[1];
rz(-1.4380102) q[1];
sx q[1];
rz(-1.8627235) q[1];
rz(-0.92613756) q[3];
sx q[3];
rz(-1.411866) q[3];
sx q[3];
rz(1.8862806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58012086) q[2];
sx q[2];
rz(-0.46130195) q[2];
sx q[2];
rz(-1.0599773) q[2];
rz(-0.75016841) q[3];
sx q[3];
rz(-1.3786022) q[3];
sx q[3];
rz(1.2257956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2724514) q[0];
sx q[0];
rz(-1.5831818) q[0];
sx q[0];
rz(-2.9172752) q[0];
rz(-1.51651) q[1];
sx q[1];
rz(-0.61619157) q[1];
sx q[1];
rz(-0.075627653) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13091732) q[0];
sx q[0];
rz(-1.0300127) q[0];
sx q[0];
rz(-2.6877139) q[0];
rz(-pi) q[1];
rz(-0.096616726) q[2];
sx q[2];
rz(-1.268309) q[2];
sx q[2];
rz(0.58524489) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4757753) q[1];
sx q[1];
rz(-2.4815637) q[1];
sx q[1];
rz(-0.81190656) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0527961) q[3];
sx q[3];
rz(-0.78556934) q[3];
sx q[3];
rz(-1.1256684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6077967) q[2];
sx q[2];
rz(-2.3888612) q[2];
sx q[2];
rz(-0.047018615) q[2];
rz(-2.4646344) q[3];
sx q[3];
rz(-1.8432063) q[3];
sx q[3];
rz(-3.1162139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9852279) q[0];
sx q[0];
rz(-1.6522836) q[0];
sx q[0];
rz(1.697668) q[0];
rz(-2.2185183) q[1];
sx q[1];
rz(-2.1363027) q[1];
sx q[1];
rz(2.9583171) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8240487) q[0];
sx q[0];
rz(-2.5185761) q[0];
sx q[0];
rz(-2.7035294) q[0];
x q[1];
rz(1.1475816) q[2];
sx q[2];
rz(-1.1626889) q[2];
sx q[2];
rz(-2.0988879) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.35872719) q[1];
sx q[1];
rz(-2.065383) q[1];
sx q[1];
rz(-0.371277) q[1];
rz(-2.1499726) q[3];
sx q[3];
rz(-1.8133111) q[3];
sx q[3];
rz(-0.33580599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.034059374) q[2];
sx q[2];
rz(-1.4564161) q[2];
sx q[2];
rz(-3.122186) q[2];
rz(2.9803045) q[3];
sx q[3];
rz(-2.1262157) q[3];
sx q[3];
rz(-1.9136782) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1041383) q[0];
sx q[0];
rz(-0.4087953) q[0];
sx q[0];
rz(-0.092305146) q[0];
rz(-1.1760938) q[1];
sx q[1];
rz(-2.8832925) q[1];
sx q[1];
rz(1.4091122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69337849) q[0];
sx q[0];
rz(-1.4657768) q[0];
sx q[0];
rz(0.054570507) q[0];
x q[1];
rz(2.1966366) q[2];
sx q[2];
rz(-2.4279478) q[2];
sx q[2];
rz(-1.2882441) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17615023) q[1];
sx q[1];
rz(-1.9561844) q[1];
sx q[1];
rz(1.2809281) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73958379) q[3];
sx q[3];
rz(-0.7639262) q[3];
sx q[3];
rz(-0.52810625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.041542355) q[2];
sx q[2];
rz(-1.5644194) q[2];
sx q[2];
rz(1.9680295) q[2];
rz(0.38343492) q[3];
sx q[3];
rz(-1.8337367) q[3];
sx q[3];
rz(-2.0508544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0122871) q[0];
sx q[0];
rz(-1.6467935) q[0];
sx q[0];
rz(-0.20558414) q[0];
rz(-0.81941098) q[1];
sx q[1];
rz(-1.9267547) q[1];
sx q[1];
rz(2.6507071) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48614472) q[0];
sx q[0];
rz(-1.1974338) q[0];
sx q[0];
rz(-2.4400737) q[0];
x q[1];
rz(0.41342469) q[2];
sx q[2];
rz(-2.4707762) q[2];
sx q[2];
rz(2.6476423) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1236022) q[1];
sx q[1];
rz(-1.523087) q[1];
sx q[1];
rz(-1.0022267) q[1];
rz(0.098693178) q[3];
sx q[3];
rz(-2.0798363) q[3];
sx q[3];
rz(-2.4157897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4129591) q[2];
sx q[2];
rz(-2.1042991) q[2];
sx q[2];
rz(-1.884985) q[2];
rz(0.43577731) q[3];
sx q[3];
rz(-1.1083009) q[3];
sx q[3];
rz(-2.2573439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7213781) q[0];
sx q[0];
rz(-0.90737897) q[0];
sx q[0];
rz(-0.23289982) q[0];
rz(0.13889343) q[1];
sx q[1];
rz(-1.9878309) q[1];
sx q[1];
rz(2.6331666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.229436) q[0];
sx q[0];
rz(-0.60750738) q[0];
sx q[0];
rz(-1.1077393) q[0];
x q[1];
rz(-0.2744197) q[2];
sx q[2];
rz(-1.0792582) q[2];
sx q[2];
rz(-1.4948782) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.57181035) q[1];
sx q[1];
rz(-1.5787193) q[1];
sx q[1];
rz(-1.7047764) q[1];
rz(-pi) q[2];
rz(0.94975779) q[3];
sx q[3];
rz(-2.3425079) q[3];
sx q[3];
rz(-1.0423555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7306708) q[2];
sx q[2];
rz(-1.6795009) q[2];
sx q[2];
rz(0.81653583) q[2];
rz(-0.38980347) q[3];
sx q[3];
rz(-1.0023508) q[3];
sx q[3];
rz(-1.3780814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.447406) q[0];
sx q[0];
rz(-2.2185855) q[0];
sx q[0];
rz(2.9119281) q[0];
rz(1.6387088) q[1];
sx q[1];
rz(-1.3353744) q[1];
sx q[1];
rz(-2.2257805) q[1];
rz(-0.38328485) q[2];
sx q[2];
rz(-1.4105759) q[2];
sx q[2];
rz(-2.069266) q[2];
rz(-0.16666804) q[3];
sx q[3];
rz(-0.87571908) q[3];
sx q[3];
rz(-2.2805211) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

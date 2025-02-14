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
rz(-2.6925447) q[0];
sx q[0];
rz(-1.1392051) q[0];
sx q[0];
rz(2.6989302) q[0];
rz(0.17065419) q[1];
sx q[1];
rz(-2.3499188) q[1];
sx q[1];
rz(0.72905529) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1281888) q[0];
sx q[0];
rz(-2.3263651) q[0];
sx q[0];
rz(0.56437738) q[0];
rz(0.38972028) q[2];
sx q[2];
rz(-1.5699631) q[2];
sx q[2];
rz(-2.4454947) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2393136) q[1];
sx q[1];
rz(-0.89717509) q[1];
sx q[1];
rz(2.8301653) q[1];
rz(0.9425052) q[3];
sx q[3];
rz(-1.2344591) q[3];
sx q[3];
rz(2.9021341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2810104) q[2];
sx q[2];
rz(-3.0594825) q[2];
sx q[2];
rz(0.78544593) q[2];
rz(-2.1752518) q[3];
sx q[3];
rz(-1.0680501) q[3];
sx q[3];
rz(2.5016224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.53168374) q[0];
sx q[0];
rz(-2.5737679) q[0];
sx q[0];
rz(-2.5058643) q[0];
rz(0.6334148) q[1];
sx q[1];
rz(-0.76786357) q[1];
sx q[1];
rz(2.8107218) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71459717) q[0];
sx q[0];
rz(-1.2408371) q[0];
sx q[0];
rz(-0.64100463) q[0];
rz(-pi) q[1];
rz(1.517968) q[2];
sx q[2];
rz(-2.6739648) q[2];
sx q[2];
rz(0.29948452) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0111661) q[1];
sx q[1];
rz(-0.067421801) q[1];
sx q[1];
rz(-1.0481337) q[1];
rz(-pi) q[2];
rz(-0.14563609) q[3];
sx q[3];
rz(-3.0560827) q[3];
sx q[3];
rz(-0.95216361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.50515455) q[2];
sx q[2];
rz(-2.324489) q[2];
sx q[2];
rz(2.6804697) q[2];
rz(-0.72502208) q[3];
sx q[3];
rz(-0.45261639) q[3];
sx q[3];
rz(0.66864526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4918936) q[0];
sx q[0];
rz(-1.7576341) q[0];
sx q[0];
rz(2.6032676) q[0];
rz(0.0093983924) q[1];
sx q[1];
rz(-2.5098269) q[1];
sx q[1];
rz(1.1993154) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56012404) q[0];
sx q[0];
rz(-2.631979) q[0];
sx q[0];
rz(-1.6625924) q[0];
x q[1];
rz(-1.7117768) q[2];
sx q[2];
rz(-2.2491697) q[2];
sx q[2];
rz(1.3305261) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.44502434) q[1];
sx q[1];
rz(-2.4557132) q[1];
sx q[1];
rz(-2.1522151) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84405293) q[3];
sx q[3];
rz(-1.7622602) q[3];
sx q[3];
rz(1.307631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7678396) q[2];
sx q[2];
rz(-1.491863) q[2];
sx q[2];
rz(-2.4191432) q[2];
rz(-2.5433698) q[3];
sx q[3];
rz(-0.0349508) q[3];
sx q[3];
rz(1.9237062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5825321) q[0];
sx q[0];
rz(-1.8836972) q[0];
sx q[0];
rz(0.019056888) q[0];
rz(-1.4590774) q[1];
sx q[1];
rz(-2.9190813) q[1];
sx q[1];
rz(-0.51024514) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048240926) q[0];
sx q[0];
rz(-0.66931319) q[0];
sx q[0];
rz(-2.0055683) q[0];
rz(-pi) q[1];
rz(-1.8809659) q[2];
sx q[2];
rz(-0.87291832) q[2];
sx q[2];
rz(0.073402799) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7170808) q[1];
sx q[1];
rz(-1.1649017) q[1];
sx q[1];
rz(-0.67274977) q[1];
x q[2];
rz(2.8250349) q[3];
sx q[3];
rz(-2.7084841) q[3];
sx q[3];
rz(0.18726997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9973008) q[2];
sx q[2];
rz(-1.4520626) q[2];
sx q[2];
rz(-0.22225456) q[2];
rz(0.079515919) q[3];
sx q[3];
rz(-2.5955718) q[3];
sx q[3];
rz(0.63681805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52614373) q[0];
sx q[0];
rz(-2.0455102) q[0];
sx q[0];
rz(1.8002321) q[0];
rz(0.51097956) q[1];
sx q[1];
rz(-1.7781517) q[1];
sx q[1];
rz(-2.708639) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2596672) q[0];
sx q[0];
rz(-2.8779463) q[0];
sx q[0];
rz(-0.22048283) q[0];
x q[1];
rz(2.0413002) q[2];
sx q[2];
rz(-2.6950201) q[2];
sx q[2];
rz(0.64556015) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4652327) q[1];
sx q[1];
rz(-1.3238878) q[1];
sx q[1];
rz(1.2311155) q[1];
rz(-pi) q[2];
rz(1.8515685) q[3];
sx q[3];
rz(-2.3237356) q[3];
sx q[3];
rz(-0.75033891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0114835) q[2];
sx q[2];
rz(-0.92769647) q[2];
sx q[2];
rz(-0.9978655) q[2];
rz(-2.4326371) q[3];
sx q[3];
rz(-2.7537789) q[3];
sx q[3];
rz(-2.1628105) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82584941) q[0];
sx q[0];
rz(-0.5873) q[0];
sx q[0];
rz(-2.24776) q[0];
rz(-1.029344) q[1];
sx q[1];
rz(-0.7411595) q[1];
sx q[1];
rz(-3.0267874) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.478677) q[0];
sx q[0];
rz(-1.3548242) q[0];
sx q[0];
rz(1.8195282) q[0];
rz(-pi) q[1];
rz(-1.9283781) q[2];
sx q[2];
rz(-2.1809357) q[2];
sx q[2];
rz(-1.4550191) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0521212) q[1];
sx q[1];
rz(-1.2942181) q[1];
sx q[1];
rz(2.8514423) q[1];
rz(2.8607058) q[3];
sx q[3];
rz(-1.6320737) q[3];
sx q[3];
rz(0.83524246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43742314) q[2];
sx q[2];
rz(-0.88714522) q[2];
sx q[2];
rz(0.21370055) q[2];
rz(2.5050488) q[3];
sx q[3];
rz(-2.611219) q[3];
sx q[3];
rz(2.4056733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32695025) q[0];
sx q[0];
rz(-2.3982168) q[0];
sx q[0];
rz(-3.0315234) q[0];
rz(0.07668177) q[1];
sx q[1];
rz(-0.86014599) q[1];
sx q[1];
rz(-2.0032047) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6784043) q[0];
sx q[0];
rz(-0.41751719) q[0];
sx q[0];
rz(0.56708401) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3537972) q[2];
sx q[2];
rz(-1.5681364) q[2];
sx q[2];
rz(0.12622368) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0323022) q[1];
sx q[1];
rz(-1.3869973) q[1];
sx q[1];
rz(1.079783) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8951169) q[3];
sx q[3];
rz(-1.3630056) q[3];
sx q[3];
rz(-2.7026619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.25973213) q[2];
sx q[2];
rz(-1.0819819) q[2];
sx q[2];
rz(-2.760375) q[2];
rz(-3.0242053) q[3];
sx q[3];
rz(-2.6365247) q[3];
sx q[3];
rz(-0.086732619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36670244) q[0];
sx q[0];
rz(-2.369304) q[0];
sx q[0];
rz(-0.49852398) q[0];
rz(2.3122834) q[1];
sx q[1];
rz(-2.3259951) q[1];
sx q[1];
rz(-0.097749762) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5437336) q[0];
sx q[0];
rz(-1.2654788) q[0];
sx q[0];
rz(-1.1174563) q[0];
rz(2.7404064) q[2];
sx q[2];
rz(-2.7747535) q[2];
sx q[2];
rz(-0.40415472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6548582) q[1];
sx q[1];
rz(-2.2796101) q[1];
sx q[1];
rz(2.4177119) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5427715) q[3];
sx q[3];
rz(-2.8198979) q[3];
sx q[3];
rz(-3.0672586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8349614) q[2];
sx q[2];
rz(-1.1502879) q[2];
sx q[2];
rz(-2.5508896) q[2];
rz(3.0513884) q[3];
sx q[3];
rz(-0.43961757) q[3];
sx q[3];
rz(2.2265767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0331994) q[0];
sx q[0];
rz(-2.2985701) q[0];
sx q[0];
rz(-0.75476187) q[0];
rz(-2.8056878) q[1];
sx q[1];
rz(-2.5867808) q[1];
sx q[1];
rz(-2.1051443) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51506804) q[0];
sx q[0];
rz(-1.4534411) q[0];
sx q[0];
rz(0.68674725) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48123863) q[2];
sx q[2];
rz(-1.5526774) q[2];
sx q[2];
rz(-0.5258403) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9687313) q[1];
sx q[1];
rz(-1.5318499) q[1];
sx q[1];
rz(2.624818) q[1];
rz(1.3410946) q[3];
sx q[3];
rz(-1.5840845) q[3];
sx q[3];
rz(-0.54315996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.094940946) q[2];
sx q[2];
rz(-1.4163821) q[2];
sx q[2];
rz(-0.45647344) q[2];
rz(-0.60232317) q[3];
sx q[3];
rz(-2.9538302) q[3];
sx q[3];
rz(2.203234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25191054) q[0];
sx q[0];
rz(-1.2924117) q[0];
sx q[0];
rz(-0.46345261) q[0];
rz(-3.1262596) q[1];
sx q[1];
rz(-0.066611193) q[1];
sx q[1];
rz(-0.32036805) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.326691) q[0];
sx q[0];
rz(-1.4938271) q[0];
sx q[0];
rz(0.14369731) q[0];
rz(-2.8276075) q[2];
sx q[2];
rz(-1.3330402) q[2];
sx q[2];
rz(-1.2842922) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55904222) q[1];
sx q[1];
rz(-2.1782785) q[1];
sx q[1];
rz(-2.5431125) q[1];
x q[2];
rz(1.8611365) q[3];
sx q[3];
rz(-1.8986033) q[3];
sx q[3];
rz(-0.39738712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1128803) q[2];
sx q[2];
rz(-2.4148603) q[2];
sx q[2];
rz(1.0090562) q[2];
rz(-2.3446337) q[3];
sx q[3];
rz(-1.2538486) q[3];
sx q[3];
rz(2.8033281) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0491895) q[0];
sx q[0];
rz(-1.6272463) q[0];
sx q[0];
rz(1.8820681) q[0];
rz(0.084820329) q[1];
sx q[1];
rz(-1.6592574) q[1];
sx q[1];
rz(1.4316373) q[1];
rz(-0.063719393) q[2];
sx q[2];
rz(-1.812171) q[2];
sx q[2];
rz(-2.4952427) q[2];
rz(1.3884801) q[3];
sx q[3];
rz(-1.8686466) q[3];
sx q[3];
rz(2.2212088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

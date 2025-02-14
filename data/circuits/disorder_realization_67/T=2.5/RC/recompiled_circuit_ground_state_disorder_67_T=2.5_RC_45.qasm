OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.32970348) q[0];
sx q[0];
rz(-2.831037) q[0];
sx q[0];
rz(1.9429053) q[0];
rz(-2.0074453) q[1];
sx q[1];
rz(-0.79678798) q[1];
sx q[1];
rz(-2.8400583) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26623959) q[0];
sx q[0];
rz(-1.4090898) q[0];
sx q[0];
rz(-1.7577359) q[0];
rz(2.1211336) q[2];
sx q[2];
rz(-0.34780234) q[2];
sx q[2];
rz(-2.2344799) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.49820013) q[1];
sx q[1];
rz(-0.82005703) q[1];
sx q[1];
rz(-3.0489075) q[1];
rz(-1.1687247) q[3];
sx q[3];
rz(-1.3969867) q[3];
sx q[3];
rz(0.92648348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.6887168) q[2];
sx q[2];
rz(-2.0417002) q[2];
sx q[2];
rz(-1.1348881) q[2];
rz(0.12198837) q[3];
sx q[3];
rz(-2.205409) q[3];
sx q[3];
rz(-3.0854935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4695456) q[0];
sx q[0];
rz(-2.8944954) q[0];
sx q[0];
rz(3.1392414) q[0];
rz(2.7408842) q[1];
sx q[1];
rz(-1.7727163) q[1];
sx q[1];
rz(2.2854663) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8979748) q[0];
sx q[0];
rz(-1.8254733) q[0];
sx q[0];
rz(-2.2092878) q[0];
rz(-2.8555774) q[2];
sx q[2];
rz(-1.0723503) q[2];
sx q[2];
rz(-1.9749383) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0217738) q[1];
sx q[1];
rz(-2.0975862) q[1];
sx q[1];
rz(-3.0837545) q[1];
rz(0.54243954) q[3];
sx q[3];
rz(-0.63839165) q[3];
sx q[3];
rz(0.29104376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.037228435) q[2];
sx q[2];
rz(-1.9375485) q[2];
sx q[2];
rz(-0.47002235) q[2];
rz(-2.5445599) q[3];
sx q[3];
rz(-0.11490122) q[3];
sx q[3];
rz(1.443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3700579) q[0];
sx q[0];
rz(-0.4929339) q[0];
sx q[0];
rz(-0.22802995) q[0];
rz(-0.44949624) q[1];
sx q[1];
rz(-2.1411965) q[1];
sx q[1];
rz(-2.1953348) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7957009) q[0];
sx q[0];
rz(-1.3145224) q[0];
sx q[0];
rz(-1.9836224) q[0];
rz(-pi) q[1];
rz(-2.3870638) q[2];
sx q[2];
rz(-2.3268394) q[2];
sx q[2];
rz(-0.94902869) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6201278) q[1];
sx q[1];
rz(-1.2343654) q[1];
sx q[1];
rz(1.7640339) q[1];
x q[2];
rz(1.854293) q[3];
sx q[3];
rz(-1.4186267) q[3];
sx q[3];
rz(-2.7622472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9049282) q[2];
sx q[2];
rz(-1.6764574) q[2];
sx q[2];
rz(1.6620592) q[2];
rz(-2.9605401) q[3];
sx q[3];
rz(-0.79020399) q[3];
sx q[3];
rz(2.2553867) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36412305) q[0];
sx q[0];
rz(-1.5793707) q[0];
sx q[0];
rz(2.2430578) q[0];
rz(-2.6354375) q[1];
sx q[1];
rz(-1.2944784) q[1];
sx q[1];
rz(1.6815965) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0943992) q[0];
sx q[0];
rz(-1.6364105) q[0];
sx q[0];
rz(2.608247) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3982749) q[2];
sx q[2];
rz(-1.2365336) q[2];
sx q[2];
rz(1.0570132) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0467028) q[1];
sx q[1];
rz(-3.0071435) q[1];
sx q[1];
rz(1.0195169) q[1];
x q[2];
rz(-1.6494895) q[3];
sx q[3];
rz(-0.65641145) q[3];
sx q[3];
rz(2.479913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26607457) q[2];
sx q[2];
rz(-0.47553277) q[2];
sx q[2];
rz(1.6391222) q[2];
rz(1.255704) q[3];
sx q[3];
rz(-1.4016822) q[3];
sx q[3];
rz(0.046253117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2332377) q[0];
sx q[0];
rz(-2.4112356) q[0];
sx q[0];
rz(2.9610942) q[0];
rz(1.7620697) q[1];
sx q[1];
rz(-2.5738398) q[1];
sx q[1];
rz(-1.0951805) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22893208) q[0];
sx q[0];
rz(-1.7062418) q[0];
sx q[0];
rz(1.6333461) q[0];
x q[1];
rz(-3.0600431) q[2];
sx q[2];
rz(-2.3429541) q[2];
sx q[2];
rz(-2.7872686) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3928581) q[1];
sx q[1];
rz(-1.7104516) q[1];
sx q[1];
rz(0.77844324) q[1];
rz(-pi) q[2];
rz(3.0083382) q[3];
sx q[3];
rz(-2.2550003) q[3];
sx q[3];
rz(-1.2017991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77202648) q[2];
sx q[2];
rz(-0.84731805) q[2];
sx q[2];
rz(-1.0687211) q[2];
rz(-2.7145743) q[3];
sx q[3];
rz(-0.87978274) q[3];
sx q[3];
rz(-1.8900185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0291979) q[0];
sx q[0];
rz(-2.2258832) q[0];
sx q[0];
rz(-2.8503964) q[0];
rz(1.683782) q[1];
sx q[1];
rz(-2.1144512) q[1];
sx q[1];
rz(-2.2529032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5989953) q[0];
sx q[0];
rz(-2.3610736) q[0];
sx q[0];
rz(0.88135834) q[0];
x q[1];
rz(2.3843147) q[2];
sx q[2];
rz(-1.5895529) q[2];
sx q[2];
rz(-0.64964408) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7900438) q[1];
sx q[1];
rz(-2.1218532) q[1];
sx q[1];
rz(-0.95981564) q[1];
rz(-pi) q[2];
rz(1.8848041) q[3];
sx q[3];
rz(-1.2824138) q[3];
sx q[3];
rz(0.42321229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.39151829) q[2];
sx q[2];
rz(-1.8374279) q[2];
sx q[2];
rz(0.94775003) q[2];
rz(-0.1968955) q[3];
sx q[3];
rz(-0.24053776) q[3];
sx q[3];
rz(0.31015629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31203684) q[0];
sx q[0];
rz(-2.0289679) q[0];
sx q[0];
rz(2.6716676) q[0];
rz(1.6795109) q[1];
sx q[1];
rz(-1.3009678) q[1];
sx q[1];
rz(-1.7376815) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71588281) q[0];
sx q[0];
rz(-0.68001594) q[0];
sx q[0];
rz(-2.143857) q[0];
rz(-pi) q[1];
x q[1];
rz(2.933104) q[2];
sx q[2];
rz(-1.1013774) q[2];
sx q[2];
rz(-2.4698348) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.95937377) q[1];
sx q[1];
rz(-1.2101908) q[1];
sx q[1];
rz(2.1350767) q[1];
rz(-pi) q[2];
rz(-0.58340729) q[3];
sx q[3];
rz(-1.2800763) q[3];
sx q[3];
rz(0.099778508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5002084) q[2];
sx q[2];
rz(-2.310014) q[2];
sx q[2];
rz(-0.10759648) q[2];
rz(-2.522116) q[3];
sx q[3];
rz(-1.8457125) q[3];
sx q[3];
rz(-1.3639601) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7463995) q[0];
sx q[0];
rz(-1.0522319) q[0];
sx q[0];
rz(-2.8885544) q[0];
rz(3.053275) q[1];
sx q[1];
rz(-2.2679236) q[1];
sx q[1];
rz(-0.94892445) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2491566) q[0];
sx q[0];
rz(-1.8674441) q[0];
sx q[0];
rz(-0.15693024) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1429092) q[2];
sx q[2];
rz(-0.44369953) q[2];
sx q[2];
rz(0.86715172) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9946549) q[1];
sx q[1];
rz(-2.6587464) q[1];
sx q[1];
rz(-1.6864683) q[1];
rz(-pi) q[2];
rz(-2.3036495) q[3];
sx q[3];
rz(-2.4334014) q[3];
sx q[3];
rz(1.4670682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.39763149) q[2];
sx q[2];
rz(-0.88699114) q[2];
sx q[2];
rz(-0.31069791) q[2];
rz(-2.9001696) q[3];
sx q[3];
rz(-1.9390691) q[3];
sx q[3];
rz(-2.0588622) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0017515) q[0];
sx q[0];
rz(-2.7353291) q[0];
sx q[0];
rz(1.9061506) q[0];
rz(-2.6605117) q[1];
sx q[1];
rz(-0.95293871) q[1];
sx q[1];
rz(-1.7135886) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9426711) q[0];
sx q[0];
rz(-2.3301972) q[0];
sx q[0];
rz(-1.5862203) q[0];
x q[1];
rz(-0.88793869) q[2];
sx q[2];
rz(-0.35338923) q[2];
sx q[2];
rz(1.2438174) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.075632) q[1];
sx q[1];
rz(-1.8143225) q[1];
sx q[1];
rz(1.8357518) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8974182) q[3];
sx q[3];
rz(-2.2664968) q[3];
sx q[3];
rz(-0.43750924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1741751) q[2];
sx q[2];
rz(-2.316541) q[2];
sx q[2];
rz(-0.027912557) q[2];
rz(-0.86999718) q[3];
sx q[3];
rz(-2.3614707) q[3];
sx q[3];
rz(1.9546668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1160527) q[0];
sx q[0];
rz(-1.5929796) q[0];
sx q[0];
rz(1.0593587) q[0];
rz(1.4147883) q[1];
sx q[1];
rz(-1.6635514) q[1];
sx q[1];
rz(2.6639604) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85738289) q[0];
sx q[0];
rz(-2.3071831) q[0];
sx q[0];
rz(1.7718023) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5651836) q[2];
sx q[2];
rz(-2.8057348) q[2];
sx q[2];
rz(-2.4120021) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.79523008) q[1];
sx q[1];
rz(-1.6246965) q[1];
sx q[1];
rz(-1.4583711) q[1];
rz(2.2077363) q[3];
sx q[3];
rz(-2.2837486) q[3];
sx q[3];
rz(-1.0757004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16964218) q[2];
sx q[2];
rz(-2.1596238) q[2];
sx q[2];
rz(2.7726445) q[2];
rz(1.2667027) q[3];
sx q[3];
rz(-2.4167175) q[3];
sx q[3];
rz(2.1784311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0030768) q[0];
sx q[0];
rz(-1.930548) q[0];
sx q[0];
rz(1.7532274) q[0];
rz(2.6970462) q[1];
sx q[1];
rz(-0.81041705) q[1];
sx q[1];
rz(3.0141426) q[1];
rz(0.44345345) q[2];
sx q[2];
rz(-1.8691861) q[2];
sx q[2];
rz(-0.791823) q[2];
rz(-0.038158554) q[3];
sx q[3];
rz(-2.6737619) q[3];
sx q[3];
rz(3.0554042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-2.0353844) q[0];
sx q[0];
rz(-0.68618965) q[0];
sx q[0];
rz(1.9880779) q[0];
rz(-1.8312307) q[1];
sx q[1];
rz(-2.6481833) q[1];
sx q[1];
rz(0.44841132) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36487406) q[0];
sx q[0];
rz(-1.6239415) q[0];
sx q[0];
rz(2.4928635) q[0];
x q[1];
rz(1.6719867) q[2];
sx q[2];
rz(-1.3162344) q[2];
sx q[2];
rz(-2.2806666) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5858305) q[1];
sx q[1];
rz(-2.1735648) q[1];
sx q[1];
rz(0.65048154) q[1];
x q[2];
rz(0.5011933) q[3];
sx q[3];
rz(-1.5287011) q[3];
sx q[3];
rz(2.834109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51556921) q[2];
sx q[2];
rz(-1.3884576) q[2];
sx q[2];
rz(-2.5173729) q[2];
rz(0.37912399) q[3];
sx q[3];
rz(-0.99025327) q[3];
sx q[3];
rz(2.8930801) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9666331) q[0];
sx q[0];
rz(-0.81460726) q[0];
sx q[0];
rz(2.1061184) q[0];
rz(1.8677208) q[1];
sx q[1];
rz(-0.73718166) q[1];
sx q[1];
rz(0.48286352) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.14014) q[0];
sx q[0];
rz(-1.4743544) q[0];
sx q[0];
rz(-1.524855) q[0];
rz(-pi) q[1];
rz(-1.4522885) q[2];
sx q[2];
rz(-1.7359455) q[2];
sx q[2];
rz(-1.120795) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9804364) q[1];
sx q[1];
rz(-2.709678) q[1];
sx q[1];
rz(-1.5483787) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6522406) q[3];
sx q[3];
rz(-1.8907428) q[3];
sx q[3];
rz(-0.68405747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1172993) q[2];
sx q[2];
rz(-1.6397986) q[2];
sx q[2];
rz(-1.2705605) q[2];
rz(-2.471586) q[3];
sx q[3];
rz(-0.62205258) q[3];
sx q[3];
rz(2.2445934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73827493) q[0];
sx q[0];
rz(-0.44602317) q[0];
sx q[0];
rz(-0.73905149) q[0];
rz(0.96145472) q[1];
sx q[1];
rz(-2.8196204) q[1];
sx q[1];
rz(2.3846073) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2852427) q[0];
sx q[0];
rz(-0.77613554) q[0];
sx q[0];
rz(-1.1936643) q[0];
rz(-pi) q[1];
rz(-2.4503373) q[2];
sx q[2];
rz(-1.8658085) q[2];
sx q[2];
rz(-1.6694836) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0464335) q[1];
sx q[1];
rz(-2.1951402) q[1];
sx q[1];
rz(1.7364362) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2631031) q[3];
sx q[3];
rz(-1.6374) q[3];
sx q[3];
rz(0.25552017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5028533) q[2];
sx q[2];
rz(-2.2245378) q[2];
sx q[2];
rz(2.4364831) q[2];
rz(-0.32143587) q[3];
sx q[3];
rz(-1.0175846) q[3];
sx q[3];
rz(-2.7307935) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.131677) q[0];
sx q[0];
rz(-0.83808815) q[0];
sx q[0];
rz(-1.9158844) q[0];
rz(1.2340087) q[1];
sx q[1];
rz(-1.0154513) q[1];
sx q[1];
rz(-0.066224901) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4915337) q[0];
sx q[0];
rz(-0.26187944) q[0];
sx q[0];
rz(1.9240379) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8784857) q[2];
sx q[2];
rz(-2.3480847) q[2];
sx q[2];
rz(-2.6621278) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6844422) q[1];
sx q[1];
rz(-1.8891786) q[1];
sx q[1];
rz(1.063698) q[1];
x q[2];
rz(-1.6965167) q[3];
sx q[3];
rz(-1.2074405) q[3];
sx q[3];
rz(-1.1896127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0951198) q[2];
sx q[2];
rz(-0.9943704) q[2];
sx q[2];
rz(-2.7009916) q[2];
rz(-1.0062086) q[3];
sx q[3];
rz(-1.3738084) q[3];
sx q[3];
rz(-0.83427507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42136583) q[0];
sx q[0];
rz(-0.81291968) q[0];
sx q[0];
rz(1.5059858) q[0];
rz(1.6141363) q[1];
sx q[1];
rz(-1.9215877) q[1];
sx q[1];
rz(-1.45586) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9855749) q[0];
sx q[0];
rz(-0.7472207) q[0];
sx q[0];
rz(1.3559627) q[0];
rz(1.4541995) q[2];
sx q[2];
rz(-1.7966088) q[2];
sx q[2];
rz(1.9436398) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8511635) q[1];
sx q[1];
rz(-2.4896087) q[1];
sx q[1];
rz(0.033520582) q[1];
rz(-pi) q[2];
rz(0.60815706) q[3];
sx q[3];
rz(-1.6975612) q[3];
sx q[3];
rz(-0.34846445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27318925) q[2];
sx q[2];
rz(-1.5634147) q[2];
sx q[2];
rz(-2.2301162) q[2];
rz(0.8738001) q[3];
sx q[3];
rz(-0.74207145) q[3];
sx q[3];
rz(0.0066643683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28354302) q[0];
sx q[0];
rz(-2.8908505) q[0];
sx q[0];
rz(-0.13033303) q[0];
rz(2.6538972) q[1];
sx q[1];
rz(-2.6002488) q[1];
sx q[1];
rz(2.3977051) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.374342) q[0];
sx q[0];
rz(-3.0453186) q[0];
sx q[0];
rz(1.0262579) q[0];
x q[1];
rz(1.3792453) q[2];
sx q[2];
rz(-2.5073176) q[2];
sx q[2];
rz(0.77855643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6106657) q[1];
sx q[1];
rz(-1.0333038) q[1];
sx q[1];
rz(-2.8100138) q[1];
rz(2.8091431) q[3];
sx q[3];
rz(-0.5088734) q[3];
sx q[3];
rz(2.8116941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.1193715) q[2];
sx q[2];
rz(-0.85249844) q[2];
sx q[2];
rz(-2.1807097) q[2];
rz(-0.39572257) q[3];
sx q[3];
rz(-1.3362249) q[3];
sx q[3];
rz(3.0254288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6845067) q[0];
sx q[0];
rz(-1.1154024) q[0];
sx q[0];
rz(2.0945666) q[0];
rz(-0.60415769) q[1];
sx q[1];
rz(-1.2774757) q[1];
sx q[1];
rz(-2.7488757) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34907863) q[0];
sx q[0];
rz(-2.431708) q[0];
sx q[0];
rz(-1.841808) q[0];
rz(1.4899859) q[2];
sx q[2];
rz(-0.62056345) q[2];
sx q[2];
rz(-1.1309689) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8128374) q[1];
sx q[1];
rz(-0.27662524) q[1];
sx q[1];
rz(0.098621086) q[1];
rz(-pi) q[2];
rz(1.8363399) q[3];
sx q[3];
rz(-1.8887541) q[3];
sx q[3];
rz(-1.6443271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.087223209) q[2];
sx q[2];
rz(-1.1944218) q[2];
sx q[2];
rz(-1.3090022) q[2];
rz(-1.3537815) q[3];
sx q[3];
rz(-1.196922) q[3];
sx q[3];
rz(2.8712809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29276174) q[0];
sx q[0];
rz(-1.6852385) q[0];
sx q[0];
rz(-0.30717474) q[0];
rz(1.3245026) q[1];
sx q[1];
rz(-0.38506404) q[1];
sx q[1];
rz(2.2408392) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7867635) q[0];
sx q[0];
rz(-2.8788924) q[0];
sx q[0];
rz(3.0765615) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2239752) q[2];
sx q[2];
rz(-1.2003044) q[2];
sx q[2];
rz(2.8965184) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14054243) q[1];
sx q[1];
rz(-1.2549572) q[1];
sx q[1];
rz(1.2211826) q[1];
x q[2];
rz(1.4291271) q[3];
sx q[3];
rz(-0.9441388) q[3];
sx q[3];
rz(0.10654813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8073392) q[2];
sx q[2];
rz(-0.05143493) q[2];
sx q[2];
rz(2.5267498) q[2];
rz(1.0963415) q[3];
sx q[3];
rz(-0.59211007) q[3];
sx q[3];
rz(0.69826564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39356247) q[0];
sx q[0];
rz(-2.5108971) q[0];
sx q[0];
rz(-2.6322741) q[0];
rz(-0.18109426) q[1];
sx q[1];
rz(-1.4053922) q[1];
sx q[1];
rz(2.1720355) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90617243) q[0];
sx q[0];
rz(-2.3114822) q[0];
sx q[0];
rz(1.7801911) q[0];
x q[1];
rz(1.0636341) q[2];
sx q[2];
rz(-1.615287) q[2];
sx q[2];
rz(2.8756093) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6020376) q[1];
sx q[1];
rz(-1.4629211) q[1];
sx q[1];
rz(2.8007568) q[1];
rz(2.2990555) q[3];
sx q[3];
rz(-1.8368097) q[3];
sx q[3];
rz(2.0913948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8934882) q[2];
sx q[2];
rz(-0.23427811) q[2];
sx q[2];
rz(-1.0815557) q[2];
rz(2.9099416) q[3];
sx q[3];
rz(-1.7678429) q[3];
sx q[3];
rz(-2.933568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3807826) q[0];
sx q[0];
rz(-0.2247227) q[0];
sx q[0];
rz(2.494452) q[0];
rz(-2.6864247) q[1];
sx q[1];
rz(-1.7166694) q[1];
sx q[1];
rz(0.761935) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3477701) q[0];
sx q[0];
rz(-2.2552654) q[0];
sx q[0];
rz(-1.9467627) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64024957) q[2];
sx q[2];
rz(-0.92776042) q[2];
sx q[2];
rz(-1.6563479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0207723) q[1];
sx q[1];
rz(-1.6050395) q[1];
sx q[1];
rz(-0.65837752) q[1];
x q[2];
rz(-1.5860709) q[3];
sx q[3];
rz(-2.6741283) q[3];
sx q[3];
rz(-1.338442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.053293856) q[2];
sx q[2];
rz(-2.9247354) q[2];
sx q[2];
rz(-2.2376412) q[2];
rz(1.6186742) q[3];
sx q[3];
rz(-2.121033) q[3];
sx q[3];
rz(-1.1618377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91697964) q[0];
sx q[0];
rz(-1.1876748) q[0];
sx q[0];
rz(-1.351958) q[0];
rz(-3.0147973) q[1];
sx q[1];
rz(-2.3190111) q[1];
sx q[1];
rz(0.93228985) q[1];
rz(2.8895072) q[2];
sx q[2];
rz(-1.132195) q[2];
sx q[2];
rz(1.8762527) q[2];
rz(2.7115962) q[3];
sx q[3];
rz(-1.5734133) q[3];
sx q[3];
rz(-1.6071241) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

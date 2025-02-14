OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.10915169) q[0];
sx q[0];
rz(-1.9885795) q[0];
sx q[0];
rz(0.074540019) q[0];
rz(1.6115161) q[1];
sx q[1];
rz(3.2009701) q[1];
sx q[1];
rz(8.9006807) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68428129) q[0];
sx q[0];
rz(-1.0012828) q[0];
sx q[0];
rz(2.7049541) q[0];
rz(-pi) q[1];
rz(1.8139548) q[2];
sx q[2];
rz(-1.8983525) q[2];
sx q[2];
rz(2.067468) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9896696) q[1];
sx q[1];
rz(-1.0115342) q[1];
sx q[1];
rz(-0.86821809) q[1];
x q[2];
rz(-1.4788116) q[3];
sx q[3];
rz(-1.6114724) q[3];
sx q[3];
rz(-2.3856861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3268299) q[2];
sx q[2];
rz(-2.6430898) q[2];
sx q[2];
rz(-2.2461183) q[2];
rz(0.26252663) q[3];
sx q[3];
rz(-2.7956876) q[3];
sx q[3];
rz(3.053022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9240016) q[0];
sx q[0];
rz(-1.1955248) q[0];
sx q[0];
rz(-0.1161282) q[0];
rz(2.5092292) q[1];
sx q[1];
rz(-2.875681) q[1];
sx q[1];
rz(-1.4858474) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6054942) q[0];
sx q[0];
rz(-1.5935531) q[0];
sx q[0];
rz(3.0748585) q[0];
x q[1];
rz(2.9048278) q[2];
sx q[2];
rz(-2.1000266) q[2];
sx q[2];
rz(-2.6108612) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1097488) q[1];
sx q[1];
rz(-0.1662152) q[1];
sx q[1];
rz(-0.69352652) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13074517) q[3];
sx q[3];
rz(-0.8746038) q[3];
sx q[3];
rz(-1.8896904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.31132013) q[2];
sx q[2];
rz(-0.98036426) q[2];
sx q[2];
rz(0.92973989) q[2];
rz(2.7021507) q[3];
sx q[3];
rz(-1.5403055) q[3];
sx q[3];
rz(-0.75696993) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0759401) q[0];
sx q[0];
rz(-2.3891698) q[0];
sx q[0];
rz(-2.2546076) q[0];
rz(1.2607964) q[1];
sx q[1];
rz(-2.0340684) q[1];
sx q[1];
rz(-1.4874123) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5652147) q[0];
sx q[0];
rz(-1.6883255) q[0];
sx q[0];
rz(3.1410724) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6662798) q[2];
sx q[2];
rz(-1.1794568) q[2];
sx q[2];
rz(-3.0961159) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9650813) q[1];
sx q[1];
rz(-1.2413032) q[1];
sx q[1];
rz(-1.5894057) q[1];
x q[2];
rz(2.5599407) q[3];
sx q[3];
rz(-1.2949139) q[3];
sx q[3];
rz(1.1946071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4999353) q[2];
sx q[2];
rz(-0.96330088) q[2];
sx q[2];
rz(-0.43528834) q[2];
rz(0.38517243) q[3];
sx q[3];
rz(-1.9305072) q[3];
sx q[3];
rz(-0.96863532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4061072) q[0];
sx q[0];
rz(-3.0574953) q[0];
sx q[0];
rz(-0.054585833) q[0];
rz(-1.6361884) q[1];
sx q[1];
rz(-1.9278229) q[1];
sx q[1];
rz(-0.41753599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5180001) q[0];
sx q[0];
rz(-1.5125649) q[0];
sx q[0];
rz(-2.7043155) q[0];
rz(0.436385) q[2];
sx q[2];
rz(-2.4616995) q[2];
sx q[2];
rz(-1.5117548) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7723267) q[1];
sx q[1];
rz(-1.7290745) q[1];
sx q[1];
rz(-1.5637828) q[1];
rz(-1.3261111) q[3];
sx q[3];
rz(-0.26991329) q[3];
sx q[3];
rz(-3.1343366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2893082) q[2];
sx q[2];
rz(-2.4384629) q[2];
sx q[2];
rz(3.1404148) q[2];
rz(-1.3834472) q[3];
sx q[3];
rz(-0.26418424) q[3];
sx q[3];
rz(1.0435102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99172878) q[0];
sx q[0];
rz(-2.9809451) q[0];
sx q[0];
rz(0.37586656) q[0];
rz(2.1916892) q[1];
sx q[1];
rz(-1.4774731) q[1];
sx q[1];
rz(-2.4163767) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0501971) q[0];
sx q[0];
rz(-1.689347) q[0];
sx q[0];
rz(-1.4539745) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64230772) q[2];
sx q[2];
rz(-2.1644582) q[2];
sx q[2];
rz(0.35333179) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29595536) q[1];
sx q[1];
rz(-1.6682677) q[1];
sx q[1];
rz(-1.5649765) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8309202) q[3];
sx q[3];
rz(-0.88444607) q[3];
sx q[3];
rz(-2.7174866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3788562) q[2];
sx q[2];
rz(-1.6894268) q[2];
sx q[2];
rz(-0.47259304) q[2];
rz(0.11911123) q[3];
sx q[3];
rz(-2.5178858) q[3];
sx q[3];
rz(-2.3314893) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3177719) q[0];
sx q[0];
rz(-2.9404984) q[0];
sx q[0];
rz(-0.066545181) q[0];
rz(-2.9464974) q[1];
sx q[1];
rz(-2.0654443) q[1];
sx q[1];
rz(-1.9871575) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6789952) q[0];
sx q[0];
rz(-0.53946153) q[0];
sx q[0];
rz(-2.0360002) q[0];
rz(0.17627861) q[2];
sx q[2];
rz(-0.46277324) q[2];
sx q[2];
rz(1.1413934) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6007541) q[1];
sx q[1];
rz(-2.7620865) q[1];
sx q[1];
rz(1.862816) q[1];
rz(-1.5905753) q[3];
sx q[3];
rz(-0.58073509) q[3];
sx q[3];
rz(-1.6466717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.92516148) q[2];
sx q[2];
rz(-1.5387646) q[2];
sx q[2];
rz(2.4161941) q[2];
rz(1.409449) q[3];
sx q[3];
rz(-1.9603445) q[3];
sx q[3];
rz(-0.60666549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055939097) q[0];
sx q[0];
rz(-0.5760718) q[0];
sx q[0];
rz(2.2357909) q[0];
rz(-0.55465758) q[1];
sx q[1];
rz(-0.83811086) q[1];
sx q[1];
rz(2.99756) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5921123) q[0];
sx q[0];
rz(-2.0680111) q[0];
sx q[0];
rz(-0.29659941) q[0];
rz(2.0833182) q[2];
sx q[2];
rz(-1.4965759) q[2];
sx q[2];
rz(0.74736258) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1220458) q[1];
sx q[1];
rz(-0.48312369) q[1];
sx q[1];
rz(-0.93576099) q[1];
rz(2.1062076) q[3];
sx q[3];
rz(-0.28942063) q[3];
sx q[3];
rz(-1.1076224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4732699) q[2];
sx q[2];
rz(-0.39445764) q[2];
sx q[2];
rz(-2.1486166) q[2];
rz(2.9571577) q[3];
sx q[3];
rz(-0.95576972) q[3];
sx q[3];
rz(-2.1485645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011878012) q[0];
sx q[0];
rz(-2.5373902) q[0];
sx q[0];
rz(-2.7469444) q[0];
rz(-2.6036085) q[1];
sx q[1];
rz(-0.6901651) q[1];
sx q[1];
rz(-1.1916377) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79121548) q[0];
sx q[0];
rz(-2.1095304) q[0];
sx q[0];
rz(1.2722434) q[0];
rz(-1.3778119) q[2];
sx q[2];
rz(-1.9375083) q[2];
sx q[2];
rz(-0.16736469) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.059947) q[1];
sx q[1];
rz(-0.73127247) q[1];
sx q[1];
rz(-0.88986963) q[1];
x q[2];
rz(-1.1147145) q[3];
sx q[3];
rz(-0.9852162) q[3];
sx q[3];
rz(-2.8600313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3523606) q[2];
sx q[2];
rz(-1.5287986) q[2];
sx q[2];
rz(-0.94190502) q[2];
rz(2.3801129) q[3];
sx q[3];
rz(-2.0165636) q[3];
sx q[3];
rz(-0.040508125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86904675) q[0];
sx q[0];
rz(-1.3665552) q[0];
sx q[0];
rz(-1.0719365) q[0];
rz(2.5275224) q[1];
sx q[1];
rz(-2.2449988) q[1];
sx q[1];
rz(0.44642064) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20782411) q[0];
sx q[0];
rz(-0.75382262) q[0];
sx q[0];
rz(-1.1541744) q[0];
x q[1];
rz(-0.69009366) q[2];
sx q[2];
rz(-1.141667) q[2];
sx q[2];
rz(-0.60455017) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3371967) q[1];
sx q[1];
rz(-2.2257518) q[1];
sx q[1];
rz(1.7095997) q[1];
x q[2];
rz(-1.4037973) q[3];
sx q[3];
rz(-1.0713959) q[3];
sx q[3];
rz(1.7399519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.3093695) q[2];
sx q[2];
rz(-2.1840405) q[2];
sx q[2];
rz(1.2584125) q[2];
rz(1.9084515) q[3];
sx q[3];
rz(-2.5150053) q[3];
sx q[3];
rz(1.3174177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42245361) q[0];
sx q[0];
rz(-2.3261676) q[0];
sx q[0];
rz(-1.7183787) q[0];
rz(-2.8990959) q[1];
sx q[1];
rz(-2.6624694) q[1];
sx q[1];
rz(-1.2538145) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72287382) q[0];
sx q[0];
rz(-0.44298816) q[0];
sx q[0];
rz(1.1305869) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0331633) q[2];
sx q[2];
rz(-0.067010894) q[2];
sx q[2];
rz(-2.5052414) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4776626) q[1];
sx q[1];
rz(-0.60873181) q[1];
sx q[1];
rz(-2.1109646) q[1];
rz(-2.9217072) q[3];
sx q[3];
rz(-1.9337092) q[3];
sx q[3];
rz(0.21838926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8157876) q[2];
sx q[2];
rz(-0.82745224) q[2];
sx q[2];
rz(3.117756) q[2];
rz(1.2787974) q[3];
sx q[3];
rz(-0.33134225) q[3];
sx q[3];
rz(2.3648025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9780289) q[0];
sx q[0];
rz(-1.704287) q[0];
sx q[0];
rz(-1.1821672) q[0];
rz(0.72721807) q[1];
sx q[1];
rz(-1.4677508) q[1];
sx q[1];
rz(-0.8263091) q[1];
rz(-1.4197311) q[2];
sx q[2];
rz(-1.551566) q[2];
sx q[2];
rz(-2.8788064) q[2];
rz(0.047688382) q[3];
sx q[3];
rz(-1.449122) q[3];
sx q[3];
rz(2.4301152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

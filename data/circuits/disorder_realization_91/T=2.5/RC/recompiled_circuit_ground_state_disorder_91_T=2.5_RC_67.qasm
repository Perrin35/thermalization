OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8102326) q[0];
sx q[0];
rz(-1.8868089) q[0];
sx q[0];
rz(-0.64594185) q[0];
rz(1.6826001) q[1];
sx q[1];
rz(3.269722) q[1];
sx q[1];
rz(10.092957) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0872535) q[0];
sx q[0];
rz(-0.26832661) q[0];
sx q[0];
rz(-2.7032204) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92360381) q[2];
sx q[2];
rz(-1.3062813) q[2];
sx q[2];
rz(2.2326927) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1937374) q[1];
sx q[1];
rz(-1.1754218) q[1];
sx q[1];
rz(-2.1389524) q[1];
rz(-0.48623881) q[3];
sx q[3];
rz(-2.231153) q[3];
sx q[3];
rz(0.59729353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6231923) q[2];
sx q[2];
rz(-0.65541583) q[2];
sx q[2];
rz(2.0584959) q[2];
rz(-2.8862503) q[3];
sx q[3];
rz(-2.2446003) q[3];
sx q[3];
rz(-1.9820836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.979368) q[0];
sx q[0];
rz(-0.14142445) q[0];
sx q[0];
rz(-0.033578385) q[0];
rz(1.1383188) q[1];
sx q[1];
rz(-0.99457026) q[1];
sx q[1];
rz(-0.23695645) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19582307) q[0];
sx q[0];
rz(-1.7993712) q[0];
sx q[0];
rz(1.4287455) q[0];
x q[1];
rz(0.18743022) q[2];
sx q[2];
rz(-2.0399562) q[2];
sx q[2];
rz(-2.3694866) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.33920501) q[1];
sx q[1];
rz(-1.5453813) q[1];
sx q[1];
rz(-1.7632381) q[1];
rz(-pi) q[2];
rz(-0.44741607) q[3];
sx q[3];
rz(-2.1037115) q[3];
sx q[3];
rz(2.3428903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.529155) q[2];
sx q[2];
rz(-2.0945022) q[2];
sx q[2];
rz(-3.0212413) q[2];
rz(2.555661) q[3];
sx q[3];
rz(-1.7610995) q[3];
sx q[3];
rz(-1.2683292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.878207) q[0];
sx q[0];
rz(-2.1757941) q[0];
sx q[0];
rz(2.1534488) q[0];
rz(-1.490961) q[1];
sx q[1];
rz(-0.71690503) q[1];
sx q[1];
rz(1.5879226) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40000382) q[0];
sx q[0];
rz(-0.015946139) q[0];
sx q[0];
rz(3.0724597) q[0];
rz(-1.6397255) q[2];
sx q[2];
rz(-1.7048827) q[2];
sx q[2];
rz(0.20251911) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.37582477) q[1];
sx q[1];
rz(-1.8319523) q[1];
sx q[1];
rz(3.1118666) q[1];
x q[2];
rz(-0.6122784) q[3];
sx q[3];
rz(-2.4796072) q[3];
sx q[3];
rz(-0.16890165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3042018) q[2];
sx q[2];
rz(-1.3434255) q[2];
sx q[2];
rz(-0.070579441) q[2];
rz(2.6523759) q[3];
sx q[3];
rz(-0.990812) q[3];
sx q[3];
rz(0.14744559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1864784) q[0];
sx q[0];
rz(-0.69545737) q[0];
sx q[0];
rz(1.4008993) q[0];
rz(-2.9715624) q[1];
sx q[1];
rz(-1.731512) q[1];
sx q[1];
rz(-1.2331351) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9233386) q[0];
sx q[0];
rz(-1.6955351) q[0];
sx q[0];
rz(-1.1472923) q[0];
rz(2.4958359) q[2];
sx q[2];
rz(-2.3120572) q[2];
sx q[2];
rz(0.84596201) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.94205663) q[1];
sx q[1];
rz(-0.28008533) q[1];
sx q[1];
rz(-1.5166609) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5828269) q[3];
sx q[3];
rz(-2.0352484) q[3];
sx q[3];
rz(-2.8248276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3418545) q[2];
sx q[2];
rz(-1.7065455) q[2];
sx q[2];
rz(-0.62961659) q[2];
rz(-3.0047505) q[3];
sx q[3];
rz(-0.62961951) q[3];
sx q[3];
rz(1.7262044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28302309) q[0];
sx q[0];
rz(-0.50823277) q[0];
sx q[0];
rz(0.10230219) q[0];
rz(-1.6434068) q[1];
sx q[1];
rz(-1.3003277) q[1];
sx q[1];
rz(2.7844875) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84182318) q[0];
sx q[0];
rz(-1.7938611) q[0];
sx q[0];
rz(1.7905495) q[0];
rz(-3.018385) q[2];
sx q[2];
rz(-1.5438617) q[2];
sx q[2];
rz(-1.7719442) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6372924) q[1];
sx q[1];
rz(-1.8691113) q[1];
sx q[1];
rz(-0.53108414) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0856241) q[3];
sx q[3];
rz(-1.9686804) q[3];
sx q[3];
rz(1.1345989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.47306481) q[2];
sx q[2];
rz(-1.7230956) q[2];
sx q[2];
rz(-1.3225887) q[2];
rz(-1.4332917) q[3];
sx q[3];
rz(-1.8247484) q[3];
sx q[3];
rz(0.1639666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.598031) q[0];
sx q[0];
rz(-0.99995166) q[0];
sx q[0];
rz(-1.7033956) q[0];
rz(-1.2403129) q[1];
sx q[1];
rz(-0.65483171) q[1];
sx q[1];
rz(1.3528489) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.529218) q[0];
sx q[0];
rz(-1.1657526) q[0];
sx q[0];
rz(-0.10251001) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8565494) q[2];
sx q[2];
rz(-0.54423287) q[2];
sx q[2];
rz(-2.2533992) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.50033113) q[1];
sx q[1];
rz(-1.2029543) q[1];
sx q[1];
rz(-1.0235051) q[1];
rz(0.39522533) q[3];
sx q[3];
rz(-1.6754284) q[3];
sx q[3];
rz(2.6992309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.214434) q[2];
sx q[2];
rz(-0.36403251) q[2];
sx q[2];
rz(-0.39043179) q[2];
rz(-1.3123784) q[3];
sx q[3];
rz(-2.3305011) q[3];
sx q[3];
rz(2.3217679) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.696233) q[0];
sx q[0];
rz(-2.2177028) q[0];
sx q[0];
rz(-1.4129289) q[0];
rz(1.6784809) q[1];
sx q[1];
rz(-1.5556346) q[1];
sx q[1];
rz(2.5020592) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7752658) q[0];
sx q[0];
rz(-1.9444939) q[0];
sx q[0];
rz(2.213272) q[0];
x q[1];
rz(-2.9908871) q[2];
sx q[2];
rz(-1.6455212) q[2];
sx q[2];
rz(-3.0091803) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23375854) q[1];
sx q[1];
rz(-2.4423264) q[1];
sx q[1];
rz(-2.722867) q[1];
x q[2];
rz(-1.9652548) q[3];
sx q[3];
rz(-1.6812075) q[3];
sx q[3];
rz(-0.73964707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.21706906) q[2];
sx q[2];
rz(-1.8127433) q[2];
sx q[2];
rz(-2.6893943) q[2];
rz(-2.9186115) q[3];
sx q[3];
rz(-2.860234) q[3];
sx q[3];
rz(-2.9862459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0105791) q[0];
sx q[0];
rz(-0.92229811) q[0];
sx q[0];
rz(0.16192326) q[0];
rz(2.7936392) q[1];
sx q[1];
rz(-1.8674928) q[1];
sx q[1];
rz(-0.23439342) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9471301) q[0];
sx q[0];
rz(-0.85887733) q[0];
sx q[0];
rz(-2.3361337) q[0];
rz(1.2599808) q[2];
sx q[2];
rz(-0.87500415) q[2];
sx q[2];
rz(2.0134848) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1474851) q[1];
sx q[1];
rz(-0.6692769) q[1];
sx q[1];
rz(1.2662751) q[1];
x q[2];
rz(0.23080821) q[3];
sx q[3];
rz(-0.2838906) q[3];
sx q[3];
rz(0.47286716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.28635412) q[2];
sx q[2];
rz(-2.0704634) q[2];
sx q[2];
rz(0.62412778) q[2];
rz(0.74328077) q[3];
sx q[3];
rz(-0.99370876) q[3];
sx q[3];
rz(0.68964094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2766992) q[0];
sx q[0];
rz(-1.6962637) q[0];
sx q[0];
rz(-2.0546761) q[0];
rz(-1.2326321) q[1];
sx q[1];
rz(-1.107736) q[1];
sx q[1];
rz(1.3023652) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24827458) q[0];
sx q[0];
rz(-2.9030307) q[0];
sx q[0];
rz(-0.49805157) q[0];
x q[1];
rz(2.0696569) q[2];
sx q[2];
rz(-2.4833224) q[2];
sx q[2];
rz(1.2128304) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2757146) q[1];
sx q[1];
rz(-0.68531407) q[1];
sx q[1];
rz(-2.4908743) q[1];
rz(2.5185561) q[3];
sx q[3];
rz(-1.5554232) q[3];
sx q[3];
rz(-0.29237177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6752601) q[2];
sx q[2];
rz(-1.2559428) q[2];
sx q[2];
rz(2.7272398) q[2];
rz(0.77053344) q[3];
sx q[3];
rz(-1.4221752) q[3];
sx q[3];
rz(0.93366247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3122124) q[0];
sx q[0];
rz(-2.3214564) q[0];
sx q[0];
rz(2.6729551) q[0];
rz(-1.8679484) q[1];
sx q[1];
rz(-1.2761152) q[1];
sx q[1];
rz(0.31148568) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2149406) q[0];
sx q[0];
rz(-2.4209341) q[0];
sx q[0];
rz(-1.65833) q[0];
rz(-1.6577254) q[2];
sx q[2];
rz(-1.6875899) q[2];
sx q[2];
rz(-0.96125666) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2356687) q[1];
sx q[1];
rz(-2.0579268) q[1];
sx q[1];
rz(3.0887927) q[1];
rz(2.4021637) q[3];
sx q[3];
rz(-0.66777705) q[3];
sx q[3];
rz(-1.2583994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1062539) q[2];
sx q[2];
rz(-1.1941348) q[2];
sx q[2];
rz(-1.7485471) q[2];
rz(-0.92750183) q[3];
sx q[3];
rz(-1.5774957) q[3];
sx q[3];
rz(-0.86021304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8504234) q[0];
sx q[0];
rz(-2.8130154) q[0];
sx q[0];
rz(-1.6541506) q[0];
rz(-0.18118478) q[1];
sx q[1];
rz(-2.1962427) q[1];
sx q[1];
rz(-0.93999351) q[1];
rz(2.9857582) q[2];
sx q[2];
rz(-2.0991785) q[2];
sx q[2];
rz(3.0002158) q[2];
rz(0.90418935) q[3];
sx q[3];
rz(-2.0695569) q[3];
sx q[3];
rz(0.29026779) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

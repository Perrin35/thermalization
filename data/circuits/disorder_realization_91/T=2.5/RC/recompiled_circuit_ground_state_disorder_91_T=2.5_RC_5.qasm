OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.33136) q[0];
sx q[0];
rz(-1.2547837) q[0];
sx q[0];
rz(0.64594185) q[0];
rz(-1.4589925) q[1];
sx q[1];
rz(-0.12812935) q[1];
sx q[1];
rz(2.473414) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9409882) q[0];
sx q[0];
rz(-1.683569) q[0];
sx q[0];
rz(-0.24399816) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8142848) q[2];
sx q[2];
rz(-2.1919554) q[2];
sx q[2];
rz(0.46681625) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5234452) q[1];
sx q[1];
rz(-1.0511569) q[1];
sx q[1];
rz(2.6818399) q[1];
x q[2];
rz(0.84998895) q[3];
sx q[3];
rz(-1.1927943) q[3];
sx q[3];
rz(1.8545618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6231923) q[2];
sx q[2];
rz(-0.65541583) q[2];
sx q[2];
rz(-2.0584959) q[2];
rz(0.25534233) q[3];
sx q[3];
rz(-2.2446003) q[3];
sx q[3];
rz(1.159509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1622247) q[0];
sx q[0];
rz(-3.0001682) q[0];
sx q[0];
rz(0.033578385) q[0];
rz(-1.1383188) q[1];
sx q[1];
rz(-0.99457026) q[1];
sx q[1];
rz(-2.9046362) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9457696) q[0];
sx q[0];
rz(-1.3422215) q[0];
sx q[0];
rz(-1.4287455) q[0];
rz(1.218538) q[2];
sx q[2];
rz(-0.50261231) q[2];
sx q[2];
rz(-1.1692804) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9050497) q[1];
sx q[1];
rz(-1.3784175) q[1];
sx q[1];
rz(3.1156999) q[1];
rz(-2.2036425) q[3];
sx q[3];
rz(-0.68162912) q[3];
sx q[3];
rz(-1.5555895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.529155) q[2];
sx q[2];
rz(-2.0945022) q[2];
sx q[2];
rz(3.0212413) q[2];
rz(-0.58593166) q[3];
sx q[3];
rz(-1.3804932) q[3];
sx q[3];
rz(-1.8732635) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.878207) q[0];
sx q[0];
rz(-0.96579856) q[0];
sx q[0];
rz(-2.1534488) q[0];
rz(1.6506317) q[1];
sx q[1];
rz(-2.4246876) q[1];
sx q[1];
rz(-1.5879226) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33086209) q[0];
sx q[0];
rz(-1.5867044) q[0];
sx q[0];
rz(-1.5718979) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13440172) q[2];
sx q[2];
rz(-1.5024868) q[2];
sx q[2];
rz(-1.7640863) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26116672) q[1];
sx q[1];
rz(-0.26280394) q[1];
sx q[1];
rz(-1.4600423) q[1];
rz(-pi) q[2];
rz(1.1497028) q[3];
sx q[3];
rz(-2.097887) q[3];
sx q[3];
rz(0.89638174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3042018) q[2];
sx q[2];
rz(-1.3434255) q[2];
sx q[2];
rz(-3.0710132) q[2];
rz(0.48921674) q[3];
sx q[3];
rz(-2.1507806) q[3];
sx q[3];
rz(-2.9941471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9551142) q[0];
sx q[0];
rz(-0.69545737) q[0];
sx q[0];
rz(-1.4008993) q[0];
rz(2.9715624) q[1];
sx q[1];
rz(-1.4100807) q[1];
sx q[1];
rz(-1.2331351) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9233386) q[0];
sx q[0];
rz(-1.6955351) q[0];
sx q[0];
rz(1.1472923) q[0];
rz(-pi) q[1];
rz(-0.64575671) q[2];
sx q[2];
rz(-2.3120572) q[2];
sx q[2];
rz(0.84596201) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68076949) q[1];
sx q[1];
rz(-1.5857547) q[1];
sx q[1];
rz(1.2911002) q[1];
rz(-pi) q[2];
rz(-1.5828269) q[3];
sx q[3];
rz(-2.0352484) q[3];
sx q[3];
rz(0.31676502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.79973811) q[2];
sx q[2];
rz(-1.4350472) q[2];
sx q[2];
rz(2.5119761) q[2];
rz(-3.0047505) q[3];
sx q[3];
rz(-0.62961951) q[3];
sx q[3];
rz(-1.4153882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28302309) q[0];
sx q[0];
rz(-0.50823277) q[0];
sx q[0];
rz(-0.10230219) q[0];
rz(-1.4981859) q[1];
sx q[1];
rz(-1.3003277) q[1];
sx q[1];
rz(-2.7844875) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67960167) q[0];
sx q[0];
rz(-1.7850189) q[0];
sx q[0];
rz(2.9132183) q[0];
rz(-pi) q[1];
rz(-0.12320766) q[2];
sx q[2];
rz(-1.5977309) q[2];
sx q[2];
rz(-1.7719442) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.60254492) q[1];
sx q[1];
rz(-2.5395406) q[1];
sx q[1];
rz(-2.5959488) q[1];
x q[2];
rz(2.0559685) q[3];
sx q[3];
rz(-1.1729122) q[3];
sx q[3];
rz(-1.1345989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.47306481) q[2];
sx q[2];
rz(-1.7230956) q[2];
sx q[2];
rz(-1.8190039) q[2];
rz(1.4332917) q[3];
sx q[3];
rz(-1.8247484) q[3];
sx q[3];
rz(2.9776261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54356164) q[0];
sx q[0];
rz(-2.141641) q[0];
sx q[0];
rz(-1.438197) q[0];
rz(-1.9012798) q[1];
sx q[1];
rz(-0.65483171) q[1];
sx q[1];
rz(1.7887438) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677296) q[0];
sx q[0];
rz(-2.7244748) q[0];
sx q[0];
rz(-1.8050844) q[0];
x q[1];
rz(-2.8565494) q[2];
sx q[2];
rz(-2.5973598) q[2];
sx q[2];
rz(0.88819347) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2861915) q[1];
sx q[1];
rz(-2.0778065) q[1];
sx q[1];
rz(-0.42393522) q[1];
rz(-2.7463673) q[3];
sx q[3];
rz(-1.4661643) q[3];
sx q[3];
rz(0.44236174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9271586) q[2];
sx q[2];
rz(-0.36403251) q[2];
sx q[2];
rz(0.39043179) q[2];
rz(-1.8292142) q[3];
sx q[3];
rz(-0.81109154) q[3];
sx q[3];
rz(-0.81982476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.696233) q[0];
sx q[0];
rz(-0.92388988) q[0];
sx q[0];
rz(1.7286638) q[0];
rz(-1.4631118) q[1];
sx q[1];
rz(-1.585958) q[1];
sx q[1];
rz(0.63953343) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93776751) q[0];
sx q[0];
rz(-0.97903189) q[0];
sx q[0];
rz(-0.45543619) q[0];
rz(-pi) q[1];
rz(-0.15070559) q[2];
sx q[2];
rz(-1.4960714) q[2];
sx q[2];
rz(0.13241235) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8486316) q[1];
sx q[1];
rz(-2.1994414) q[1];
sx q[1];
rz(1.2413003) q[1];
rz(-1.1763379) q[3];
sx q[3];
rz(-1.6812075) q[3];
sx q[3];
rz(0.73964707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.21706906) q[2];
sx q[2];
rz(-1.8127433) q[2];
sx q[2];
rz(-2.6893943) q[2];
rz(-0.2229812) q[3];
sx q[3];
rz(-0.28135869) q[3];
sx q[3];
rz(0.15534672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0105791) q[0];
sx q[0];
rz(-0.92229811) q[0];
sx q[0];
rz(-0.16192326) q[0];
rz(-2.7936392) q[1];
sx q[1];
rz(-1.8674928) q[1];
sx q[1];
rz(-2.9071992) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9471301) q[0];
sx q[0];
rz(-2.2827153) q[0];
sx q[0];
rz(-0.80545896) q[0];
x q[1];
rz(-2.7905445) q[2];
sx q[2];
rz(-2.3902438) q[2];
sx q[2];
rz(0.66354942) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.99410759) q[1];
sx q[1];
rz(-2.4723158) q[1];
sx q[1];
rz(1.2662751) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23080821) q[3];
sx q[3];
rz(-0.2838906) q[3];
sx q[3];
rz(-0.47286716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28635412) q[2];
sx q[2];
rz(-2.0704634) q[2];
sx q[2];
rz(-2.5174649) q[2];
rz(-0.74328077) q[3];
sx q[3];
rz(-0.99370876) q[3];
sx q[3];
rz(2.4519517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8648935) q[0];
sx q[0];
rz(-1.6962637) q[0];
sx q[0];
rz(1.0869166) q[0];
rz(1.2326321) q[1];
sx q[1];
rz(-2.0338567) q[1];
sx q[1];
rz(-1.8392275) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3051599) q[0];
sx q[0];
rz(-1.4576685) q[0];
sx q[0];
rz(0.21048429) q[0];
rz(-1.0719358) q[2];
sx q[2];
rz(-0.65827024) q[2];
sx q[2];
rz(-1.2128304) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2375396) q[1];
sx q[1];
rz(-1.1773279) q[1];
sx q[1];
rz(-2.5649125) q[1];
rz(-pi) q[2];
rz(-2.5185561) q[3];
sx q[3];
rz(-1.5554232) q[3];
sx q[3];
rz(0.29237177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46633259) q[2];
sx q[2];
rz(-1.8856498) q[2];
sx q[2];
rz(-2.7272398) q[2];
rz(-0.77053344) q[3];
sx q[3];
rz(-1.4221752) q[3];
sx q[3];
rz(-0.93366247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3122124) q[0];
sx q[0];
rz(-0.82013622) q[0];
sx q[0];
rz(0.46863753) q[0];
rz(-1.8679484) q[1];
sx q[1];
rz(-1.8654774) q[1];
sx q[1];
rz(2.830107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92665206) q[0];
sx q[0];
rz(-2.4209341) q[0];
sx q[0];
rz(1.4832627) q[0];
rz(-0.63705541) q[2];
sx q[2];
rz(-2.9961176) q[2];
sx q[2];
rz(1.6033974) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1232417) q[1];
sx q[1];
rz(-0.48975485) q[1];
sx q[1];
rz(1.4714929) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0592779) q[3];
sx q[3];
rz(-2.0460108) q[3];
sx q[3];
rz(2.7431874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0353388) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2911693) q[0];
sx q[0];
rz(-2.8130154) q[0];
sx q[0];
rz(-1.6541506) q[0];
rz(-2.9604079) q[1];
sx q[1];
rz(-0.94534992) q[1];
sx q[1];
rz(2.2015991) q[1];
rz(-0.15583444) q[2];
sx q[2];
rz(-2.0991785) q[2];
sx q[2];
rz(3.0002158) q[2];
rz(-0.84862205) q[3];
sx q[3];
rz(-0.80905882) q[3];
sx q[3];
rz(-0.73425135) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

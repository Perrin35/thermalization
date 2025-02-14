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
rz(-2.4956508) q[0];
rz(-1.4589925) q[1];
sx q[1];
rz(-0.12812935) q[1];
sx q[1];
rz(2.473414) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7433919) q[0];
sx q[0];
rz(-1.3283793) q[0];
sx q[0];
rz(1.4546118) q[0];
rz(1.1485841) q[2];
sx q[2];
rz(-0.6919043) q[2];
sx q[2];
rz(-0.9949323) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.61814743) q[1];
sx q[1];
rz(-2.0904358) q[1];
sx q[1];
rz(0.45975273) q[1];
rz(-2.6553538) q[3];
sx q[3];
rz(-0.91043962) q[3];
sx q[3];
rz(-2.5442991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6231923) q[2];
sx q[2];
rz(-2.4861768) q[2];
sx q[2];
rz(2.0584959) q[2];
rz(-2.8862503) q[3];
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
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.979368) q[0];
sx q[0];
rz(-0.14142445) q[0];
sx q[0];
rz(0.033578385) q[0];
rz(-1.1383188) q[1];
sx q[1];
rz(-0.99457026) q[1];
sx q[1];
rz(0.23695645) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7743999) q[0];
sx q[0];
rz(-0.26846186) q[0];
sx q[0];
rz(2.5949097) q[0];
rz(-pi) q[1];
rz(1.0944738) q[2];
sx q[2];
rz(-1.4038205) q[2];
sx q[2];
rz(0.71314883) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33920501) q[1];
sx q[1];
rz(-1.5453813) q[1];
sx q[1];
rz(-1.3783546) q[1];
rz(-2.2036425) q[3];
sx q[3];
rz(-0.68162912) q[3];
sx q[3];
rz(1.5860032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.61243764) q[2];
sx q[2];
rz(-2.0945022) q[2];
sx q[2];
rz(-0.1203514) q[2];
rz(2.555661) q[3];
sx q[3];
rz(-1.7610995) q[3];
sx q[3];
rz(-1.2683292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2633857) q[0];
sx q[0];
rz(-2.1757941) q[0];
sx q[0];
rz(0.98814386) q[0];
rz(1.490961) q[1];
sx q[1];
rz(-2.4246876) q[1];
sx q[1];
rz(1.5879226) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2399167) q[0];
sx q[0];
rz(-1.5718978) q[0];
sx q[0];
rz(3.1256846) q[0];
x q[1];
rz(-1.5018671) q[2];
sx q[2];
rz(-1.7048827) q[2];
sx q[2];
rz(2.9390735) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7657679) q[1];
sx q[1];
rz(-1.8319523) q[1];
sx q[1];
rz(3.1118666) q[1];
rz(-1.1497028) q[3];
sx q[3];
rz(-1.0437056) q[3];
sx q[3];
rz(-2.2452109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1864784) q[0];
sx q[0];
rz(-2.4461353) q[0];
sx q[0];
rz(1.4008993) q[0];
rz(-0.1700302) q[1];
sx q[1];
rz(-1.4100807) q[1];
sx q[1];
rz(-1.2331351) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2182541) q[0];
sx q[0];
rz(-1.6955351) q[0];
sx q[0];
rz(1.9943004) q[0];
rz(-2.1523684) q[2];
sx q[2];
rz(-2.2007341) q[2];
sx q[2];
rz(3.1357855) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.68076949) q[1];
sx q[1];
rz(-1.5857547) q[1];
sx q[1];
rz(-1.2911002) q[1];
rz(-1.5828269) q[3];
sx q[3];
rz(-1.1063442) q[3];
sx q[3];
rz(2.8248276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3418545) q[2];
sx q[2];
rz(-1.4350472) q[2];
sx q[2];
rz(2.5119761) q[2];
rz(-0.13684212) q[3];
sx q[3];
rz(-0.62961951) q[3];
sx q[3];
rz(1.4153882) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28302309) q[0];
sx q[0];
rz(-2.6333599) q[0];
sx q[0];
rz(-0.10230219) q[0];
rz(1.6434068) q[1];
sx q[1];
rz(-1.841265) q[1];
sx q[1];
rz(-0.35710517) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.461991) q[0];
sx q[0];
rz(-1.3565738) q[0];
sx q[0];
rz(-0.22837436) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9257878) q[2];
sx q[2];
rz(-0.12610283) q[2];
sx q[2];
rz(-0.012995586) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50430027) q[1];
sx q[1];
rz(-1.8691113) q[1];
sx q[1];
rz(2.6105085) q[1];
rz(-pi) q[2];
rz(0.83729845) q[3];
sx q[3];
rz(-2.5243763) q[3];
sx q[3];
rz(-2.0715947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6685278) q[2];
sx q[2];
rz(-1.7230956) q[2];
sx q[2];
rz(-1.8190039) q[2];
rz(-1.708301) q[3];
sx q[3];
rz(-1.8247484) q[3];
sx q[3];
rz(2.9776261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54356164) q[0];
sx q[0];
rz(-2.141641) q[0];
sx q[0];
rz(-1.7033956) q[0];
rz(1.9012798) q[1];
sx q[1];
rz(-0.65483171) q[1];
sx q[1];
rz(-1.7887438) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0594994) q[0];
sx q[0];
rz(-1.4766066) q[0];
sx q[0];
rz(-1.1638428) q[0];
rz(-0.2850432) q[2];
sx q[2];
rz(-2.5973598) q[2];
sx q[2];
rz(-0.88819347) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2861915) q[1];
sx q[1];
rz(-2.0778065) q[1];
sx q[1];
rz(0.42393522) q[1];
x q[2];
rz(2.7463673) q[3];
sx q[3];
rz(-1.6754284) q[3];
sx q[3];
rz(0.44236174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.214434) q[2];
sx q[2];
rz(-2.7775601) q[2];
sx q[2];
rz(-0.39043179) q[2];
rz(1.3123784) q[3];
sx q[3];
rz(-0.81109154) q[3];
sx q[3];
rz(2.3217679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
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
rz(-1.585958) q[1];
sx q[1];
rz(0.63953343) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36632682) q[0];
sx q[0];
rz(-1.9444939) q[0];
sx q[0];
rz(0.92832066) q[0];
rz(-0.15070559) q[2];
sx q[2];
rz(-1.6455212) q[2];
sx q[2];
rz(3.0091803) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.29296103) q[1];
sx q[1];
rz(-0.94215122) q[1];
sx q[1];
rz(1.2413003) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11951167) q[3];
sx q[3];
rz(-1.1788713) q[3];
sx q[3];
rz(2.3562795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9245236) q[2];
sx q[2];
rz(-1.8127433) q[2];
sx q[2];
rz(0.45219839) q[2];
rz(-0.2229812) q[3];
sx q[3];
rz(-0.28135869) q[3];
sx q[3];
rz(0.15534672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-3.0105791) q[0];
sx q[0];
rz(-0.92229811) q[0];
sx q[0];
rz(2.9796694) q[0];
rz(2.7936392) q[1];
sx q[1];
rz(-1.8674928) q[1];
sx q[1];
rz(2.9071992) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1944626) q[0];
sx q[0];
rz(-0.85887733) q[0];
sx q[0];
rz(-0.80545896) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72004628) q[2];
sx q[2];
rz(-1.8077501) q[2];
sx q[2];
rz(-2.4958378) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99410759) q[1];
sx q[1];
rz(-0.6692769) q[1];
sx q[1];
rz(1.8753176) q[1];
rz(-0.23080821) q[3];
sx q[3];
rz(-0.2838906) q[3];
sx q[3];
rz(-0.47286716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8552385) q[2];
sx q[2];
rz(-2.0704634) q[2];
sx q[2];
rz(0.62412778) q[2];
rz(2.3983119) q[3];
sx q[3];
rz(-0.99370876) q[3];
sx q[3];
rz(2.4519517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2766992) q[0];
sx q[0];
rz(-1.445329) q[0];
sx q[0];
rz(1.0869166) q[0];
rz(1.9089606) q[1];
sx q[1];
rz(-1.107736) q[1];
sx q[1];
rz(1.3023652) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83643276) q[0];
sx q[0];
rz(-1.4576685) q[0];
sx q[0];
rz(0.21048429) q[0];
rz(-0.97424284) q[2];
sx q[2];
rz(-1.2737717) q[2];
sx q[2];
rz(3.0926306) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2757146) q[1];
sx q[1];
rz(-0.68531407) q[1];
sx q[1];
rz(0.65071836) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5518673) q[3];
sx q[3];
rz(-2.193748) q[3];
sx q[3];
rz(-1.8521233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6752601) q[2];
sx q[2];
rz(-1.8856498) q[2];
sx q[2];
rz(-2.7272398) q[2];
rz(0.77053344) q[3];
sx q[3];
rz(-1.7194175) q[3];
sx q[3];
rz(-0.93366247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82938021) q[0];
sx q[0];
rz(-2.3214564) q[0];
sx q[0];
rz(-2.6729551) q[0];
rz(-1.2736443) q[1];
sx q[1];
rz(-1.2761152) q[1];
sx q[1];
rz(2.830107) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5783008) q[0];
sx q[0];
rz(-1.5130763) q[0];
sx q[0];
rz(2.2895534) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1172322) q[2];
sx q[2];
rz(-1.4844609) q[2];
sx q[2];
rz(-2.5422079) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.905924) q[1];
sx q[1];
rz(-1.0836658) q[1];
sx q[1];
rz(0.052799932) q[1];
rz(-pi) q[2];
rz(1.0823147) q[3];
sx q[3];
rz(-2.0460108) q[3];
sx q[3];
rz(2.7431874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1062539) q[2];
sx q[2];
rz(-1.1941348) q[2];
sx q[2];
rz(-1.7485471) q[2];
rz(-0.92750183) q[3];
sx q[3];
rz(-1.564097) q[3];
sx q[3];
rz(-2.2813796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.8504234) q[0];
sx q[0];
rz(-0.32857729) q[0];
sx q[0];
rz(1.487442) q[0];
rz(-2.9604079) q[1];
sx q[1];
rz(-0.94534992) q[1];
sx q[1];
rz(2.2015991) q[1];
rz(-2.9857582) q[2];
sx q[2];
rz(-1.0424141) q[2];
sx q[2];
rz(-0.14137683) q[2];
rz(2.5355382) q[3];
sx q[3];
rz(-0.99679508) q[3];
sx q[3];
rz(1.5011277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

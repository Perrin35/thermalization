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
rz(1.6826001) q[1];
sx q[1];
rz(3.269722) q[1];
sx q[1];
rz(10.092957) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9409882) q[0];
sx q[0];
rz(-1.683569) q[0];
sx q[0];
rz(0.24399816) q[0];
rz(-pi) q[1];
rz(0.92360381) q[2];
sx q[2];
rz(-1.8353113) q[2];
sx q[2];
rz(-0.90889999) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5234452) q[1];
sx q[1];
rz(-2.0904358) q[1];
sx q[1];
rz(-2.6818399) q[1];
rz(-pi) q[2];
rz(-2.6553538) q[3];
sx q[3];
rz(-0.91043962) q[3];
sx q[3];
rz(0.59729353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6231923) q[2];
sx q[2];
rz(-0.65541583) q[2];
sx q[2];
rz(-1.0830967) q[2];
rz(-2.8862503) q[3];
sx q[3];
rz(-0.89699236) q[3];
sx q[3];
rz(-1.159509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1622247) q[0];
sx q[0];
rz(-0.14142445) q[0];
sx q[0];
rz(-0.033578385) q[0];
rz(1.1383188) q[1];
sx q[1];
rz(-0.99457026) q[1];
sx q[1];
rz(-0.23695645) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4073674) q[0];
sx q[0];
rz(-1.7091284) q[0];
sx q[0];
rz(0.23081918) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0471188) q[2];
sx q[2];
rz(-1.4038205) q[2];
sx q[2];
rz(2.4284438) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8023876) q[1];
sx q[1];
rz(-1.5453813) q[1];
sx q[1];
rz(-1.7632381) q[1];
rz(-pi) q[2];
rz(0.99144793) q[3];
sx q[3];
rz(-1.1889403) q[3];
sx q[3];
rz(0.53296158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61243764) q[2];
sx q[2];
rz(-2.0945022) q[2];
sx q[2];
rz(0.1203514) q[2];
rz(-0.58593166) q[3];
sx q[3];
rz(-1.3804932) q[3];
sx q[3];
rz(1.2683292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.878207) q[0];
sx q[0];
rz(-0.96579856) q[0];
sx q[0];
rz(-2.1534488) q[0];
rz(1.490961) q[1];
sx q[1];
rz(-2.4246876) q[1];
sx q[1];
rz(-1.5536701) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33086209) q[0];
sx q[0];
rz(-1.5548883) q[0];
sx q[0];
rz(-1.5718979) q[0];
rz(-pi) q[1];
rz(-0.13440172) q[2];
sx q[2];
rz(-1.5024868) q[2];
sx q[2];
rz(1.7640863) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7657679) q[1];
sx q[1];
rz(-1.8319523) q[1];
sx q[1];
rz(3.1118666) q[1];
rz(-pi) q[2];
rz(-0.6122784) q[3];
sx q[3];
rz(-2.4796072) q[3];
sx q[3];
rz(-0.16890165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.83739088) q[2];
sx q[2];
rz(-1.3434255) q[2];
sx q[2];
rz(-3.0710132) q[2];
rz(-0.48921674) q[3];
sx q[3];
rz(-2.1507806) q[3];
sx q[3];
rz(-0.14744559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1864784) q[0];
sx q[0];
rz(-2.4461353) q[0];
sx q[0];
rz(-1.4008993) q[0];
rz(-2.9715624) q[1];
sx q[1];
rz(-1.4100807) q[1];
sx q[1];
rz(1.2331351) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2182541) q[0];
sx q[0];
rz(-1.4460576) q[0];
sx q[0];
rz(1.1472923) q[0];
x q[1];
rz(0.64575671) q[2];
sx q[2];
rz(-2.3120572) q[2];
sx q[2];
rz(2.2956306) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.255862) q[1];
sx q[1];
rz(-1.8504603) q[1];
sx q[1];
rz(-0.015563029) q[1];
rz(-pi) q[2];
rz(-0.024007576) q[3];
sx q[3];
rz(-2.6769961) q[3];
sx q[3];
rz(-0.28991227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79973811) q[2];
sx q[2];
rz(-1.7065455) q[2];
sx q[2];
rz(0.62961659) q[2];
rz(-0.13684212) q[3];
sx q[3];
rz(-2.5119731) q[3];
sx q[3];
rz(-1.4153882) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8585696) q[0];
sx q[0];
rz(-0.50823277) q[0];
sx q[0];
rz(-0.10230219) q[0];
rz(1.4981859) q[1];
sx q[1];
rz(-1.3003277) q[1];
sx q[1];
rz(-0.35710517) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5095561) q[0];
sx q[0];
rz(-2.8297544) q[0];
sx q[0];
rz(-2.3760892) q[0];
x q[1];
rz(3.018385) q[2];
sx q[2];
rz(-1.5438617) q[2];
sx q[2];
rz(-1.3696485) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5390477) q[1];
sx q[1];
rz(-0.60205205) q[1];
sx q[1];
rz(0.54564387) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44355794) q[3];
sx q[3];
rz(-1.1263811) q[3];
sx q[3];
rz(-0.23469521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6685278) q[2];
sx q[2];
rz(-1.7230956) q[2];
sx q[2];
rz(-1.8190039) q[2];
rz(1.708301) q[3];
sx q[3];
rz(-1.3168443) q[3];
sx q[3];
rz(2.9776261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54356164) q[0];
sx q[0];
rz(-0.99995166) q[0];
sx q[0];
rz(-1.438197) q[0];
rz(1.9012798) q[1];
sx q[1];
rz(-2.4867609) q[1];
sx q[1];
rz(-1.3528489) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2738631) q[0];
sx q[0];
rz(-2.7244748) q[0];
sx q[0];
rz(1.3365082) q[0];
x q[1];
rz(-1.7393624) q[2];
sx q[2];
rz(-1.050808) q[2];
sx q[2];
rz(-1.2181768) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6042418) q[1];
sx q[1];
rz(-2.4927995) q[1];
sx q[1];
rz(-2.2082445) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8753136) q[3];
sx q[3];
rz(-2.7334573) q[3];
sx q[3];
rz(0.88312393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9271586) q[2];
sx q[2];
rz(-0.36403251) q[2];
sx q[2];
rz(2.7511609) q[2];
rz(-1.3123784) q[3];
sx q[3];
rz(-0.81109154) q[3];
sx q[3];
rz(-2.3217679) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44535962) q[0];
sx q[0];
rz(-0.92388988) q[0];
sx q[0];
rz(-1.4129289) q[0];
rz(-1.6784809) q[1];
sx q[1];
rz(-1.585958) q[1];
sx q[1];
rz(-0.63953343) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7752658) q[0];
sx q[0];
rz(-1.1970988) q[0];
sx q[0];
rz(0.92832066) q[0];
rz(2.9908871) q[2];
sx q[2];
rz(-1.4960714) q[2];
sx q[2];
rz(-3.0091803) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9078341) q[1];
sx q[1];
rz(-0.69926622) q[1];
sx q[1];
rz(-2.722867) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11951167) q[3];
sx q[3];
rz(-1.1788713) q[3];
sx q[3];
rz(-0.78531314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9245236) q[2];
sx q[2];
rz(-1.3288493) q[2];
sx q[2];
rz(0.45219839) q[2];
rz(2.9186115) q[3];
sx q[3];
rz(-0.28135869) q[3];
sx q[3];
rz(0.15534672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0105791) q[0];
sx q[0];
rz(-0.92229811) q[0];
sx q[0];
rz(-0.16192326) q[0];
rz(-0.34795347) q[1];
sx q[1];
rz(-1.8674928) q[1];
sx q[1];
rz(2.9071992) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9574675) q[0];
sx q[0];
rz(-2.1229366) q[0];
sx q[0];
rz(-2.2669621) q[0];
rz(-pi) q[1];
rz(-0.72004628) q[2];
sx q[2];
rz(-1.3338425) q[2];
sx q[2];
rz(-0.64575486) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.99410759) q[1];
sx q[1];
rz(-2.4723158) q[1];
sx q[1];
rz(1.2662751) q[1];
rz(-0.23080821) q[3];
sx q[3];
rz(-0.2838906) q[3];
sx q[3];
rz(-0.47286716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.28635412) q[2];
sx q[2];
rz(-2.0704634) q[2];
sx q[2];
rz(-0.62412778) q[2];
rz(-0.74328077) q[3];
sx q[3];
rz(-2.1478839) q[3];
sx q[3];
rz(0.68964094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2766992) q[0];
sx q[0];
rz(-1.445329) q[0];
sx q[0];
rz(-2.0546761) q[0];
rz(1.2326321) q[1];
sx q[1];
rz(-1.107736) q[1];
sx q[1];
rz(-1.3023652) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83643276) q[0];
sx q[0];
rz(-1.4576685) q[0];
sx q[0];
rz(-0.21048429) q[0];
x q[1];
rz(2.0696569) q[2];
sx q[2];
rz(-0.65827024) q[2];
sx q[2];
rz(1.9287623) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2375396) q[1];
sx q[1];
rz(-1.9642648) q[1];
sx q[1];
rz(2.5649125) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62303657) q[3];
sx q[3];
rz(-1.5554232) q[3];
sx q[3];
rz(-2.8492209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6752601) q[2];
sx q[2];
rz(-1.8856498) q[2];
sx q[2];
rz(-0.41435286) q[2];
rz(2.3710592) q[3];
sx q[3];
rz(-1.4221752) q[3];
sx q[3];
rz(-0.93366247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.82938021) q[0];
sx q[0];
rz(-0.82013622) q[0];
sx q[0];
rz(-2.6729551) q[0];
rz(1.8679484) q[1];
sx q[1];
rz(-1.2761152) q[1];
sx q[1];
rz(-0.31148568) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0429223) q[0];
sx q[0];
rz(-2.2880974) q[0];
sx q[0];
rz(0.076626549) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63705541) q[2];
sx q[2];
rz(-2.9961176) q[2];
sx q[2];
rz(1.5381952) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35986082) q[1];
sx q[1];
rz(-1.6174498) q[1];
sx q[1];
rz(2.0585039) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6139975) q[3];
sx q[3];
rz(-1.1403392) q[3];
sx q[3];
rz(-2.2077219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0353388) q[2];
sx q[2];
rz(-1.9474578) q[2];
sx q[2];
rz(-1.7485471) q[2];
rz(2.2140908) q[3];
sx q[3];
rz(-1.564097) q[3];
sx q[3];
rz(-2.2813796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2911693) q[0];
sx q[0];
rz(-2.8130154) q[0];
sx q[0];
rz(-1.6541506) q[0];
rz(0.18118478) q[1];
sx q[1];
rz(-0.94534992) q[1];
sx q[1];
rz(2.2015991) q[1];
rz(1.0370902) q[2];
sx q[2];
rz(-1.4363534) q[2];
sx q[2];
rz(-1.7912122) q[2];
rz(0.84862205) q[3];
sx q[3];
rz(-2.3325338) q[3];
sx q[3];
rz(2.4073413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

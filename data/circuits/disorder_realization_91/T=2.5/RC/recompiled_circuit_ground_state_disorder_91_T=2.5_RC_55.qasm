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
rz(4.3963764) q[0];
sx q[0];
rz(8.7788361) q[0];
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
rz(-2.7433919) q[0];
sx q[0];
rz(-1.8132134) q[0];
sx q[0];
rz(-1.6869808) q[0];
x q[1];
rz(-1.9930085) q[2];
sx q[2];
rz(-2.4496884) q[2];
sx q[2];
rz(0.9949323) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5234452) q[1];
sx q[1];
rz(-1.0511569) q[1];
sx q[1];
rz(-0.45975273) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6553538) q[3];
sx q[3];
rz(-2.231153) q[3];
sx q[3];
rz(2.5442991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6231923) q[2];
sx q[2];
rz(-2.4861768) q[2];
sx q[2];
rz(1.0830967) q[2];
rz(2.8862503) q[3];
sx q[3];
rz(-2.2446003) q[3];
sx q[3];
rz(-1.159509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.979368) q[0];
sx q[0];
rz(-3.0001682) q[0];
sx q[0];
rz(3.1080143) q[0];
rz(2.0032739) q[1];
sx q[1];
rz(-0.99457026) q[1];
sx q[1];
rz(0.23695645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4073674) q[0];
sx q[0];
rz(-1.4324643) q[0];
sx q[0];
rz(-2.9107735) q[0];
rz(-pi) q[1];
rz(-1.9230546) q[2];
sx q[2];
rz(-2.6389803) q[2];
sx q[2];
rz(1.1692804) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33920501) q[1];
sx q[1];
rz(-1.5453813) q[1];
sx q[1];
rz(-1.7632381) q[1];
rz(-pi) q[2];
rz(0.99144793) q[3];
sx q[3];
rz(-1.9526523) q[3];
sx q[3];
rz(2.6086311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.61243764) q[2];
sx q[2];
rz(-2.0945022) q[2];
sx q[2];
rz(3.0212413) q[2];
rz(0.58593166) q[3];
sx q[3];
rz(-1.7610995) q[3];
sx q[3];
rz(-1.8732635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(1.5536701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7415888) q[0];
sx q[0];
rz(-3.1256465) q[0];
sx q[0];
rz(0.069132968) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47207295) q[2];
sx q[2];
rz(-2.9909212) q[2];
sx q[2];
rz(-0.27418384) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7657679) q[1];
sx q[1];
rz(-1.3096403) q[1];
sx q[1];
rz(-3.1118666) q[1];
x q[2];
rz(2.5293143) q[3];
sx q[3];
rz(-0.66198549) q[3];
sx q[3];
rz(0.16890165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.83739088) q[2];
sx q[2];
rz(-1.7981671) q[2];
sx q[2];
rz(-0.070579441) q[2];
rz(-2.6523759) q[3];
sx q[3];
rz(-0.990812) q[3];
sx q[3];
rz(-0.14744559) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(-0.69545737) q[0];
sx q[0];
rz(1.7406933) q[0];
rz(2.9715624) q[1];
sx q[1];
rz(-1.731512) q[1];
sx q[1];
rz(-1.9084575) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2182541) q[0];
sx q[0];
rz(-1.4460576) q[0];
sx q[0];
rz(-1.9943004) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4958359) q[2];
sx q[2];
rz(-2.3120572) q[2];
sx q[2];
rz(-2.2956306) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.255862) q[1];
sx q[1];
rz(-1.8504603) q[1];
sx q[1];
rz(-3.1260296) q[1];
rz(3.1175851) q[3];
sx q[3];
rz(-0.46459651) q[3];
sx q[3];
rz(0.28991227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3418545) q[2];
sx q[2];
rz(-1.7065455) q[2];
sx q[2];
rz(0.62961659) q[2];
rz(0.13684212) q[3];
sx q[3];
rz(-2.5119731) q[3];
sx q[3];
rz(1.4153882) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28302309) q[0];
sx q[0];
rz(-2.6333599) q[0];
sx q[0];
rz(3.0392905) q[0];
rz(-1.6434068) q[1];
sx q[1];
rz(-1.841265) q[1];
sx q[1];
rz(-2.7844875) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.461991) q[0];
sx q[0];
rz(-1.3565738) q[0];
sx q[0];
rz(-2.9132183) q[0];
rz(-pi) q[1];
rz(-1.5979366) q[2];
sx q[2];
rz(-1.4476336) q[2];
sx q[2];
rz(2.9371098) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2374463) q[1];
sx q[1];
rz(-2.0761298) q[1];
sx q[1];
rz(1.2282441) q[1];
rz(-pi) q[2];
rz(-2.3042942) q[3];
sx q[3];
rz(-2.5243763) q[3];
sx q[3];
rz(-2.0715947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6685278) q[2];
sx q[2];
rz(-1.4184971) q[2];
sx q[2];
rz(1.3225887) q[2];
rz(-1.708301) q[3];
sx q[3];
rz(-1.8247484) q[3];
sx q[3];
rz(-0.1639666) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.598031) q[0];
sx q[0];
rz(-2.141641) q[0];
sx q[0];
rz(1.7033956) q[0];
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
rz(3.0594994) q[0];
sx q[0];
rz(-1.4766066) q[0];
sx q[0];
rz(1.1638428) q[0];
rz(1.7393624) q[2];
sx q[2];
rz(-2.0907846) q[2];
sx q[2];
rz(-1.2181768) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8554012) q[1];
sx q[1];
rz(-2.0778065) q[1];
sx q[1];
rz(-2.7176574) q[1];
rz(-0.26627906) q[3];
sx q[3];
rz(-2.7334573) q[3];
sx q[3];
rz(-0.88312393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.214434) q[2];
sx q[2];
rz(-0.36403251) q[2];
sx q[2];
rz(2.7511609) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.696233) q[0];
sx q[0];
rz(-2.2177028) q[0];
sx q[0];
rz(-1.7286638) q[0];
rz(-1.6784809) q[1];
sx q[1];
rz(-1.585958) q[1];
sx q[1];
rz(-0.63953343) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4832925) q[0];
sx q[0];
rz(-0.72972238) q[0];
sx q[0];
rz(2.150282) q[0];
rz(-2.6790304) q[2];
sx q[2];
rz(-2.9735045) q[2];
sx q[2];
rz(-1.8953022) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8486316) q[1];
sx q[1];
rz(-0.94215122) q[1];
sx q[1];
rz(-1.9002923) q[1];
x q[2];
rz(-0.11951167) q[3];
sx q[3];
rz(-1.9627213) q[3];
sx q[3];
rz(2.3562795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9245236) q[2];
sx q[2];
rz(-1.3288493) q[2];
sx q[2];
rz(0.45219839) q[2];
rz(-0.2229812) q[3];
sx q[3];
rz(-2.860234) q[3];
sx q[3];
rz(2.9862459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(3.0105791) q[0];
sx q[0];
rz(-2.2192945) q[0];
sx q[0];
rz(-0.16192326) q[0];
rz(0.34795347) q[1];
sx q[1];
rz(-1.8674928) q[1];
sx q[1];
rz(-2.9071992) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1944626) q[0];
sx q[0];
rz(-2.2827153) q[0];
sx q[0];
rz(-0.80545896) q[0];
rz(-pi) q[1];
rz(-0.35104819) q[2];
sx q[2];
rz(-0.75134885) q[2];
sx q[2];
rz(-2.4780432) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33501262) q[1];
sx q[1];
rz(-1.7579105) q[1];
sx q[1];
rz(2.2172865) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.864847) q[3];
sx q[3];
rz(-1.5066772) q[3];
sx q[3];
rz(-1.319805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.28635412) q[2];
sx q[2];
rz(-1.0711292) q[2];
sx q[2];
rz(-2.5174649) q[2];
rz(0.74328077) q[3];
sx q[3];
rz(-0.99370876) q[3];
sx q[3];
rz(-2.4519517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2766992) q[0];
sx q[0];
rz(-1.445329) q[0];
sx q[0];
rz(-1.0869166) q[0];
rz(1.9089606) q[1];
sx q[1];
rz(-2.0338567) q[1];
sx q[1];
rz(-1.3023652) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24827458) q[0];
sx q[0];
rz(-0.23856197) q[0];
sx q[0];
rz(-2.6435411) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0719358) q[2];
sx q[2];
rz(-2.4833224) q[2];
sx q[2];
rz(1.2128304) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.088868695) q[1];
sx q[1];
rz(-2.0985328) q[1];
sx q[1];
rz(-2.0305968) q[1];
rz(0.026342197) q[3];
sx q[3];
rz(-0.62320101) q[3];
sx q[3];
rz(-1.2570326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6752601) q[2];
sx q[2];
rz(-1.8856498) q[2];
sx q[2];
rz(2.7272398) q[2];
rz(-0.77053344) q[3];
sx q[3];
rz(-1.7194175) q[3];
sx q[3];
rz(-2.2079302) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82938021) q[0];
sx q[0];
rz(-2.3214564) q[0];
sx q[0];
rz(0.46863753) q[0];
rz(1.2736443) q[1];
sx q[1];
rz(-1.2761152) q[1];
sx q[1];
rz(-2.830107) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92665206) q[0];
sx q[0];
rz(-2.4209341) q[0];
sx q[0];
rz(-1.4832627) q[0];
x q[1];
rz(2.5045372) q[2];
sx q[2];
rz(-2.9961176) q[2];
sx q[2];
rz(1.6033974) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35986082) q[1];
sx q[1];
rz(-1.5241429) q[1];
sx q[1];
rz(1.0830888) q[1];
x q[2];
rz(-2.4021637) q[3];
sx q[3];
rz(-0.66777705) q[3];
sx q[3];
rz(-1.8831933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1062539) q[2];
sx q[2];
rz(-1.1941348) q[2];
sx q[2];
rz(-1.3930456) q[2];
rz(-2.2140908) q[3];
sx q[3];
rz(-1.564097) q[3];
sx q[3];
rz(2.2813796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2911693) q[0];
sx q[0];
rz(-2.8130154) q[0];
sx q[0];
rz(-1.6541506) q[0];
rz(-0.18118478) q[1];
sx q[1];
rz(-2.1962427) q[1];
sx q[1];
rz(-0.93999351) q[1];
rz(2.1045024) q[2];
sx q[2];
rz(-1.7052393) q[2];
sx q[2];
rz(1.3503804) q[2];
rz(-2.2929706) q[3];
sx q[3];
rz(-2.3325338) q[3];
sx q[3];
rz(2.4073413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

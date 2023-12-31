OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55943638) q[0];
sx q[0];
rz(3.732582) q[0];
sx q[0];
rz(8.8413722) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(-0.89259994) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8569782) q[0];
sx q[0];
rz(-2.1224788) q[0];
sx q[0];
rz(0.3509699) q[0];
rz(-0.66944389) q[2];
sx q[2];
rz(-1.9813683) q[2];
sx q[2];
rz(0.66561156) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8436369) q[1];
sx q[1];
rz(-1.2601818) q[1];
sx q[1];
rz(-0.85308869) q[1];
x q[2];
rz(1.6139908) q[3];
sx q[3];
rz(-0.82026635) q[3];
sx q[3];
rz(0.84212069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2261752) q[2];
sx q[2];
rz(-1.6025275) q[2];
sx q[2];
rz(-1.1606476) q[2];
rz(0.21696572) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(-1.0552361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0579257) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(2.5657186) q[0];
rz(1.2469762) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(-1.1670246) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407699) q[0];
sx q[0];
rz(-1.5123899) q[0];
sx q[0];
rz(1.3129243) q[0];
rz(1.0300693) q[2];
sx q[2];
rz(-1.3582555) q[2];
sx q[2];
rz(-2.9589047) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6845219) q[1];
sx q[1];
rz(-1.6324537) q[1];
sx q[1];
rz(-2.3350299) q[1];
x q[2];
rz(-2.4715273) q[3];
sx q[3];
rz(-1.9668005) q[3];
sx q[3];
rz(2.0585287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.033826753) q[2];
sx q[2];
rz(-1.1788538) q[2];
sx q[2];
rz(-2.9272184) q[2];
rz(0.073444627) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630994) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(1.9146772) q[0];
rz(0.40027174) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(-2.1267557) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9432362) q[0];
sx q[0];
rz(-2.495932) q[0];
sx q[0];
rz(-0.093430324) q[0];
x q[1];
rz(2.8995908) q[2];
sx q[2];
rz(-1.9654462) q[2];
sx q[2];
rz(-1.637527) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6940569) q[1];
sx q[1];
rz(-2.5390365) q[1];
sx q[1];
rz(1.7130997) q[1];
x q[2];
rz(0.26168163) q[3];
sx q[3];
rz(-0.90173429) q[3];
sx q[3];
rz(0.79479587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0959452) q[2];
sx q[2];
rz(-2.084705) q[2];
sx q[2];
rz(0.8992368) q[2];
rz(0.69747654) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(0.60788679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4863481) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(-2.7096601) q[0];
rz(2.5090384) q[1];
sx q[1];
rz(-0.41699854) q[1];
sx q[1];
rz(-0.63582173) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.302127) q[0];
sx q[0];
rz(-1.3913904) q[0];
sx q[0];
rz(-1.0666215) q[0];
rz(-pi) q[1];
rz(1.3859149) q[2];
sx q[2];
rz(-0.73542483) q[2];
sx q[2];
rz(-2.3787969) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1078474) q[1];
sx q[1];
rz(-2.3305232) q[1];
sx q[1];
rz(-1.8268405) q[1];
x q[2];
rz(-1.3577348) q[3];
sx q[3];
rz(-1.2902998) q[3];
sx q[3];
rz(-0.76663843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2137961) q[2];
sx q[2];
rz(-0.14363229) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(-0.33411807) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(2.9978602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14389811) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(1.0744263) q[0];
rz(-0.74514666) q[1];
sx q[1];
rz(-1.4829758) q[1];
sx q[1];
rz(-2.863046) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4771839) q[0];
sx q[0];
rz(-1.9946949) q[0];
sx q[0];
rz(2.3197078) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4611545) q[2];
sx q[2];
rz(-0.96193681) q[2];
sx q[2];
rz(1.7313752) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46036938) q[1];
sx q[1];
rz(-1.2037828) q[1];
sx q[1];
rz(-2.2296434) q[1];
rz(-pi) q[2];
rz(-0.55322247) q[3];
sx q[3];
rz(-2.2461938) q[3];
sx q[3];
rz(1.355987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6986065) q[2];
sx q[2];
rz(-1.3977945) q[2];
sx q[2];
rz(-2.8473575) q[2];
rz(3.0596628) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(-0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1086403) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(-0.0090573514) q[0];
rz(2.5065705) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(-0.10805282) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8979643) q[0];
sx q[0];
rz(-1.122323) q[0];
sx q[0];
rz(-1.2178221) q[0];
rz(-pi) q[1];
rz(0.33424218) q[2];
sx q[2];
rz(-2.7396482) q[2];
sx q[2];
rz(-2.3964756) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6335771) q[1];
sx q[1];
rz(-1.0883691) q[1];
sx q[1];
rz(1.1303933) q[1];
rz(-pi) q[2];
rz(-1.024527) q[3];
sx q[3];
rz(-1.7987935) q[3];
sx q[3];
rz(2.8834164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.99284995) q[2];
sx q[2];
rz(-1.7603346) q[2];
sx q[2];
rz(-1.2711058) q[2];
rz(-3.0631915) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(1.1318077) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25061297) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(-2.4839731) q[0];
rz(1.3972067) q[1];
sx q[1];
rz(-1.1746635) q[1];
sx q[1];
rz(2.2479642) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0890869) q[0];
sx q[0];
rz(-0.94852175) q[0];
sx q[0];
rz(-1.9786579) q[0];
rz(0.29427476) q[2];
sx q[2];
rz(-1.1106967) q[2];
sx q[2];
rz(-1.0690881) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5205198) q[1];
sx q[1];
rz(-2.2398661) q[1];
sx q[1];
rz(1.3135765) q[1];
x q[2];
rz(-0.75974792) q[3];
sx q[3];
rz(-1.500251) q[3];
sx q[3];
rz(-2.7618046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7509193) q[2];
sx q[2];
rz(-2.1942287) q[2];
sx q[2];
rz(-2.612109) q[2];
rz(0.47618619) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(2.7990394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7628409) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(-2.0513127) q[0];
rz(3.0293363) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(1.9876678) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33015108) q[0];
sx q[0];
rz(-2.220287) q[0];
sx q[0];
rz(1.5280456) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23837337) q[2];
sx q[2];
rz(-2.5775238) q[2];
sx q[2];
rz(-2.7114045) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6653319) q[1];
sx q[1];
rz(-1.6325103) q[1];
sx q[1];
rz(0.76121059) q[1];
x q[2];
rz(1.0941986) q[3];
sx q[3];
rz(-2.4123441) q[3];
sx q[3];
rz(0.13343982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.551679) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(-2.7611458) q[2];
rz(-2.0137265) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(-1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8001051) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(0.63968101) q[0];
rz(-1.2387964) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(1.9715086) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026924883) q[0];
sx q[0];
rz(-1.9244734) q[0];
sx q[0];
rz(0.35004079) q[0];
rz(-pi) q[1];
rz(1.8413999) q[2];
sx q[2];
rz(-0.81297183) q[2];
sx q[2];
rz(1.1635309) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.77800345) q[1];
sx q[1];
rz(-0.52007857) q[1];
sx q[1];
rz(-0.37104721) q[1];
rz(-pi) q[2];
rz(-0.13799237) q[3];
sx q[3];
rz(-0.90413168) q[3];
sx q[3];
rz(-0.43825144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8273948) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(-2.7588552) q[2];
rz(-0.9283723) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(0.66463566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2255573) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(-2.8826707) q[0];
rz(0.71031538) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(0.47992596) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5537162) q[0];
sx q[0];
rz(-1.9037316) q[0];
sx q[0];
rz(-2.2399708) q[0];
rz(-pi) q[1];
rz(-2.2357113) q[2];
sx q[2];
rz(-1.7582298) q[2];
sx q[2];
rz(-2.4311709) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.48206115) q[1];
sx q[1];
rz(-1.2188984) q[1];
sx q[1];
rz(2.5955276) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89703538) q[3];
sx q[3];
rz(-1.1437136) q[3];
sx q[3];
rz(0.87891146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2991128) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(1.8245565) q[2];
rz(-1.2420098) q[3];
sx q[3];
rz(-2.1879991) q[3];
sx q[3];
rz(-0.73808134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15923545) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(-2.1622529) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(2.813415) q[2];
sx q[2];
rz(-2.6371418) q[2];
sx q[2];
rz(-0.20124659) q[2];
rz(-2.8638774) q[3];
sx q[3];
rz(-1.82901) q[3];
sx q[3];
rz(1.3808586) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

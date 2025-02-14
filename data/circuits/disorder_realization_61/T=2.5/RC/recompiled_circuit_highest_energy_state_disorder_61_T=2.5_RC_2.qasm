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
rz(0.50245589) q[0];
sx q[0];
rz(4.2491663) q[0];
sx q[0];
rz(7.6585461) q[0];
rz(-0.31887588) q[1];
sx q[1];
rz(-1.640929) q[1];
sx q[1];
rz(-2.4102416) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1494182) q[0];
sx q[0];
rz(-0.96620027) q[0];
sx q[0];
rz(0.3913619) q[0];
x q[1];
rz(-1.574975) q[2];
sx q[2];
rz(-2.2327645) q[2];
sx q[2];
rz(2.3062492) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13818422) q[1];
sx q[1];
rz(-2.3560682) q[1];
sx q[1];
rz(-1.5064606) q[1];
rz(-pi) q[2];
rz(-0.5238138) q[3];
sx q[3];
rz(-0.18693811) q[3];
sx q[3];
rz(-2.4863249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2607164) q[2];
sx q[2];
rz(-1.2141576) q[2];
sx q[2];
rz(-2.1991849) q[2];
rz(0.85855329) q[3];
sx q[3];
rz(-2.5311354) q[3];
sx q[3];
rz(2.4220991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.079916) q[0];
sx q[0];
rz(-2.583857) q[0];
sx q[0];
rz(3.0829561) q[0];
rz(-3.0851641) q[1];
sx q[1];
rz(-1.7516878) q[1];
sx q[1];
rz(-1.1172392) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65172136) q[0];
sx q[0];
rz(-0.84017032) q[0];
sx q[0];
rz(-2.1108759) q[0];
x q[1];
rz(-1.8705759) q[2];
sx q[2];
rz(-1.6938871) q[2];
sx q[2];
rz(-1.276615) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2721418) q[1];
sx q[1];
rz(-1.488416) q[1];
sx q[1];
rz(-1.8100965) q[1];
rz(-0.47567077) q[3];
sx q[3];
rz(-1.6783617) q[3];
sx q[3];
rz(1.8799262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3650018) q[2];
sx q[2];
rz(-2.8114909) q[2];
sx q[2];
rz(-2.7351725) q[2];
rz(2.8255919) q[3];
sx q[3];
rz(-1.4603115) q[3];
sx q[3];
rz(-2.3587904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86421788) q[0];
sx q[0];
rz(-2.6835231) q[0];
sx q[0];
rz(-1.3360485) q[0];
rz(2.863073) q[1];
sx q[1];
rz(-1.3986992) q[1];
sx q[1];
rz(-2.1727402) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012798158) q[0];
sx q[0];
rz(-1.9417106) q[0];
sx q[0];
rz(2.4910035) q[0];
rz(0.79040852) q[2];
sx q[2];
rz(-1.9876754) q[2];
sx q[2];
rz(-1.1000772) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.90889764) q[1];
sx q[1];
rz(-1.2746342) q[1];
sx q[1];
rz(3.0975268) q[1];
rz(-2.8519451) q[3];
sx q[3];
rz(-2.395219) q[3];
sx q[3];
rz(1.0371072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99732533) q[2];
sx q[2];
rz(-0.48339016) q[2];
sx q[2];
rz(-2.9664795) q[2];
rz(-1.8363606) q[3];
sx q[3];
rz(-2.1988018) q[3];
sx q[3];
rz(-1.85873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1433732) q[0];
sx q[0];
rz(-1.136919) q[0];
sx q[0];
rz(-2.112222) q[0];
rz(2.3221305) q[1];
sx q[1];
rz(-1.1623323) q[1];
sx q[1];
rz(-0.79684657) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40596698) q[0];
sx q[0];
rz(-1.7628364) q[0];
sx q[0];
rz(-2.6940932) q[0];
x q[1];
rz(2.098747) q[2];
sx q[2];
rz(-1.5656131) q[2];
sx q[2];
rz(0.56765616) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0681035) q[1];
sx q[1];
rz(-1.7760381) q[1];
sx q[1];
rz(-1.2935335) q[1];
x q[2];
rz(0.53189532) q[3];
sx q[3];
rz(-0.73557702) q[3];
sx q[3];
rz(0.07594219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.13021079) q[2];
sx q[2];
rz(-3.104976) q[2];
sx q[2];
rz(-1.2288564) q[2];
rz(2.7025488) q[3];
sx q[3];
rz(-1.565758) q[3];
sx q[3];
rz(1.9108093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.48290408) q[0];
sx q[0];
rz(-1.6084325) q[0];
sx q[0];
rz(-0.80832344) q[0];
rz(-1.3574379) q[1];
sx q[1];
rz(-0.789855) q[1];
sx q[1];
rz(-2.435991) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43967908) q[0];
sx q[0];
rz(-1.8712303) q[0];
sx q[0];
rz(-1.6070532) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86807495) q[2];
sx q[2];
rz(-1.3931615) q[2];
sx q[2];
rz(2.9736819) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29632211) q[1];
sx q[1];
rz(-2.8883838) q[1];
sx q[1];
rz(1.5814387) q[1];
rz(-0.79885317) q[3];
sx q[3];
rz(-2.0435413) q[3];
sx q[3];
rz(-0.96351876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6610873) q[2];
sx q[2];
rz(-1.8578119) q[2];
sx q[2];
rz(2.535848) q[2];
rz(3.0211966) q[3];
sx q[3];
rz(-1.8431289) q[3];
sx q[3];
rz(-2.3275183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26009387) q[0];
sx q[0];
rz(-2.3277178) q[0];
sx q[0];
rz(3.1392198) q[0];
rz(-2.782605) q[1];
sx q[1];
rz(-1.2009883) q[1];
sx q[1];
rz(-2.733309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6946163) q[0];
sx q[0];
rz(-0.13152619) q[0];
sx q[0];
rz(1.294554) q[0];
x q[1];
rz(-1.194656) q[2];
sx q[2];
rz(-0.6413817) q[2];
sx q[2];
rz(0.94571404) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1241345) q[1];
sx q[1];
rz(-1.5038905) q[1];
sx q[1];
rz(0.23048877) q[1];
rz(-pi) q[2];
rz(-2.1903135) q[3];
sx q[3];
rz(-0.86182098) q[3];
sx q[3];
rz(-0.10892841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4957054) q[2];
sx q[2];
rz(-2.3453823) q[2];
sx q[2];
rz(0.1296981) q[2];
rz(0.4194704) q[3];
sx q[3];
rz(-1.9286112) q[3];
sx q[3];
rz(-0.81370846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7523338) q[0];
sx q[0];
rz(-0.77521721) q[0];
sx q[0];
rz(2.4580521) q[0];
rz(0.93841249) q[1];
sx q[1];
rz(-1.7529731) q[1];
sx q[1];
rz(-1.9155115) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3049604) q[0];
sx q[0];
rz(-1.9373405) q[0];
sx q[0];
rz(2.5185597) q[0];
rz(-pi) q[1];
rz(-0.60003321) q[2];
sx q[2];
rz(-1.4156315) q[2];
sx q[2];
rz(-0.70916353) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9270806) q[1];
sx q[1];
rz(-1.1225554) q[1];
sx q[1];
rz(-0.31772504) q[1];
rz(-pi) q[2];
rz(-0.88507401) q[3];
sx q[3];
rz(-0.7739017) q[3];
sx q[3];
rz(0.89088551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.73416) q[2];
sx q[2];
rz(-1.7075044) q[2];
sx q[2];
rz(-0.70719353) q[2];
rz(1.6678984) q[3];
sx q[3];
rz(-2.4652822) q[3];
sx q[3];
rz(-1.7128806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4555175) q[0];
sx q[0];
rz(-0.34611836) q[0];
sx q[0];
rz(2.2177875) q[0];
rz(1.0651945) q[1];
sx q[1];
rz(-1.1895475) q[1];
sx q[1];
rz(-2.0069897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4489733) q[0];
sx q[0];
rz(-1.53063) q[0];
sx q[0];
rz(1.9461826) q[0];
x q[1];
rz(-0.8834817) q[2];
sx q[2];
rz(-1.0687912) q[2];
sx q[2];
rz(-2.4740681) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2345786) q[1];
sx q[1];
rz(-1.5173787) q[1];
sx q[1];
rz(2.1018886) q[1];
rz(-pi) q[2];
rz(1.9530961) q[3];
sx q[3];
rz(-0.99671066) q[3];
sx q[3];
rz(2.7758383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.0040466641) q[2];
sx q[2];
rz(-1.7377487) q[2];
sx q[2];
rz(-0.64296067) q[2];
rz(0.17063394) q[3];
sx q[3];
rz(-2.4891487) q[3];
sx q[3];
rz(1.0954866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74649015) q[0];
sx q[0];
rz(-3.009142) q[0];
sx q[0];
rz(0.80628959) q[0];
rz(-0.0099446615) q[1];
sx q[1];
rz(-2.5237623) q[1];
sx q[1];
rz(-2.8795805) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0341971) q[0];
sx q[0];
rz(-0.55343628) q[0];
sx q[0];
rz(-1.2518) q[0];
rz(-pi) q[1];
rz(-2.9295983) q[2];
sx q[2];
rz(-0.65644453) q[2];
sx q[2];
rz(-2.580263) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6422185) q[1];
sx q[1];
rz(-1.1406349) q[1];
sx q[1];
rz(-2.9146934) q[1];
rz(-1.069182) q[3];
sx q[3];
rz(-1.8610238) q[3];
sx q[3];
rz(-0.25267664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3917196) q[2];
sx q[2];
rz(-2.5842857) q[2];
sx q[2];
rz(0.14287512) q[2];
rz(-2.2376132) q[3];
sx q[3];
rz(-1.4986135) q[3];
sx q[3];
rz(-2.2601295) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3751462) q[0];
sx q[0];
rz(-1.2297933) q[0];
sx q[0];
rz(-2.4850856) q[0];
rz(-0.87908602) q[1];
sx q[1];
rz(-1.9706985) q[1];
sx q[1];
rz(-0.62969977) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64163369) q[0];
sx q[0];
rz(-2.2370506) q[0];
sx q[0];
rz(-2.0883661) q[0];
x q[1];
rz(-1.1348502) q[2];
sx q[2];
rz(-2.5647047) q[2];
sx q[2];
rz(-2.1132121) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.39552626) q[1];
sx q[1];
rz(-1.6159355) q[1];
sx q[1];
rz(-3.1137115) q[1];
rz(-1.0694169) q[3];
sx q[3];
rz(-0.58749858) q[3];
sx q[3];
rz(-2.0794433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.62275824) q[2];
sx q[2];
rz(-0.82509416) q[2];
sx q[2];
rz(2.4614914) q[2];
rz(2.425219) q[3];
sx q[3];
rz(-3.0383737) q[3];
sx q[3];
rz(1.1187925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685365) q[0];
sx q[0];
rz(-1.2328883) q[0];
sx q[0];
rz(0.20986025) q[0];
rz(-0.14280351) q[1];
sx q[1];
rz(-2.0695984) q[1];
sx q[1];
rz(-2.4048068) q[1];
rz(0.057351107) q[2];
sx q[2];
rz(-1.8125368) q[2];
sx q[2];
rz(-2.9619458) q[2];
rz(-2.18022) q[3];
sx q[3];
rz(-1.537804) q[3];
sx q[3];
rz(0.15450879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

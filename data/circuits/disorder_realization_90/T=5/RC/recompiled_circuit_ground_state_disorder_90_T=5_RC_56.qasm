OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2990155) q[0];
sx q[0];
rz(-2.9691073) q[0];
sx q[0];
rz(-0.01699288) q[0];
rz(0.51141557) q[1];
sx q[1];
rz(-0.49632448) q[1];
sx q[1];
rz(0.48547784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6926138) q[0];
sx q[0];
rz(-2.9831121) q[0];
sx q[0];
rz(0.7787718) q[0];
x q[1];
rz(0.83037776) q[2];
sx q[2];
rz(-0.86566209) q[2];
sx q[2];
rz(-0.0041088897) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37066165) q[1];
sx q[1];
rz(-2.0699796) q[1];
sx q[1];
rz(-0.73727495) q[1];
x q[2];
rz(1.2258244) q[3];
sx q[3];
rz(-0.32805035) q[3];
sx q[3];
rz(-1.5313182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39712507) q[2];
sx q[2];
rz(-3.0391389) q[2];
sx q[2];
rz(2.3490119) q[2];
rz(0.98627311) q[3];
sx q[3];
rz(-1.4873742) q[3];
sx q[3];
rz(3.0084394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2546286) q[0];
sx q[0];
rz(-1.5381085) q[0];
sx q[0];
rz(-0.1057374) q[0];
rz(-0.75791439) q[1];
sx q[1];
rz(-0.71280232) q[1];
sx q[1];
rz(2.6656718) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71005586) q[0];
sx q[0];
rz(-1.8049631) q[0];
sx q[0];
rz(1.6706927) q[0];
rz(1.7143102) q[2];
sx q[2];
rz(-1.0821618) q[2];
sx q[2];
rz(1.2554864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9502752) q[1];
sx q[1];
rz(-1.4498267) q[1];
sx q[1];
rz(1.0005534) q[1];
x q[2];
rz(1.9571689) q[3];
sx q[3];
rz(-0.68102057) q[3];
sx q[3];
rz(1.0134361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9953352) q[2];
sx q[2];
rz(-2.6671851) q[2];
sx q[2];
rz(1.8246626) q[2];
rz(-2.3138192) q[3];
sx q[3];
rz(-2.0483978) q[3];
sx q[3];
rz(0.83724666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4227609) q[0];
sx q[0];
rz(-1.3439002) q[0];
sx q[0];
rz(1.9047009) q[0];
rz(1.749136) q[1];
sx q[1];
rz(-1.720865) q[1];
sx q[1];
rz(-0.69033355) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25094098) q[0];
sx q[0];
rz(-1.1162045) q[0];
sx q[0];
rz(0.62418749) q[0];
x q[1];
rz(-0.0039991556) q[2];
sx q[2];
rz(-2.0341113) q[2];
sx q[2];
rz(-0.057553854) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.54385932) q[1];
sx q[1];
rz(-1.3784474) q[1];
sx q[1];
rz(-1.400977) q[1];
rz(-pi) q[2];
rz(-0.9564444) q[3];
sx q[3];
rz(-2.2057475) q[3];
sx q[3];
rz(0.81723833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51848015) q[2];
sx q[2];
rz(-1.5223576) q[2];
sx q[2];
rz(3.0677262) q[2];
rz(1.9833924) q[3];
sx q[3];
rz(-1.9246212) q[3];
sx q[3];
rz(1.7346252) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2288007) q[0];
sx q[0];
rz(-3.085629) q[0];
sx q[0];
rz(2.9486616) q[0];
rz(1.8979161) q[1];
sx q[1];
rz(-1.396215) q[1];
sx q[1];
rz(-0.23304932) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0351999) q[0];
sx q[0];
rz(-1.4241204) q[0];
sx q[0];
rz(1.4752409) q[0];
rz(-pi) q[1];
rz(-1.6150171) q[2];
sx q[2];
rz(-2.4323175) q[2];
sx q[2];
rz(1.0929293) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0600507) q[1];
sx q[1];
rz(-0.31585187) q[1];
sx q[1];
rz(-2.9595023) q[1];
rz(2.1904864) q[3];
sx q[3];
rz(-2.4220805) q[3];
sx q[3];
rz(0.48966416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9243246) q[2];
sx q[2];
rz(-1.3793719) q[2];
sx q[2];
rz(3.0775089) q[2];
rz(-0.25121769) q[3];
sx q[3];
rz(-2.678674) q[3];
sx q[3];
rz(-0.29138756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083704405) q[0];
sx q[0];
rz(-1.8479481) q[0];
sx q[0];
rz(1.8573014) q[0];
rz(-0.0074370782) q[1];
sx q[1];
rz(-1.9556655) q[1];
sx q[1];
rz(-2.1844905) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5610661) q[0];
sx q[0];
rz(-1.2532803) q[0];
sx q[0];
rz(2.0317715) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61422698) q[2];
sx q[2];
rz(-2.0822656) q[2];
sx q[2];
rz(1.5661256) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9451198) q[1];
sx q[1];
rz(-0.34356657) q[1];
sx q[1];
rz(-2.8445704) q[1];
rz(-0.65517218) q[3];
sx q[3];
rz(-1.7196894) q[3];
sx q[3];
rz(-2.6096491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1288422) q[2];
sx q[2];
rz(-2.4350171) q[2];
sx q[2];
rz(1.0303222) q[2];
rz(-0.71850592) q[3];
sx q[3];
rz(-2.0593144) q[3];
sx q[3];
rz(-2.4350186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4265823) q[0];
sx q[0];
rz(-2.3937245) q[0];
sx q[0];
rz(-0.36710516) q[0];
rz(-0.35588613) q[1];
sx q[1];
rz(-2.0242736) q[1];
sx q[1];
rz(1.5497367) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7283495) q[0];
sx q[0];
rz(-1.6093045) q[0];
sx q[0];
rz(2.4000077) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0075304) q[2];
sx q[2];
rz(-1.6730089) q[2];
sx q[2];
rz(1.619316) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.983596) q[1];
sx q[1];
rz(-1.776773) q[1];
sx q[1];
rz(2.7218444) q[1];
x q[2];
rz(-2.995302) q[3];
sx q[3];
rz(-2.6669569) q[3];
sx q[3];
rz(-0.65505799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1648569) q[2];
sx q[2];
rz(-1.848449) q[2];
sx q[2];
rz(0.0034275835) q[2];
rz(-2.2554743) q[3];
sx q[3];
rz(-2.0485853) q[3];
sx q[3];
rz(-1.4440943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.7057328) q[0];
sx q[0];
rz(-1.4737361) q[0];
sx q[0];
rz(-2.4263897) q[0];
rz(0.12229478) q[1];
sx q[1];
rz(-1.0194174) q[1];
sx q[1];
rz(2.9939647) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0334761) q[0];
sx q[0];
rz(-1.8980935) q[0];
sx q[0];
rz(1.9294338) q[0];
x q[1];
rz(-1.9609902) q[2];
sx q[2];
rz(-1.1653333) q[2];
sx q[2];
rz(-2.6801339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6138654) q[1];
sx q[1];
rz(-2.4151049) q[1];
sx q[1];
rz(-1.8469091) q[1];
rz(-1.857599) q[3];
sx q[3];
rz(-1.5032603) q[3];
sx q[3];
rz(0.90112858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3132396) q[2];
sx q[2];
rz(-0.30892631) q[2];
sx q[2];
rz(-0.1114791) q[2];
rz(0.76505032) q[3];
sx q[3];
rz(-1.4234411) q[3];
sx q[3];
rz(-3.1077207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2400804) q[0];
sx q[0];
rz(-2.1724048) q[0];
sx q[0];
rz(2.1387658) q[0];
rz(0.78549939) q[1];
sx q[1];
rz(-1.5248884) q[1];
sx q[1];
rz(-1.2202107) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5883137) q[0];
sx q[0];
rz(-1.0568585) q[0];
sx q[0];
rz(1.5945503) q[0];
rz(-pi) q[1];
rz(3.0790824) q[2];
sx q[2];
rz(-1.4803807) q[2];
sx q[2];
rz(-0.16363283) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.073382811) q[1];
sx q[1];
rz(-1.0676976) q[1];
sx q[1];
rz(2.7791609) q[1];
rz(-pi) q[2];
rz(1.7986913) q[3];
sx q[3];
rz(-0.64809767) q[3];
sx q[3];
rz(2.9886887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9057374) q[2];
sx q[2];
rz(-1.7451655) q[2];
sx q[2];
rz(3.1077969) q[2];
rz(-2.3679521) q[3];
sx q[3];
rz(-1.6126817) q[3];
sx q[3];
rz(-1.0627559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3646024) q[0];
sx q[0];
rz(-2.5211625) q[0];
sx q[0];
rz(-0.34238368) q[0];
rz(1.4211593) q[1];
sx q[1];
rz(-2.3183289) q[1];
sx q[1];
rz(1.8470496) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0598885) q[0];
sx q[0];
rz(-1.9317288) q[0];
sx q[0];
rz(1.4128708) q[0];
rz(-pi) q[1];
x q[1];
rz(2.420408) q[2];
sx q[2];
rz(-0.9936665) q[2];
sx q[2];
rz(0.78109988) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8322014) q[1];
sx q[1];
rz(-2.3497255) q[1];
sx q[1];
rz(1.743508) q[1];
rz(-2.6813981) q[3];
sx q[3];
rz(-0.37307533) q[3];
sx q[3];
rz(2.4378547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2932059) q[2];
sx q[2];
rz(-1.3853955) q[2];
sx q[2];
rz(0.53655857) q[2];
rz(2.7287591) q[3];
sx q[3];
rz(-0.53240132) q[3];
sx q[3];
rz(-1.1134953) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8263016) q[0];
sx q[0];
rz(-1.8719801) q[0];
sx q[0];
rz(1.6444561) q[0];
rz(-2.3313088) q[1];
sx q[1];
rz(-1.5565926) q[1];
sx q[1];
rz(-2.2946045) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4627914) q[0];
sx q[0];
rz(-1.7399995) q[0];
sx q[0];
rz(-2.9786318) q[0];
rz(-pi) q[1];
rz(-2.2184847) q[2];
sx q[2];
rz(-1.1121684) q[2];
sx q[2];
rz(3.0748526) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3851953) q[1];
sx q[1];
rz(-2.1902731) q[1];
sx q[1];
rz(-1.9036915) q[1];
rz(0.096135898) q[3];
sx q[3];
rz(-1.8435974) q[3];
sx q[3];
rz(-0.48855272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.102313) q[2];
sx q[2];
rz(-1.8259093) q[2];
sx q[2];
rz(-0.63423356) q[2];
rz(2.5041194) q[3];
sx q[3];
rz(-2.0135148) q[3];
sx q[3];
rz(-2.9243961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98916003) q[0];
sx q[0];
rz(-1.9518873) q[0];
sx q[0];
rz(1.963203) q[0];
rz(2.4329026) q[1];
sx q[1];
rz(-1.7955753) q[1];
sx q[1];
rz(1.2830455) q[1];
rz(1.218956) q[2];
sx q[2];
rz(-1.2165804) q[2];
sx q[2];
rz(2.658398) q[2];
rz(-1.1905963) q[3];
sx q[3];
rz(-2.1722542) q[3];
sx q[3];
rz(0.62361591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

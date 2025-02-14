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
rz(0.99474466) q[0];
sx q[0];
rz(-2.5543307) q[0];
sx q[0];
rz(0.62556148) q[0];
rz(-2.5211531) q[1];
sx q[1];
rz(-0.99045366) q[1];
sx q[1];
rz(2.3576696) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44255689) q[0];
sx q[0];
rz(-0.43414206) q[0];
sx q[0];
rz(3.1331557) q[0];
x q[1];
rz(2.6033635) q[2];
sx q[2];
rz(-2.9799649) q[2];
sx q[2];
rz(1.5618351) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82134889) q[1];
sx q[1];
rz(-0.87961266) q[1];
sx q[1];
rz(2.3752179) q[1];
rz(-pi) q[2];
rz(-1.4830071) q[3];
sx q[3];
rz(-0.79500073) q[3];
sx q[3];
rz(0.82124293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1938532) q[2];
sx q[2];
rz(-1.0449907) q[2];
sx q[2];
rz(0.6673153) q[2];
rz(-0.31467485) q[3];
sx q[3];
rz(-2.5421725) q[3];
sx q[3];
rz(-2.2144894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84745234) q[0];
sx q[0];
rz(-0.84646928) q[0];
sx q[0];
rz(2.8417929) q[0];
rz(-1.0046129) q[1];
sx q[1];
rz(-2.424898) q[1];
sx q[1];
rz(2.2915548) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3687821) q[0];
sx q[0];
rz(-1.4748992) q[0];
sx q[0];
rz(-0.96536388) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9710567) q[2];
sx q[2];
rz(-1.6493622) q[2];
sx q[2];
rz(-2.6002392) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0176702) q[1];
sx q[1];
rz(-2.7233989) q[1];
sx q[1];
rz(-1.0355509) q[1];
rz(-2.7369954) q[3];
sx q[3];
rz(-1.5834152) q[3];
sx q[3];
rz(0.22443988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5726418) q[2];
sx q[2];
rz(-0.53223842) q[2];
sx q[2];
rz(0.66298318) q[2];
rz(-2.3891383) q[3];
sx q[3];
rz(-1.821725) q[3];
sx q[3];
rz(-2.9451356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74333423) q[0];
sx q[0];
rz(-0.4011811) q[0];
sx q[0];
rz(0.67984003) q[0];
rz(2.4196449) q[1];
sx q[1];
rz(-2.5847021) q[1];
sx q[1];
rz(-2.360875) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3762622) q[0];
sx q[0];
rz(-1.9714545) q[0];
sx q[0];
rz(-0.0082992502) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6569312) q[2];
sx q[2];
rz(-1.5670098) q[2];
sx q[2];
rz(0.45391339) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5131849) q[1];
sx q[1];
rz(-0.76090136) q[1];
sx q[1];
rz(-0.55879467) q[1];
rz(-0.2188245) q[3];
sx q[3];
rz(-2.2940562) q[3];
sx q[3];
rz(1.3506324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9146933) q[2];
sx q[2];
rz(-1.9136027) q[2];
sx q[2];
rz(-2.4719888) q[2];
rz(-2.2412444) q[3];
sx q[3];
rz(-1.0122296) q[3];
sx q[3];
rz(-2.3646234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25903073) q[0];
sx q[0];
rz(-0.43864033) q[0];
sx q[0];
rz(-0.90748179) q[0];
rz(0.59440815) q[1];
sx q[1];
rz(-0.92955697) q[1];
sx q[1];
rz(-1.1297191) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27760273) q[0];
sx q[0];
rz(-2.4765795) q[0];
sx q[0];
rz(-1.2590842) q[0];
rz(-2.1024804) q[2];
sx q[2];
rz(-1.6999107) q[2];
sx q[2];
rz(1.521284) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0320624) q[1];
sx q[1];
rz(-2.1550165) q[1];
sx q[1];
rz(1.24448) q[1];
rz(-1.2775189) q[3];
sx q[3];
rz(-1.8403111) q[3];
sx q[3];
rz(0.5130918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.10924673) q[2];
sx q[2];
rz(-1.3999908) q[2];
sx q[2];
rz(0.67231154) q[2];
rz(-3.1224871) q[3];
sx q[3];
rz(-0.34135434) q[3];
sx q[3];
rz(0.13489558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8371209) q[0];
sx q[0];
rz(-0.79455513) q[0];
sx q[0];
rz(-0.16803148) q[0];
rz(1.9145603) q[1];
sx q[1];
rz(-1.9214168) q[1];
sx q[1];
rz(2.8438445) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7452324) q[0];
sx q[0];
rz(-2.0848635) q[0];
sx q[0];
rz(2.5395509) q[0];
rz(-2.5644659) q[2];
sx q[2];
rz(-1.071605) q[2];
sx q[2];
rz(-1.2291723) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7682957) q[1];
sx q[1];
rz(-1.776564) q[1];
sx q[1];
rz(-1.5453669) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5750507) q[3];
sx q[3];
rz(-2.9828072) q[3];
sx q[3];
rz(2.6095853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1905404) q[2];
sx q[2];
rz(-1.0995883) q[2];
sx q[2];
rz(2.5229048) q[2];
rz(2.3330073) q[3];
sx q[3];
rz(-0.4777258) q[3];
sx q[3];
rz(1.8019069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.488778) q[0];
sx q[0];
rz(-1.5050911) q[0];
sx q[0];
rz(-2.6763647) q[0];
rz(-0.48509994) q[1];
sx q[1];
rz(-2.7310889) q[1];
sx q[1];
rz(3.0533275) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6104077) q[0];
sx q[0];
rz(-1.0779128) q[0];
sx q[0];
rz(2.6238074) q[0];
rz(-0.012862269) q[2];
sx q[2];
rz(-0.84190997) q[2];
sx q[2];
rz(0.60573214) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3151206) q[1];
sx q[1];
rz(-1.7355647) q[1];
sx q[1];
rz(0.75854782) q[1];
rz(0.28782423) q[3];
sx q[3];
rz(-1.6555641) q[3];
sx q[3];
rz(-0.945795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7469067) q[2];
sx q[2];
rz(-0.8553108) q[2];
sx q[2];
rz(2.6439903) q[2];
rz(-0.48007128) q[3];
sx q[3];
rz(-0.92104715) q[3];
sx q[3];
rz(3.0729496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8467872) q[0];
sx q[0];
rz(-2.6755896) q[0];
sx q[0];
rz(-1.71126) q[0];
rz(1.826674) q[1];
sx q[1];
rz(-2.7480835) q[1];
sx q[1];
rz(1.5865145) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2908132) q[0];
sx q[0];
rz(-1.6899791) q[0];
sx q[0];
rz(1.1759773) q[0];
x q[1];
rz(0.85679634) q[2];
sx q[2];
rz(-1.1044377) q[2];
sx q[2];
rz(3.015632) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79826971) q[1];
sx q[1];
rz(-2.1499691) q[1];
sx q[1];
rz(-0.12076853) q[1];
rz(-pi) q[2];
rz(2.2022543) q[3];
sx q[3];
rz(-1.5087481) q[3];
sx q[3];
rz(-2.5417059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.84307182) q[2];
sx q[2];
rz(-2.4303747) q[2];
sx q[2];
rz(2.2141875) q[2];
rz(1.1585506) q[3];
sx q[3];
rz(-0.68803334) q[3];
sx q[3];
rz(3.0454175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4226828) q[0];
sx q[0];
rz(-2.2490608) q[0];
sx q[0];
rz(2.6655777) q[0];
rz(-0.99126518) q[1];
sx q[1];
rz(-0.74949336) q[1];
sx q[1];
rz(3.057726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059444025) q[0];
sx q[0];
rz(-2.5184439) q[0];
sx q[0];
rz(-2.2150458) q[0];
x q[1];
rz(-2.4886683) q[2];
sx q[2];
rz(-2.5819375) q[2];
sx q[2];
rz(1.4687572) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0202581) q[1];
sx q[1];
rz(-0.42829207) q[1];
sx q[1];
rz(-1.7767724) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9355385) q[3];
sx q[3];
rz(-1.4901154) q[3];
sx q[3];
rz(1.9239359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91741651) q[2];
sx q[2];
rz(-0.26599628) q[2];
sx q[2];
rz(2.4931397) q[2];
rz(-2.2367541) q[3];
sx q[3];
rz(-1.6831393) q[3];
sx q[3];
rz(2.8308433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7998578) q[0];
sx q[0];
rz(-0.71001995) q[0];
sx q[0];
rz(2.572686) q[0];
rz(2.3078602) q[1];
sx q[1];
rz(-1.8582452) q[1];
sx q[1];
rz(-1.0047097) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92623989) q[0];
sx q[0];
rz(-2.4808501) q[0];
sx q[0];
rz(2.8244429) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9453853) q[2];
sx q[2];
rz(-0.56972144) q[2];
sx q[2];
rz(-1.8853055) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6697996) q[1];
sx q[1];
rz(-2.2679015) q[1];
sx q[1];
rz(-3.0473188) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58752693) q[3];
sx q[3];
rz(-1.6391132) q[3];
sx q[3];
rz(-2.5053566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7052762) q[2];
sx q[2];
rz(-2.269561) q[2];
sx q[2];
rz(-0.3479859) q[2];
rz(0.82585382) q[3];
sx q[3];
rz(-2.9233942) q[3];
sx q[3];
rz(-2.5364449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0572877) q[0];
sx q[0];
rz(-2.075752) q[0];
sx q[0];
rz(-0.2070981) q[0];
rz(2.6026978) q[1];
sx q[1];
rz(-0.62108827) q[1];
sx q[1];
rz(-2.5394687) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3054863) q[0];
sx q[0];
rz(-2.0417111) q[0];
sx q[0];
rz(-0.54022809) q[0];
x q[1];
rz(-1.978157) q[2];
sx q[2];
rz(-1.0216139) q[2];
sx q[2];
rz(2.801535) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4796743) q[1];
sx q[1];
rz(-1.9768324) q[1];
sx q[1];
rz(3.0597294) q[1];
x q[2];
rz(2.3942457) q[3];
sx q[3];
rz(-1.3266194) q[3];
sx q[3];
rz(-1.3507193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.78486097) q[2];
sx q[2];
rz(-2.1897903) q[2];
sx q[2];
rz(-1.3447364) q[2];
rz(0.52062672) q[3];
sx q[3];
rz(-2.8713363) q[3];
sx q[3];
rz(2.3394913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6945334) q[0];
sx q[0];
rz(-0.72493989) q[0];
sx q[0];
rz(-1.3865393) q[0];
rz(0.78320349) q[1];
sx q[1];
rz(-1.422311) q[1];
sx q[1];
rz(1.5625988) q[1];
rz(0.34306768) q[2];
sx q[2];
rz(-1.3160327) q[2];
sx q[2];
rz(2.119488) q[2];
rz(0.047719638) q[3];
sx q[3];
rz(-2.9012091) q[3];
sx q[3];
rz(-0.11500959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

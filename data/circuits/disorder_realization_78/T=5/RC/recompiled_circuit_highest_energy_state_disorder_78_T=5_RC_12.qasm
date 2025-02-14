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
rz(-1.6260835) q[0];
sx q[0];
rz(-0.27684394) q[0];
sx q[0];
rz(-3.0247363) q[0];
rz(2.6729743) q[1];
sx q[1];
rz(2.7873971) q[1];
sx q[1];
rz(8.9335557) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9899333) q[0];
sx q[0];
rz(-2.058968) q[0];
sx q[0];
rz(3.0753646) q[0];
x q[1];
rz(-2.6543219) q[2];
sx q[2];
rz(-1.3830162) q[2];
sx q[2];
rz(-2.1405381) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.74581214) q[1];
sx q[1];
rz(-0.94008259) q[1];
sx q[1];
rz(1.4480026) q[1];
rz(-pi) q[2];
rz(3.1163636) q[3];
sx q[3];
rz(-1.9140035) q[3];
sx q[3];
rz(2.0398839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7626071) q[2];
sx q[2];
rz(-2.8540322) q[2];
sx q[2];
rz(2.911705) q[2];
rz(3.0843132) q[3];
sx q[3];
rz(-0.66904896) q[3];
sx q[3];
rz(0.91276401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463992) q[0];
sx q[0];
rz(-2.754358) q[0];
sx q[0];
rz(-0.21695319) q[0];
rz(0.03288658) q[1];
sx q[1];
rz(-0.93371987) q[1];
sx q[1];
rz(-0.65509534) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0467619) q[0];
sx q[0];
rz(-2.3944986) q[0];
sx q[0];
rz(-1.2741035) q[0];
rz(-pi) q[1];
rz(-2.7883165) q[2];
sx q[2];
rz(-0.52918079) q[2];
sx q[2];
rz(0.25111094) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.078121878) q[1];
sx q[1];
rz(-0.95032108) q[1];
sx q[1];
rz(-0.93777754) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3052196) q[3];
sx q[3];
rz(-1.4640088) q[3];
sx q[3];
rz(-1.9143145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32404262) q[2];
sx q[2];
rz(-2.5174759) q[2];
sx q[2];
rz(1.8435271) q[2];
rz(1.921418) q[3];
sx q[3];
rz(-1.7036567) q[3];
sx q[3];
rz(2.4841681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.87797457) q[0];
sx q[0];
rz(-2.2442696) q[0];
sx q[0];
rz(2.5275912) q[0];
rz(1.3677431) q[1];
sx q[1];
rz(-0.65679336) q[1];
sx q[1];
rz(0.49555379) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14392631) q[0];
sx q[0];
rz(-0.94519061) q[0];
sx q[0];
rz(2.1493874) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96769963) q[2];
sx q[2];
rz(-2.357238) q[2];
sx q[2];
rz(-1.3788144) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6537955) q[1];
sx q[1];
rz(-1.3843517) q[1];
sx q[1];
rz(1.814278) q[1];
x q[2];
rz(0.63276799) q[3];
sx q[3];
rz(-2.583021) q[3];
sx q[3];
rz(-2.4726923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.76501781) q[2];
sx q[2];
rz(-2.0707264) q[2];
sx q[2];
rz(-2.5490254) q[2];
rz(-0.23558922) q[3];
sx q[3];
rz(-2.3841136) q[3];
sx q[3];
rz(-0.45430115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.4950824) q[0];
sx q[0];
rz(-2.3401234) q[0];
sx q[0];
rz(0.0058280514) q[0];
rz(2.5616772) q[1];
sx q[1];
rz(-0.35570759) q[1];
sx q[1];
rz(-2.6204806) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65893764) q[0];
sx q[0];
rz(-1.5636496) q[0];
sx q[0];
rz(3.1369476) q[0];
x q[1];
rz(0.02946202) q[2];
sx q[2];
rz(-0.97004393) q[2];
sx q[2];
rz(3.1371869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.24758928) q[1];
sx q[1];
rz(-2.517068) q[1];
sx q[1];
rz(-2.4127059) q[1];
rz(-2.8981325) q[3];
sx q[3];
rz(-0.32103466) q[3];
sx q[3];
rz(2.3784172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8899322) q[2];
sx q[2];
rz(-2.6106788) q[2];
sx q[2];
rz(0.31822515) q[2];
rz(-2.4288154) q[3];
sx q[3];
rz(-2.8409499) q[3];
sx q[3];
rz(-1.4819063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7707959) q[0];
sx q[0];
rz(-2.9279121) q[0];
sx q[0];
rz(0.036291432) q[0];
rz(-2.6151784) q[1];
sx q[1];
rz(-1.6830091) q[1];
sx q[1];
rz(-2.859419) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9045893) q[0];
sx q[0];
rz(-2.1891199) q[0];
sx q[0];
rz(2.8008921) q[0];
rz(0.93499489) q[2];
sx q[2];
rz(-0.80342573) q[2];
sx q[2];
rz(-0.73551169) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5314922) q[1];
sx q[1];
rz(-0.86750114) q[1];
sx q[1];
rz(1.478757) q[1];
x q[2];
rz(2.3822278) q[3];
sx q[3];
rz(-0.96624422) q[3];
sx q[3];
rz(0.60455632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7669749) q[2];
sx q[2];
rz(-0.66156113) q[2];
sx q[2];
rz(-3.0401163) q[2];
rz(-2.1756419) q[3];
sx q[3];
rz(-0.92065293) q[3];
sx q[3];
rz(-3.0275893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58979708) q[0];
sx q[0];
rz(-2.9647201) q[0];
sx q[0];
rz(-1.1740603) q[0];
rz(0.38408285) q[1];
sx q[1];
rz(-1.9575155) q[1];
sx q[1];
rz(0.039965872) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9703107) q[0];
sx q[0];
rz(-1.2501646) q[0];
sx q[0];
rz(-0.21954222) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1388117) q[2];
sx q[2];
rz(-0.68772763) q[2];
sx q[2];
rz(-2.1124396) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3921622) q[1];
sx q[1];
rz(-2.431173) q[1];
sx q[1];
rz(-2.602052) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82620718) q[3];
sx q[3];
rz(-1.4109269) q[3];
sx q[3];
rz(1.9217971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2018955) q[2];
sx q[2];
rz(-0.59212089) q[2];
sx q[2];
rz(-0.8197909) q[2];
rz(-0.78153265) q[3];
sx q[3];
rz(-0.90618366) q[3];
sx q[3];
rz(1.3801105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49067295) q[0];
sx q[0];
rz(-2.1020205) q[0];
sx q[0];
rz(1.8248722) q[0];
rz(-1.3680178) q[1];
sx q[1];
rz(-1.9098234) q[1];
sx q[1];
rz(0.43693158) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7765771) q[0];
sx q[0];
rz(-1.3804589) q[0];
sx q[0];
rz(-0.037419407) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75114366) q[2];
sx q[2];
rz(-0.1607543) q[2];
sx q[2];
rz(3.1243665) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1679042) q[1];
sx q[1];
rz(-1.2120795) q[1];
sx q[1];
rz(0.78679797) q[1];
x q[2];
rz(2.715492) q[3];
sx q[3];
rz(-2.0668732) q[3];
sx q[3];
rz(-0.29447281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74066073) q[2];
sx q[2];
rz(-1.0793945) q[2];
sx q[2];
rz(0.13460049) q[2];
rz(1.2234737) q[3];
sx q[3];
rz(-1.7569907) q[3];
sx q[3];
rz(-1.9756165) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80158919) q[0];
sx q[0];
rz(-0.29186258) q[0];
sx q[0];
rz(2.959751) q[0];
rz(-2.8002411) q[1];
sx q[1];
rz(-1.1266484) q[1];
sx q[1];
rz(-0.15886074) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5349992) q[0];
sx q[0];
rz(-0.45289055) q[0];
sx q[0];
rz(0.15720982) q[0];
rz(2.5652021) q[2];
sx q[2];
rz(-0.49242556) q[2];
sx q[2];
rz(-2.3362291) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6312069) q[1];
sx q[1];
rz(-2.800731) q[1];
sx q[1];
rz(2.4426961) q[1];
rz(2.587593) q[3];
sx q[3];
rz(-1.7007728) q[3];
sx q[3];
rz(2.8838571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2813256) q[2];
sx q[2];
rz(-2.211326) q[2];
sx q[2];
rz(2.2567828) q[2];
rz(1.7224711) q[3];
sx q[3];
rz(-2.116674) q[3];
sx q[3];
rz(0.44986808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.849843) q[0];
sx q[0];
rz(-1.2717286) q[0];
sx q[0];
rz(2.1480609) q[0];
rz(1.3604856) q[1];
sx q[1];
rz(-0.88881701) q[1];
sx q[1];
rz(2.5885168) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4489514) q[0];
sx q[0];
rz(-1.7990489) q[0];
sx q[0];
rz(-0.3897764) q[0];
rz(-pi) q[1];
rz(-0.73846784) q[2];
sx q[2];
rz(-1.6203801) q[2];
sx q[2];
rz(2.2253583) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7375664) q[1];
sx q[1];
rz(-1.8282923) q[1];
sx q[1];
rz(-2.5016258) q[1];
rz(-2.870375) q[3];
sx q[3];
rz(-0.94908774) q[3];
sx q[3];
rz(0.78628899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20455827) q[2];
sx q[2];
rz(-2.5356346) q[2];
sx q[2];
rz(0.23703144) q[2];
rz(1.6009181) q[3];
sx q[3];
rz(-2.377244) q[3];
sx q[3];
rz(1.1698394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071526214) q[0];
sx q[0];
rz(-1.6306174) q[0];
sx q[0];
rz(-3.0226829) q[0];
rz(-0.43926829) q[1];
sx q[1];
rz(-1.8426789) q[1];
sx q[1];
rz(0.063442245) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6574888) q[0];
sx q[0];
rz(-2.2229337) q[0];
sx q[0];
rz(0.87638292) q[0];
x q[1];
rz(2.5173152) q[2];
sx q[2];
rz(-1.7221525) q[2];
sx q[2];
rz(0.69619149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1181142) q[1];
sx q[1];
rz(-2.8148068) q[1];
sx q[1];
rz(-1.5972196) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2584553) q[3];
sx q[3];
rz(-2.4534907) q[3];
sx q[3];
rz(-2.7795252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.1305609) q[2];
sx q[2];
rz(-1.0728027) q[2];
sx q[2];
rz(1.4981184) q[2];
rz(-0.9324075) q[3];
sx q[3];
rz(-2.0216209) q[3];
sx q[3];
rz(-1.2137265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3907923) q[0];
sx q[0];
rz(-1.9600497) q[0];
sx q[0];
rz(2.0586769) q[0];
rz(2.2687268) q[1];
sx q[1];
rz(-2.0582336) q[1];
sx q[1];
rz(2.2179926) q[1];
rz(2.4058061) q[2];
sx q[2];
rz(-0.83333165) q[2];
sx q[2];
rz(0.32902645) q[2];
rz(-3.031293) q[3];
sx q[3];
rz(-0.58021373) q[3];
sx q[3];
rz(1.8662682) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

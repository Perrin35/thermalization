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
rz(1.5155091) q[0];
sx q[0];
rz(3.4184366) q[0];
sx q[0];
rz(9.3079216) q[0];
rz(-0.46861831) q[1];
sx q[1];
rz(-2.7873971) q[1];
sx q[1];
rz(2.6503704) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4502342) q[0];
sx q[0];
rz(-1.629279) q[0];
sx q[0];
rz(-2.0598777) q[0];
x q[1];
rz(-2.6543219) q[2];
sx q[2];
rz(-1.3830162) q[2];
sx q[2];
rz(-2.1405381) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.95211103) q[1];
sx q[1];
rz(-0.64095381) q[1];
sx q[1];
rz(0.16619311) q[1];
x q[2];
rz(-1.9141044) q[3];
sx q[3];
rz(-1.5945537) q[3];
sx q[3];
rz(-0.46059617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.7626071) q[2];
sx q[2];
rz(-2.8540322) q[2];
sx q[2];
rz(-0.22988764) q[2];
rz(0.057279438) q[3];
sx q[3];
rz(-2.4725437) q[3];
sx q[3];
rz(-2.2288286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463992) q[0];
sx q[0];
rz(-0.38723463) q[0];
sx q[0];
rz(-2.9246395) q[0];
rz(0.03288658) q[1];
sx q[1];
rz(-2.2078728) q[1];
sx q[1];
rz(0.65509534) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65194535) q[0];
sx q[0];
rz(-2.2781452) q[0];
sx q[0];
rz(2.8771557) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3711602) q[2];
sx q[2];
rz(-2.0642274) q[2];
sx q[2];
rz(-0.1525998) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8959143) q[1];
sx q[1];
rz(-2.0729471) q[1];
sx q[1];
rz(-2.4163626) q[1];
x q[2];
rz(2.9981587) q[3];
sx q[3];
rz(-2.3000882) q[3];
sx q[3];
rz(0.4394596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32404262) q[2];
sx q[2];
rz(-0.62411672) q[2];
sx q[2];
rz(1.8435271) q[2];
rz(1.2201747) q[3];
sx q[3];
rz(-1.4379359) q[3];
sx q[3];
rz(-0.65742457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2636181) q[0];
sx q[0];
rz(-2.2442696) q[0];
sx q[0];
rz(-0.61400145) q[0];
rz(1.7738495) q[1];
sx q[1];
rz(-2.4847993) q[1];
sx q[1];
rz(0.49555379) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7921831) q[0];
sx q[0];
rz(-2.0300477) q[0];
sx q[0];
rz(0.71190603) q[0];
rz(-2.2587288) q[2];
sx q[2];
rz(-1.9830215) q[2];
sx q[2];
rz(-0.26160535) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55807251) q[1];
sx q[1];
rz(-0.30553833) q[1];
sx q[1];
rz(-2.2347441) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5088247) q[3];
sx q[3];
rz(-0.55857165) q[3];
sx q[3];
rz(0.66890034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76501781) q[2];
sx q[2];
rz(-1.0708662) q[2];
sx q[2];
rz(2.5490254) q[2];
rz(-0.23558922) q[3];
sx q[3];
rz(-0.75747907) q[3];
sx q[3];
rz(0.45430115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6465103) q[0];
sx q[0];
rz(-2.3401234) q[0];
sx q[0];
rz(0.0058280514) q[0];
rz(2.5616772) q[1];
sx q[1];
rz(-0.35570759) q[1];
sx q[1];
rz(0.52111202) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.482655) q[0];
sx q[0];
rz(-1.577943) q[0];
sx q[0];
rz(0.004645017) q[0];
rz(-pi) q[1];
rz(1.5278339) q[2];
sx q[2];
rz(-0.60138541) q[2];
sx q[2];
rz(-0.056494519) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.24758928) q[1];
sx q[1];
rz(-0.62452468) q[1];
sx q[1];
rz(0.7288868) q[1];
rz(-pi) q[2];
x q[2];
rz(1.490805) q[3];
sx q[3];
rz(-1.2595525) q[3];
sx q[3];
rz(2.1224006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2516605) q[2];
sx q[2];
rz(-0.53091383) q[2];
sx q[2];
rz(0.31822515) q[2];
rz(0.71277726) q[3];
sx q[3];
rz(-0.30064279) q[3];
sx q[3];
rz(1.4819063) q[3];
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
rz(2.7707959) q[0];
sx q[0];
rz(-0.21368055) q[0];
sx q[0];
rz(0.036291432) q[0];
rz(0.52641422) q[1];
sx q[1];
rz(-1.6830091) q[1];
sx q[1];
rz(-2.859419) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.605114) q[0];
sx q[0];
rz(-1.295034) q[0];
sx q[0];
rz(-2.2173397) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55183175) q[2];
sx q[2];
rz(-0.95317344) q[2];
sx q[2];
rz(-3.0611401) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4683405) q[1];
sx q[1];
rz(-2.4333204) q[1];
sx q[1];
rz(0.10797031) q[1];
rz(2.3544694) q[3];
sx q[3];
rz(-0.93138444) q[3];
sx q[3];
rz(-1.5057664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7669749) q[2];
sx q[2];
rz(-2.4800315) q[2];
sx q[2];
rz(0.10147632) q[2];
rz(-2.1756419) q[3];
sx q[3];
rz(-0.92065293) q[3];
sx q[3];
rz(0.11400338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5517956) q[0];
sx q[0];
rz(-0.17687251) q[0];
sx q[0];
rz(-1.1740603) q[0];
rz(2.7575098) q[1];
sx q[1];
rz(-1.1840772) q[1];
sx q[1];
rz(0.039965872) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.171282) q[0];
sx q[0];
rz(-1.2501646) q[0];
sx q[0];
rz(2.9220504) q[0];
x q[1];
rz(-2.2117679) q[2];
sx q[2];
rz(-1.3017968) q[2];
sx q[2];
rz(-0.19942927) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.080885305) q[1];
sx q[1];
rz(-0.97700143) q[1];
sx q[1];
rz(1.9869366) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82620718) q[3];
sx q[3];
rz(-1.7306657) q[3];
sx q[3];
rz(1.2197956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2018955) q[2];
sx q[2];
rz(-0.59212089) q[2];
sx q[2];
rz(-2.3218018) q[2];
rz(0.78153265) q[3];
sx q[3];
rz(-2.235409) q[3];
sx q[3];
rz(-1.7614822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6509197) q[0];
sx q[0];
rz(-2.1020205) q[0];
sx q[0];
rz(1.8248722) q[0];
rz(1.7735749) q[1];
sx q[1];
rz(-1.9098234) q[1];
sx q[1];
rz(-2.7046611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56037134) q[0];
sx q[0];
rz(-0.19393714) q[0];
sx q[0];
rz(1.7625798) q[0];
rz(2.390449) q[2];
sx q[2];
rz(-0.1607543) q[2];
sx q[2];
rz(-3.1243665) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26035151) q[1];
sx q[1];
rz(-0.8484183) q[1];
sx q[1];
rz(-2.6546225) q[1];
rz(-0.42610069) q[3];
sx q[3];
rz(-2.0668732) q[3];
sx q[3];
rz(-0.29447281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74066073) q[2];
sx q[2];
rz(-2.0621982) q[2];
sx q[2];
rz(-0.13460049) q[2];
rz(-1.918119) q[3];
sx q[3];
rz(-1.7569907) q[3];
sx q[3];
rz(1.1659762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3400035) q[0];
sx q[0];
rz(-0.29186258) q[0];
sx q[0];
rz(2.959751) q[0];
rz(-0.34135154) q[1];
sx q[1];
rz(-1.1266484) q[1];
sx q[1];
rz(-2.9827319) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0358064) q[0];
sx q[0];
rz(-1.6393568) q[0];
sx q[0];
rz(0.44802702) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42271013) q[2];
sx q[2];
rz(-1.8313932) q[2];
sx q[2];
rz(2.8963331) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5103858) q[1];
sx q[1];
rz(-2.800731) q[1];
sx q[1];
rz(2.4426961) q[1];
rz(-pi) q[2];
rz(-1.7233061) q[3];
sx q[3];
rz(-1.0220064) q[3];
sx q[3];
rz(-1.9085409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.86026704) q[2];
sx q[2];
rz(-0.93026668) q[2];
sx q[2];
rz(-2.2567828) q[2];
rz(1.7224711) q[3];
sx q[3];
rz(-2.116674) q[3];
sx q[3];
rz(-2.6917246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.849843) q[0];
sx q[0];
rz(-1.869864) q[0];
sx q[0];
rz(2.1480609) q[0];
rz(-1.3604856) q[1];
sx q[1];
rz(-0.88881701) q[1];
sx q[1];
rz(-2.5885168) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6926413) q[0];
sx q[0];
rz(-1.3425437) q[0];
sx q[0];
rz(-0.3897764) q[0];
rz(2.4031248) q[2];
sx q[2];
rz(-1.6203801) q[2];
sx q[2];
rz(-0.91623437) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7375664) q[1];
sx q[1];
rz(-1.8282923) q[1];
sx q[1];
rz(-0.63996686) q[1];
rz(-pi) q[2];
x q[2];
rz(2.210238) q[3];
sx q[3];
rz(-1.7903312) q[3];
sx q[3];
rz(-0.94506782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9370344) q[2];
sx q[2];
rz(-0.60595804) q[2];
sx q[2];
rz(-2.9045612) q[2];
rz(1.6009181) q[3];
sx q[3];
rz(-2.377244) q[3];
sx q[3];
rz(-1.9717533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(3.0700664) q[0];
sx q[0];
rz(-1.6306174) q[0];
sx q[0];
rz(-3.0226829) q[0];
rz(0.43926829) q[1];
sx q[1];
rz(-1.8426789) q[1];
sx q[1];
rz(-0.063442245) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.760323) q[0];
sx q[0];
rz(-2.1043692) q[0];
sx q[0];
rz(-0.78223451) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5173152) q[2];
sx q[2];
rz(-1.4194402) q[2];
sx q[2];
rz(2.4454012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0513775) q[1];
sx q[1];
rz(-1.8974638) q[1];
sx q[1];
rz(0.0089545687) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88313734) q[3];
sx q[3];
rz(-0.68810191) q[3];
sx q[3];
rz(2.7795252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1305609) q[2];
sx q[2];
rz(-2.06879) q[2];
sx q[2];
rz(1.6434742) q[2];
rz(0.9324075) q[3];
sx q[3];
rz(-2.0216209) q[3];
sx q[3];
rz(1.2137265) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7508004) q[0];
sx q[0];
rz(-1.9600497) q[0];
sx q[0];
rz(2.0586769) q[0];
rz(-0.87286585) q[1];
sx q[1];
rz(-2.0582336) q[1];
sx q[1];
rz(2.2179926) q[1];
rz(0.93449705) q[2];
sx q[2];
rz(-2.1515982) q[2];
sx q[2];
rz(1.2610255) q[2];
rz(0.11029966) q[3];
sx q[3];
rz(-0.58021373) q[3];
sx q[3];
rz(1.8662682) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

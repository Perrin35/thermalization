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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15165937) q[0];
sx q[0];
rz(-2.058968) q[0];
sx q[0];
rz(-3.0753646) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7826177) q[2];
sx q[2];
rz(-1.0928177) q[2];
sx q[2];
rz(2.4732531) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95211103) q[1];
sx q[1];
rz(-2.5006388) q[1];
sx q[1];
rz(2.9753995) q[1];
rz(-pi) q[2];
rz(-3.1163636) q[3];
sx q[3];
rz(-1.2275891) q[3];
sx q[3];
rz(2.0398839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.7626071) q[2];
sx q[2];
rz(-2.8540322) q[2];
sx q[2];
rz(-0.22988764) q[2];
rz(3.0843132) q[3];
sx q[3];
rz(-2.4725437) q[3];
sx q[3];
rz(2.2288286) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463992) q[0];
sx q[0];
rz(-0.38723463) q[0];
sx q[0];
rz(-2.9246395) q[0];
rz(3.1087061) q[1];
sx q[1];
rz(-0.93371987) q[1];
sx q[1];
rz(-2.4864973) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74467662) q[0];
sx q[0];
rz(-1.370805) q[0];
sx q[0];
rz(-2.295664) q[0];
rz(-pi) q[1];
rz(1.7704324) q[2];
sx q[2];
rz(-2.0642274) q[2];
sx q[2];
rz(-0.1525998) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.82258525) q[1];
sx q[1];
rz(-2.2862541) q[1];
sx q[1];
rz(2.4501178) q[1];
rz(-pi) q[2];
rz(-2.3052196) q[3];
sx q[3];
rz(-1.6775838) q[3];
sx q[3];
rz(-1.9143145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.81755) q[2];
sx q[2];
rz(-0.62411672) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87797457) q[0];
sx q[0];
rz(-2.2442696) q[0];
sx q[0];
rz(0.61400145) q[0];
rz(1.3677431) q[1];
sx q[1];
rz(-2.4847993) q[1];
sx q[1];
rz(-0.49555379) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9976663) q[0];
sx q[0];
rz(-2.196402) q[0];
sx q[0];
rz(0.99220522) q[0];
x q[1];
rz(-2.2587288) q[2];
sx q[2];
rz(-1.1585711) q[2];
sx q[2];
rz(-2.8799873) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6537955) q[1];
sx q[1];
rz(-1.7572409) q[1];
sx q[1];
rz(-1.814278) q[1];
x q[2];
rz(-1.2167778) q[3];
sx q[3];
rz(-2.0123768) q[3];
sx q[3];
rz(1.7596678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3765748) q[2];
sx q[2];
rz(-2.0707264) q[2];
sx q[2];
rz(-2.5490254) q[2];
rz(-2.9060034) q[3];
sx q[3];
rz(-0.75747907) q[3];
sx q[3];
rz(-0.45430115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6465103) q[0];
sx q[0];
rz(-0.80146924) q[0];
sx q[0];
rz(3.1357646) q[0];
rz(-2.5616772) q[1];
sx q[1];
rz(-0.35570759) q[1];
sx q[1];
rz(-0.52111202) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0590044) q[0];
sx q[0];
rz(-3.1330691) q[0];
sx q[0];
rz(0.99446358) q[0];
rz(0.96984152) q[2];
sx q[2];
rz(-1.5464939) q[2];
sx q[2];
rz(1.5497335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24758928) q[1];
sx q[1];
rz(-0.62452468) q[1];
sx q[1];
rz(-0.7288868) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.490805) q[3];
sx q[3];
rz(-1.8820401) q[3];
sx q[3];
rz(-1.019192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8899322) q[2];
sx q[2];
rz(-0.53091383) q[2];
sx q[2];
rz(0.31822515) q[2];
rz(0.71277726) q[3];
sx q[3];
rz(-2.8409499) q[3];
sx q[3];
rz(-1.4819063) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37079674) q[0];
sx q[0];
rz(-2.9279121) q[0];
sx q[0];
rz(0.036291432) q[0];
rz(-0.52641422) q[1];
sx q[1];
rz(-1.6830091) q[1];
sx q[1];
rz(-0.2821736) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68809915) q[0];
sx q[0];
rz(-0.69506139) q[0];
sx q[0];
rz(-1.1316677) q[0];
rz(0.93499489) q[2];
sx q[2];
rz(-0.80342573) q[2];
sx q[2];
rz(-0.73551169) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6101004) q[1];
sx q[1];
rz(-0.86750114) q[1];
sx q[1];
rz(1.6628357) q[1];
rz(-pi) q[2];
rz(0.7871233) q[3];
sx q[3];
rz(-2.2102082) q[3];
sx q[3];
rz(-1.5057664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37461773) q[2];
sx q[2];
rz(-2.4800315) q[2];
sx q[2];
rz(0.10147632) q[2];
rz(-2.1756419) q[3];
sx q[3];
rz(-0.92065293) q[3];
sx q[3];
rz(-3.0275893) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5517956) q[0];
sx q[0];
rz(-2.9647201) q[0];
sx q[0];
rz(-1.9675323) q[0];
rz(2.7575098) q[1];
sx q[1];
rz(-1.9575155) q[1];
sx q[1];
rz(-0.039965872) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46972457) q[0];
sx q[0];
rz(-1.7789808) q[0];
sx q[0];
rz(1.2428268) q[0];
x q[1];
rz(-2.8103175) q[2];
sx q[2];
rz(-2.185198) q[2];
sx q[2];
rz(-1.5671052) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.080885305) q[1];
sx q[1];
rz(-2.1645912) q[1];
sx q[1];
rz(-1.9869366) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9257366) q[3];
sx q[3];
rz(-0.83789589) q[3];
sx q[3];
rz(0.49666109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.93969718) q[2];
sx q[2];
rz(-0.59212089) q[2];
sx q[2];
rz(-2.3218018) q[2];
rz(-0.78153265) q[3];
sx q[3];
rz(-2.235409) q[3];
sx q[3];
rz(-1.3801105) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49067295) q[0];
sx q[0];
rz(-2.1020205) q[0];
sx q[0];
rz(1.3167205) q[0];
rz(1.7735749) q[1];
sx q[1];
rz(-1.9098234) q[1];
sx q[1];
rz(-2.7046611) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1986982) q[0];
sx q[0];
rz(-1.6075396) q[0];
sx q[0];
rz(1.7612639) q[0];
rz(-pi) q[1];
rz(-1.4605791) q[2];
sx q[2];
rz(-1.6880562) q[2];
sx q[2];
rz(2.3667468) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97368846) q[1];
sx q[1];
rz(-1.9295132) q[1];
sx q[1];
rz(0.78679797) q[1];
rz(0.42610069) q[3];
sx q[3];
rz(-2.0668732) q[3];
sx q[3];
rz(-2.8471198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.74066073) q[2];
sx q[2];
rz(-2.0621982) q[2];
sx q[2];
rz(3.0069922) q[2];
rz(1.918119) q[3];
sx q[3];
rz(-1.384602) q[3];
sx q[3];
rz(1.1659762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3400035) q[0];
sx q[0];
rz(-2.8497301) q[0];
sx q[0];
rz(0.18184161) q[0];
rz(-2.8002411) q[1];
sx q[1];
rz(-2.0149442) q[1];
sx q[1];
rz(-2.9827319) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0358064) q[0];
sx q[0];
rz(-1.6393568) q[0];
sx q[0];
rz(2.6935656) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7188825) q[2];
sx q[2];
rz(-1.8313932) q[2];
sx q[2];
rz(-2.8963331) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6312069) q[1];
sx q[1];
rz(-2.800731) q[1];
sx q[1];
rz(0.69889654) q[1];
rz(-pi) q[2];
rz(2.587593) q[3];
sx q[3];
rz(-1.4408198) q[3];
sx q[3];
rz(-2.8838571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2813256) q[2];
sx q[2];
rz(-0.93026668) q[2];
sx q[2];
rz(-2.2567828) q[2];
rz(1.4191215) q[3];
sx q[3];
rz(-1.0249187) q[3];
sx q[3];
rz(-2.6917246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.849843) q[0];
sx q[0];
rz(-1.2717286) q[0];
sx q[0];
rz(0.9935317) q[0];
rz(-1.3604856) q[1];
sx q[1];
rz(-0.88881701) q[1];
sx q[1];
rz(0.55307585) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9270614) q[0];
sx q[0];
rz(-1.1916516) q[0];
sx q[0];
rz(1.3247471) q[0];
x q[1];
rz(3.0680066) q[2];
sx q[2];
rz(-2.4017757) q[2];
sx q[2];
rz(0.70895665) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.83718761) q[1];
sx q[1];
rz(-0.68301979) q[1];
sx q[1];
rz(-0.4153312) q[1];
rz(-pi) q[2];
rz(2.870375) q[3];
sx q[3];
rz(-2.1925049) q[3];
sx q[3];
rz(0.78628899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9370344) q[2];
sx q[2];
rz(-2.5356346) q[2];
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
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0700664) q[0];
sx q[0];
rz(-1.6306174) q[0];
sx q[0];
rz(-0.1189098) q[0];
rz(2.7023244) q[1];
sx q[1];
rz(-1.2989137) q[1];
sx q[1];
rz(-0.063442245) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3812696) q[0];
sx q[0];
rz(-2.1043692) q[0];
sx q[0];
rz(0.78223451) q[0];
x q[1];
rz(0.25524885) q[2];
sx q[2];
rz(-0.63997696) q[2];
sx q[2];
rz(2.4733123) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0234785) q[1];
sx q[1];
rz(-2.8148068) q[1];
sx q[1];
rz(1.5972196) q[1];
rz(-2.2584553) q[3];
sx q[3];
rz(-0.68810191) q[3];
sx q[3];
rz(-0.36206743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0110317) q[2];
sx q[2];
rz(-2.06879) q[2];
sx q[2];
rz(-1.4981184) q[2];
rz(-2.2091852) q[3];
sx q[3];
rz(-1.1199718) q[3];
sx q[3];
rz(-1.2137265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-0.93449705) q[2];
sx q[2];
rz(-0.9899944) q[2];
sx q[2];
rz(-1.8805671) q[2];
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

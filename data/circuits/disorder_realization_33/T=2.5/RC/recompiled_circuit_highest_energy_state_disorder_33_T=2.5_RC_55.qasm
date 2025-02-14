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
rz(-2.9774732) q[0];
sx q[0];
rz(-1.0241221) q[0];
sx q[0];
rz(-0.48409387) q[0];
rz(2.3789499) q[1];
sx q[1];
rz(-1.5634544) q[1];
sx q[1];
rz(-0.054952316) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8380175) q[0];
sx q[0];
rz(-0.9009046) q[0];
sx q[0];
rz(-0.039867) q[0];
rz(3.0829551) q[2];
sx q[2];
rz(-1.8203041) q[2];
sx q[2];
rz(-0.80896689) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.311885) q[1];
sx q[1];
rz(-1.4730318) q[1];
sx q[1];
rz(3.0086161) q[1];
x q[2];
rz(2.1922429) q[3];
sx q[3];
rz(-1.2934522) q[3];
sx q[3];
rz(0.30888939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.83076465) q[2];
sx q[2];
rz(-2.1711633) q[2];
sx q[2];
rz(-0.27466276) q[2];
rz(2.2667387) q[3];
sx q[3];
rz(-2.0503876) q[3];
sx q[3];
rz(-2.8680475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.4817568) q[0];
sx q[0];
rz(-0.6898703) q[0];
sx q[0];
rz(2.7428395) q[0];
rz(-2.8975471) q[1];
sx q[1];
rz(-1.2363385) q[1];
sx q[1];
rz(-2.0358548) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2139521) q[0];
sx q[0];
rz(-2.0114779) q[0];
sx q[0];
rz(-1.8233612) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91715889) q[2];
sx q[2];
rz(-1.0984525) q[2];
sx q[2];
rz(1.2231959) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.51197375) q[1];
sx q[1];
rz(-1.5675768) q[1];
sx q[1];
rz(-1.5734488) q[1];
rz(-pi) q[2];
rz(3.0209474) q[3];
sx q[3];
rz(-1.835726) q[3];
sx q[3];
rz(0.93075965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4411321) q[2];
sx q[2];
rz(-1.9596142) q[2];
sx q[2];
rz(2.0665118) q[2];
rz(-2.4454146) q[3];
sx q[3];
rz(-1.8680365) q[3];
sx q[3];
rz(0.16304326) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27580801) q[0];
sx q[0];
rz(-1.2514665) q[0];
sx q[0];
rz(-0.21670093) q[0];
rz(-1.7614583) q[1];
sx q[1];
rz(-2.7237027) q[1];
sx q[1];
rz(-0.61947852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9443837) q[0];
sx q[0];
rz(-1.608662) q[0];
sx q[0];
rz(0.029649563) q[0];
rz(1.6939075) q[2];
sx q[2];
rz(-1.6173714) q[2];
sx q[2];
rz(-1.2706626) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6782376) q[1];
sx q[1];
rz(-1.2690407) q[1];
sx q[1];
rz(-0.22308087) q[1];
x q[2];
rz(0.61010078) q[3];
sx q[3];
rz(-1.0786295) q[3];
sx q[3];
rz(-0.51612332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2565903) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(-0.090506434) q[2];
rz(-2.6835486) q[3];
sx q[3];
rz(-2.6543255) q[3];
sx q[3];
rz(-1.1080144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3312382) q[0];
sx q[0];
rz(-2.2585223) q[0];
sx q[0];
rz(2.2099387) q[0];
rz(-0.16904198) q[1];
sx q[1];
rz(-0.5642429) q[1];
sx q[1];
rz(-1.0394675) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66327205) q[0];
sx q[0];
rz(-2.2815858) q[0];
sx q[0];
rz(1.8904314) q[0];
rz(-pi) q[1];
rz(1.3170856) q[2];
sx q[2];
rz(-1.3519962) q[2];
sx q[2];
rz(-1.7148866) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9586981) q[1];
sx q[1];
rz(-2.9258203) q[1];
sx q[1];
rz(1.853626) q[1];
rz(2.5040197) q[3];
sx q[3];
rz(-2.7927783) q[3];
sx q[3];
rz(1.82774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52183759) q[2];
sx q[2];
rz(-2.0401185) q[2];
sx q[2];
rz(-2.6452765) q[2];
rz(-0.79658341) q[3];
sx q[3];
rz(-0.28592548) q[3];
sx q[3];
rz(-1.0930141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.8747044) q[0];
sx q[0];
rz(-0.58458352) q[0];
sx q[0];
rz(-2.1697178) q[0];
rz(3.019849) q[1];
sx q[1];
rz(-1.5309445) q[1];
sx q[1];
rz(-1.8121388) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0861391) q[0];
sx q[0];
rz(-0.3124961) q[0];
sx q[0];
rz(-1.7715095) q[0];
rz(-pi) q[1];
rz(2.0636286) q[2];
sx q[2];
rz(-1.9415641) q[2];
sx q[2];
rz(2.855122) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5962474) q[1];
sx q[1];
rz(-1.9831428) q[1];
sx q[1];
rz(-0.72280563) q[1];
rz(-1.6498783) q[3];
sx q[3];
rz(-1.5800161) q[3];
sx q[3];
rz(-0.10445933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.04008) q[2];
sx q[2];
rz(-1.2206581) q[2];
sx q[2];
rz(0.15138781) q[2];
rz(0.76465145) q[3];
sx q[3];
rz(-2.305856) q[3];
sx q[3];
rz(-1.5449272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0720035) q[0];
sx q[0];
rz(-1.66865) q[0];
sx q[0];
rz(-2.0701011) q[0];
rz(0.53120652) q[1];
sx q[1];
rz(-1.4899645) q[1];
sx q[1];
rz(2.5154617) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1086639) q[0];
sx q[0];
rz(-1.5877921) q[0];
sx q[0];
rz(-0.0047608382) q[0];
rz(-pi) q[1];
rz(1.8161723) q[2];
sx q[2];
rz(-0.98773709) q[2];
sx q[2];
rz(2.8193605) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6408566) q[1];
sx q[1];
rz(-1.8759512) q[1];
sx q[1];
rz(-0.63668294) q[1];
x q[2];
rz(-0.49895309) q[3];
sx q[3];
rz(-1.1756983) q[3];
sx q[3];
rz(-0.56624352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7834187) q[2];
sx q[2];
rz(-2.5219315) q[2];
sx q[2];
rz(-2.1273071) q[2];
rz(1.8861534) q[3];
sx q[3];
rz(-1.7858601) q[3];
sx q[3];
rz(2.0881418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1805873) q[0];
sx q[0];
rz(-0.87562457) q[0];
sx q[0];
rz(-2.2210806) q[0];
rz(3.0275184) q[1];
sx q[1];
rz(-1.736085) q[1];
sx q[1];
rz(2.0106409) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3217472) q[0];
sx q[0];
rz(-1.5764569) q[0];
sx q[0];
rz(-0.59771718) q[0];
x q[1];
rz(-2.1423856) q[2];
sx q[2];
rz(-2.0519369) q[2];
sx q[2];
rz(0.66758093) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1753833) q[1];
sx q[1];
rz(-1.4384603) q[1];
sx q[1];
rz(0.90465178) q[1];
rz(1.4907964) q[3];
sx q[3];
rz(-1.0558075) q[3];
sx q[3];
rz(1.1752216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7941234) q[2];
sx q[2];
rz(-2.4994714) q[2];
sx q[2];
rz(-2.7327909) q[2];
rz(-2.9583904) q[3];
sx q[3];
rz(-1.5513523) q[3];
sx q[3];
rz(2.0028152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0278397) q[0];
sx q[0];
rz(-3.0945859) q[0];
sx q[0];
rz(0.29079944) q[0];
rz(-0.94789061) q[1];
sx q[1];
rz(-0.46698505) q[1];
sx q[1];
rz(1.1999757) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24058293) q[0];
sx q[0];
rz(-1.9685317) q[0];
sx q[0];
rz(-0.57698864) q[0];
rz(0.84917111) q[2];
sx q[2];
rz(-0.90540041) q[2];
sx q[2];
rz(2.9862491) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51344752) q[1];
sx q[1];
rz(-1.2772868) q[1];
sx q[1];
rz(-0.048444466) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7539976) q[3];
sx q[3];
rz(-1.4716513) q[3];
sx q[3];
rz(0.20618901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7802508) q[2];
sx q[2];
rz(-1.4486518) q[2];
sx q[2];
rz(-1.7928436) q[2];
rz(1.5032984) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(-0.1851113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0272442) q[0];
sx q[0];
rz(-2.769727) q[0];
sx q[0];
rz(1.0741023) q[0];
rz(0.4298003) q[1];
sx q[1];
rz(-1.0440412) q[1];
sx q[1];
rz(2.6780186) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6493452) q[0];
sx q[0];
rz(-1.6657636) q[0];
sx q[0];
rz(-1.7596721) q[0];
rz(-1.3691291) q[2];
sx q[2];
rz(-1.8275765) q[2];
sx q[2];
rz(-0.078291206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8817655) q[1];
sx q[1];
rz(-2.1912088) q[1];
sx q[1];
rz(0.46807351) q[1];
x q[2];
rz(-0.55863278) q[3];
sx q[3];
rz(-1.6470634) q[3];
sx q[3];
rz(2.7879451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.27434906) q[2];
sx q[2];
rz(-2.5090019) q[2];
sx q[2];
rz(-2.3331433) q[2];
rz(-2.6570053) q[3];
sx q[3];
rz(-0.37823585) q[3];
sx q[3];
rz(0.43420473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7782068) q[0];
sx q[0];
rz(-0.45553842) q[0];
sx q[0];
rz(-1.1325915) q[0];
rz(-3.1032108) q[1];
sx q[1];
rz(-1.8033586) q[1];
sx q[1];
rz(-2.8448232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7420121) q[0];
sx q[0];
rz(-1.0431847) q[0];
sx q[0];
rz(2.8199838) q[0];
rz(2.0149258) q[2];
sx q[2];
rz(-1.9812366) q[2];
sx q[2];
rz(-3.0164312) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3546875) q[1];
sx q[1];
rz(-2.9070435) q[1];
sx q[1];
rz(-0.3292747) q[1];
rz(-1.1007916) q[3];
sx q[3];
rz(-0.74277011) q[3];
sx q[3];
rz(-2.7239885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9639637) q[2];
sx q[2];
rz(-2.0202426) q[2];
sx q[2];
rz(-0.56619823) q[2];
rz(1.1244134) q[3];
sx q[3];
rz(-3.0163613) q[3];
sx q[3];
rz(1.9521149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.88846702) q[0];
sx q[0];
rz(-0.69857004) q[0];
sx q[0];
rz(0.10733124) q[0];
rz(2.879907) q[1];
sx q[1];
rz(-1.9410004) q[1];
sx q[1];
rz(0.95071361) q[1];
rz(0.08200866) q[2];
sx q[2];
rz(-1.5463943) q[2];
sx q[2];
rz(0.80873185) q[2];
rz(-2.2812609) q[3];
sx q[3];
rz(-2.0900149) q[3];
sx q[3];
rz(2.7505977) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

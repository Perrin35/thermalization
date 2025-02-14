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
rz(1.264313) q[0];
sx q[0];
rz(-2.8180583) q[0];
sx q[0];
rz(-2.6927595) q[0];
rz(-1.1558865) q[1];
sx q[1];
rz(-1.7855676) q[1];
sx q[1];
rz(1.9670271) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4455152) q[0];
sx q[0];
rz(-2.4598274) q[0];
sx q[0];
rz(-2.5755139) q[0];
x q[1];
rz(1.6392133) q[2];
sx q[2];
rz(-1.8937773) q[2];
sx q[2];
rz(-0.078106192) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5331363) q[1];
sx q[1];
rz(-1.605833) q[1];
sx q[1];
rz(-0.94989454) q[1];
rz(-pi) q[2];
x q[2];
rz(1.470119) q[3];
sx q[3];
rz(-2.0824471) q[3];
sx q[3];
rz(2.4477521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.55149469) q[2];
sx q[2];
rz(-1.3000725) q[2];
sx q[2];
rz(-0.91280118) q[2];
rz(1.270795) q[3];
sx q[3];
rz(-2.1015034) q[3];
sx q[3];
rz(-0.84996581) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0047282334) q[0];
sx q[0];
rz(-0.57924634) q[0];
sx q[0];
rz(2.7194523) q[0];
rz(-0.36960754) q[1];
sx q[1];
rz(-1.0969176) q[1];
sx q[1];
rz(0.94013989) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62806118) q[0];
sx q[0];
rz(-1.2413176) q[0];
sx q[0];
rz(0.0026436289) q[0];
x q[1];
rz(1.0375848) q[2];
sx q[2];
rz(-1.4681446) q[2];
sx q[2];
rz(-0.29666049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9782422) q[1];
sx q[1];
rz(-0.85127318) q[1];
sx q[1];
rz(-2.4656328) q[1];
rz(-pi) q[2];
rz(0.89408447) q[3];
sx q[3];
rz(-1.7853123) q[3];
sx q[3];
rz(-1.0971951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9422354) q[2];
sx q[2];
rz(-0.67023674) q[2];
sx q[2];
rz(0.90636903) q[2];
rz(2.1622315) q[3];
sx q[3];
rz(-2.7370079) q[3];
sx q[3];
rz(-1.8766859) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1449428) q[0];
sx q[0];
rz(-3.0746089) q[0];
sx q[0];
rz(1.2445194) q[0];
rz(2.7527346) q[1];
sx q[1];
rz(-1.2797979) q[1];
sx q[1];
rz(0.32274524) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9661315) q[0];
sx q[0];
rz(-1.4585988) q[0];
sx q[0];
rz(-2.4787122) q[0];
rz(-pi) q[1];
rz(2.0409731) q[2];
sx q[2];
rz(-2.5289446) q[2];
sx q[2];
rz(2.8620811) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9936064) q[1];
sx q[1];
rz(-1.8022416) q[1];
sx q[1];
rz(2.3416421) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1348008) q[3];
sx q[3];
rz(-0.35405891) q[3];
sx q[3];
rz(1.1969138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2286223) q[2];
sx q[2];
rz(-1.7614438) q[2];
sx q[2];
rz(1.7886394) q[2];
rz(2.8690423) q[3];
sx q[3];
rz(-1.8504146) q[3];
sx q[3];
rz(-1.4636309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54101855) q[0];
sx q[0];
rz(-2.8193642) q[0];
sx q[0];
rz(-1.2433276) q[0];
rz(-0.86878949) q[1];
sx q[1];
rz(-0.96005762) q[1];
sx q[1];
rz(2.0713846) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7420827) q[0];
sx q[0];
rz(-0.94219452) q[0];
sx q[0];
rz(0.6357155) q[0];
rz(-2.0212428) q[2];
sx q[2];
rz(-2.3131158) q[2];
sx q[2];
rz(2.7591005) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94302109) q[1];
sx q[1];
rz(-2.7987438) q[1];
sx q[1];
rz(1.2238316) q[1];
rz(0.34262212) q[3];
sx q[3];
rz(-1.7512158) q[3];
sx q[3];
rz(-2.8460549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.59226817) q[2];
sx q[2];
rz(-0.78846875) q[2];
sx q[2];
rz(-2.2198524) q[2];
rz(2.8978469) q[3];
sx q[3];
rz(-1.7071525) q[3];
sx q[3];
rz(-0.29269472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9859966) q[0];
sx q[0];
rz(-1.3256185) q[0];
sx q[0];
rz(-0.50088125) q[0];
rz(2.0241375) q[1];
sx q[1];
rz(-1.3331579) q[1];
sx q[1];
rz(0.31455988) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.019032) q[0];
sx q[0];
rz(-1.9261596) q[0];
sx q[0];
rz(1.2073231) q[0];
x q[1];
rz(1.4681513) q[2];
sx q[2];
rz(-0.043641239) q[2];
sx q[2];
rz(-0.99422821) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3479473) q[1];
sx q[1];
rz(-2.1160908) q[1];
sx q[1];
rz(-2.6752363) q[1];
rz(-pi) q[2];
rz(1.3077478) q[3];
sx q[3];
rz(-2.7125053) q[3];
sx q[3];
rz(-1.3671631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8353117) q[2];
sx q[2];
rz(-1.1113144) q[2];
sx q[2];
rz(1.8446946) q[2];
rz(2.6563307) q[3];
sx q[3];
rz(-1.5596215) q[3];
sx q[3];
rz(1.1135134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6473963) q[0];
sx q[0];
rz(-1.2911456) q[0];
sx q[0];
rz(-0.099763481) q[0];
rz(1.2707204) q[1];
sx q[1];
rz(-2.2427509) q[1];
sx q[1];
rz(-1.7521923) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4609328) q[0];
sx q[0];
rz(-1.4062728) q[0];
sx q[0];
rz(2.8712832) q[0];
x q[1];
rz(1.4657137) q[2];
sx q[2];
rz(-1.3977524) q[2];
sx q[2];
rz(-1.6478754) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8956611) q[1];
sx q[1];
rz(-0.97890039) q[1];
sx q[1];
rz(-0.81879692) q[1];
rz(-3.0148195) q[3];
sx q[3];
rz(-0.78713464) q[3];
sx q[3];
rz(-1.8925557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.49326593) q[2];
sx q[2];
rz(-0.47097012) q[2];
sx q[2];
rz(0.89034447) q[2];
rz(-1.1912311) q[3];
sx q[3];
rz(-1.3072562) q[3];
sx q[3];
rz(2.2395535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4297727) q[0];
sx q[0];
rz(-0.044782488) q[0];
sx q[0];
rz(-0.032715948) q[0];
rz(0.50321594) q[1];
sx q[1];
rz(-1.0992173) q[1];
sx q[1];
rz(0.12282664) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5635934) q[0];
sx q[0];
rz(-1.2482572) q[0];
sx q[0];
rz(0.75046993) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5678504) q[2];
sx q[2];
rz(-0.57663585) q[2];
sx q[2];
rz(1.6078311) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87007755) q[1];
sx q[1];
rz(-2.6598499) q[1];
sx q[1];
rz(1.7863356) q[1];
x q[2];
rz(-1.7006247) q[3];
sx q[3];
rz(-1.9205838) q[3];
sx q[3];
rz(1.2799124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.34933019) q[2];
sx q[2];
rz(-2.1749039) q[2];
sx q[2];
rz(1.9897423) q[2];
rz(2.2564015) q[3];
sx q[3];
rz(-1.3693634) q[3];
sx q[3];
rz(-1.7159897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.952878) q[0];
sx q[0];
rz(-0.9335683) q[0];
sx q[0];
rz(-2.2123912) q[0];
rz(0.62905351) q[1];
sx q[1];
rz(-1.7864497) q[1];
sx q[1];
rz(-2.3984875) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0068664) q[0];
sx q[0];
rz(-1.6496302) q[0];
sx q[0];
rz(-2.0187442) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9628108) q[2];
sx q[2];
rz(-2.6465073) q[2];
sx q[2];
rz(-2.924233) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.23452246) q[1];
sx q[1];
rz(-1.6369695) q[1];
sx q[1];
rz(-2.727134) q[1];
rz(-0.943103) q[3];
sx q[3];
rz(-0.73929683) q[3];
sx q[3];
rz(0.044667808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5762081) q[2];
sx q[2];
rz(-0.78018633) q[2];
sx q[2];
rz(0.72594491) q[2];
rz(-2.7706326) q[3];
sx q[3];
rz(-1.5981263) q[3];
sx q[3];
rz(0.36271873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87042701) q[0];
sx q[0];
rz(-2.3112264) q[0];
sx q[0];
rz(0.56394947) q[0];
rz(1.14934) q[1];
sx q[1];
rz(-2.1795858) q[1];
sx q[1];
rz(2.3990778) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19631736) q[0];
sx q[0];
rz(-2.4565897) q[0];
sx q[0];
rz(2.3514868) q[0];
rz(-1.276539) q[2];
sx q[2];
rz(-1.9694491) q[2];
sx q[2];
rz(-0.63766232) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0418675) q[1];
sx q[1];
rz(-2.0085287) q[1];
sx q[1];
rz(1.5220286) q[1];
x q[2];
rz(2.2854076) q[3];
sx q[3];
rz(-2.0171229) q[3];
sx q[3];
rz(0.7577589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.62770647) q[2];
sx q[2];
rz(-1.0739001) q[2];
sx q[2];
rz(-0.34454301) q[2];
rz(-2.6978317) q[3];
sx q[3];
rz(-2.0659476) q[3];
sx q[3];
rz(-1.3226604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3798856) q[0];
sx q[0];
rz(-2.8590617) q[0];
sx q[0];
rz(0.66201061) q[0];
rz(2.8061891) q[1];
sx q[1];
rz(-1.6796203) q[1];
sx q[1];
rz(2.2047156) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5479845) q[0];
sx q[0];
rz(-2.4146705) q[0];
sx q[0];
rz(1.0127823) q[0];
x q[1];
rz(3.1179948) q[2];
sx q[2];
rz(-1.3006217) q[2];
sx q[2];
rz(-2.9058356) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7422694) q[1];
sx q[1];
rz(-1.284467) q[1];
sx q[1];
rz(2.7469357) q[1];
rz(-pi) q[2];
rz(0.85031894) q[3];
sx q[3];
rz(-0.26278824) q[3];
sx q[3];
rz(-0.14498728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4235437) q[2];
sx q[2];
rz(-0.23423883) q[2];
sx q[2];
rz(2.6424109) q[2];
rz(1.4623803) q[3];
sx q[3];
rz(-1.9951818) q[3];
sx q[3];
rz(-1.4596938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.725631) q[0];
sx q[0];
rz(-1.9209296) q[0];
sx q[0];
rz(-0.85335535) q[0];
rz(0.60494963) q[1];
sx q[1];
rz(-1.9974983) q[1];
sx q[1];
rz(-2.6430184) q[1];
rz(0.5729948) q[2];
sx q[2];
rz(-1.4481432) q[2];
sx q[2];
rz(-1.0021423) q[2];
rz(2.8145335) q[3];
sx q[3];
rz(-0.52486692) q[3];
sx q[3];
rz(-0.4384144) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

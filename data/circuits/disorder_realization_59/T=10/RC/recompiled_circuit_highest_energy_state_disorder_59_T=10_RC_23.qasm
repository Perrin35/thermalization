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
rz(-1.1745656) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6960775) q[0];
sx q[0];
rz(-0.68176523) q[0];
sx q[0];
rz(2.5755139) q[0];
rz(-pi) q[1];
rz(0.20148142) q[2];
sx q[2];
rz(-2.8116932) q[2];
sx q[2];
rz(3.0070674) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5331363) q[1];
sx q[1];
rz(-1.605833) q[1];
sx q[1];
rz(-0.94989454) q[1];
x q[2];
rz(1.6714736) q[3];
sx q[3];
rz(-1.0591456) q[3];
sx q[3];
rz(2.4477521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.55149469) q[2];
sx q[2];
rz(-1.8415201) q[2];
sx q[2];
rz(-0.91280118) q[2];
rz(-1.8707976) q[3];
sx q[3];
rz(-2.1015034) q[3];
sx q[3];
rz(2.2916268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1368644) q[0];
sx q[0];
rz(-2.5623463) q[0];
sx q[0];
rz(-0.42214033) q[0];
rz(2.7719851) q[1];
sx q[1];
rz(-2.044675) q[1];
sx q[1];
rz(2.2014528) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62806118) q[0];
sx q[0];
rz(-1.9002751) q[0];
sx q[0];
rz(0.0026436289) q[0];
rz(-pi) q[1];
rz(-1.3708417) q[2];
sx q[2];
rz(-2.5995289) q[2];
sx q[2];
rz(1.102239) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.247929) q[1];
sx q[1];
rz(-1.0808696) q[1];
sx q[1];
rz(0.72743261) q[1];
x q[2];
rz(0.89408447) q[3];
sx q[3];
rz(-1.7853123) q[3];
sx q[3];
rz(-1.0971951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9422354) q[2];
sx q[2];
rz(-0.67023674) q[2];
sx q[2];
rz(2.2352236) q[2];
rz(2.1622315) q[3];
sx q[3];
rz(-0.40458471) q[3];
sx q[3];
rz(-1.2649068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1449428) q[0];
sx q[0];
rz(-3.0746089) q[0];
sx q[0];
rz(1.8970733) q[0];
rz(0.38885802) q[1];
sx q[1];
rz(-1.2797979) q[1];
sx q[1];
rz(2.8188474) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9661315) q[0];
sx q[0];
rz(-1.4585988) q[0];
sx q[0];
rz(0.66288046) q[0];
x q[1];
rz(-1.0110485) q[2];
sx q[2];
rz(-1.8343535) q[2];
sx q[2];
rz(-1.4562869) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9936064) q[1];
sx q[1];
rz(-1.3393511) q[1];
sx q[1];
rz(0.79995059) q[1];
x q[2];
rz(3.0067919) q[3];
sx q[3];
rz(-2.7875337) q[3];
sx q[3];
rz(-1.1969138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2286223) q[2];
sx q[2];
rz(-1.3801489) q[2];
sx q[2];
rz(1.7886394) q[2];
rz(0.27255034) q[3];
sx q[3];
rz(-1.8504146) q[3];
sx q[3];
rz(-1.6779617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54101855) q[0];
sx q[0];
rz(-2.8193642) q[0];
sx q[0];
rz(1.2433276) q[0];
rz(0.86878949) q[1];
sx q[1];
rz(-2.181535) q[1];
sx q[1];
rz(-1.0702081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5015429) q[0];
sx q[0];
rz(-0.86210712) q[0];
sx q[0];
rz(-0.88592822) q[0];
rz(1.1203499) q[2];
sx q[2];
rz(-0.82847682) q[2];
sx q[2];
rz(-2.7591005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8319885) q[1];
sx q[1];
rz(-1.8924531) q[1];
sx q[1];
rz(3.020806) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7989705) q[3];
sx q[3];
rz(-1.3903769) q[3];
sx q[3];
rz(-0.29553771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59226817) q[2];
sx q[2];
rz(-2.3531239) q[2];
sx q[2];
rz(-2.2198524) q[2];
rz(0.24374572) q[3];
sx q[3];
rz(-1.4344401) q[3];
sx q[3];
rz(-0.29269472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9859966) q[0];
sx q[0];
rz(-1.3256185) q[0];
sx q[0];
rz(-2.6407114) q[0];
rz(-1.1174551) q[1];
sx q[1];
rz(-1.8084348) q[1];
sx q[1];
rz(-0.31455988) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18911823) q[0];
sx q[0];
rz(-0.50273147) q[0];
sx q[0];
rz(-0.76393868) q[0];
x q[1];
rz(-3.1371181) q[2];
sx q[2];
rz(-1.5273849) q[2];
sx q[2];
rz(-0.89148607) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.52175888) q[1];
sx q[1];
rz(-1.17621) q[1];
sx q[1];
rz(0.97418262) q[1];
rz(-pi) q[2];
rz(-0.1184095) q[3];
sx q[3];
rz(-1.1573912) q[3];
sx q[3];
rz(-1.6550696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.306281) q[2];
sx q[2];
rz(-1.1113144) q[2];
sx q[2];
rz(1.8446946) q[2];
rz(2.6563307) q[3];
sx q[3];
rz(-1.5819712) q[3];
sx q[3];
rz(-1.1135134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6473963) q[0];
sx q[0];
rz(-1.8504471) q[0];
sx q[0];
rz(3.0418292) q[0];
rz(-1.2707204) q[1];
sx q[1];
rz(-2.2427509) q[1];
sx q[1];
rz(-1.3894003) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6806599) q[0];
sx q[0];
rz(-1.4062728) q[0];
sx q[0];
rz(-2.8712832) q[0];
rz(-2.9676081) q[2];
sx q[2];
rz(-1.4672888) q[2];
sx q[2];
rz(-0.058920842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.278881) q[1];
sx q[1];
rz(-0.91966719) q[1];
sx q[1];
rz(-0.79336262) q[1];
x q[2];
rz(3.0148195) q[3];
sx q[3];
rz(-0.78713464) q[3];
sx q[3];
rz(1.8925557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6483267) q[2];
sx q[2];
rz(-2.6706225) q[2];
sx q[2];
rz(-0.89034447) q[2];
rz(1.1912311) q[3];
sx q[3];
rz(-1.3072562) q[3];
sx q[3];
rz(-2.2395535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.7118199) q[0];
sx q[0];
rz(-0.044782488) q[0];
sx q[0];
rz(3.1088767) q[0];
rz(-2.6383767) q[1];
sx q[1];
rz(-2.0423753) q[1];
sx q[1];
rz(3.018766) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5779992) q[0];
sx q[0];
rz(-1.8933354) q[0];
sx q[0];
rz(0.75046993) q[0];
rz(-pi) q[1];
rz(-1.2314447) q[2];
sx q[2];
rz(-1.095158) q[2];
sx q[2];
rz(2.1905157) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.87007755) q[1];
sx q[1];
rz(-0.48174274) q[1];
sx q[1];
rz(1.7863356) q[1];
rz(-1.7006247) q[3];
sx q[3];
rz(-1.2210088) q[3];
sx q[3];
rz(-1.2799124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7922625) q[2];
sx q[2];
rz(-0.96668875) q[2];
sx q[2];
rz(-1.1518504) q[2];
rz(-2.2564015) q[3];
sx q[3];
rz(-1.3693634) q[3];
sx q[3];
rz(1.7159897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18871466) q[0];
sx q[0];
rz(-2.2080244) q[0];
sx q[0];
rz(-0.92920148) q[0];
rz(0.62905351) q[1];
sx q[1];
rz(-1.355143) q[1];
sx q[1];
rz(-0.74310511) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0068664) q[0];
sx q[0];
rz(-1.4919625) q[0];
sx q[0];
rz(-1.1228485) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6531947) q[2];
sx q[2];
rz(-1.6553859) q[2];
sx q[2];
rz(-1.9458488) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3071909) q[1];
sx q[1];
rz(-1.1573004) q[1];
sx q[1];
rz(-1.4985227) q[1];
x q[2];
rz(-2.1984897) q[3];
sx q[3];
rz(-0.73929683) q[3];
sx q[3];
rz(3.0969248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5653845) q[2];
sx q[2];
rz(-0.78018633) q[2];
sx q[2];
rz(-0.72594491) q[2];
rz(-0.37096008) q[3];
sx q[3];
rz(-1.5981263) q[3];
sx q[3];
rz(2.7788739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2711656) q[0];
sx q[0];
rz(-2.3112264) q[0];
sx q[0];
rz(-0.56394947) q[0];
rz(1.14934) q[1];
sx q[1];
rz(-2.1795858) q[1];
sx q[1];
rz(-0.74251485) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19631736) q[0];
sx q[0];
rz(-2.4565897) q[0];
sx q[0];
rz(-2.3514868) q[0];
rz(-2.5385802) q[2];
sx q[2];
rz(-2.6508287) q[2];
sx q[2];
rz(3.1163851) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2143605) q[1];
sx q[1];
rz(-0.4402658) q[1];
sx q[1];
rz(0.10378598) q[1];
rz(-pi) q[2];
rz(0.85618505) q[3];
sx q[3];
rz(-1.1244698) q[3];
sx q[3];
rz(-2.3838338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5138862) q[2];
sx q[2];
rz(-2.0676925) q[2];
sx q[2];
rz(0.34454301) q[2];
rz(0.44376093) q[3];
sx q[3];
rz(-1.0756451) q[3];
sx q[3];
rz(1.3226604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76170707) q[0];
sx q[0];
rz(-0.2825309) q[0];
sx q[0];
rz(2.479582) q[0];
rz(0.33540353) q[1];
sx q[1];
rz(-1.6796203) q[1];
sx q[1];
rz(-2.2047156) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7279908) q[0];
sx q[0];
rz(-1.2112036) q[0];
sx q[0];
rz(2.2171564) q[0];
rz(-1.4858021) q[2];
sx q[2];
rz(-2.8704145) q[2];
sx q[2];
rz(-0.32395872) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7422694) q[1];
sx q[1];
rz(-1.8571257) q[1];
sx q[1];
rz(2.7469357) q[1];
rz(1.7702661) q[3];
sx q[3];
rz(-1.7430309) q[3];
sx q[3];
rz(-1.0125835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.718049) q[2];
sx q[2];
rz(-2.9073538) q[2];
sx q[2];
rz(2.6424109) q[2];
rz(-1.6792123) q[3];
sx q[3];
rz(-1.9951818) q[3];
sx q[3];
rz(-1.4596938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.725631) q[0];
sx q[0];
rz(-1.9209296) q[0];
sx q[0];
rz(-0.85335535) q[0];
rz(2.536643) q[1];
sx q[1];
rz(-1.1440944) q[1];
sx q[1];
rz(0.49857421) q[1];
rz(-2.9180183) q[2];
sx q[2];
rz(-0.58453544) q[2];
sx q[2];
rz(0.75605308) q[2];
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

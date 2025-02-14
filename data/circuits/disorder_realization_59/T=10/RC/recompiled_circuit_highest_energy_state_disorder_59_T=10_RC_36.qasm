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
rz(1.356025) q[1];
sx q[1];
rz(7.4577509) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0102219) q[0];
sx q[0];
rz(-1.0099942) q[0];
sx q[0];
rz(-1.9813374) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9401112) q[2];
sx q[2];
rz(-0.32989943) q[2];
sx q[2];
rz(3.0070674) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2043031) q[1];
sx q[1];
rz(-2.1912592) q[1];
sx q[1];
rz(3.0985249) q[1];
x q[2];
rz(1.470119) q[3];
sx q[3];
rz(-1.0591456) q[3];
sx q[3];
rz(0.6938405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.55149469) q[2];
sx q[2];
rz(-1.3000725) q[2];
sx q[2];
rz(2.2287915) q[2];
rz(1.8707976) q[3];
sx q[3];
rz(-1.0400892) q[3];
sx q[3];
rz(-0.84996581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0047282334) q[0];
sx q[0];
rz(-2.5623463) q[0];
sx q[0];
rz(-2.7194523) q[0];
rz(0.36960754) q[1];
sx q[1];
rz(-2.044675) q[1];
sx q[1];
rz(0.94013989) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5053609) q[0];
sx q[0];
rz(-2.8121037) q[0];
sx q[0];
rz(-1.5785274) q[0];
rz(2.1040078) q[2];
sx q[2];
rz(-1.4681446) q[2];
sx q[2];
rz(0.29666049) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.2804037) q[1];
sx q[1];
rz(-2.197816) q[1];
sx q[1];
rz(2.1908733) q[1];
rz(-pi) q[2];
rz(1.2359911) q[3];
sx q[3];
rz(-2.4368317) q[3];
sx q[3];
rz(-2.4089264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1993572) q[2];
sx q[2];
rz(-0.67023674) q[2];
sx q[2];
rz(-0.90636903) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9966499) q[0];
sx q[0];
rz(-0.066983797) q[0];
sx q[0];
rz(-1.8970733) q[0];
rz(2.7527346) q[1];
sx q[1];
rz(-1.2797979) q[1];
sx q[1];
rz(0.32274524) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8334482) q[0];
sx q[0];
rz(-2.2287773) q[0];
sx q[0];
rz(1.7127772) q[0];
x q[1];
rz(-1.1006196) q[2];
sx q[2];
rz(-0.61264804) q[2];
sx q[2];
rz(-2.8620811) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6420329) q[1];
sx q[1];
rz(-2.316058) q[1];
sx q[1];
rz(2.8241629) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1348008) q[3];
sx q[3];
rz(-2.7875337) q[3];
sx q[3];
rz(1.9446789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9129703) q[2];
sx q[2];
rz(-1.7614438) q[2];
sx q[2];
rz(-1.3529533) q[2];
rz(0.27255034) q[3];
sx q[3];
rz(-1.8504146) q[3];
sx q[3];
rz(1.4636309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54101855) q[0];
sx q[0];
rz(-0.32222846) q[0];
sx q[0];
rz(1.2433276) q[0];
rz(-2.2728032) q[1];
sx q[1];
rz(-2.181535) q[1];
sx q[1];
rz(-1.0702081) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6400498) q[0];
sx q[0];
rz(-0.86210712) q[0];
sx q[0];
rz(2.2556644) q[0];
x q[1];
rz(-2.0212428) q[2];
sx q[2];
rz(-0.82847682) q[2];
sx q[2];
rz(0.38249215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1985716) q[1];
sx q[1];
rz(-0.34284886) q[1];
sx q[1];
rz(-1.2238316) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49740187) q[3];
sx q[3];
rz(-0.38555749) q[3];
sx q[3];
rz(-1.7413643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5493245) q[2];
sx q[2];
rz(-2.3531239) q[2];
sx q[2];
rz(-0.92174021) q[2];
rz(2.8978469) q[3];
sx q[3];
rz(-1.7071525) q[3];
sx q[3];
rz(2.8488979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9859966) q[0];
sx q[0];
rz(-1.8159741) q[0];
sx q[0];
rz(-0.50088125) q[0];
rz(1.1174551) q[1];
sx q[1];
rz(-1.3331579) q[1];
sx q[1];
rz(2.8270328) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4582493) q[0];
sx q[0];
rz(-1.2309845) q[0];
sx q[0];
rz(0.37796867) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4681513) q[2];
sx q[2];
rz(-3.0979514) q[2];
sx q[2];
rz(0.99422821) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52175888) q[1];
sx q[1];
rz(-1.9653826) q[1];
sx q[1];
rz(-2.16741) q[1];
rz(-pi) q[2];
rz(-1.9867927) q[3];
sx q[3];
rz(-1.6791897) q[3];
sx q[3];
rz(3.1050753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8353117) q[2];
sx q[2];
rz(-2.0302782) q[2];
sx q[2];
rz(1.8446946) q[2];
rz(-0.48526192) q[3];
sx q[3];
rz(-1.5819712) q[3];
sx q[3];
rz(-1.1135134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49419633) q[0];
sx q[0];
rz(-1.2911456) q[0];
sx q[0];
rz(0.099763481) q[0];
rz(-1.2707204) q[1];
sx q[1];
rz(-0.89884177) q[1];
sx q[1];
rz(-1.7521923) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064512028) q[0];
sx q[0];
rz(-1.837366) q[0];
sx q[0];
rz(-1.4001911) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9676081) q[2];
sx q[2];
rz(-1.4672888) q[2];
sx q[2];
rz(0.058920842) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8956611) q[1];
sx q[1];
rz(-0.97890039) q[1];
sx q[1];
rz(-2.3227957) q[1];
rz(-pi) q[2];
rz(-1.6969958) q[3];
sx q[3];
rz(-0.79168237) q[3];
sx q[3];
rz(1.070397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6483267) q[2];
sx q[2];
rz(-2.6706225) q[2];
sx q[2];
rz(-2.2512482) q[2];
rz(-1.1912311) q[3];
sx q[3];
rz(-1.8343364) q[3];
sx q[3];
rz(-2.2395535) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4297727) q[0];
sx q[0];
rz(-0.044782488) q[0];
sx q[0];
rz(0.032715948) q[0];
rz(2.6383767) q[1];
sx q[1];
rz(-2.0423753) q[1];
sx q[1];
rz(0.12282664) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5635934) q[0];
sx q[0];
rz(-1.2482572) q[0];
sx q[0];
rz(-0.75046993) q[0];
rz(0.57374222) q[2];
sx q[2];
rz(-2.5649568) q[2];
sx q[2];
rz(1.5337616) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.029307) q[1];
sx q[1];
rz(-1.1011135) q[1];
sx q[1];
rz(0.11135688) q[1];
rz(1.7006247) q[3];
sx q[3];
rz(-1.9205838) q[3];
sx q[3];
rz(-1.2799124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7922625) q[2];
sx q[2];
rz(-2.1749039) q[2];
sx q[2];
rz(1.1518504) q[2];
rz(-2.2564015) q[3];
sx q[3];
rz(-1.3693634) q[3];
sx q[3];
rz(-1.4256029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(2.952878) q[0];
sx q[0];
rz(-2.2080244) q[0];
sx q[0];
rz(-0.92920148) q[0];
rz(-2.5125391) q[1];
sx q[1];
rz(-1.355143) q[1];
sx q[1];
rz(-0.74310511) q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(-1.4862068) q[2];
sx q[2];
rz(-1.1957439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3071909) q[1];
sx q[1];
rz(-1.9842923) q[1];
sx q[1];
rz(1.6430699) q[1];
x q[2];
rz(0.93502126) q[3];
sx q[3];
rz(-1.9776157) q[3];
sx q[3];
rz(-2.0182899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5653845) q[2];
sx q[2];
rz(-2.3614063) q[2];
sx q[2];
rz(-2.4156477) q[2];
rz(0.37096008) q[3];
sx q[3];
rz(-1.5434664) q[3];
sx q[3];
rz(-0.36271873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87042701) q[0];
sx q[0];
rz(-2.3112264) q[0];
sx q[0];
rz(-0.56394947) q[0];
rz(1.9922527) q[1];
sx q[1];
rz(-2.1795858) q[1];
sx q[1];
rz(0.74251485) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71096984) q[0];
sx q[0];
rz(-2.0369684) q[0];
sx q[0];
rz(2.6197893) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41461035) q[2];
sx q[2];
rz(-1.8413723) q[2];
sx q[2];
rz(2.0913578) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0997252) q[1];
sx q[1];
rz(-2.0085287) q[1];
sx q[1];
rz(1.619564) q[1];
rz(-pi) q[2];
rz(2.2015195) q[3];
sx q[3];
rz(-2.3204062) q[3];
sx q[3];
rz(1.8668777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5138862) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76170707) q[0];
sx q[0];
rz(-0.2825309) q[0];
sx q[0];
rz(0.66201061) q[0];
rz(2.8061891) q[1];
sx q[1];
rz(-1.4619724) q[1];
sx q[1];
rz(-2.2047156) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59360817) q[0];
sx q[0];
rz(-0.72692211) q[0];
sx q[0];
rz(1.0127823) q[0];
rz(3.1179948) q[2];
sx q[2];
rz(-1.3006217) q[2];
sx q[2];
rz(0.23575704) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7422694) q[1];
sx q[1];
rz(-1.8571257) q[1];
sx q[1];
rz(-2.7469357) q[1];
x q[2];
rz(-2.2912737) q[3];
sx q[3];
rz(-0.26278824) q[3];
sx q[3];
rz(2.9966054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4235437) q[2];
sx q[2];
rz(-0.23423883) q[2];
sx q[2];
rz(0.49918175) q[2];
rz(1.4623803) q[3];
sx q[3];
rz(-1.9951818) q[3];
sx q[3];
rz(1.6818989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
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
rz(-1.7164604) q[2];
sx q[2];
rz(-1.002641) q[2];
sx q[2];
rz(-2.6517131) q[2];
rz(-1.3868757) q[3];
sx q[3];
rz(-1.0763604) q[3];
sx q[3];
rz(3.0767783) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

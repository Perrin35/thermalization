OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.225086) q[0];
sx q[0];
rz(-0.14296159) q[0];
sx q[0];
rz(11.077865) q[0];
rz(-2.0090964) q[1];
sx q[1];
rz(-0.43192616) q[1];
sx q[1];
rz(-2.5052524) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36447517) q[0];
sx q[0];
rz(-0.72491966) q[0];
sx q[0];
rz(-2.1687647) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9540399) q[2];
sx q[2];
rz(-1.2072717) q[2];
sx q[2];
rz(0.55628796) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.591037) q[1];
sx q[1];
rz(-1.6948957) q[1];
sx q[1];
rz(-1.9684588) q[1];
x q[2];
rz(-0.79145517) q[3];
sx q[3];
rz(-0.88578445) q[3];
sx q[3];
rz(2.5572957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1119614) q[2];
sx q[2];
rz(-1.1250857) q[2];
sx q[2];
rz(-2.4911528) q[2];
rz(1.8247617) q[3];
sx q[3];
rz(-1.7951169) q[3];
sx q[3];
rz(3.0397547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0874852) q[0];
sx q[0];
rz(-0.58687812) q[0];
sx q[0];
rz(2.6869539) q[0];
rz(-3.1184323) q[1];
sx q[1];
rz(-1.4440447) q[1];
sx q[1];
rz(-0.32618162) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5989203) q[0];
sx q[0];
rz(-0.70264757) q[0];
sx q[0];
rz(-0.62851466) q[0];
rz(-0.97450476) q[2];
sx q[2];
rz(-2.6428004) q[2];
sx q[2];
rz(1.5128262) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7506894) q[1];
sx q[1];
rz(-1.1613966) q[1];
sx q[1];
rz(-0.55623325) q[1];
x q[2];
rz(-3.0241859) q[3];
sx q[3];
rz(-2.2281584) q[3];
sx q[3];
rz(2.9598992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.11639103) q[2];
sx q[2];
rz(-1.9390743) q[2];
sx q[2];
rz(-0.95428673) q[2];
rz(1.8918234) q[3];
sx q[3];
rz(-2.818675) q[3];
sx q[3];
rz(3.1402816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2713476) q[0];
sx q[0];
rz(-1.1529237) q[0];
sx q[0];
rz(-1.1767607) q[0];
rz(0.36508834) q[1];
sx q[1];
rz(-1.3026214) q[1];
sx q[1];
rz(1.5286068) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9823526) q[0];
sx q[0];
rz(-1.591658) q[0];
sx q[0];
rz(0.013289159) q[0];
rz(-pi) q[1];
rz(0.079054376) q[2];
sx q[2];
rz(-2.0131265) q[2];
sx q[2];
rz(2.4238264) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4499719) q[1];
sx q[1];
rz(-0.1547389) q[1];
sx q[1];
rz(1.7697033) q[1];
rz(-0.068644329) q[3];
sx q[3];
rz(-1.6283096) q[3];
sx q[3];
rz(-0.36348331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24885808) q[2];
sx q[2];
rz(-2.5881519) q[2];
sx q[2];
rz(-1.152732) q[2];
rz(-1.24498) q[3];
sx q[3];
rz(-0.98074073) q[3];
sx q[3];
rz(0.28321701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98274851) q[0];
sx q[0];
rz(-pi/12) q[0];
sx q[0];
rz(-0.68956476) q[0];
rz(-3.116563) q[1];
sx q[1];
rz(-1.5182779) q[1];
sx q[1];
rz(-1.4208581) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8129028) q[0];
sx q[0];
rz(-2.5089009) q[0];
sx q[0];
rz(-2.9660874) q[0];
x q[1];
rz(2.000932) q[2];
sx q[2];
rz(-1.7342415) q[2];
sx q[2];
rz(-2.8317833) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5031227) q[1];
sx q[1];
rz(-1.5984189) q[1];
sx q[1];
rz(0.50325985) q[1];
rz(2.8042364) q[3];
sx q[3];
rz(-1.8280262) q[3];
sx q[3];
rz(1.6599865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.932852) q[2];
sx q[2];
rz(-0.17540652) q[2];
sx q[2];
rz(-1.9746732) q[2];
rz(-2.6973727) q[3];
sx q[3];
rz(-1.9053713) q[3];
sx q[3];
rz(-2.8319632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2883478) q[0];
sx q[0];
rz(-1.3695559) q[0];
sx q[0];
rz(-2.7254768) q[0];
rz(0.5101282) q[1];
sx q[1];
rz(-2.3704539) q[1];
sx q[1];
rz(-0.77521926) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8213971) q[0];
sx q[0];
rz(-2.1761083) q[0];
sx q[0];
rz(1.3277675) q[0];
rz(-1.3890319) q[2];
sx q[2];
rz(-1.9700288) q[2];
sx q[2];
rz(2.4295074) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0803804) q[1];
sx q[1];
rz(-2.9728372) q[1];
sx q[1];
rz(-2.3024998) q[1];
rz(2.6121077) q[3];
sx q[3];
rz(-0.66505243) q[3];
sx q[3];
rz(2.3163026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9267209) q[2];
sx q[2];
rz(-0.99433172) q[2];
sx q[2];
rz(1.035824) q[2];
rz(2.9521247) q[3];
sx q[3];
rz(-2.6697956) q[3];
sx q[3];
rz(-2.6389879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.254461) q[0];
sx q[0];
rz(-3.0939565) q[0];
sx q[0];
rz(2.7857842) q[0];
rz(1.4954781) q[1];
sx q[1];
rz(-2.1864086) q[1];
sx q[1];
rz(1.1411508) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3095207) q[0];
sx q[0];
rz(-2.4395925) q[0];
sx q[0];
rz(-1.2418447) q[0];
rz(1.5836493) q[2];
sx q[2];
rz(-0.61159482) q[2];
sx q[2];
rz(-2.0919959) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6427152) q[1];
sx q[1];
rz(-1.4659473) q[1];
sx q[1];
rz(3.0003606) q[1];
rz(-pi) q[2];
rz(-3.1294501) q[3];
sx q[3];
rz(-2.0046765) q[3];
sx q[3];
rz(-1.7858684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.63153875) q[2];
sx q[2];
rz(-0.71439356) q[2];
sx q[2];
rz(1.7047403) q[2];
rz(-1.0531744) q[3];
sx q[3];
rz(-1.6097693) q[3];
sx q[3];
rz(0.50301445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.961504) q[0];
sx q[0];
rz(-0.29619521) q[0];
sx q[0];
rz(-0.12251138) q[0];
rz(-0.65613121) q[1];
sx q[1];
rz(-1.8729112) q[1];
sx q[1];
rz(1.8001385) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3229834) q[0];
sx q[0];
rz(-3.083813) q[0];
sx q[0];
rz(0.66814439) q[0];
x q[1];
rz(-1.6788922) q[2];
sx q[2];
rz(-1.2891532) q[2];
sx q[2];
rz(0.070527129) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.16640284) q[1];
sx q[1];
rz(-0.67968785) q[1];
sx q[1];
rz(-1.1389334) q[1];
x q[2];
rz(2.6977067) q[3];
sx q[3];
rz(-1.4680913) q[3];
sx q[3];
rz(-3.1004124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.573367) q[2];
sx q[2];
rz(-1.2218916) q[2];
sx q[2];
rz(1.4637671) q[2];
rz(-0.73417869) q[3];
sx q[3];
rz(-0.28407431) q[3];
sx q[3];
rz(-1.8370139) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7952591) q[0];
sx q[0];
rz(-1.3672375) q[0];
sx q[0];
rz(-2.6829868) q[0];
rz(2.1268225) q[1];
sx q[1];
rz(-1.1943123) q[1];
sx q[1];
rz(2.0789304) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.748007) q[0];
sx q[0];
rz(-0.17608041) q[0];
sx q[0];
rz(2.8066471) q[0];
x q[1];
rz(0.52824654) q[2];
sx q[2];
rz(-2.2604803) q[2];
sx q[2];
rz(1.6843585) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.96982376) q[1];
sx q[1];
rz(-0.14482982) q[1];
sx q[1];
rz(-2.0098196) q[1];
rz(-1.0414391) q[3];
sx q[3];
rz(-2.5564407) q[3];
sx q[3];
rz(-0.063466788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8899272) q[2];
sx q[2];
rz(-1.3883611) q[2];
sx q[2];
rz(-1.3512705) q[2];
rz(-1.2190367) q[3];
sx q[3];
rz(-0.42870298) q[3];
sx q[3];
rz(-0.8333227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079000533) q[0];
sx q[0];
rz(-1.7575678) q[0];
sx q[0];
rz(0.90743995) q[0];
rz(-2.9145248) q[1];
sx q[1];
rz(-2.0957969) q[1];
sx q[1];
rz(-2.2423832) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19455043) q[0];
sx q[0];
rz(-0.83861807) q[0];
sx q[0];
rz(-2.28449) q[0];
rz(2.3311798) q[2];
sx q[2];
rz(-0.11356662) q[2];
sx q[2];
rz(-1.0134987) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6688706) q[1];
sx q[1];
rz(-0.54520118) q[1];
sx q[1];
rz(1.2073869) q[1];
rz(-pi) q[2];
rz(1.5158922) q[3];
sx q[3];
rz(-2.3927352) q[3];
sx q[3];
rz(-2.1532173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1850618) q[2];
sx q[2];
rz(-1.1659634) q[2];
sx q[2];
rz(-0.60066191) q[2];
rz(-1.0968084) q[3];
sx q[3];
rz(-0.59022248) q[3];
sx q[3];
rz(-2.2685331) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3908865) q[0];
sx q[0];
rz(-1.0989256) q[0];
sx q[0];
rz(-2.1571958) q[0];
rz(-0.4920494) q[1];
sx q[1];
rz(-1.4130519) q[1];
sx q[1];
rz(1.1955998) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0427143) q[0];
sx q[0];
rz(-0.98473583) q[0];
sx q[0];
rz(-1.9978961) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2480505) q[2];
sx q[2];
rz(-1.9697426) q[2];
sx q[2];
rz(-0.91541399) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8859133) q[1];
sx q[1];
rz(-1.9530086) q[1];
sx q[1];
rz(-3.0542733) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3687184) q[3];
sx q[3];
rz(-1.2307271) q[3];
sx q[3];
rz(-3.056385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.055858) q[2];
sx q[2];
rz(-1.6207638) q[2];
sx q[2];
rz(-2.1585507) q[2];
rz(0.25162697) q[3];
sx q[3];
rz(-2.7139137) q[3];
sx q[3];
rz(1.0035286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6428103) q[0];
sx q[0];
rz(-2.1693873) q[0];
sx q[0];
rz(0.65709773) q[0];
rz(1.2596754) q[1];
sx q[1];
rz(-1.4114264) q[1];
sx q[1];
rz(-2.2687601) q[1];
rz(0.85763422) q[2];
sx q[2];
rz(-1.4076172) q[2];
sx q[2];
rz(2.313688) q[2];
rz(-1.891741) q[3];
sx q[3];
rz(-1.2529916) q[3];
sx q[3];
rz(2.8734251) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

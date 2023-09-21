OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.94937593) q[0];
sx q[0];
rz(-1.047171) q[0];
sx q[0];
rz(0.068724364) q[0];
rz(-1.3955431) q[1];
sx q[1];
rz(-1.5323324) q[1];
sx q[1];
rz(1.2083763) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8954593) q[0];
sx q[0];
rz(-1.5243422) q[0];
sx q[0];
rz(1.4746656) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5586833) q[2];
sx q[2];
rz(-1.4450057) q[2];
sx q[2];
rz(-2.8319401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.87885107) q[1];
sx q[1];
rz(-0.95572119) q[1];
sx q[1];
rz(1.4351074) q[1];
rz(-2.6192059) q[3];
sx q[3];
rz(-1.1484227) q[3];
sx q[3];
rz(1.9463604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3258813) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(1.4665843) q[2];
rz(2.4438434) q[3];
sx q[3];
rz(-1.1013228) q[3];
sx q[3];
rz(-2.3944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0242457) q[0];
sx q[0];
rz(-1.1276561) q[0];
sx q[0];
rz(-1.9673989) q[0];
rz(-2.970447) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(-2.8443964) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3127808) q[0];
sx q[0];
rz(-0.096185616) q[0];
sx q[0];
rz(-2.0975153) q[0];
rz(-pi) q[1];
rz(-2.4291123) q[2];
sx q[2];
rz(-1.9409632) q[2];
sx q[2];
rz(-0.21288255) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.17644037) q[1];
sx q[1];
rz(-2.0808176) q[1];
sx q[1];
rz(-2.7494207) q[1];
rz(-pi) q[2];
rz(2.3791749) q[3];
sx q[3];
rz(-1.8652417) q[3];
sx q[3];
rz(0.03014119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.16195665) q[2];
sx q[2];
rz(-1.9571783) q[2];
sx q[2];
rz(-0.61398181) q[2];
rz(-2.2654514) q[3];
sx q[3];
rz(-0.66771475) q[3];
sx q[3];
rz(2.5045625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0959594) q[0];
sx q[0];
rz(-1.8563844) q[0];
sx q[0];
rz(-2.8339548) q[0];
rz(0.74854198) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(0.83980733) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15159431) q[0];
sx q[0];
rz(-2.6322106) q[0];
sx q[0];
rz(0.61181061) q[0];
x q[1];
rz(0.26250458) q[2];
sx q[2];
rz(-0.35048198) q[2];
sx q[2];
rz(1.0002491) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0452022) q[1];
sx q[1];
rz(-2.5844816) q[1];
sx q[1];
rz(2.5379009) q[1];
rz(-pi) q[2];
rz(-0.17383667) q[3];
sx q[3];
rz(-0.94508119) q[3];
sx q[3];
rz(-3.0415149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86747375) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(-1.8236558) q[2];
rz(1.9258202) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(1.5095476) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76386219) q[0];
sx q[0];
rz(-1.3803991) q[0];
sx q[0];
rz(0.44556251) q[0];
rz(-0.62082779) q[1];
sx q[1];
rz(-1.8955684) q[1];
sx q[1];
rz(-0.96558085) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5774028) q[0];
sx q[0];
rz(-1.8468542) q[0];
sx q[0];
rz(0.77703707) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.031164073) q[2];
sx q[2];
rz(-0.84261299) q[2];
sx q[2];
rz(-0.63809168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.671215) q[1];
sx q[1];
rz(-0.2556076) q[1];
sx q[1];
rz(-1.6971991) q[1];
rz(-pi) q[2];
rz(-3.0442763) q[3];
sx q[3];
rz(-1.4574377) q[3];
sx q[3];
rz(1.3219584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.821637) q[2];
sx q[2];
rz(-2.379202) q[2];
sx q[2];
rz(2.2122673) q[2];
rz(-2.4980513) q[3];
sx q[3];
rz(-2.1054335) q[3];
sx q[3];
rz(-2.1381366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4716361) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(-1.2840282) q[0];
rz(-0.28981003) q[1];
sx q[1];
rz(-2.402014) q[1];
sx q[1];
rz(1.0481542) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7832344) q[0];
sx q[0];
rz(-0.90792197) q[0];
sx q[0];
rz(2.5213084) q[0];
rz(-pi) q[1];
rz(0.83277793) q[2];
sx q[2];
rz(-1.1156429) q[2];
sx q[2];
rz(-0.27311329) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3874665) q[1];
sx q[1];
rz(-1.4782463) q[1];
sx q[1];
rz(2.8184163) q[1];
rz(-pi) q[2];
rz(-2.6606779) q[3];
sx q[3];
rz(-2.1926583) q[3];
sx q[3];
rz(-2.8180168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8706878) q[2];
sx q[2];
rz(-1.015816) q[2];
sx q[2];
rz(-2.2407545) q[2];
rz(2.0488996) q[3];
sx q[3];
rz(-2.1381502) q[3];
sx q[3];
rz(1.9074915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822405) q[0];
sx q[0];
rz(-1.2077967) q[0];
sx q[0];
rz(0.59610468) q[0];
rz(1.4959363) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(-1.2449107) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0743474) q[0];
sx q[0];
rz(-2.0739569) q[0];
sx q[0];
rz(2.7748681) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5858324) q[2];
sx q[2];
rz(-0.55609497) q[2];
sx q[2];
rz(0.13242002) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4391172) q[1];
sx q[1];
rz(-0.87279746) q[1];
sx q[1];
rz(2.4707787) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8935029) q[3];
sx q[3];
rz(-1.7956942) q[3];
sx q[3];
rz(3.0450862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.231455) q[2];
sx q[2];
rz(-0.69880501) q[2];
sx q[2];
rz(2.9166252) q[2];
rz(-3.0531626) q[3];
sx q[3];
rz(-1.6902573) q[3];
sx q[3];
rz(0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46463075) q[0];
sx q[0];
rz(-1.9043652) q[0];
sx q[0];
rz(0.27994573) q[0];
rz(1.6784558) q[1];
sx q[1];
rz(-1.8755553) q[1];
sx q[1];
rz(2.8889012) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42851105) q[0];
sx q[0];
rz(-1.5230852) q[0];
sx q[0];
rz(-1.8639355) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93637034) q[2];
sx q[2];
rz(-1.5555824) q[2];
sx q[2];
rz(2.0903367) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.686324) q[1];
sx q[1];
rz(-2.2043921) q[1];
sx q[1];
rz(-0.92647657) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2884376) q[3];
sx q[3];
rz(-0.95082885) q[3];
sx q[3];
rz(-0.48228797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8213886) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(1.5318711) q[2];
rz(1.1931217) q[3];
sx q[3];
rz(-0.90819287) q[3];
sx q[3];
rz(2.0549205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10483345) q[0];
sx q[0];
rz(-1.3830673) q[0];
sx q[0];
rz(-1.6554792) q[0];
rz(-0.2688109) q[1];
sx q[1];
rz(-1.1299645) q[1];
sx q[1];
rz(-0.2789467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3725961) q[0];
sx q[0];
rz(-2.0645752) q[0];
sx q[0];
rz(0.84817024) q[0];
x q[1];
rz(-2.0729162) q[2];
sx q[2];
rz(-1.7218105) q[2];
sx q[2];
rz(-2.6627024) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8262417) q[1];
sx q[1];
rz(-2.3126174) q[1];
sx q[1];
rz(2.6226603) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99854462) q[3];
sx q[3];
rz(-1.4456985) q[3];
sx q[3];
rz(2.3541114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.802861) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(1.5489244) q[2];
rz(-1.0507978) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(1.9416434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779697) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(-1.2532225) q[0];
rz(2.5190917) q[1];
sx q[1];
rz(-1.4600735) q[1];
sx q[1];
rz(-1.1463096) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63908731) q[0];
sx q[0];
rz(-2.3168132) q[0];
sx q[0];
rz(-1.6665002) q[0];
rz(2.4856604) q[2];
sx q[2];
rz(-2.1428875) q[2];
sx q[2];
rz(2.7149372) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4790926) q[1];
sx q[1];
rz(-0.88516419) q[1];
sx q[1];
rz(0.094866026) q[1];
rz(-pi) q[2];
rz(0.54543145) q[3];
sx q[3];
rz(-2.430393) q[3];
sx q[3];
rz(1.4382854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15381947) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(-1.0158319) q[2];
rz(-0.9097957) q[3];
sx q[3];
rz(-2.4715021) q[3];
sx q[3];
rz(-0.36809665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0125473) q[0];
sx q[0];
rz(-0.61360252) q[0];
sx q[0];
rz(0.01424271) q[0];
rz(-2.2968538) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(1.7600118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2505014) q[0];
sx q[0];
rz(-1.7133461) q[0];
sx q[0];
rz(0.010682627) q[0];
x q[1];
rz(-1.472677) q[2];
sx q[2];
rz(-2.5604904) q[2];
sx q[2];
rz(-0.89400089) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4631008) q[1];
sx q[1];
rz(-1.4271724) q[1];
sx q[1];
rz(-1.0113869) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36079447) q[3];
sx q[3];
rz(-0.56305712) q[3];
sx q[3];
rz(0.35239708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0649197) q[2];
sx q[2];
rz(-1.939247) q[2];
sx q[2];
rz(-0.98999611) q[2];
rz(2.2475217) q[3];
sx q[3];
rz(-0.48864135) q[3];
sx q[3];
rz(-0.22542424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0513231) q[0];
sx q[0];
rz(-0.6482424) q[0];
sx q[0];
rz(2.0363664) q[0];
rz(1.3399711) q[1];
sx q[1];
rz(-0.62146386) q[1];
sx q[1];
rz(0.38846831) q[1];
rz(-2.4253035) q[2];
sx q[2];
rz(-1.204797) q[2];
sx q[2];
rz(2.7439678) q[2];
rz(2.7064825) q[3];
sx q[3];
rz(-0.995244) q[3];
sx q[3];
rz(1.4959195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

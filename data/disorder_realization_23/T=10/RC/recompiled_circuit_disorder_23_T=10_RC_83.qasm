OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(7.6927778) q[0];
sx q[0];
rz(11.132244) q[0];
rz(-2.5073476) q[1];
sx q[1];
rz(-0.60159644) q[1];
sx q[1];
rz(-2.7231725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.943676) q[0];
sx q[0];
rz(-2.0352053) q[0];
sx q[0];
rz(2.8465413) q[0];
x q[1];
rz(1.7069874) q[2];
sx q[2];
rz(-1.4601267) q[2];
sx q[2];
rz(1.9905123) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99416713) q[1];
sx q[1];
rz(-2.6174183) q[1];
sx q[1];
rz(-1.0679354) q[1];
rz(-pi) q[2];
rz(-2.0099785) q[3];
sx q[3];
rz(-1.3703128) q[3];
sx q[3];
rz(0.85103121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.819954) q[2];
sx q[2];
rz(-1.7724089) q[2];
sx q[2];
rz(-0.83797541) q[2];
rz(2.6485802) q[3];
sx q[3];
rz(-0.27291441) q[3];
sx q[3];
rz(3.0626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.74719602) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(-1.863742) q[0];
rz(2.9648119) q[1];
sx q[1];
rz(-1.3143833) q[1];
sx q[1];
rz(-2.7094254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1457739) q[0];
sx q[0];
rz(-0.60129014) q[0];
sx q[0];
rz(-2.6390618) q[0];
rz(-pi) q[1];
rz(2.7414397) q[2];
sx q[2];
rz(-1.1740985) q[2];
sx q[2];
rz(0.19043365) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4914815) q[1];
sx q[1];
rz(-1.1307798) q[1];
sx q[1];
rz(0.12932175) q[1];
rz(-pi) q[2];
rz(1.0057955) q[3];
sx q[3];
rz(-1.5917935) q[3];
sx q[3];
rz(-2.3137623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.54923487) q[2];
sx q[2];
rz(-1.2640415) q[2];
sx q[2];
rz(-0.48669997) q[2];
rz(-1.3782079) q[3];
sx q[3];
rz(-1.8816201) q[3];
sx q[3];
rz(-0.53282213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.10087) q[0];
sx q[0];
rz(-0.7464872) q[0];
sx q[0];
rz(-2.7242463) q[0];
rz(1.6529282) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(-2.6352077) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8538118) q[0];
sx q[0];
rz(-1.3940485) q[0];
sx q[0];
rz(-1.9275097) q[0];
rz(-pi) q[1];
rz(-1.8110397) q[2];
sx q[2];
rz(-0.85075399) q[2];
sx q[2];
rz(-2.4400997) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0326929) q[1];
sx q[1];
rz(-1.1316205) q[1];
sx q[1];
rz(-1.2296618) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18249986) q[3];
sx q[3];
rz(-1.9347408) q[3];
sx q[3];
rz(2.5516627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.69904077) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(-2.55012) q[2];
rz(2.5555723) q[3];
sx q[3];
rz(-1.2089217) q[3];
sx q[3];
rz(-1.7104141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9716924) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(-2.3024094) q[0];
rz(3.1160141) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(-1.5485839) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73491053) q[0];
sx q[0];
rz(-1.6618177) q[0];
sx q[0];
rz(-2.0216366) q[0];
x q[1];
rz(0.84393878) q[2];
sx q[2];
rz(-2.2288449) q[2];
sx q[2];
rz(-0.3447926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81739391) q[1];
sx q[1];
rz(-2.3014268) q[1];
sx q[1];
rz(2.6170931) q[1];
rz(-1.8539092) q[3];
sx q[3];
rz(-1.9861756) q[3];
sx q[3];
rz(1.1113885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8362391) q[2];
sx q[2];
rz(-2.2576015) q[2];
sx q[2];
rz(-3.0419066) q[2];
rz(2.1827407) q[3];
sx q[3];
rz(-1.3189664) q[3];
sx q[3];
rz(-1.8166186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754159) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(-2.8856522) q[0];
rz(0.4610962) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(2.3815313) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5998659) q[0];
sx q[0];
rz(-1.4540298) q[0];
sx q[0];
rz(-1.6745425) q[0];
rz(-3.0779482) q[2];
sx q[2];
rz(-1.1568937) q[2];
sx q[2];
rz(0.49583437) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.45436817) q[1];
sx q[1];
rz(-1.4441274) q[1];
sx q[1];
rz(-1.6997937) q[1];
rz(-pi) q[2];
rz(2.2171668) q[3];
sx q[3];
rz(-1.8674208) q[3];
sx q[3];
rz(-2.7588206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.57050675) q[2];
sx q[2];
rz(-1.3723624) q[2];
sx q[2];
rz(2.467353) q[2];
rz(-0.21480602) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(-3.1242483) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6102585) q[0];
sx q[0];
rz(-1.4704309) q[0];
sx q[0];
rz(1.9624788) q[0];
rz(-0.20482652) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(1.0669605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63283352) q[0];
sx q[0];
rz(-1.5501223) q[0];
sx q[0];
rz(-1.5563006) q[0];
x q[1];
rz(-1.1322137) q[2];
sx q[2];
rz(-1.9595993) q[2];
sx q[2];
rz(2.4593381) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8300007) q[1];
sx q[1];
rz(-2.316906) q[1];
sx q[1];
rz(3.0768865) q[1];
rz(2.9818929) q[3];
sx q[3];
rz(-2.3387863) q[3];
sx q[3];
rz(-1.6646977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.431488) q[2];
sx q[2];
rz(-1.2892712) q[2];
sx q[2];
rz(2.5816494) q[2];
rz(-0.7263178) q[3];
sx q[3];
rz(-2.8328219) q[3];
sx q[3];
rz(0.30549756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2840435) q[0];
sx q[0];
rz(-2.5950268) q[0];
sx q[0];
rz(1.7204826) q[0];
rz(-0.20206085) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(2.2834159) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9439529) q[0];
sx q[0];
rz(-1.7174935) q[0];
sx q[0];
rz(-3.0118045) q[0];
x q[1];
rz(1.981295) q[2];
sx q[2];
rz(-1.6121284) q[2];
sx q[2];
rz(-1.873204) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.757829) q[1];
sx q[1];
rz(-1.541829) q[1];
sx q[1];
rz(-1.5555192) q[1];
rz(-pi) q[2];
rz(2.6303597) q[3];
sx q[3];
rz(-0.63496642) q[3];
sx q[3];
rz(-0.2583897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.28785607) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(-1.8161592) q[2];
rz(-2.251513) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(-2.1530698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4637852) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(2.7668787) q[0];
rz(2.162714) q[1];
sx q[1];
rz(-0.68190494) q[1];
sx q[1];
rz(1.3495061) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5193072) q[0];
sx q[0];
rz(-1.1847727) q[0];
sx q[0];
rz(-0.65667721) q[0];
rz(-pi) q[1];
x q[1];
rz(0.026272341) q[2];
sx q[2];
rz(-0.71825829) q[2];
sx q[2];
rz(0.19933137) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.3186036) q[1];
sx q[1];
rz(-2.9530596) q[1];
sx q[1];
rz(1.3033426) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7262444) q[3];
sx q[3];
rz(-0.58598622) q[3];
sx q[3];
rz(-1.2522445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.65138856) q[2];
sx q[2];
rz(-0.75275246) q[2];
sx q[2];
rz(-0.46869579) q[2];
rz(1.9474585) q[3];
sx q[3];
rz(-1.8374551) q[3];
sx q[3];
rz(1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5114708) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(-1.8796896) q[0];
rz(2.966554) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(1.5375686) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8309098) q[0];
sx q[0];
rz(-1.3498107) q[0];
sx q[0];
rz(1.2587147) q[0];
rz(-pi) q[1];
rz(1.4666918) q[2];
sx q[2];
rz(-1.1541919) q[2];
sx q[2];
rz(-0.57550752) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.21737145) q[1];
sx q[1];
rz(-1.6842168) q[1];
sx q[1];
rz(-0.70593112) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1101301) q[3];
sx q[3];
rz(-1.2445645) q[3];
sx q[3];
rz(-2.5903451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1016772) q[2];
sx q[2];
rz(-0.95255178) q[2];
sx q[2];
rz(2.3802479) q[2];
rz(2.2411761) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(3.0925687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.802357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(2.9246869) q[0];
rz(0.63198173) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(-0.95473081) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023708658) q[0];
sx q[0];
rz(-2.3291991) q[0];
sx q[0];
rz(-2.6931767) q[0];
rz(-pi) q[1];
rz(1.9913313) q[2];
sx q[2];
rz(-2.0805801) q[2];
sx q[2];
rz(1.2812986) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2659457) q[1];
sx q[1];
rz(-1.6469643) q[1];
sx q[1];
rz(-2.6850558) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6984152) q[3];
sx q[3];
rz(-0.9842397) q[3];
sx q[3];
rz(0.12936684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0163991) q[2];
sx q[2];
rz(-1.23896) q[2];
sx q[2];
rz(0.94474244) q[2];
rz(2.7567806) q[3];
sx q[3];
rz(-2.0231569) q[3];
sx q[3];
rz(0.95705664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.64086296) q[0];
sx q[0];
rz(-0.50518112) q[0];
sx q[0];
rz(1.5541979) q[0];
rz(-0.8846994) q[1];
sx q[1];
rz(-2.2365166) q[1];
sx q[1];
rz(2.8832163) q[1];
rz(0.89818556) q[2];
sx q[2];
rz(-2.2834416) q[2];
sx q[2];
rz(-1.7590547) q[2];
rz(1.1273884) q[3];
sx q[3];
rz(-1.756712) q[3];
sx q[3];
rz(-2.0851019) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
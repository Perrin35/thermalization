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
rz(-1.7320002) q[0];
sx q[0];
rz(-1.707466) q[0];
rz(-2.5073476) q[1];
sx q[1];
rz(-0.60159644) q[1];
sx q[1];
rz(-2.7231725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9040065) q[0];
sx q[0];
rz(-1.8338086) q[0];
sx q[0];
rz(1.0884652) q[0];
rz(1.7069874) q[2];
sx q[2];
rz(-1.681466) q[2];
sx q[2];
rz(1.1510804) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.99416713) q[1];
sx q[1];
rz(-2.6174183) q[1];
sx q[1];
rz(2.0736573) q[1];
x q[2];
rz(2.9207346) q[3];
sx q[3];
rz(-2.0005895) q[3];
sx q[3];
rz(2.5151099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3216386) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(-2.3036172) q[2];
rz(2.6485802) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(-3.0626007) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3943966) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(1.863742) q[0];
rz(2.9648119) q[1];
sx q[1];
rz(-1.8272094) q[1];
sx q[1];
rz(-0.4321672) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85053274) q[0];
sx q[0];
rz(-1.8467554) q[0];
sx q[0];
rz(-0.54131298) q[0];
x q[1];
rz(-1.997666) q[2];
sx q[2];
rz(-1.2032713) q[2];
sx q[2];
rz(1.2183684) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4914815) q[1];
sx q[1];
rz(-1.1307798) q[1];
sx q[1];
rz(-0.12932175) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1167332) q[3];
sx q[3];
rz(-2.1356574) q[3];
sx q[3];
rz(-0.75627518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54923487) q[2];
sx q[2];
rz(-1.8775512) q[2];
sx q[2];
rz(0.48669997) q[2];
rz(-1.3782079) q[3];
sx q[3];
rz(-1.8816201) q[3];
sx q[3];
rz(2.6087705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040722672) q[0];
sx q[0];
rz(-0.7464872) q[0];
sx q[0];
rz(2.7242463) q[0];
rz(-1.4886645) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(0.506385) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9837512) q[0];
sx q[0];
rz(-2.7451773) q[0];
sx q[0];
rz(-2.0435964) q[0];
rz(-pi) q[1];
rz(-1.330553) q[2];
sx q[2];
rz(-2.2908387) q[2];
sx q[2];
rz(0.70149295) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3370034) q[1];
sx q[1];
rz(-2.5924006) q[1];
sx q[1];
rz(-2.5227491) q[1];
rz(-2.9590928) q[3];
sx q[3];
rz(-1.2068519) q[3];
sx q[3];
rz(2.5516627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.69904077) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(2.55012) q[2];
rz(2.5555723) q[3];
sx q[3];
rz(-1.2089217) q[3];
sx q[3];
rz(1.4311786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9716924) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(-0.83918321) q[0];
rz(-3.1160141) q[1];
sx q[1];
rz(-2.4459116) q[1];
sx q[1];
rz(1.5930088) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4066821) q[0];
sx q[0];
rz(-1.6618177) q[0];
sx q[0];
rz(2.0216366) q[0];
rz(-2.2976539) q[2];
sx q[2];
rz(-0.91274777) q[2];
sx q[2];
rz(-2.7968) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.81739391) q[1];
sx q[1];
rz(-2.3014268) q[1];
sx q[1];
rz(2.6170931) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7110093) q[3];
sx q[3];
rz(-1.3123371) q[3];
sx q[3];
rz(0.57627288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30535355) q[2];
sx q[2];
rz(-2.2576015) q[2];
sx q[2];
rz(-0.099686064) q[2];
rz(0.95885197) q[3];
sx q[3];
rz(-1.3189664) q[3];
sx q[3];
rz(-1.3249741) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5754159) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(-0.25594041) q[0];
rz(-2.6804965) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(2.3815313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1003935) q[0];
sx q[0];
rz(-1.6738335) q[0];
sx q[0];
rz(-3.0242007) q[0];
rz(-pi) q[1];
rz(-3.0779482) q[2];
sx q[2];
rz(-1.1568937) q[2];
sx q[2];
rz(0.49583437) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45436817) q[1];
sx q[1];
rz(-1.6974653) q[1];
sx q[1];
rz(-1.441799) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7759336) q[3];
sx q[3];
rz(-2.1846111) q[3];
sx q[3];
rz(-0.97096503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5710859) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(-0.6742397) q[2];
rz(2.9267866) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(-3.1242483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53133416) q[0];
sx q[0];
rz(-1.6711618) q[0];
sx q[0];
rz(1.1791139) q[0];
rz(-0.20482652) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(-2.0746322) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63283352) q[0];
sx q[0];
rz(-1.5501223) q[0];
sx q[0];
rz(-1.585292) q[0];
rz(0.42491575) q[2];
sx q[2];
rz(-1.9746466) q[2];
sx q[2];
rz(-1.0645107) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9251717) q[1];
sx q[1];
rz(-0.74837084) q[1];
sx q[1];
rz(-1.640635) q[1];
rz(-pi) q[2];
rz(1.7339891) q[3];
sx q[3];
rz(-0.7810775) q[3];
sx q[3];
rz(1.8925325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.431488) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(-2.5816494) q[2];
rz(0.7263178) q[3];
sx q[3];
rz(-0.30877078) q[3];
sx q[3];
rz(0.30549756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.2840435) q[0];
sx q[0];
rz(-2.5950268) q[0];
sx q[0];
rz(-1.42111) q[0];
rz(-0.20206085) q[1];
sx q[1];
rz(-1.7077363) q[1];
sx q[1];
rz(-0.85817671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7875123) q[0];
sx q[0];
rz(-1.6991827) q[0];
sx q[0];
rz(1.4228729) q[0];
rz(-1.1602976) q[2];
sx q[2];
rz(-1.6121284) q[2];
sx q[2];
rz(-1.873204) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.757829) q[1];
sx q[1];
rz(-1.5997636) q[1];
sx q[1];
rz(1.5860735) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9167561) q[3];
sx q[3];
rz(-2.1145027) q[3];
sx q[3];
rz(0.35017761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8537366) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(-1.3254335) q[2];
rz(0.89007968) q[3];
sx q[3];
rz(-1.9944913) q[3];
sx q[3];
rz(-0.98852283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6778075) q[0];
sx q[0];
rz(-0.57254922) q[0];
sx q[0];
rz(-2.7668787) q[0];
rz(0.97887865) q[1];
sx q[1];
rz(-0.68190494) q[1];
sx q[1];
rz(1.7920866) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5193072) q[0];
sx q[0];
rz(-1.1847727) q[0];
sx q[0];
rz(2.4849154) q[0];
x q[1];
rz(0.71808727) q[2];
sx q[2];
rz(-1.5535083) q[2];
sx q[2];
rz(1.7503439) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.046557758) q[1];
sx q[1];
rz(-1.3890424) q[1];
sx q[1];
rz(0.05038105) q[1];
rz(-pi) q[2];
rz(-1.8324864) q[3];
sx q[3];
rz(-2.1013386) q[3];
sx q[3];
rz(-1.4025276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4902041) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(0.46869579) q[2];
rz(-1.9474585) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63012183) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(-1.8796896) q[0];
rz(0.17503861) q[1];
sx q[1];
rz(-1.1418399) q[1];
sx q[1];
rz(1.5375686) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81048548) q[0];
sx q[0];
rz(-1.2665505) q[0];
sx q[0];
rz(0.23181339) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23065718) q[2];
sx q[2];
rz(-0.42867491) q[2];
sx q[2];
rz(-0.32282695) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2214339) q[1];
sx q[1];
rz(-0.71343525) q[1];
sx q[1];
rz(2.9677797) q[1];
rz(-0.92026199) q[3];
sx q[3];
rz(-0.55763054) q[3];
sx q[3];
rz(2.695431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1016772) q[2];
sx q[2];
rz(-0.95255178) q[2];
sx q[2];
rz(-2.3802479) q[2];
rz(2.2411761) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(-0.049023978) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.802357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(-2.9246869) q[0];
rz(-2.5096109) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(2.1868618) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5548984) q[0];
sx q[0];
rz(-2.2838755) q[0];
sx q[0];
rz(-1.1416392) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54951349) q[2];
sx q[2];
rz(-1.2064484) q[2];
sx q[2];
rz(0.074631045) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.875647) q[1];
sx q[1];
rz(-1.4946283) q[1];
sx q[1];
rz(2.6850558) q[1];
rz(-0.44317742) q[3];
sx q[3];
rz(-0.9842397) q[3];
sx q[3];
rz(3.0122258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0163991) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(-2.1968502) q[2];
rz(-2.7567806) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(0.95705664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64086296) q[0];
sx q[0];
rz(-0.50518112) q[0];
sx q[0];
rz(1.5541979) q[0];
rz(0.8846994) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(-2.2434071) q[2];
sx q[2];
rz(-2.2834416) q[2];
sx q[2];
rz(-1.7590547) q[2];
rz(0.20529071) q[3];
sx q[3];
rz(-2.0060354) q[3];
sx q[3];
rz(2.5397186) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

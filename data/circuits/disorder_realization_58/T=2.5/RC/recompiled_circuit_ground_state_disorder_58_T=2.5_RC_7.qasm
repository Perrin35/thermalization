OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6405606) q[0];
sx q[0];
rz(8.6241047) q[0];
sx q[0];
rz(9.5700349) q[0];
rz(-0.89291209) q[1];
sx q[1];
rz(3.555759) q[1];
sx q[1];
rz(10.531737) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021028886) q[0];
sx q[0];
rz(-1.0146838) q[0];
sx q[0];
rz(-2.4377941) q[0];
rz(-2.8584004) q[2];
sx q[2];
rz(-1.257466) q[2];
sx q[2];
rz(1.6746132) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9913276) q[1];
sx q[1];
rz(-1.8561761) q[1];
sx q[1];
rz(-2.1726514) q[1];
rz(2.0756196) q[3];
sx q[3];
rz(-2.7144066) q[3];
sx q[3];
rz(1.8452725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5408111) q[2];
sx q[2];
rz(-1.3593707) q[2];
sx q[2];
rz(-0.092279807) q[2];
rz(-1.7021092) q[3];
sx q[3];
rz(-1.1441792) q[3];
sx q[3];
rz(0.93851844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1321024) q[0];
sx q[0];
rz(-1.1630031) q[0];
sx q[0];
rz(2.5376885) q[0];
rz(-1.5314792) q[1];
sx q[1];
rz(-1.8096626) q[1];
sx q[1];
rz(1.5276705) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85717843) q[0];
sx q[0];
rz(-0.99161327) q[0];
sx q[0];
rz(0.37838899) q[0];
rz(2.1555734) q[2];
sx q[2];
rz(-0.9169609) q[2];
sx q[2];
rz(-2.3399835) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3221047) q[1];
sx q[1];
rz(-2.8675962) q[1];
sx q[1];
rz(-0.042983965) q[1];
rz(1.5714721) q[3];
sx q[3];
rz(-1.4726032) q[3];
sx q[3];
rz(2.6718966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4332726) q[2];
sx q[2];
rz(-1.6927787) q[2];
sx q[2];
rz(-2.0937008) q[2];
rz(0.6692872) q[3];
sx q[3];
rz(-0.71988121) q[3];
sx q[3];
rz(1.9765114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3480551) q[0];
sx q[0];
rz(-2.3568643) q[0];
sx q[0];
rz(-2.134557) q[0];
rz(-0.29193613) q[1];
sx q[1];
rz(-2.597229) q[1];
sx q[1];
rz(-2.4709591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29771996) q[0];
sx q[0];
rz(-2.5150635) q[0];
sx q[0];
rz(1.421375) q[0];
x q[1];
rz(-0.88460716) q[2];
sx q[2];
rz(-2.2401407) q[2];
sx q[2];
rz(1.455292) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.933316) q[1];
sx q[1];
rz(-1.9066992) q[1];
sx q[1];
rz(-0.64888727) q[1];
x q[2];
rz(3.0257053) q[3];
sx q[3];
rz(-1.8975583) q[3];
sx q[3];
rz(0.22921697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2900419) q[2];
sx q[2];
rz(-2.2883577) q[2];
sx q[2];
rz(0.89149371) q[2];
rz(-2.0391035) q[3];
sx q[3];
rz(-2.1474371) q[3];
sx q[3];
rz(-1.4055143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82146984) q[0];
sx q[0];
rz(-1.7550884) q[0];
sx q[0];
rz(-1.0572877) q[0];
rz(3.140246) q[1];
sx q[1];
rz(-2.3206382) q[1];
sx q[1];
rz(-2.3473306) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47424537) q[0];
sx q[0];
rz(-0.66215958) q[0];
sx q[0];
rz(-0.029900877) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7873331) q[2];
sx q[2];
rz(-1.2984167) q[2];
sx q[2];
rz(0.59857063) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0483305) q[1];
sx q[1];
rz(-0.89846134) q[1];
sx q[1];
rz(0.59497084) q[1];
rz(-pi) q[2];
rz(2.7806588) q[3];
sx q[3];
rz(-2.8246701) q[3];
sx q[3];
rz(-1.288687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1069676) q[2];
sx q[2];
rz(-2.5204973) q[2];
sx q[2];
rz(1.0221488) q[2];
rz(1.243783) q[3];
sx q[3];
rz(-2.1918178) q[3];
sx q[3];
rz(1.7826084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6720402) q[0];
sx q[0];
rz(-1.38009) q[0];
sx q[0];
rz(-0.45355466) q[0];
rz(1.0376616) q[1];
sx q[1];
rz(-2.0062168) q[1];
sx q[1];
rz(0.79197788) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4290473) q[0];
sx q[0];
rz(-2.8600116) q[0];
sx q[0];
rz(1.4123807) q[0];
x q[1];
rz(2.1486695) q[2];
sx q[2];
rz(-1.2348334) q[2];
sx q[2];
rz(3.1022037) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8199948) q[1];
sx q[1];
rz(-1.454108) q[1];
sx q[1];
rz(-0.34081809) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1230917) q[3];
sx q[3];
rz(-1.7089239) q[3];
sx q[3];
rz(-0.02397315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8301293) q[2];
sx q[2];
rz(-2.4757803) q[2];
sx q[2];
rz(-2.8361481) q[2];
rz(-2.9597802) q[3];
sx q[3];
rz(-1.5287377) q[3];
sx q[3];
rz(2.4055433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55564725) q[0];
sx q[0];
rz(-0.82252994) q[0];
sx q[0];
rz(3.0162051) q[0];
rz(1.5665945) q[1];
sx q[1];
rz(-1.4488723) q[1];
sx q[1];
rz(0.032141846) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1458932) q[0];
sx q[0];
rz(-0.75135485) q[0];
sx q[0];
rz(2.8180647) q[0];
x q[1];
rz(-1.6400385) q[2];
sx q[2];
rz(-0.73261315) q[2];
sx q[2];
rz(0.68147269) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6010203) q[1];
sx q[1];
rz(-0.8972315) q[1];
sx q[1];
rz(-0.62255287) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4549082) q[3];
sx q[3];
rz(-0.82852302) q[3];
sx q[3];
rz(-0.39388408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0142168) q[2];
sx q[2];
rz(-1.2140423) q[2];
sx q[2];
rz(2.0188913) q[2];
rz(-0.073401062) q[3];
sx q[3];
rz(-0.95728907) q[3];
sx q[3];
rz(-0.61298031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4323394) q[0];
sx q[0];
rz(-1.657635) q[0];
sx q[0];
rz(0.51026979) q[0];
rz(2.7773652) q[1];
sx q[1];
rz(-2.7204456) q[1];
sx q[1];
rz(-1.4998923) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7059874) q[0];
sx q[0];
rz(-1.4650808) q[0];
sx q[0];
rz(1.1054429) q[0];
rz(-0.11830637) q[2];
sx q[2];
rz(-1.7276754) q[2];
sx q[2];
rz(-1.1821234) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8270696) q[1];
sx q[1];
rz(-1.0770969) q[1];
sx q[1];
rz(1.3160454) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33320547) q[3];
sx q[3];
rz(-2.3764046) q[3];
sx q[3];
rz(3.0195587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.69786543) q[2];
sx q[2];
rz(-1.5550193) q[2];
sx q[2];
rz(2.2743684) q[2];
rz(1.5602268) q[3];
sx q[3];
rz(-2.9428704) q[3];
sx q[3];
rz(-1.3377415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4174058) q[0];
sx q[0];
rz(-3.0751808) q[0];
sx q[0];
rz(2.7253286) q[0];
rz(-1.1626214) q[1];
sx q[1];
rz(-1.1888209) q[1];
sx q[1];
rz(-2.3847041) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74706057) q[0];
sx q[0];
rz(-1.146527) q[0];
sx q[0];
rz(0.31655689) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95332884) q[2];
sx q[2];
rz(-0.67425113) q[2];
sx q[2];
rz(1.9566388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50584201) q[1];
sx q[1];
rz(-1.7598087) q[1];
sx q[1];
rz(-2.0057949) q[1];
x q[2];
rz(-2.5060095) q[3];
sx q[3];
rz(-1.9548237) q[3];
sx q[3];
rz(1.7049198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68652144) q[2];
sx q[2];
rz(-1.4342118) q[2];
sx q[2];
rz(-1.3580458) q[2];
rz(2.3250735) q[3];
sx q[3];
rz(-1.2314545) q[3];
sx q[3];
rz(0.16286287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74129504) q[0];
sx q[0];
rz(-1.6930027) q[0];
sx q[0];
rz(-1.4073538) q[0];
rz(-0.74527144) q[1];
sx q[1];
rz(-2.407357) q[1];
sx q[1];
rz(-1.5819246) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0228221) q[0];
sx q[0];
rz(-1.8466936) q[0];
sx q[0];
rz(0.89275443) q[0];
rz(-3.09868) q[2];
sx q[2];
rz(-2.4933698) q[2];
sx q[2];
rz(-0.15951482) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7551387) q[1];
sx q[1];
rz(-1.6999131) q[1];
sx q[1];
rz(2.0172202) q[1];
rz(-pi) q[2];
rz(0.060272597) q[3];
sx q[3];
rz(-2.1456686) q[3];
sx q[3];
rz(-2.4337492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.72426307) q[2];
sx q[2];
rz(-1.686692) q[2];
sx q[2];
rz(2.0261185) q[2];
rz(-2.8880902) q[3];
sx q[3];
rz(-1.2616254) q[3];
sx q[3];
rz(0.0089664627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.020092) q[0];
sx q[0];
rz(-1.2637063) q[0];
sx q[0];
rz(-2.3790835) q[0];
rz(-2.1992042) q[1];
sx q[1];
rz(-1.0647048) q[1];
sx q[1];
rz(1.2368088) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.002703) q[0];
sx q[0];
rz(-1.5093439) q[0];
sx q[0];
rz(-0.49646722) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57970826) q[2];
sx q[2];
rz(-0.97203185) q[2];
sx q[2];
rz(2.6331462) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2357465) q[1];
sx q[1];
rz(-1.5457834) q[1];
sx q[1];
rz(0.025084875) q[1];
rz(1.7364301) q[3];
sx q[3];
rz(-0.40963033) q[3];
sx q[3];
rz(1.7543242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2463871) q[2];
sx q[2];
rz(-0.095194101) q[2];
sx q[2];
rz(3.134356) q[2];
rz(-0.17035189) q[3];
sx q[3];
rz(-1.4879613) q[3];
sx q[3];
rz(-0.65792221) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8728747) q[0];
sx q[0];
rz(-1.7418516) q[0];
sx q[0];
rz(-2.0280784) q[0];
rz(0.089182236) q[1];
sx q[1];
rz(-1.5166278) q[1];
sx q[1];
rz(1.2022432) q[1];
rz(-2.4287672) q[2];
sx q[2];
rz(-0.20607866) q[2];
sx q[2];
rz(-1.1417749) q[2];
rz(-1.3842267) q[3];
sx q[3];
rz(-1.2751147) q[3];
sx q[3];
rz(-0.57300344) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(0.060184181) q[0];
sx q[0];
rz(-2.0826075) q[0];
sx q[0];
rz(-1.0101779) q[0];
rz(-0.8085568) q[1];
sx q[1];
rz(2.8454236) q[1];
sx q[1];
rz(12.256395) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5910019) q[0];
sx q[0];
rz(-0.43568107) q[0];
sx q[0];
rz(0.26623078) q[0];
rz(0.10711889) q[2];
sx q[2];
rz(-2.9351882) q[2];
sx q[2];
rz(0.53910461) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9604608) q[1];
sx q[1];
rz(-0.38552654) q[1];
sx q[1];
rz(1.0268282) q[1];
rz(-pi) q[2];
rz(-2.9397291) q[3];
sx q[3];
rz(-0.84175693) q[3];
sx q[3];
rz(2.1191459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6628722) q[2];
sx q[2];
rz(-2.8933849) q[2];
sx q[2];
rz(-3.0459246) q[2];
rz(3.0293448) q[3];
sx q[3];
rz(-0.87533689) q[3];
sx q[3];
rz(1.7101425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2410759) q[0];
sx q[0];
rz(-0.88923419) q[0];
sx q[0];
rz(-2.526793) q[0];
rz(2.2531033) q[1];
sx q[1];
rz(-1.588984) q[1];
sx q[1];
rz(2.6420171) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7730656) q[0];
sx q[0];
rz(-1.245226) q[0];
sx q[0];
rz(-0.21296176) q[0];
rz(-pi) q[1];
rz(-0.52421661) q[2];
sx q[2];
rz(-1.7078064) q[2];
sx q[2];
rz(1.7643339) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72909468) q[1];
sx q[1];
rz(-2.4654268) q[1];
sx q[1];
rz(-1.2926284) q[1];
rz(2.8039475) q[3];
sx q[3];
rz(-2.1269375) q[3];
sx q[3];
rz(-0.82586702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2483612) q[2];
sx q[2];
rz(-1.2373368) q[2];
sx q[2];
rz(0.020817967) q[2];
rz(-1.9013532) q[3];
sx q[3];
rz(-1.6720684) q[3];
sx q[3];
rz(-1.1851236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5027387) q[0];
sx q[0];
rz(-0.26832142) q[0];
sx q[0];
rz(-0.33367208) q[0];
rz(-2.4036713) q[1];
sx q[1];
rz(-1.9155904) q[1];
sx q[1];
rz(-3.0121682) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0282818) q[0];
sx q[0];
rz(-2.2419562) q[0];
sx q[0];
rz(-2.9379803) q[0];
rz(-2.4366662) q[2];
sx q[2];
rz(-0.80653301) q[2];
sx q[2];
rz(0.091900983) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.94679615) q[1];
sx q[1];
rz(-1.4490286) q[1];
sx q[1];
rz(0.090784723) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6629642) q[3];
sx q[3];
rz(-0.84175292) q[3];
sx q[3];
rz(0.42985502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1239329) q[2];
sx q[2];
rz(-1.61597) q[2];
sx q[2];
rz(-1.7714436) q[2];
rz(0.74639368) q[3];
sx q[3];
rz(-1.8236225) q[3];
sx q[3];
rz(0.082775041) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.133404) q[0];
sx q[0];
rz(-1.2097825) q[0];
sx q[0];
rz(-2.5764537) q[0];
rz(0.42916974) q[1];
sx q[1];
rz(-0.64650911) q[1];
sx q[1];
rz(-2.3133004) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9681184) q[0];
sx q[0];
rz(-2.3096363) q[0];
sx q[0];
rz(1.678874) q[0];
rz(-1.0847237) q[2];
sx q[2];
rz(-1.3157433) q[2];
sx q[2];
rz(1.6549095) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8737303) q[1];
sx q[1];
rz(-0.56247382) q[1];
sx q[1];
rz(2.6270694) q[1];
rz(-pi) q[2];
rz(2.0873484) q[3];
sx q[3];
rz(-1.787546) q[3];
sx q[3];
rz(-2.3473489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68507489) q[2];
sx q[2];
rz(-1.061331) q[2];
sx q[2];
rz(1.7100517) q[2];
rz(-0.93959129) q[3];
sx q[3];
rz(-2.6288433) q[3];
sx q[3];
rz(0.00069869839) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13570304) q[0];
sx q[0];
rz(-0.26517427) q[0];
sx q[0];
rz(-1.8121207) q[0];
rz(1.9895408) q[1];
sx q[1];
rz(-1.0877129) q[1];
sx q[1];
rz(-3.1413445) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5329929) q[0];
sx q[0];
rz(-1.8155351) q[0];
sx q[0];
rz(-0.18969638) q[0];
rz(-pi) q[1];
x q[1];
rz(1.423944) q[2];
sx q[2];
rz(-0.53895437) q[2];
sx q[2];
rz(0.94648628) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2622288) q[1];
sx q[1];
rz(-1.9161703) q[1];
sx q[1];
rz(-1.0443347) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6432006) q[3];
sx q[3];
rz(-0.84057099) q[3];
sx q[3];
rz(2.8080468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1382711) q[2];
sx q[2];
rz(-1.891581) q[2];
sx q[2];
rz(3.1357583) q[2];
rz(2.2516069) q[3];
sx q[3];
rz(-2.4599288) q[3];
sx q[3];
rz(-1.1267004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936845) q[0];
sx q[0];
rz(-1.6981145) q[0];
sx q[0];
rz(-1.6538612) q[0];
rz(0.030390175) q[1];
sx q[1];
rz(-2.0149714) q[1];
sx q[1];
rz(-2.1991275) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4543991) q[0];
sx q[0];
rz(-1.8158578) q[0];
sx q[0];
rz(1.3769727) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29877383) q[2];
sx q[2];
rz(-1.818383) q[2];
sx q[2];
rz(-0.58591671) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.16557517) q[1];
sx q[1];
rz(-1.8569131) q[1];
sx q[1];
rz(-1.0817492) q[1];
rz(-pi) q[2];
rz(-1.0932572) q[3];
sx q[3];
rz(-0.52553383) q[3];
sx q[3];
rz(0.20697396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5421062) q[2];
sx q[2];
rz(-0.78080559) q[2];
sx q[2];
rz(-1.3055118) q[2];
rz(0.54715884) q[3];
sx q[3];
rz(-1.9360417) q[3];
sx q[3];
rz(-2.3356596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.9611573) q[0];
sx q[0];
rz(-2.1077709) q[0];
sx q[0];
rz(0.10051522) q[0];
rz(2.6590977) q[1];
sx q[1];
rz(-1.1588187) q[1];
sx q[1];
rz(-0.85711342) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8554358) q[0];
sx q[0];
rz(-1.8925987) q[0];
sx q[0];
rz(-2.4706173) q[0];
x q[1];
rz(0.64684644) q[2];
sx q[2];
rz(-1.610272) q[2];
sx q[2];
rz(-2.0496862) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.313011) q[1];
sx q[1];
rz(-1.9691113) q[1];
sx q[1];
rz(1.7969153) q[1];
rz(-pi) q[2];
rz(0.25730142) q[3];
sx q[3];
rz(-2.1273779) q[3];
sx q[3];
rz(-2.0428994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7414005) q[2];
sx q[2];
rz(-2.1705748) q[2];
sx q[2];
rz(-2.8391489) q[2];
rz(0.97964573) q[3];
sx q[3];
rz(-1.0023578) q[3];
sx q[3];
rz(1.8625331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7804467) q[0];
sx q[0];
rz(-3.0028711) q[0];
sx q[0];
rz(-0.030990344) q[0];
rz(0.57394761) q[1];
sx q[1];
rz(-1.4731044) q[1];
sx q[1];
rz(-1.2773638) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4784067) q[0];
sx q[0];
rz(-2.660523) q[0];
sx q[0];
rz(-2.3548954) q[0];
rz(-0.24279109) q[2];
sx q[2];
rz(-2.0567679) q[2];
sx q[2];
rz(0.15835855) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.053169202) q[1];
sx q[1];
rz(-0.49234566) q[1];
sx q[1];
rz(-1.4150934) q[1];
rz(-pi) q[2];
rz(0.061789024) q[3];
sx q[3];
rz(-1.8553858) q[3];
sx q[3];
rz(-2.9348843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0370827) q[2];
sx q[2];
rz(-3.0057378) q[2];
sx q[2];
rz(-0.48745298) q[2];
rz(-0.69665748) q[3];
sx q[3];
rz(-2.2127377) q[3];
sx q[3];
rz(-3.0776183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071851991) q[0];
sx q[0];
rz(-0.47436473) q[0];
sx q[0];
rz(-0.35097861) q[0];
rz(0.17414302) q[1];
sx q[1];
rz(-1.6330279) q[1];
sx q[1];
rz(-1.4016271) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1202209) q[0];
sx q[0];
rz(-1.8184156) q[0];
sx q[0];
rz(1.6970474) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0171595) q[2];
sx q[2];
rz(-0.87136641) q[2];
sx q[2];
rz(2.1289189) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3199215) q[1];
sx q[1];
rz(-0.18583365) q[1];
sx q[1];
rz(-1.7286848) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4003567) q[3];
sx q[3];
rz(-1.110807) q[3];
sx q[3];
rz(2.7335087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.031124) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(-0.6366716) q[2];
rz(-2.0311671) q[3];
sx q[3];
rz(-1.597155) q[3];
sx q[3];
rz(-2.6231664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7745895) q[0];
sx q[0];
rz(-0.29357266) q[0];
sx q[0];
rz(1.0585744) q[0];
rz(0.038657945) q[1];
sx q[1];
rz(-1.6148753) q[1];
sx q[1];
rz(-2.0775332) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3389726) q[0];
sx q[0];
rz(-1.5741328) q[0];
sx q[0];
rz(0.0044857684) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8644807) q[2];
sx q[2];
rz(-1.2925576) q[2];
sx q[2];
rz(0.81467512) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7239889) q[1];
sx q[1];
rz(-0.072712459) q[1];
sx q[1];
rz(-2.5515208) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5359155) q[3];
sx q[3];
rz(-2.0300755) q[3];
sx q[3];
rz(1.1920795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4260063) q[2];
sx q[2];
rz(-2.6435489) q[2];
sx q[2];
rz(-3.0675724) q[2];
rz(2.7600539) q[3];
sx q[3];
rz(-1.8266725) q[3];
sx q[3];
rz(0.35416245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.030180177) q[0];
sx q[0];
rz(-1.3127865) q[0];
sx q[0];
rz(2.4319613) q[0];
rz(0.61182712) q[1];
sx q[1];
rz(-2.430293) q[1];
sx q[1];
rz(-1.5878955) q[1];
rz(0.62803531) q[2];
sx q[2];
rz(-0.59878329) q[2];
sx q[2];
rz(-1.5046635) q[2];
rz(2.6388219) q[3];
sx q[3];
rz(-2.3270861) q[3];
sx q[3];
rz(2.6181639) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

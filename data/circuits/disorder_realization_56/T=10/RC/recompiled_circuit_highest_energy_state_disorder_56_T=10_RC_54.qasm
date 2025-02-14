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
rz(-0.40060488) q[0];
sx q[0];
rz(-2.7317943) q[0];
sx q[0];
rz(2.0706489) q[0];
rz(0.61869705) q[1];
sx q[1];
rz(3.7096042) q[1];
sx q[1];
rz(8.5811442) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7535665) q[0];
sx q[0];
rz(-0.26615289) q[0];
sx q[0];
rz(-1.4583336) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6981612) q[2];
sx q[2];
rz(-0.3249692) q[2];
sx q[2];
rz(-2.5814499) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8275717) q[1];
sx q[1];
rz(-2.2306475) q[1];
sx q[1];
rz(-2.72675) q[1];
rz(-pi) q[2];
rz(1.0759962) q[3];
sx q[3];
rz(-1.9952979) q[3];
sx q[3];
rz(0.44874661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.3598651) q[2];
sx q[2];
rz(-0.83911506) q[2];
sx q[2];
rz(3.1245933) q[2];
rz(-1.7403691) q[3];
sx q[3];
rz(-0.56178105) q[3];
sx q[3];
rz(1.2167654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5725937) q[0];
sx q[0];
rz(-1.3324791) q[0];
sx q[0];
rz(0.78729415) q[0];
rz(0.17838082) q[1];
sx q[1];
rz(-1.9352103) q[1];
sx q[1];
rz(1.9040727) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2923861) q[0];
sx q[0];
rz(-0.9885177) q[0];
sx q[0];
rz(1.4968064) q[0];
x q[1];
rz(-1.9745047) q[2];
sx q[2];
rz(-2.4879527) q[2];
sx q[2];
rz(0.47932759) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.54531714) q[1];
sx q[1];
rz(-1.8870237) q[1];
sx q[1];
rz(1.113722) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54333289) q[3];
sx q[3];
rz(-0.50341922) q[3];
sx q[3];
rz(-1.0375298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7201207) q[2];
sx q[2];
rz(-1.6933491) q[2];
sx q[2];
rz(-0.063701542) q[2];
rz(1.2740159) q[3];
sx q[3];
rz(-2.2985022) q[3];
sx q[3];
rz(-1.5773076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.098123) q[0];
sx q[0];
rz(-1.948057) q[0];
sx q[0];
rz(-0.7435588) q[0];
rz(-2.6877563) q[1];
sx q[1];
rz(-1.791714) q[1];
sx q[1];
rz(0.41890621) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.009368816) q[0];
sx q[0];
rz(-2.0133349) q[0];
sx q[0];
rz(2.924463) q[0];
x q[1];
rz(-0.30861993) q[2];
sx q[2];
rz(-2.2815124) q[2];
sx q[2];
rz(-2.1580704) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3943399) q[1];
sx q[1];
rz(-1.2868007) q[1];
sx q[1];
rz(-0.92628493) q[1];
rz(-pi) q[2];
rz(0.13678582) q[3];
sx q[3];
rz(-0.89623755) q[3];
sx q[3];
rz(-3.0709895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0539187) q[2];
sx q[2];
rz(-0.81834617) q[2];
sx q[2];
rz(-1.9107001) q[2];
rz(-2.6026717) q[3];
sx q[3];
rz(-1.6659707) q[3];
sx q[3];
rz(1.4628791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.8777799) q[0];
sx q[0];
rz(-2.3872264) q[0];
sx q[0];
rz(-0.60212773) q[0];
rz(1.6944132) q[1];
sx q[1];
rz(-2.4592631) q[1];
sx q[1];
rz(0.86300659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38463889) q[0];
sx q[0];
rz(-2.2881329) q[0];
sx q[0];
rz(1.4199281) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30089897) q[2];
sx q[2];
rz(-1.2105807) q[2];
sx q[2];
rz(0.27733251) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.26559184) q[1];
sx q[1];
rz(-2.1419898) q[1];
sx q[1];
rz(2.8286562) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3087048) q[3];
sx q[3];
rz(-1.8084267) q[3];
sx q[3];
rz(0.72494635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0756695) q[2];
sx q[2];
rz(-1.1752335) q[2];
sx q[2];
rz(2.2389331) q[2];
rz(2.3265808) q[3];
sx q[3];
rz(-2.5383526) q[3];
sx q[3];
rz(-0.41316113) q[3];
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
rz(1.0937423) q[0];
sx q[0];
rz(-0.047690064) q[0];
sx q[0];
rz(-0.93511859) q[0];
rz(-1.5330261) q[1];
sx q[1];
rz(-0.32564274) q[1];
sx q[1];
rz(-2.8757222) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2926798) q[0];
sx q[0];
rz(-0.32912185) q[0];
sx q[0];
rz(1.2364682) q[0];
rz(0.75030783) q[2];
sx q[2];
rz(-1.3641832) q[2];
sx q[2];
rz(2.7624102) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.26755493) q[1];
sx q[1];
rz(-1.7312164) q[1];
sx q[1];
rz(-2.7880413) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3060027) q[3];
sx q[3];
rz(-1.6700498) q[3];
sx q[3];
rz(-2.3276727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.461146) q[2];
sx q[2];
rz(-1.5004044) q[2];
sx q[2];
rz(-0.45073304) q[2];
rz(1.1905253) q[3];
sx q[3];
rz(-0.7014941) q[3];
sx q[3];
rz(0.43959555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.180069) q[0];
sx q[0];
rz(-2.838205) q[0];
sx q[0];
rz(-3.0971089) q[0];
rz(-0.26875177) q[1];
sx q[1];
rz(-2.2046397) q[1];
sx q[1];
rz(-2.7106947) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.160153) q[0];
sx q[0];
rz(-0.68228506) q[0];
sx q[0];
rz(-0.42295608) q[0];
rz(-2.2399432) q[2];
sx q[2];
rz(-1.2001654) q[2];
sx q[2];
rz(-1.8805875) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9932672) q[1];
sx q[1];
rz(-1.1701411) q[1];
sx q[1];
rz(2.2877778) q[1];
x q[2];
rz(-1.5650125) q[3];
sx q[3];
rz(-0.91051379) q[3];
sx q[3];
rz(-2.3920812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0253133) q[2];
sx q[2];
rz(-2.7677324) q[2];
sx q[2];
rz(-2.5171793) q[2];
rz(-2.0404909) q[3];
sx q[3];
rz(-0.81239429) q[3];
sx q[3];
rz(-2.1541434) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2295912) q[0];
sx q[0];
rz(-1.0593375) q[0];
sx q[0];
rz(-0.68369317) q[0];
rz(1.4423485) q[1];
sx q[1];
rz(-2.3577299) q[1];
sx q[1];
rz(-0.20194617) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5366906) q[0];
sx q[0];
rz(-2.8134007) q[0];
sx q[0];
rz(-1.5668014) q[0];
rz(2.3181386) q[2];
sx q[2];
rz(-1.5353893) q[2];
sx q[2];
rz(3.0454783) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3298885) q[1];
sx q[1];
rz(-2.0913908) q[1];
sx q[1];
rz(-1.64774) q[1];
rz(-pi) q[2];
rz(1.6138329) q[3];
sx q[3];
rz(-2.0282054) q[3];
sx q[3];
rz(0.73107536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0553637) q[2];
sx q[2];
rz(-1.7217041) q[2];
sx q[2];
rz(-1.281338) q[2];
rz(-0.66644871) q[3];
sx q[3];
rz(-1.0969578) q[3];
sx q[3];
rz(2.5950529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1619038) q[0];
sx q[0];
rz(-2.5894916) q[0];
sx q[0];
rz(-2.6276278) q[0];
rz(-2.1886096) q[1];
sx q[1];
rz(-0.62478462) q[1];
sx q[1];
rz(2.0228588) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1361439) q[0];
sx q[0];
rz(-1.4892254) q[0];
sx q[0];
rz(2.7177627) q[0];
rz(-pi) q[1];
rz(0.4301901) q[2];
sx q[2];
rz(-2.3025844) q[2];
sx q[2];
rz(-2.6553287) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7756164) q[1];
sx q[1];
rz(-1.9785641) q[1];
sx q[1];
rz(-2.9080176) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5947078) q[3];
sx q[3];
rz(-0.62084475) q[3];
sx q[3];
rz(1.3399233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3861367) q[2];
sx q[2];
rz(-0.35217199) q[2];
sx q[2];
rz(0.90292162) q[2];
rz(-1.6821945) q[3];
sx q[3];
rz(-1.3420339) q[3];
sx q[3];
rz(0.82227796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.67895401) q[0];
sx q[0];
rz(-0.32005388) q[0];
sx q[0];
rz(1.951304) q[0];
rz(0.32084385) q[1];
sx q[1];
rz(-1.6879993) q[1];
sx q[1];
rz(-1.2215337) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7090593) q[0];
sx q[0];
rz(-1.4343059) q[0];
sx q[0];
rz(0.20705072) q[0];
x q[1];
rz(-2.8655445) q[2];
sx q[2];
rz(-0.85535565) q[2];
sx q[2];
rz(-2.9710342) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22449422) q[1];
sx q[1];
rz(-2.3116744) q[1];
sx q[1];
rz(1.8486449) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3873877) q[3];
sx q[3];
rz(-0.26421031) q[3];
sx q[3];
rz(-2.9809706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90423501) q[2];
sx q[2];
rz(-2.9866437) q[2];
sx q[2];
rz(0.3332738) q[2];
rz(-1.0171558) q[3];
sx q[3];
rz(-2.1867496) q[3];
sx q[3];
rz(2.8216383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79621133) q[0];
sx q[0];
rz(-2.6295202) q[0];
sx q[0];
rz(0.91878015) q[0];
rz(-2.9504919) q[1];
sx q[1];
rz(-0.42593503) q[1];
sx q[1];
rz(0.79413116) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92768663) q[0];
sx q[0];
rz(-2.0964668) q[0];
sx q[0];
rz(-0.77362432) q[0];
rz(1.3599858) q[2];
sx q[2];
rz(-1.2173869) q[2];
sx q[2];
rz(-1.4251874) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5672553) q[1];
sx q[1];
rz(-2.2681288) q[1];
sx q[1];
rz(2.789807) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32148949) q[3];
sx q[3];
rz(-1.7712542) q[3];
sx q[3];
rz(3.0908302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11067757) q[2];
sx q[2];
rz(-0.6157178) q[2];
sx q[2];
rz(-0.75882971) q[2];
rz(-2.2714254) q[3];
sx q[3];
rz(-2.3535959) q[3];
sx q[3];
rz(-2.0738535) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1398685) q[0];
sx q[0];
rz(-1.4586466) q[0];
sx q[0];
rz(1.9139105) q[0];
rz(-2.7036746) q[1];
sx q[1];
rz(-1.4358078) q[1];
sx q[1];
rz(-1.5461071) q[1];
rz(0.48702892) q[2];
sx q[2];
rz(-1.9430046) q[2];
sx q[2];
rz(-0.76443048) q[2];
rz(0.78183635) q[3];
sx q[3];
rz(-1.1734097) q[3];
sx q[3];
rz(2.9076037) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

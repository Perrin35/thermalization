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
rz(0.33557284) q[0];
sx q[0];
rz(4.8290401) q[0];
sx q[0];
rz(7.3168559) q[0];
rz(1.4312862) q[1];
sx q[1];
rz(-0.55748504) q[1];
sx q[1];
rz(-1.0097591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.570734) q[0];
sx q[0];
rz(-1.6964127) q[0];
sx q[0];
rz(-1.3799589) q[0];
rz(-pi) q[1];
rz(0.0037438914) q[2];
sx q[2];
rz(-1.6532678) q[2];
sx q[2];
rz(-1.0325026) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7518471) q[1];
sx q[1];
rz(-1.455014) q[1];
sx q[1];
rz(2.9751361) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14222957) q[3];
sx q[3];
rz(-1.5664219) q[3];
sx q[3];
rz(1.7612918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7734163) q[2];
sx q[2];
rz(-2.2645617) q[2];
sx q[2];
rz(-0.82484335) q[2];
rz(-2.8090737) q[3];
sx q[3];
rz(-0.85086346) q[3];
sx q[3];
rz(0.48377812) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9124209) q[0];
sx q[0];
rz(-0.084675463) q[0];
sx q[0];
rz(0.53572768) q[0];
rz(-0.76672211) q[1];
sx q[1];
rz(-1.7886536) q[1];
sx q[1];
rz(-1.3923233) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9014382) q[0];
sx q[0];
rz(-1.2291698) q[0];
sx q[0];
rz(2.065395) q[0];
x q[1];
rz(2.8004592) q[2];
sx q[2];
rz(-1.5798414) q[2];
sx q[2];
rz(1.9655399) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5114047) q[1];
sx q[1];
rz(-1.5923481) q[1];
sx q[1];
rz(-0.81540108) q[1];
rz(0.68269888) q[3];
sx q[3];
rz(-2.9761752) q[3];
sx q[3];
rz(1.9660275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9048189) q[2];
sx q[2];
rz(-0.83319131) q[2];
sx q[2];
rz(-2.0429677) q[2];
rz(2.3084579) q[3];
sx q[3];
rz(-1.5980709) q[3];
sx q[3];
rz(2.0875077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5611834) q[0];
sx q[0];
rz(-1.1495178) q[0];
sx q[0];
rz(2.4533601) q[0];
rz(1.2511823) q[1];
sx q[1];
rz(-0.65443188) q[1];
sx q[1];
rz(1.2122663) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1935342) q[0];
sx q[0];
rz(-2.2431264) q[0];
sx q[0];
rz(-0.95114077) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86264865) q[2];
sx q[2];
rz(-2.7852163) q[2];
sx q[2];
rz(1.732812) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7408091) q[1];
sx q[1];
rz(-2.2297528) q[1];
sx q[1];
rz(1.5779005) q[1];
rz(-3.1367407) q[3];
sx q[3];
rz(-1.650164) q[3];
sx q[3];
rz(0.30747372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4091829) q[2];
sx q[2];
rz(-0.99486351) q[2];
sx q[2];
rz(-2.691972) q[2];
rz(2.2879587) q[3];
sx q[3];
rz(-2.4946404) q[3];
sx q[3];
rz(0.70931119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7483826) q[0];
sx q[0];
rz(-1.0164096) q[0];
sx q[0];
rz(-0.40312314) q[0];
rz(-0.77278167) q[1];
sx q[1];
rz(-1.9085596) q[1];
sx q[1];
rz(2.3162139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7220424) q[0];
sx q[0];
rz(-2.3539676) q[0];
sx q[0];
rz(2.032729) q[0];
rz(-pi) q[1];
rz(3.1146997) q[2];
sx q[2];
rz(-1.5468569) q[2];
sx q[2];
rz(-2.8223512) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5063012) q[1];
sx q[1];
rz(-1.4951733) q[1];
sx q[1];
rz(2.1525394) q[1];
rz(-pi) q[2];
rz(2.230497) q[3];
sx q[3];
rz(-2.0269217) q[3];
sx q[3];
rz(-2.5191189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1021425) q[2];
sx q[2];
rz(-1.2692229) q[2];
sx q[2];
rz(0.5086745) q[2];
rz(-1.8968808) q[3];
sx q[3];
rz(-2.0679943) q[3];
sx q[3];
rz(-1.8083474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73388571) q[0];
sx q[0];
rz(-1.5835967) q[0];
sx q[0];
rz(-0.35633126) q[0];
rz(1.412926) q[1];
sx q[1];
rz(-1.8316385) q[1];
sx q[1];
rz(-1.8467356) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.750941) q[0];
sx q[0];
rz(-1.7897507) q[0];
sx q[0];
rz(0.44415565) q[0];
rz(-pi) q[1];
rz(1.7412498) q[2];
sx q[2];
rz(-0.91636412) q[2];
sx q[2];
rz(2.7026388) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7562189) q[1];
sx q[1];
rz(-1.9665779) q[1];
sx q[1];
rz(-1.8204921) q[1];
rz(0.44900198) q[3];
sx q[3];
rz(-1.0706361) q[3];
sx q[3];
rz(-1.7412108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20092189) q[2];
sx q[2];
rz(-1.7565497) q[2];
sx q[2];
rz(-0.6400288) q[2];
rz(-2.707543) q[3];
sx q[3];
rz(-2.1899352) q[3];
sx q[3];
rz(-2.0210338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.5459179) q[0];
sx q[0];
rz(-0.43788236) q[0];
sx q[0];
rz(0.75089279) q[0];
rz(0.76260507) q[1];
sx q[1];
rz(-1.7060988) q[1];
sx q[1];
rz(-2.6079752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5845244) q[0];
sx q[0];
rz(-1.2831339) q[0];
sx q[0];
rz(0.30867048) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4328753) q[2];
sx q[2];
rz(-2.2895669) q[2];
sx q[2];
rz(1.5611888) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68249629) q[1];
sx q[1];
rz(-0.67855103) q[1];
sx q[1];
rz(1.7595923) q[1];
rz(2.5848006) q[3];
sx q[3];
rz(-2.0813312) q[3];
sx q[3];
rz(-0.42312121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0642455) q[2];
sx q[2];
rz(-0.40234819) q[2];
sx q[2];
rz(0.8858436) q[2];
rz(-1.3028076) q[3];
sx q[3];
rz(-1.9582483) q[3];
sx q[3];
rz(-2.6634789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.2122022) q[0];
sx q[0];
rz(-0.92340702) q[0];
sx q[0];
rz(-0.59342629) q[0];
rz(-1.628283) q[1];
sx q[1];
rz(-1.9624886) q[1];
sx q[1];
rz(2.5679307) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79686576) q[0];
sx q[0];
rz(-1.0930976) q[0];
sx q[0];
rz(-1.1487238) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5638177) q[2];
sx q[2];
rz(-1.7794357) q[2];
sx q[2];
rz(-1.8633757) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1036914) q[1];
sx q[1];
rz(-1.9677487) q[1];
sx q[1];
rz(1.8381702) q[1];
x q[2];
rz(-0.57145561) q[3];
sx q[3];
rz(-0.09342362) q[3];
sx q[3];
rz(-2.759552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4297428) q[2];
sx q[2];
rz(-1.9387127) q[2];
sx q[2];
rz(1.708606) q[2];
rz(1.1687219) q[3];
sx q[3];
rz(-0.43787268) q[3];
sx q[3];
rz(-1.5168813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.2699921) q[0];
sx q[0];
rz(-2.2347436) q[0];
sx q[0];
rz(-2.9378743) q[0];
rz(-1.3153971) q[1];
sx q[1];
rz(-1.8351277) q[1];
sx q[1];
rz(1.3320097) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7527134) q[0];
sx q[0];
rz(-2.1915771) q[0];
sx q[0];
rz(-0.78691532) q[0];
rz(-pi) q[1];
rz(-0.90317995) q[2];
sx q[2];
rz(-1.1051264) q[2];
sx q[2];
rz(-2.6671034) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89573025) q[1];
sx q[1];
rz(-1.9612211) q[1];
sx q[1];
rz(0.67710442) q[1];
rz(-1.733794) q[3];
sx q[3];
rz(-1.6987112) q[3];
sx q[3];
rz(-0.35234141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77067644) q[2];
sx q[2];
rz(-3.0597661) q[2];
sx q[2];
rz(1.7601684) q[2];
rz(2.5904739) q[3];
sx q[3];
rz(-1.8128606) q[3];
sx q[3];
rz(2.3962925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4296752) q[0];
sx q[0];
rz(-1.6918007) q[0];
sx q[0];
rz(1.2592738) q[0];
rz(1.7970386) q[1];
sx q[1];
rz(-0.63251907) q[1];
sx q[1];
rz(1.22619) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2215458) q[0];
sx q[0];
rz(-1.4355005) q[0];
sx q[0];
rz(1.8578803) q[0];
x q[1];
rz(-2.2974469) q[2];
sx q[2];
rz(-2.591989) q[2];
sx q[2];
rz(1.5908102) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3873432) q[1];
sx q[1];
rz(-2.3918536) q[1];
sx q[1];
rz(0.096258817) q[1];
x q[2];
rz(2.0894946) q[3];
sx q[3];
rz(-1.8648793) q[3];
sx q[3];
rz(0.034386612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58174497) q[2];
sx q[2];
rz(-0.86220828) q[2];
sx q[2];
rz(1.5398514) q[2];
rz(1.7848484) q[3];
sx q[3];
rz(-0.53250766) q[3];
sx q[3];
rz(-1.9073568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27720472) q[0];
sx q[0];
rz(-1.5437523) q[0];
sx q[0];
rz(0.10449115) q[0];
rz(-1.3970207) q[1];
sx q[1];
rz(-1.3518159) q[1];
sx q[1];
rz(-1.5207965) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45986555) q[0];
sx q[0];
rz(-0.84375536) q[0];
sx q[0];
rz(1.8521502) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42664032) q[2];
sx q[2];
rz(-0.37084118) q[2];
sx q[2];
rz(-1.9538572) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7193774) q[1];
sx q[1];
rz(-1.6746192) q[1];
sx q[1];
rz(0.28909282) q[1];
x q[2];
rz(-1.4827221) q[3];
sx q[3];
rz(-2.5582426) q[3];
sx q[3];
rz(-2.6627735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4349159) q[2];
sx q[2];
rz(-0.84182635) q[2];
sx q[2];
rz(-1.8527156) q[2];
rz(-1.0051109) q[3];
sx q[3];
rz(-1.0594599) q[3];
sx q[3];
rz(0.10678261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2224251) q[0];
sx q[0];
rz(-1.5821624) q[0];
sx q[0];
rz(-0.17056175) q[0];
rz(0.87986058) q[1];
sx q[1];
rz(-2.6251371) q[1];
sx q[1];
rz(0.62216204) q[1];
rz(3.0237985) q[2];
sx q[2];
rz(-2.6391023) q[2];
sx q[2];
rz(-0.53822623) q[2];
rz(-0.60012695) q[3];
sx q[3];
rz(-1.6356346) q[3];
sx q[3];
rz(0.076250565) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

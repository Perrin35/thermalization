OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.29326987) q[0];
sx q[0];
rz(-2.8889416) q[0];
sx q[0];
rz(1.2055612) q[0];
rz(2.0134917) q[1];
sx q[1];
rz(-1.3265346) q[1];
sx q[1];
rz(0.74572745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099961258) q[0];
sx q[0];
rz(-1.3486762) q[0];
sx q[0];
rz(-0.061454031) q[0];
rz(-pi) q[1];
rz(-1.5775852) q[2];
sx q[2];
rz(-1.3464084) q[2];
sx q[2];
rz(-2.9780088) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.89167833) q[1];
sx q[1];
rz(-2.7418156) q[1];
sx q[1];
rz(2.1348025) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9336003) q[3];
sx q[3];
rz(-1.1054512) q[3];
sx q[3];
rz(0.057010827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91743177) q[2];
sx q[2];
rz(-1.5956343) q[2];
sx q[2];
rz(1.4707627) q[2];
rz(-1.5077695) q[3];
sx q[3];
rz(-0.016851146) q[3];
sx q[3];
rz(-2.2182218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5690174) q[0];
sx q[0];
rz(-1.9439531) q[0];
sx q[0];
rz(-1.5684599) q[0];
rz(2.9720427) q[1];
sx q[1];
rz(-0.1145656) q[1];
sx q[1];
rz(3.0068908) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3189663) q[0];
sx q[0];
rz(-1.9958226) q[0];
sx q[0];
rz(-1.951773) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5021654) q[2];
sx q[2];
rz(-1.5830212) q[2];
sx q[2];
rz(3.1032012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.70054322) q[1];
sx q[1];
rz(-1.8722539) q[1];
sx q[1];
rz(0.45977199) q[1];
rz(1.7487583) q[3];
sx q[3];
rz(-1.6497532) q[3];
sx q[3];
rz(0.7782026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0160825) q[2];
sx q[2];
rz(-1.6566015) q[2];
sx q[2];
rz(0.15277319) q[2];
rz(1.7759391) q[3];
sx q[3];
rz(-0.036866166) q[3];
sx q[3];
rz(2.9735145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.633054) q[0];
sx q[0];
rz(-0.80798739) q[0];
sx q[0];
rz(-0.48164865) q[0];
rz(-2.956849) q[1];
sx q[1];
rz(-1.7748723) q[1];
sx q[1];
rz(2.1462323) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48993123) q[0];
sx q[0];
rz(-1.7999987) q[0];
sx q[0];
rz(1.2152753) q[0];
rz(1.6312509) q[2];
sx q[2];
rz(-1.5222856) q[2];
sx q[2];
rz(-1.2918351) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0936443) q[1];
sx q[1];
rz(-1.8694832) q[1];
sx q[1];
rz(0.68409749) q[1];
x q[2];
rz(1.5257201) q[3];
sx q[3];
rz(-2.5451676) q[3];
sx q[3];
rz(-3.0395122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2153726) q[2];
sx q[2];
rz(-0.054457713) q[2];
sx q[2];
rz(0.064662956) q[2];
rz(2.0926545) q[3];
sx q[3];
rz(-0.026853042) q[3];
sx q[3];
rz(-1.8612727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30508405) q[0];
sx q[0];
rz(-3.0236112) q[0];
sx q[0];
rz(2.3205561) q[0];
rz(-3.0631284) q[1];
sx q[1];
rz(-1.4469701) q[1];
sx q[1];
rz(2.1441114) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4560735) q[0];
sx q[0];
rz(-2.4302097) q[0];
sx q[0];
rz(2.6472241) q[0];
x q[1];
rz(1.6299963) q[2];
sx q[2];
rz(-1.6069876) q[2];
sx q[2];
rz(2.9179887) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8197937) q[1];
sx q[1];
rz(-2.5410278) q[1];
sx q[1];
rz(-1.0933881) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6281582) q[3];
sx q[3];
rz(-2.0485176) q[3];
sx q[3];
rz(2.3582332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0492101) q[2];
sx q[2];
rz(-3.1123078) q[2];
sx q[2];
rz(-1.1537665) q[2];
rz(2.9554101) q[3];
sx q[3];
rz(-3.0598873) q[3];
sx q[3];
rz(2.8103099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1368197) q[0];
sx q[0];
rz(-2.3240219) q[0];
sx q[0];
rz(1.9983043) q[0];
rz(1.1413057) q[1];
sx q[1];
rz(-2.3316796) q[1];
sx q[1];
rz(2.6236261) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2835613) q[0];
sx q[0];
rz(-1.0378583) q[0];
sx q[0];
rz(-2.3951247) q[0];
rz(-pi) q[1];
rz(-1.5557655) q[2];
sx q[2];
rz(-1.5717197) q[2];
sx q[2];
rz(-2.8944588) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8570366) q[1];
sx q[1];
rz(-1.9003881) q[1];
sx q[1];
rz(-0.28917851) q[1];
x q[2];
rz(0.2980026) q[3];
sx q[3];
rz(-1.7116705) q[3];
sx q[3];
rz(-0.44671392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8225857) q[2];
sx q[2];
rz(-0.95390445) q[2];
sx q[2];
rz(2.8138568) q[2];
rz(1.1945126) q[3];
sx q[3];
rz(-0.12735282) q[3];
sx q[3];
rz(2.2849042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1389393) q[0];
sx q[0];
rz(-2.8728573) q[0];
sx q[0];
rz(2.5797381) q[0];
rz(1.4909164) q[1];
sx q[1];
rz(-1.5166538) q[1];
sx q[1];
rz(-3.0432826) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5033103) q[0];
sx q[0];
rz(-1.9925963) q[0];
sx q[0];
rz(-3.1027334) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1409522) q[2];
sx q[2];
rz(-1.5709236) q[2];
sx q[2];
rz(1.887111) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.505762) q[1];
sx q[1];
rz(-1.1019076) q[1];
sx q[1];
rz(-1.6230349) q[1];
x q[2];
rz(-2.0344072) q[3];
sx q[3];
rz(-2.0258198) q[3];
sx q[3];
rz(2.1912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0564698) q[2];
sx q[2];
rz(-0.19530185) q[2];
sx q[2];
rz(-1.0657715) q[2];
rz(-0.39984518) q[3];
sx q[3];
rz(-2.6078434) q[3];
sx q[3];
rz(1.8736418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14168508) q[0];
sx q[0];
rz(-0.081900224) q[0];
sx q[0];
rz(-1.6837233) q[0];
rz(-2.0211925) q[1];
sx q[1];
rz(-3.0034062) q[1];
sx q[1];
rz(0.33946005) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65050478) q[0];
sx q[0];
rz(-1.5840713) q[0];
sx q[0];
rz(-0.076939452) q[0];
x q[1];
rz(0.55468126) q[2];
sx q[2];
rz(-1.5590057) q[2];
sx q[2];
rz(1.5761216) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.016638674) q[1];
sx q[1];
rz(-1.5724239) q[1];
sx q[1];
rz(1.6563936) q[1];
rz(-1.6704329) q[3];
sx q[3];
rz(-1.2335586) q[3];
sx q[3];
rz(1.8189614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.3642984) q[2];
sx q[2];
rz(-3.1396781) q[2];
sx q[2];
rz(-2.778229) q[2];
rz(2.0511138) q[3];
sx q[3];
rz(-2.5616779) q[3];
sx q[3];
rz(2.0232078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5217487) q[0];
sx q[0];
rz(-2.813297) q[0];
sx q[0];
rz(-0.92292619) q[0];
rz(-1.6587616) q[1];
sx q[1];
rz(-0.62983477) q[1];
sx q[1];
rz(3.1373851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18989604) q[0];
sx q[0];
rz(-2.0810662) q[0];
sx q[0];
rz(2.4416591) q[0];
rz(-1.1599067) q[2];
sx q[2];
rz(-1.5762323) q[2];
sx q[2];
rz(1.5559352) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.002122238) q[1];
sx q[1];
rz(-1.0957901) q[1];
sx q[1];
rz(3.0691181) q[1];
rz(-pi) q[2];
rz(-3.0757853) q[3];
sx q[3];
rz(-2.5281457) q[3];
sx q[3];
rz(2.891401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5620455) q[2];
sx q[2];
rz(-1.6041218) q[2];
sx q[2];
rz(-1.9407678) q[2];
rz(-1.3935401) q[3];
sx q[3];
rz(-0.0037007185) q[3];
sx q[3];
rz(-2.4414731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30018184) q[0];
sx q[0];
rz(-2.6926079) q[0];
sx q[0];
rz(-1.9218943) q[0];
rz(-1.8118743) q[1];
sx q[1];
rz(-2.0025573) q[1];
sx q[1];
rz(0.16389287) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9326646) q[0];
sx q[0];
rz(-2.744356) q[0];
sx q[0];
rz(-1.1139289) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58456011) q[2];
sx q[2];
rz(-1.341408) q[2];
sx q[2];
rz(3.0363415) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.039197415) q[1];
sx q[1];
rz(-1.6886687) q[1];
sx q[1];
rz(0.060023569) q[1];
rz(2.4285942) q[3];
sx q[3];
rz(-0.87942356) q[3];
sx q[3];
rz(2.612243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.71358744) q[2];
sx q[2];
rz(-0.089944936) q[2];
sx q[2];
rz(2.5186727) q[2];
rz(3.0981787) q[3];
sx q[3];
rz(-0.89689887) q[3];
sx q[3];
rz(2.3626732) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6503705) q[0];
sx q[0];
rz(-0.0064042052) q[0];
sx q[0];
rz(-0.48625913) q[0];
rz(2.4468415) q[1];
sx q[1];
rz(-0.34725747) q[1];
sx q[1];
rz(-2.6659226) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18494517) q[0];
sx q[0];
rz(-1.6066736) q[0];
sx q[0];
rz(0.033680276) q[0];
x q[1];
rz(0.034157201) q[2];
sx q[2];
rz(-2.2766621) q[2];
sx q[2];
rz(-1.5547084) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.039363843) q[1];
sx q[1];
rz(-0.96323955) q[1];
sx q[1];
rz(-1.5786912) q[1];
rz(-pi) q[2];
rz(0.31020152) q[3];
sx q[3];
rz(-1.7473012) q[3];
sx q[3];
rz(-1.6506529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67705578) q[2];
sx q[2];
rz(-3.0812283) q[2];
sx q[2];
rz(1.2976868) q[2];
rz(-1.248598) q[3];
sx q[3];
rz(-0.55515754) q[3];
sx q[3];
rz(-0.36778522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0155335) q[0];
sx q[0];
rz(-1.5801237) q[0];
sx q[0];
rz(1.7042241) q[0];
rz(-2.2592648) q[1];
sx q[1];
rz(-3.075141) q[1];
sx q[1];
rz(-2.3022423) q[1];
rz(-1.5346943) q[2];
sx q[2];
rz(-1.4786199) q[2];
sx q[2];
rz(0.27603966) q[2];
rz(1.5658548) q[3];
sx q[3];
rz(-1.8099804) q[3];
sx q[3];
rz(-3.119131) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

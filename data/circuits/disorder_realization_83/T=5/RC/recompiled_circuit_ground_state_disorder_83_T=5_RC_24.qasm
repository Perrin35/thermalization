OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2228425) q[0];
sx q[0];
rz(-1.9434682) q[0];
sx q[0];
rz(-0.54984468) q[0];
rz(-0.066601872) q[1];
sx q[1];
rz(5.0636518) q[1];
sx q[1];
rz(10.650462) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9860814) q[0];
sx q[0];
rz(-1.5442463) q[0];
sx q[0];
rz(-0.20602137) q[0];
rz(-pi) q[1];
rz(-3.0751905) q[2];
sx q[2];
rz(-1.4208671) q[2];
sx q[2];
rz(0.87519803) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5388214) q[1];
sx q[1];
rz(-0.32646561) q[1];
sx q[1];
rz(-1.5646255) q[1];
rz(-2.5494895) q[3];
sx q[3];
rz(-1.0463011) q[3];
sx q[3];
rz(1.8703574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77039069) q[2];
sx q[2];
rz(-1.0545701) q[2];
sx q[2];
rz(0.87898177) q[2];
rz(-0.061035872) q[3];
sx q[3];
rz(-1.6998484) q[3];
sx q[3];
rz(0.92476168) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9691129) q[0];
sx q[0];
rz(-1.1107439) q[0];
sx q[0];
rz(1.9222395) q[0];
rz(-0.99766937) q[1];
sx q[1];
rz(-1.4439293) q[1];
sx q[1];
rz(-2.3692621) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0270841) q[0];
sx q[0];
rz(-2.3128384) q[0];
sx q[0];
rz(2.7161612) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6222811) q[2];
sx q[2];
rz(-0.33308187) q[2];
sx q[2];
rz(2.5903828) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5520305) q[1];
sx q[1];
rz(-2.8300474) q[1];
sx q[1];
rz(-3.0208842) q[1];
rz(0.62735438) q[3];
sx q[3];
rz(-1.5511906) q[3];
sx q[3];
rz(-0.62999187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4274009) q[2];
sx q[2];
rz(-1.6855468) q[2];
sx q[2];
rz(-1.5217155) q[2];
rz(-1.4766258) q[3];
sx q[3];
rz(-1.0790389) q[3];
sx q[3];
rz(0.21292201) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6427226) q[0];
sx q[0];
rz(-1.1551789) q[0];
sx q[0];
rz(-2.2464519) q[0];
rz(0.25998947) q[1];
sx q[1];
rz(-1.3464059) q[1];
sx q[1];
rz(-2.3666429) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4473559) q[0];
sx q[0];
rz(-2.1292344) q[0];
sx q[0];
rz(1.1456699) q[0];
x q[1];
rz(-0.16037143) q[2];
sx q[2];
rz(-2.7087492) q[2];
sx q[2];
rz(2.3059887) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.65371014) q[1];
sx q[1];
rz(-1.0431494) q[1];
sx q[1];
rz(2.5608173) q[1];
rz(-pi) q[2];
rz(0.633274) q[3];
sx q[3];
rz(-0.97304854) q[3];
sx q[3];
rz(-2.5839295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8746752) q[2];
sx q[2];
rz(-0.83628925) q[2];
sx q[2];
rz(-0.8379035) q[2];
rz(-0.72495929) q[3];
sx q[3];
rz(-0.91491142) q[3];
sx q[3];
rz(2.8673577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.2570141) q[0];
sx q[0];
rz(-3.0138636) q[0];
sx q[0];
rz(1.9789486) q[0];
rz(-2.0024025) q[1];
sx q[1];
rz(-1.2290686) q[1];
sx q[1];
rz(-0.3831648) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59012198) q[0];
sx q[0];
rz(-1.5801727) q[0];
sx q[0];
rz(-2.7290415) q[0];
rz(-1.2960984) q[2];
sx q[2];
rz(-0.75680671) q[2];
sx q[2];
rz(1.5001378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.63190118) q[1];
sx q[1];
rz(-1.5576524) q[1];
sx q[1];
rz(0.34850328) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64997079) q[3];
sx q[3];
rz(-1.293217) q[3];
sx q[3];
rz(1.3560719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.42134103) q[2];
sx q[2];
rz(-2.4946419) q[2];
sx q[2];
rz(1.7163537) q[2];
rz(-2.0187812) q[3];
sx q[3];
rz(-0.69901005) q[3];
sx q[3];
rz(-0.56306806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4238033) q[0];
sx q[0];
rz(-0.20861067) q[0];
sx q[0];
rz(-0.31722379) q[0];
rz(-2.9605561) q[1];
sx q[1];
rz(-1.92417) q[1];
sx q[1];
rz(-0.6212298) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8375299) q[0];
sx q[0];
rz(-2.3637536) q[0];
sx q[0];
rz(-1.2647948) q[0];
x q[1];
rz(-1.6124837) q[2];
sx q[2];
rz(-1.0752077) q[2];
sx q[2];
rz(2.1245196) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.23074977) q[1];
sx q[1];
rz(-0.98346868) q[1];
sx q[1];
rz(3.0452646) q[1];
rz(0.43790169) q[3];
sx q[3];
rz(-2.2931778) q[3];
sx q[3];
rz(3.0423328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8165596) q[2];
sx q[2];
rz(-2.4123522) q[2];
sx q[2];
rz(2.0751591) q[2];
rz(2.1224497) q[3];
sx q[3];
rz(-2.416553) q[3];
sx q[3];
rz(1.2044005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2524734) q[0];
sx q[0];
rz(-2.7121565) q[0];
sx q[0];
rz(-0.47527894) q[0];
rz(0.94398445) q[1];
sx q[1];
rz(-1.9119268) q[1];
sx q[1];
rz(2.7164187) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19744148) q[0];
sx q[0];
rz(-1.7867435) q[0];
sx q[0];
rz(-0.39255377) q[0];
rz(-pi) q[1];
x q[1];
rz(1.556301) q[2];
sx q[2];
rz(-1.0767848) q[2];
sx q[2];
rz(1.838889) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7092998) q[1];
sx q[1];
rz(-1.4008351) q[1];
sx q[1];
rz(-1.0372838) q[1];
rz(-pi) q[2];
rz(-0.3495174) q[3];
sx q[3];
rz(-1.8807285) q[3];
sx q[3];
rz(-0.30924451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7298594) q[2];
sx q[2];
rz(-1.4889762) q[2];
sx q[2];
rz(-0.64424166) q[2];
rz(-1.980137) q[3];
sx q[3];
rz(-2.1803653) q[3];
sx q[3];
rz(-0.78288356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.2162868) q[0];
sx q[0];
rz(-1.1290978) q[0];
sx q[0];
rz(-2.570545) q[0];
rz(-1.1043999) q[1];
sx q[1];
rz(-2.0490502) q[1];
sx q[1];
rz(-0.33448321) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5307823) q[0];
sx q[0];
rz(-0.40662262) q[0];
sx q[0];
rz(1.6645891) q[0];
x q[1];
rz(1.1715631) q[2];
sx q[2];
rz(-2.3314771) q[2];
sx q[2];
rz(-0.1619815) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7240136) q[1];
sx q[1];
rz(-1.3162606) q[1];
sx q[1];
rz(1.2271787) q[1];
rz(-pi) q[2];
rz(-0.35558139) q[3];
sx q[3];
rz(-2.4556841) q[3];
sx q[3];
rz(2.4816328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0337446) q[2];
sx q[2];
rz(-1.1128384) q[2];
sx q[2];
rz(3.0372078) q[2];
rz(-2.8153822) q[3];
sx q[3];
rz(-2.2741337) q[3];
sx q[3];
rz(0.67702684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0831083) q[0];
sx q[0];
rz(-1.9309738) q[0];
sx q[0];
rz(2.6404358) q[0];
rz(1.1309364) q[1];
sx q[1];
rz(-0.94291818) q[1];
sx q[1];
rz(-1.2859734) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20984291) q[0];
sx q[0];
rz(-0.49389687) q[0];
sx q[0];
rz(-1.2487869) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3589972) q[2];
sx q[2];
rz(-2.411826) q[2];
sx q[2];
rz(1.9791918) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0803538) q[1];
sx q[1];
rz(-0.29268943) q[1];
sx q[1];
rz(2.5848542) q[1];
x q[2];
rz(0.85190947) q[3];
sx q[3];
rz(-0.63792568) q[3];
sx q[3];
rz(-1.7450116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8274294) q[2];
sx q[2];
rz(-2.8848727) q[2];
sx q[2];
rz(0.87744212) q[2];
rz(2.1432803) q[3];
sx q[3];
rz(-2.5613997) q[3];
sx q[3];
rz(-1.72054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5494635) q[0];
sx q[0];
rz(-1.9100459) q[0];
sx q[0];
rz(0.25729427) q[0];
rz(2.4843702) q[1];
sx q[1];
rz(-2.6723599) q[1];
sx q[1];
rz(1.8990272) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94958828) q[0];
sx q[0];
rz(-0.87118252) q[0];
sx q[0];
rz(-2.5618895) q[0];
rz(-pi) q[1];
rz(-1.7153344) q[2];
sx q[2];
rz(-0.78510127) q[2];
sx q[2];
rz(2.7018765) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5661652) q[1];
sx q[1];
rz(-1.0707483) q[1];
sx q[1];
rz(1.0779672) q[1];
rz(2.6493862) q[3];
sx q[3];
rz(-2.7653381) q[3];
sx q[3];
rz(0.87615651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.88951748) q[2];
sx q[2];
rz(-1.7097079) q[2];
sx q[2];
rz(0.70402181) q[2];
rz(1.4975123) q[3];
sx q[3];
rz(-1.5181395) q[3];
sx q[3];
rz(-1.2161072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4154476) q[0];
sx q[0];
rz(-2.4496267) q[0];
sx q[0];
rz(2.6494001) q[0];
rz(-0.65912143) q[1];
sx q[1];
rz(-0.4159795) q[1];
sx q[1];
rz(-0.98512828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27714849) q[0];
sx q[0];
rz(-1.5168958) q[0];
sx q[0];
rz(3.0646851) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3756204) q[2];
sx q[2];
rz(-1.6723126) q[2];
sx q[2];
rz(-0.40292172) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6669504) q[1];
sx q[1];
rz(-1.1792151) q[1];
sx q[1];
rz(1.8285059) q[1];
rz(-pi) q[2];
rz(3.1019347) q[3];
sx q[3];
rz(-0.83740202) q[3];
sx q[3];
rz(-1.7379675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27348406) q[2];
sx q[2];
rz(-2.2874139) q[2];
sx q[2];
rz(-1.9197397) q[2];
rz(-0.46011225) q[3];
sx q[3];
rz(-1.8729788) q[3];
sx q[3];
rz(2.4251895) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39501326) q[0];
sx q[0];
rz(-1.4656675) q[0];
sx q[0];
rz(-1.898484) q[0];
rz(1.5383491) q[1];
sx q[1];
rz(-1.0335045) q[1];
sx q[1];
rz(-0.43362591) q[1];
rz(-1.4547841) q[2];
sx q[2];
rz(-2.3064936) q[2];
sx q[2];
rz(2.4207122) q[2];
rz(0.015751377) q[3];
sx q[3];
rz(-1.2917662) q[3];
sx q[3];
rz(0.36242604) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(5.0603428) q[0];
sx q[0];
rz(4.3397171) q[0];
sx q[0];
rz(8.8749333) q[0];
rz(-0.066601872) q[1];
sx q[1];
rz(-1.2195335) q[1];
sx q[1];
rz(1.2256844) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.420833) q[0];
sx q[0];
rz(-1.3648486) q[0];
sx q[0];
rz(1.5979196) q[0];
rz(3.0751905) q[2];
sx q[2];
rz(-1.4208671) q[2];
sx q[2];
rz(-0.87519803) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5962564) q[1];
sx q[1];
rz(-1.2443372) q[1];
sx q[1];
rz(0.0020892987) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96192653) q[3];
sx q[3];
rz(-1.0666218) q[3];
sx q[3];
rz(0.025328764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.371202) q[2];
sx q[2];
rz(-2.0870225) q[2];
sx q[2];
rz(0.87898177) q[2];
rz(-3.0805568) q[3];
sx q[3];
rz(-1.6998484) q[3];
sx q[3];
rz(-0.92476168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9691129) q[0];
sx q[0];
rz(-1.1107439) q[0];
sx q[0];
rz(1.2193532) q[0];
rz(-0.99766937) q[1];
sx q[1];
rz(-1.6976633) q[1];
sx q[1];
rz(2.3692621) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7051429) q[0];
sx q[0];
rz(-2.3068743) q[0];
sx q[0];
rz(1.993773) q[0];
rz(1.2381239) q[2];
sx q[2];
rz(-1.5539697) q[2];
sx q[2];
rz(2.1706659) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0961842) q[1];
sx q[1];
rz(-1.533877) q[1];
sx q[1];
rz(0.30942076) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1082013) q[3];
sx q[3];
rz(-0.62761939) q[3];
sx q[3];
rz(-0.96783584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4274009) q[2];
sx q[2];
rz(-1.4560459) q[2];
sx q[2];
rz(1.6198772) q[2];
rz(-1.6649668) q[3];
sx q[3];
rz(-2.0625538) q[3];
sx q[3];
rz(0.21292201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6427226) q[0];
sx q[0];
rz(-1.9864137) q[0];
sx q[0];
rz(0.89514071) q[0];
rz(0.25998947) q[1];
sx q[1];
rz(-1.7951868) q[1];
sx q[1];
rz(2.3666429) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74028984) q[0];
sx q[0];
rz(-2.4537114) q[0];
sx q[0];
rz(2.5581261) q[0];
rz(1.4971447) q[2];
sx q[2];
rz(-1.9977187) q[2];
sx q[2];
rz(-2.1296453) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.59794489) q[1];
sx q[1];
rz(-1.0769118) q[1];
sx q[1];
rz(-2.1795401) q[1];
rz(-2.5083187) q[3];
sx q[3];
rz(-0.97304854) q[3];
sx q[3];
rz(0.55766314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8746752) q[2];
sx q[2];
rz(-0.83628925) q[2];
sx q[2];
rz(0.8379035) q[2];
rz(-0.72495929) q[3];
sx q[3];
rz(-0.91491142) q[3];
sx q[3];
rz(-0.27423492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.2570141) q[0];
sx q[0];
rz(-3.0138636) q[0];
sx q[0];
rz(-1.1626441) q[0];
rz(1.1391901) q[1];
sx q[1];
rz(-1.9125241) q[1];
sx q[1];
rz(-2.7584279) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1823381) q[0];
sx q[0];
rz(-2.7289411) q[0];
sx q[0];
rz(3.1182108) q[0];
rz(0.25077925) q[2];
sx q[2];
rz(-0.84883142) q[2];
sx q[2];
rz(2.0112558) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.93411968) q[1];
sx q[1];
rz(-1.9192682) q[1];
sx q[1];
rz(-1.5568118) q[1];
rz(-2.7015436) q[3];
sx q[3];
rz(-2.4428133) q[3];
sx q[3];
rz(-0.13126365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.42134103) q[2];
sx q[2];
rz(-0.64695078) q[2];
sx q[2];
rz(1.7163537) q[2];
rz(2.0187812) q[3];
sx q[3];
rz(-0.69901005) q[3];
sx q[3];
rz(-2.5785246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7177893) q[0];
sx q[0];
rz(-0.20861067) q[0];
sx q[0];
rz(-0.31722379) q[0];
rz(-0.18103655) q[1];
sx q[1];
rz(-1.92417) q[1];
sx q[1];
rz(0.6212298) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7214523) q[0];
sx q[0];
rz(-0.83774191) q[0];
sx q[0];
rz(2.8531403) q[0];
rz(-0.076940342) q[2];
sx q[2];
rz(-0.49719337) q[2];
sx q[2];
rz(-2.212008) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7382051) q[1];
sx q[1];
rz(-0.59425577) q[1];
sx q[1];
rz(-1.7142833) q[1];
x q[2];
rz(-2.3425668) q[3];
sx q[3];
rz(-1.247041) q[3];
sx q[3];
rz(1.369839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.3250331) q[2];
sx q[2];
rz(-0.72924048) q[2];
sx q[2];
rz(2.0751591) q[2];
rz(-2.1224497) q[3];
sx q[3];
rz(-2.416553) q[3];
sx q[3];
rz(-1.2044005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8891193) q[0];
sx q[0];
rz(-0.42943615) q[0];
sx q[0];
rz(-0.47527894) q[0];
rz(-0.94398445) q[1];
sx q[1];
rz(-1.2296659) q[1];
sx q[1];
rz(-0.42517391) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89579534) q[0];
sx q[0];
rz(-2.6962792) q[0];
sx q[0];
rz(2.6209339) q[0];
x q[1];
rz(-0.026907909) q[2];
sx q[2];
rz(-0.49420658) q[2];
sx q[2];
rz(-1.869452) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2380772) q[1];
sx q[1];
rz(-1.0457731) q[1];
sx q[1];
rz(0.19673853) q[1];
x q[2];
rz(-1.2422855) q[3];
sx q[3];
rz(-1.2385912) q[3];
sx q[3];
rz(-1.1508416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7298594) q[2];
sx q[2];
rz(-1.4889762) q[2];
sx q[2];
rz(0.64424166) q[2];
rz(-1.1614557) q[3];
sx q[3];
rz(-2.1803653) q[3];
sx q[3];
rz(-2.3587091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2162868) q[0];
sx q[0];
rz(-2.0124948) q[0];
sx q[0];
rz(2.570545) q[0];
rz(-2.0371927) q[1];
sx q[1];
rz(-2.0490502) q[1];
sx q[1];
rz(0.33448321) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5307823) q[0];
sx q[0];
rz(-0.40662262) q[0];
sx q[0];
rz(-1.4770035) q[0];
rz(2.7538497) q[2];
sx q[2];
rz(-2.3015112) q[2];
sx q[2];
rz(-2.4304488) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7664288) q[1];
sx q[1];
rz(-2.7169982) q[1];
sx q[1];
rz(0.91318513) q[1];
rz(-1.2932106) q[3];
sx q[3];
rz(-0.93507877) q[3];
sx q[3];
rz(1.1073974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.107848) q[2];
sx q[2];
rz(-1.1128384) q[2];
sx q[2];
rz(-0.10438485) q[2];
rz(-0.32621041) q[3];
sx q[3];
rz(-0.86745894) q[3];
sx q[3];
rz(0.67702684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0584843) q[0];
sx q[0];
rz(-1.9309738) q[0];
sx q[0];
rz(-0.50115681) q[0];
rz(-1.1309364) q[1];
sx q[1];
rz(-0.94291818) q[1];
sx q[1];
rz(-1.8556192) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5720309) q[0];
sx q[0];
rz(-2.0372197) q[0];
sx q[0];
rz(2.9728229) q[0];
x q[1];
rz(-2.4444163) q[2];
sx q[2];
rz(-1.3343658) q[2];
sx q[2];
rz(-0.6811179) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0612388) q[1];
sx q[1];
rz(-2.8489032) q[1];
sx q[1];
rz(-2.5848542) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0619265) q[3];
sx q[3];
rz(-1.1677907) q[3];
sx q[3];
rz(2.7030735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8274294) q[2];
sx q[2];
rz(-2.8848727) q[2];
sx q[2];
rz(0.87744212) q[2];
rz(0.99831239) q[3];
sx q[3];
rz(-2.5613997) q[3];
sx q[3];
rz(1.72054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5921291) q[0];
sx q[0];
rz(-1.2315467) q[0];
sx q[0];
rz(-0.25729427) q[0];
rz(-2.4843702) q[1];
sx q[1];
rz(-2.6723599) q[1];
sx q[1];
rz(-1.8990272) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9193756) q[0];
sx q[0];
rz(-1.1383445) q[0];
sx q[0];
rz(-2.3591756) q[0];
x q[1];
rz(-0.79093604) q[2];
sx q[2];
rz(-1.6727912) q[2];
sx q[2];
rz(-1.9079218) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4081289) q[1];
sx q[1];
rz(-2.4545547) q[1];
sx q[1];
rz(-0.71367674) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3862196) q[3];
sx q[3];
rz(-1.2410302) q[3];
sx q[3];
rz(0.3531526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2520752) q[2];
sx q[2];
rz(-1.7097079) q[2];
sx q[2];
rz(-0.70402181) q[2];
rz(-1.6440803) q[3];
sx q[3];
rz(-1.6234532) q[3];
sx q[3];
rz(1.2161072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7261451) q[0];
sx q[0];
rz(-0.69196597) q[0];
sx q[0];
rz(0.49219254) q[0];
rz(-0.65912143) q[1];
sx q[1];
rz(-0.4159795) q[1];
sx q[1];
rz(-0.98512828) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8644442) q[0];
sx q[0];
rz(-1.6246968) q[0];
sx q[0];
rz(3.0646851) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8707377) q[2];
sx q[2];
rz(-0.388467) q[2];
sx q[2];
rz(-2.2252639) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6669504) q[1];
sx q[1];
rz(-1.9623775) q[1];
sx q[1];
rz(1.3130867) q[1];
x q[2];
rz(-2.3045818) q[3];
sx q[3];
rz(-1.5413377) q[3];
sx q[3];
rz(0.1406167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27348406) q[2];
sx q[2];
rz(-0.85417875) q[2];
sx q[2];
rz(-1.2218529) q[2];
rz(-2.6814804) q[3];
sx q[3];
rz(-1.8729788) q[3];
sx q[3];
rz(-2.4251895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7465794) q[0];
sx q[0];
rz(-1.4656675) q[0];
sx q[0];
rz(-1.898484) q[0];
rz(-1.5383491) q[1];
sx q[1];
rz(-2.1080882) q[1];
sx q[1];
rz(2.7079667) q[1];
rz(-1.6868085) q[2];
sx q[2];
rz(-0.83509904) q[2];
sx q[2];
rz(-0.72088045) q[2];
rz(1.8498593) q[3];
sx q[3];
rz(-1.5859384) q[3];
sx q[3];
rz(1.9288837) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

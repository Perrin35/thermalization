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
rz(3.0749908) q[1];
sx q[1];
rz(-1.9220592) q[1];
sx q[1];
rz(1.9159082) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8526575) q[0];
sx q[0];
rz(-2.9338917) q[0];
sx q[0];
rz(3.0124979) q[0];
rz(-pi) q[1];
rz(3.0751905) q[2];
sx q[2];
rz(-1.7207256) q[2];
sx q[2];
rz(-2.2663946) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5388214) q[1];
sx q[1];
rz(-0.32646561) q[1];
sx q[1];
rz(1.5646255) q[1];
rz(-pi) q[2];
rz(0.59210316) q[3];
sx q[3];
rz(-2.0952916) q[3];
sx q[3];
rz(-1.8703574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.371202) q[2];
sx q[2];
rz(-1.0545701) q[2];
sx q[2];
rz(2.2626109) q[2];
rz(-3.0805568) q[3];
sx q[3];
rz(-1.6998484) q[3];
sx q[3];
rz(-0.92476168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9691129) q[0];
sx q[0];
rz(-2.0308487) q[0];
sx q[0];
rz(-1.2193532) q[0];
rz(0.99766937) q[1];
sx q[1];
rz(-1.6976633) q[1];
sx q[1];
rz(-2.3692621) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1145086) q[0];
sx q[0];
rz(-0.82875427) q[0];
sx q[0];
rz(-2.7161612) q[0];
rz(1.6222811) q[2];
sx q[2];
rz(-0.33308187) q[2];
sx q[2];
rz(-0.55120984) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.46281262) q[1];
sx q[1];
rz(-1.2615934) q[1];
sx q[1];
rz(1.5320381) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1082013) q[3];
sx q[3];
rz(-0.62761939) q[3];
sx q[3];
rz(-2.1737568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.71419174) q[2];
sx q[2];
rz(-1.6855468) q[2];
sx q[2];
rz(-1.6198772) q[2];
rz(1.6649668) q[3];
sx q[3];
rz(-1.0790389) q[3];
sx q[3];
rz(-2.9286706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6427226) q[0];
sx q[0];
rz(-1.1551789) q[0];
sx q[0];
rz(2.2464519) q[0];
rz(2.8816032) q[1];
sx q[1];
rz(-1.3464059) q[1];
sx q[1];
rz(2.3666429) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7827137) q[0];
sx q[0];
rz(-1.9281328) q[0];
sx q[0];
rz(0.60114791) q[0];
x q[1];
rz(1.644448) q[2];
sx q[2];
rz(-1.143874) q[2];
sx q[2];
rz(1.0119473) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59794489) q[1];
sx q[1];
rz(-2.0646808) q[1];
sx q[1];
rz(0.96205254) q[1];
x q[2];
rz(-0.85525671) q[3];
sx q[3];
rz(-0.84153131) q[3];
sx q[3];
rz(-2.7824912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8746752) q[2];
sx q[2];
rz(-0.83628925) q[2];
sx q[2];
rz(0.8379035) q[2];
rz(-2.4166334) q[3];
sx q[3];
rz(-2.2266812) q[3];
sx q[3];
rz(2.8673577) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2570141) q[0];
sx q[0];
rz(-3.0138636) q[0];
sx q[0];
rz(1.9789486) q[0];
rz(1.1391901) q[1];
sx q[1];
rz(-1.2290686) q[1];
sx q[1];
rz(-0.3831648) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1823381) q[0];
sx q[0];
rz(-2.7289411) q[0];
sx q[0];
rz(0.023381845) q[0];
rz(-pi) q[1];
rz(2.3085528) q[2];
sx q[2];
rz(-1.7581356) q[2];
sx q[2];
rz(0.27275547) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.207473) q[1];
sx q[1];
rz(-1.2223244) q[1];
sx q[1];
rz(1.5847809) q[1];
rz(1.9145033) q[3];
sx q[3];
rz(-0.94961221) q[3];
sx q[3];
rz(-2.7215001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7202516) q[2];
sx q[2];
rz(-2.4946419) q[2];
sx q[2];
rz(-1.7163537) q[2];
rz(-1.1228115) q[3];
sx q[3];
rz(-0.69901005) q[3];
sx q[3];
rz(0.56306806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.4238033) q[0];
sx q[0];
rz(-2.932982) q[0];
sx q[0];
rz(-0.31722379) q[0];
rz(-2.9605561) q[1];
sx q[1];
rz(-1.92417) q[1];
sx q[1];
rz(-0.6212298) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7214523) q[0];
sx q[0];
rz(-0.83774191) q[0];
sx q[0];
rz(0.28845235) q[0];
x q[1];
rz(-1.5291089) q[2];
sx q[2];
rz(-2.0663849) q[2];
sx q[2];
rz(-1.017073) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40338755) q[1];
sx q[1];
rz(-0.59425577) q[1];
sx q[1];
rz(-1.4273093) q[1];
rz(-2.703691) q[3];
sx q[3];
rz(-0.84841484) q[3];
sx q[3];
rz(-3.0423328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8165596) q[2];
sx q[2];
rz(-2.4123522) q[2];
sx q[2];
rz(-2.0751591) q[2];
rz(-2.1224497) q[3];
sx q[3];
rz(-2.416553) q[3];
sx q[3];
rz(1.9371921) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8891193) q[0];
sx q[0];
rz(-0.42943615) q[0];
sx q[0];
rz(0.47527894) q[0];
rz(0.94398445) q[1];
sx q[1];
rz(-1.2296659) q[1];
sx q[1];
rz(-2.7164187) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9441512) q[0];
sx q[0];
rz(-1.3548491) q[0];
sx q[0];
rz(0.39255377) q[0];
rz(-pi) q[1];
rz(0.026907909) q[2];
sx q[2];
rz(-0.49420658) q[2];
sx q[2];
rz(-1.2721407) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2819967) q[1];
sx q[1];
rz(-2.5841671) q[1];
sx q[1];
rz(1.2453399) q[1];
x q[2];
rz(-2.7920753) q[3];
sx q[3];
rz(-1.2608642) q[3];
sx q[3];
rz(-0.30924451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7298594) q[2];
sx q[2];
rz(-1.4889762) q[2];
sx q[2];
rz(-2.497351) q[2];
rz(-1.980137) q[3];
sx q[3];
rz(-2.1803653) q[3];
sx q[3];
rz(2.3587091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92530584) q[0];
sx q[0];
rz(-2.0124948) q[0];
sx q[0];
rz(-0.5710477) q[0];
rz(-2.0371927) q[1];
sx q[1];
rz(-1.0925424) q[1];
sx q[1];
rz(2.8071094) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0153941) q[0];
sx q[0];
rz(-1.6078464) q[0];
sx q[0];
rz(-1.1657715) q[0];
rz(2.3399722) q[2];
sx q[2];
rz(-1.2853664) q[2];
sx q[2];
rz(1.1257671) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4175791) q[1];
sx q[1];
rz(-1.8253321) q[1];
sx q[1];
rz(-1.9144139) q[1];
rz(-pi) q[2];
rz(-1.848382) q[3];
sx q[3];
rz(-0.93507877) q[3];
sx q[3];
rz(-1.1073974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0337446) q[2];
sx q[2];
rz(-2.0287543) q[2];
sx q[2];
rz(-3.0372078) q[2];
rz(0.32621041) q[3];
sx q[3];
rz(-0.86745894) q[3];
sx q[3];
rz(2.4645658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0584843) q[0];
sx q[0];
rz(-1.9309738) q[0];
sx q[0];
rz(0.50115681) q[0];
rz(2.0106563) q[1];
sx q[1];
rz(-2.1986745) q[1];
sx q[1];
rz(-1.2859734) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5720309) q[0];
sx q[0];
rz(-1.104373) q[0];
sx q[0];
rz(0.16876975) q[0];
x q[1];
rz(-1.8752918) q[2];
sx q[2];
rz(-2.2448953) q[2];
sx q[2];
rz(-0.69597352) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0612388) q[1];
sx q[1];
rz(-2.8489032) q[1];
sx q[1];
rz(0.55673843) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6874316) q[3];
sx q[3];
rz(-1.1060904) q[3];
sx q[3];
rz(0.91688076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3141632) q[2];
sx q[2];
rz(-0.25671998) q[2];
sx q[2];
rz(-2.2641505) q[2];
rz(-2.1432803) q[3];
sx q[3];
rz(-2.5613997) q[3];
sx q[3];
rz(-1.4210526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5494635) q[0];
sx q[0];
rz(-1.2315467) q[0];
sx q[0];
rz(2.8842984) q[0];
rz(-2.4843702) q[1];
sx q[1];
rz(-0.46923271) q[1];
sx q[1];
rz(-1.2425655) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1920044) q[0];
sx q[0];
rz(-2.2704101) q[0];
sx q[0];
rz(-2.5618895) q[0];
rz(-pi) q[1];
rz(1.7153344) q[2];
sx q[2];
rz(-2.3564914) q[2];
sx q[2];
rz(-0.43971616) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8942106) q[1];
sx q[1];
rz(-1.9989493) q[1];
sx q[1];
rz(2.5864771) q[1];
rz(2.8065422) q[3];
sx q[3];
rz(-1.7453225) q[3];
sx q[3];
rz(-1.1572591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2520752) q[2];
sx q[2];
rz(-1.7097079) q[2];
sx q[2];
rz(0.70402181) q[2];
rz(1.4975123) q[3];
sx q[3];
rz(-1.6234532) q[3];
sx q[3];
rz(1.2161072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4154476) q[0];
sx q[0];
rz(-0.69196597) q[0];
sx q[0];
rz(-2.6494001) q[0];
rz(-2.4824712) q[1];
sx q[1];
rz(-2.7256131) q[1];
sx q[1];
rz(2.1564644) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2894963) q[0];
sx q[0];
rz(-1.4940007) q[0];
sx q[0];
rz(-1.6248563) q[0];
x q[1];
rz(1.4617301) q[2];
sx q[2];
rz(-1.1972053) q[2];
sx q[2];
rz(-1.2078169) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4746423) q[1];
sx q[1];
rz(-1.9623775) q[1];
sx q[1];
rz(1.8285059) q[1];
x q[2];
rz(-3.1019347) q[3];
sx q[3];
rz(-0.83740202) q[3];
sx q[3];
rz(1.7379675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8681086) q[2];
sx q[2];
rz(-2.2874139) q[2];
sx q[2];
rz(-1.2218529) q[2];
rz(0.46011225) q[3];
sx q[3];
rz(-1.8729788) q[3];
sx q[3];
rz(0.71640316) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7465794) q[0];
sx q[0];
rz(-1.6759251) q[0];
sx q[0];
rz(1.2431086) q[0];
rz(-1.5383491) q[1];
sx q[1];
rz(-2.1080882) q[1];
sx q[1];
rz(2.7079667) q[1];
rz(1.6868085) q[2];
sx q[2];
rz(-2.3064936) q[2];
sx q[2];
rz(2.4207122) q[2];
rz(-3.1258413) q[3];
sx q[3];
rz(-1.2917662) q[3];
sx q[3];
rz(0.36242604) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

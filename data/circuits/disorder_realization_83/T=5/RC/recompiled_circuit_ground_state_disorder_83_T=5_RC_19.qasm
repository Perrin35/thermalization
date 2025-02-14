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
rz(-1.2256844) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9860814) q[0];
sx q[0];
rz(-1.5442463) q[0];
sx q[0];
rz(0.20602137) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7210518) q[2];
sx q[2];
rz(-1.5051402) q[2];
sx q[2];
rz(0.70553095) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5453363) q[1];
sx q[1];
rz(-1.2443372) q[1];
sx q[1];
rz(-3.1395034) q[1];
x q[2];
rz(-2.5494895) q[3];
sx q[3];
rz(-1.0463011) q[3];
sx q[3];
rz(1.8703574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.371202) q[2];
sx q[2];
rz(-1.0545701) q[2];
sx q[2];
rz(-2.2626109) q[2];
rz(3.0805568) q[3];
sx q[3];
rz(-1.6998484) q[3];
sx q[3];
rz(0.92476168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1724797) q[0];
sx q[0];
rz(-1.1107439) q[0];
sx q[0];
rz(1.9222395) q[0];
rz(-0.99766937) q[1];
sx q[1];
rz(-1.6976633) q[1];
sx q[1];
rz(2.3692621) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0270841) q[0];
sx q[0];
rz(-0.82875427) q[0];
sx q[0];
rz(2.7161612) q[0];
rz(-pi) q[1];
rz(-1.2381239) q[2];
sx q[2];
rz(-1.587623) q[2];
sx q[2];
rz(-0.97092674) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0961842) q[1];
sx q[1];
rz(-1.6077157) q[1];
sx q[1];
rz(2.8321719) q[1];
rz(-pi) q[2];
rz(1.5950117) q[3];
sx q[3];
rz(-2.1980114) q[3];
sx q[3];
rz(-2.2150018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.71419174) q[2];
sx q[2];
rz(-1.4560459) q[2];
sx q[2];
rz(1.6198772) q[2];
rz(1.6649668) q[3];
sx q[3];
rz(-2.0625538) q[3];
sx q[3];
rz(2.9286706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6427226) q[0];
sx q[0];
rz(-1.9864137) q[0];
sx q[0];
rz(2.2464519) q[0];
rz(2.8816032) q[1];
sx q[1];
rz(-1.7951868) q[1];
sx q[1];
rz(-2.3666429) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74028984) q[0];
sx q[0];
rz(-2.4537114) q[0];
sx q[0];
rz(2.5581261) q[0];
x q[1];
rz(-0.16037143) q[2];
sx q[2];
rz(-2.7087492) q[2];
sx q[2];
rz(-0.835604) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4878825) q[1];
sx q[1];
rz(-1.0431494) q[1];
sx q[1];
rz(2.5608173) q[1];
x q[2];
rz(2.5083187) q[3];
sx q[3];
rz(-0.97304854) q[3];
sx q[3];
rz(-0.55766314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.2669175) q[2];
sx q[2];
rz(-0.83628925) q[2];
sx q[2];
rz(2.3036892) q[2];
rz(-2.4166334) q[3];
sx q[3];
rz(-2.2266812) q[3];
sx q[3];
rz(-0.27423492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2570141) q[0];
sx q[0];
rz(-0.12772904) q[0];
sx q[0];
rz(-1.1626441) q[0];
rz(2.0024025) q[1];
sx q[1];
rz(-1.2290686) q[1];
sx q[1];
rz(-2.7584279) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95925453) q[0];
sx q[0];
rz(-2.7289411) q[0];
sx q[0];
rz(-0.023381845) q[0];
x q[1];
rz(2.3085528) q[2];
sx q[2];
rz(-1.383457) q[2];
sx q[2];
rz(2.8688372) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.97505442) q[1];
sx q[1];
rz(-0.34874094) q[1];
sx q[1];
rz(0.038473126) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2270893) q[3];
sx q[3];
rz(-2.1919804) q[3];
sx q[3];
rz(0.4200926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42134103) q[2];
sx q[2];
rz(-2.4946419) q[2];
sx q[2];
rz(-1.7163537) q[2];
rz(-1.1228115) q[3];
sx q[3];
rz(-2.4425826) q[3];
sx q[3];
rz(2.5785246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4238033) q[0];
sx q[0];
rz(-2.932982) q[0];
sx q[0];
rz(0.31722379) q[0];
rz(2.9605561) q[1];
sx q[1];
rz(-1.92417) q[1];
sx q[1];
rz(0.6212298) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4201403) q[0];
sx q[0];
rz(-0.83774191) q[0];
sx q[0];
rz(0.28845235) q[0];
rz(-pi) q[1];
rz(1.5291089) q[2];
sx q[2];
rz(-1.0752077) q[2];
sx q[2];
rz(-1.017073) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2865528) q[1];
sx q[1];
rz(-1.650944) q[1];
sx q[1];
rz(0.98132332) q[1];
rz(-pi) q[2];
rz(-0.79902584) q[3];
sx q[3];
rz(-1.8945516) q[3];
sx q[3];
rz(-1.7717537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8165596) q[2];
sx q[2];
rz(-0.72924048) q[2];
sx q[2];
rz(2.0751591) q[2];
rz(1.019143) q[3];
sx q[3];
rz(-2.416553) q[3];
sx q[3];
rz(1.9371921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2524734) q[0];
sx q[0];
rz(-0.42943615) q[0];
sx q[0];
rz(-0.47527894) q[0];
rz(2.1976082) q[1];
sx q[1];
rz(-1.9119268) q[1];
sx q[1];
rz(-2.7164187) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4618414) q[0];
sx q[0];
rz(-1.9537524) q[0];
sx q[0];
rz(1.337685) q[0];
rz(-pi) q[1];
rz(0.4940554) q[2];
sx q[2];
rz(-1.5835585) q[2];
sx q[2];
rz(-0.27496613) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.43229285) q[1];
sx q[1];
rz(-1.4008351) q[1];
sx q[1];
rz(-2.1043088) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7920753) q[3];
sx q[3];
rz(-1.8807285) q[3];
sx q[3];
rz(-0.30924451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7298594) q[2];
sx q[2];
rz(-1.6526165) q[2];
sx q[2];
rz(-2.497351) q[2];
rz(-1.980137) q[3];
sx q[3];
rz(-0.96122733) q[3];
sx q[3];
rz(0.78288356) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
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
rz(0.5710477) q[0];
rz(2.0371927) q[1];
sx q[1];
rz(-1.0925424) q[1];
sx q[1];
rz(-2.8071094) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7128744) q[0];
sx q[0];
rz(-1.1660657) q[0];
sx q[0];
rz(-3.1012845) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38774298) q[2];
sx q[2];
rz(-2.3015112) q[2];
sx q[2];
rz(2.4304488) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4175791) q[1];
sx q[1];
rz(-1.8253321) q[1];
sx q[1];
rz(1.2271787) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2932106) q[3];
sx q[3];
rz(-0.93507877) q[3];
sx q[3];
rz(-1.1073974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.107848) q[2];
sx q[2];
rz(-2.0287543) q[2];
sx q[2];
rz(-0.10438485) q[2];
rz(-2.8153822) q[3];
sx q[3];
rz(-2.2741337) q[3];
sx q[3];
rz(-2.4645658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0584843) q[0];
sx q[0];
rz(-1.9309738) q[0];
sx q[0];
rz(-0.50115681) q[0];
rz(2.0106563) q[1];
sx q[1];
rz(-0.94291818) q[1];
sx q[1];
rz(-1.8556192) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9317497) q[0];
sx q[0];
rz(-2.6476958) q[0];
sx q[0];
rz(-1.2487869) q[0];
x q[1];
rz(-0.69717631) q[2];
sx q[2];
rz(-1.8072268) q[2];
sx q[2];
rz(2.4604748) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.046989881) q[1];
sx q[1];
rz(-1.7238574) q[1];
sx q[1];
rz(-2.8911289) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.454161) q[3];
sx q[3];
rz(-1.1060904) q[3];
sx q[3];
rz(0.91688076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8274294) q[2];
sx q[2];
rz(-2.8848727) q[2];
sx q[2];
rz(0.87744212) q[2];
rz(0.99831239) q[3];
sx q[3];
rz(-2.5613997) q[3];
sx q[3];
rz(-1.4210526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5494635) q[0];
sx q[0];
rz(-1.9100459) q[0];
sx q[0];
rz(2.8842984) q[0];
rz(-2.4843702) q[1];
sx q[1];
rz(-2.6723599) q[1];
sx q[1];
rz(-1.8990272) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22221702) q[0];
sx q[0];
rz(-2.0032481) q[0];
sx q[0];
rz(2.3591756) q[0];
rz(0.79093604) q[2];
sx q[2];
rz(-1.4688014) q[2];
sx q[2];
rz(-1.9079218) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.73346371) q[1];
sx q[1];
rz(-2.4545547) q[1];
sx q[1];
rz(-2.4279159) q[1];
x q[2];
rz(0.33505043) q[3];
sx q[3];
rz(-1.3962702) q[3];
sx q[3];
rz(-1.1572591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.88951748) q[2];
sx q[2];
rz(-1.4318848) q[2];
sx q[2];
rz(0.70402181) q[2];
rz(1.4975123) q[3];
sx q[3];
rz(-1.6234532) q[3];
sx q[3];
rz(-1.9254855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7261451) q[0];
sx q[0];
rz(-2.4496267) q[0];
sx q[0];
rz(-0.49219254) q[0];
rz(-2.4824712) q[1];
sx q[1];
rz(-2.7256131) q[1];
sx q[1];
rz(-0.98512828) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27714849) q[0];
sx q[0];
rz(-1.6246968) q[0];
sx q[0];
rz(-3.0646851) q[0];
rz(-pi) q[1];
rz(-2.8707377) q[2];
sx q[2];
rz(-2.7531257) q[2];
sx q[2];
rz(0.91632878) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6669504) q[1];
sx q[1];
rz(-1.9623775) q[1];
sx q[1];
rz(-1.3130867) q[1];
rz(2.3045818) q[3];
sx q[3];
rz(-1.5413377) q[3];
sx q[3];
rz(3.000976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8681086) q[2];
sx q[2];
rz(-0.85417875) q[2];
sx q[2];
rz(-1.2218529) q[2];
rz(2.6814804) q[3];
sx q[3];
rz(-1.8729788) q[3];
sx q[3];
rz(-0.71640316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7465794) q[0];
sx q[0];
rz(-1.6759251) q[0];
sx q[0];
rz(1.2431086) q[0];
rz(1.6032435) q[1];
sx q[1];
rz(-2.1080882) q[1];
sx q[1];
rz(2.7079667) q[1];
rz(2.4025386) q[2];
sx q[2];
rz(-1.6567163) q[2];
sx q[2];
rz(-2.2136282) q[2];
rz(3.1258413) q[3];
sx q[3];
rz(-1.8498265) q[3];
sx q[3];
rz(-2.7791666) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

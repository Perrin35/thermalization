OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.76685846) q[0];
sx q[0];
rz(-0.91349608) q[0];
sx q[0];
rz(1.7537533) q[0];
rz(0.18985441) q[1];
sx q[1];
rz(4.4603577) q[1];
sx q[1];
rz(9.7108884) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7095735) q[0];
sx q[0];
rz(-0.46105024) q[0];
sx q[0];
rz(-0.053857164) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0593023) q[2];
sx q[2];
rz(-0.9613885) q[2];
sx q[2];
rz(-1.0179969) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3158947) q[1];
sx q[1];
rz(-1.8301536) q[1];
sx q[1];
rz(-2.5729235) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0864618) q[3];
sx q[3];
rz(-1.1305439) q[3];
sx q[3];
rz(2.3326479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0949969) q[2];
sx q[2];
rz(-1.4904212) q[2];
sx q[2];
rz(0.68520927) q[2];
rz(1.4677216) q[3];
sx q[3];
rz(-1.0823715) q[3];
sx q[3];
rz(2.8732324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0005223) q[0];
sx q[0];
rz(-1.594161) q[0];
sx q[0];
rz(-2.1970314) q[0];
rz(-3.0682849) q[1];
sx q[1];
rz(-0.83275515) q[1];
sx q[1];
rz(-1.2578957) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1441378) q[0];
sx q[0];
rz(-1.0187444) q[0];
sx q[0];
rz(2.6004418) q[0];
rz(-pi) q[1];
rz(0.33690764) q[2];
sx q[2];
rz(-0.55977976) q[2];
sx q[2];
rz(2.155456) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91852616) q[1];
sx q[1];
rz(-2.7246973) q[1];
sx q[1];
rz(-2.7781163) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89711758) q[3];
sx q[3];
rz(-0.059487933) q[3];
sx q[3];
rz(-0.18732303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8153317) q[2];
sx q[2];
rz(-1.7504642) q[2];
sx q[2];
rz(2.1168671) q[2];
rz(2.4552086) q[3];
sx q[3];
rz(-2.4095583) q[3];
sx q[3];
rz(-2.2288403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.82655418) q[0];
sx q[0];
rz(-1.2737561) q[0];
sx q[0];
rz(2.3231373) q[0];
rz(1.2089027) q[1];
sx q[1];
rz(-1.5070288) q[1];
sx q[1];
rz(2.2817629) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61469528) q[0];
sx q[0];
rz(-0.69577587) q[0];
sx q[0];
rz(-2.1012563) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6099714) q[2];
sx q[2];
rz(-0.79382703) q[2];
sx q[2];
rz(-0.17228157) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2261994) q[1];
sx q[1];
rz(-2.2186154) q[1];
sx q[1];
rz(-0.97572216) q[1];
rz(-1.1016261) q[3];
sx q[3];
rz(-2.8896324) q[3];
sx q[3];
rz(0.89706883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7439421) q[2];
sx q[2];
rz(-0.33471477) q[2];
sx q[2];
rz(0.115455) q[2];
rz(0.3913106) q[3];
sx q[3];
rz(-1.0262998) q[3];
sx q[3];
rz(-1.9976043) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7855969) q[0];
sx q[0];
rz(-1.9556671) q[0];
sx q[0];
rz(2.9010229) q[0];
rz(-1.7754414) q[1];
sx q[1];
rz(-1.4931449) q[1];
sx q[1];
rz(-1.2496525) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8842953) q[0];
sx q[0];
rz(-2.4438071) q[0];
sx q[0];
rz(1.9724413) q[0];
rz(-2.7848148) q[2];
sx q[2];
rz(-0.98519221) q[2];
sx q[2];
rz(-2.9139589) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27423672) q[1];
sx q[1];
rz(-2.4411267) q[1];
sx q[1];
rz(2.7718616) q[1];
x q[2];
rz(0.68393882) q[3];
sx q[3];
rz(-1.1653644) q[3];
sx q[3];
rz(2.8059174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.25632855) q[2];
sx q[2];
rz(-0.54577959) q[2];
sx q[2];
rz(1.3679999) q[2];
rz(2.8374425) q[3];
sx q[3];
rz(-1.5757898) q[3];
sx q[3];
rz(-0.69653851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2279385) q[0];
sx q[0];
rz(-2.9253687) q[0];
sx q[0];
rz(-0.53713334) q[0];
rz(-2.3917603) q[1];
sx q[1];
rz(-2.0593819) q[1];
sx q[1];
rz(-2.4044663) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43057399) q[0];
sx q[0];
rz(-1.6665742) q[0];
sx q[0];
rz(-1.6349313) q[0];
rz(-1.7456878) q[2];
sx q[2];
rz(-1.6564257) q[2];
sx q[2];
rz(1.081664) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29038844) q[1];
sx q[1];
rz(-1.3121288) q[1];
sx q[1];
rz(1.438739) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26925663) q[3];
sx q[3];
rz(-0.33447166) q[3];
sx q[3];
rz(1.0448898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6970814) q[2];
sx q[2];
rz(-0.35327521) q[2];
sx q[2];
rz(-0.14661655) q[2];
rz(1.4491436) q[3];
sx q[3];
rz(-1.861898) q[3];
sx q[3];
rz(-0.52391887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.523282) q[0];
sx q[0];
rz(-1.016541) q[0];
sx q[0];
rz(-1.4215533) q[0];
rz(1.2591741) q[1];
sx q[1];
rz(-2.4651395) q[1];
sx q[1];
rz(0.4037942) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97966563) q[0];
sx q[0];
rz(-2.296519) q[0];
sx q[0];
rz(2.7620074) q[0];
x q[1];
rz(2.5228259) q[2];
sx q[2];
rz(-0.31224373) q[2];
sx q[2];
rz(-2.00756) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48133367) q[1];
sx q[1];
rz(-1.0254045) q[1];
sx q[1];
rz(2.27687) q[1];
rz(-1.8122702) q[3];
sx q[3];
rz(-0.89697402) q[3];
sx q[3];
rz(-0.50844565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8518565) q[2];
sx q[2];
rz(-1.7278262) q[2];
sx q[2];
rz(-2.9750321) q[2];
rz(1.5038495) q[3];
sx q[3];
rz(-0.47347355) q[3];
sx q[3];
rz(0.09854266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.903776) q[0];
sx q[0];
rz(-2.7523968) q[0];
sx q[0];
rz(-2.26407) q[0];
rz(-2.1794043) q[1];
sx q[1];
rz(-2.018237) q[1];
sx q[1];
rz(-2.7872564) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7276132) q[0];
sx q[0];
rz(-0.2383735) q[0];
sx q[0];
rz(0.075043126) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.118204) q[2];
sx q[2];
rz(-1.4260458) q[2];
sx q[2];
rz(-2.5770066) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6207032) q[1];
sx q[1];
rz(-0.9714618) q[1];
sx q[1];
rz(-0.66461892) q[1];
rz(-pi) q[2];
rz(-1.6463424) q[3];
sx q[3];
rz(-2.0795238) q[3];
sx q[3];
rz(-1.4010324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.80498901) q[2];
sx q[2];
rz(-0.73770928) q[2];
sx q[2];
rz(1.0450411) q[2];
rz(2.7392144) q[3];
sx q[3];
rz(-1.1674403) q[3];
sx q[3];
rz(-1.6654525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43110111) q[0];
sx q[0];
rz(-3.0248108) q[0];
sx q[0];
rz(0.85103273) q[0];
rz(0.35375133) q[1];
sx q[1];
rz(-1.4101135) q[1];
sx q[1];
rz(0.85465777) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25682043) q[0];
sx q[0];
rz(-1.2899634) q[0];
sx q[0];
rz(1.7600842) q[0];
rz(-pi) q[1];
rz(1.4209916) q[2];
sx q[2];
rz(-1.8179107) q[2];
sx q[2];
rz(2.832156) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.51616299) q[1];
sx q[1];
rz(-0.99402797) q[1];
sx q[1];
rz(1.7712084) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1885957) q[3];
sx q[3];
rz(-2.7321987) q[3];
sx q[3];
rz(-0.90367095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19055584) q[2];
sx q[2];
rz(-2.4943116) q[2];
sx q[2];
rz(0.62038842) q[2];
rz(-0.22647151) q[3];
sx q[3];
rz(-1.3969996) q[3];
sx q[3];
rz(-0.58102077) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9613551) q[0];
sx q[0];
rz(-1.9397475) q[0];
sx q[0];
rz(-0.015137976) q[0];
rz(2.119428) q[1];
sx q[1];
rz(-0.41003761) q[1];
sx q[1];
rz(1.2000363) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2101636) q[0];
sx q[0];
rz(-2.4151122) q[0];
sx q[0];
rz(1.8379712) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5483918) q[2];
sx q[2];
rz(-1.6199379) q[2];
sx q[2];
rz(-1.0047439) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.92654681) q[1];
sx q[1];
rz(-2.6290942) q[1];
sx q[1];
rz(1.5096751) q[1];
rz(-pi) q[2];
rz(-1.0608114) q[3];
sx q[3];
rz(-1.1378764) q[3];
sx q[3];
rz(2.1149389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.52593652) q[2];
sx q[2];
rz(-2.2286712) q[2];
sx q[2];
rz(-0.23942648) q[2];
rz(-0.058852363) q[3];
sx q[3];
rz(-2.3299496) q[3];
sx q[3];
rz(-2.8656901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9489768) q[0];
sx q[0];
rz(-1.2139576) q[0];
sx q[0];
rz(-0.97491997) q[0];
rz(2.4661567) q[1];
sx q[1];
rz(-2.2223739) q[1];
sx q[1];
rz(2.1598037) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47768053) q[0];
sx q[0];
rz(-0.64723158) q[0];
sx q[0];
rz(-2.0114698) q[0];
rz(-pi) q[1];
rz(-1.6390332) q[2];
sx q[2];
rz(-1.1600798) q[2];
sx q[2];
rz(1.0194743) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8999381) q[1];
sx q[1];
rz(-1.3119427) q[1];
sx q[1];
rz(2.3090825) q[1];
rz(1.6493787) q[3];
sx q[3];
rz(-2.693697) q[3];
sx q[3];
rz(0.092002951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5741817) q[2];
sx q[2];
rz(-2.1705706) q[2];
sx q[2];
rz(1.9166454) q[2];
rz(0.36568668) q[3];
sx q[3];
rz(-0.94958011) q[3];
sx q[3];
rz(1.1398116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0443307) q[0];
sx q[0];
rz(-1.2428357) q[0];
sx q[0];
rz(-1.6999929) q[0];
rz(3.0615831) q[1];
sx q[1];
rz(-1.7623822) q[1];
sx q[1];
rz(-0.0026127876) q[1];
rz(1.710113) q[2];
sx q[2];
rz(-1.1291383) q[2];
sx q[2];
rz(0.51869803) q[2];
rz(0.45251485) q[3];
sx q[3];
rz(-1.5502676) q[3];
sx q[3];
rz(1.4905139) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

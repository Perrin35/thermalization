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
rz(-2.4177457) q[0];
sx q[0];
rz(-1.9404193) q[0];
sx q[0];
rz(-0.49402753) q[0];
rz(1.5552893) q[1];
sx q[1];
rz(-0.92702335) q[1];
sx q[1];
rz(0.85049367) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5828398) q[0];
sx q[0];
rz(-0.90769288) q[0];
sx q[0];
rz(2.9302674) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57588864) q[2];
sx q[2];
rz(-2.401899) q[2];
sx q[2];
rz(1.0583371) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9697992) q[1];
sx q[1];
rz(-1.7928908) q[1];
sx q[1];
rz(-0.33116688) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1876202) q[3];
sx q[3];
rz(-0.37196649) q[3];
sx q[3];
rz(0.99727977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4832619) q[2];
sx q[2];
rz(-0.61507812) q[2];
sx q[2];
rz(-0.94493803) q[2];
rz(-0.0065217892) q[3];
sx q[3];
rz(-2.3801453) q[3];
sx q[3];
rz(2.884088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1061123) q[0];
sx q[0];
rz(-2.1529614) q[0];
sx q[0];
rz(2.535787) q[0];
rz(2.2094191) q[1];
sx q[1];
rz(-1.4235539) q[1];
sx q[1];
rz(0.1056284) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86565844) q[0];
sx q[0];
rz(-1.7248885) q[0];
sx q[0];
rz(0.6886512) q[0];
x q[1];
rz(-0.35657652) q[2];
sx q[2];
rz(-2.0322213) q[2];
sx q[2];
rz(-0.8666641) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.16766891) q[1];
sx q[1];
rz(-0.55096704) q[1];
sx q[1];
rz(-0.58631222) q[1];
x q[2];
rz(-0.37704269) q[3];
sx q[3];
rz(-1.1867513) q[3];
sx q[3];
rz(2.3689601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2848844) q[2];
sx q[2];
rz(-1.7785037) q[2];
sx q[2];
rz(-1.6178097) q[2];
rz(2.3273322) q[3];
sx q[3];
rz(-1.7819449) q[3];
sx q[3];
rz(0.8849357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098323671) q[0];
sx q[0];
rz(-0.79541484) q[0];
sx q[0];
rz(-1.5650308) q[0];
rz(2.1463429) q[1];
sx q[1];
rz(-2.1627656) q[1];
sx q[1];
rz(-0.78688041) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3580619) q[0];
sx q[0];
rz(-1.3710877) q[0];
sx q[0];
rz(-3.003503) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22898211) q[2];
sx q[2];
rz(-0.35576421) q[2];
sx q[2];
rz(-0.34504978) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39549024) q[1];
sx q[1];
rz(-1.8999892) q[1];
sx q[1];
rz(-1.2389061) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0819503) q[3];
sx q[3];
rz(-1.5111877) q[3];
sx q[3];
rz(-2.8622421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3028822) q[2];
sx q[2];
rz(-1.4883214) q[2];
sx q[2];
rz(-2.5035456) q[2];
rz(-2.6436515) q[3];
sx q[3];
rz(-1.0718071) q[3];
sx q[3];
rz(2.0518484) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8968673) q[0];
sx q[0];
rz(-1.3571955) q[0];
sx q[0];
rz(-3.0778399) q[0];
rz(-1.6925192) q[1];
sx q[1];
rz(-1.8359102) q[1];
sx q[1];
rz(1.8316899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5638014) q[0];
sx q[0];
rz(-3.0534857) q[0];
sx q[0];
rz(-2.7035575) q[0];
rz(-1.3856349) q[2];
sx q[2];
rz(-1.7660391) q[2];
sx q[2];
rz(0.38530063) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4938072) q[1];
sx q[1];
rz(-1.4469379) q[1];
sx q[1];
rz(0.74614831) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2579262) q[3];
sx q[3];
rz(-1.2901297) q[3];
sx q[3];
rz(-1.7629881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.070907585) q[2];
sx q[2];
rz(-2.0347774) q[2];
sx q[2];
rz(3.0359388) q[2];
rz(1.1834772) q[3];
sx q[3];
rz(-2.2453997) q[3];
sx q[3];
rz(0.67160523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79528177) q[0];
sx q[0];
rz(-0.91049894) q[0];
sx q[0];
rz(1.7972535) q[0];
rz(-0.67289871) q[1];
sx q[1];
rz(-1.3105323) q[1];
sx q[1];
rz(0.013462822) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1707662) q[0];
sx q[0];
rz(-1.544569) q[0];
sx q[0];
rz(0.88084014) q[0];
x q[1];
rz(1.5384664) q[2];
sx q[2];
rz(-2.7271977) q[2];
sx q[2];
rz(2.3737645) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9078482) q[1];
sx q[1];
rz(-1.1554759) q[1];
sx q[1];
rz(-1.585258) q[1];
x q[2];
rz(2.0125403) q[3];
sx q[3];
rz(-2.3886282) q[3];
sx q[3];
rz(-2.4747965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5992735) q[2];
sx q[2];
rz(-0.66212526) q[2];
sx q[2];
rz(1.2467965) q[2];
rz(0.072619297) q[3];
sx q[3];
rz(-0.19425546) q[3];
sx q[3];
rz(0.3092002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4319864) q[0];
sx q[0];
rz(-2.015634) q[0];
sx q[0];
rz(-3.0288938) q[0];
rz(-0.20507774) q[1];
sx q[1];
rz(-1.8094742) q[1];
sx q[1];
rz(0.90788666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8255709) q[0];
sx q[0];
rz(-2.1793723) q[0];
sx q[0];
rz(0.84644239) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61300437) q[2];
sx q[2];
rz(-1.2132436) q[2];
sx q[2];
rz(1.845971) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4389613) q[1];
sx q[1];
rz(-1.4794032) q[1];
sx q[1];
rz(-2.5083739) q[1];
rz(-1.7849633) q[3];
sx q[3];
rz(-0.91717824) q[3];
sx q[3];
rz(0.71398338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1418566) q[2];
sx q[2];
rz(-0.33679589) q[2];
sx q[2];
rz(-0.067961819) q[2];
rz(1.9939907) q[3];
sx q[3];
rz(-1.8736898) q[3];
sx q[3];
rz(1.2609153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1934018) q[0];
sx q[0];
rz(-1.8302487) q[0];
sx q[0];
rz(-0.46352682) q[0];
rz(2.9947128) q[1];
sx q[1];
rz(-1.2356707) q[1];
sx q[1];
rz(0.99259496) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50336058) q[0];
sx q[0];
rz(-1.8192756) q[0];
sx q[0];
rz(2.4186224) q[0];
rz(-pi) q[1];
rz(1.2847177) q[2];
sx q[2];
rz(-1.5062638) q[2];
sx q[2];
rz(2.4308824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9149962) q[1];
sx q[1];
rz(-2.2961756) q[1];
sx q[1];
rz(0.62208773) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73294183) q[3];
sx q[3];
rz(-1.9187201) q[3];
sx q[3];
rz(-0.80845736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2831882) q[2];
sx q[2];
rz(-1.4294727) q[2];
sx q[2];
rz(-0.44537133) q[2];
rz(1.8177659) q[3];
sx q[3];
rz(-2.0711074) q[3];
sx q[3];
rz(-2.0557192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.104798) q[0];
sx q[0];
rz(-3.1322271) q[0];
sx q[0];
rz(-0.15583663) q[0];
rz(1.2771295) q[1];
sx q[1];
rz(-1.3469478) q[1];
sx q[1];
rz(-3.1256622) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4919419) q[0];
sx q[0];
rz(-1.8109522) q[0];
sx q[0];
rz(2.35046) q[0];
x q[1];
rz(2.7819949) q[2];
sx q[2];
rz(-0.82530515) q[2];
sx q[2];
rz(-1.8458837) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4602892) q[1];
sx q[1];
rz(-2.2573514) q[1];
sx q[1];
rz(2.249447) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5032477) q[3];
sx q[3];
rz(-1.2223635) q[3];
sx q[3];
rz(-1.3909457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.62898034) q[2];
sx q[2];
rz(-0.11233687) q[2];
sx q[2];
rz(0.12574276) q[2];
rz(-2.1851152) q[3];
sx q[3];
rz(-1.7185017) q[3];
sx q[3];
rz(-2.8023348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15928957) q[0];
sx q[0];
rz(-0.24022261) q[0];
sx q[0];
rz(-0.45641986) q[0];
rz(1.5727111) q[1];
sx q[1];
rz(-1.0036889) q[1];
sx q[1];
rz(-0.68663866) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68701247) q[0];
sx q[0];
rz(-2.0622776) q[0];
sx q[0];
rz(-2.4647852) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1988611) q[2];
sx q[2];
rz(-0.21272993) q[2];
sx q[2];
rz(-0.63791927) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9514044) q[1];
sx q[1];
rz(-1.4238289) q[1];
sx q[1];
rz(1.4495871) q[1];
x q[2];
rz(-1.6410646) q[3];
sx q[3];
rz(-1.7450376) q[3];
sx q[3];
rz(-1.4605923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.75866428) q[2];
sx q[2];
rz(-1.1229346) q[2];
sx q[2];
rz(-1.135896) q[2];
rz(0.9797594) q[3];
sx q[3];
rz(-1.3330678) q[3];
sx q[3];
rz(-0.69825828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6887688) q[0];
sx q[0];
rz(-1.3249506) q[0];
sx q[0];
rz(-0.45595566) q[0];
rz(3.0988354) q[1];
sx q[1];
rz(-1.2084992) q[1];
sx q[1];
rz(2.4483689) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9589577) q[0];
sx q[0];
rz(-1.6696249) q[0];
sx q[0];
rz(1.1846428) q[0];
x q[1];
rz(1.1666388) q[2];
sx q[2];
rz(-1.1956788) q[2];
sx q[2];
rz(0.98304316) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.4148414) q[1];
sx q[1];
rz(-1.5778687) q[1];
sx q[1];
rz(1.4727791) q[1];
rz(-pi) q[2];
rz(1.8878804) q[3];
sx q[3];
rz(-1.1389995) q[3];
sx q[3];
rz(-0.82796873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2263055) q[2];
sx q[2];
rz(-1.8733571) q[2];
sx q[2];
rz(0.66429663) q[2];
rz(1.9281467) q[3];
sx q[3];
rz(-0.72187859) q[3];
sx q[3];
rz(0.87219316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.450347) q[0];
sx q[0];
rz(-2.3007614) q[0];
sx q[0];
rz(1.7578516) q[0];
rz(-1.4585523) q[1];
sx q[1];
rz(-2.8589307) q[1];
sx q[1];
rz(0.9160441) q[1];
rz(-0.37049313) q[2];
sx q[2];
rz(-1.5975614) q[2];
sx q[2];
rz(2.3536828) q[2];
rz(-2.4356859) q[3];
sx q[3];
rz(-2.0173895) q[3];
sx q[3];
rz(2.8593393) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

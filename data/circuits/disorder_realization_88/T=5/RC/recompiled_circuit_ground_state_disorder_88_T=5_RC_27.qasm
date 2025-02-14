OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4360566) q[0];
sx q[0];
rz(-0.85890618) q[0];
sx q[0];
rz(0.98281759) q[0];
rz(4.3217826) q[1];
sx q[1];
rz(5.5569841) q[1];
sx q[1];
rz(10.32294) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9010278) q[0];
sx q[0];
rz(-1.254558) q[0];
sx q[0];
rz(-1.680512) q[0];
rz(-0.54141374) q[2];
sx q[2];
rz(-1.8820418) q[2];
sx q[2];
rz(2.1212861) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8232728) q[1];
sx q[1];
rz(-1.987769) q[1];
sx q[1];
rz(-0.61064536) q[1];
x q[2];
rz(1.4709365) q[3];
sx q[3];
rz(-1.7152399) q[3];
sx q[3];
rz(0.58510548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1232221) q[2];
sx q[2];
rz(-1.3192588) q[2];
sx q[2];
rz(2.3940864) q[2];
rz(-1.8788762) q[3];
sx q[3];
rz(-1.5371753) q[3];
sx q[3];
rz(3.1178764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30775192) q[0];
sx q[0];
rz(-3.096014) q[0];
sx q[0];
rz(2.0763092) q[0];
rz(0.65159687) q[1];
sx q[1];
rz(-2.7862796) q[1];
sx q[1];
rz(-1.6337055) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7972854) q[0];
sx q[0];
rz(-2.4977614) q[0];
sx q[0];
rz(1.0539529) q[0];
rz(-2.5152757) q[2];
sx q[2];
rz(-2.5765315) q[2];
sx q[2];
rz(0.14742278) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3976685) q[1];
sx q[1];
rz(-2.5062291) q[1];
sx q[1];
rz(0.85224908) q[1];
x q[2];
rz(-1.5298858) q[3];
sx q[3];
rz(-2.5299151) q[3];
sx q[3];
rz(-0.024723147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4773341) q[2];
sx q[2];
rz(-1.4724255) q[2];
sx q[2];
rz(0.092546917) q[2];
rz(-2.7009098) q[3];
sx q[3];
rz(-2.7350072) q[3];
sx q[3];
rz(-1.996076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4785081) q[0];
sx q[0];
rz(-0.19190754) q[0];
sx q[0];
rz(1.9051911) q[0];
rz(-2.9036486) q[1];
sx q[1];
rz(-1.6704208) q[1];
sx q[1];
rz(-0.96756378) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0210812) q[0];
sx q[0];
rz(-0.97645611) q[0];
sx q[0];
rz(-0.8404151) q[0];
rz(2.601449) q[2];
sx q[2];
rz(-0.61068084) q[2];
sx q[2];
rz(1.3474828) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0261126) q[1];
sx q[1];
rz(-2.7819595) q[1];
sx q[1];
rz(-2.1277277) q[1];
x q[2];
rz(1.4400515) q[3];
sx q[3];
rz(-2.7249911) q[3];
sx q[3];
rz(2.1611913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48270109) q[2];
sx q[2];
rz(-0.37223688) q[2];
sx q[2];
rz(2.2230395) q[2];
rz(2.3679768) q[3];
sx q[3];
rz(-2.2548803) q[3];
sx q[3];
rz(-2.1913967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2660265) q[0];
sx q[0];
rz(-2.4022864) q[0];
sx q[0];
rz(-2.8493122) q[0];
rz(-2.267011) q[1];
sx q[1];
rz(-2.5129109) q[1];
sx q[1];
rz(-0.94327092) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2349512) q[0];
sx q[0];
rz(-0.72777339) q[0];
sx q[0];
rz(-1.1207188) q[0];
rz(-2.9246632) q[2];
sx q[2];
rz(-1.1116905) q[2];
sx q[2];
rz(-0.87068671) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0318289) q[1];
sx q[1];
rz(-2.2110143) q[1];
sx q[1];
rz(-0.39217453) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1499829) q[3];
sx q[3];
rz(-2.5217223) q[3];
sx q[3];
rz(0.068178328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5248519) q[2];
sx q[2];
rz(-2.0276232) q[2];
sx q[2];
rz(-2.2576766) q[2];
rz(1.135745) q[3];
sx q[3];
rz(-0.92848778) q[3];
sx q[3];
rz(-2.4511724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.32490548) q[0];
sx q[0];
rz(-0.08520928) q[0];
sx q[0];
rz(-2.9499522) q[0];
rz(1.2958255) q[1];
sx q[1];
rz(-1.192966) q[1];
sx q[1];
rz(0.4506909) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2795581) q[0];
sx q[0];
rz(-2.4992895) q[0];
sx q[0];
rz(-1.7234283) q[0];
rz(-pi) q[1];
rz(-0.24709267) q[2];
sx q[2];
rz(-0.28045248) q[2];
sx q[2];
rz(0.3874028) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1289336) q[1];
sx q[1];
rz(-1.4988572) q[1];
sx q[1];
rz(-2.1484445) q[1];
x q[2];
rz(-2.8487318) q[3];
sx q[3];
rz(-1.5080875) q[3];
sx q[3];
rz(0.69578275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.4638764) q[2];
sx q[2];
rz(-1.2748984) q[2];
sx q[2];
rz(-0.95787588) q[2];
rz(3.0887582) q[3];
sx q[3];
rz(-2.5962679) q[3];
sx q[3];
rz(0.42832819) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65474725) q[0];
sx q[0];
rz(-2.0992278) q[0];
sx q[0];
rz(-2.6690707) q[0];
rz(-0.76332244) q[1];
sx q[1];
rz(-0.7205874) q[1];
sx q[1];
rz(0.29019132) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71745738) q[0];
sx q[0];
rz(-0.73834921) q[0];
sx q[0];
rz(0.33070578) q[0];
rz(-pi) q[1];
rz(1.0792208) q[2];
sx q[2];
rz(-0.60814684) q[2];
sx q[2];
rz(-1.7660994) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5033333) q[1];
sx q[1];
rz(-2.4730254) q[1];
sx q[1];
rz(0.16600538) q[1];
rz(-pi) q[2];
rz(1.3715586) q[3];
sx q[3];
rz(-1.7459948) q[3];
sx q[3];
rz(0.74149473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7758238) q[2];
sx q[2];
rz(-2.2799728) q[2];
sx q[2];
rz(0.59930581) q[2];
rz(1.722909) q[3];
sx q[3];
rz(-1.8214046) q[3];
sx q[3];
rz(-2.1883709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80239427) q[0];
sx q[0];
rz(-1.9559487) q[0];
sx q[0];
rz(-1.3108569) q[0];
rz(-1.6400379) q[1];
sx q[1];
rz(-2.7412667) q[1];
sx q[1];
rz(0.60360533) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30234114) q[0];
sx q[0];
rz(-0.44284824) q[0];
sx q[0];
rz(-1.2542972) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3519751) q[2];
sx q[2];
rz(-2.1116774) q[2];
sx q[2];
rz(-0.70798028) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.43638602) q[1];
sx q[1];
rz(-1.5028186) q[1];
sx q[1];
rz(-1.8118565) q[1];
x q[2];
rz(2.4000136) q[3];
sx q[3];
rz(-1.7125579) q[3];
sx q[3];
rz(-0.89270702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32460585) q[2];
sx q[2];
rz(-2.2548455) q[2];
sx q[2];
rz(2.9465607) q[2];
rz(-2.4845691) q[3];
sx q[3];
rz(-1.7200836) q[3];
sx q[3];
rz(-3.0278964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3312155) q[0];
sx q[0];
rz(-1.4137784) q[0];
sx q[0];
rz(-0.072176607) q[0];
rz(0.78308925) q[1];
sx q[1];
rz(-2.5091645) q[1];
sx q[1];
rz(2.2236688) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1872301) q[0];
sx q[0];
rz(-1.3636806) q[0];
sx q[0];
rz(2.6595479) q[0];
rz(-pi) q[1];
x q[1];
rz(2.520944) q[2];
sx q[2];
rz(-1.3507412) q[2];
sx q[2];
rz(-1.4043851) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0705441) q[1];
sx q[1];
rz(-1.8080825) q[1];
sx q[1];
rz(-0.29233934) q[1];
rz(-pi) q[2];
rz(1.1456756) q[3];
sx q[3];
rz(-1.3367869) q[3];
sx q[3];
rz(-1.0581512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5808631) q[2];
sx q[2];
rz(-1.2141289) q[2];
sx q[2];
rz(1.6403991) q[2];
rz(-2.255127) q[3];
sx q[3];
rz(-1.3366046) q[3];
sx q[3];
rz(-0.10543536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60370541) q[0];
sx q[0];
rz(-2.7925346) q[0];
sx q[0];
rz(-2.6339997) q[0];
rz(0.72788584) q[1];
sx q[1];
rz(-1.4898841) q[1];
sx q[1];
rz(-1.1061888) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9637601) q[0];
sx q[0];
rz(-2.7179657) q[0];
sx q[0];
rz(-2.2209441) q[0];
rz(0.44194965) q[2];
sx q[2];
rz(-2.7486251) q[2];
sx q[2];
rz(3.0073187) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4708097) q[1];
sx q[1];
rz(-1.2074665) q[1];
sx q[1];
rz(-2.7881507) q[1];
x q[2];
rz(0.98420268) q[3];
sx q[3];
rz(-1.9846669) q[3];
sx q[3];
rz(-2.4936287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9375406) q[2];
sx q[2];
rz(-2.0378518) q[2];
sx q[2];
rz(-2.6943915) q[2];
rz(-1.4372545) q[3];
sx q[3];
rz(-2.3279326) q[3];
sx q[3];
rz(-1.4849496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60780418) q[0];
sx q[0];
rz(-2.0543126) q[0];
sx q[0];
rz(2.6265889) q[0];
rz(-1.1732514) q[1];
sx q[1];
rz(-1.2898022) q[1];
sx q[1];
rz(0.44949964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0093548) q[0];
sx q[0];
rz(-1.5482246) q[0];
sx q[0];
rz(1.2448553) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1859288) q[2];
sx q[2];
rz(-1.1329831) q[2];
sx q[2];
rz(0.73305886) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.044552) q[1];
sx q[1];
rz(-2.1637193) q[1];
sx q[1];
rz(1.0059936) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8490542) q[3];
sx q[3];
rz(-2.4026767) q[3];
sx q[3];
rz(-2.7310039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0476734) q[2];
sx q[2];
rz(-0.27457044) q[2];
sx q[2];
rz(0.74335113) q[2];
rz(-0.36176935) q[3];
sx q[3];
rz(-2.1105364) q[3];
sx q[3];
rz(0.37775347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.34306985) q[0];
sx q[0];
rz(-1.1953851) q[0];
sx q[0];
rz(0.60085798) q[0];
rz(-2.9877904) q[1];
sx q[1];
rz(-1.8228795) q[1];
sx q[1];
rz(0.22620329) q[1];
rz(-1.1719731) q[2];
sx q[2];
rz(-0.33179596) q[2];
sx q[2];
rz(-1.1992762) q[2];
rz(2.5220925) q[3];
sx q[3];
rz(-1.8826857) q[3];
sx q[3];
rz(2.7553325) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

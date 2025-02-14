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
rz(-2.7231554) q[0];
sx q[0];
rz(-2.1783481) q[0];
sx q[0];
rz(-0.20382398) q[0];
rz(-2.0889497) q[1];
sx q[1];
rz(-0.79336762) q[1];
sx q[1];
rz(-1.5348943) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8255768) q[0];
sx q[0];
rz(-1.1281779) q[0];
sx q[0];
rz(2.8809263) q[0];
rz(-pi) q[1];
rz(-2.7947172) q[2];
sx q[2];
rz(-0.87793186) q[2];
sx q[2];
rz(1.3077867) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1032627) q[1];
sx q[1];
rz(-0.94417773) q[1];
sx q[1];
rz(-2.7019397) q[1];
rz(2.5584869) q[3];
sx q[3];
rz(-0.2421549) q[3];
sx q[3];
rz(2.9271169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4172998) q[2];
sx q[2];
rz(-2.5827926) q[2];
sx q[2];
rz(2.144045) q[2];
rz(-0.7558465) q[3];
sx q[3];
rz(-1.7852802) q[3];
sx q[3];
rz(1.0409748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34870979) q[0];
sx q[0];
rz(-0.097276874) q[0];
sx q[0];
rz(-1.2874228) q[0];
rz(-2.6584794) q[1];
sx q[1];
rz(-2.099642) q[1];
sx q[1];
rz(-1.9226673) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12139509) q[0];
sx q[0];
rz(-1.4123865) q[0];
sx q[0];
rz(1.7686469) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8701209) q[2];
sx q[2];
rz(-0.69005943) q[2];
sx q[2];
rz(-0.59218237) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.375735) q[1];
sx q[1];
rz(-1.2925783) q[1];
sx q[1];
rz(-0.15032676) q[1];
rz(2.2628701) q[3];
sx q[3];
rz(-1.2936397) q[3];
sx q[3];
rz(2.4038278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7634742) q[2];
sx q[2];
rz(-3.0749622) q[2];
sx q[2];
rz(2.5197869) q[2];
rz(-2.5032737) q[3];
sx q[3];
rz(-0.74898762) q[3];
sx q[3];
rz(-2.2973255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8353552) q[0];
sx q[0];
rz(-2.4995646) q[0];
sx q[0];
rz(-2.9252692) q[0];
rz(-2.5430039) q[1];
sx q[1];
rz(-1.3052156) q[1];
sx q[1];
rz(-0.52282202) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3591839) q[0];
sx q[0];
rz(-2.6893927) q[0];
sx q[0];
rz(2.3308999) q[0];
x q[1];
rz(-1.4141358) q[2];
sx q[2];
rz(-1.9576503) q[2];
sx q[2];
rz(-2.000282) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8717958) q[1];
sx q[1];
rz(-1.7933473) q[1];
sx q[1];
rz(-2.6128164) q[1];
x q[2];
rz(1.2402291) q[3];
sx q[3];
rz(-1.2798314) q[3];
sx q[3];
rz(-1.7592848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9598976) q[2];
sx q[2];
rz(-1.2744224) q[2];
sx q[2];
rz(1.6240906) q[2];
rz(-2.6767139) q[3];
sx q[3];
rz(-2.2113694) q[3];
sx q[3];
rz(1.5215065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14438039) q[0];
sx q[0];
rz(-2.7535487) q[0];
sx q[0];
rz(-1.602518) q[0];
rz(-1.0333565) q[1];
sx q[1];
rz(-0.76834232) q[1];
sx q[1];
rz(1.8446406) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.813445) q[0];
sx q[0];
rz(-1.7792276) q[0];
sx q[0];
rz(0.84020241) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4960852) q[2];
sx q[2];
rz(-0.69715188) q[2];
sx q[2];
rz(0.32767228) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9819543) q[1];
sx q[1];
rz(-0.57483638) q[1];
sx q[1];
rz(2.3153618) q[1];
rz(-pi) q[2];
rz(0.4850895) q[3];
sx q[3];
rz(-0.61044932) q[3];
sx q[3];
rz(-1.4097253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79232717) q[2];
sx q[2];
rz(-2.5040864) q[2];
sx q[2];
rz(1.0456592) q[2];
rz(0.21444923) q[3];
sx q[3];
rz(-2.4172473) q[3];
sx q[3];
rz(-1.1469871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8056718) q[0];
sx q[0];
rz(-1.8593973) q[0];
sx q[0];
rz(-0.51934284) q[0];
rz(0.94201159) q[1];
sx q[1];
rz(-0.28128925) q[1];
sx q[1];
rz(1.6090144) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9056775) q[0];
sx q[0];
rz(-2.2122243) q[0];
sx q[0];
rz(1.1113313) q[0];
rz(-1.7132883) q[2];
sx q[2];
rz(-1.3918557) q[2];
sx q[2];
rz(-1.6662836) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6841507) q[1];
sx q[1];
rz(-1.5760211) q[1];
sx q[1];
rz(-2.0713776) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.665409) q[3];
sx q[3];
rz(-1.4615371) q[3];
sx q[3];
rz(-3.008568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64451009) q[2];
sx q[2];
rz(-0.77748674) q[2];
sx q[2];
rz(-0.28437781) q[2];
rz(2.583368) q[3];
sx q[3];
rz(-1.7773209) q[3];
sx q[3];
rz(0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.072642) q[0];
sx q[0];
rz(-0.41256368) q[0];
sx q[0];
rz(-2.5010342) q[0];
rz(1.6259646) q[1];
sx q[1];
rz(-1.1311572) q[1];
sx q[1];
rz(2.9277149) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2757077) q[0];
sx q[0];
rz(-1.8731706) q[0];
sx q[0];
rz(2.9024966) q[0];
rz(0.40753813) q[2];
sx q[2];
rz(-2.354051) q[2];
sx q[2];
rz(-2.119273) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1024903) q[1];
sx q[1];
rz(-2.5081345) q[1];
sx q[1];
rz(1.4536812) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6951896) q[3];
sx q[3];
rz(-2.6181707) q[3];
sx q[3];
rz(1.0833797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3512257) q[2];
sx q[2];
rz(-0.48206097) q[2];
sx q[2];
rz(-0.17769979) q[2];
rz(-1.3972345) q[3];
sx q[3];
rz(-1.9652365) q[3];
sx q[3];
rz(-2.7241404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099667065) q[0];
sx q[0];
rz(-2.5302027) q[0];
sx q[0];
rz(2.0360816) q[0];
rz(-0.28327495) q[1];
sx q[1];
rz(-0.53422821) q[1];
sx q[1];
rz(-1.6800605) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4257129) q[0];
sx q[0];
rz(-2.113803) q[0];
sx q[0];
rz(0.79820658) q[0];
rz(2.2093532) q[2];
sx q[2];
rz(-1.9912212) q[2];
sx q[2];
rz(2.1115542) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6532306) q[1];
sx q[1];
rz(-1.4547336) q[1];
sx q[1];
rz(-0.56452063) q[1];
x q[2];
rz(2.6631132) q[3];
sx q[3];
rz(-0.63204403) q[3];
sx q[3];
rz(1.5173591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0854411) q[2];
sx q[2];
rz(-1.5957007) q[2];
sx q[2];
rz(2.936787) q[2];
rz(-1.0497931) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(-2.7753196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7751223) q[0];
sx q[0];
rz(-2.6752052) q[0];
sx q[0];
rz(-3.108016) q[0];
rz(1.6665943) q[1];
sx q[1];
rz(-0.74445236) q[1];
sx q[1];
rz(-1.9974476) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0559674) q[0];
sx q[0];
rz(-0.48729839) q[0];
sx q[0];
rz(3.0874599) q[0];
rz(-0.099191908) q[2];
sx q[2];
rz(-1.8934403) q[2];
sx q[2];
rz(1.1485554) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.28539666) q[1];
sx q[1];
rz(-2.2335529) q[1];
sx q[1];
rz(-1.0876571) q[1];
rz(-pi) q[2];
rz(1.2253292) q[3];
sx q[3];
rz(-1.8356712) q[3];
sx q[3];
rz(0.027907413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7555776) q[2];
sx q[2];
rz(-2.3365946) q[2];
sx q[2];
rz(-1.8899274) q[2];
rz(1.3517316) q[3];
sx q[3];
rz(-0.91910619) q[3];
sx q[3];
rz(0.98021093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0018472483) q[0];
sx q[0];
rz(-2.503105) q[0];
sx q[0];
rz(-1.8852604) q[0];
rz(3.046335) q[1];
sx q[1];
rz(-0.88905159) q[1];
sx q[1];
rz(0.5307861) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7266709) q[0];
sx q[0];
rz(-1.4106531) q[0];
sx q[0];
rz(-1.3199575) q[0];
x q[1];
rz(2.8235648) q[2];
sx q[2];
rz(-1.5340759) q[2];
sx q[2];
rz(-1.6317692) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2797045) q[1];
sx q[1];
rz(-1.3315531) q[1];
sx q[1];
rz(-1.7025392) q[1];
rz(2.6952031) q[3];
sx q[3];
rz(-1.4671031) q[3];
sx q[3];
rz(-2.5959792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.329616) q[2];
sx q[2];
rz(-1.5831542) q[2];
sx q[2];
rz(-0.54083332) q[2];
rz(-0.81816188) q[3];
sx q[3];
rz(-1.4427789) q[3];
sx q[3];
rz(0.63280672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.624991) q[0];
sx q[0];
rz(-1.7778968) q[0];
sx q[0];
rz(-0.3987819) q[0];
rz(1.7575691) q[1];
sx q[1];
rz(-0.75474352) q[1];
sx q[1];
rz(-3.0922281) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20959768) q[0];
sx q[0];
rz(-1.5100129) q[0];
sx q[0];
rz(-1.6796735) q[0];
rz(0.10714679) q[2];
sx q[2];
rz(-1.309589) q[2];
sx q[2];
rz(3.0334453) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4409742) q[1];
sx q[1];
rz(-1.6153533) q[1];
sx q[1];
rz(-2.150321) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7743763) q[3];
sx q[3];
rz(-0.94881159) q[3];
sx q[3];
rz(-2.3443215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0177239) q[2];
sx q[2];
rz(-1.3088635) q[2];
sx q[2];
rz(1.6309942) q[2];
rz(-0.92783582) q[3];
sx q[3];
rz(-1.3414914) q[3];
sx q[3];
rz(-1.235289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8265726) q[0];
sx q[0];
rz(-2.6980504) q[0];
sx q[0];
rz(2.0499688) q[0];
rz(-2.0768968) q[1];
sx q[1];
rz(-2.6152492) q[1];
sx q[1];
rz(-1.1660887) q[1];
rz(2.3774316) q[2];
sx q[2];
rz(-2.4780826) q[2];
sx q[2];
rz(1.6586951) q[2];
rz(1.0743027) q[3];
sx q[3];
rz(-0.60582325) q[3];
sx q[3];
rz(-0.7614991) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

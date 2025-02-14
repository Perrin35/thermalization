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
rz(2.9377687) q[0];
rz(1.0526429) q[1];
sx q[1];
rz(3.9349603) q[1];
sx q[1];
rz(10.959672) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8255768) q[0];
sx q[0];
rz(-1.1281779) q[0];
sx q[0];
rz(-2.8809263) q[0];
rz(-pi) q[1];
rz(-2.7947172) q[2];
sx q[2];
rz(-0.87793186) q[2];
sx q[2];
rz(1.3077867) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.878256) q[1];
sx q[1];
rz(-1.2188101) q[1];
sx q[1];
rz(2.2455567) q[1];
rz(-0.20333692) q[3];
sx q[3];
rz(-1.4383738) q[3];
sx q[3];
rz(1.9258969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4172998) q[2];
sx q[2];
rz(-2.5827926) q[2];
sx q[2];
rz(2.144045) q[2];
rz(-0.7558465) q[3];
sx q[3];
rz(-1.3563124) q[3];
sx q[3];
rz(2.1006179) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34870979) q[0];
sx q[0];
rz(-0.097276874) q[0];
sx q[0];
rz(1.8541699) q[0];
rz(-0.48311326) q[1];
sx q[1];
rz(-2.099642) q[1];
sx q[1];
rz(1.9226673) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12139509) q[0];
sx q[0];
rz(-1.4123865) q[0];
sx q[0];
rz(-1.7686469) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67181113) q[2];
sx q[2];
rz(-1.7423358) q[2];
sx q[2];
rz(1.1900657) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.23663737) q[1];
sx q[1];
rz(-1.7153011) q[1];
sx q[1];
rz(-1.2895686) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2628701) q[3];
sx q[3];
rz(-1.2936397) q[3];
sx q[3];
rz(2.4038278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7634742) q[2];
sx q[2];
rz(-0.066630445) q[2];
sx q[2];
rz(2.5197869) q[2];
rz(0.63831896) q[3];
sx q[3];
rz(-2.392605) q[3];
sx q[3];
rz(-0.84426713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
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
rz(2.6187706) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7824088) q[0];
sx q[0];
rz(-0.45219996) q[0];
sx q[0];
rz(-2.3308999) q[0];
rz(-0.36575138) q[2];
sx q[2];
rz(-2.7257082) q[2];
sx q[2];
rz(-1.6037841) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8717958) q[1];
sx q[1];
rz(-1.7933473) q[1];
sx q[1];
rz(2.6128164) q[1];
rz(1.2402291) q[3];
sx q[3];
rz(-1.2798314) q[3];
sx q[3];
rz(1.3823079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9598976) q[2];
sx q[2];
rz(-1.8671702) q[2];
sx q[2];
rz(1.6240906) q[2];
rz(-2.6767139) q[3];
sx q[3];
rz(-0.93022323) q[3];
sx q[3];
rz(-1.5215065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14438039) q[0];
sx q[0];
rz(-0.388044) q[0];
sx q[0];
rz(-1.5390747) q[0];
rz(-2.1082361) q[1];
sx q[1];
rz(-2.3732503) q[1];
sx q[1];
rz(-1.296952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3281477) q[0];
sx q[0];
rz(-1.362365) q[0];
sx q[0];
rz(0.84020241) q[0];
rz(-pi) q[1];
rz(2.2665735) q[2];
sx q[2];
rz(-1.6187374) q[2];
sx q[2];
rz(1.9557916) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.32758157) q[1];
sx q[1];
rz(-1.1594698) q[1];
sx q[1];
rz(-2.7279305) q[1];
rz(-pi) q[2];
rz(2.5874073) q[3];
sx q[3];
rz(-1.300214) q[3];
sx q[3];
rz(-2.8949646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79232717) q[2];
sx q[2];
rz(-0.63750625) q[2];
sx q[2];
rz(1.0456592) q[2];
rz(2.9271434) q[3];
sx q[3];
rz(-0.72434536) q[3];
sx q[3];
rz(1.9946056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3359208) q[0];
sx q[0];
rz(-1.2821953) q[0];
sx q[0];
rz(-2.6222498) q[0];
rz(-2.1995811) q[1];
sx q[1];
rz(-2.8603034) q[1];
sx q[1];
rz(-1.6090144) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.518884) q[0];
sx q[0];
rz(-1.2075338) q[0];
sx q[0];
rz(2.4469482) q[0];
x q[1];
rz(-0.66560676) q[2];
sx q[2];
rz(-2.9133248) q[2];
sx q[2];
rz(2.1537202) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.457442) q[1];
sx q[1];
rz(-1.5760211) q[1];
sx q[1];
rz(-2.0713776) q[1];
rz(-pi) q[2];
rz(-3.0318465) q[3];
sx q[3];
rz(-1.4767495) q[3];
sx q[3];
rz(-1.4274244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4970826) q[2];
sx q[2];
rz(-0.77748674) q[2];
sx q[2];
rz(0.28437781) q[2];
rz(2.583368) q[3];
sx q[3];
rz(-1.7773209) q[3];
sx q[3];
rz(0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068950653) q[0];
sx q[0];
rz(-0.41256368) q[0];
sx q[0];
rz(0.64055842) q[0];
rz(1.6259646) q[1];
sx q[1];
rz(-2.0104355) q[1];
sx q[1];
rz(-2.9277149) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17992526) q[0];
sx q[0];
rz(-2.7583987) q[0];
sx q[0];
rz(-2.2201594) q[0];
rz(-pi) q[1];
rz(1.9496232) q[2];
sx q[2];
rz(-2.2791499) q[2];
sx q[2];
rz(1.5701936) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7678309) q[1];
sx q[1];
rz(-1.6400178) q[1];
sx q[1];
rz(0.94061416) q[1];
rz(-pi) q[2];
rz(1.446403) q[3];
sx q[3];
rz(-0.52342192) q[3];
sx q[3];
rz(1.0833797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3512257) q[2];
sx q[2];
rz(-0.48206097) q[2];
sx q[2];
rz(2.9638929) q[2];
rz(1.7443582) q[3];
sx q[3];
rz(-1.9652365) q[3];
sx q[3];
rz(0.41745225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099667065) q[0];
sx q[0];
rz(-2.5302027) q[0];
sx q[0];
rz(1.1055111) q[0];
rz(0.28327495) q[1];
sx q[1];
rz(-0.53422821) q[1];
sx q[1];
rz(-1.4615321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71587976) q[0];
sx q[0];
rz(-1.0277896) q[0];
sx q[0];
rz(2.3433861) q[0];
x q[1];
rz(-2.6335476) q[2];
sx q[2];
rz(-0.99544243) q[2];
sx q[2];
rz(-2.3066556) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9859559) q[1];
sx q[1];
rz(-1.0105304) q[1];
sx q[1];
rz(1.7079279) q[1];
rz(-pi) q[2];
rz(-2.6631132) q[3];
sx q[3];
rz(-0.63204403) q[3];
sx q[3];
rz(1.6242336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0561515) q[2];
sx q[2];
rz(-1.5458919) q[2];
sx q[2];
rz(-2.936787) q[2];
rz(-2.0917995) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(-0.36627305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3664704) q[0];
sx q[0];
rz(-2.6752052) q[0];
sx q[0];
rz(-3.108016) q[0];
rz(1.4749984) q[1];
sx q[1];
rz(-2.3971403) q[1];
sx q[1];
rz(1.144145) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9947203) q[0];
sx q[0];
rz(-2.0573186) q[0];
sx q[0];
rz(-1.5994607) q[0];
rz(-pi) q[1];
rz(3.0424007) q[2];
sx q[2];
rz(-1.8934403) q[2];
sx q[2];
rz(1.1485554) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.543964) q[1];
sx q[1];
rz(-1.9457327) q[1];
sx q[1];
rz(2.4191394) q[1];
rz(0.89544483) q[3];
sx q[3];
rz(-2.709528) q[3];
sx q[3];
rz(-0.91401446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7555776) q[2];
sx q[2];
rz(-2.3365946) q[2];
sx q[2];
rz(1.2516652) q[2];
rz(1.3517316) q[3];
sx q[3];
rz(-2.2224865) q[3];
sx q[3];
rz(2.1613817) q[3];
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
rz(-pi/2) q[0];
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
rz(-2.2525411) q[1];
sx q[1];
rz(-0.5307861) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7408001) q[0];
sx q[0];
rz(-0.2966899) q[0];
sx q[0];
rz(-2.147697) q[0];
rz(-0.1169493) q[2];
sx q[2];
rz(-2.8215234) q[2];
sx q[2];
rz(0.17203278) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.8618882) q[1];
sx q[1];
rz(-1.3315531) q[1];
sx q[1];
rz(-1.4390535) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6952031) q[3];
sx q[3];
rz(-1.4671031) q[3];
sx q[3];
rz(-0.54561347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.329616) q[2];
sx q[2];
rz(-1.5584385) q[2];
sx q[2];
rz(-0.54083332) q[2];
rz(0.81816188) q[3];
sx q[3];
rz(-1.6988138) q[3];
sx q[3];
rz(0.63280672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
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
rz(0.049364518) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.931995) q[0];
sx q[0];
rz(-1.6315797) q[0];
sx q[0];
rz(-1.4619191) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3081506) q[2];
sx q[2];
rz(-1.4672973) q[2];
sx q[2];
rz(1.4348794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2422707) q[1];
sx q[1];
rz(-0.99192109) q[1];
sx q[1];
rz(0.053236628) q[1];
rz(0.63189854) q[3];
sx q[3];
rz(-1.4057341) q[3];
sx q[3];
rz(0.89323211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0177239) q[2];
sx q[2];
rz(-1.8327291) q[2];
sx q[2];
rz(1.5105985) q[2];
rz(0.92783582) q[3];
sx q[3];
rz(-1.3414914) q[3];
sx q[3];
rz(1.235289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8265726) q[0];
sx q[0];
rz(-0.44354225) q[0];
sx q[0];
rz(-1.0916239) q[0];
rz(2.0768968) q[1];
sx q[1];
rz(-0.5263435) q[1];
sx q[1];
rz(1.975504) q[1];
rz(-2.3774316) q[2];
sx q[2];
rz(-0.66351009) q[2];
sx q[2];
rz(-1.4828975) q[2];
rz(2.06729) q[3];
sx q[3];
rz(-2.5357694) q[3];
sx q[3];
rz(2.3800935) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

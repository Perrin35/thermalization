OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4322296) q[0];
sx q[0];
rz(-0.95786434) q[0];
sx q[0];
rz(0.14444484) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(-0.52539879) q[1];
sx q[1];
rz(0.97775835) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35740556) q[0];
sx q[0];
rz(-1.6561243) q[0];
sx q[0];
rz(0.25934319) q[0];
rz(-pi) q[1];
rz(1.0563645) q[2];
sx q[2];
rz(-1.8407514) q[2];
sx q[2];
rz(1.1411238) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2797151) q[1];
sx q[1];
rz(-1.4289083) q[1];
sx q[1];
rz(0.14076294) q[1];
x q[2];
rz(1.3917189) q[3];
sx q[3];
rz(-0.79154166) q[3];
sx q[3];
rz(-2.4812428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4686761) q[2];
sx q[2];
rz(-1.9925995) q[2];
sx q[2];
rz(-2.2093175) q[2];
rz(2.9428234) q[3];
sx q[3];
rz(-1.1105744) q[3];
sx q[3];
rz(-2.1762302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4345877) q[0];
sx q[0];
rz(-2.2362237) q[0];
sx q[0];
rz(2.7804651) q[0];
rz(-1.4350285) q[1];
sx q[1];
rz(-1.7838493) q[1];
sx q[1];
rz(-0.8180058) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0571787) q[0];
sx q[0];
rz(-1.1857496) q[0];
sx q[0];
rz(-1.8684698) q[0];
x q[1];
rz(1.1992707) q[2];
sx q[2];
rz(-1.4563592) q[2];
sx q[2];
rz(-2.3105846) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4600196) q[1];
sx q[1];
rz(-0.96833723) q[1];
sx q[1];
rz(0.54089344) q[1];
x q[2];
rz(-3.0844968) q[3];
sx q[3];
rz(-1.9340056) q[3];
sx q[3];
rz(2.925194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.25257418) q[2];
sx q[2];
rz(-0.44712862) q[2];
sx q[2];
rz(1.7209631) q[2];
rz(1.3160926) q[3];
sx q[3];
rz(-2.3836453) q[3];
sx q[3];
rz(0.38823286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4017568) q[0];
sx q[0];
rz(-2.6494884) q[0];
sx q[0];
rz(2.2706568) q[0];
rz(0.31618205) q[1];
sx q[1];
rz(-2.8600287) q[1];
sx q[1];
rz(-0.2972163) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1145828) q[0];
sx q[0];
rz(-1.5073338) q[0];
sx q[0];
rz(-0.40513904) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4386114) q[2];
sx q[2];
rz(-2.2361122) q[2];
sx q[2];
rz(-0.92232982) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3618468) q[1];
sx q[1];
rz(-1.6184813) q[1];
sx q[1];
rz(-1.3839128) q[1];
x q[2];
rz(-0.86621427) q[3];
sx q[3];
rz(-0.88170393) q[3];
sx q[3];
rz(0.44486526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3729942) q[2];
sx q[2];
rz(-0.82295376) q[2];
sx q[2];
rz(-2.7139943) q[2];
rz(-1.9528495) q[3];
sx q[3];
rz(-2.5116428) q[3];
sx q[3];
rz(-2.6141613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7052085) q[0];
sx q[0];
rz(-1.5642865) q[0];
sx q[0];
rz(-2.3676681) q[0];
rz(-2.4286843) q[1];
sx q[1];
rz(-1.0293101) q[1];
sx q[1];
rz(-2.4598222) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1958774) q[0];
sx q[0];
rz(-1.8316852) q[0];
sx q[0];
rz(2.9257141) q[0];
rz(0.96659987) q[2];
sx q[2];
rz(-1.7503498) q[2];
sx q[2];
rz(2.0008848) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75979739) q[1];
sx q[1];
rz(-1.1361546) q[1];
sx q[1];
rz(-0.97310193) q[1];
x q[2];
rz(2.79822) q[3];
sx q[3];
rz(-1.2265172) q[3];
sx q[3];
rz(-1.4624649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.008808) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(0.95820367) q[2];
rz(1.0664252) q[3];
sx q[3];
rz(-1.8582148) q[3];
sx q[3];
rz(2.7619894) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.585007) q[0];
sx q[0];
rz(-1.3469232) q[0];
sx q[0];
rz(2.5812896) q[0];
rz(-0.99984461) q[1];
sx q[1];
rz(-2.9381349) q[1];
sx q[1];
rz(1.6220185) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.863127) q[0];
sx q[0];
rz(-2.0533877) q[0];
sx q[0];
rz(2.0762073) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1278473) q[2];
sx q[2];
rz(-2.1079014) q[2];
sx q[2];
rz(-2.277166) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49415627) q[1];
sx q[1];
rz(-1.7687706) q[1];
sx q[1];
rz(0.76311771) q[1];
x q[2];
rz(0.95789692) q[3];
sx q[3];
rz(-1.4688204) q[3];
sx q[3];
rz(2.0436055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6891629) q[2];
sx q[2];
rz(-2.5897557) q[2];
sx q[2];
rz(2.9186644) q[2];
rz(-3.1068504) q[3];
sx q[3];
rz(-1.7545173) q[3];
sx q[3];
rz(3.070014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8054304) q[0];
sx q[0];
rz(-0.31328377) q[0];
sx q[0];
rz(-1.0700595) q[0];
rz(-1.7806212) q[1];
sx q[1];
rz(-0.37934163) q[1];
sx q[1];
rz(1.8575352) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4053455) q[0];
sx q[0];
rz(-0.9760455) q[0];
sx q[0];
rz(2.2527762) q[0];
rz(-pi) q[1];
rz(-0.27988866) q[2];
sx q[2];
rz(-1.1942689) q[2];
sx q[2];
rz(2.1044452) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2594086) q[1];
sx q[1];
rz(-0.55667294) q[1];
sx q[1];
rz(2.5747719) q[1];
rz(1.0809903) q[3];
sx q[3];
rz(-0.59270699) q[3];
sx q[3];
rz(2.5268775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6665035) q[2];
sx q[2];
rz(-1.44375) q[2];
sx q[2];
rz(0.63759032) q[2];
rz(-2.3049138) q[3];
sx q[3];
rz(-2.0357318) q[3];
sx q[3];
rz(-1.3440514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-2.5752983) q[0];
sx q[0];
rz(-0.42995444) q[0];
sx q[0];
rz(-0.56754011) q[0];
rz(-2.7138846) q[1];
sx q[1];
rz(-1.5274915) q[1];
sx q[1];
rz(2.2033851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3283549) q[0];
sx q[0];
rz(-2.6211446) q[0];
sx q[0];
rz(-0.92207272) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6193879) q[2];
sx q[2];
rz(-2.3084547) q[2];
sx q[2];
rz(-0.39820652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.92749121) q[1];
sx q[1];
rz(-0.92988211) q[1];
sx q[1];
rz(0.44968857) q[1];
x q[2];
rz(1.2660965) q[3];
sx q[3];
rz(-2.2548171) q[3];
sx q[3];
rz(-1.5414433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2618835) q[2];
sx q[2];
rz(-1.1738913) q[2];
sx q[2];
rz(-1.7555457) q[2];
rz(1.322768) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(-3.1125606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8742074) q[0];
sx q[0];
rz(-0.33339849) q[0];
sx q[0];
rz(-1.4338795) q[0];
rz(1.8677615) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(-0.83126718) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2682876) q[0];
sx q[0];
rz(-1.9718861) q[0];
sx q[0];
rz(-1.4108765) q[0];
rz(1.7449042) q[2];
sx q[2];
rz(-2.499352) q[2];
sx q[2];
rz(1.0345392) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.23849328) q[1];
sx q[1];
rz(-2.0598754) q[1];
sx q[1];
rz(0.78519435) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23037489) q[3];
sx q[3];
rz(-1.7883693) q[3];
sx q[3];
rz(0.93254596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5780118) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(-2.1772299) q[2];
rz(-1.9780805) q[3];
sx q[3];
rz(-0.61146277) q[3];
sx q[3];
rz(-2.4826629) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37373856) q[0];
sx q[0];
rz(-0.73497325) q[0];
sx q[0];
rz(-0.97737616) q[0];
rz(1.3865698) q[1];
sx q[1];
rz(-1.3061378) q[1];
sx q[1];
rz(1.1057373) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2290105) q[0];
sx q[0];
rz(-1.6874896) q[0];
sx q[0];
rz(2.6936561) q[0];
x q[1];
rz(-2.2374723) q[2];
sx q[2];
rz(-1.6534001) q[2];
sx q[2];
rz(0.77043515) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9165009) q[1];
sx q[1];
rz(-1.3396001) q[1];
sx q[1];
rz(2.9514312) q[1];
rz(2.4615199) q[3];
sx q[3];
rz(-0.44049965) q[3];
sx q[3];
rz(0.77914731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.608312) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(2.9917955) q[2];
rz(1.3730565) q[3];
sx q[3];
rz(-1.7356197) q[3];
sx q[3];
rz(-1.8201374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49366429) q[0];
sx q[0];
rz(-0.51769185) q[0];
sx q[0];
rz(-3.0143484) q[0];
rz(-1.6607025) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(-0.1677992) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084899336) q[0];
sx q[0];
rz(-3.1057682) q[0];
sx q[0];
rz(-0.56532677) q[0];
x q[1];
rz(0.60815732) q[2];
sx q[2];
rz(-1.4434575) q[2];
sx q[2];
rz(-1.5978158) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92229453) q[1];
sx q[1];
rz(-2.1592327) q[1];
sx q[1];
rz(1.5790308) q[1];
rz(0.12331788) q[3];
sx q[3];
rz(-1.4535558) q[3];
sx q[3];
rz(0.58837147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.026713513) q[2];
sx q[2];
rz(-2.202704) q[2];
sx q[2];
rz(-2.3804469) q[2];
rz(0.090027697) q[3];
sx q[3];
rz(-1.0031676) q[3];
sx q[3];
rz(2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5933843) q[0];
sx q[0];
rz(-1.1593288) q[0];
sx q[0];
rz(2.819084) q[0];
rz(2.7535915) q[1];
sx q[1];
rz(-1.3996268) q[1];
sx q[1];
rz(-0.7849801) q[1];
rz(-0.81007304) q[2];
sx q[2];
rz(-2.0740866) q[2];
sx q[2];
rz(-1.496051) q[2];
rz(1.2896982) q[3];
sx q[3];
rz(-0.55681183) q[3];
sx q[3];
rz(-1.1435215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

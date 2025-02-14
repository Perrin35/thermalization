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
rz(2.9003484) q[0];
sx q[0];
rz(-2.9304507) q[0];
sx q[0];
rz(-0.40786064) q[0];
rz(1.3762228) q[1];
sx q[1];
rz(-1.9150182) q[1];
sx q[1];
rz(3.0768375) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5714343) q[0];
sx q[0];
rz(-1.6862094) q[0];
sx q[0];
rz(3.136271) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32081713) q[2];
sx q[2];
rz(-1.0091558) q[2];
sx q[2];
rz(-0.80036847) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9038856) q[1];
sx q[1];
rz(-0.14169417) q[1];
sx q[1];
rz(-3.0631439) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14459855) q[3];
sx q[3];
rz(-0.62718348) q[3];
sx q[3];
rz(-2.7759564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18225081) q[2];
sx q[2];
rz(-2.0472517) q[2];
sx q[2];
rz(-2.0056637) q[2];
rz(2.411339) q[3];
sx q[3];
rz(-1.7466931) q[3];
sx q[3];
rz(1.5203389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0486883) q[0];
sx q[0];
rz(-2.1850259) q[0];
sx q[0];
rz(0.055305716) q[0];
rz(-2.7703908) q[1];
sx q[1];
rz(-0.35938811) q[1];
sx q[1];
rz(-0.41195437) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5980976) q[0];
sx q[0];
rz(-2.6109222) q[0];
sx q[0];
rz(-0.12904386) q[0];
rz(-1.2367593) q[2];
sx q[2];
rz(-1.8530117) q[2];
sx q[2];
rz(-1.7441234) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4186379) q[1];
sx q[1];
rz(-1.4297863) q[1];
sx q[1];
rz(-0.73608558) q[1];
rz(-0.37825122) q[3];
sx q[3];
rz(-0.95767271) q[3];
sx q[3];
rz(2.48191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5281795) q[2];
sx q[2];
rz(-1.4643022) q[2];
sx q[2];
rz(-2.2589653) q[2];
rz(1.4334076) q[3];
sx q[3];
rz(-1.2127168) q[3];
sx q[3];
rz(-0.17361704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16647896) q[0];
sx q[0];
rz(-0.011703165) q[0];
sx q[0];
rz(2.574918) q[0];
rz(0.4176248) q[1];
sx q[1];
rz(-1.6028701) q[1];
sx q[1];
rz(1.8142987) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2921995) q[0];
sx q[0];
rz(-2.4367094) q[0];
sx q[0];
rz(-1.7017176) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94608091) q[2];
sx q[2];
rz(-1.6951121) q[2];
sx q[2];
rz(0.92783463) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23073634) q[1];
sx q[1];
rz(-1.5720782) q[1];
sx q[1];
rz(1.2627453) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2977975) q[3];
sx q[3];
rz(-1.2183172) q[3];
sx q[3];
rz(-1.6683287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7495482) q[2];
sx q[2];
rz(-1.0907402) q[2];
sx q[2];
rz(-2.2331179) q[2];
rz(-2.3467017) q[3];
sx q[3];
rz(-2.0386233) q[3];
sx q[3];
rz(3.095043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6598776) q[0];
sx q[0];
rz(-2.2292723) q[0];
sx q[0];
rz(-0.6657486) q[0];
rz(0.84716973) q[1];
sx q[1];
rz(-2.8916292) q[1];
sx q[1];
rz(-0.4096823) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85532641) q[0];
sx q[0];
rz(-0.5200035) q[0];
sx q[0];
rz(2.8306243) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80862019) q[2];
sx q[2];
rz(-1.3043772) q[2];
sx q[2];
rz(-3.1297822) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.078881513) q[1];
sx q[1];
rz(-2.1797175) q[1];
sx q[1];
rz(2.8669861) q[1];
rz(-2.014854) q[3];
sx q[3];
rz(-1.5448625) q[3];
sx q[3];
rz(-1.8166722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72982558) q[2];
sx q[2];
rz(-2.4335786) q[2];
sx q[2];
rz(-1.3789619) q[2];
rz(-2.0046558) q[3];
sx q[3];
rz(-1.7863019) q[3];
sx q[3];
rz(-0.8374477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1788504) q[0];
sx q[0];
rz(-2.0806291) q[0];
sx q[0];
rz(-2.7233997) q[0];
rz(1.7182982) q[1];
sx q[1];
rz(-1.1885252) q[1];
sx q[1];
rz(1.4994015) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20883501) q[0];
sx q[0];
rz(-1.3262188) q[0];
sx q[0];
rz(0.45758684) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33376212) q[2];
sx q[2];
rz(-1.859772) q[2];
sx q[2];
rz(-3.1372276) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5562999) q[1];
sx q[1];
rz(-1.628581) q[1];
sx q[1];
rz(-1.789566) q[1];
rz(-pi) q[2];
rz(-3.0762227) q[3];
sx q[3];
rz(-1.7002859) q[3];
sx q[3];
rz(1.1060028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2749918) q[2];
sx q[2];
rz(-1.9860705) q[2];
sx q[2];
rz(0.17821136) q[2];
rz(2.7675659) q[3];
sx q[3];
rz(-2.1686797) q[3];
sx q[3];
rz(1.1427574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48443925) q[0];
sx q[0];
rz(-1.8171808) q[0];
sx q[0];
rz(-1.665218) q[0];
rz(-0.58552512) q[1];
sx q[1];
rz(-0.89591566) q[1];
sx q[1];
rz(2.4077328) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39291362) q[0];
sx q[0];
rz(-2.3322912) q[0];
sx q[0];
rz(-2.5038263) q[0];
x q[1];
rz(-2.8346905) q[2];
sx q[2];
rz(-1.2595121) q[2];
sx q[2];
rz(-2.5800623) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1230525) q[1];
sx q[1];
rz(-2.5830659) q[1];
sx q[1];
rz(1.9473507) q[1];
x q[2];
rz(-1.9820342) q[3];
sx q[3];
rz(-0.48732584) q[3];
sx q[3];
rz(-0.05035487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0840941) q[2];
sx q[2];
rz(-1.7666662) q[2];
sx q[2];
rz(-1.6424087) q[2];
rz(-1.5205787) q[3];
sx q[3];
rz(-2.4937544) q[3];
sx q[3];
rz(2.9852941) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75477377) q[0];
sx q[0];
rz(-2.2235625) q[0];
sx q[0];
rz(1.3936438) q[0];
rz(-2.5539894) q[1];
sx q[1];
rz(-1.5005451) q[1];
sx q[1];
rz(0.29744068) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7845532) q[0];
sx q[0];
rz(-1.926657) q[0];
sx q[0];
rz(-2.3054992) q[0];
x q[1];
rz(-1.3002505) q[2];
sx q[2];
rz(-1.5890117) q[2];
sx q[2];
rz(1.8730989) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1042034) q[1];
sx q[1];
rz(-2.2504077) q[1];
sx q[1];
rz(0.012455539) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2138811) q[3];
sx q[3];
rz(-1.3731908) q[3];
sx q[3];
rz(-1.9324945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8224767) q[2];
sx q[2];
rz(-1.3603223) q[2];
sx q[2];
rz(0.3817257) q[2];
rz(-0.60025674) q[3];
sx q[3];
rz(-2.468942) q[3];
sx q[3];
rz(2.8762347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6622019) q[0];
sx q[0];
rz(-1.5150161) q[0];
sx q[0];
rz(1.8040682) q[0];
rz(-0.0094825347) q[1];
sx q[1];
rz(-0.9597221) q[1];
sx q[1];
rz(1.7705852) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0616546) q[0];
sx q[0];
rz(-1.9091055) q[0];
sx q[0];
rz(1.7163971) q[0];
rz(-0.90223899) q[2];
sx q[2];
rz(-1.1663365) q[2];
sx q[2];
rz(2.8079833) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4552729) q[1];
sx q[1];
rz(-2.3523221) q[1];
sx q[1];
rz(-1.7474921) q[1];
x q[2];
rz(-1.6608606) q[3];
sx q[3];
rz(-1.1582631) q[3];
sx q[3];
rz(1.5145258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3235772) q[2];
sx q[2];
rz(-2.24021) q[2];
sx q[2];
rz(-0.14249194) q[2];
rz(0.88322181) q[3];
sx q[3];
rz(-1.764069) q[3];
sx q[3];
rz(-1.4628598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0505117) q[0];
sx q[0];
rz(-0.091878042) q[0];
sx q[0];
rz(-1.8930513) q[0];
rz(-0.40731353) q[1];
sx q[1];
rz(-1.7946323) q[1];
sx q[1];
rz(1.9479082) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4536205) q[0];
sx q[0];
rz(-1.9000475) q[0];
sx q[0];
rz(-0.27197522) q[0];
x q[1];
rz(-2.3825109) q[2];
sx q[2];
rz(-2.1164701) q[2];
sx q[2];
rz(-2.1496043) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5919009) q[1];
sx q[1];
rz(-0.61335269) q[1];
sx q[1];
rz(-2.5332622) q[1];
rz(-pi) q[2];
x q[2];
rz(3.3832559e-05) q[3];
sx q[3];
rz(-1.8305998) q[3];
sx q[3];
rz(-2.3855748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0672037) q[2];
sx q[2];
rz(-0.95778242) q[2];
sx q[2];
rz(-1.1129145) q[2];
rz(-0.042996081) q[3];
sx q[3];
rz(-1.4837416) q[3];
sx q[3];
rz(-2.903741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0766895) q[0];
sx q[0];
rz(-1.7273128) q[0];
sx q[0];
rz(2.9458556) q[0];
rz(-1.9211357) q[1];
sx q[1];
rz(-0.38206044) q[1];
sx q[1];
rz(-2.1641796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61924926) q[0];
sx q[0];
rz(-0.27740208) q[0];
sx q[0];
rz(-2.7707556) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8711817) q[2];
sx q[2];
rz(-1.7160048) q[2];
sx q[2];
rz(-1.5086439) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4949242) q[1];
sx q[1];
rz(-2.1416602) q[1];
sx q[1];
rz(2.5350556) q[1];
rz(-0.39267003) q[3];
sx q[3];
rz(-1.304198) q[3];
sx q[3];
rz(2.6368898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0105373) q[2];
sx q[2];
rz(-0.26630339) q[2];
sx q[2];
rz(-1.6923426) q[2];
rz(0.87861711) q[3];
sx q[3];
rz(-2.0802458) q[3];
sx q[3];
rz(-1.7884458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5129678) q[0];
sx q[0];
rz(-1.1911363) q[0];
sx q[0];
rz(-1.3918899) q[0];
rz(-1.3799008) q[1];
sx q[1];
rz(-1.5329783) q[1];
sx q[1];
rz(3.0402532) q[1];
rz(-2.352667) q[2];
sx q[2];
rz(-2.8876705) q[2];
sx q[2];
rz(-1.3063274) q[2];
rz(2.8158549) q[3];
sx q[3];
rz(-0.67254638) q[3];
sx q[3];
rz(-1.2253958) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

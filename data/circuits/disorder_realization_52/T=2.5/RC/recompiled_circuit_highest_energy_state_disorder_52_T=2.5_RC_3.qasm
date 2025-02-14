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
rz(-1.988451) q[0];
sx q[0];
rz(-2.3260131) q[0];
sx q[0];
rz(-2.3834035) q[0];
rz(-0.012501333) q[1];
sx q[1];
rz(-1.3036417) q[1];
sx q[1];
rz(1.5703896) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0051828) q[0];
sx q[0];
rz(-1.9290975) q[0];
sx q[0];
rz(-0.65922064) q[0];
rz(-pi) q[1];
rz(1.9704029) q[2];
sx q[2];
rz(-1.6328703) q[2];
sx q[2];
rz(0.8555748) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.058152288) q[1];
sx q[1];
rz(-1.6279164) q[1];
sx q[1];
rz(-1.9316462) q[1];
rz(1.38405) q[3];
sx q[3];
rz(-2.8272326) q[3];
sx q[3];
rz(-1.7160742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5300753) q[2];
sx q[2];
rz(-3.1334183) q[2];
sx q[2];
rz(2.7008936) q[2];
rz(3.0618482) q[3];
sx q[3];
rz(-0.00010448797) q[3];
sx q[3];
rz(-2.0174446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1752862) q[0];
sx q[0];
rz(-3.0847302) q[0];
sx q[0];
rz(0.16800198) q[0];
rz(0.02027823) q[1];
sx q[1];
rz(-2.8333277) q[1];
sx q[1];
rz(1.5365938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.080606) q[0];
sx q[0];
rz(-1.3362552) q[0];
sx q[0];
rz(-1.3479665) q[0];
rz(-pi) q[1];
rz(-3.1390983) q[2];
sx q[2];
rz(-1.5809142) q[2];
sx q[2];
rz(1.5455139) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.78129301) q[1];
sx q[1];
rz(-1.5661245) q[1];
sx q[1];
rz(-3.1404936) q[1];
rz(-1.5797938) q[3];
sx q[3];
rz(-1.5285042) q[3];
sx q[3];
rz(-2.472714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7961879) q[2];
sx q[2];
rz(-2.2296843) q[2];
sx q[2];
rz(-1.7339285) q[2];
rz(-2.0965072) q[3];
sx q[3];
rz(-3.0920691) q[3];
sx q[3];
rz(2.869587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3097836) q[0];
sx q[0];
rz(-0.97667664) q[0];
sx q[0];
rz(-2.580544) q[0];
rz(-0.27684119) q[1];
sx q[1];
rz(-3.1287153) q[1];
sx q[1];
rz(-1.3078088) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9914068) q[0];
sx q[0];
rz(-1.5298109) q[0];
sx q[0];
rz(-1.8677554) q[0];
rz(-0.00066999992) q[2];
sx q[2];
rz(-3.1342034) q[2];
sx q[2];
rz(2.0787897) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6456994) q[1];
sx q[1];
rz(-0.99826854) q[1];
sx q[1];
rz(3.0719724) q[1];
rz(-pi) q[2];
rz(-0.8933634) q[3];
sx q[3];
rz(-1.9382538) q[3];
sx q[3];
rz(-1.8830255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7967367) q[2];
sx q[2];
rz(-0.00011809706) q[2];
sx q[2];
rz(2.5522088) q[2];
rz(1.2423337) q[3];
sx q[3];
rz(-0.012367736) q[3];
sx q[3];
rz(1.3602268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0068483343) q[0];
sx q[0];
rz(-2.6297748) q[0];
sx q[0];
rz(1.3486598) q[0];
rz(3.1349365) q[1];
sx q[1];
rz(-1.8204047) q[1];
sx q[1];
rz(-3.1080918) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.63147) q[0];
sx q[0];
rz(-2.2167248) q[0];
sx q[0];
rz(0.21696802) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5662996) q[2];
sx q[2];
rz(-1.6904313) q[2];
sx q[2];
rz(0.41882354) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9814947) q[1];
sx q[1];
rz(-1.3019053) q[1];
sx q[1];
rz(0.010543028) q[1];
x q[2];
rz(1.7351522) q[3];
sx q[3];
rz(-1.6911611) q[3];
sx q[3];
rz(0.09935483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.032430705) q[2];
sx q[2];
rz(-0.0065294821) q[2];
sx q[2];
rz(0.29864857) q[2];
rz(-1.8715035) q[3];
sx q[3];
rz(-0.01615571) q[3];
sx q[3];
rz(0.0531918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3700767) q[0];
sx q[0];
rz(-1.6271485) q[0];
sx q[0];
rz(2.694743) q[0];
rz(2.9554548) q[1];
sx q[1];
rz(-0.061807241) q[1];
sx q[1];
rz(-1.7284547) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9139149) q[0];
sx q[0];
rz(-2.9083038) q[0];
sx q[0];
rz(1.4248217) q[0];
x q[1];
rz(2.9252626) q[2];
sx q[2];
rz(-1.1573175) q[2];
sx q[2];
rz(-1.0461996) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4889989) q[1];
sx q[1];
rz(-1.6203383) q[1];
sx q[1];
rz(-0.045157305) q[1];
rz(-pi) q[2];
rz(-0.90425332) q[3];
sx q[3];
rz(-2.7891141) q[3];
sx q[3];
rz(0.74732399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3073005) q[2];
sx q[2];
rz(-1.5520381) q[2];
sx q[2];
rz(-2.6429122) q[2];
rz(0.57791609) q[3];
sx q[3];
rz(-0.48360616) q[3];
sx q[3];
rz(-0.59670603) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96108288) q[0];
sx q[0];
rz(-2.0208277) q[0];
sx q[0];
rz(2.7951796) q[0];
rz(2.5397884) q[1];
sx q[1];
rz(-1.5608984) q[1];
sx q[1];
rz(-0.7535038) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2997871) q[0];
sx q[0];
rz(-2.8827169) q[0];
sx q[0];
rz(-0.73200925) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97917231) q[2];
sx q[2];
rz(-0.1340999) q[2];
sx q[2];
rz(0.025452415) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0776163) q[1];
sx q[1];
rz(-1.9850296) q[1];
sx q[1];
rz(2.2760681) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4077134) q[3];
sx q[3];
rz(-1.4797416) q[3];
sx q[3];
rz(1.0331105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.57009131) q[2];
sx q[2];
rz(-0.0034905958) q[2];
sx q[2];
rz(-1.6174779) q[2];
rz(-0.12024719) q[3];
sx q[3];
rz(-3.1382939) q[3];
sx q[3];
rz(-0.53774589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26871249) q[0];
sx q[0];
rz(-2.217642) q[0];
sx q[0];
rz(0.22126108) q[0];
rz(1.4606754) q[1];
sx q[1];
rz(-0.9333846) q[1];
sx q[1];
rz(3.0642919) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0780405) q[0];
sx q[0];
rz(-1.5722599) q[0];
sx q[0];
rz(1.5290676) q[0];
x q[1];
rz(0.004782025) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(0.21172571) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.12880023) q[1];
sx q[1];
rz(-0.18928738) q[1];
sx q[1];
rz(0.38893338) q[1];
x q[2];
rz(0.12760136) q[3];
sx q[3];
rz(-1.3584879) q[3];
sx q[3];
rz(1.177976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7920502) q[2];
sx q[2];
rz(-3.1304066) q[2];
sx q[2];
rz(-0.95996094) q[2];
rz(-2.8121484) q[3];
sx q[3];
rz(-3.1335148) q[3];
sx q[3];
rz(0.85028696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.195381) q[0];
sx q[0];
rz(-2.52849) q[0];
sx q[0];
rz(-0.10928133) q[0];
rz(-0.37336135) q[1];
sx q[1];
rz(-2.3318718) q[1];
sx q[1];
rz(-1.2304617) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62723535) q[0];
sx q[0];
rz(-0.97313303) q[0];
sx q[0];
rz(2.3742832) q[0];
rz(2.3745499) q[2];
sx q[2];
rz(-0.27077507) q[2];
sx q[2];
rz(0.8052288) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.44430393) q[1];
sx q[1];
rz(-3.0737801) q[1];
sx q[1];
rz(-1.4013002) q[1];
rz(-pi) q[2];
rz(-0.91992141) q[3];
sx q[3];
rz(-1.3215142) q[3];
sx q[3];
rz(-1.9712308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5665148) q[2];
sx q[2];
rz(-1.9058303) q[2];
sx q[2];
rz(-1.3285948) q[2];
rz(-1.7447507) q[3];
sx q[3];
rz(-3.1378523) q[3];
sx q[3];
rz(-2.1197135) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0777271) q[0];
sx q[0];
rz(-1.6985748) q[0];
sx q[0];
rz(2.5685837) q[0];
rz(-0.30814463) q[1];
sx q[1];
rz(-2.7317218) q[1];
sx q[1];
rz(2.134197) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6992189) q[0];
sx q[0];
rz(-1.464252) q[0];
sx q[0];
rz(-0.005597896) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.043906004) q[2];
sx q[2];
rz(-1.708235) q[2];
sx q[2];
rz(3.1071747) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9646109) q[1];
sx q[1];
rz(-1.4876502) q[1];
sx q[1];
rz(0.040772922) q[1];
rz(1.6009523) q[3];
sx q[3];
rz(-1.5573676) q[3];
sx q[3];
rz(-1.7755938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8241626) q[2];
sx q[2];
rz(-2.514826) q[2];
sx q[2];
rz(-2.7516348) q[2];
rz(-0.069084875) q[3];
sx q[3];
rz(-0.0091113541) q[3];
sx q[3];
rz(0.33153427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-3.1214509) q[0];
sx q[0];
rz(-2.3905601) q[0];
sx q[0];
rz(-2.6556515) q[0];
rz(0.87156975) q[1];
sx q[1];
rz(-1.3078682) q[1];
sx q[1];
rz(-1.4942687) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93372351) q[0];
sx q[0];
rz(-1.2656801) q[0];
sx q[0];
rz(-0.56201571) q[0];
rz(-pi) q[1];
rz(1.5266085) q[2];
sx q[2];
rz(-0.95643294) q[2];
sx q[2];
rz(3.0725266) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.20032665) q[1];
sx q[1];
rz(-1.8706053) q[1];
sx q[1];
rz(-0.34383641) q[1];
x q[2];
rz(-1.5288197) q[3];
sx q[3];
rz(-1.4577565) q[3];
sx q[3];
rz(0.50769317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5668874) q[2];
sx q[2];
rz(-0.042782728) q[2];
sx q[2];
rz(-3.1097143) q[2];
rz(0.77830642) q[3];
sx q[3];
rz(-0.0068155546) q[3];
sx q[3];
rz(-0.2955029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7185709) q[0];
sx q[0];
rz(-1.6091249) q[0];
sx q[0];
rz(-1.3269497) q[0];
rz(-0.12693916) q[1];
sx q[1];
rz(-2.9025684) q[1];
sx q[1];
rz(-2.92166) q[1];
rz(1.4311287) q[2];
sx q[2];
rz(-1.5652547) q[2];
sx q[2];
rz(1.8143285) q[2];
rz(-0.89613468) q[3];
sx q[3];
rz(-1.5265161) q[3];
sx q[3];
rz(2.0731887) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

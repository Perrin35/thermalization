OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(2.2979484) q[0];
sx q[0];
rz(9.2568682) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(-2.8462703) q[1];
sx q[1];
rz(0.056161031) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6457155) q[0];
sx q[0];
rz(-1.0951395) q[0];
sx q[0];
rz(0.23947421) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21284717) q[2];
sx q[2];
rz(-2.2058862) q[2];
sx q[2];
rz(1.1320621) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4178884) q[1];
sx q[1];
rz(-1.0350409) q[1];
sx q[1];
rz(0.31236155) q[1];
x q[2];
rz(2.8817301) q[3];
sx q[3];
rz(-1.6136323) q[3];
sx q[3];
rz(0.44997893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.37796676) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(-0.4326694) q[2];
rz(1.9487322) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(-0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52779657) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(-1.8288076) q[0];
rz(0.20547543) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(1.9899433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5713455) q[0];
sx q[0];
rz(-2.3803108) q[0];
sx q[0];
rz(-2.0061357) q[0];
x q[1];
rz(-2.0552911) q[2];
sx q[2];
rz(-1.0420348) q[2];
sx q[2];
rz(-0.27258401) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0484867) q[1];
sx q[1];
rz(-1.152532) q[1];
sx q[1];
rz(2.6622245) q[1];
x q[2];
rz(-1.1851951) q[3];
sx q[3];
rz(-1.10023) q[3];
sx q[3];
rz(3.1009931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0097222086) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(0.37718537) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.8310228) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(-3.1047399) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(3.085014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7650334) q[0];
sx q[0];
rz(-0.74900904) q[0];
sx q[0];
rz(0.88699938) q[0];
x q[1];
rz(-2.8851606) q[2];
sx q[2];
rz(-1.4415381) q[2];
sx q[2];
rz(-1.3916707) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6371582) q[1];
sx q[1];
rz(-1.27099) q[1];
sx q[1];
rz(1.7791041) q[1];
x q[2];
rz(-2.4967381) q[3];
sx q[3];
rz(-1.8861024) q[3];
sx q[3];
rz(0.23526084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8686707) q[2];
sx q[2];
rz(-1.4977095) q[2];
sx q[2];
rz(-0.92612129) q[2];
rz(-2.5849294) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91519231) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(0.74209374) q[0];
rz(2.0023951) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(0.46359584) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84309834) q[0];
sx q[0];
rz(-1.2756186) q[0];
sx q[0];
rz(0.85423268) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5092588) q[2];
sx q[2];
rz(-1.2386285) q[2];
sx q[2];
rz(0.57519826) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3568748) q[1];
sx q[1];
rz(-2.3349635) q[1];
sx q[1];
rz(-2.9033317) q[1];
x q[2];
rz(-0.47834088) q[3];
sx q[3];
rz(-2.0739177) q[3];
sx q[3];
rz(0.92418811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7745557) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(0.049499361) q[2];
rz(0.1285304) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(-0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13609919) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(-0.29770011) q[0];
rz(-2.659335) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(-0.94435) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56086841) q[0];
sx q[0];
rz(-1.3618016) q[0];
sx q[0];
rz(3.0955549) q[0];
rz(0.77563939) q[2];
sx q[2];
rz(-1.2795942) q[2];
sx q[2];
rz(1.8679384) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.797232) q[1];
sx q[1];
rz(-1.5096944) q[1];
sx q[1];
rz(-3.1394715) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3145507) q[3];
sx q[3];
rz(-1.1503997) q[3];
sx q[3];
rz(-0.30403578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8828316) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(-3.0333701) q[2];
rz(3.1392858) q[3];
sx q[3];
rz(-1.6121515) q[3];
sx q[3];
rz(-0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7047983) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(2.5571402) q[0];
rz(-0.8862409) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(-3.086673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7391101) q[0];
sx q[0];
rz(-1.4882898) q[0];
sx q[0];
rz(1.301618) q[0];
rz(3.1228742) q[2];
sx q[2];
rz(-1.2476377) q[2];
sx q[2];
rz(-1.9402372) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7683148) q[1];
sx q[1];
rz(-1.0265961) q[1];
sx q[1];
rz(1.1276223) q[1];
x q[2];
rz(2.1543598) q[3];
sx q[3];
rz(-1.1757441) q[3];
sx q[3];
rz(2.6810255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75446689) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(1.0167271) q[2];
rz(2.5975442) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6761557) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(0.28453919) q[0];
rz(-0.94447213) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(-0.91032666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1032573) q[0];
sx q[0];
rz(-0.084780134) q[0];
sx q[0];
rz(0.8707365) q[0];
rz(-0.89332135) q[2];
sx q[2];
rz(-2.6258694) q[2];
sx q[2];
rz(2.233778) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0525166) q[1];
sx q[1];
rz(-2.0564582) q[1];
sx q[1];
rz(-1.3658701) q[1];
rz(-pi) q[2];
rz(0.63776871) q[3];
sx q[3];
rz(-0.83003269) q[3];
sx q[3];
rz(2.4142767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(-2.2195623) q[2];
rz(-2.5743124) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(-2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2475125) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(0.12938736) q[0];
rz(2.5091876) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(2.8410889) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26586543) q[0];
sx q[0];
rz(-0.89571307) q[0];
sx q[0];
rz(0.48203326) q[0];
x q[1];
rz(2.1883165) q[2];
sx q[2];
rz(-0.91070181) q[2];
sx q[2];
rz(2.2613139) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.955519) q[1];
sx q[1];
rz(-1.3961853) q[1];
sx q[1];
rz(0.34592918) q[1];
x q[2];
rz(-1.9037876) q[3];
sx q[3];
rz(-2.3559642) q[3];
sx q[3];
rz(1.1803407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5552716) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(0.78197455) q[2];
rz(0.55109763) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5683811) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(3.0138299) q[0];
rz(2.5993775) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(-2.382747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42620537) q[0];
sx q[0];
rz(-1.1466768) q[0];
sx q[0];
rz(-0.018652648) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27300948) q[2];
sx q[2];
rz(-2.5148279) q[2];
sx q[2];
rz(0.11944709) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2821741) q[1];
sx q[1];
rz(-2.4255883) q[1];
sx q[1];
rz(0.23846682) q[1];
rz(-pi) q[2];
rz(2.2364053) q[3];
sx q[3];
rz(-1.5822516) q[3];
sx q[3];
rz(-0.36597914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1252497) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(-2.6515567) q[2];
rz(-1.7193433) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(-1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816417) q[0];
sx q[0];
rz(-0.61976969) q[0];
sx q[0];
rz(3.066257) q[0];
rz(-2.244859) q[1];
sx q[1];
rz(-1.9672085) q[1];
sx q[1];
rz(-2.5316701) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2595554) q[0];
sx q[0];
rz(-1.8543188) q[0];
sx q[0];
rz(-2.3638704) q[0];
x q[1];
rz(2.5209849) q[2];
sx q[2];
rz(-0.45244103) q[2];
sx q[2];
rz(1.1021745) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6590609) q[1];
sx q[1];
rz(-0.7809124) q[1];
sx q[1];
rz(-2.798583) q[1];
rz(-2.6033953) q[3];
sx q[3];
rz(-0.81848577) q[3];
sx q[3];
rz(-1.8173816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.909409) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(-0.71371901) q[2];
rz(-2.7632726) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(-0.87987125) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80355766) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(-2.4907885) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(1.4177633) q[2];
sx q[2];
rz(-2.9433708) q[2];
sx q[2];
rz(2.3797258) q[2];
rz(1.0999023) q[3];
sx q[3];
rz(-1.3445911) q[3];
sx q[3];
rz(0.54683987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

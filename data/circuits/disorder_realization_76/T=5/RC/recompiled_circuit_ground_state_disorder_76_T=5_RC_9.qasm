OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4467093) q[0];
sx q[0];
rz(-2.0599685) q[0];
sx q[0];
rz(-2.4021436) q[0];
rz(2.3820355) q[1];
sx q[1];
rz(1.8269202) q[1];
sx q[1];
rz(4.8575525) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6524871) q[0];
sx q[0];
rz(-1.28363) q[0];
sx q[0];
rz(1.0190359) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6447634) q[2];
sx q[2];
rz(-1.8112) q[2];
sx q[2];
rz(-0.59076004) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.4233746) q[1];
sx q[1];
rz(-1.4932695) q[1];
sx q[1];
rz(1.422387) q[1];
x q[2];
rz(1.6207148) q[3];
sx q[3];
rz(-2.422159) q[3];
sx q[3];
rz(1.0843383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2322959) q[2];
sx q[2];
rz(-0.84315073) q[2];
sx q[2];
rz(2.9006531) q[2];
rz(-3.1124034) q[3];
sx q[3];
rz(-1.8029282) q[3];
sx q[3];
rz(1.1791112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9772684) q[0];
sx q[0];
rz(-1.203953) q[0];
sx q[0];
rz(2.6569195) q[0];
rz(2.4618705) q[1];
sx q[1];
rz(-1.8599963) q[1];
sx q[1];
rz(1.9814804) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74014464) q[0];
sx q[0];
rz(-2.8370246) q[0];
sx q[0];
rz(0.061621678) q[0];
x q[1];
rz(-0.76267879) q[2];
sx q[2];
rz(-2.2644223) q[2];
sx q[2];
rz(2.9528303) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1801123) q[1];
sx q[1];
rz(-2.4870076) q[1];
sx q[1];
rz(-2.727319) q[1];
rz(-pi) q[2];
rz(2.6187702) q[3];
sx q[3];
rz(-1.2827875) q[3];
sx q[3];
rz(0.097974591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.844187) q[2];
sx q[2];
rz(-2.83941) q[2];
sx q[2];
rz(0.45787946) q[2];
rz(1.1229905) q[3];
sx q[3];
rz(-2.1419958) q[3];
sx q[3];
rz(-0.56639731) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26894012) q[0];
sx q[0];
rz(-0.70916969) q[0];
sx q[0];
rz(0.031524468) q[0];
rz(0.28745502) q[1];
sx q[1];
rz(-2.2645686) q[1];
sx q[1];
rz(-1.9256176) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4327964) q[0];
sx q[0];
rz(-2.2026718) q[0];
sx q[0];
rz(0.89544501) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38162614) q[2];
sx q[2];
rz(-0.51593057) q[2];
sx q[2];
rz(-2.6037773) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.79920125) q[1];
sx q[1];
rz(-1.2974129) q[1];
sx q[1];
rz(-2.8552516) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76283703) q[3];
sx q[3];
rz(-0.82218542) q[3];
sx q[3];
rz(1.959182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1027801) q[2];
sx q[2];
rz(-1.5422042) q[2];
sx q[2];
rz(0.83596027) q[2];
rz(-1.7838259) q[3];
sx q[3];
rz(-2.416555) q[3];
sx q[3];
rz(2.3939705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8589856) q[0];
sx q[0];
rz(-1.7361807) q[0];
sx q[0];
rz(0.080168515) q[0];
rz(2.0591586) q[1];
sx q[1];
rz(-2.9083462) q[1];
sx q[1];
rz(-1.8128043) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6850615) q[0];
sx q[0];
rz(-2.079112) q[0];
sx q[0];
rz(1.9446745) q[0];
rz(1.7432418) q[2];
sx q[2];
rz(-1.2524464) q[2];
sx q[2];
rz(-2.6686252) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7693682) q[1];
sx q[1];
rz(-1.9980248) q[1];
sx q[1];
rz(0.15924304) q[1];
x q[2];
rz(-2.8058047) q[3];
sx q[3];
rz(-1.3779252) q[3];
sx q[3];
rz(-1.646281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3749915) q[2];
sx q[2];
rz(-2.1583755) q[2];
sx q[2];
rz(0.90744606) q[2];
rz(-2.4713016) q[3];
sx q[3];
rz(-0.82796103) q[3];
sx q[3];
rz(-1.7214187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9012673) q[0];
sx q[0];
rz(-0.87258744) q[0];
sx q[0];
rz(2.7101044) q[0];
rz(1.1314499) q[1];
sx q[1];
rz(-1.1791469) q[1];
sx q[1];
rz(-0.44050899) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90135306) q[0];
sx q[0];
rz(-0.47751891) q[0];
sx q[0];
rz(-1.8285455) q[0];
x q[1];
rz(-1.155987) q[2];
sx q[2];
rz(-1.7602663) q[2];
sx q[2];
rz(-2.1304325) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1738834) q[1];
sx q[1];
rz(-1.6102127) q[1];
sx q[1];
rz(-2.8110678) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2045361) q[3];
sx q[3];
rz(-1.0770633) q[3];
sx q[3];
rz(0.99668324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12209192) q[2];
sx q[2];
rz(-2.2152405) q[2];
sx q[2];
rz(-0.41395536) q[2];
rz(1.7533938) q[3];
sx q[3];
rz(-1.5475169) q[3];
sx q[3];
rz(-3.0164914) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51467657) q[0];
sx q[0];
rz(-1.4056982) q[0];
sx q[0];
rz(-2.2077014) q[0];
rz(2.4712708) q[1];
sx q[1];
rz(-1.8972081) q[1];
sx q[1];
rz(1.8062887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2349699) q[0];
sx q[0];
rz(-0.9714533) q[0];
sx q[0];
rz(-0.12110981) q[0];
rz(-pi) q[1];
x q[1];
rz(2.983903) q[2];
sx q[2];
rz(-1.5957498) q[2];
sx q[2];
rz(1.4098997) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0508479) q[1];
sx q[1];
rz(-0.11493348) q[1];
sx q[1];
rz(1.7687709) q[1];
rz(-pi) q[2];
x q[2];
rz(0.011739775) q[3];
sx q[3];
rz(-0.63705963) q[3];
sx q[3];
rz(-2.5979192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2495217) q[2];
sx q[2];
rz(-2.8688909) q[2];
sx q[2];
rz(-1.8708694) q[2];
rz(-1.7719841) q[3];
sx q[3];
rz(-1.2415875) q[3];
sx q[3];
rz(2.6197267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524566) q[0];
sx q[0];
rz(-0.58681762) q[0];
sx q[0];
rz(2.9879046) q[0];
rz(-2.9391089) q[1];
sx q[1];
rz(-1.2362365) q[1];
sx q[1];
rz(2.1930146) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.089854) q[0];
sx q[0];
rz(-1.3696545) q[0];
sx q[0];
rz(2.5178943) q[0];
x q[1];
rz(1.1467298) q[2];
sx q[2];
rz(-2.4575007) q[2];
sx q[2];
rz(0.37307326) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2413509) q[1];
sx q[1];
rz(-0.35120041) q[1];
sx q[1];
rz(-2.3022149) q[1];
rz(-3.0074869) q[3];
sx q[3];
rz(-0.92536345) q[3];
sx q[3];
rz(-0.2303309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.004868) q[2];
sx q[2];
rz(-2.0437045) q[2];
sx q[2];
rz(3.1401805) q[2];
rz(3.0794365) q[3];
sx q[3];
rz(-1.4096189) q[3];
sx q[3];
rz(-0.96127659) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.38178) q[0];
sx q[0];
rz(-2.3842922) q[0];
sx q[0];
rz(1.9267474) q[0];
rz(-0.75792056) q[1];
sx q[1];
rz(-2.6140723) q[1];
sx q[1];
rz(-2.5835999) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76282952) q[0];
sx q[0];
rz(-0.42233322) q[0];
sx q[0];
rz(0.93393737) q[0];
rz(0.2336913) q[2];
sx q[2];
rz(-0.56392852) q[2];
sx q[2];
rz(-0.83966161) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.26436603) q[1];
sx q[1];
rz(-2.5747882) q[1];
sx q[1];
rz(2.1071438) q[1];
rz(-pi) q[2];
rz(0.026838035) q[3];
sx q[3];
rz(-1.8478113) q[3];
sx q[3];
rz(2.7153496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5743635) q[2];
sx q[2];
rz(-1.1939253) q[2];
sx q[2];
rz(0.26926678) q[2];
rz(-1.1632129) q[3];
sx q[3];
rz(-2.5389157) q[3];
sx q[3];
rz(2.9009624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6729386) q[0];
sx q[0];
rz(-2.1042295) q[0];
sx q[0];
rz(-2.0682251) q[0];
rz(2.0138373) q[1];
sx q[1];
rz(-0.22722166) q[1];
sx q[1];
rz(-0.55666298) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4024324) q[0];
sx q[0];
rz(-2.7910829) q[0];
sx q[0];
rz(1.9950161) q[0];
x q[1];
rz(2.1416592) q[2];
sx q[2];
rz(-2.0682671) q[2];
sx q[2];
rz(-1.9299049) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5334415) q[1];
sx q[1];
rz(-1.9284298) q[1];
sx q[1];
rz(-2.3945892) q[1];
x q[2];
rz(-2.7438874) q[3];
sx q[3];
rz(-2.8880062) q[3];
sx q[3];
rz(1.2802779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.90159455) q[2];
sx q[2];
rz(-2.0549213) q[2];
sx q[2];
rz(-1.1617917) q[2];
rz(2.6084172) q[3];
sx q[3];
rz(-2.6170001) q[3];
sx q[3];
rz(-2.9840792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6530782) q[0];
sx q[0];
rz(-0.97452679) q[0];
sx q[0];
rz(2.5467806) q[0];
rz(-0.57890233) q[1];
sx q[1];
rz(-1.6659104) q[1];
sx q[1];
rz(0.65931177) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3357617) q[0];
sx q[0];
rz(-1.7044401) q[0];
sx q[0];
rz(-2.5138096) q[0];
rz(-2.012678) q[2];
sx q[2];
rz(-2.127248) q[2];
sx q[2];
rz(-0.2054727) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54592268) q[1];
sx q[1];
rz(-2.8575142) q[1];
sx q[1];
rz(1.0417263) q[1];
rz(-pi) q[2];
rz(-2.3385184) q[3];
sx q[3];
rz(-2.7702906) q[3];
sx q[3];
rz(-2.7367531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.351563) q[2];
sx q[2];
rz(-1.943925) q[2];
sx q[2];
rz(0.053038049) q[2];
rz(1.3430345) q[3];
sx q[3];
rz(-1.2767295) q[3];
sx q[3];
rz(-3.067335) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2863083) q[0];
sx q[0];
rz(-1.3636148) q[0];
sx q[0];
rz(-1.6028945) q[0];
rz(-3.042649) q[1];
sx q[1];
rz(-1.3700486) q[1];
sx q[1];
rz(0.055421967) q[1];
rz(-1.5794051) q[2];
sx q[2];
rz(-1.1653524) q[2];
sx q[2];
rz(-0.63962519) q[2];
rz(2.9704499) q[3];
sx q[3];
rz(-0.74872331) q[3];
sx q[3];
rz(-3.1414422) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

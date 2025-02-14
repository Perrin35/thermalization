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
rz(1.1531416) q[0];
sx q[0];
rz(2.3260131) q[0];
sx q[0];
rz(10.182967) q[0];
rz(3.1290913) q[1];
sx q[1];
rz(-1.8379509) q[1];
sx q[1];
rz(-1.5703896) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2822097) q[0];
sx q[0];
rz(-2.4042685) q[0];
sx q[0];
rz(0.54873772) q[0];
rz(-pi) q[1];
rz(-3.0742253) q[2];
sx q[2];
rz(-1.172003) q[2];
sx q[2];
rz(-0.68902868) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4788073) q[1];
sx q[1];
rz(-0.36514716) q[1];
sx q[1];
rz(-1.7313596) q[1];
rz(-pi) q[2];
rz(0.060293555) q[3];
sx q[3];
rz(-1.8795085) q[3];
sx q[3];
rz(1.2293881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9663064) q[0];
sx q[0];
rz(-3.0847302) q[0];
sx q[0];
rz(-2.9735907) q[0];
rz(-3.1213144) q[1];
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
rz(2.5791721) q[0];
sx q[0];
rz(-1.7874266) q[0];
sx q[0];
rz(-2.9013293) q[0];
rz(-pi) q[1];
x q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3602996) q[1];
sx q[1];
rz(-1.5754682) q[1];
sx q[1];
rz(0.0010990573) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.042293799) q[3];
sx q[3];
rz(-1.5618069) q[3];
sx q[3];
rz(-2.2400554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.34540471) q[2];
sx q[2];
rz(-0.91190839) q[2];
sx q[2];
rz(1.7339285) q[2];
rz(2.0965072) q[3];
sx q[3];
rz(-0.049523517) q[3];
sx q[3];
rz(-0.27200562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8318091) q[0];
sx q[0];
rz(-2.164916) q[0];
sx q[0];
rz(-2.580544) q[0];
rz(0.27684119) q[1];
sx q[1];
rz(-3.1287153) q[1];
sx q[1];
rz(-1.8337839) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55371743) q[0];
sx q[0];
rz(-2.8419015) q[0];
sx q[0];
rz(-1.7100348) q[0];
rz(-pi) q[1];
x q[1];
rz(0.00066999992) q[2];
sx q[2];
rz(-0.0073892842) q[2];
sx q[2];
rz(-1.0628029) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3678811) q[1];
sx q[1];
rz(-0.57627541) q[1];
sx q[1];
rz(-1.4632844) q[1];
rz(-pi) q[2];
rz(-0.4588608) q[3];
sx q[3];
rz(-0.94597497) q[3];
sx q[3];
rz(-3.1106839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3448559) q[2];
sx q[2];
rz(-0.00011809706) q[2];
sx q[2];
rz(0.5893839) q[2];
rz(-1.2423337) q[3];
sx q[3];
rz(-0.012367736) q[3];
sx q[3];
rz(-1.3602268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(0.0068483343) q[0];
sx q[0];
rz(-2.6297748) q[0];
sx q[0];
rz(-1.7929329) q[0];
rz(0.0066561247) q[1];
sx q[1];
rz(-1.8204047) q[1];
sx q[1];
rz(-0.033500813) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.63147) q[0];
sx q[0];
rz(-2.2167248) q[0];
sx q[0];
rz(2.9246246) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0219565) q[2];
sx q[2];
rz(-1.575261) q[2];
sx q[2];
rz(1.9890832) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0211612) q[1];
sx q[1];
rz(-2.8725) q[1];
sx q[1];
rz(-1.6090367) q[1];
rz(-pi) q[2];
rz(2.2073645) q[3];
sx q[3];
rz(-2.9381972) q[3];
sx q[3];
rz(1.0435085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1091619) q[2];
sx q[2];
rz(-3.1350632) q[2];
sx q[2];
rz(-2.8429441) q[2];
rz(1.2700891) q[3];
sx q[3];
rz(-3.1254369) q[3];
sx q[3];
rz(-0.0531918) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3700767) q[0];
sx q[0];
rz(-1.5144441) q[0];
sx q[0];
rz(2.694743) q[0];
rz(2.9554548) q[1];
sx q[1];
rz(-3.0797854) q[1];
sx q[1];
rz(-1.4131379) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37765682) q[0];
sx q[0];
rz(-1.8015588) q[0];
sx q[0];
rz(-0.034548977) q[0];
rz(-pi) q[1];
rz(2.0257641) q[2];
sx q[2];
rz(-0.46374292) q[2];
sx q[2];
rz(-0.5456897) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4889989) q[1];
sx q[1];
rz(-1.5212544) q[1];
sx q[1];
rz(0.045157305) q[1];
rz(-pi) q[2];
rz(-0.90425332) q[3];
sx q[3];
rz(-0.35247856) q[3];
sx q[3];
rz(-0.74732399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3073005) q[2];
sx q[2];
rz(-1.5895546) q[2];
sx q[2];
rz(0.49868047) q[2];
rz(-0.57791609) q[3];
sx q[3];
rz(-2.6579865) q[3];
sx q[3];
rz(2.5448866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805098) q[0];
sx q[0];
rz(-1.120765) q[0];
sx q[0];
rz(0.34641308) q[0];
rz(-0.60180426) q[1];
sx q[1];
rz(-1.5806942) q[1];
sx q[1];
rz(0.7535038) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44418535) q[0];
sx q[0];
rz(-1.7427398) q[0];
sx q[0];
rz(-0.19449046) q[0];
x q[1];
rz(3.0664938) q[2];
sx q[2];
rz(-1.6819998) q[2];
sx q[2];
rz(-2.5203343) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9484472) q[1];
sx q[1];
rz(-2.3420534) q[1];
sx q[1];
rz(0.97481291) q[1];
rz(0.092272225) q[3];
sx q[3];
rz(-1.408395) q[3];
sx q[3];
rz(0.55264651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.57009131) q[2];
sx q[2];
rz(-0.0034905958) q[2];
sx q[2];
rz(-1.6174779) q[2];
rz(0.12024719) q[3];
sx q[3];
rz(-3.1382939) q[3];
sx q[3];
rz(0.53774589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26871249) q[0];
sx q[0];
rz(-2.217642) q[0];
sx q[0];
rz(-0.22126108) q[0];
rz(-1.6809173) q[1];
sx q[1];
rz(-2.2082081) q[1];
sx q[1];
rz(-3.0642919) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0635522) q[0];
sx q[0];
rz(-1.5722599) q[0];
sx q[0];
rz(-1.6125251) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0650778) q[2];
sx q[2];
rz(-3.1315127) q[2];
sx q[2];
rz(-2.8590157) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.52409808) q[1];
sx q[1];
rz(-1.7457944) q[1];
sx q[1];
rz(-1.4982759) q[1];
x q[2];
rz(1.3568001) q[3];
sx q[3];
rz(-1.6955175) q[3];
sx q[3];
rz(0.41984841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3495425) q[2];
sx q[2];
rz(-0.011186102) q[2];
sx q[2];
rz(-2.1816317) q[2];
rz(-0.32944426) q[3];
sx q[3];
rz(-3.1335148) q[3];
sx q[3];
rz(-0.85028696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.195381) q[0];
sx q[0];
rz(-0.61310261) q[0];
sx q[0];
rz(-3.0323113) q[0];
rz(-2.7682313) q[1];
sx q[1];
rz(-0.80972087) q[1];
sx q[1];
rz(-1.2304617) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5143573) q[0];
sx q[0];
rz(-2.1684596) q[0];
sx q[0];
rz(2.3742832) q[0];
x q[1];
rz(2.9443342) q[2];
sx q[2];
rz(-1.3840809) q[2];
sx q[2];
rz(3.124539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8671678) q[1];
sx q[1];
rz(-1.503957) q[1];
sx q[1];
rz(3.1301366) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8318895) q[3];
sx q[3];
rz(-2.1983357) q[3];
sx q[3];
rz(-2.5554339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5665148) q[2];
sx q[2];
rz(-1.9058303) q[2];
sx q[2];
rz(-1.3285948) q[2];
rz(1.7447507) q[3];
sx q[3];
rz(-0.003740398) q[3];
sx q[3];
rz(-2.1197135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063865572) q[0];
sx q[0];
rz(-1.4430178) q[0];
sx q[0];
rz(-2.5685837) q[0];
rz(0.30814463) q[1];
sx q[1];
rz(-2.7317218) q[1];
sx q[1];
rz(-2.134197) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38978168) q[0];
sx q[0];
rz(-0.10669076) q[0];
sx q[0];
rz(-1.5185028) q[0];
rz(-pi) q[1];
rz(-0.043906004) q[2];
sx q[2];
rz(-1.708235) q[2];
sx q[2];
rz(3.1071747) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.17698174) q[1];
sx q[1];
rz(-1.4876502) q[1];
sx q[1];
rz(0.040772922) q[1];
x q[2];
rz(0.013434826) q[3];
sx q[3];
rz(-1.6009496) q[3];
sx q[3];
rz(2.9363901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.31743) q[2];
sx q[2];
rz(-0.62676668) q[2];
sx q[2];
rz(0.38995788) q[2];
rz(-0.069084875) q[3];
sx q[3];
rz(-0.0091113541) q[3];
sx q[3];
rz(0.33153427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020141715) q[0];
sx q[0];
rz(-0.75103253) q[0];
sx q[0];
rz(2.6556515) q[0];
rz(0.87156975) q[1];
sx q[1];
rz(-1.3078682) q[1];
sx q[1];
rz(1.647324) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6914929) q[0];
sx q[0];
rz(-2.1039824) q[0];
sx q[0];
rz(1.9271118) q[0];
rz(-pi) q[1];
rz(2.5267692) q[2];
sx q[2];
rz(-1.6069001) q[2];
sx q[2];
rz(-1.6143798) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4758265) q[1];
sx q[1];
rz(-1.2428871) q[1];
sx q[1];
rz(1.2535415) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.028454) q[3];
sx q[3];
rz(-1.5290878) q[3];
sx q[3];
rz(-2.0737518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5668874) q[2];
sx q[2];
rz(-0.042782728) q[2];
sx q[2];
rz(-0.031878397) q[2];
rz(2.3632862) q[3];
sx q[3];
rz(-3.1347771) q[3];
sx q[3];
rz(2.8460898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42302172) q[0];
sx q[0];
rz(-1.6091249) q[0];
sx q[0];
rz(-1.3269497) q[0];
rz(-3.0146535) q[1];
sx q[1];
rz(-0.23902421) q[1];
sx q[1];
rz(0.21993266) q[1];
rz(1.610582) q[2];
sx q[2];
rz(-3.0018158) q[2];
sx q[2];
rz(-2.9374585) q[2];
rz(-1.4999785) q[3];
sx q[3];
rz(-0.67588617) q[3];
sx q[3];
rz(-2.5839154) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

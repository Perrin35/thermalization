OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.863707) q[0];
sx q[0];
rz(-0.42930332) q[0];
sx q[0];
rz(-2.8306146) q[0];
rz(-1.3485981) q[1];
sx q[1];
rz(-2.0452979) q[1];
sx q[1];
rz(-0.57408339) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9756843) q[0];
sx q[0];
rz(-0.57281369) q[0];
sx q[0];
rz(-2.2982135) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98784222) q[2];
sx q[2];
rz(-1.5997837) q[2];
sx q[2];
rz(1.0430973) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8771389) q[1];
sx q[1];
rz(-1.2065556) q[1];
sx q[1];
rz(3.0458782) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9472136) q[3];
sx q[3];
rz(-1.6302846) q[3];
sx q[3];
rz(0.74947442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8671888) q[2];
sx q[2];
rz(-1.8827266) q[2];
sx q[2];
rz(-1.3192419) q[2];
rz(0.073171767) q[3];
sx q[3];
rz(-1.3061482) q[3];
sx q[3];
rz(2.3225972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0935593) q[0];
sx q[0];
rz(-1.9608542) q[0];
sx q[0];
rz(0.67980415) q[0];
rz(-0.80348429) q[1];
sx q[1];
rz(-1.7796703) q[1];
sx q[1];
rz(0.63978535) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6006192) q[0];
sx q[0];
rz(-2.2007211) q[0];
sx q[0];
rz(-0.39430228) q[0];
rz(2.4117208) q[2];
sx q[2];
rz(-1.6449071) q[2];
sx q[2];
rz(2.6593936) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3099932) q[1];
sx q[1];
rz(-0.91504708) q[1];
sx q[1];
rz(0.80226957) q[1];
rz(2.5211447) q[3];
sx q[3];
rz(-2.7683966) q[3];
sx q[3];
rz(-2.9377112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5378319) q[2];
sx q[2];
rz(-2.64309) q[2];
sx q[2];
rz(-0.70408386) q[2];
rz(-3.0294561) q[3];
sx q[3];
rz(-1.6928558) q[3];
sx q[3];
rz(-1.0224379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10244399) q[0];
sx q[0];
rz(-2.0344489) q[0];
sx q[0];
rz(-0.33945864) q[0];
rz(-2.0231694) q[1];
sx q[1];
rz(-2.2148841) q[1];
sx q[1];
rz(0.035645398) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6528856) q[0];
sx q[0];
rz(-1.4849326) q[0];
sx q[0];
rz(-2.9339132) q[0];
x q[1];
rz(-0.33907922) q[2];
sx q[2];
rz(-1.2847752) q[2];
sx q[2];
rz(1.3468483) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0633723) q[1];
sx q[1];
rz(-2.2148892) q[1];
sx q[1];
rz(-1.3011342) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6142637) q[3];
sx q[3];
rz(-1.8395367) q[3];
sx q[3];
rz(0.72878557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.79225916) q[2];
sx q[2];
rz(-0.68816853) q[2];
sx q[2];
rz(-2.2815857) q[2];
rz(-2.3490014) q[3];
sx q[3];
rz(-2.9383797) q[3];
sx q[3];
rz(1.0205166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.129313) q[0];
sx q[0];
rz(-1.7565933) q[0];
sx q[0];
rz(0.63283515) q[0];
rz(1.1526147) q[1];
sx q[1];
rz(-2.2668138) q[1];
sx q[1];
rz(1.6713743) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8058384) q[0];
sx q[0];
rz(-0.64941657) q[0];
sx q[0];
rz(1.592167) q[0];
x q[1];
rz(0.15369065) q[2];
sx q[2];
rz(-0.77506232) q[2];
sx q[2];
rz(-2.3714921) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8006261) q[1];
sx q[1];
rz(-1.3142085) q[1];
sx q[1];
rz(2.2881195) q[1];
rz(-pi) q[2];
rz(-0.95872236) q[3];
sx q[3];
rz(-2.1321724) q[3];
sx q[3];
rz(1.4179729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4319438) q[2];
sx q[2];
rz(-2.180763) q[2];
sx q[2];
rz(0.43080899) q[2];
rz(-2.4591947) q[3];
sx q[3];
rz(-1.4902481) q[3];
sx q[3];
rz(-2.1916913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34669852) q[0];
sx q[0];
rz(-1.2656724) q[0];
sx q[0];
rz(1.9516113) q[0];
rz(2.8311912) q[1];
sx q[1];
rz(-1.0810532) q[1];
sx q[1];
rz(2.6365872) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81722084) q[0];
sx q[0];
rz(-1.1341057) q[0];
sx q[0];
rz(-0.94193108) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1558542) q[2];
sx q[2];
rz(-1.4507074) q[2];
sx q[2];
rz(-1.1121554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68143916) q[1];
sx q[1];
rz(-2.6362754) q[1];
sx q[1];
rz(2.9605335) q[1];
rz(-pi) q[2];
rz(1.4851863) q[3];
sx q[3];
rz(-2.2946649) q[3];
sx q[3];
rz(1.6023265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40671047) q[2];
sx q[2];
rz(-2.4411185) q[2];
sx q[2];
rz(-0.90275466) q[2];
rz(1.0550176) q[3];
sx q[3];
rz(-0.86692923) q[3];
sx q[3];
rz(0.46190754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2257253) q[0];
sx q[0];
rz(-0.71284717) q[0];
sx q[0];
rz(2.0676887) q[0];
rz(1.7621) q[1];
sx q[1];
rz(-2.4122489) q[1];
sx q[1];
rz(-0.097537907) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3341971) q[0];
sx q[0];
rz(-2.1046035) q[0];
sx q[0];
rz(-1.7814723) q[0];
x q[1];
rz(0.16867192) q[2];
sx q[2];
rz(-1.9430117) q[2];
sx q[2];
rz(1.8255359) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8733858) q[1];
sx q[1];
rz(-1.6318351) q[1];
sx q[1];
rz(3.0931506) q[1];
rz(-pi) q[2];
x q[2];
rz(2.318368) q[3];
sx q[3];
rz(-0.71079094) q[3];
sx q[3];
rz(-1.2430134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0292616) q[2];
sx q[2];
rz(-1.4533726) q[2];
sx q[2];
rz(-2.7304999) q[2];
rz(1.7279651) q[3];
sx q[3];
rz(-0.77815431) q[3];
sx q[3];
rz(1.5948064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3635062) q[0];
sx q[0];
rz(-3.111105) q[0];
sx q[0];
rz(1.7331069) q[0];
rz(2.1259437) q[1];
sx q[1];
rz(-1.3280832) q[1];
sx q[1];
rz(-1.4782864) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41970872) q[0];
sx q[0];
rz(-0.96707487) q[0];
sx q[0];
rz(-1.2092071) q[0];
rz(2.4934216) q[2];
sx q[2];
rz(-1.0888466) q[2];
sx q[2];
rz(-0.75343695) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1455545) q[1];
sx q[1];
rz(-1.1281456) q[1];
sx q[1];
rz(3.0575322) q[1];
rz(-pi) q[2];
rz(1.3331378) q[3];
sx q[3];
rz(-0.51222748) q[3];
sx q[3];
rz(-2.81306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0279072) q[2];
sx q[2];
rz(-2.5137641) q[2];
sx q[2];
rz(2.9713463) q[2];
rz(1.2804735) q[3];
sx q[3];
rz(-1.6672641) q[3];
sx q[3];
rz(1.1580275) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4729507) q[0];
sx q[0];
rz(-2.1773715) q[0];
sx q[0];
rz(2.6131795) q[0];
rz(2.2083185) q[1];
sx q[1];
rz(-2.2163138) q[1];
sx q[1];
rz(2.5859213) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9708722) q[0];
sx q[0];
rz(-1.375388) q[0];
sx q[0];
rz(0.94265818) q[0];
rz(-pi) q[1];
rz(-1.6711177) q[2];
sx q[2];
rz(-2.360376) q[2];
sx q[2];
rz(2.2364716) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.96413999) q[1];
sx q[1];
rz(-2.429306) q[1];
sx q[1];
rz(-3.0925372) q[1];
x q[2];
rz(-2.7714588) q[3];
sx q[3];
rz(-2.4151093) q[3];
sx q[3];
rz(-0.47011061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.8133424) q[2];
sx q[2];
rz(-2.3395061) q[2];
sx q[2];
rz(-0.30725202) q[2];
rz(-2.3452289) q[3];
sx q[3];
rz(-2.4107404) q[3];
sx q[3];
rz(-3.0439607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55427134) q[0];
sx q[0];
rz(-2.0926496) q[0];
sx q[0];
rz(-1.3294504) q[0];
rz(0.21996552) q[1];
sx q[1];
rz(-2.797762) q[1];
sx q[1];
rz(1.8230009) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7466965) q[0];
sx q[0];
rz(-1.3550955) q[0];
sx q[0];
rz(-1.8442276) q[0];
rz(-pi) q[1];
rz(-0.49974738) q[2];
sx q[2];
rz(-1.6732273) q[2];
sx q[2];
rz(0.62291716) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.66178759) q[1];
sx q[1];
rz(-2.0347682) q[1];
sx q[1];
rz(-2.7590092) q[1];
rz(-0.72033003) q[3];
sx q[3];
rz(-1.2033278) q[3];
sx q[3];
rz(-3.0331963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69507504) q[2];
sx q[2];
rz(-1.4277642) q[2];
sx q[2];
rz(1.8565149) q[2];
rz(2.255693) q[3];
sx q[3];
rz(-1.748184) q[3];
sx q[3];
rz(-1.5281965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1000243) q[0];
sx q[0];
rz(-2.4055241) q[0];
sx q[0];
rz(2.8247483) q[0];
rz(-0.44081229) q[1];
sx q[1];
rz(-2.3471954) q[1];
sx q[1];
rz(-0.37030927) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49476162) q[0];
sx q[0];
rz(-1.6017822) q[0];
sx q[0];
rz(1.9520743) q[0];
x q[1];
rz(-2.4914527) q[2];
sx q[2];
rz(-0.95528379) q[2];
sx q[2];
rz(2.5742449) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1052484) q[1];
sx q[1];
rz(-0.40648983) q[1];
sx q[1];
rz(2.2986733) q[1];
rz(-0.050906128) q[3];
sx q[3];
rz(-1.6744782) q[3];
sx q[3];
rz(0.10242505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73654282) q[2];
sx q[2];
rz(-0.49095792) q[2];
sx q[2];
rz(-1.6323818) q[2];
rz(-0.80471188) q[3];
sx q[3];
rz(-1.5038306) q[3];
sx q[3];
rz(-3.0338083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1152773) q[0];
sx q[0];
rz(-1.7290709) q[0];
sx q[0];
rz(1.2363634) q[0];
rz(-2.9112877) q[1];
sx q[1];
rz(-0.31619148) q[1];
sx q[1];
rz(-2.2264623) q[1];
rz(-0.98501803) q[2];
sx q[2];
rz(-2.14902) q[2];
sx q[2];
rz(2.4075748) q[2];
rz(-0.95190081) q[3];
sx q[3];
rz(-1.8397844) q[3];
sx q[3];
rz(2.7927273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

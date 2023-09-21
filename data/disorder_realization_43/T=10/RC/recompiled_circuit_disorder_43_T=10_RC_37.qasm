OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(1.8068846) q[0];
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(-2.1656353) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87603509) q[0];
sx q[0];
rz(-1.3296488) q[0];
sx q[0];
rz(0.59224706) q[0];
rz(-1.0693597) q[2];
sx q[2];
rz(-0.62383365) q[2];
sx q[2];
rz(1.8471579) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9589899) q[1];
sx q[1];
rz(-0.62287736) q[1];
sx q[1];
rz(1.3401003) q[1];
rz(2.7482412) q[3];
sx q[3];
rz(-0.8439807) q[3];
sx q[3];
rz(-0.8918744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7261937) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(0.78380084) q[2];
rz(0.33256724) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(1.3403085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6587104) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(0.17856199) q[0];
rz(-1.3372955) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(3.1352502) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3865249) q[0];
sx q[0];
rz(-1.7641983) q[0];
sx q[0];
rz(-0.13271876) q[0];
rz(-0.86654051) q[2];
sx q[2];
rz(-2.4079977) q[2];
sx q[2];
rz(1.3726335) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9176863) q[1];
sx q[1];
rz(-2.8250541) q[1];
sx q[1];
rz(0.19885893) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7553495) q[3];
sx q[3];
rz(-0.59730232) q[3];
sx q[3];
rz(2.951705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1938842) q[2];
sx q[2];
rz(-0.42703736) q[2];
sx q[2];
rz(2.2634899) q[2];
rz(2.2795423) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(-3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5851615) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(0.01097824) q[0];
rz(0.36704656) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(3.045851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11287963) q[0];
sx q[0];
rz(-2.9071147) q[0];
sx q[0];
rz(-2.2857091) q[0];
rz(-pi) q[1];
rz(2.0289621) q[2];
sx q[2];
rz(-1.1095699) q[2];
sx q[2];
rz(-2.0756276) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2925551) q[1];
sx q[1];
rz(-1.109086) q[1];
sx q[1];
rz(-2.1828116) q[1];
rz(-pi) q[2];
rz(-1.5355574) q[3];
sx q[3];
rz(-1.0198776) q[3];
sx q[3];
rz(0.86722022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0720955) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(0.83479184) q[2];
rz(-2.9299724) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(0.37477469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19462207) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(0.34657493) q[0];
rz(-2.6158781) q[1];
sx q[1];
rz(-2.321967) q[1];
sx q[1];
rz(-2.1077164) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0948807) q[0];
sx q[0];
rz(-2.046642) q[0];
sx q[0];
rz(-2.4701719) q[0];
x q[1];
rz(0.526555) q[2];
sx q[2];
rz(-1.0609846) q[2];
sx q[2];
rz(1.9499792) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.66400601) q[1];
sx q[1];
rz(-0.2020745) q[1];
sx q[1];
rz(-1.2408153) q[1];
rz(2.8598966) q[3];
sx q[3];
rz(-1.2558189) q[3];
sx q[3];
rz(-1.7236934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6870849) q[2];
sx q[2];
rz(-1.4322174) q[2];
sx q[2];
rz(1.3467849) q[2];
rz(2.7205617) q[3];
sx q[3];
rz(-1.0166758) q[3];
sx q[3];
rz(-2.4436061) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7088257) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(2.72686) q[0];
rz(-1.3955836) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(2.5674852) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610883) q[0];
sx q[0];
rz(-1.442712) q[0];
sx q[0];
rz(-2.1696521) q[0];
rz(0.026002361) q[2];
sx q[2];
rz(-0.71547316) q[2];
sx q[2];
rz(2.0040087) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9088604) q[1];
sx q[1];
rz(-2.6870011) q[1];
sx q[1];
rz(0.67635398) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19458171) q[3];
sx q[3];
rz(-2.2268725) q[3];
sx q[3];
rz(1.7723099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.85270143) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(1.5138907) q[2];
rz(-3.1001575) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067327499) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(-2.4434027) q[0];
rz(-0.16695887) q[1];
sx q[1];
rz(-2.0623465) q[1];
sx q[1];
rz(2.4093157) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0064272881) q[0];
sx q[0];
rz(-0.26627243) q[0];
sx q[0];
rz(-0.64384319) q[0];
x q[1];
rz(-0.68533021) q[2];
sx q[2];
rz(-2.2517859) q[2];
sx q[2];
rz(2.6076917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.60019833) q[1];
sx q[1];
rz(-1.3458034) q[1];
sx q[1];
rz(0.66585983) q[1];
rz(0.5703339) q[3];
sx q[3];
rz(-1.7551646) q[3];
sx q[3];
rz(1.0384699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7531062) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(-1.9801271) q[2];
rz(0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(-1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56918615) q[0];
sx q[0];
rz(-0.82018667) q[0];
sx q[0];
rz(-2.9851595) q[0];
rz(0.45174831) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(2.3715473) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2609445) q[0];
sx q[0];
rz(-2.5136607) q[0];
sx q[0];
rz(-0.29675608) q[0];
rz(-pi) q[1];
rz(-2.7355843) q[2];
sx q[2];
rz(-1.0205262) q[2];
sx q[2];
rz(-1.6183491) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5321977) q[1];
sx q[1];
rz(-0.97071338) q[1];
sx q[1];
rz(2.1983912) q[1];
rz(1.2392063) q[3];
sx q[3];
rz(-1.0911897) q[3];
sx q[3];
rz(0.94707205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6440755) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(-2.1264123) q[2];
rz(1.4366359) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(-0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5148233) q[0];
sx q[0];
rz(-0.55861449) q[0];
sx q[0];
rz(2.8998937) q[0];
rz(2.4027951) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(-2.7752005) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.217427) q[0];
sx q[0];
rz(-2.1806742) q[0];
sx q[0];
rz(0.92745552) q[0];
rz(-pi) q[1];
rz(2.0051458) q[2];
sx q[2];
rz(-0.98698101) q[2];
sx q[2];
rz(0.39603147) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0802676) q[1];
sx q[1];
rz(-0.72239164) q[1];
sx q[1];
rz(0.25218833) q[1];
rz(-1.0189692) q[3];
sx q[3];
rz(-1.5295715) q[3];
sx q[3];
rz(3.105203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1239132) q[2];
sx q[2];
rz(-1.9423449) q[2];
sx q[2];
rz(0.29385847) q[2];
rz(3.1270694) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83207399) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(2.2576387) q[0];
rz(-0.66028315) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(-2.3502137) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6521097) q[0];
sx q[0];
rz(-1.8882505) q[0];
sx q[0];
rz(2.2788458) q[0];
rz(-pi) q[1];
rz(2.3320138) q[2];
sx q[2];
rz(-2.5123345) q[2];
sx q[2];
rz(1.4776243) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8225122) q[1];
sx q[1];
rz(-1.984239) q[1];
sx q[1];
rz(2.6507069) q[1];
x q[2];
rz(2.1095554) q[3];
sx q[3];
rz(-1.4269392) q[3];
sx q[3];
rz(-2.0335576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0749977) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(-1.0212612) q[2];
rz(-2.7630473) q[3];
sx q[3];
rz(-2.5773541) q[3];
sx q[3];
rz(2.5436201) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026697712) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(0.60780203) q[0];
rz(-0.23884493) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(-1.3806608) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5138807) q[0];
sx q[0];
rz(-2.4988334) q[0];
sx q[0];
rz(1.7895262) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1157687) q[2];
sx q[2];
rz(-1.1354453) q[2];
sx q[2];
rz(-3.0961852) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.41439357) q[1];
sx q[1];
rz(-0.52800035) q[1];
sx q[1];
rz(-2.2382733) q[1];
rz(-pi) q[2];
rz(-0.18908432) q[3];
sx q[3];
rz(-1.4049279) q[3];
sx q[3];
rz(-1.3454816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3537139) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(2.805368) q[2];
rz(-2.0119038) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(-2.0000134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1375785) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(-0.54429383) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(2.8019194) q[2];
sx q[2];
rz(-2.2876231) q[2];
sx q[2];
rz(-3.0291578) q[2];
rz(1.0014793) q[3];
sx q[3];
rz(-1.3712728) q[3];
sx q[3];
rz(-1.2029592) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

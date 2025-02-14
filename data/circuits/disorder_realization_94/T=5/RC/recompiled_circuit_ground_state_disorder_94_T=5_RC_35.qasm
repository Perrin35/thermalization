OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2151467) q[0];
sx q[0];
rz(-2.4617221) q[0];
sx q[0];
rz(0.1006861) q[0];
rz(-2.6868532) q[1];
sx q[1];
rz(-2.4603031) q[1];
sx q[1];
rz(-2.1159621) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68368841) q[0];
sx q[0];
rz(-1.6207028) q[0];
sx q[0];
rz(-1.3752974) q[0];
rz(0.16029393) q[2];
sx q[2];
rz(-2.829814) q[2];
sx q[2];
rz(-0.82138731) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4512349) q[1];
sx q[1];
rz(-0.75690311) q[1];
sx q[1];
rz(-0.1166269) q[1];
rz(2.1061534) q[3];
sx q[3];
rz(-2.1637754) q[3];
sx q[3];
rz(1.6361332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92496282) q[2];
sx q[2];
rz(-1.7731885) q[2];
sx q[2];
rz(-1.487757) q[2];
rz(1.9803068) q[3];
sx q[3];
rz(-1.0043283) q[3];
sx q[3];
rz(1.0950834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0682812) q[0];
sx q[0];
rz(-2.735205) q[0];
sx q[0];
rz(0.43346369) q[0];
rz(1.06217) q[1];
sx q[1];
rz(-1.559161) q[1];
sx q[1];
rz(2.4210222) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67340785) q[0];
sx q[0];
rz(-1.5338729) q[0];
sx q[0];
rz(1.4440749) q[0];
x q[1];
rz(2.8938204) q[2];
sx q[2];
rz(-1.0508014) q[2];
sx q[2];
rz(1.8682075) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51746619) q[1];
sx q[1];
rz(-1.4223012) q[1];
sx q[1];
rz(-0.35617574) q[1];
rz(-pi) q[2];
rz(-2.8470542) q[3];
sx q[3];
rz(-1.1014767) q[3];
sx q[3];
rz(2.81463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7945107) q[2];
sx q[2];
rz(-1.9158659) q[2];
sx q[2];
rz(0.98428717) q[2];
rz(2.2847564) q[3];
sx q[3];
rz(-2.555116) q[3];
sx q[3];
rz(1.2543359) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03265753) q[0];
sx q[0];
rz(-0.28626838) q[0];
sx q[0];
rz(2.3987067) q[0];
rz(-0.0097097857) q[1];
sx q[1];
rz(-1.5747036) q[1];
sx q[1];
rz(0.28365338) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8946202) q[0];
sx q[0];
rz(-2.5067911) q[0];
sx q[0];
rz(-2.532124) q[0];
rz(-2.9417324) q[2];
sx q[2];
rz(-1.7283354) q[2];
sx q[2];
rz(-3.1197366) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4176826) q[1];
sx q[1];
rz(-1.7387448) q[1];
sx q[1];
rz(-1.8803094) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1183692) q[3];
sx q[3];
rz(-1.4693714) q[3];
sx q[3];
rz(-1.2161127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.11188406) q[2];
sx q[2];
rz(-1.4446898) q[2];
sx q[2];
rz(-2.4857944) q[2];
rz(-0.2119952) q[3];
sx q[3];
rz(-2.5097804) q[3];
sx q[3];
rz(1.1727715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1126605) q[0];
sx q[0];
rz(-0.49429587) q[0];
sx q[0];
rz(-1.5869045) q[0];
rz(-2.9117865) q[1];
sx q[1];
rz(-2.4202012) q[1];
sx q[1];
rz(-1.302964) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6617216) q[0];
sx q[0];
rz(-0.73307836) q[0];
sx q[0];
rz(-2.4443611) q[0];
rz(-1.8660061) q[2];
sx q[2];
rz(-1.8285995) q[2];
sx q[2];
rz(2.757673) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11889549) q[1];
sx q[1];
rz(-0.52403421) q[1];
sx q[1];
rz(-0.043020821) q[1];
rz(1.0423562) q[3];
sx q[3];
rz(-2.5109249) q[3];
sx q[3];
rz(-1.3321277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0064696781) q[2];
sx q[2];
rz(-2.3389356) q[2];
sx q[2];
rz(0.50372493) q[2];
rz(1.6040365) q[3];
sx q[3];
rz(-1.9683014) q[3];
sx q[3];
rz(1.5776207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2876494) q[0];
sx q[0];
rz(-0.33677736) q[0];
sx q[0];
rz(2.7119998) q[0];
rz(-2.8751539) q[1];
sx q[1];
rz(-2.3315505) q[1];
sx q[1];
rz(1.0688759) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42035139) q[0];
sx q[0];
rz(-1.1614033) q[0];
sx q[0];
rz(-0.73914011) q[0];
rz(-1.2524302) q[2];
sx q[2];
rz(-2.3471544) q[2];
sx q[2];
rz(-0.81046644) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4539448) q[1];
sx q[1];
rz(-1.3140716) q[1];
sx q[1];
rz(-2.686108) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3232628) q[3];
sx q[3];
rz(-1.8283852) q[3];
sx q[3];
rz(-0.12572337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2499007) q[2];
sx q[2];
rz(-0.72600681) q[2];
sx q[2];
rz(-0.1740087) q[2];
rz(-1.8788667) q[3];
sx q[3];
rz(-1.5401253) q[3];
sx q[3];
rz(-1.5866097) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1332909) q[0];
sx q[0];
rz(-1.0187848) q[0];
sx q[0];
rz(-2.2513576) q[0];
rz(1.7091735) q[1];
sx q[1];
rz(-2.1232257) q[1];
sx q[1];
rz(-3.0216246) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1748808) q[0];
sx q[0];
rz(-2.565756) q[0];
sx q[0];
rz(-2.8737646) q[0];
x q[1];
rz(-0.34949686) q[2];
sx q[2];
rz(-1.1634524) q[2];
sx q[2];
rz(1.1098338) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.240472) q[1];
sx q[1];
rz(-1.8163306) q[1];
sx q[1];
rz(-2.4811109) q[1];
rz(0.31474416) q[3];
sx q[3];
rz(-2.5786244) q[3];
sx q[3];
rz(1.1551577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1098108) q[2];
sx q[2];
rz(-1.6620771) q[2];
sx q[2];
rz(-0.39153448) q[2];
rz(-1.7560962) q[3];
sx q[3];
rz(-0.23844312) q[3];
sx q[3];
rz(2.5029206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2077797) q[0];
sx q[0];
rz(-2.5901828) q[0];
sx q[0];
rz(-1.3054003) q[0];
rz(-2.6370559) q[1];
sx q[1];
rz(-1.7326109) q[1];
sx q[1];
rz(-0.22444589) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87418694) q[0];
sx q[0];
rz(-1.5874369) q[0];
sx q[0];
rz(1.5550625) q[0];
rz(-pi) q[1];
rz(2.6487892) q[2];
sx q[2];
rz(-1.7597919) q[2];
sx q[2];
rz(-2.1616621) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1149772) q[1];
sx q[1];
rz(-1.3040286) q[1];
sx q[1];
rz(1.90735) q[1];
rz(-pi) q[2];
x q[2];
rz(2.283758) q[3];
sx q[3];
rz(-2.668475) q[3];
sx q[3];
rz(2.3046062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9597783) q[2];
sx q[2];
rz(-1.10428) q[2];
sx q[2];
rz(-2.1941198) q[2];
rz(0.13253658) q[3];
sx q[3];
rz(-2.4819076) q[3];
sx q[3];
rz(0.98602492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9082311) q[0];
sx q[0];
rz(-0.34603226) q[0];
sx q[0];
rz(-0.37286266) q[0];
rz(2.8064959) q[1];
sx q[1];
rz(-1.9357888) q[1];
sx q[1];
rz(1.0710517) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6823387) q[0];
sx q[0];
rz(-1.5821348) q[0];
sx q[0];
rz(3.0724694) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55166371) q[2];
sx q[2];
rz(-1.7529738) q[2];
sx q[2];
rz(-2.7523486) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4651686) q[1];
sx q[1];
rz(-1.81552) q[1];
sx q[1];
rz(1.3307491) q[1];
rz(-3.0464977) q[3];
sx q[3];
rz(-1.250306) q[3];
sx q[3];
rz(1.2456428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4272473) q[2];
sx q[2];
rz(-1.599396) q[2];
sx q[2];
rz(-2.074312) q[2];
rz(-2.769477) q[3];
sx q[3];
rz(-2.1477063) q[3];
sx q[3];
rz(-0.20017643) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2315955) q[0];
sx q[0];
rz(-2.5439926) q[0];
sx q[0];
rz(1.4071314) q[0];
rz(-2.097997) q[1];
sx q[1];
rz(-2.2035972) q[1];
sx q[1];
rz(-0.49120206) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5564767) q[0];
sx q[0];
rz(-2.399754) q[0];
sx q[0];
rz(1.1354574) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7411072) q[2];
sx q[2];
rz(-1.7372157) q[2];
sx q[2];
rz(-0.8569878) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8206827) q[1];
sx q[1];
rz(-1.9100128) q[1];
sx q[1];
rz(1.7103819) q[1];
x q[2];
rz(1.3902499) q[3];
sx q[3];
rz(-1.9420764) q[3];
sx q[3];
rz(-2.2425687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4444943) q[2];
sx q[2];
rz(-1.5105057) q[2];
sx q[2];
rz(2.0972882) q[2];
rz(0.77110428) q[3];
sx q[3];
rz(-1.5325357) q[3];
sx q[3];
rz(-0.3573629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9412909) q[0];
sx q[0];
rz(-1.2933949) q[0];
sx q[0];
rz(1.6218761) q[0];
rz(1.3007851) q[1];
sx q[1];
rz(-1.3011353) q[1];
sx q[1];
rz(-2.4470952) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8845399) q[0];
sx q[0];
rz(-1.5844795) q[0];
sx q[0];
rz(-1.3428997) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6171309) q[2];
sx q[2];
rz(-1.6860644) q[2];
sx q[2];
rz(1.0712717) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.84811454) q[1];
sx q[1];
rz(-1.9047996) q[1];
sx q[1];
rz(-3.1298742) q[1];
rz(-pi) q[2];
rz(2.7312134) q[3];
sx q[3];
rz(-0.79550064) q[3];
sx q[3];
rz(-2.8727227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0684315) q[2];
sx q[2];
rz(-1.7794926) q[2];
sx q[2];
rz(2.6279602) q[2];
rz(-2.8237776) q[3];
sx q[3];
rz(-1.5376667) q[3];
sx q[3];
rz(-2.1224799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1401055) q[0];
sx q[0];
rz(-1.6248063) q[0];
sx q[0];
rz(-0.90909062) q[0];
rz(-1.9152676) q[1];
sx q[1];
rz(-1.5227804) q[1];
sx q[1];
rz(-2.7979122) q[1];
rz(-2.4438159) q[2];
sx q[2];
rz(-1.1502167) q[2];
sx q[2];
rz(0.95982734) q[2];
rz(-1.5671586) q[3];
sx q[3];
rz(-0.52593492) q[3];
sx q[3];
rz(1.8686663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

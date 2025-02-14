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
rz(-3.0946331) q[0];
sx q[0];
rz(-1.5927915) q[0];
sx q[0];
rz(-1.3943075) q[0];
rz(2.6810763) q[1];
sx q[1];
rz(-1.2682275) q[1];
sx q[1];
rz(0.4761129) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44350177) q[0];
sx q[0];
rz(-2.4558668) q[0];
sx q[0];
rz(-2.1954407) q[0];
rz(-pi) q[1];
rz(-1.3018487) q[2];
sx q[2];
rz(-1.7978872) q[2];
sx q[2];
rz(-0.64193945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75501498) q[1];
sx q[1];
rz(-1.5881032) q[1];
sx q[1];
rz(2.9933418) q[1];
rz(1.5067817) q[3];
sx q[3];
rz(-1.5435757) q[3];
sx q[3];
rz(-0.79827362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.32690471) q[2];
sx q[2];
rz(-1.2219595) q[2];
sx q[2];
rz(-1.8587221) q[2];
rz(3.0021216) q[3];
sx q[3];
rz(-1.8255511) q[3];
sx q[3];
rz(-2.4563346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4271456) q[0];
sx q[0];
rz(-2.7834263) q[0];
sx q[0];
rz(0.32592475) q[0];
rz(0.70949316) q[1];
sx q[1];
rz(-0.26115099) q[1];
sx q[1];
rz(2.479898) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70764953) q[0];
sx q[0];
rz(-3.1295332) q[0];
sx q[0];
rz(-0.90977408) q[0];
rz(-pi) q[1];
rz(2.6432132) q[2];
sx q[2];
rz(-2.1415258) q[2];
sx q[2];
rz(-0.98266593) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0584512) q[1];
sx q[1];
rz(-1.3399235) q[1];
sx q[1];
rz(-1.9564232) q[1];
rz(-pi) q[2];
rz(-2.4567408) q[3];
sx q[3];
rz(-1.0954482) q[3];
sx q[3];
rz(0.97739391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1208531) q[2];
sx q[2];
rz(-1.8031305) q[2];
sx q[2];
rz(0.2600812) q[2];
rz(-2.2777879) q[3];
sx q[3];
rz(-1.6983906) q[3];
sx q[3];
rz(1.9728194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9048555) q[0];
sx q[0];
rz(-2.3598292) q[0];
sx q[0];
rz(-2.8535063) q[0];
rz(1.9347363) q[1];
sx q[1];
rz(-2.0286045) q[1];
sx q[1];
rz(-0.77817121) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93302762) q[0];
sx q[0];
rz(-1.3817245) q[0];
sx q[0];
rz(-2.2651432) q[0];
rz(-pi) q[1];
rz(0.7581832) q[2];
sx q[2];
rz(-1.6155525) q[2];
sx q[2];
rz(0.54162177) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.4716123) q[1];
sx q[1];
rz(-1.1410574) q[1];
sx q[1];
rz(0.55008908) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2178389) q[3];
sx q[3];
rz(-1.2525932) q[3];
sx q[3];
rz(1.7585839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.244016) q[2];
sx q[2];
rz(-1.1268758) q[2];
sx q[2];
rz(0.24813949) q[2];
rz(2.1085619) q[3];
sx q[3];
rz(-1.6525729) q[3];
sx q[3];
rz(-1.8831133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7018062) q[0];
sx q[0];
rz(-0.93455625) q[0];
sx q[0];
rz(2.5820861) q[0];
rz(-1.7350381) q[1];
sx q[1];
rz(-2.1969257) q[1];
sx q[1];
rz(1.14538) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6481768) q[0];
sx q[0];
rz(-2.6828565) q[0];
sx q[0];
rz(-0.10347314) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4907095) q[2];
sx q[2];
rz(-2.0325629) q[2];
sx q[2];
rz(-2.8956763) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5694847) q[1];
sx q[1];
rz(-0.83781257) q[1];
sx q[1];
rz(-1.9721287) q[1];
rz(-3.1265852) q[3];
sx q[3];
rz(-1.0010442) q[3];
sx q[3];
rz(-0.80591312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.8765325) q[2];
sx q[2];
rz(-1.0141076) q[2];
sx q[2];
rz(1.8939135) q[2];
rz(1.0809336) q[3];
sx q[3];
rz(-1.388988) q[3];
sx q[3];
rz(-1.8814253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0697698) q[0];
sx q[0];
rz(-2.6520196) q[0];
sx q[0];
rz(2.3967632) q[0];
rz(1.7597594) q[1];
sx q[1];
rz(-1.1188743) q[1];
sx q[1];
rz(-0.098085731) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8902405) q[0];
sx q[0];
rz(-2.6664263) q[0];
sx q[0];
rz(2.422228) q[0];
rz(0.44386835) q[2];
sx q[2];
rz(-1.8763246) q[2];
sx q[2];
rz(-1.6084087) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.67691117) q[1];
sx q[1];
rz(-1.5891074) q[1];
sx q[1];
rz(0.69377884) q[1];
x q[2];
rz(1.7616055) q[3];
sx q[3];
rz(-2.6529783) q[3];
sx q[3];
rz(-3.0861175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.44094008) q[2];
sx q[2];
rz(-2.7589189) q[2];
sx q[2];
rz(-1.0901701) q[2];
rz(-2.8160461) q[3];
sx q[3];
rz(-1.5109589) q[3];
sx q[3];
rz(2.8946099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.9149822) q[0];
sx q[0];
rz(-0.53737265) q[0];
sx q[0];
rz(1.9052624) q[0];
rz(-2.5044598) q[1];
sx q[1];
rz(-0.56012154) q[1];
sx q[1];
rz(2.2083652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7819311) q[0];
sx q[0];
rz(-2.1225516) q[0];
sx q[0];
rz(-2.7130068) q[0];
rz(-1.1652214) q[2];
sx q[2];
rz(-0.45171192) q[2];
sx q[2];
rz(-3.0690985) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.98249032) q[1];
sx q[1];
rz(-1.4878927) q[1];
sx q[1];
rz(-0.29183951) q[1];
rz(2.0108443) q[3];
sx q[3];
rz(-0.55956105) q[3];
sx q[3];
rz(-0.51792972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8895662) q[2];
sx q[2];
rz(-1.7539975) q[2];
sx q[2];
rz(-1.8760366) q[2];
rz(1.7172074) q[3];
sx q[3];
rz(-2.6129183) q[3];
sx q[3];
rz(0.85844794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76765656) q[0];
sx q[0];
rz(-2.7601384) q[0];
sx q[0];
rz(2.9071627) q[0];
rz(-0.43342057) q[1];
sx q[1];
rz(-1.0973009) q[1];
sx q[1];
rz(-2.0030599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9397126) q[0];
sx q[0];
rz(-2.1576886) q[0];
sx q[0];
rz(0.71016772) q[0];
rz(-pi) q[1];
rz(1.578733) q[2];
sx q[2];
rz(-2.0614377) q[2];
sx q[2];
rz(-0.68665169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1416152) q[1];
sx q[1];
rz(-2.6385698) q[1];
sx q[1];
rz(-1.5161689) q[1];
rz(-pi) q[2];
rz(-0.23536162) q[3];
sx q[3];
rz(-0.70066626) q[3];
sx q[3];
rz(-1.9266025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66702691) q[2];
sx q[2];
rz(-0.67956769) q[2];
sx q[2];
rz(-0.95728528) q[2];
rz(2.3840617) q[3];
sx q[3];
rz(-1.4742955) q[3];
sx q[3];
rz(2.9981414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677419) q[0];
sx q[0];
rz(-2.0389281) q[0];
sx q[0];
rz(2.962501) q[0];
rz(2.9415019) q[1];
sx q[1];
rz(-2.4724019) q[1];
sx q[1];
rz(-0.77654138) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0514798) q[0];
sx q[0];
rz(-2.5830373) q[0];
sx q[0];
rz(1.9165975) q[0];
rz(-0.29751038) q[2];
sx q[2];
rz(-0.97504598) q[2];
sx q[2];
rz(1.7864625) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69570615) q[1];
sx q[1];
rz(-1.4974471) q[1];
sx q[1];
rz(0.35584764) q[1];
rz(-pi) q[2];
rz(0.051335585) q[3];
sx q[3];
rz(-1.253479) q[3];
sx q[3];
rz(2.7961344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17579235) q[2];
sx q[2];
rz(-2.0696023) q[2];
sx q[2];
rz(-1.3624066) q[2];
rz(-1.3029178) q[3];
sx q[3];
rz(-1.4785654) q[3];
sx q[3];
rz(-1.038507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4565444) q[0];
sx q[0];
rz(-1.1948723) q[0];
sx q[0];
rz(0.1142647) q[0];
rz(-1.9396797) q[1];
sx q[1];
rz(-0.54113954) q[1];
sx q[1];
rz(3.0466383) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33126011) q[0];
sx q[0];
rz(-0.55821943) q[0];
sx q[0];
rz(-1.2329739) q[0];
x q[1];
rz(0.46038119) q[2];
sx q[2];
rz(-2.2308439) q[2];
sx q[2];
rz(-0.4649064) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.11647955) q[1];
sx q[1];
rz(-2.1201061) q[1];
sx q[1];
rz(1.6366556) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41422959) q[3];
sx q[3];
rz(-1.0981907) q[3];
sx q[3];
rz(0.35652682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.3007043) q[2];
sx q[2];
rz(-1.7965094) q[2];
sx q[2];
rz(0.67187205) q[2];
rz(-2.4554409) q[3];
sx q[3];
rz(-0.66303623) q[3];
sx q[3];
rz(1.2175951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9753863) q[0];
sx q[0];
rz(-0.19151846) q[0];
sx q[0];
rz(1.6084877) q[0];
rz(0.32471049) q[1];
sx q[1];
rz(-1.8638116) q[1];
sx q[1];
rz(2.8285573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7164118) q[0];
sx q[0];
rz(-1.6508174) q[0];
sx q[0];
rz(-0.99786134) q[0];
x q[1];
rz(0.061278419) q[2];
sx q[2];
rz(-1.3512502) q[2];
sx q[2];
rz(-0.82876182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5165007) q[1];
sx q[1];
rz(-2.2269669) q[1];
sx q[1];
rz(0.8582127) q[1];
x q[2];
rz(0.5996941) q[3];
sx q[3];
rz(-0.71772493) q[3];
sx q[3];
rz(2.5120171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93840557) q[2];
sx q[2];
rz(-1.9878191) q[2];
sx q[2];
rz(-1.3333092) q[2];
rz(0.61776727) q[3];
sx q[3];
rz(-0.94527644) q[3];
sx q[3];
rz(1.1355404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4278605) q[0];
sx q[0];
rz(-1.8129616) q[0];
sx q[0];
rz(0.34995361) q[0];
rz(2.2433544) q[1];
sx q[1];
rz(-1.0844834) q[1];
sx q[1];
rz(-0.35379298) q[1];
rz(0.66302115) q[2];
sx q[2];
rz(-0.86227476) q[2];
sx q[2];
rz(-3.0515565) q[2];
rz(-3.114936) q[3];
sx q[3];
rz(-2.4116357) q[3];
sx q[3];
rz(0.38105376) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

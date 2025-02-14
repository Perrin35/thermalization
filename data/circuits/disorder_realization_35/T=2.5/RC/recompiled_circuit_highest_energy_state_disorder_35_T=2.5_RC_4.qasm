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
rz(-0.91312042) q[0];
sx q[0];
rz(-1.0364391) q[0];
sx q[0];
rz(1.6057462) q[0];
rz(-2.4802471) q[1];
sx q[1];
rz(-1.9547434) q[1];
sx q[1];
rz(-0.43651906) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52659863) q[0];
sx q[0];
rz(-0.33373555) q[0];
sx q[0];
rz(-1.2646535) q[0];
rz(2.5769325) q[2];
sx q[2];
rz(-1.1810978) q[2];
sx q[2];
rz(-2.1291358) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84387842) q[1];
sx q[1];
rz(-1.6116287) q[1];
sx q[1];
rz(0.86754129) q[1];
rz(-0.8869041) q[3];
sx q[3];
rz(-1.7085848) q[3];
sx q[3];
rz(-2.6215214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46301699) q[2];
sx q[2];
rz(-0.43872321) q[2];
sx q[2];
rz(-2.2185745) q[2];
rz(0.065741278) q[3];
sx q[3];
rz(-1.8140503) q[3];
sx q[3];
rz(-1.8008697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098669212) q[0];
sx q[0];
rz(-1.0634402) q[0];
sx q[0];
rz(0.03751066) q[0];
rz(-0.58049479) q[1];
sx q[1];
rz(-0.66941222) q[1];
sx q[1];
rz(0.69211778) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0789675) q[0];
sx q[0];
rz(-0.37165549) q[0];
sx q[0];
rz(-1.070648) q[0];
x q[1];
rz(1.0893482) q[2];
sx q[2];
rz(-2.5855067) q[2];
sx q[2];
rz(0.5769481) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.9034391) q[1];
sx q[1];
rz(-0.80793706) q[1];
sx q[1];
rz(2.417516) q[1];
rz(0.18499891) q[3];
sx q[3];
rz(-2.3810412) q[3];
sx q[3];
rz(0.669125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6633501) q[2];
sx q[2];
rz(-1.1967836) q[2];
sx q[2];
rz(-1.7459858) q[2];
rz(2.9546402) q[3];
sx q[3];
rz(-1.9713277) q[3];
sx q[3];
rz(-2.601534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4726987) q[0];
sx q[0];
rz(-2.5085594) q[0];
sx q[0];
rz(-2.582666) q[0];
rz(0.82866296) q[1];
sx q[1];
rz(-1.5907954) q[1];
sx q[1];
rz(2.8899946) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54788816) q[0];
sx q[0];
rz(-1.3798837) q[0];
sx q[0];
rz(-0.44661354) q[0];
rz(-pi) q[1];
rz(0.34464406) q[2];
sx q[2];
rz(-2.0953669) q[2];
sx q[2];
rz(-1.8328431) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.025574) q[1];
sx q[1];
rz(-1.3085828) q[1];
sx q[1];
rz(-0.88974726) q[1];
x q[2];
rz(2.4647242) q[3];
sx q[3];
rz(-2.1123611) q[3];
sx q[3];
rz(1.3071909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.805213) q[2];
sx q[2];
rz(-0.42625913) q[2];
sx q[2];
rz(-1.4317929) q[2];
rz(-0.21670565) q[3];
sx q[3];
rz(-1.6101937) q[3];
sx q[3];
rz(1.3592892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25877229) q[0];
sx q[0];
rz(-2.6348305) q[0];
sx q[0];
rz(-1.850542) q[0];
rz(0.82921118) q[1];
sx q[1];
rz(-2.0334838) q[1];
sx q[1];
rz(2.1270027) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9982581) q[0];
sx q[0];
rz(-0.59276544) q[0];
sx q[0];
rz(2.696585) q[0];
rz(-pi) q[1];
rz(-2.5780771) q[2];
sx q[2];
rz(-1.2402099) q[2];
sx q[2];
rz(-0.78621582) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9660849) q[1];
sx q[1];
rz(-2.1455742) q[1];
sx q[1];
rz(-0.0055343363) q[1];
rz(-pi) q[2];
rz(0.75506702) q[3];
sx q[3];
rz(-1.4658548) q[3];
sx q[3];
rz(-1.5615095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0097051) q[2];
sx q[2];
rz(-1.2010801) q[2];
sx q[2];
rz(0.76060549) q[2];
rz(2.6724114) q[3];
sx q[3];
rz(-1.6719336) q[3];
sx q[3];
rz(0.73452264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0483265) q[0];
sx q[0];
rz(-2.6272197) q[0];
sx q[0];
rz(2.560428) q[0];
rz(-1.5515074) q[1];
sx q[1];
rz(-2.4947512) q[1];
sx q[1];
rz(-2.4466628) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7030554) q[0];
sx q[0];
rz(-1.8037272) q[0];
sx q[0];
rz(1.2715879) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29232358) q[2];
sx q[2];
rz(-1.4806804) q[2];
sx q[2];
rz(-2.9499049) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0241922) q[1];
sx q[1];
rz(-1.0276737) q[1];
sx q[1];
rz(1.5195816) q[1];
rz(-0.23628037) q[3];
sx q[3];
rz(-1.7485011) q[3];
sx q[3];
rz(-1.1760933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.86357) q[2];
sx q[2];
rz(-1.5333971) q[2];
sx q[2];
rz(0.30280534) q[2];
rz(1.4263724) q[3];
sx q[3];
rz(-0.36077603) q[3];
sx q[3];
rz(0.58437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67702174) q[0];
sx q[0];
rz(-2.7320221) q[0];
sx q[0];
rz(-2.4466178) q[0];
rz(0.28680828) q[1];
sx q[1];
rz(-0.60568714) q[1];
sx q[1];
rz(1.8668176) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5999278) q[0];
sx q[0];
rz(-2.554691) q[0];
sx q[0];
rz(2.15564) q[0];
x q[1];
rz(-2.5996501) q[2];
sx q[2];
rz(-1.9795609) q[2];
sx q[2];
rz(-0.33318502) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.66836897) q[1];
sx q[1];
rz(-2.372207) q[1];
sx q[1];
rz(-1.046087) q[1];
x q[2];
rz(-0.58198317) q[3];
sx q[3];
rz(-1.3701539) q[3];
sx q[3];
rz(2.0394005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0810762) q[2];
sx q[2];
rz(-0.21848564) q[2];
sx q[2];
rz(2.1742353) q[2];
rz(-0.71990144) q[3];
sx q[3];
rz(-2.1981809) q[3];
sx q[3];
rz(1.8657743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8571781) q[0];
sx q[0];
rz(-0.030817742) q[0];
sx q[0];
rz(1.4469294) q[0];
rz(-2.6743496) q[1];
sx q[1];
rz(-1.2930361) q[1];
sx q[1];
rz(-1.1955059) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.104983) q[0];
sx q[0];
rz(-2.4901142) q[0];
sx q[0];
rz(2.3250513) q[0];
rz(2.325752) q[2];
sx q[2];
rz(-1.8058449) q[2];
sx q[2];
rz(1.5790958) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6536234) q[1];
sx q[1];
rz(-1.5088827) q[1];
sx q[1];
rz(2.0701522) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2310813) q[3];
sx q[3];
rz(-0.038148316) q[3];
sx q[3];
rz(-1.3919786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3168827) q[2];
sx q[2];
rz(-0.83493817) q[2];
sx q[2];
rz(-2.6173124) q[2];
rz(-2.8892062) q[3];
sx q[3];
rz(-0.10656825) q[3];
sx q[3];
rz(-3.0805123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97279945) q[0];
sx q[0];
rz(-2.0357098) q[0];
sx q[0];
rz(1.9148781) q[0];
rz(-2.3854158) q[1];
sx q[1];
rz(-0.94804472) q[1];
sx q[1];
rz(2.7511168) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7602063) q[0];
sx q[0];
rz(-1.1293656) q[0];
sx q[0];
rz(1.8618097) q[0];
rz(-3.0608589) q[2];
sx q[2];
rz(-2.6962523) q[2];
sx q[2];
rz(0.75917086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9041216) q[1];
sx q[1];
rz(-1.7842222) q[1];
sx q[1];
rz(-1.2280812) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.195695) q[3];
sx q[3];
rz(-2.3970553) q[3];
sx q[3];
rz(2.8685399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8172424) q[2];
sx q[2];
rz(-0.70049006) q[2];
sx q[2];
rz(-1.3520757) q[2];
rz(-2.9663441) q[3];
sx q[3];
rz(-1.8027571) q[3];
sx q[3];
rz(-1.2043183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7290203) q[0];
sx q[0];
rz(-0.46171284) q[0];
sx q[0];
rz(2.4793258) q[0];
rz(-2.5114255) q[1];
sx q[1];
rz(-1.8616734) q[1];
sx q[1];
rz(0.23552775) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2199687) q[0];
sx q[0];
rz(-1.7818091) q[0];
sx q[0];
rz(1.7485662) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9959683) q[2];
sx q[2];
rz(-2.5674985) q[2];
sx q[2];
rz(0.51806322) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3709066) q[1];
sx q[1];
rz(-1.5471493) q[1];
sx q[1];
rz(1.5093944) q[1];
x q[2];
rz(3.0074545) q[3];
sx q[3];
rz(-2.2629177) q[3];
sx q[3];
rz(1.1617171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0887289) q[2];
sx q[2];
rz(-0.95053089) q[2];
sx q[2];
rz(-2.5843184) q[2];
rz(0.82734621) q[3];
sx q[3];
rz(-1.6297623) q[3];
sx q[3];
rz(1.6415589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.92852229) q[0];
sx q[0];
rz(-0.7970354) q[0];
sx q[0];
rz(2.573977) q[0];
rz(-2.7396743) q[1];
sx q[1];
rz(-2.4979976) q[1];
sx q[1];
rz(-1.3622805) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6797736) q[0];
sx q[0];
rz(-1.1704233) q[0];
sx q[0];
rz(2.7072565) q[0];
rz(-pi) q[1];
rz(0.2593597) q[2];
sx q[2];
rz(-0.96732124) q[2];
sx q[2];
rz(0.063613907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.030008121) q[1];
sx q[1];
rz(-1.6492731) q[1];
sx q[1];
rz(2.2021738) q[1];
x q[2];
rz(2.3188842) q[3];
sx q[3];
rz(-0.97859446) q[3];
sx q[3];
rz(1.6975192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7207429) q[2];
sx q[2];
rz(-1.2030615) q[2];
sx q[2];
rz(2.8813349) q[2];
rz(2.4506954) q[3];
sx q[3];
rz(-0.8553718) q[3];
sx q[3];
rz(-1.5232085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33547587) q[0];
sx q[0];
rz(-1.5806883) q[0];
sx q[0];
rz(-1.6089532) q[0];
rz(-2.7083022) q[1];
sx q[1];
rz(-1.3229803) q[1];
sx q[1];
rz(1.9569474) q[1];
rz(-2.4337089) q[2];
sx q[2];
rz(-2.08692) q[2];
sx q[2];
rz(-0.89649123) q[2];
rz(-0.79037249) q[3];
sx q[3];
rz(-1.6675259) q[3];
sx q[3];
rz(-0.70601757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

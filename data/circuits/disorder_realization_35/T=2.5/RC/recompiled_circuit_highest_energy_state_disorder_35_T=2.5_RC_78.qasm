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
rz(0.66134557) q[1];
sx q[1];
rz(-1.1868492) q[1];
sx q[1];
rz(-2.7050736) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2921599) q[0];
sx q[0];
rz(-1.8884553) q[0];
sx q[0];
rz(3.0374797) q[0];
x q[1];
rz(1.1182734) q[2];
sx q[2];
rz(-1.0529537) q[2];
sx q[2];
rz(-2.81942) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4492717) q[1];
sx q[1];
rz(-2.2733449) q[1];
sx q[1];
rz(-0.0535123) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7868452) q[3];
sx q[3];
rz(-2.4461544) q[3];
sx q[3];
rz(1.2176633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6785757) q[2];
sx q[2];
rz(-0.43872321) q[2];
sx q[2];
rz(2.2185745) q[2];
rz(-0.065741278) q[3];
sx q[3];
rz(-1.3275423) q[3];
sx q[3];
rz(1.340723) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098669212) q[0];
sx q[0];
rz(-1.0634402) q[0];
sx q[0];
rz(0.03751066) q[0];
rz(-2.5610979) q[1];
sx q[1];
rz(-0.66941222) q[1];
sx q[1];
rz(2.4494749) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6094309) q[0];
sx q[0];
rz(-1.2464644) q[0];
sx q[0];
rz(-2.956809) q[0];
rz(2.0742998) q[2];
sx q[2];
rz(-1.3238591) q[2];
sx q[2];
rz(-1.4114789) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.925732) q[1];
sx q[1];
rz(-2.070148) q[1];
sx q[1];
rz(-2.4768922) q[1];
rz(-pi) q[2];
rz(1.7440657) q[3];
sx q[3];
rz(-0.82635802) q[3];
sx q[3];
rz(-2.2196774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.47824255) q[2];
sx q[2];
rz(-1.9448091) q[2];
sx q[2];
rz(1.7459858) q[2];
rz(-0.1869525) q[3];
sx q[3];
rz(-1.9713277) q[3];
sx q[3];
rz(-2.601534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66889399) q[0];
sx q[0];
rz(-0.63303328) q[0];
sx q[0];
rz(-0.5589267) q[0];
rz(0.82866296) q[1];
sx q[1];
rz(-1.5507973) q[1];
sx q[1];
rz(-2.8899946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64564451) q[0];
sx q[0];
rz(-0.48316607) q[0];
sx q[0];
rz(-0.42074411) q[0];
rz(-0.34464406) q[2];
sx q[2];
rz(-1.0462257) q[2];
sx q[2];
rz(-1.8328431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.14519599) q[1];
sx q[1];
rz(-2.4193978) q[1];
sx q[1];
rz(-1.1678371) q[1];
x q[2];
rz(-2.4647242) q[3];
sx q[3];
rz(-1.0292316) q[3];
sx q[3];
rz(1.3071909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.805213) q[2];
sx q[2];
rz(-2.7153335) q[2];
sx q[2];
rz(1.4317929) q[2];
rz(2.924887) q[3];
sx q[3];
rz(-1.5313989) q[3];
sx q[3];
rz(1.7823035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25877229) q[0];
sx q[0];
rz(-2.6348305) q[0];
sx q[0];
rz(-1.850542) q[0];
rz(-0.82921118) q[1];
sx q[1];
rz(-1.1081089) q[1];
sx q[1];
rz(2.1270027) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0907959) q[0];
sx q[0];
rz(-1.8136588) q[0];
sx q[0];
rz(-2.5953351) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9564187) q[2];
sx q[2];
rz(-2.1004538) q[2];
sx q[2];
rz(-0.5822863) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.17550771) q[1];
sx q[1];
rz(-2.1455742) q[1];
sx q[1];
rz(-0.0055343363) q[1];
rz(-pi) q[2];
rz(-0.15249522) q[3];
sx q[3];
rz(-0.76089459) q[3];
sx q[3];
rz(3.0214579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0097051) q[2];
sx q[2];
rz(-1.2010801) q[2];
sx q[2];
rz(-2.3809872) q[2];
rz(0.46918121) q[3];
sx q[3];
rz(-1.6719336) q[3];
sx q[3];
rz(-0.73452264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0483265) q[0];
sx q[0];
rz(-0.51437298) q[0];
sx q[0];
rz(0.58116466) q[0];
rz(-1.5900853) q[1];
sx q[1];
rz(-0.64684144) q[1];
sx q[1];
rz(0.69492984) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0804176) q[0];
sx q[0];
rz(-1.2799078) q[0];
sx q[0];
rz(0.24334749) q[0];
rz(2.8377436) q[2];
sx q[2];
rz(-2.8360747) q[2];
sx q[2];
rz(1.6696827) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1174005) q[1];
sx q[1];
rz(-2.1139189) q[1];
sx q[1];
rz(-1.5195816) q[1];
x q[2];
rz(-2.9053123) q[3];
sx q[3];
rz(-1.7485011) q[3];
sx q[3];
rz(1.1760933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.86357) q[2];
sx q[2];
rz(-1.5333971) q[2];
sx q[2];
rz(2.8387873) q[2];
rz(-1.4263724) q[3];
sx q[3];
rz(-2.7808166) q[3];
sx q[3];
rz(0.58437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.4645709) q[0];
sx q[0];
rz(-2.7320221) q[0];
sx q[0];
rz(2.4466178) q[0];
rz(0.28680828) q[1];
sx q[1];
rz(-0.60568714) q[1];
sx q[1];
rz(-1.274775) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6086084) q[0];
sx q[0];
rz(-1.8814981) q[0];
sx q[0];
rz(-2.0771189) q[0];
x q[1];
rz(0.69853492) q[2];
sx q[2];
rz(-0.66627995) q[2];
sx q[2];
rz(2.4874788) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4732237) q[1];
sx q[1];
rz(-2.372207) q[1];
sx q[1];
rz(1.046087) q[1];
rz(2.5596095) q[3];
sx q[3];
rz(-1.7714388) q[3];
sx q[3];
rz(-2.0394005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0810762) q[2];
sx q[2];
rz(-2.923107) q[2];
sx q[2];
rz(2.1742353) q[2];
rz(-2.4216912) q[3];
sx q[3];
rz(-2.1981809) q[3];
sx q[3];
rz(-1.8657743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8571781) q[0];
sx q[0];
rz(-3.1107749) q[0];
sx q[0];
rz(-1.6946633) q[0];
rz(-2.6743496) q[1];
sx q[1];
rz(-1.8485565) q[1];
sx q[1];
rz(-1.9460868) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9050483) q[0];
sx q[0];
rz(-1.1130739) q[0];
sx q[0];
rz(-0.48120705) q[0];
rz(-pi) q[1];
x q[1];
rz(2.325752) q[2];
sx q[2];
rz(-1.8058449) q[2];
sx q[2];
rz(1.5790958) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.48796921) q[1];
sx q[1];
rz(-1.6327099) q[1];
sx q[1];
rz(-2.0701522) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6067664) q[3];
sx q[3];
rz(-1.5580873) q[3];
sx q[3];
rz(-0.1606687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82470992) q[2];
sx q[2];
rz(-2.3066545) q[2];
sx q[2];
rz(2.6173124) q[2];
rz(-2.8892062) q[3];
sx q[3];
rz(-3.0350244) q[3];
sx q[3];
rz(3.0805123) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1687932) q[0];
sx q[0];
rz(-2.0357098) q[0];
sx q[0];
rz(-1.2267145) q[0];
rz(0.75617689) q[1];
sx q[1];
rz(-2.1935479) q[1];
sx q[1];
rz(-2.7511168) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0794546) q[0];
sx q[0];
rz(-1.833217) q[0];
sx q[0];
rz(2.6833437) q[0];
rz(0.080733732) q[2];
sx q[2];
rz(-2.6962523) q[2];
sx q[2];
rz(0.75917086) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23747102) q[1];
sx q[1];
rz(-1.3573705) q[1];
sx q[1];
rz(1.9135114) q[1];
rz(-0.94589767) q[3];
sx q[3];
rz(-2.3970553) q[3];
sx q[3];
rz(0.27305279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8172424) q[2];
sx q[2];
rz(-0.70049006) q[2];
sx q[2];
rz(1.7895169) q[2];
rz(2.9663441) q[3];
sx q[3];
rz(-1.3388355) q[3];
sx q[3];
rz(-1.2043183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7290203) q[0];
sx q[0];
rz(-2.6798798) q[0];
sx q[0];
rz(2.4793258) q[0];
rz(-2.5114255) q[1];
sx q[1];
rz(-1.8616734) q[1];
sx q[1];
rz(-2.9060649) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2199687) q[0];
sx q[0];
rz(-1.3597836) q[0];
sx q[0];
rz(-1.7485662) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6643737) q[2];
sx q[2];
rz(-2.1380599) q[2];
sx q[2];
rz(0.34514375) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9743226) q[1];
sx q[1];
rz(-0.065792699) q[1];
sx q[1];
rz(1.2029103) q[1];
rz(-pi) q[2];
rz(3.0074545) q[3];
sx q[3];
rz(-0.87867498) q[3];
sx q[3];
rz(1.9798756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0528637) q[2];
sx q[2];
rz(-2.1910618) q[2];
sx q[2];
rz(-0.55727422) q[2];
rz(2.3142464) q[3];
sx q[3];
rz(-1.5118303) q[3];
sx q[3];
rz(-1.5000337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2130704) q[0];
sx q[0];
rz(-0.7970354) q[0];
sx q[0];
rz(-2.573977) q[0];
rz(-2.7396743) q[1];
sx q[1];
rz(-2.4979976) q[1];
sx q[1];
rz(1.7793122) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80776725) q[0];
sx q[0];
rz(-2.5596715) q[0];
sx q[0];
rz(2.3533217) q[0];
x q[1];
rz(0.2593597) q[2];
sx q[2];
rz(-2.1742714) q[2];
sx q[2];
rz(3.0779787) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4940004) q[1];
sx q[1];
rz(-2.5060182) q[1];
sx q[1];
rz(1.7032436) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79093334) q[3];
sx q[3];
rz(-0.91703992) q[3];
sx q[3];
rz(2.4733651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.42084971) q[2];
sx q[2];
rz(-1.2030615) q[2];
sx q[2];
rz(2.8813349) q[2];
rz(-2.4506954) q[3];
sx q[3];
rz(-2.2862209) q[3];
sx q[3];
rz(-1.5232085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8061168) q[0];
sx q[0];
rz(-1.5609043) q[0];
sx q[0];
rz(1.5326395) q[0];
rz(-0.43329049) q[1];
sx q[1];
rz(-1.8186124) q[1];
sx q[1];
rz(-1.1846452) q[1];
rz(0.70788371) q[2];
sx q[2];
rz(-2.08692) q[2];
sx q[2];
rz(-0.89649123) q[2];
rz(2.3512202) q[3];
sx q[3];
rz(-1.6675259) q[3];
sx q[3];
rz(-0.70601757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

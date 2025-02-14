OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4545257) q[0];
sx q[0];
rz(-1.2626167) q[0];
sx q[0];
rz(0.73448056) q[0];
rz(1.6425411) q[1];
sx q[1];
rz(-1.2999111) q[1];
sx q[1];
rz(-2.1631961) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7534417) q[0];
sx q[0];
rz(-1.7072258) q[0];
sx q[0];
rz(-0.0078972422) q[0];
rz(-2.5402563) q[2];
sx q[2];
rz(-0.80500717) q[2];
sx q[2];
rz(0.43078255) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.731175) q[1];
sx q[1];
rz(-1.2520619) q[1];
sx q[1];
rz(-2.9916863) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0722854) q[3];
sx q[3];
rz(-1.5007334) q[3];
sx q[3];
rz(-2.4481776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.80028278) q[2];
sx q[2];
rz(-1.492123) q[2];
sx q[2];
rz(-0.71551234) q[2];
rz(-1.4095151) q[3];
sx q[3];
rz(-0.87604299) q[3];
sx q[3];
rz(-1.5040262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8228804) q[0];
sx q[0];
rz(-0.57682288) q[0];
sx q[0];
rz(0.034828287) q[0];
rz(2.4233129) q[1];
sx q[1];
rz(-1.7363345) q[1];
sx q[1];
rz(1.1911596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.98004) q[0];
sx q[0];
rz(-1.411952) q[0];
sx q[0];
rz(-0.94670403) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4592752) q[2];
sx q[2];
rz(-1.6860262) q[2];
sx q[2];
rz(0.29733411) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9893966) q[1];
sx q[1];
rz(-0.65534822) q[1];
sx q[1];
rz(-2.6776777) q[1];
rz(-pi) q[2];
rz(2.9510849) q[3];
sx q[3];
rz(-1.1393875) q[3];
sx q[3];
rz(2.0634758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3388153) q[2];
sx q[2];
rz(-1.2884527) q[2];
sx q[2];
rz(3.0954933) q[2];
rz(-2.3378546) q[3];
sx q[3];
rz(-1.9019144) q[3];
sx q[3];
rz(-0.55155915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7130724) q[0];
sx q[0];
rz(-1.5621194) q[0];
sx q[0];
rz(-2.5285316) q[0];
rz(0.26519457) q[1];
sx q[1];
rz(-1.3399597) q[1];
sx q[1];
rz(0.048096098) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0408024) q[0];
sx q[0];
rz(-1.5682966) q[0];
sx q[0];
rz(-1.0889755) q[0];
x q[1];
rz(0.71718054) q[2];
sx q[2];
rz(-2.4240026) q[2];
sx q[2];
rz(-2.8924475) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9768757) q[1];
sx q[1];
rz(-1.6304468) q[1];
sx q[1];
rz(-0.05257923) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7907324) q[3];
sx q[3];
rz(-1.8985629) q[3];
sx q[3];
rz(-1.8192489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.6951822) q[2];
sx q[2];
rz(-1.9789275) q[2];
sx q[2];
rz(-0.41333684) q[2];
rz(0.6684331) q[3];
sx q[3];
rz(-1.8139402) q[3];
sx q[3];
rz(-1.7641164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
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
rz(0.17871019) q[0];
sx q[0];
rz(-2.9229735) q[0];
sx q[0];
rz(-2.9982153) q[0];
rz(-0.42477056) q[1];
sx q[1];
rz(-1.9353341) q[1];
sx q[1];
rz(-0.43407789) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33725564) q[0];
sx q[0];
rz(-1.6640264) q[0];
sx q[0];
rz(0.038196724) q[0];
rz(0.6521449) q[2];
sx q[2];
rz(-1.6639478) q[2];
sx q[2];
rz(0.31983122) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.078316704) q[1];
sx q[1];
rz(-2.6022291) q[1];
sx q[1];
rz(-1.2961948) q[1];
rz(1.3531209) q[3];
sx q[3];
rz(-0.77131144) q[3];
sx q[3];
rz(1.2319437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8726337) q[2];
sx q[2];
rz(-2.3838398) q[2];
sx q[2];
rz(0.41716519) q[2];
rz(-0.68552351) q[3];
sx q[3];
rz(-1.439582) q[3];
sx q[3];
rz(3.0912257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79614574) q[0];
sx q[0];
rz(-0.31679994) q[0];
sx q[0];
rz(-1.5078804) q[0];
rz(-0.66868526) q[1];
sx q[1];
rz(-1.94328) q[1];
sx q[1];
rz(1.1315469) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9848081) q[0];
sx q[0];
rz(-2.1685792) q[0];
sx q[0];
rz(-3.1020107) q[0];
rz(-pi) q[1];
rz(1.0102398) q[2];
sx q[2];
rz(-1.8565289) q[2];
sx q[2];
rz(0.066372685) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3858429) q[1];
sx q[1];
rz(-1.3846701) q[1];
sx q[1];
rz(-2.666074) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.442836) q[3];
sx q[3];
rz(-1.4528414) q[3];
sx q[3];
rz(-1.02533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23182997) q[2];
sx q[2];
rz(-2.6015687) q[2];
sx q[2];
rz(0.46169272) q[2];
rz(-0.26990226) q[3];
sx q[3];
rz(-1.88068) q[3];
sx q[3];
rz(-0.65129483) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57399026) q[0];
sx q[0];
rz(-0.42581588) q[0];
sx q[0];
rz(2.0800166) q[0];
rz(-0.65879446) q[1];
sx q[1];
rz(-1.7196722) q[1];
sx q[1];
rz(-2.7366791) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.957307) q[0];
sx q[0];
rz(-2.4290963) q[0];
sx q[0];
rz(0.16100968) q[0];
rz(-pi) q[1];
rz(-0.54203029) q[2];
sx q[2];
rz(-1.932229) q[2];
sx q[2];
rz(0.76914495) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4872849) q[1];
sx q[1];
rz(-0.52738076) q[1];
sx q[1];
rz(1.2913778) q[1];
rz(-1.177321) q[3];
sx q[3];
rz(-2.4654536) q[3];
sx q[3];
rz(-0.87489389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0460661) q[2];
sx q[2];
rz(-2.1370856) q[2];
sx q[2];
rz(1.3859762) q[2];
rz(-0.69685495) q[3];
sx q[3];
rz(-0.98074061) q[3];
sx q[3];
rz(-0.81010404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54414576) q[0];
sx q[0];
rz(-0.61328855) q[0];
sx q[0];
rz(2.5469653) q[0];
rz(1.1320629) q[1];
sx q[1];
rz(-1.4040399) q[1];
sx q[1];
rz(-2.6304257) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5385732) q[0];
sx q[0];
rz(-1.5165189) q[0];
sx q[0];
rz(0.57247573) q[0];
rz(-pi) q[1];
rz(0.86894247) q[2];
sx q[2];
rz(-2.7232183) q[2];
sx q[2];
rz(1.4792031) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0446544) q[1];
sx q[1];
rz(-1.6676098) q[1];
sx q[1];
rz(1.3083878) q[1];
rz(1.0812382) q[3];
sx q[3];
rz(-2.2415906) q[3];
sx q[3];
rz(1.2257934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1749997) q[2];
sx q[2];
rz(-2.9006557) q[2];
sx q[2];
rz(2.8540376) q[2];
rz(-1.1047085) q[3];
sx q[3];
rz(-1.46773) q[3];
sx q[3];
rz(-3.0159359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6519046) q[0];
sx q[0];
rz(-0.37864417) q[0];
sx q[0];
rz(-2.4878159) q[0];
rz(-2.6405624) q[1];
sx q[1];
rz(-2.4153695) q[1];
sx q[1];
rz(-0.78530606) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0806431) q[0];
sx q[0];
rz(-1.4333581) q[0];
sx q[0];
rz(-2.5713443) q[0];
x q[1];
rz(-2.5748524) q[2];
sx q[2];
rz(-1.9575685) q[2];
sx q[2];
rz(0.45505986) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.056902) q[1];
sx q[1];
rz(-2.4246019) q[1];
sx q[1];
rz(-2.8193974) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4887962) q[3];
sx q[3];
rz(-1.7991711) q[3];
sx q[3];
rz(2.5395951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1886957) q[2];
sx q[2];
rz(-1.5263824) q[2];
sx q[2];
rz(2.7719899) q[2];
rz(0.26436198) q[3];
sx q[3];
rz(-1.2814458) q[3];
sx q[3];
rz(3.1407691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13067506) q[0];
sx q[0];
rz(-2.4322746) q[0];
sx q[0];
rz(-1.5189019) q[0];
rz(1.2394637) q[1];
sx q[1];
rz(-2.0294956) q[1];
sx q[1];
rz(0.94737238) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1993123) q[0];
sx q[0];
rz(-2.7373642) q[0];
sx q[0];
rz(1.72615) q[0];
x q[1];
rz(2.7804393) q[2];
sx q[2];
rz(-1.5120107) q[2];
sx q[2];
rz(1.0967983) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.67657802) q[1];
sx q[1];
rz(-1.4099219) q[1];
sx q[1];
rz(-1.4032928) q[1];
x q[2];
rz(-1.2227433) q[3];
sx q[3];
rz(-2.9677941) q[3];
sx q[3];
rz(2.9617975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47306791) q[2];
sx q[2];
rz(-0.99978414) q[2];
sx q[2];
rz(2.0295985) q[2];
rz(0.2019349) q[3];
sx q[3];
rz(-0.58378059) q[3];
sx q[3];
rz(2.7659265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428381) q[0];
sx q[0];
rz(-0.18809479) q[0];
sx q[0];
rz(1.4659708) q[0];
rz(1.6000043) q[1];
sx q[1];
rz(-1.3009289) q[1];
sx q[1];
rz(-2.8864554) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8153471) q[0];
sx q[0];
rz(-1.7044984) q[0];
sx q[0];
rz(3.1241425) q[0];
x q[1];
rz(-2.3681247) q[2];
sx q[2];
rz(-1.8102526) q[2];
sx q[2];
rz(2.4635893) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3465187) q[1];
sx q[1];
rz(-2.8523905) q[1];
sx q[1];
rz(1.2917915) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3794243) q[3];
sx q[3];
rz(-0.81052033) q[3];
sx q[3];
rz(2.9178491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5280891) q[2];
sx q[2];
rz(-1.1897503) q[2];
sx q[2];
rz(-0.36821723) q[2];
rz(-0.28299371) q[3];
sx q[3];
rz(-2.8734983) q[3];
sx q[3];
rz(0.32576573) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1267283) q[0];
sx q[0];
rz(-0.8664425) q[0];
sx q[0];
rz(2.0436825) q[0];
rz(2.6425843) q[1];
sx q[1];
rz(-1.9022763) q[1];
sx q[1];
rz(0.3872445) q[1];
rz(1.6989742) q[2];
sx q[2];
rz(-1.9009931) q[2];
sx q[2];
rz(-0.92585678) q[2];
rz(1.2756497) q[3];
sx q[3];
rz(-2.9316918) q[3];
sx q[3];
rz(1.5422921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

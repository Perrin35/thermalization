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
rz(5.0205686) q[0];
sx q[0];
rz(10.159259) q[0];
rz(1.6425411) q[1];
sx q[1];
rz(1.8416815) q[1];
sx q[1];
rz(11.587974) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1837195) q[0];
sx q[0];
rz(-1.5629725) q[0];
sx q[0];
rz(1.70723) q[0];
rz(-2.5402563) q[2];
sx q[2];
rz(-2.3365855) q[2];
sx q[2];
rz(2.7108101) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4104176) q[1];
sx q[1];
rz(-1.2520619) q[1];
sx q[1];
rz(0.14990633) q[1];
x q[2];
rz(-0.079732883) q[3];
sx q[3];
rz(-1.0736205) q[3];
sx q[3];
rz(2.2261208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3413099) q[2];
sx q[2];
rz(-1.6494696) q[2];
sx q[2];
rz(-0.71551234) q[2];
rz(1.7320775) q[3];
sx q[3];
rz(-0.87604299) q[3];
sx q[3];
rz(1.6375665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3187123) q[0];
sx q[0];
rz(-2.5647698) q[0];
sx q[0];
rz(3.1067644) q[0];
rz(-0.71827978) q[1];
sx q[1];
rz(-1.4052582) q[1];
sx q[1];
rz(-1.1911596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.98004) q[0];
sx q[0];
rz(-1.7296407) q[0];
sx q[0];
rz(0.94670403) q[0];
rz(-pi) q[1];
rz(-1.7188363) q[2];
sx q[2];
rz(-0.89385539) q[2];
sx q[2];
rz(-1.9612775) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.71514788) q[1];
sx q[1];
rz(-0.99438018) q[1];
sx q[1];
rz(-1.2395241) q[1];
rz(-2.9510849) q[3];
sx q[3];
rz(-1.1393875) q[3];
sx q[3];
rz(-2.0634758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8027773) q[2];
sx q[2];
rz(-1.85314) q[2];
sx q[2];
rz(-3.0954933) q[2];
rz(-0.80373803) q[3];
sx q[3];
rz(-1.2396783) q[3];
sx q[3];
rz(-0.55155915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7130724) q[0];
sx q[0];
rz(-1.5621194) q[0];
sx q[0];
rz(0.61306104) q[0];
rz(-2.8763981) q[1];
sx q[1];
rz(-1.801633) q[1];
sx q[1];
rz(-0.048096098) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4652258) q[0];
sx q[0];
rz(-0.48182677) q[0];
sx q[0];
rz(-1.5761907) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4244121) q[2];
sx q[2];
rz(-0.71759008) q[2];
sx q[2];
rz(-2.8924475) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4029419) q[1];
sx q[1];
rz(-1.5183107) q[1];
sx q[1];
rz(-1.5110635) q[1];
x q[2];
rz(0.35086029) q[3];
sx q[3];
rz(-1.8985629) q[3];
sx q[3];
rz(-1.8192489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4464104) q[2];
sx q[2];
rz(-1.1626652) q[2];
sx q[2];
rz(0.41333684) q[2];
rz(-2.4731596) q[3];
sx q[3];
rz(-1.3276525) q[3];
sx q[3];
rz(-1.3774762) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9628825) q[0];
sx q[0];
rz(-2.9229735) q[0];
sx q[0];
rz(0.14337732) q[0];
rz(0.42477056) q[1];
sx q[1];
rz(-1.9353341) q[1];
sx q[1];
rz(-2.7075148) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052267678) q[0];
sx q[0];
rz(-0.10073034) q[0];
sx q[0];
rz(-1.9585376) q[0];
x q[1];
rz(-1.4537895) q[2];
sx q[2];
rz(-0.9219578) q[2];
sx q[2];
rz(-1.1800571) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0632759) q[1];
sx q[1];
rz(-2.6022291) q[1];
sx q[1];
rz(-1.2961948) q[1];
x q[2];
rz(1.7884718) q[3];
sx q[3];
rz(-0.77131144) q[3];
sx q[3];
rz(-1.2319437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.268959) q[2];
sx q[2];
rz(-2.3838398) q[2];
sx q[2];
rz(0.41716519) q[2];
rz(0.68552351) q[3];
sx q[3];
rz(-1.7020107) q[3];
sx q[3];
rz(-0.050366966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.79614574) q[0];
sx q[0];
rz(-0.31679994) q[0];
sx q[0];
rz(1.5078804) q[0];
rz(-2.4729074) q[1];
sx q[1];
rz(-1.94328) q[1];
sx q[1];
rz(-1.1315469) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.749866) q[0];
sx q[0];
rz(-1.6035115) q[0];
sx q[0];
rz(-2.1689438) q[0];
rz(-pi) q[1];
rz(-2.8077233) q[2];
sx q[2];
rz(-2.1060995) q[2];
sx q[2];
rz(1.4620505) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.9100489) q[1];
sx q[1];
rz(-1.1041512) q[1];
sx q[1];
rz(-1.3620767) q[1];
rz(-0.11891811) q[3];
sx q[3];
rz(-1.6978627) q[3];
sx q[3];
rz(-0.53032622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9097627) q[2];
sx q[2];
rz(-0.54002395) q[2];
sx q[2];
rz(-0.46169272) q[2];
rz(0.26990226) q[3];
sx q[3];
rz(-1.88068) q[3];
sx q[3];
rz(0.65129483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57399026) q[0];
sx q[0];
rz(-0.42581588) q[0];
sx q[0];
rz(1.061576) q[0];
rz(-0.65879446) q[1];
sx q[1];
rz(-1.7196722) q[1];
sx q[1];
rz(-2.7366791) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8773738) q[0];
sx q[0];
rz(-1.675791) q[0];
sx q[0];
rz(-0.70606249) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63240856) q[2];
sx q[2];
rz(-0.64130613) q[2];
sx q[2];
rz(0.270688) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8078708) q[1];
sx q[1];
rz(-2.0757339) q[1];
sx q[1];
rz(-0.15926475) q[1];
rz(-pi) q[2];
x q[2];
rz(1.177321) q[3];
sx q[3];
rz(-2.4654536) q[3];
sx q[3];
rz(-2.2666988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.095526516) q[2];
sx q[2];
rz(-2.1370856) q[2];
sx q[2];
rz(1.7556165) q[2];
rz(0.69685495) q[3];
sx q[3];
rz(-2.160852) q[3];
sx q[3];
rz(-0.81010404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.54414576) q[0];
sx q[0];
rz(-0.61328855) q[0];
sx q[0];
rz(2.5469653) q[0];
rz(-2.0095297) q[1];
sx q[1];
rz(-1.7375528) q[1];
sx q[1];
rz(-0.51116699) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1388672) q[0];
sx q[0];
rz(-2.1423233) q[0];
sx q[0];
rz(-1.6353428) q[0];
rz(0.2795477) q[2];
sx q[2];
rz(-1.8862533) q[2];
sx q[2];
rz(-2.4089898) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6936924) q[1];
sx q[1];
rz(-1.8319472) q[1];
sx q[1];
rz(3.0413701) q[1];
rz(-pi) q[2];
rz(2.0603544) q[3];
sx q[3];
rz(-2.2415906) q[3];
sx q[3];
rz(1.9157992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1749997) q[2];
sx q[2];
rz(-2.9006557) q[2];
sx q[2];
rz(-2.8540376) q[2];
rz(-2.0368841) q[3];
sx q[3];
rz(-1.6738626) q[3];
sx q[3];
rz(-3.0159359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6519046) q[0];
sx q[0];
rz(-0.37864417) q[0];
sx q[0];
rz(2.4878159) q[0];
rz(-2.6405624) q[1];
sx q[1];
rz(-0.72622314) q[1];
sx q[1];
rz(-2.3562866) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5637999) q[0];
sx q[0];
rz(-1.0065838) q[0];
sx q[0];
rz(1.4079421) q[0];
x q[1];
rz(-2.5748524) q[2];
sx q[2];
rz(-1.9575685) q[2];
sx q[2];
rz(-2.6865328) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.056902) q[1];
sx q[1];
rz(-0.71699079) q[1];
sx q[1];
rz(0.3221953) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2291181) q[3];
sx q[3];
rz(-1.6506628) q[3];
sx q[3];
rz(0.95019482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1886957) q[2];
sx q[2];
rz(-1.6152103) q[2];
sx q[2];
rz(0.3696028) q[2];
rz(0.26436198) q[3];
sx q[3];
rz(-1.8601469) q[3];
sx q[3];
rz(-3.1407691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13067506) q[0];
sx q[0];
rz(-2.4322746) q[0];
sx q[0];
rz(-1.5189019) q[0];
rz(1.9021289) q[1];
sx q[1];
rz(-2.0294956) q[1];
sx q[1];
rz(2.1942203) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1993123) q[0];
sx q[0];
rz(-0.40422842) q[0];
sx q[0];
rz(-1.4154427) q[0];
rz(-pi) q[1];
rz(-0.36115335) q[2];
sx q[2];
rz(-1.5120107) q[2];
sx q[2];
rz(1.0967983) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4650146) q[1];
sx q[1];
rz(-1.4099219) q[1];
sx q[1];
rz(-1.4032928) q[1];
rz(-1.4072284) q[3];
sx q[3];
rz(-1.6298098) q[3];
sx q[3];
rz(-1.407377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.47306791) q[2];
sx q[2];
rz(-0.99978414) q[2];
sx q[2];
rz(-2.0295985) q[2];
rz(-2.9396577) q[3];
sx q[3];
rz(-0.58378059) q[3];
sx q[3];
rz(-0.37566617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(2.9428381) q[0];
sx q[0];
rz(-2.9534979) q[0];
sx q[0];
rz(1.4659708) q[0];
rz(-1.5415883) q[1];
sx q[1];
rz(-1.3009289) q[1];
sx q[1];
rz(0.25513729) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45642325) q[0];
sx q[0];
rz(-0.13482929) q[0];
sx q[0];
rz(1.441787) q[0];
rz(-pi) q[1];
rz(-1.2419644) q[2];
sx q[2];
rz(-0.82482289) q[2];
sx q[2];
rz(1.1203631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0855116) q[1];
sx q[1];
rz(-1.8485118) q[1];
sx q[1];
rz(3.059832) q[1];
x q[2];
rz(2.3721123) q[3];
sx q[3];
rz(-1.4325241) q[3];
sx q[3];
rz(-1.4797803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5280891) q[2];
sx q[2];
rz(-1.9518423) q[2];
sx q[2];
rz(2.7733754) q[2];
rz(0.28299371) q[3];
sx q[3];
rz(-0.26809433) q[3];
sx q[3];
rz(0.32576573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1267283) q[0];
sx q[0];
rz(-0.8664425) q[0];
sx q[0];
rz(2.0436825) q[0];
rz(0.49900838) q[1];
sx q[1];
rz(-1.2393163) q[1];
sx q[1];
rz(-2.7543482) q[1];
rz(-0.3327315) q[2];
sx q[2];
rz(-1.6920148) q[2];
sx q[2];
rz(0.68670338) q[2];
rz(1.8659429) q[3];
sx q[3];
rz(-0.20990089) q[3];
sx q[3];
rz(-1.5993006) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

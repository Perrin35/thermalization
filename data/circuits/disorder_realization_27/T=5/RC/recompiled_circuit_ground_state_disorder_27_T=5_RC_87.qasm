OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34484997) q[0];
sx q[0];
rz(-0.27422187) q[0];
sx q[0];
rz(0.56871498) q[0];
rz(-5.0721726) q[1];
sx q[1];
rz(0.99744263) q[1];
sx q[1];
rz(9.6923516) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28963213) q[0];
sx q[0];
rz(-1.7718162) q[0];
sx q[0];
rz(-2.5876849) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11315364) q[2];
sx q[2];
rz(-1.9027766) q[2];
sx q[2];
rz(-1.6210213) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.888511) q[1];
sx q[1];
rz(-0.0086697658) q[1];
sx q[1];
rz(-1.8523097) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1564134) q[3];
sx q[3];
rz(-2.9880045) q[3];
sx q[3];
rz(-1.6382662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8895175) q[2];
sx q[2];
rz(-1.0675665) q[2];
sx q[2];
rz(-2.0102823) q[2];
rz(0.45025292) q[3];
sx q[3];
rz(-2.4501652) q[3];
sx q[3];
rz(0.94436193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64304072) q[0];
sx q[0];
rz(-2.2096071) q[0];
sx q[0];
rz(-1.407628) q[0];
rz(-2.6990926) q[1];
sx q[1];
rz(-1.718037) q[1];
sx q[1];
rz(-0.59534591) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.731643) q[0];
sx q[0];
rz(-2.7945257) q[0];
sx q[0];
rz(-0.76526977) q[0];
rz(-pi) q[1];
rz(2.1791502) q[2];
sx q[2];
rz(-2.3909937) q[2];
sx q[2];
rz(-0.19718328) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7253984) q[1];
sx q[1];
rz(-0.42896118) q[1];
sx q[1];
rz(-0.012417656) q[1];
rz(-0.54259681) q[3];
sx q[3];
rz(-2.4568181) q[3];
sx q[3];
rz(-2.7921576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2596316) q[2];
sx q[2];
rz(-2.5233848) q[2];
sx q[2];
rz(-2.372443) q[2];
rz(-1.7800356) q[3];
sx q[3];
rz(-0.96746126) q[3];
sx q[3];
rz(0.59534016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.896987) q[0];
sx q[0];
rz(-2.2266882) q[0];
sx q[0];
rz(-0.33682522) q[0];
rz(-1.4312076) q[1];
sx q[1];
rz(-2.2957048) q[1];
sx q[1];
rz(-0.097188458) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4405085) q[0];
sx q[0];
rz(-1.3680662) q[0];
sx q[0];
rz(1.0725783) q[0];
rz(0.23938208) q[2];
sx q[2];
rz(-2.7500543) q[2];
sx q[2];
rz(1.7042314) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5775902) q[1];
sx q[1];
rz(-1.9158984) q[1];
sx q[1];
rz(1.7175098) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95366565) q[3];
sx q[3];
rz(-0.78804555) q[3];
sx q[3];
rz(-0.55824454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1691957) q[2];
sx q[2];
rz(-1.376386) q[2];
sx q[2];
rz(0.81673679) q[2];
rz(2.2190602) q[3];
sx q[3];
rz(-0.43729344) q[3];
sx q[3];
rz(-1.1923265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794491) q[0];
sx q[0];
rz(-1.3379931) q[0];
sx q[0];
rz(2.2747967) q[0];
rz(1.9056994) q[1];
sx q[1];
rz(-1.1150603) q[1];
sx q[1];
rz(-2.8996276) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2529255) q[0];
sx q[0];
rz(-2.2641716) q[0];
sx q[0];
rz(2.7916551) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.073512065) q[2];
sx q[2];
rz(-1.4097555) q[2];
sx q[2];
rz(1.1390151) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.82034012) q[1];
sx q[1];
rz(-2.7124321) q[1];
sx q[1];
rz(-0.33351516) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6446043) q[3];
sx q[3];
rz(-1.6145633) q[3];
sx q[3];
rz(0.62326335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8786826) q[2];
sx q[2];
rz(-1.6308558) q[2];
sx q[2];
rz(-0.53544694) q[2];
rz(0.49324909) q[3];
sx q[3];
rz(-2.2507164) q[3];
sx q[3];
rz(0.44805995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.066147476) q[0];
sx q[0];
rz(-0.15114052) q[0];
sx q[0];
rz(-2.2976663) q[0];
rz(-1.6663724) q[1];
sx q[1];
rz(-1.6553144) q[1];
sx q[1];
rz(-0.39316887) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61680865) q[0];
sx q[0];
rz(-2.3032585) q[0];
sx q[0];
rz(-2.7969267) q[0];
rz(-0.63661036) q[2];
sx q[2];
rz(-1.3946956) q[2];
sx q[2];
rz(-2.2409093) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1816493) q[1];
sx q[1];
rz(-1.8487329) q[1];
sx q[1];
rz(-0.47680579) q[1];
rz(-2.353998) q[3];
sx q[3];
rz(-1.2114085) q[3];
sx q[3];
rz(1.5415292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7273442) q[2];
sx q[2];
rz(-2.9372637) q[2];
sx q[2];
rz(-0.3981398) q[2];
rz(2.5967755) q[3];
sx q[3];
rz(-0.78854338) q[3];
sx q[3];
rz(-1.3032234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2274461) q[0];
sx q[0];
rz(-1.0108203) q[0];
sx q[0];
rz(2.2858802) q[0];
rz(0.69264597) q[1];
sx q[1];
rz(-0.99450642) q[1];
sx q[1];
rz(-2.8725502) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1097668) q[0];
sx q[0];
rz(-2.9485011) q[0];
sx q[0];
rz(1.6982128) q[0];
x q[1];
rz(-1.6330209) q[2];
sx q[2];
rz(-2.6341558) q[2];
sx q[2];
rz(-0.98558805) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0647886) q[1];
sx q[1];
rz(-2.2858983) q[1];
sx q[1];
rz(0.24311693) q[1];
rz(-2.2773197) q[3];
sx q[3];
rz(-2.8916997) q[3];
sx q[3];
rz(-2.8003729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0272224) q[2];
sx q[2];
rz(-1.3193069) q[2];
sx q[2];
rz(3.031292) q[2];
rz(0.86841622) q[3];
sx q[3];
rz(-1.1766368) q[3];
sx q[3];
rz(1.4364012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39349839) q[0];
sx q[0];
rz(-0.49999923) q[0];
sx q[0];
rz(-0.86135832) q[0];
rz(1.512108) q[1];
sx q[1];
rz(-1.6810828) q[1];
sx q[1];
rz(0.82383627) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31164661) q[0];
sx q[0];
rz(-0.59893805) q[0];
sx q[0];
rz(-1.7700559) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4365224) q[2];
sx q[2];
rz(-2.3425205) q[2];
sx q[2];
rz(2.7288849) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41125248) q[1];
sx q[1];
rz(-1.7837875) q[1];
sx q[1];
rz(0.61720444) q[1];
rz(-0.40315513) q[3];
sx q[3];
rz(-0.69270999) q[3];
sx q[3];
rz(1.5628536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5065826) q[2];
sx q[2];
rz(-0.51598769) q[2];
sx q[2];
rz(0.79279509) q[2];
rz(-0.005216287) q[3];
sx q[3];
rz(-0.79013932) q[3];
sx q[3];
rz(1.4504356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1714627) q[0];
sx q[0];
rz(-1.9487533) q[0];
sx q[0];
rz(-0.18950732) q[0];
rz(2.3686523) q[1];
sx q[1];
rz(-2.645292) q[1];
sx q[1];
rz(-2.5453087) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4204105) q[0];
sx q[0];
rz(-2.2503958) q[0];
sx q[0];
rz(-1.3534989) q[0];
rz(-2.3438966) q[2];
sx q[2];
rz(-2.0828649) q[2];
sx q[2];
rz(1.860581) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.50508037) q[1];
sx q[1];
rz(-0.13992913) q[1];
sx q[1];
rz(-2.416978) q[1];
rz(-pi) q[2];
rz(2.1607481) q[3];
sx q[3];
rz(-2.337237) q[3];
sx q[3];
rz(2.435702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.72558713) q[2];
sx q[2];
rz(-3.1276939) q[2];
sx q[2];
rz(-1.2131946) q[2];
rz(-0.98617918) q[3];
sx q[3];
rz(-1.4158019) q[3];
sx q[3];
rz(1.2590316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0202494) q[0];
sx q[0];
rz(-1.5556524) q[0];
sx q[0];
rz(1.1100618) q[0];
rz(2.8129261) q[1];
sx q[1];
rz(-1.5866491) q[1];
sx q[1];
rz(1.2967671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27889869) q[0];
sx q[0];
rz(-1.6644913) q[0];
sx q[0];
rz(-2.2298032) q[0];
x q[1];
rz(-1.718156) q[2];
sx q[2];
rz(-2.7412716) q[2];
sx q[2];
rz(1.6860698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3039141) q[1];
sx q[1];
rz(-2.0635567) q[1];
sx q[1];
rz(-1.4999092) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4761183) q[3];
sx q[3];
rz(-1.7413119) q[3];
sx q[3];
rz(-2.1297034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0729596) q[2];
sx q[2];
rz(-2.1559842) q[2];
sx q[2];
rz(-0.10406058) q[2];
rz(1.1348628) q[3];
sx q[3];
rz(-1.3645423) q[3];
sx q[3];
rz(-0.22741905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.65570152) q[0];
sx q[0];
rz(-0.98485297) q[0];
sx q[0];
rz(2.4110598) q[0];
rz(0.55745521) q[1];
sx q[1];
rz(-1.1703706) q[1];
sx q[1];
rz(0.4298068) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8107501) q[0];
sx q[0];
rz(-1.5702899) q[0];
sx q[0];
rz(2.3798126) q[0];
x q[1];
rz(2.7676959) q[2];
sx q[2];
rz(-1.8086401) q[2];
sx q[2];
rz(-1.6872981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.046172) q[1];
sx q[1];
rz(-0.37994994) q[1];
sx q[1];
rz(-1.7871961) q[1];
rz(-0.12390512) q[3];
sx q[3];
rz(-1.9384111) q[3];
sx q[3];
rz(-2.1218561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.80959117) q[2];
sx q[2];
rz(-1.0217228) q[2];
sx q[2];
rz(2.8487955) q[2];
rz(3.0012567) q[3];
sx q[3];
rz(-1.0293181) q[3];
sx q[3];
rz(0.50104195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(0.3705227) q[0];
sx q[0];
rz(-1.6765544) q[0];
sx q[0];
rz(0.21677207) q[0];
rz(-2.4222005) q[1];
sx q[1];
rz(-1.2750625) q[1];
sx q[1];
rz(0.19663179) q[1];
rz(1.8875296) q[2];
sx q[2];
rz(-2.6120196) q[2];
sx q[2];
rz(2.0092464) q[2];
rz(1.8932395) q[3];
sx q[3];
rz(-1.697968) q[3];
sx q[3];
rz(-1.9597114) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

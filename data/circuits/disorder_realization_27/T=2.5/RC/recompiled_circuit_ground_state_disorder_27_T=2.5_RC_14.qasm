OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5840924) q[0];
sx q[0];
rz(3.1182365) q[0];
sx q[0];
rz(11.630848) q[0];
rz(-1.6159396) q[1];
sx q[1];
rz(-1.5204117) q[1];
sx q[1];
rz(2.867155) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68831681) q[0];
sx q[0];
rz(-0.098628086) q[0];
sx q[0];
rz(0.80699705) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1054613) q[2];
sx q[2];
rz(-1.0792152) q[2];
sx q[2];
rz(-0.21869379) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.5861862) q[1];
sx q[1];
rz(-1.6073391) q[1];
sx q[1];
rz(0.010374109) q[1];
x q[2];
rz(-1.7414344) q[3];
sx q[3];
rz(-1.6434165) q[3];
sx q[3];
rz(-2.448248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9258257) q[2];
sx q[2];
rz(-3.1302858) q[2];
sx q[2];
rz(2.0634148) q[2];
rz(-0.82873851) q[3];
sx q[3];
rz(-1.6308035) q[3];
sx q[3];
rz(-0.80882788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082212903) q[0];
sx q[0];
rz(-1.8635211) q[0];
sx q[0];
rz(-1.7510121) q[0];
rz(1.4311721) q[1];
sx q[1];
rz(-0.0043967604) q[1];
sx q[1];
rz(-1.7079401) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2947049) q[0];
sx q[0];
rz(-1.4661015) q[0];
sx q[0];
rz(-2.2082735) q[0];
x q[1];
rz(1.1270056) q[2];
sx q[2];
rz(-1.5856885) q[2];
sx q[2];
rz(-3.1010951) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.28273496) q[1];
sx q[1];
rz(-1.5703859) q[1];
sx q[1];
rz(1.5617227) q[1];
rz(-pi) q[2];
rz(0.065304718) q[3];
sx q[3];
rz(-2.3529144) q[3];
sx q[3];
rz(-2.1826133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7626875) q[2];
sx q[2];
rz(-1.5451558) q[2];
sx q[2];
rz(-1.572466) q[2];
rz(2.3804741) q[3];
sx q[3];
rz(-0.053839024) q[3];
sx q[3];
rz(-1.9332473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91956562) q[0];
sx q[0];
rz(-2.5313105) q[0];
sx q[0];
rz(-2.5797504) q[0];
rz(-1.5737083) q[1];
sx q[1];
rz(-1.578873) q[1];
sx q[1];
rz(3.124253) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0578831) q[0];
sx q[0];
rz(-2.0590933) q[0];
sx q[0];
rz(2.6668307) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1465535) q[2];
sx q[2];
rz(-2.3684566) q[2];
sx q[2];
rz(0.38627689) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8120809) q[1];
sx q[1];
rz(-1.2083665) q[1];
sx q[1];
rz(0.017048841) q[1];
rz(-1.5478327) q[3];
sx q[3];
rz(-2.0809789) q[3];
sx q[3];
rz(1.2018454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5138381) q[2];
sx q[2];
rz(-1.4189812) q[2];
sx q[2];
rz(2.5886152) q[2];
rz(-1.2051469) q[3];
sx q[3];
rz(-1.5670992) q[3];
sx q[3];
rz(1.5945826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39128458) q[0];
sx q[0];
rz(-1.0306083) q[0];
sx q[0];
rz(1.3543825) q[0];
rz(-1.7693819) q[1];
sx q[1];
rz(-0.0028227614) q[1];
sx q[1];
rz(1.3815968) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1005122) q[0];
sx q[0];
rz(-1.1100475) q[0];
sx q[0];
rz(2.6897088) q[0];
rz(-pi) q[1];
rz(-0.69475485) q[2];
sx q[2];
rz(-3.1379897) q[2];
sx q[2];
rz(2.3191593) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.83643944) q[1];
sx q[1];
rz(-0.90739319) q[1];
sx q[1];
rz(-2.7022578) q[1];
rz(-pi) q[2];
rz(-2.3390807) q[3];
sx q[3];
rz(-1.7798692) q[3];
sx q[3];
rz(2.0627562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0733205) q[2];
sx q[2];
rz(-3.1219411) q[2];
sx q[2];
rz(1.9015296) q[2];
rz(-0.23010075) q[3];
sx q[3];
rz(-0.0041882526) q[3];
sx q[3];
rz(-0.41964644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0536026) q[0];
sx q[0];
rz(-2.3542861) q[0];
sx q[0];
rz(-1.2552274) q[0];
rz(3.1322196) q[1];
sx q[1];
rz(-1.7724937) q[1];
sx q[1];
rz(0.031551687) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2173715) q[0];
sx q[0];
rz(-1.517202) q[0];
sx q[0];
rz(1.0384667) q[0];
x q[1];
rz(1.2826233) q[2];
sx q[2];
rz(-2.8409578) q[2];
sx q[2];
rz(-2.2733781) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0134558) q[1];
sx q[1];
rz(-2.0990879) q[1];
sx q[1];
rz(-0.25685132) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58720354) q[3];
sx q[3];
rz(-1.0343583) q[3];
sx q[3];
rz(-1.0565384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8072529) q[2];
sx q[2];
rz(-3.1355317) q[2];
sx q[2];
rz(0.80017153) q[2];
rz(-0.79450327) q[3];
sx q[3];
rz(-0.032363351) q[3];
sx q[3];
rz(-0.86875027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.054258) q[0];
sx q[0];
rz(-2.9825409) q[0];
sx q[0];
rz(1.6416838) q[0];
rz(0.17290393) q[1];
sx q[1];
rz(-0.042363107) q[1];
sx q[1];
rz(-3.0917047) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.848441) q[0];
sx q[0];
rz(-1.3209136) q[0];
sx q[0];
rz(-0.94012733) q[0];
rz(-pi) q[1];
rz(0.11325963) q[2];
sx q[2];
rz(-2.335603) q[2];
sx q[2];
rz(1.2330556) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3059002) q[1];
sx q[1];
rz(-2.3819807) q[1];
sx q[1];
rz(1.9363957) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6088151) q[3];
sx q[3];
rz(-2.2190337) q[3];
sx q[3];
rz(1.6870354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.738203) q[2];
sx q[2];
rz(-3.093779) q[2];
sx q[2];
rz(1.8955463) q[2];
rz(-1.3561148) q[3];
sx q[3];
rz(-0.035364371) q[3];
sx q[3];
rz(-1.7252007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4288915) q[0];
sx q[0];
rz(-2.3172947) q[0];
sx q[0];
rz(-1.4225381) q[0];
rz(1.2369583) q[1];
sx q[1];
rz(-0.037038602) q[1];
sx q[1];
rz(0.2027771) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13064676) q[0];
sx q[0];
rz(-1.8955064) q[0];
sx q[0];
rz(0.79239158) q[0];
rz(-pi) q[1];
rz(-0.76272623) q[2];
sx q[2];
rz(-2.2710851) q[2];
sx q[2];
rz(-1.3415847) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7084658) q[1];
sx q[1];
rz(-2.7964292) q[1];
sx q[1];
rz(2.9776447) q[1];
rz(-1.0596541) q[3];
sx q[3];
rz(-1.7082001) q[3];
sx q[3];
rz(2.5554772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.612959) q[2];
sx q[2];
rz(-3.0411159) q[2];
sx q[2];
rz(0.67451492) q[2];
rz(-1.3450735) q[3];
sx q[3];
rz(-0.14480545) q[3];
sx q[3];
rz(1.323918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8737008) q[0];
sx q[0];
rz(-2.38509) q[0];
sx q[0];
rz(0.85195136) q[0];
rz(-2.9497228) q[1];
sx q[1];
rz(-3.1286897) q[1];
sx q[1];
rz(0.26564863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.083355) q[0];
sx q[0];
rz(-1.7449656) q[0];
sx q[0];
rz(0.41634286) q[0];
rz(-0.4531817) q[2];
sx q[2];
rz(-0.66614775) q[2];
sx q[2];
rz(1.11048) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77924624) q[1];
sx q[1];
rz(-1.6542477) q[1];
sx q[1];
rz(1.5117743) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5784078) q[3];
sx q[3];
rz(-2.5068589) q[3];
sx q[3];
rz(-2.4358482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3702281) q[2];
sx q[2];
rz(-3.0493272) q[2];
sx q[2];
rz(-2.6293758) q[2];
rz(-0.15906119) q[3];
sx q[3];
rz(-0.03511196) q[3];
sx q[3];
rz(-1.4353282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8856186) q[0];
sx q[0];
rz(-1.4775448) q[0];
sx q[0];
rz(-2.0632099) q[0];
rz(1.6400317) q[1];
sx q[1];
rz(-2.9402132) q[1];
sx q[1];
rz(-1.5837502) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5659065) q[0];
sx q[0];
rz(-0.68725902) q[0];
sx q[0];
rz(-0.024373011) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7388849) q[2];
sx q[2];
rz(-2.2103849) q[2];
sx q[2];
rz(1.0693897) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8858151) q[1];
sx q[1];
rz(-3.1193135) q[1];
sx q[1];
rz(-2.9905969) q[1];
rz(-pi) q[2];
rz(2.0191865) q[3];
sx q[3];
rz(-0.70647722) q[3];
sx q[3];
rz(-2.024533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.25643361) q[2];
sx q[2];
rz(-0.014545518) q[2];
sx q[2];
rz(-0.069615901) q[2];
rz(-3.0749248) q[3];
sx q[3];
rz(-1.0144517) q[3];
sx q[3];
rz(0.67924172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2081864) q[0];
sx q[0];
rz(-1.8649768) q[0];
sx q[0];
rz(2.8896914) q[0];
rz(1.4725641) q[1];
sx q[1];
rz(-2.9258969) q[1];
sx q[1];
rz(0.073898166) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2011482) q[0];
sx q[0];
rz(-1.5000888) q[0];
sx q[0];
rz(1.2051677) q[0];
x q[1];
rz(0.017357512) q[2];
sx q[2];
rz(-1.6021172) q[2];
sx q[2];
rz(2.2407414) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.628129) q[1];
sx q[1];
rz(-2.3766915) q[1];
sx q[1];
rz(-0.35180636) q[1];
x q[2];
rz(2.4611887) q[3];
sx q[3];
rz(-2.8814253) q[3];
sx q[3];
rz(2.5870067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4613688) q[2];
sx q[2];
rz(-0.0071439925) q[2];
sx q[2];
rz(2.3738677) q[2];
rz(-1.7420306) q[3];
sx q[3];
rz(-3.1407686) q[3];
sx q[3];
rz(0.50423938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5117699) q[0];
sx q[0];
rz(-0.98291021) q[0];
sx q[0];
rz(1.7194189) q[0];
rz(-3.1172251) q[1];
sx q[1];
rz(-0.15932803) q[1];
sx q[1];
rz(-2.9111964) q[1];
rz(1.1047614) q[2];
sx q[2];
rz(-0.73095041) q[2];
sx q[2];
rz(-1.6292844) q[2];
rz(1.4735994) q[3];
sx q[3];
rz(-1.9862277) q[3];
sx q[3];
rz(1.338892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

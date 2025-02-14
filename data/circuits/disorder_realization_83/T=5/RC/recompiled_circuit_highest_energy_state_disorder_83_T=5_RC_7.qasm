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
rz(-1.8981847) q[0];
sx q[0];
rz(-1.564448) q[0];
sx q[0];
rz(-0.91140437) q[0];
rz(2.4721594) q[1];
sx q[1];
rz(-2.8652006) q[1];
sx q[1];
rz(-0.82958329) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88032915) q[0];
sx q[0];
rz(-1.1733049) q[0];
sx q[0];
rz(2.7842159) q[0];
x q[1];
rz(1.6638512) q[2];
sx q[2];
rz(-1.3246127) q[2];
sx q[2];
rz(3.1174146) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1948283) q[1];
sx q[1];
rz(-0.84673893) q[1];
sx q[1];
rz(1.9763234) q[1];
x q[2];
rz(-0.31949701) q[3];
sx q[3];
rz(-2.6852859) q[3];
sx q[3];
rz(0.97684492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36898819) q[2];
sx q[2];
rz(-0.96869865) q[2];
sx q[2];
rz(-1.1645092) q[2];
rz(-0.72689593) q[3];
sx q[3];
rz(-1.8098857) q[3];
sx q[3];
rz(-3.0708142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1233391) q[0];
sx q[0];
rz(-2.4788719) q[0];
sx q[0];
rz(-0.29139274) q[0];
rz(-2.5413051) q[1];
sx q[1];
rz(-2.184506) q[1];
sx q[1];
rz(-2.9765863) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8931029) q[0];
sx q[0];
rz(-1.4997673) q[0];
sx q[0];
rz(1.8336589) q[0];
rz(-pi) q[1];
rz(1.1218614) q[2];
sx q[2];
rz(-1.672272) q[2];
sx q[2];
rz(-1.5047274) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9143888) q[1];
sx q[1];
rz(-1.640135) q[1];
sx q[1];
rz(0.68080618) q[1];
x q[2];
rz(0.24834541) q[3];
sx q[3];
rz(-1.6237061) q[3];
sx q[3];
rz(1.2179045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0642285) q[2];
sx q[2];
rz(-0.95244971) q[2];
sx q[2];
rz(2.9020818) q[2];
rz(0.32660487) q[3];
sx q[3];
rz(-0.17290792) q[3];
sx q[3];
rz(-0.094836205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8851682) q[0];
sx q[0];
rz(-0.435985) q[0];
sx q[0];
rz(-1.5727795) q[0];
rz(0.39372152) q[1];
sx q[1];
rz(-2.5865793) q[1];
sx q[1];
rz(0.47600019) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5167907) q[0];
sx q[0];
rz(-1.3627407) q[0];
sx q[0];
rz(2.1387908) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3838686) q[2];
sx q[2];
rz(-1.6003967) q[2];
sx q[2];
rz(1.5160949) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90895528) q[1];
sx q[1];
rz(-2.401315) q[1];
sx q[1];
rz(-0.93772447) q[1];
rz(-pi) q[2];
rz(1.908444) q[3];
sx q[3];
rz(-1.3364571) q[3];
sx q[3];
rz(0.657814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8126882) q[2];
sx q[2];
rz(-0.71949553) q[2];
sx q[2];
rz(2.1997814) q[2];
rz(-1.4224667) q[3];
sx q[3];
rz(-0.37426451) q[3];
sx q[3];
rz(-2.4678738) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0734237) q[0];
sx q[0];
rz(-0.24003679) q[0];
sx q[0];
rz(1.6085251) q[0];
rz(2.7948921) q[1];
sx q[1];
rz(-2.0997212) q[1];
sx q[1];
rz(-2.9516721) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6081734) q[0];
sx q[0];
rz(-0.77866879) q[0];
sx q[0];
rz(0.37663583) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0713351) q[2];
sx q[2];
rz(-0.53218666) q[2];
sx q[2];
rz(2.4376873) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.080481361) q[1];
sx q[1];
rz(-2.5330392) q[1];
sx q[1];
rz(-1.670865) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40855405) q[3];
sx q[3];
rz(-1.1748649) q[3];
sx q[3];
rz(0.2847288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7803663) q[2];
sx q[2];
rz(-0.85005886) q[2];
sx q[2];
rz(-0.38193199) q[2];
rz(0.10968883) q[3];
sx q[3];
rz(-2.1083125) q[3];
sx q[3];
rz(0.1194574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0383976) q[0];
sx q[0];
rz(-1.1829475) q[0];
sx q[0];
rz(-0.84684816) q[0];
rz(3.0021744) q[1];
sx q[1];
rz(-2.3250695) q[1];
sx q[1];
rz(0.74904186) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0160834) q[0];
sx q[0];
rz(-0.88403832) q[0];
sx q[0];
rz(1.3796877) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48495082) q[2];
sx q[2];
rz(-1.4451988) q[2];
sx q[2];
rz(-1.4664949) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7498503) q[1];
sx q[1];
rz(-0.23255177) q[1];
sx q[1];
rz(2.560023) q[1];
rz(-pi) q[2];
rz(-2.0723125) q[3];
sx q[3];
rz(-1.4677047) q[3];
sx q[3];
rz(-0.13971288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2060812) q[2];
sx q[2];
rz(-1.2568018) q[2];
sx q[2];
rz(1.067266) q[2];
rz(-0.71077985) q[3];
sx q[3];
rz(-1.658541) q[3];
sx q[3];
rz(1.7162286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34273219) q[0];
sx q[0];
rz(-1.6628168) q[0];
sx q[0];
rz(2.1730098) q[0];
rz(-0.69106483) q[1];
sx q[1];
rz(-1.7993118) q[1];
sx q[1];
rz(1.8341281) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.867315) q[0];
sx q[0];
rz(-2.4558892) q[0];
sx q[0];
rz(-2.7327442) q[0];
x q[1];
rz(0.2910462) q[2];
sx q[2];
rz(-1.7564764) q[2];
sx q[2];
rz(-3.1094375) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.56388748) q[1];
sx q[1];
rz(-1.1839377) q[1];
sx q[1];
rz(3.0847286) q[1];
rz(-0.18344603) q[3];
sx q[3];
rz(-1.8357092) q[3];
sx q[3];
rz(-1.9062756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4751733) q[2];
sx q[2];
rz(-0.23491493) q[2];
sx q[2];
rz(1.2972181) q[2];
rz(1.3951067) q[3];
sx q[3];
rz(-1.3336983) q[3];
sx q[3];
rz(1.5313799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.34201) q[0];
sx q[0];
rz(-0.81556773) q[0];
sx q[0];
rz(2.2014501) q[0];
rz(2.4932585) q[1];
sx q[1];
rz(-1.2038566) q[1];
sx q[1];
rz(-0.26888332) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20737442) q[0];
sx q[0];
rz(-1.3257799) q[0];
sx q[0];
rz(-0.55565683) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36264561) q[2];
sx q[2];
rz(-0.73604903) q[2];
sx q[2];
rz(0.48378746) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9098879) q[1];
sx q[1];
rz(-0.99867601) q[1];
sx q[1];
rz(-2.8673299) q[1];
rz(-0.42952764) q[3];
sx q[3];
rz(-1.2614354) q[3];
sx q[3];
rz(-1.7842899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0073283422) q[2];
sx q[2];
rz(-1.2790044) q[2];
sx q[2];
rz(-2.4556665) q[2];
rz(-2.4053597) q[3];
sx q[3];
rz(-1.0400306) q[3];
sx q[3];
rz(-2.9161016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69970423) q[0];
sx q[0];
rz(-1.2208953) q[0];
sx q[0];
rz(-2.3928483) q[0];
rz(-3.0486095) q[1];
sx q[1];
rz(-0.75983202) q[1];
sx q[1];
rz(-2.0704796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94959983) q[0];
sx q[0];
rz(-1.1566391) q[0];
sx q[0];
rz(0.55121514) q[0];
rz(-pi) q[1];
rz(1.1125715) q[2];
sx q[2];
rz(-2.1734218) q[2];
sx q[2];
rz(-2.9650709) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0784356) q[1];
sx q[1];
rz(-0.74045762) q[1];
sx q[1];
rz(-2.034212) q[1];
x q[2];
rz(-1.981826) q[3];
sx q[3];
rz(-1.618652) q[3];
sx q[3];
rz(-1.6896629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.45212713) q[2];
sx q[2];
rz(-0.63673821) q[2];
sx q[2];
rz(-2.9316736) q[2];
rz(-0.12957761) q[3];
sx q[3];
rz(-0.94730535) q[3];
sx q[3];
rz(1.5642222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0196446) q[0];
sx q[0];
rz(-0.26109281) q[0];
sx q[0];
rz(0.013539465) q[0];
rz(0.53567046) q[1];
sx q[1];
rz(-2.5745013) q[1];
sx q[1];
rz(-0.013669107) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5608198) q[0];
sx q[0];
rz(-0.81308621) q[0];
sx q[0];
rz(1.8845425) q[0];
rz(-2.057836) q[2];
sx q[2];
rz(-1.9197316) q[2];
sx q[2];
rz(-1.2020122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0108384) q[1];
sx q[1];
rz(-1.2474244) q[1];
sx q[1];
rz(2.4107928) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7444856) q[3];
sx q[3];
rz(-2.350961) q[3];
sx q[3];
rz(1.1943898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51836625) q[2];
sx q[2];
rz(-1.7609111) q[2];
sx q[2];
rz(1.7784485) q[2];
rz(0.81950435) q[3];
sx q[3];
rz(-0.47911152) q[3];
sx q[3];
rz(-0.68377408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.118947) q[0];
sx q[0];
rz(-0.85856694) q[0];
sx q[0];
rz(-1.6534506) q[0];
rz(-2.7986774) q[1];
sx q[1];
rz(-1.5838793) q[1];
sx q[1];
rz(0.75103474) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2011178) q[0];
sx q[0];
rz(-1.507326) q[0];
sx q[0];
rz(-1.2918143) q[0];
rz(1.1840759) q[2];
sx q[2];
rz(-2.6474075) q[2];
sx q[2];
rz(1.7750211) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6760867) q[1];
sx q[1];
rz(-1.4526396) q[1];
sx q[1];
rz(0.60176577) q[1];
rz(-pi) q[2];
rz(0.033174935) q[3];
sx q[3];
rz(-1.6609332) q[3];
sx q[3];
rz(-1.4226923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.47120961) q[2];
sx q[2];
rz(-2.1977916) q[2];
sx q[2];
rz(-0.38718265) q[2];
rz(-2.6662628) q[3];
sx q[3];
rz(-1.9742222) q[3];
sx q[3];
rz(-2.3502684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.249007) q[0];
sx q[0];
rz(-2.1795166) q[0];
sx q[0];
rz(-2.4695061) q[0];
rz(-2.7723906) q[1];
sx q[1];
rz(-0.75405706) q[1];
sx q[1];
rz(-1.3934607) q[1];
rz(1.3849707) q[2];
sx q[2];
rz(-2.2792049) q[2];
sx q[2];
rz(-0.21273108) q[2];
rz(1.8440856) q[3];
sx q[3];
rz(-1.4617625) q[3];
sx q[3];
rz(1.1980496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

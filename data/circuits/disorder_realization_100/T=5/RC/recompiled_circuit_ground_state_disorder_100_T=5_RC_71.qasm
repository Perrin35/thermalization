OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72803175) q[0];
sx q[0];
rz(-0.95946884) q[0];
sx q[0];
rz(2.5897107) q[0];
rz(2.7565487) q[1];
sx q[1];
rz(-1.7793964) q[1];
sx q[1];
rz(2.722932) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075567632) q[0];
sx q[0];
rz(-0.93187823) q[0];
sx q[0];
rz(-0.35982168) q[0];
rz(-pi) q[1];
rz(0.44960449) q[2];
sx q[2];
rz(-1.7309963) q[2];
sx q[2];
rz(2.2267603) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4837515) q[1];
sx q[1];
rz(-2.9194909) q[1];
sx q[1];
rz(-2.7276843) q[1];
rz(-pi) q[2];
rz(2.0830293) q[3];
sx q[3];
rz(-0.034907428) q[3];
sx q[3];
rz(1.3024746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.60575134) q[2];
sx q[2];
rz(-1.8391106) q[2];
sx q[2];
rz(0.38899404) q[2];
rz(-2.244106) q[3];
sx q[3];
rz(-0.56098452) q[3];
sx q[3];
rz(1.2787904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2314583) q[0];
sx q[0];
rz(-0.95160216) q[0];
sx q[0];
rz(0.32904539) q[0];
rz(-1.344205) q[1];
sx q[1];
rz(-2.3842594) q[1];
sx q[1];
rz(-0.12408852) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1162001) q[0];
sx q[0];
rz(-1.7804464) q[0];
sx q[0];
rz(-1.8127182) q[0];
rz(0.28606881) q[2];
sx q[2];
rz(-1.4382312) q[2];
sx q[2];
rz(2.2372467) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8777496) q[1];
sx q[1];
rz(-1.9274492) q[1];
sx q[1];
rz(1.0984196) q[1];
rz(1.8120519) q[3];
sx q[3];
rz(-0.71863758) q[3];
sx q[3];
rz(2.1332316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8590392) q[2];
sx q[2];
rz(-2.6586847) q[2];
sx q[2];
rz(-0.60749751) q[2];
rz(-2.2533158) q[3];
sx q[3];
rz(-1.5253303) q[3];
sx q[3];
rz(-1.7150778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4597976) q[0];
sx q[0];
rz(-1.8495704) q[0];
sx q[0];
rz(2.9070396) q[0];
rz(1.0598496) q[1];
sx q[1];
rz(-1.9748961) q[1];
sx q[1];
rz(-0.92811981) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7674196) q[0];
sx q[0];
rz(-2.1819127) q[0];
sx q[0];
rz(1.060063) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4305827) q[2];
sx q[2];
rz(-0.59525049) q[2];
sx q[2];
rz(-0.75314116) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.69423543) q[1];
sx q[1];
rz(-0.35590812) q[1];
sx q[1];
rz(0.39519542) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9506128) q[3];
sx q[3];
rz(-0.39576021) q[3];
sx q[3];
rz(0.8518962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5633391) q[2];
sx q[2];
rz(-1.4207062) q[2];
sx q[2];
rz(1.2467747) q[2];
rz(2.9747544) q[3];
sx q[3];
rz(-0.92883795) q[3];
sx q[3];
rz(-1.5290574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.700915) q[0];
sx q[0];
rz(-3.0755141) q[0];
sx q[0];
rz(-0.75827688) q[0];
rz(1.5658763) q[1];
sx q[1];
rz(-0.58454746) q[1];
sx q[1];
rz(-0.87055269) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93773952) q[0];
sx q[0];
rz(-1.1273307) q[0];
sx q[0];
rz(-0.70758836) q[0];
rz(-pi) q[1];
rz(2.5640268) q[2];
sx q[2];
rz(-1.1223167) q[2];
sx q[2];
rz(0.73671337) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4728001) q[1];
sx q[1];
rz(-1.1684946) q[1];
sx q[1];
rz(-0.16568664) q[1];
rz(-pi) q[2];
rz(2.5136886) q[3];
sx q[3];
rz(-1.6293173) q[3];
sx q[3];
rz(-2.1128775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42252758) q[2];
sx q[2];
rz(-1.0079404) q[2];
sx q[2];
rz(-0.75971216) q[2];
rz(-2.4921913) q[3];
sx q[3];
rz(-1.6067182) q[3];
sx q[3];
rz(-1.5627129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11288697) q[0];
sx q[0];
rz(-0.6539456) q[0];
sx q[0];
rz(-1.9586067) q[0];
rz(0.68823367) q[1];
sx q[1];
rz(-1.3166683) q[1];
sx q[1];
rz(-2.7722955) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0205905) q[0];
sx q[0];
rz(-2.2341616) q[0];
sx q[0];
rz(-1.6055907) q[0];
rz(-pi) q[1];
rz(2.4299116) q[2];
sx q[2];
rz(-1.7717012) q[2];
sx q[2];
rz(0.37362675) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0329513) q[1];
sx q[1];
rz(-0.88283112) q[1];
sx q[1];
rz(2.9017519) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7820205) q[3];
sx q[3];
rz(-0.45392515) q[3];
sx q[3];
rz(-0.63185121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1987622) q[2];
sx q[2];
rz(-1.8161215) q[2];
sx q[2];
rz(-2.14373) q[2];
rz(-2.5241847) q[3];
sx q[3];
rz(-0.95881763) q[3];
sx q[3];
rz(-0.82120419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1029516) q[0];
sx q[0];
rz(-1.8674253) q[0];
sx q[0];
rz(-1.8983023) q[0];
rz(0.78168166) q[1];
sx q[1];
rz(-1.8676753) q[1];
sx q[1];
rz(2.8807358) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92453018) q[0];
sx q[0];
rz(-2.6017761) q[0];
sx q[0];
rz(-2.7227719) q[0];
rz(-pi) q[1];
rz(0.86791188) q[2];
sx q[2];
rz(-2.5943668) q[2];
sx q[2];
rz(-2.1457714) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8159224) q[1];
sx q[1];
rz(-1.7223893) q[1];
sx q[1];
rz(1.1071015) q[1];
rz(-pi) q[2];
rz(-1.8571768) q[3];
sx q[3];
rz(-2.4900511) q[3];
sx q[3];
rz(-2.1809354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6692052) q[2];
sx q[2];
rz(-0.86296764) q[2];
sx q[2];
rz(2.6589987) q[2];
rz(-0.11416301) q[3];
sx q[3];
rz(-2.0518905) q[3];
sx q[3];
rz(-0.94858661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77940762) q[0];
sx q[0];
rz(-0.055883378) q[0];
sx q[0];
rz(2.8756397) q[0];
rz(-2.7159122) q[1];
sx q[1];
rz(-1.2454147) q[1];
sx q[1];
rz(-0.46806213) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45586203) q[0];
sx q[0];
rz(-2.1121368) q[0];
sx q[0];
rz(-0.0079519072) q[0];
x q[1];
rz(-1.8276026) q[2];
sx q[2];
rz(-2.7903284) q[2];
sx q[2];
rz(2.6333269) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.148273) q[1];
sx q[1];
rz(-0.56374218) q[1];
sx q[1];
rz(-1.0757273) q[1];
rz(-pi) q[2];
rz(-3.1105177) q[3];
sx q[3];
rz(-1.0508176) q[3];
sx q[3];
rz(-1.7620373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.88428664) q[2];
sx q[2];
rz(-1.851119) q[2];
sx q[2];
rz(0.5536983) q[2];
rz(0.92343679) q[3];
sx q[3];
rz(-0.90673509) q[3];
sx q[3];
rz(-0.060801774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4629352) q[0];
sx q[0];
rz(-2.904992) q[0];
sx q[0];
rz(-0.67805725) q[0];
rz(0.97738114) q[1];
sx q[1];
rz(-1.1481608) q[1];
sx q[1];
rz(1.4378907) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17409731) q[0];
sx q[0];
rz(-0.50104841) q[0];
sx q[0];
rz(1.9809841) q[0];
rz(-3.1372141) q[2];
sx q[2];
rz(-2.5776049) q[2];
sx q[2];
rz(2.0019238) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4297144) q[1];
sx q[1];
rz(-2.7949882) q[1];
sx q[1];
rz(2.4019353) q[1];
rz(-pi) q[2];
rz(0.57084031) q[3];
sx q[3];
rz(-2.3085706) q[3];
sx q[3];
rz(0.86467165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.73072726) q[2];
sx q[2];
rz(-2.5812456) q[2];
sx q[2];
rz(0.66696683) q[2];
rz(-3.0242331) q[3];
sx q[3];
rz(-1.8662235) q[3];
sx q[3];
rz(-1.3807152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19408195) q[0];
sx q[0];
rz(-1.5912594) q[0];
sx q[0];
rz(2.7981753) q[0];
rz(1.6471479) q[1];
sx q[1];
rz(-0.98973715) q[1];
sx q[1];
rz(2.6572878) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13499943) q[0];
sx q[0];
rz(-1.3065225) q[0];
sx q[0];
rz(2.5359383) q[0];
x q[1];
rz(-0.42895173) q[2];
sx q[2];
rz(-1.9520734) q[2];
sx q[2];
rz(-2.5457053) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8645939) q[1];
sx q[1];
rz(-0.62886695) q[1];
sx q[1];
rz(0.38284812) q[1];
rz(-pi) q[2];
rz(-2.8922594) q[3];
sx q[3];
rz(-1.3568546) q[3];
sx q[3];
rz(0.46505022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1962428) q[2];
sx q[2];
rz(-1.1002772) q[2];
sx q[2];
rz(-0.84570447) q[2];
rz(-1.2196994) q[3];
sx q[3];
rz(-1.2518576) q[3];
sx q[3];
rz(0.20488258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4591111) q[0];
sx q[0];
rz(-2.8963608) q[0];
sx q[0];
rz(-2.5526175) q[0];
rz(-2.466195) q[1];
sx q[1];
rz(-2.1928936) q[1];
sx q[1];
rz(1.677547) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9927514) q[0];
sx q[0];
rz(-2.0939079) q[0];
sx q[0];
rz(2.6639725) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6590933) q[2];
sx q[2];
rz(-0.73511926) q[2];
sx q[2];
rz(0.59150782) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3879536) q[1];
sx q[1];
rz(-1.5245226) q[1];
sx q[1];
rz(-0.58677499) q[1];
rz(-pi) q[2];
rz(-2.0328224) q[3];
sx q[3];
rz(-2.2115876) q[3];
sx q[3];
rz(-2.6907211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7318763) q[2];
sx q[2];
rz(-1.3940553) q[2];
sx q[2];
rz(-2.0173006) q[2];
rz(-0.8479979) q[3];
sx q[3];
rz(-0.8046937) q[3];
sx q[3];
rz(0.030357411) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.547434) q[0];
sx q[0];
rz(-2.0068598) q[0];
sx q[0];
rz(1.6225847) q[0];
rz(2.2183954) q[1];
sx q[1];
rz(-2.1122439) q[1];
sx q[1];
rz(-1.5069638) q[1];
rz(-1.3849003) q[2];
sx q[2];
rz(-1.8757314) q[2];
sx q[2];
rz(1.3396946) q[2];
rz(2.2062929) q[3];
sx q[3];
rz(-0.37105303) q[3];
sx q[3];
rz(1.1340352) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(0.49864545) q[0];
sx q[0];
rz(3.6874229) q[0];
sx q[0];
rz(9.3064718) q[0];
rz(0.91953295) q[1];
sx q[1];
rz(-1.8383205) q[1];
sx q[1];
rz(-1.2143171) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8144381) q[0];
sx q[0];
rz(-2.0025578) q[0];
sx q[0];
rz(-0.23349725) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75096424) q[2];
sx q[2];
rz(-2.3625824) q[2];
sx q[2];
rz(2.2476803) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.37049261) q[1];
sx q[1];
rz(-2.1861861) q[1];
sx q[1];
rz(2.9740488) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2976152) q[3];
sx q[3];
rz(-1.5500549) q[3];
sx q[3];
rz(2.1368049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.168657) q[2];
sx q[2];
rz(-1.1309705) q[2];
sx q[2];
rz(-1.0966148) q[2];
rz(-0.24122572) q[3];
sx q[3];
rz(-1.3696407) q[3];
sx q[3];
rz(0.84478861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089652561) q[0];
sx q[0];
rz(-2.0515428) q[0];
sx q[0];
rz(0.3845149) q[0];
rz(-2.8334726) q[1];
sx q[1];
rz(-0.48079753) q[1];
sx q[1];
rz(-0.62517977) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2848672) q[0];
sx q[0];
rz(-1.2846886) q[0];
sx q[0];
rz(-1.6038143) q[0];
rz(-pi) q[1];
rz(2.9956498) q[2];
sx q[2];
rz(-2.4193442) q[2];
sx q[2];
rz(-0.099522245) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2634266) q[1];
sx q[1];
rz(-1.485865) q[1];
sx q[1];
rz(-2.4919266) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6513807) q[3];
sx q[3];
rz(-1.0102838) q[3];
sx q[3];
rz(2.0106767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.008931) q[2];
sx q[2];
rz(-2.2098358) q[2];
sx q[2];
rz(-2.8399732) q[2];
rz(-0.53145069) q[3];
sx q[3];
rz(-0.88574946) q[3];
sx q[3];
rz(2.0187995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5136435) q[0];
sx q[0];
rz(-2.1301837) q[0];
sx q[0];
rz(-1.1847786) q[0];
rz(-0.14492497) q[1];
sx q[1];
rz(-1.8791608) q[1];
sx q[1];
rz(-1.5941934) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14159053) q[0];
sx q[0];
rz(-0.98391082) q[0];
sx q[0];
rz(-0.97372719) q[0];
x q[1];
rz(0.1833488) q[2];
sx q[2];
rz(-1.9446704) q[2];
sx q[2];
rz(-3.0484859) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6559534) q[1];
sx q[1];
rz(-1.1879293) q[1];
sx q[1];
rz(-2.1583771) q[1];
rz(-1.8933442) q[3];
sx q[3];
rz(-1.9574021) q[3];
sx q[3];
rz(0.69535386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3323815) q[2];
sx q[2];
rz(-0.96497649) q[2];
sx q[2];
rz(-1.8154209) q[2];
rz(2.2281846) q[3];
sx q[3];
rz(-1.2181506) q[3];
sx q[3];
rz(2.6083045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7186385) q[0];
sx q[0];
rz(-2.664743) q[0];
sx q[0];
rz(1.3612716) q[0];
rz(2.4286229) q[1];
sx q[1];
rz(-2.0013335) q[1];
sx q[1];
rz(0.68351173) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26091247) q[0];
sx q[0];
rz(-2.144038) q[0];
sx q[0];
rz(-3.1295958) q[0];
rz(-pi) q[1];
rz(-2.9320549) q[2];
sx q[2];
rz(-1.7895797) q[2];
sx q[2];
rz(-2.0444586) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5665413) q[1];
sx q[1];
rz(-1.05541) q[1];
sx q[1];
rz(-0.75830663) q[1];
rz(-pi) q[2];
rz(1.470742) q[3];
sx q[3];
rz(-1.3276427) q[3];
sx q[3];
rz(0.2167165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6354562) q[2];
sx q[2];
rz(-2.7244302) q[2];
sx q[2];
rz(-0.48640856) q[2];
rz(0.5168612) q[3];
sx q[3];
rz(-1.9603399) q[3];
sx q[3];
rz(1.0335056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085623398) q[0];
sx q[0];
rz(-1.509868) q[0];
sx q[0];
rz(-1.7686718) q[0];
rz(1.4068475) q[1];
sx q[1];
rz(-0.56766784) q[1];
sx q[1];
rz(-0.47703201) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5236899) q[0];
sx q[0];
rz(-1.7672024) q[0];
sx q[0];
rz(1.1454885) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7482175) q[2];
sx q[2];
rz(-2.057339) q[2];
sx q[2];
rz(-1.7560619) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.43852311) q[1];
sx q[1];
rz(-2.065899) q[1];
sx q[1];
rz(0.34214051) q[1];
rz(2.3588657) q[3];
sx q[3];
rz(-2.1313792) q[3];
sx q[3];
rz(2.4334139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.10739747) q[2];
sx q[2];
rz(-0.98661462) q[2];
sx q[2];
rz(0.46727115) q[2];
rz(-0.16921903) q[3];
sx q[3];
rz(-1.7305814) q[3];
sx q[3];
rz(-0.30538487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53448236) q[0];
sx q[0];
rz(-0.87962532) q[0];
sx q[0];
rz(0.18950263) q[0];
rz(-2.8684008) q[1];
sx q[1];
rz(-1.0739645) q[1];
sx q[1];
rz(2.3409519) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55274945) q[0];
sx q[0];
rz(-0.23736033) q[0];
sx q[0];
rz(1.3861603) q[0];
rz(-pi) q[1];
rz(0.34363644) q[2];
sx q[2];
rz(-2.8784907) q[2];
sx q[2];
rz(-3.0575036) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8756518) q[1];
sx q[1];
rz(-1.3752626) q[1];
sx q[1];
rz(0.13130782) q[1];
x q[2];
rz(0.54038911) q[3];
sx q[3];
rz(-2.0107993) q[3];
sx q[3];
rz(2.7203181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.93929401) q[2];
sx q[2];
rz(-2.9677128) q[2];
sx q[2];
rz(2.6046216) q[2];
rz(1.2823074) q[3];
sx q[3];
rz(-1.8351646) q[3];
sx q[3];
rz(-1.4793226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.031161664) q[0];
sx q[0];
rz(-0.81542504) q[0];
sx q[0];
rz(1.8102616) q[0];
rz(1.7000465) q[1];
sx q[1];
rz(-1.6621637) q[1];
sx q[1];
rz(1.2225245) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13200296) q[0];
sx q[0];
rz(-1.9870269) q[0];
sx q[0];
rz(2.8819109) q[0];
rz(-pi) q[1];
rz(-1.6454846) q[2];
sx q[2];
rz(-2.1213581) q[2];
sx q[2];
rz(-2.2979743) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8795506) q[1];
sx q[1];
rz(-0.73644222) q[1];
sx q[1];
rz(1.8950736) q[1];
x q[2];
rz(3.0223373) q[3];
sx q[3];
rz(-1.1290765) q[3];
sx q[3];
rz(2.2252803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20087251) q[2];
sx q[2];
rz(-0.49786374) q[2];
sx q[2];
rz(-2.8087924) q[2];
rz(2.8424272) q[3];
sx q[3];
rz(-1.2406415) q[3];
sx q[3];
rz(3.1201194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12896319) q[0];
sx q[0];
rz(-0.7681995) q[0];
sx q[0];
rz(0.27498883) q[0];
rz(-3.0601652) q[1];
sx q[1];
rz(-0.80137253) q[1];
sx q[1];
rz(-1.8191232) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0596432) q[0];
sx q[0];
rz(-2.9026051) q[0];
sx q[0];
rz(-1.3168524) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1312348) q[2];
sx q[2];
rz(-0.73652525) q[2];
sx q[2];
rz(0.16939746) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7863955) q[1];
sx q[1];
rz(-1.9405559) q[1];
sx q[1];
rz(-0.74798297) q[1];
rz(-pi) q[2];
rz(-2.4422936) q[3];
sx q[3];
rz(-0.85996503) q[3];
sx q[3];
rz(0.028472326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.72121173) q[2];
sx q[2];
rz(-3.0215441) q[2];
sx q[2];
rz(-1.4018641) q[2];
rz(-1.2841355) q[3];
sx q[3];
rz(-1.9522342) q[3];
sx q[3];
rz(3.106626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5818802) q[0];
sx q[0];
rz(-1.548883) q[0];
sx q[0];
rz(-0.48932073) q[0];
rz(3.1210461) q[1];
sx q[1];
rz(-1.9053883) q[1];
sx q[1];
rz(2.4403341) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72557455) q[0];
sx q[0];
rz(-1.3301992) q[0];
sx q[0];
rz(-0.98576905) q[0];
rz(1.5145153) q[2];
sx q[2];
rz(-2.4372134) q[2];
sx q[2];
rz(-3.0075108) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.012246) q[1];
sx q[1];
rz(-1.3353026) q[1];
sx q[1];
rz(0.78462623) q[1];
rz(-pi) q[2];
rz(0.68655007) q[3];
sx q[3];
rz(-1.1995763) q[3];
sx q[3];
rz(0.29774259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7265085) q[2];
sx q[2];
rz(-2.5231611) q[2];
sx q[2];
rz(2.6611967) q[2];
rz(-1.7324309) q[3];
sx q[3];
rz(-1.3536072) q[3];
sx q[3];
rz(2.8044146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3163863) q[0];
sx q[0];
rz(-1.9012863) q[0];
sx q[0];
rz(0.41148841) q[0];
rz(1.2443789) q[1];
sx q[1];
rz(-2.0020516) q[1];
sx q[1];
rz(-2.8172475) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1169057) q[0];
sx q[0];
rz(-1.9431264) q[0];
sx q[0];
rz(0.056554746) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3408808) q[2];
sx q[2];
rz(-0.47242785) q[2];
sx q[2];
rz(2.1683482) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0592546) q[1];
sx q[1];
rz(-1.6538701) q[1];
sx q[1];
rz(2.1359813) q[1];
x q[2];
rz(-3.0203811) q[3];
sx q[3];
rz(-1.5577661) q[3];
sx q[3];
rz(0.37076326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9559429) q[2];
sx q[2];
rz(-2.0273384) q[2];
sx q[2];
rz(-3.0246217) q[2];
rz(2.1864435) q[3];
sx q[3];
rz(-0.17894608) q[3];
sx q[3];
rz(0.54168934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87880001) q[0];
sx q[0];
rz(-2.265082) q[0];
sx q[0];
rz(-2.2646917) q[0];
rz(-2.0211438) q[1];
sx q[1];
rz(-0.81605492) q[1];
sx q[1];
rz(1.1647404) q[1];
rz(-2.0777069) q[2];
sx q[2];
rz(-1.3851266) q[2];
sx q[2];
rz(1.2141488) q[2];
rz(-1.4799812) q[3];
sx q[3];
rz(-0.27169966) q[3];
sx q[3];
rz(-3.0063831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(2.6440808) q[0];
sx q[0];
rz(-1.3389791) q[0];
sx q[0];
rz(-0.31565491) q[0];
rz(1.1701801) q[1];
sx q[1];
rz(-0.38771954) q[1];
sx q[1];
rz(0.60400909) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3919298) q[0];
sx q[0];
rz(-1.6214341) q[0];
sx q[0];
rz(-1.2681566) q[0];
rz(-pi) q[1];
rz(2.6409615) q[2];
sx q[2];
rz(-0.96304578) q[2];
sx q[2];
rz(1.0333824) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4872514) q[1];
sx q[1];
rz(-2.5120334) q[1];
sx q[1];
rz(0.97626026) q[1];
rz(-0.1084566) q[3];
sx q[3];
rz(-1.8780969) q[3];
sx q[3];
rz(1.9815597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.151256) q[2];
sx q[2];
rz(-1.9467111) q[2];
sx q[2];
rz(1.7830431) q[2];
rz(-1.547706) q[3];
sx q[3];
rz(-2.2113776) q[3];
sx q[3];
rz(1.5096629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58594054) q[0];
sx q[0];
rz(-2.5037615) q[0];
sx q[0];
rz(1.9441388) q[0];
rz(-2.6400631) q[1];
sx q[1];
rz(-1.5781559) q[1];
sx q[1];
rz(1.4362358) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1284416) q[0];
sx q[0];
rz(-1.3396946) q[0];
sx q[0];
rz(2.1559155) q[0];
rz(-pi) q[1];
rz(-0.70296944) q[2];
sx q[2];
rz(-1.6466993) q[2];
sx q[2];
rz(0.10026201) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2982218) q[1];
sx q[1];
rz(-2.9195115) q[1];
sx q[1];
rz(2.6574617) q[1];
rz(-0.40879876) q[3];
sx q[3];
rz(-0.49064562) q[3];
sx q[3];
rz(1.7237494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2531835) q[2];
sx q[2];
rz(-1.9056355) q[2];
sx q[2];
rz(0.02296981) q[2];
rz(-2.2394771) q[3];
sx q[3];
rz(-1.2227367) q[3];
sx q[3];
rz(-2.7149916) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4334634) q[0];
sx q[0];
rz(-1.3490278) q[0];
sx q[0];
rz(-2.1500812) q[0];
rz(-1.0333215) q[1];
sx q[1];
rz(-0.28305498) q[1];
sx q[1];
rz(1.2352157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69312364) q[0];
sx q[0];
rz(-2.4938994) q[0];
sx q[0];
rz(-0.41843398) q[0];
rz(-pi) q[1];
rz(-0.60977817) q[2];
sx q[2];
rz(-1.9710566) q[2];
sx q[2];
rz(-0.61851172) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24656235) q[1];
sx q[1];
rz(-0.8324648) q[1];
sx q[1];
rz(0.33729302) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3091812) q[3];
sx q[3];
rz(-2.4042272) q[3];
sx q[3];
rz(-2.5016145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.027355) q[2];
sx q[2];
rz(-2.3894252) q[2];
sx q[2];
rz(-0.99119622) q[2];
rz(1.8441955) q[3];
sx q[3];
rz(-1.2605366) q[3];
sx q[3];
rz(-1.4170925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7436413) q[0];
sx q[0];
rz(-2.209111) q[0];
sx q[0];
rz(0.81746307) q[0];
rz(-0.3262597) q[1];
sx q[1];
rz(-1.9605109) q[1];
sx q[1];
rz(0.78525966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0164675) q[0];
sx q[0];
rz(-0.40832357) q[0];
sx q[0];
rz(-0.84618469) q[0];
x q[1];
rz(-2.8815124) q[2];
sx q[2];
rz(-0.65979119) q[2];
sx q[2];
rz(1.545797) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37267673) q[1];
sx q[1];
rz(-2.0119889) q[1];
sx q[1];
rz(-0.4695453) q[1];
x q[2];
rz(2.7866632) q[3];
sx q[3];
rz(-1.4961492) q[3];
sx q[3];
rz(1.9982893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.12104812) q[2];
sx q[2];
rz(-0.64457568) q[2];
sx q[2];
rz(2.5784967) q[2];
rz(2.2480887) q[3];
sx q[3];
rz(-1.1734791) q[3];
sx q[3];
rz(1.2347429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.23984443) q[0];
sx q[0];
rz(-2.8185066) q[0];
sx q[0];
rz(1.595994) q[0];
rz(-2.1196938) q[1];
sx q[1];
rz(-1.1232168) q[1];
sx q[1];
rz(-0.18128577) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2342628) q[0];
sx q[0];
rz(-2.0600304) q[0];
sx q[0];
rz(0.2184043) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0490169) q[2];
sx q[2];
rz(-0.47738722) q[2];
sx q[2];
rz(0.30711781) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.1039155) q[1];
sx q[1];
rz(-1.5427151) q[1];
sx q[1];
rz(1.3485391) q[1];
x q[2];
rz(-0.44293176) q[3];
sx q[3];
rz(-1.4287523) q[3];
sx q[3];
rz(-2.1403109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.01866092) q[2];
sx q[2];
rz(-2.4666726) q[2];
sx q[2];
rz(1.8956511) q[2];
rz(0.59257007) q[3];
sx q[3];
rz(-1.6962681) q[3];
sx q[3];
rz(2.3004801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5331921) q[0];
sx q[0];
rz(-2.1493981) q[0];
sx q[0];
rz(0.30884185) q[0];
rz(-2.8187075) q[1];
sx q[1];
rz(-0.37728089) q[1];
sx q[1];
rz(3.0016532) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2353781) q[0];
sx q[0];
rz(-1.3227292) q[0];
sx q[0];
rz(-0.95309044) q[0];
x q[1];
rz(0.094801589) q[2];
sx q[2];
rz(-2.0249512) q[2];
sx q[2];
rz(-0.73341441) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6778826) q[1];
sx q[1];
rz(-1.8720798) q[1];
sx q[1];
rz(0.82656411) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9512964) q[3];
sx q[3];
rz(-0.69529136) q[3];
sx q[3];
rz(-2.8075346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.66190019) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(-0.24197401) q[2];
rz(0.74609977) q[3];
sx q[3];
rz(-1.7753764) q[3];
sx q[3];
rz(-1.411875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11874966) q[0];
sx q[0];
rz(-2.0976837) q[0];
sx q[0];
rz(-1.2710849) q[0];
rz(-0.56308833) q[1];
sx q[1];
rz(-2.2124898) q[1];
sx q[1];
rz(1.44106) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4615501) q[0];
sx q[0];
rz(-0.6667887) q[0];
sx q[0];
rz(1.1229442) q[0];
rz(-0.38132847) q[2];
sx q[2];
rz(-1.6994119) q[2];
sx q[2];
rz(0.53024697) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.064441) q[1];
sx q[1];
rz(-1.3485326) q[1];
sx q[1];
rz(2.7138316) q[1];
rz(-0.7251803) q[3];
sx q[3];
rz(-1.4075052) q[3];
sx q[3];
rz(2.0189328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9219804) q[2];
sx q[2];
rz(-1.8461123) q[2];
sx q[2];
rz(-1.316635) q[2];
rz(0.5611788) q[3];
sx q[3];
rz(-2.1771274) q[3];
sx q[3];
rz(1.6130028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.104326) q[0];
sx q[0];
rz(-2.9236561) q[0];
sx q[0];
rz(-2.2547146) q[0];
rz(-0.2001702) q[1];
sx q[1];
rz(-1.472241) q[1];
sx q[1];
rz(2.1194469) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8093531) q[0];
sx q[0];
rz(-0.13357559) q[0];
sx q[0];
rz(-1.3788401) q[0];
rz(-pi) q[1];
rz(-0.81424197) q[2];
sx q[2];
rz(-0.98746429) q[2];
sx q[2];
rz(0.63177201) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.31769842) q[1];
sx q[1];
rz(-0.824747) q[1];
sx q[1];
rz(0.0024248799) q[1];
rz(-pi) q[2];
rz(-0.57557801) q[3];
sx q[3];
rz(-1.8691571) q[3];
sx q[3];
rz(-2.9090867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48978051) q[2];
sx q[2];
rz(-0.45280364) q[2];
sx q[2];
rz(-2.32302) q[2];
rz(0.98322785) q[3];
sx q[3];
rz(-2.1849617) q[3];
sx q[3];
rz(0.53808588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1061358) q[0];
sx q[0];
rz(-1.1470733) q[0];
sx q[0];
rz(-1.5552833) q[0];
rz(0.69681329) q[1];
sx q[1];
rz(-2.9882444) q[1];
sx q[1];
rz(0.78479016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1589543) q[0];
sx q[0];
rz(-1.0305911) q[0];
sx q[0];
rz(1.0544974) q[0];
rz(-pi) q[1];
rz(-2.2124039) q[2];
sx q[2];
rz(-2.8904473) q[2];
sx q[2];
rz(-2.8751858) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0052629) q[1];
sx q[1];
rz(-1.6044093) q[1];
sx q[1];
rz(-0.43515794) q[1];
rz(-pi) q[2];
rz(2.5303477) q[3];
sx q[3];
rz(-1.9211624) q[3];
sx q[3];
rz(-2.4823852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4697504) q[2];
sx q[2];
rz(-1.8639114) q[2];
sx q[2];
rz(-2.3308241) q[2];
rz(1.9226711) q[3];
sx q[3];
rz(-2.6007077) q[3];
sx q[3];
rz(-2.0764652) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3588381) q[0];
sx q[0];
rz(-0.29958075) q[0];
sx q[0];
rz(-1.4240356) q[0];
rz(-2.3639823) q[1];
sx q[1];
rz(-2.168226) q[1];
sx q[1];
rz(-2.3861859) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41235218) q[0];
sx q[0];
rz(-1.5662114) q[0];
sx q[0];
rz(1.7807177) q[0];
x q[1];
rz(1.5883114) q[2];
sx q[2];
rz(-1.1743769) q[2];
sx q[2];
rz(-2.6902386) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.091614231) q[1];
sx q[1];
rz(-2.0777933) q[1];
sx q[1];
rz(1.1683639) q[1];
rz(-pi) q[2];
rz(-1.8978664) q[3];
sx q[3];
rz(-0.9216412) q[3];
sx q[3];
rz(2.7130733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5648254) q[2];
sx q[2];
rz(-0.29594031) q[2];
sx q[2];
rz(-2.9887065) q[2];
rz(1.2711924) q[3];
sx q[3];
rz(-1.8891687) q[3];
sx q[3];
rz(-0.66711867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85970238) q[0];
sx q[0];
rz(-2.6489881) q[0];
sx q[0];
rz(-0.76982605) q[0];
rz(1.9181171) q[1];
sx q[1];
rz(-1.2212831) q[1];
sx q[1];
rz(2.4881359) q[1];
rz(-0.59719795) q[2];
sx q[2];
rz(-1.2033249) q[2];
sx q[2];
rz(-2.0133599) q[2];
rz(-0.35265233) q[3];
sx q[3];
rz(-1.9518387) q[3];
sx q[3];
rz(0.52458682) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

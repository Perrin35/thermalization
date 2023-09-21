OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9632602) q[0];
sx q[0];
rz(-1.652521) q[0];
sx q[0];
rz(0.89515495) q[0];
rz(2.826638) q[1];
sx q[1];
rz(-1.0840253) q[1];
sx q[1];
rz(1.6853583) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.668894) q[0];
sx q[0];
rz(-1.5905252) q[0];
sx q[0];
rz(-1.701645) q[0];
x q[1];
rz(-0.99631359) q[2];
sx q[2];
rz(-2.516054) q[2];
sx q[2];
rz(-1.1686981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.284277) q[1];
sx q[1];
rz(-0.89968649) q[1];
sx q[1];
rz(1.3778694) q[1];
rz(-pi) q[2];
x q[2];
rz(1.625657) q[3];
sx q[3];
rz(-1.8901955) q[3];
sx q[3];
rz(1.1064135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.75498092) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(-2.6888729) q[2];
rz(0.1581986) q[3];
sx q[3];
rz(-0.69163624) q[3];
sx q[3];
rz(2.246777) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92900705) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(1.1520977) q[0];
rz(-1.903803) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(-0.47098413) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5574871) q[0];
sx q[0];
rz(-0.4852681) q[0];
sx q[0];
rz(1.6886061) q[0];
rz(-0.97371308) q[2];
sx q[2];
rz(-1.7657585) q[2];
sx q[2];
rz(-2.3174469) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2742548) q[1];
sx q[1];
rz(-1.5161637) q[1];
sx q[1];
rz(2.0591303) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4906293) q[3];
sx q[3];
rz(-2.2647694) q[3];
sx q[3];
rz(2.699664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.058078893) q[2];
sx q[2];
rz(-0.60516548) q[2];
sx q[2];
rz(-2.8857968) q[2];
rz(-1.4852218) q[3];
sx q[3];
rz(-1.9519613) q[3];
sx q[3];
rz(-0.16168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8783022) q[0];
sx q[0];
rz(-2.0354164) q[0];
sx q[0];
rz(0.27134744) q[0];
rz(0.73633206) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(0.43930611) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8205748) q[0];
sx q[0];
rz(-1.6512172) q[0];
sx q[0];
rz(0.029650173) q[0];
rz(2.7810532) q[2];
sx q[2];
rz(-1.6191102) q[2];
sx q[2];
rz(2.609032) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15553741) q[1];
sx q[1];
rz(-1.7441161) q[1];
sx q[1];
rz(-0.9539414) q[1];
x q[2];
rz(-2.1268232) q[3];
sx q[3];
rz(-1.945567) q[3];
sx q[3];
rz(2.8127363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1674041) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(0.20047323) q[2];
rz(0.75508562) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(-2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.077483594) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(-1.2444929) q[0];
rz(-2.3311133) q[1];
sx q[1];
rz(-1.8119718) q[1];
sx q[1];
rz(-2.2669852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0229605) q[0];
sx q[0];
rz(-1.7264139) q[0];
sx q[0];
rz(-2.1778584) q[0];
rz(-0.34747296) q[2];
sx q[2];
rz(-2.7345737) q[2];
sx q[2];
rz(-2.4328872) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.306327) q[1];
sx q[1];
rz(-1.3439461) q[1];
sx q[1];
rz(2.7579685) q[1];
x q[2];
rz(-3.0365385) q[3];
sx q[3];
rz(-1.4964364) q[3];
sx q[3];
rz(-2.8021333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.13742927) q[2];
sx q[2];
rz(-0.88976294) q[2];
sx q[2];
rz(1.7144263) q[2];
rz(-0.066120474) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(1.6920413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604597) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(-0.57058913) q[0];
rz(-2.5866306) q[1];
sx q[1];
rz(-0.73736063) q[1];
sx q[1];
rz(0.74329174) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76737228) q[0];
sx q[0];
rz(-1.7905856) q[0];
sx q[0];
rz(0.90669294) q[0];
x q[1];
rz(1.4666124) q[2];
sx q[2];
rz(-1.0717234) q[2];
sx q[2];
rz(-2.4521329) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44504657) q[1];
sx q[1];
rz(-1.7940709) q[1];
sx q[1];
rz(2.9162507) q[1];
rz(0.27637847) q[3];
sx q[3];
rz(-0.82377269) q[3];
sx q[3];
rz(-2.7793022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9479998) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(2.0949481) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(2.823765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.627581) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(-2.6348689) q[0];
rz(2.8880033) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(-2.0862897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9133638) q[0];
sx q[0];
rz(-2.1852487) q[0];
sx q[0];
rz(2.9173304) q[0];
x q[1];
rz(-2.2248473) q[2];
sx q[2];
rz(-1.9987717) q[2];
sx q[2];
rz(-0.59631729) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9678952) q[1];
sx q[1];
rz(-2.3579862) q[1];
sx q[1];
rz(-0.99131363) q[1];
rz(-pi) q[2];
rz(1.33527) q[3];
sx q[3];
rz(-0.84170656) q[3];
sx q[3];
rz(1.0585166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6598597) q[2];
sx q[2];
rz(-2.0969756) q[2];
sx q[2];
rz(-2.2591023) q[2];
rz(2.4957538) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(1.283949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362474) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(0.77254599) q[0];
rz(1.729471) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(-2.5678182) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6360146) q[0];
sx q[0];
rz(-1.4609219) q[0];
sx q[0];
rz(-1.3789603) q[0];
x q[1];
rz(-0.022106604) q[2];
sx q[2];
rz(-1.3706285) q[2];
sx q[2];
rz(-0.18891639) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.335865) q[1];
sx q[1];
rz(-2.1153567) q[1];
sx q[1];
rz(-0.90421275) q[1];
x q[2];
rz(1.3366367) q[3];
sx q[3];
rz(-1.4244716) q[3];
sx q[3];
rz(-0.22406604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1022169) q[2];
sx q[2];
rz(-0.45414671) q[2];
sx q[2];
rz(2.3708564) q[2];
rz(2.7052774) q[3];
sx q[3];
rz(-1.2687012) q[3];
sx q[3];
rz(1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.27281) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(1.9708721) q[0];
rz(-0.51013485) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(-1.2957113) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0062795) q[0];
sx q[0];
rz(-1.2426002) q[0];
sx q[0];
rz(-0.4853863) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9240575) q[2];
sx q[2];
rz(-1.3152939) q[2];
sx q[2];
rz(3.0614292) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.46591972) q[1];
sx q[1];
rz(-0.95341668) q[1];
sx q[1];
rz(-2.8663551) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6298953) q[3];
sx q[3];
rz(-0.44145465) q[3];
sx q[3];
rz(-2.3578701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7065113) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(0.56813017) q[2];
rz(-0.98012296) q[3];
sx q[3];
rz(-1.0390493) q[3];
sx q[3];
rz(-0.19395104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4762964) q[0];
sx q[0];
rz(-0.81403533) q[0];
sx q[0];
rz(2.4639159) q[0];
rz(-2.9455345) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(2.303404) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3121376) q[0];
sx q[0];
rz(-1.2078309) q[0];
sx q[0];
rz(-1.8118993) q[0];
x q[1];
rz(2.3269862) q[2];
sx q[2];
rz(-1.5216773) q[2];
sx q[2];
rz(2.5533822) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.14324489) q[1];
sx q[1];
rz(-1.76941) q[1];
sx q[1];
rz(-2.3193588) q[1];
rz(-pi) q[2];
rz(1.4872876) q[3];
sx q[3];
rz(-1.828308) q[3];
sx q[3];
rz(0.53572342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80205408) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(1.2667123) q[2];
rz(-0.96261111) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(-0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72120136) q[0];
sx q[0];
rz(-1.2370011) q[0];
sx q[0];
rz(0.23751968) q[0];
rz(-1.0182084) q[1];
sx q[1];
rz(-2.2924481) q[1];
sx q[1];
rz(2.9097897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12122614) q[0];
sx q[0];
rz(-2.0775284) q[0];
sx q[0];
rz(0.10911848) q[0];
x q[1];
rz(-2.899029) q[2];
sx q[2];
rz(-1.3360268) q[2];
sx q[2];
rz(-2.501542) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.17813645) q[1];
sx q[1];
rz(-1.3061211) q[1];
sx q[1];
rz(1.4241649) q[1];
rz(0.2139123) q[3];
sx q[3];
rz(-0.89424664) q[3];
sx q[3];
rz(-0.032534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0314177) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(2.0533662) q[2];
rz(-0.38816372) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(0.80374074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5768455) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(2.172773) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(1.3303403) q[2];
sx q[2];
rz(-1.2546872) q[2];
sx q[2];
rz(1.6457002) q[2];
rz(-0.10235056) q[3];
sx q[3];
rz(-0.1889189) q[3];
sx q[3];
rz(1.2253996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-0.31495467) q[1];
sx q[1];
rz(-2.0575674) q[1];
sx q[1];
rz(-1.6853583) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8946978) q[0];
sx q[0];
rz(-0.1323192) q[0];
sx q[0];
rz(-1.7208862) q[0];
rz(-pi) q[1];
rz(-2.1158754) q[2];
sx q[2];
rz(-1.8946049) q[2];
sx q[2];
rz(-0.88534249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.284277) q[1];
sx q[1];
rz(-0.89968649) q[1];
sx q[1];
rz(-1.3778694) q[1];
rz(-pi) q[2];
rz(-1.625657) q[3];
sx q[3];
rz(-1.8901955) q[3];
sx q[3];
rz(-1.1064135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3866117) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(-0.45271978) q[2];
rz(-0.1581986) q[3];
sx q[3];
rz(-0.69163624) q[3];
sx q[3];
rz(0.89481568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92900705) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(-1.989495) q[0];
rz(-1.2377897) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(0.47098413) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5841056) q[0];
sx q[0];
rz(-2.6563245) q[0];
sx q[0];
rz(-1.4529865) q[0];
rz(-2.1678796) q[2];
sx q[2];
rz(-1.7657585) q[2];
sx q[2];
rz(-0.82414579) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60103154) q[1];
sx q[1];
rz(-0.49113501) q[1];
sx q[1];
rz(-1.4547552) q[1];
rz(-pi) q[2];
rz(-3.0456411) q[3];
sx q[3];
rz(-2.4437685) q[3];
sx q[3];
rz(-0.56688353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0835138) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(-0.2557959) q[2];
rz(-1.4852218) q[3];
sx q[3];
rz(-1.9519613) q[3];
sx q[3];
rz(-0.16168693) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26329041) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(-0.27134744) q[0];
rz(0.73633206) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(-2.7022865) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96734756) q[0];
sx q[0];
rz(-0.085701533) q[0];
sx q[0];
rz(1.218319) q[0];
x q[1];
rz(-0.36053948) q[2];
sx q[2];
rz(-1.6191102) q[2];
sx q[2];
rz(2.609032) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9648793) q[1];
sx q[1];
rz(-2.5039154) q[1];
sx q[1];
rz(-1.8646851) q[1];
rz(2.7078754) q[3];
sx q[3];
rz(-1.0573514) q[3];
sx q[3];
rz(1.0182667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.97418857) q[2];
sx q[2];
rz(-1.5891275) q[2];
sx q[2];
rz(0.20047323) q[2];
rz(-2.386507) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(-0.42373207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077483594) q[0];
sx q[0];
rz(-1.5492726) q[0];
sx q[0];
rz(1.8970998) q[0];
rz(0.81047932) q[1];
sx q[1];
rz(-1.8119718) q[1];
sx q[1];
rz(0.8746075) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32828242) q[0];
sx q[0];
rz(-0.62424849) q[0];
sx q[0];
rz(-1.3024131) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34747296) q[2];
sx q[2];
rz(-0.407019) q[2];
sx q[2];
rz(2.4328872) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7728459) q[1];
sx q[1];
rz(-0.44279848) q[1];
sx q[1];
rz(2.5889791) q[1];
rz(1.6455669) q[3];
sx q[3];
rz(-1.6755591) q[3];
sx q[3];
rz(-1.2391702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0041634) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(-1.7144263) q[2];
rz(-0.066120474) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(1.6920413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604597) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(0.57058913) q[0];
rz(2.5866306) q[1];
sx q[1];
rz(-0.73736063) q[1];
sx q[1];
rz(-0.74329174) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63440454) q[0];
sx q[0];
rz(-0.92538639) q[0];
sx q[0];
rz(2.8651644) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50136106) q[2];
sx q[2];
rz(-1.6622346) q[2];
sx q[2];
rz(0.83133343) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6965461) q[1];
sx q[1];
rz(-1.3475218) q[1];
sx q[1];
rz(0.22534196) q[1];
rz(-pi) q[2];
rz(-0.80446135) q[3];
sx q[3];
rz(-1.7723697) q[3];
sx q[3];
rz(1.7427012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9479998) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(0.26838475) q[2];
rz(2.0949481) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(-0.31782761) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5140117) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(0.50672379) q[0];
rz(0.2535893) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(-2.0862897) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21181606) q[0];
sx q[0];
rz(-1.388071) q[0];
sx q[0];
rz(0.94434785) q[0];
rz(-pi) q[1];
rz(2.6199117) q[2];
sx q[2];
rz(-2.1573967) q[2];
sx q[2];
rz(-0.66643836) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1046022) q[1];
sx q[1];
rz(-1.9676419) q[1];
sx q[1];
rz(-0.87581046) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8860693) q[3];
sx q[3];
rz(-0.75948411) q[3];
sx q[3];
rz(1.7373191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6598597) q[2];
sx q[2];
rz(-2.0969756) q[2];
sx q[2];
rz(-0.8824904) q[2];
rz(0.64583889) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(-1.283949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.5053453) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(2.3690467) q[0];
rz(-1.729471) q[1];
sx q[1];
rz(-1.2071143) q[1];
sx q[1];
rz(0.57377446) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6360146) q[0];
sx q[0];
rz(-1.6806707) q[0];
sx q[0];
rz(1.7626324) q[0];
rz(3.119486) q[2];
sx q[2];
rz(-1.3706285) q[2];
sx q[2];
rz(-0.18891639) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7942012) q[1];
sx q[1];
rz(-0.83354356) q[1];
sx q[1];
rz(2.3458523) q[1];
rz(-pi) q[2];
rz(0.15036924) q[3];
sx q[3];
rz(-1.8024076) q[3];
sx q[3];
rz(-1.3814955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1022169) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(0.77073628) q[2];
rz(-2.7052774) q[3];
sx q[3];
rz(-1.2687012) q[3];
sx q[3];
rz(-1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8687826) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(-1.1707206) q[0];
rz(0.51013485) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(-1.8458813) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2670691) q[0];
sx q[0];
rz(-1.1133615) q[0];
sx q[0];
rz(1.2033071) q[0];
rz(-pi) q[1];
rz(-0.2715271) q[2];
sx q[2];
rz(-1.9121133) q[2];
sx q[2];
rz(1.5835539) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6756729) q[1];
sx q[1];
rz(-0.95341668) q[1];
sx q[1];
rz(-2.8663551) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.027904228) q[3];
sx q[3];
rz(-1.1301665) q[3];
sx q[3];
rz(2.2925216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43508139) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(-2.5734625) q[2];
rz(-2.1614697) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(-0.19395104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4762964) q[0];
sx q[0];
rz(-0.81403533) q[0];
sx q[0];
rz(0.67767674) q[0];
rz(2.9455345) q[1];
sx q[1];
rz(-1.0124857) q[1];
sx q[1];
rz(-0.83818865) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8294551) q[0];
sx q[0];
rz(-1.2078309) q[0];
sx q[0];
rz(1.8118993) q[0];
rz(0.067473472) q[2];
sx q[2];
rz(-2.3258492) q[2];
sx q[2];
rz(-2.2052854) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.14324489) q[1];
sx q[1];
rz(-1.3721826) q[1];
sx q[1];
rz(0.82223383) q[1];
rz(-2.8348654) q[3];
sx q[3];
rz(-0.2704276) q[3];
sx q[3];
rz(-0.21817792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80205408) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(1.8748803) q[2];
rz(-2.1789815) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(-3.0468429) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4203913) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(2.904073) q[0];
rz(-1.0182084) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(0.231803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12122614) q[0];
sx q[0];
rz(-1.0640642) q[0];
sx q[0];
rz(-0.10911848) q[0];
rz(-pi) q[1];
rz(-1.8123774) q[2];
sx q[2];
rz(-1.3350147) q[2];
sx q[2];
rz(-0.98824046) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.787549) q[1];
sx q[1];
rz(-1.7122867) q[1];
sx q[1];
rz(-0.267412) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3123355) q[3];
sx q[3];
rz(-2.4371394) q[3];
sx q[3];
rz(2.7750912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1101749) q[2];
sx q[2];
rz(-1.2538223) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(-0.38816372) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(-2.3378519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5647472) q[0];
sx q[0];
rz(-1.7871465) q[0];
sx q[0];
rz(-0.47252895) q[0];
rz(2.172773) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(0.62933915) q[2];
sx q[2];
rz(-0.39471252) q[2];
sx q[2];
rz(-0.82804745) q[2];
rz(-1.590329) q[3];
sx q[3];
rz(-1.7587147) q[3];
sx q[3];
rz(1.3295909) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

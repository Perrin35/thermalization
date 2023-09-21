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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2468949) q[0];
sx q[0];
rz(-3.0092735) q[0];
sx q[0];
rz(1.7208862) q[0];
rz(-pi) q[1];
rz(-0.37402447) q[2];
sx q[2];
rz(-1.0569388) q[2];
sx q[2];
rz(-2.6467269) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.734182) q[1];
sx q[1];
rz(-1.721518) q[1];
sx q[1];
rz(-2.461344) q[1];
rz(-pi) q[2];
rz(0.31984826) q[3];
sx q[3];
rz(-1.5187129) q[3];
sx q[3];
rz(-0.44714123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.75498092) q[2];
sx q[2];
rz(-1.3957916) q[2];
sx q[2];
rz(-0.45271978) q[2];
rz(0.1581986) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(-2.246777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92900705) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(1.989495) q[0];
rz(1.903803) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(-0.47098413) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7171213) q[0];
sx q[0];
rz(-2.0524128) q[0];
sx q[0];
rz(3.0796914) q[0];
rz(-2.1678796) q[2];
sx q[2];
rz(-1.3758341) q[2];
sx q[2];
rz(0.82414579) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5405611) q[1];
sx q[1];
rz(-0.49113501) q[1];
sx q[1];
rz(-1.4547552) q[1];
rz(-pi) q[2];
rz(0.095951565) q[3];
sx q[3];
rz(-0.6978242) q[3];
sx q[3];
rz(0.56688353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0835138) q[2];
sx q[2];
rz(-0.60516548) q[2];
sx q[2];
rz(0.2557959) q[2];
rz(1.4852218) q[3];
sx q[3];
rz(-1.9519613) q[3];
sx q[3];
rz(0.16168693) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26329041) q[0];
sx q[0];
rz(-2.0354164) q[0];
sx q[0];
rz(0.27134744) q[0];
rz(0.73633206) q[1];
sx q[1];
rz(-1.6059748) q[1];
sx q[1];
rz(2.7022865) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8894316) q[0];
sx q[0];
rz(-1.541242) q[0];
sx q[0];
rz(1.4903402) q[0];
x q[1];
rz(2.7810532) q[2];
sx q[2];
rz(-1.6191102) q[2];
sx q[2];
rz(-0.53256065) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6046391) q[1];
sx q[1];
rz(-2.1770658) q[1];
sx q[1];
rz(2.9301675) q[1];
rz(-pi) q[2];
rz(2.1268232) q[3];
sx q[3];
rz(-1.945567) q[3];
sx q[3];
rz(-2.8127363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.97418857) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(0.20047323) q[2];
rz(2.386507) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(-0.42373207) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0641091) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(-1.8970998) q[0];
rz(-0.81047932) q[1];
sx q[1];
rz(-1.8119718) q[1];
sx q[1];
rz(-0.8746075) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4865206) q[0];
sx q[0];
rz(-0.97210303) q[0];
sx q[0];
rz(-0.18874164) q[0];
x q[1];
rz(1.7165519) q[2];
sx q[2];
rz(-1.1894023) q[2];
sx q[2];
rz(-0.33304735) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9676535) q[1];
sx q[1];
rz(-1.9441009) q[1];
sx q[1];
rz(-1.3268382) q[1];
x q[2];
rz(2.5238958) q[3];
sx q[3];
rz(-3.0129637) q[3];
sx q[3];
rz(-0.61755141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13742927) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(1.7144263) q[2];
rz(-0.066120474) q[3];
sx q[3];
rz(-2.7757006) q[3];
sx q[3];
rz(-1.6920413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-2.5710035) q[0];
rz(-2.5866306) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(2.3983009) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0750908) q[0];
sx q[0];
rz(-0.69426232) q[0];
sx q[0];
rz(-1.2230722) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6749803) q[2];
sx q[2];
rz(-1.0717234) q[2];
sx q[2];
rz(-0.68945976) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35768269) q[1];
sx q[1];
rz(-2.8257113) q[1];
sx q[1];
rz(2.3481579) q[1];
rz(2.8652142) q[3];
sx q[3];
rz(-2.31782) q[3];
sx q[3];
rz(0.36229047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1935929) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(-1.0466446) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(2.823765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.627581) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(-0.50672379) q[0];
rz(-2.8880033) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(1.0553029) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9133638) q[0];
sx q[0];
rz(-0.95634395) q[0];
sx q[0];
rz(-0.22426228) q[0];
x q[1];
rz(0.91674532) q[2];
sx q[2];
rz(-1.142821) q[2];
sx q[2];
rz(-2.5452754) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1736974) q[1];
sx q[1];
rz(-2.3579862) q[1];
sx q[1];
rz(0.99131363) q[1];
rz(-1.33527) q[3];
sx q[3];
rz(-0.84170656) q[3];
sx q[3];
rz(2.083076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6598597) q[2];
sx q[2];
rz(-2.0969756) q[2];
sx q[2];
rz(-2.2591023) q[2];
rz(-2.4957538) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(-1.8576436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362474) q[0];
sx q[0];
rz(-1.212965) q[0];
sx q[0];
rz(2.3690467) q[0];
rz(1.4121217) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(2.5678182) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.505578) q[0];
sx q[0];
rz(-1.4609219) q[0];
sx q[0];
rz(-1.3789603) q[0];
rz(-1.7710118) q[2];
sx q[2];
rz(-1.5924615) q[2];
sx q[2];
rz(1.3862762) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7942012) q[1];
sx q[1];
rz(-2.3080491) q[1];
sx q[1];
rz(-0.79574037) q[1];
x q[2];
rz(-1.8049559) q[3];
sx q[3];
rz(-1.7171211) q[3];
sx q[3];
rz(-2.9175266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1022169) q[2];
sx q[2];
rz(-0.45414671) q[2];
sx q[2];
rz(0.77073628) q[2];
rz(-2.7052774) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(-1.1707206) q[0];
rz(-2.6314578) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(-1.8458813) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8745236) q[0];
sx q[0];
rz(-2.0282312) q[0];
sx q[0];
rz(-1.2033071) q[0];
x q[1];
rz(-2.8700656) q[2];
sx q[2];
rz(-1.2294793) q[2];
sx q[2];
rz(1.5835539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8746652) q[1];
sx q[1];
rz(-1.3473359) q[1];
sx q[1];
rz(-0.93519559) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1300163) q[3];
sx q[3];
rz(-1.545558) q[3];
sx q[3];
rz(2.4079635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7065113) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(2.5734625) q[2];
rz(-0.98012296) q[3];
sx q[3];
rz(-1.0390493) q[3];
sx q[3];
rz(-0.19395104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66529626) q[0];
sx q[0];
rz(-0.81403533) q[0];
sx q[0];
rz(2.4639159) q[0];
rz(-0.19605818) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(0.83818865) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4351589) q[0];
sx q[0];
rz(-2.7088232) q[0];
sx q[0];
rz(-2.58034) q[0];
rz(0.067473472) q[2];
sx q[2];
rz(-2.3258492) q[2];
sx q[2];
rz(0.93630723) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5327685) q[1];
sx q[1];
rz(-0.84034398) q[1];
sx q[1];
rz(2.8735012) q[1];
rz(-pi) q[2];
rz(-1.4872876) q[3];
sx q[3];
rz(-1.3132846) q[3];
sx q[3];
rz(-2.6058692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3395386) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.8748803) q[2];
rz(-0.96261111) q[3];
sx q[3];
rz(-1.5714785) q[3];
sx q[3];
rz(0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.72120136) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(0.23751968) q[0];
rz(1.0182084) q[1];
sx q[1];
rz(-2.2924481) q[1];
sx q[1];
rz(0.231803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12122614) q[0];
sx q[0];
rz(-2.0775284) q[0];
sx q[0];
rz(3.0324742) q[0];
x q[1];
rz(0.24256369) q[2];
sx q[2];
rz(-1.8055658) q[2];
sx q[2];
rz(2.501542) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3540436) q[1];
sx q[1];
rz(-1.429306) q[1];
sx q[1];
rz(0.267412) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3123355) q[3];
sx q[3];
rz(-2.4371394) q[3];
sx q[3];
rz(0.36650141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0314177) q[2];
sx q[2];
rz(-1.2538223) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(-2.7534289) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(-2.3378519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5768455) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(0.96881962) q[1];
sx q[1];
rz(-0.70822721) q[1];
sx q[1];
rz(0.72475564) q[1];
rz(2.5122535) q[2];
sx q[2];
rz(-2.7468801) q[2];
sx q[2];
rz(2.3135452) q[2];
rz(0.10235056) q[3];
sx q[3];
rz(-2.9526738) q[3];
sx q[3];
rz(-1.916193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

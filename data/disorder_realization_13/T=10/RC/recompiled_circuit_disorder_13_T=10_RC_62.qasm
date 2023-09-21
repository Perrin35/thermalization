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
rz(4.6306643) q[0];
sx q[0];
rz(10.319933) q[0];
rz(2.826638) q[1];
sx q[1];
rz(-1.0840253) q[1];
sx q[1];
rz(1.6853583) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8946978) q[0];
sx q[0];
rz(-0.1323192) q[0];
sx q[0];
rz(1.7208862) q[0];
x q[1];
rz(-2.1158754) q[2];
sx q[2];
rz(-1.2469877) q[2];
sx q[2];
rz(0.88534249) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1617042) q[1];
sx q[1];
rz(-2.4474505) q[1];
sx q[1];
rz(2.90467) q[1];
rz(-pi) q[2];
rz(-2.9772894) q[3];
sx q[3];
rz(-0.32391732) q[3];
sx q[3];
rz(1.2795554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75498092) q[2];
sx q[2];
rz(-1.3957916) q[2];
sx q[2];
rz(-2.6888729) q[2];
rz(-2.9833941) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(0.89481568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92900705) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(-1.1520977) q[0];
rz(1.2377897) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(2.6706085) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5841056) q[0];
sx q[0];
rz(-0.4852681) q[0];
sx q[0];
rz(1.4529865) q[0];
x q[1];
rz(-2.9071964) q[2];
sx q[2];
rz(-2.1550551) q[2];
sx q[2];
rz(0.61569475) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.60103154) q[1];
sx q[1];
rz(-2.6504576) q[1];
sx q[1];
rz(-1.6868375) q[1];
rz(-pi) q[2];
rz(-2.4460375) q[3];
sx q[3];
rz(-1.6323946) q[3];
sx q[3];
rz(1.0775281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0835138) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(2.8857968) q[2];
rz(-1.4852218) q[3];
sx q[3];
rz(-1.9519613) q[3];
sx q[3];
rz(-0.16168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26329041) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(0.27134744) q[0];
rz(2.4052606) q[1];
sx q[1];
rz(-1.6059748) q[1];
sx q[1];
rz(-2.7022865) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8205748) q[0];
sx q[0];
rz(-1.4903755) q[0];
sx q[0];
rz(-0.029650173) q[0];
rz(-pi) q[1];
rz(-0.36053948) q[2];
sx q[2];
rz(-1.6191102) q[2];
sx q[2];
rz(2.609032) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5369536) q[1];
sx q[1];
rz(-2.1770658) q[1];
sx q[1];
rz(-0.21142516) q[1];
rz(-0.4337173) q[3];
sx q[3];
rz(-2.0842413) q[3];
sx q[3];
rz(2.123326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1674041) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(0.20047323) q[2];
rz(-0.75508562) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(-2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077483594) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(1.2444929) q[0];
rz(-0.81047932) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(0.8746075) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8133102) q[0];
sx q[0];
rz(-2.5173442) q[0];
sx q[0];
rz(-1.3024131) q[0];
x q[1];
rz(-2.7941197) q[2];
sx q[2];
rz(-0.407019) q[2];
sx q[2];
rz(0.70870542) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.7728459) q[1];
sx q[1];
rz(-0.44279848) q[1];
sx q[1];
rz(-2.5889791) q[1];
rz(2.5238958) q[3];
sx q[3];
rz(-0.12862895) q[3];
sx q[3];
rz(-2.5240412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13742927) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(1.4271663) q[2];
rz(0.066120474) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(-1.6920413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604597) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(0.57058913) q[0];
rz(0.55496201) q[1];
sx q[1];
rz(-0.73736063) q[1];
sx q[1];
rz(-2.3983009) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3742204) q[0];
sx q[0];
rz(-1.351007) q[0];
sx q[0];
rz(0.90669294) q[0];
rz(-pi) q[1];
rz(-0.18851738) q[2];
sx q[2];
rz(-0.50893116) q[2];
sx q[2];
rz(-0.9045507) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.35768269) q[1];
sx q[1];
rz(-2.8257113) q[1];
sx q[1];
rz(2.3481579) q[1];
x q[2];
rz(-2.3371313) q[3];
sx q[3];
rz(-1.7723697) q[3];
sx q[3];
rz(-1.7427012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9479998) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(1.0466446) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5140117) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(-0.50672379) q[0];
rz(-2.8880033) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(2.0862897) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9133638) q[0];
sx q[0];
rz(-0.95634395) q[0];
sx q[0];
rz(-2.9173304) q[0];
rz(2.2141586) q[2];
sx q[2];
rz(-0.76403996) q[2];
sx q[2];
rz(1.470679) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.222059) q[1];
sx q[1];
rz(-0.93898458) q[1];
sx q[1];
rz(-2.6421089) q[1];
rz(-1.8063227) q[3];
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
rz(-0.64583889) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(-1.283949) q[3];
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
rz(-1.6362474) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(0.77254599) q[0];
rz(-1.4121217) q[1];
sx q[1];
rz(-1.2071143) q[1];
sx q[1];
rz(-0.57377446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.7709641) q[2];
sx q[2];
rz(0.18891639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7942012) q[1];
sx q[1];
rz(-2.3080491) q[1];
sx q[1];
rz(-2.3458523) q[1];
x q[2];
rz(1.3366367) q[3];
sx q[3];
rz(-1.7171211) q[3];
sx q[3];
rz(-2.9175266) q[3];
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
rz(-1.8728914) q[3];
sx q[3];
rz(-1.2873945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.27281) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(1.9708721) q[0];
rz(-0.51013485) q[1];
sx q[1];
rz(-1.7850103) q[1];
sx q[1];
rz(1.2957113) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8745236) q[0];
sx q[0];
rz(-2.0282312) q[0];
sx q[0];
rz(1.9382856) q[0];
x q[1];
rz(-2.2175118) q[2];
sx q[2];
rz(-0.43281049) q[2];
sx q[2];
rz(0.88976394) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8746652) q[1];
sx q[1];
rz(-1.3473359) q[1];
sx q[1];
rz(2.2063971) q[1];
x q[2];
rz(1.6298953) q[3];
sx q[3];
rz(-0.44145465) q[3];
sx q[3];
rz(2.3578701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43508139) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(-2.5734625) q[2];
rz(0.98012296) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4762964) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(-0.67767674) q[0];
rz(-0.19605818) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(0.83818865) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7064338) q[0];
sx q[0];
rz(-0.43276946) q[0];
sx q[0];
rz(2.58034) q[0];
x q[1];
rz(2.3269862) q[2];
sx q[2];
rz(-1.6199154) q[2];
sx q[2];
rz(-2.5533822) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2182525) q[1];
sx q[1];
rz(-0.76949161) q[1];
sx q[1];
rz(-1.283265) q[1];
rz(-pi) q[2];
rz(-1.6543051) q[3];
sx q[3];
rz(-1.828308) q[3];
sx q[3];
rz(-2.6058692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3395386) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.2667123) q[2];
rz(-0.96261111) q[3];
sx q[3];
rz(-1.5714785) q[3];
sx q[3];
rz(-3.0468429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4203913) q[0];
sx q[0];
rz(-1.2370011) q[0];
sx q[0];
rz(0.23751968) q[0];
rz(-2.1233842) q[1];
sx q[1];
rz(-2.2924481) q[1];
sx q[1];
rz(-2.9097897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0203665) q[0];
sx q[0];
rz(-2.0775284) q[0];
sx q[0];
rz(3.0324742) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3292153) q[2];
sx q[2];
rz(-1.806578) q[2];
sx q[2];
rz(2.1533522) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.17813645) q[1];
sx q[1];
rz(-1.8354715) q[1];
sx q[1];
rz(1.4241649) q[1];
rz(-pi) q[2];
rz(0.2139123) q[3];
sx q[3];
rz(-0.89424664) q[3];
sx q[3];
rz(-0.032534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0314177) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(2.7534289) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(2.3378519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.5768455) q[0];
sx q[0];
rz(-1.7871465) q[0];
sx q[0];
rz(-0.47252895) q[0];
rz(2.172773) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(2.8167562) q[2];
sx q[2];
rz(-1.799121) q[2];
sx q[2];
rz(0.15098235) q[2];
rz(3.0392421) q[3];
sx q[3];
rz(-0.1889189) q[3];
sx q[3];
rz(1.2253996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
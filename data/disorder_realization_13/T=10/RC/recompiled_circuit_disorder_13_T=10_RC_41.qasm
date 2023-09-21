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
rz(2.0575674) q[1];
sx q[1];
rz(10.881012) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0460912) q[0];
sx q[0];
rz(-1.7016194) q[0];
sx q[0];
rz(3.1216937) q[0];
x q[1];
rz(-0.37402447) q[2];
sx q[2];
rz(-1.0569388) q[2];
sx q[2];
rz(0.49486578) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40741062) q[1];
sx q[1];
rz(-1.721518) q[1];
sx q[1];
rz(2.461344) q[1];
x q[2];
rz(2.8217444) q[3];
sx q[3];
rz(-1.6228798) q[3];
sx q[3];
rz(-0.44714123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3866117) q[2];
sx q[2];
rz(-1.3957916) q[2];
sx q[2];
rz(0.45271978) q[2];
rz(0.1581986) q[3];
sx q[3];
rz(-0.69163624) q[3];
sx q[3];
rz(2.246777) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92900705) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(1.989495) q[0];
rz(-1.903803) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(-2.6706085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1176227) q[0];
sx q[0];
rz(-1.6256486) q[0];
sx q[0];
rz(1.0883925) q[0];
rz(-pi) q[1];
rz(2.9071964) q[2];
sx q[2];
rz(-0.98653754) q[2];
sx q[2];
rz(-2.5258979) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4091332) q[1];
sx q[1];
rz(-2.0583378) q[1];
sx q[1];
rz(0.061846102) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4906293) q[3];
sx q[3];
rz(-0.87682322) q[3];
sx q[3];
rz(0.44192867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.058078893) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8783022) q[0];
sx q[0];
rz(-2.0354164) q[0];
sx q[0];
rz(-0.27134744) q[0];
rz(2.4052606) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(-0.43930611) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1742451) q[0];
sx q[0];
rz(-3.0558911) q[0];
sx q[0];
rz(-1.9232737) q[0];
rz(-pi) q[1];
rz(-2.7810532) q[2];
sx q[2];
rz(-1.6191102) q[2];
sx q[2];
rz(-2.609032) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9648793) q[1];
sx q[1];
rz(-0.63767725) q[1];
sx q[1];
rz(-1.8646851) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1268232) q[3];
sx q[3];
rz(-1.1960256) q[3];
sx q[3];
rz(-0.32885636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1674041) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(2.9411194) q[2];
rz(-2.386507) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(-2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0641091) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(1.8970998) q[0];
rz(2.3311133) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(0.8746075) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32828242) q[0];
sx q[0];
rz(-2.5173442) q[0];
sx q[0];
rz(1.3024131) q[0];
x q[1];
rz(-0.38509102) q[2];
sx q[2];
rz(-1.706012) q[2];
sx q[2];
rz(1.9584292) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7728459) q[1];
sx q[1];
rz(-0.44279848) q[1];
sx q[1];
rz(0.55261353) q[1];
rz(-pi) q[2];
rz(-0.10505418) q[3];
sx q[3];
rz(-1.6451562) q[3];
sx q[3];
rz(0.33945938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0041634) q[2];
sx q[2];
rz(-0.88976294) q[2];
sx q[2];
rz(1.4271663) q[2];
rz(-0.066120474) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(-1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78113294) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(0.57058913) q[0];
rz(2.5866306) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(-2.3983009) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5071881) q[0];
sx q[0];
rz(-0.92538639) q[0];
sx q[0];
rz(0.27642823) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18851738) q[2];
sx q[2];
rz(-2.6326615) q[2];
sx q[2];
rz(0.9045507) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.35768269) q[1];
sx q[1];
rz(-2.8257113) q[1];
sx q[1];
rz(-0.79343474) q[1];
rz(-pi) q[2];
rz(0.27637847) q[3];
sx q[3];
rz(-2.31782) q[3];
sx q[3];
rz(-0.36229047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9479998) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(0.26838475) q[2];
rz(1.0466446) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-1.5140117) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(-2.6348689) q[0];
rz(-2.8880033) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(-2.0862897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5366093) q[0];
sx q[0];
rz(-0.6491001) q[0];
sx q[0];
rz(-1.8761294) q[0];
rz(-pi) q[1];
rz(0.52168092) q[2];
sx q[2];
rz(-2.1573967) q[2];
sx q[2];
rz(-2.4751543) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.222059) q[1];
sx q[1];
rz(-2.2026081) q[1];
sx q[1];
rz(-2.6421089) q[1];
rz(-1.33527) q[3];
sx q[3];
rz(-2.2998861) q[3];
sx q[3];
rz(1.0585166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.48173299) q[2];
sx q[2];
rz(-2.0969756) q[2];
sx q[2];
rz(-2.2591023) q[2];
rz(-0.64583889) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(1.8576436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362474) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(-2.3690467) q[0];
rz(1.4121217) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(2.5678182) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6360146) q[0];
sx q[0];
rz(-1.4609219) q[0];
sx q[0];
rz(-1.7626324) q[0];
rz(-pi) q[1];
rz(1.4622719) q[2];
sx q[2];
rz(-2.9402241) q[2];
sx q[2];
rz(-0.078171922) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3473914) q[1];
sx q[1];
rz(-2.3080491) q[1];
sx q[1];
rz(2.3458523) q[1];
rz(-pi) q[2];
rz(-0.15036924) q[3];
sx q[3];
rz(-1.339185) q[3];
sx q[3];
rz(1.7600972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(0.77073628) q[2];
rz(0.43631521) q[3];
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
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8687826) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(1.1707206) q[0];
rz(0.51013485) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(1.2957113) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2670691) q[0];
sx q[0];
rz(-1.1133615) q[0];
sx q[0];
rz(-1.2033071) q[0];
rz(2.8700656) q[2];
sx q[2];
rz(-1.2294793) q[2];
sx q[2];
rz(1.5580387) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1295373) q[1];
sx q[1];
rz(-2.4730198) q[1];
sx q[1];
rz(-1.936391) q[1];
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
rz(pi/2) q[1];
rz(-2.7065113) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(0.56813017) q[2];
rz(-2.1614697) q[3];
sx q[3];
rz(-1.0390493) q[3];
sx q[3];
rz(-2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4351589) q[0];
sx q[0];
rz(-2.7088232) q[0];
sx q[0];
rz(-2.58034) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3269862) q[2];
sx q[2];
rz(-1.5216773) q[2];
sx q[2];
rz(-2.5533822) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9233401) q[1];
sx q[1];
rz(-2.372101) q[1];
sx q[1];
rz(-1.283265) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8348654) q[3];
sx q[3];
rz(-2.8711651) q[3];
sx q[3];
rz(0.21817792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3395386) q[2];
sx q[2];
rz(-2.4513117) q[2];
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
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72120136) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(2.904073) q[0];
rz(-2.1233842) q[1];
sx q[1];
rz(-2.2924481) q[1];
sx q[1];
rz(0.231803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7451413) q[0];
sx q[0];
rz(-1.475435) q[0];
sx q[0];
rz(-1.061529) q[0];
rz(-0.78328697) q[2];
sx q[2];
rz(-2.8056393) q[2];
sx q[2];
rz(2.96539) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4494891) q[1];
sx q[1];
rz(-2.8398501) q[1];
sx q[1];
rz(-0.49441378) q[1];
rz(-pi) q[2];
rz(-0.2139123) q[3];
sx q[3];
rz(-0.89424664) q[3];
sx q[3];
rz(0.032534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0314177) q[2];
sx q[2];
rz(-1.2538223) q[2];
sx q[2];
rz(2.0533662) q[2];
rz(-2.7534289) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(0.80374074) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5647472) q[0];
sx q[0];
rz(-1.7871465) q[0];
sx q[0];
rz(-0.47252895) q[0];
rz(-0.96881962) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(1.3303403) q[2];
sx q[2];
rz(-1.2546872) q[2];
sx q[2];
rz(1.6457002) q[2];
rz(1.5512636) q[3];
sx q[3];
rz(-1.7587147) q[3];
sx q[3];
rz(1.3295909) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
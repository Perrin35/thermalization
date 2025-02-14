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
rz(2.8259377) q[0];
rz(1.1701801) q[1];
sx q[1];
rz(-0.38771954) q[1];
sx q[1];
rz(-2.5375836) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9785288) q[0];
sx q[0];
rz(-1.2685568) q[0];
sx q[0];
rz(-3.0885484) q[0];
rz(2.2411795) q[2];
sx q[2];
rz(-1.1657823) q[2];
sx q[2];
rz(2.9069898) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4166222) q[1];
sx q[1];
rz(-1.9068826) q[1];
sx q[1];
rz(1.027847) q[1];
rz(-0.1084566) q[3];
sx q[3];
rz(-1.2634957) q[3];
sx q[3];
rz(1.1600329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.151256) q[2];
sx q[2];
rz(-1.9467111) q[2];
sx q[2];
rz(1.3585496) q[2];
rz(1.547706) q[3];
sx q[3];
rz(-0.93021506) q[3];
sx q[3];
rz(1.5096629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58594054) q[0];
sx q[0];
rz(-0.63783115) q[0];
sx q[0];
rz(1.1974539) q[0];
rz(-0.50152957) q[1];
sx q[1];
rz(-1.5781559) q[1];
sx q[1];
rz(1.7053568) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7345463) q[0];
sx q[0];
rz(-1.0031896) q[0];
sx q[0];
rz(2.8664886) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11709638) q[2];
sx q[2];
rz(-0.70636049) q[2];
sx q[2];
rz(-1.381284) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8433709) q[1];
sx q[1];
rz(-2.9195115) q[1];
sx q[1];
rz(2.6574617) q[1];
rz(-pi) q[2];
rz(2.7327939) q[3];
sx q[3];
rz(-0.49064562) q[3];
sx q[3];
rz(-1.4178432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.88840914) q[2];
sx q[2];
rz(-1.2359572) q[2];
sx q[2];
rz(3.1186228) q[2];
rz(-2.2394771) q[3];
sx q[3];
rz(-1.918856) q[3];
sx q[3];
rz(-0.42660108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4334634) q[0];
sx q[0];
rz(-1.3490278) q[0];
sx q[0];
rz(0.9915114) q[0];
rz(2.1082711) q[1];
sx q[1];
rz(-0.28305498) q[1];
sx q[1];
rz(-1.906377) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18444321) q[0];
sx q[0];
rz(-2.1547124) q[0];
sx q[0];
rz(-1.2725485) q[0];
rz(2.5052983) q[2];
sx q[2];
rz(-2.4264196) q[2];
sx q[2];
rz(-1.6805122) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24656235) q[1];
sx q[1];
rz(-2.3091279) q[1];
sx q[1];
rz(-0.33729302) q[1];
rz(-2.1622873) q[3];
sx q[3];
rz(-2.0404173) q[3];
sx q[3];
rz(-2.8036433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.027355) q[2];
sx q[2];
rz(-0.75216746) q[2];
sx q[2];
rz(-2.1503964) q[2];
rz(-1.8441955) q[3];
sx q[3];
rz(-1.2605366) q[3];
sx q[3];
rz(-1.7245002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39795136) q[0];
sx q[0];
rz(-2.209111) q[0];
sx q[0];
rz(-0.81746307) q[0];
rz(0.3262597) q[1];
sx q[1];
rz(-1.9605109) q[1];
sx q[1];
rz(-0.78525966) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2366339) q[0];
sx q[0];
rz(-1.8371305) q[0];
sx q[0];
rz(1.2575218) q[0];
rz(0.26008028) q[2];
sx q[2];
rz(-0.65979119) q[2];
sx q[2];
rz(-1.5957956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2438812) q[1];
sx q[1];
rz(-0.63279018) q[1];
sx q[1];
rz(0.80676078) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4912075) q[3];
sx q[3];
rz(-1.2168988) q[3];
sx q[3];
rz(-0.45512629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0205445) q[2];
sx q[2];
rz(-0.64457568) q[2];
sx q[2];
rz(-2.5784967) q[2];
rz(-2.2480887) q[3];
sx q[3];
rz(-1.9681135) q[3];
sx q[3];
rz(1.2347429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-2.0183759) q[1];
sx q[1];
rz(0.18128577) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3485506) q[0];
sx q[0];
rz(-0.53213813) q[0];
sx q[0];
rz(-1.9572958) q[0];
rz(1.6185804) q[2];
sx q[2];
rz(-1.0956229) q[2];
sx q[2];
rz(-0.20296861) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5904894) q[1];
sx q[1];
rz(-0.22399513) q[1];
sx q[1];
rz(1.6975387) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44293176) q[3];
sx q[3];
rz(-1.7128403) q[3];
sx q[3];
rz(1.0012817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.01866092) q[2];
sx q[2];
rz(-0.67492008) q[2];
sx q[2];
rz(1.8956511) q[2];
rz(-2.5490226) q[3];
sx q[3];
rz(-1.4453245) q[3];
sx q[3];
rz(-2.3004801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6084006) q[0];
sx q[0];
rz(-2.1493981) q[0];
sx q[0];
rz(-0.30884185) q[0];
rz(0.32288512) q[1];
sx q[1];
rz(-2.7643118) q[1];
sx q[1];
rz(-3.0016532) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4734264) q[0];
sx q[0];
rz(-0.65957171) q[0];
sx q[0];
rz(-1.1585328) q[0];
rz(-pi) q[1];
rz(1.3792737) q[2];
sx q[2];
rz(-0.46326783) q[2];
sx q[2];
rz(-2.6216216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7227712) q[1];
sx q[1];
rz(-0.79192894) q[1];
sx q[1];
rz(-2.0008816) q[1];
rz(-2.841137) q[3];
sx q[3];
rz(-2.2078035) q[3];
sx q[3];
rz(-0.14618044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4796925) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(2.8996186) q[2];
rz(0.74609977) q[3];
sx q[3];
rz(-1.7753764) q[3];
sx q[3];
rz(1.7297176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.022843) q[0];
sx q[0];
rz(-2.0976837) q[0];
sx q[0];
rz(1.8705077) q[0];
rz(-2.5785043) q[1];
sx q[1];
rz(-2.2124898) q[1];
sx q[1];
rz(-1.44106) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.010276) q[0];
sx q[0];
rz(-2.1621341) q[0];
sx q[0];
rz(-2.813126) q[0];
x q[1];
rz(-2.8071324) q[2];
sx q[2];
rz(-0.40142599) q[2];
sx q[2];
rz(-0.73094598) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.064441) q[1];
sx q[1];
rz(-1.3485326) q[1];
sx q[1];
rz(-0.42776107) q[1];
rz(-2.8981179) q[3];
sx q[3];
rz(-0.7400652) q[3];
sx q[3];
rz(-0.26672637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.21961221) q[2];
sx q[2];
rz(-1.8461123) q[2];
sx q[2];
rz(-1.8249576) q[2];
rz(2.5804139) q[3];
sx q[3];
rz(-0.96446529) q[3];
sx q[3];
rz(1.6130028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0372666) q[0];
sx q[0];
rz(-2.9236561) q[0];
sx q[0];
rz(-0.88687801) q[0];
rz(0.2001702) q[1];
sx q[1];
rz(-1.6693516) q[1];
sx q[1];
rz(-1.0221457) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048269317) q[0];
sx q[0];
rz(-1.5962068) q[0];
sx q[0];
rz(1.7019468) q[0];
rz(2.3273507) q[2];
sx q[2];
rz(-0.98746429) q[2];
sx q[2];
rz(0.63177201) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8238942) q[1];
sx q[1];
rz(-0.824747) q[1];
sx q[1];
rz(-3.1391678) q[1];
rz(0.57557801) q[3];
sx q[3];
rz(-1.8691571) q[3];
sx q[3];
rz(-0.23250599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48978051) q[2];
sx q[2];
rz(-2.688789) q[2];
sx q[2];
rz(-2.32302) q[2];
rz(2.1583648) q[3];
sx q[3];
rz(-2.1849617) q[3];
sx q[3];
rz(2.6035068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.035456903) q[0];
sx q[0];
rz(-1.1470733) q[0];
sx q[0];
rz(-1.5552833) q[0];
rz(-2.4447794) q[1];
sx q[1];
rz(-2.9882444) q[1];
sx q[1];
rz(-2.3568025) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9826384) q[0];
sx q[0];
rz(-2.1110015) q[0];
sx q[0];
rz(-1.0544974) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3680787) q[2];
sx q[2];
rz(-1.4215111) q[2];
sx q[2];
rz(1.210807) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1363298) q[1];
sx q[1];
rz(-1.5371833) q[1];
sx q[1];
rz(2.7064347) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5303477) q[3];
sx q[3];
rz(-1.9211624) q[3];
sx q[3];
rz(0.65920748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4697504) q[2];
sx q[2];
rz(-1.2776813) q[2];
sx q[2];
rz(-2.3308241) q[2];
rz(1.9226711) q[3];
sx q[3];
rz(-0.54088497) q[3];
sx q[3];
rz(2.0764652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78275457) q[0];
sx q[0];
rz(-2.8420119) q[0];
sx q[0];
rz(1.7175571) q[0];
rz(0.77761039) q[1];
sx q[1];
rz(-0.97336665) q[1];
sx q[1];
rz(2.3861859) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7292405) q[0];
sx q[0];
rz(-1.5662114) q[0];
sx q[0];
rz(-1.7807177) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.099775) q[2];
sx q[2];
rz(-2.7448069) q[2];
sx q[2];
rz(0.40601847) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3302119) q[1];
sx q[1];
rz(-2.5054058) q[1];
sx q[1];
rz(-0.61417555) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4660268) q[3];
sx q[3];
rz(-1.3119952) q[3];
sx q[3];
rz(-1.344556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5648254) q[2];
sx q[2];
rz(-0.29594031) q[2];
sx q[2];
rz(0.15288615) q[2];
rz(-1.2711924) q[3];
sx q[3];
rz(-1.252424) q[3];
sx q[3];
rz(2.474474) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2818903) q[0];
sx q[0];
rz(-0.49260456) q[0];
sx q[0];
rz(2.3717666) q[0];
rz(-1.2234756) q[1];
sx q[1];
rz(-1.2212831) q[1];
sx q[1];
rz(2.4881359) q[1];
rz(-0.59719795) q[2];
sx q[2];
rz(-1.2033249) q[2];
sx q[2];
rz(-2.0133599) q[2];
rz(-0.85930227) q[3];
sx q[3];
rz(-2.6283384) q[3];
sx q[3];
rz(2.8859861) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

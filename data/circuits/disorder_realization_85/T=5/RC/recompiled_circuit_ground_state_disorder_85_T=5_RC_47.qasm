OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.10211927) q[0];
sx q[0];
rz(4.6146225) q[0];
sx q[0];
rz(6.1518402) q[0];
rz(-3.1047473) q[1];
sx q[1];
rz(3.6689833) q[1];
sx q[1];
rz(13.875516) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9835684) q[0];
sx q[0];
rz(-0.82497549) q[0];
sx q[0];
rz(1.0452861) q[0];
x q[1];
rz(2.1151343) q[2];
sx q[2];
rz(-1.9690459) q[2];
sx q[2];
rz(1.2521225) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4507452) q[1];
sx q[1];
rz(-1.2589129) q[1];
sx q[1];
rz(-1.8212832) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9580826) q[3];
sx q[3];
rz(-1.402463) q[3];
sx q[3];
rz(0.77372293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5347791) q[2];
sx q[2];
rz(-2.7853192) q[2];
sx q[2];
rz(0.66205364) q[2];
rz(-0.74696294) q[3];
sx q[3];
rz(-2.2253939) q[3];
sx q[3];
rz(2.1257187) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5544283) q[0];
sx q[0];
rz(-0.14475188) q[0];
sx q[0];
rz(1.9594877) q[0];
rz(-2.3960522) q[1];
sx q[1];
rz(-2.5423971) q[1];
sx q[1];
rz(0.10339698) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.11907) q[0];
sx q[0];
rz(-1.3387209) q[0];
sx q[0];
rz(-2.705535) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.75704) q[2];
sx q[2];
rz(-2.1332624) q[2];
sx q[2];
rz(-0.99316521) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1704639) q[1];
sx q[1];
rz(-2.269575) q[1];
sx q[1];
rz(-2.7034195) q[1];
x q[2];
rz(-2.313226) q[3];
sx q[3];
rz(-1.6992555) q[3];
sx q[3];
rz(1.0297694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.55249247) q[2];
sx q[2];
rz(-1.8742259) q[2];
sx q[2];
rz(-2.6972771) q[2];
rz(-2.0878504) q[3];
sx q[3];
rz(-0.070662347) q[3];
sx q[3];
rz(1.1963371) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.062773) q[0];
sx q[0];
rz(-1.8906931) q[0];
sx q[0];
rz(2.479082) q[0];
rz(-1.6569116) q[1];
sx q[1];
rz(-2.1187481) q[1];
sx q[1];
rz(2.9248765) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0652567) q[0];
sx q[0];
rz(-2.2562593) q[0];
sx q[0];
rz(-0.61409593) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48061252) q[2];
sx q[2];
rz(-2.3490902) q[2];
sx q[2];
rz(0.30016826) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.413447) q[1];
sx q[1];
rz(-1.6108509) q[1];
sx q[1];
rz(-2.2437281) q[1];
x q[2];
rz(-0.49461282) q[3];
sx q[3];
rz(-1.429116) q[3];
sx q[3];
rz(-2.6050267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66389877) q[2];
sx q[2];
rz(-2.615216) q[2];
sx q[2];
rz(0.23507512) q[2];
rz(2.7663686) q[3];
sx q[3];
rz(-2.2347361) q[3];
sx q[3];
rz(2.4231329) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5594056) q[0];
sx q[0];
rz(-0.087787293) q[0];
sx q[0];
rz(-2.1413595) q[0];
rz(1.2031215) q[1];
sx q[1];
rz(-1.6327881) q[1];
sx q[1];
rz(-0.54471725) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2208015) q[0];
sx q[0];
rz(-1.615633) q[0];
sx q[0];
rz(-2.4283571) q[0];
rz(-pi) q[1];
rz(0.80766957) q[2];
sx q[2];
rz(-2.2414506) q[2];
sx q[2];
rz(-2.9602221) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1511606) q[1];
sx q[1];
rz(-0.22391388) q[1];
sx q[1];
rz(-2.1316705) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29981837) q[3];
sx q[3];
rz(-2.162419) q[3];
sx q[3];
rz(-0.064379582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7228221) q[2];
sx q[2];
rz(-2.3493769) q[2];
sx q[2];
rz(-0.79696068) q[2];
rz(0.289251) q[3];
sx q[3];
rz(-2.1713493) q[3];
sx q[3];
rz(-0.42118916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5274984) q[0];
sx q[0];
rz(-0.088936381) q[0];
sx q[0];
rz(-1.8726789) q[0];
rz(0.96619636) q[1];
sx q[1];
rz(-1.393001) q[1];
sx q[1];
rz(0.46636811) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9556676) q[0];
sx q[0];
rz(-2.4644682) q[0];
sx q[0];
rz(2.4311275) q[0];
rz(-1.2205458) q[2];
sx q[2];
rz(-2.008247) q[2];
sx q[2];
rz(-1.3910363) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.98919096) q[1];
sx q[1];
rz(-2.4394803) q[1];
sx q[1];
rz(-2.0753808) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0859231) q[3];
sx q[3];
rz(-0.72332649) q[3];
sx q[3];
rz(1.5679899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7352778) q[2];
sx q[2];
rz(-2.3390528) q[2];
sx q[2];
rz(0.80424133) q[2];
rz(-2.6546226) q[3];
sx q[3];
rz(-2.099497) q[3];
sx q[3];
rz(0.0074726661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.93953472) q[0];
sx q[0];
rz(-2.845293) q[0];
sx q[0];
rz(0.95788389) q[0];
rz(1.6288039) q[1];
sx q[1];
rz(-2.084338) q[1];
sx q[1];
rz(1.5527976) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88547546) q[0];
sx q[0];
rz(-1.6506356) q[0];
sx q[0];
rz(-0.022932963) q[0];
rz(-1.2852766) q[2];
sx q[2];
rz(-1.8381882) q[2];
sx q[2];
rz(-1.3533139) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9839638) q[1];
sx q[1];
rz(-2.2215054) q[1];
sx q[1];
rz(1.7939066) q[1];
rz(2.6090066) q[3];
sx q[3];
rz(-0.69952337) q[3];
sx q[3];
rz(-0.12040779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2019041) q[2];
sx q[2];
rz(-2.0857911) q[2];
sx q[2];
rz(-0.73516694) q[2];
rz(2.1926664) q[3];
sx q[3];
rz(-2.2831423) q[3];
sx q[3];
rz(-0.86400509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3801124) q[0];
sx q[0];
rz(-2.1367456) q[0];
sx q[0];
rz(-2.5349706) q[0];
rz(0.78423777) q[1];
sx q[1];
rz(-1.8083068) q[1];
sx q[1];
rz(-1.3444791) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30596581) q[0];
sx q[0];
rz(-0.84442645) q[0];
sx q[0];
rz(0.65610049) q[0];
x q[1];
rz(1.3188348) q[2];
sx q[2];
rz(-2.0540385) q[2];
sx q[2];
rz(1.2298707) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91025464) q[1];
sx q[1];
rz(-1.3489854) q[1];
sx q[1];
rz(2.7342058) q[1];
rz(-pi) q[2];
rz(-1.7502039) q[3];
sx q[3];
rz(-2.6699237) q[3];
sx q[3];
rz(-2.6149132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.532423) q[2];
sx q[2];
rz(-0.50749856) q[2];
sx q[2];
rz(1.3896821) q[2];
rz(2.9621647) q[3];
sx q[3];
rz(-0.98589412) q[3];
sx q[3];
rz(-2.7348886) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0824025) q[0];
sx q[0];
rz(-0.32408369) q[0];
sx q[0];
rz(-0.75505906) q[0];
rz(-0.45626196) q[1];
sx q[1];
rz(-1.4249964) q[1];
sx q[1];
rz(-1.0770575) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2307325) q[0];
sx q[0];
rz(-0.60711475) q[0];
sx q[0];
rz(-2.7348628) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78387129) q[2];
sx q[2];
rz(-1.4848564) q[2];
sx q[2];
rz(0.29168561) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7491319) q[1];
sx q[1];
rz(-1.8549671) q[1];
sx q[1];
rz(0.47537132) q[1];
rz(1.0302587) q[3];
sx q[3];
rz(-1.7071299) q[3];
sx q[3];
rz(0.34878525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4130752) q[2];
sx q[2];
rz(-1.1001526) q[2];
sx q[2];
rz(0.73823482) q[2];
rz(-1.5089367) q[3];
sx q[3];
rz(-1.3239219) q[3];
sx q[3];
rz(-1.1608605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.64410925) q[0];
sx q[0];
rz(-1.6832385) q[0];
sx q[0];
rz(-2.4926376) q[0];
rz(0.68967462) q[1];
sx q[1];
rz(-0.70976218) q[1];
sx q[1];
rz(0.41742691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6605292) q[0];
sx q[0];
rz(-1.821035) q[0];
sx q[0];
rz(-0.60926837) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6471049) q[2];
sx q[2];
rz(-1.2718437) q[2];
sx q[2];
rz(1.1079009) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3950353) q[1];
sx q[1];
rz(-0.42320028) q[1];
sx q[1];
rz(2.8854135) q[1];
x q[2];
rz(-0.77885742) q[3];
sx q[3];
rz(-0.53611272) q[3];
sx q[3];
rz(-2.6440309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4735585) q[2];
sx q[2];
rz(-1.4213976) q[2];
sx q[2];
rz(3.0958946) q[2];
rz(-2.8081196) q[3];
sx q[3];
rz(-0.57531753) q[3];
sx q[3];
rz(0.12588178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4561653) q[0];
sx q[0];
rz(-2.6021155) q[0];
sx q[0];
rz(1.217655) q[0];
rz(0.24670163) q[1];
sx q[1];
rz(-1.8496937) q[1];
sx q[1];
rz(-2.6620679) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7643023) q[0];
sx q[0];
rz(-2.3545697) q[0];
sx q[0];
rz(0.93955173) q[0];
x q[1];
rz(-2.6423655) q[2];
sx q[2];
rz(-0.93846924) q[2];
sx q[2];
rz(-1.3030753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61726924) q[1];
sx q[1];
rz(-1.8620544) q[1];
sx q[1];
rz(-2.6074383) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44879313) q[3];
sx q[3];
rz(-2.6330607) q[3];
sx q[3];
rz(-2.3859346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5852927) q[2];
sx q[2];
rz(-1.2703398) q[2];
sx q[2];
rz(1.7362107) q[2];
rz(0.65226883) q[3];
sx q[3];
rz(-2.8758719) q[3];
sx q[3];
rz(0.71776596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9115059) q[0];
sx q[0];
rz(-1.390504) q[0];
sx q[0];
rz(2.6934296) q[0];
rz(0.74408342) q[1];
sx q[1];
rz(-2.4241445) q[1];
sx q[1];
rz(-1.4345899) q[1];
rz(2.9058331) q[2];
sx q[2];
rz(-0.59638646) q[2];
sx q[2];
rz(-1.4183945) q[2];
rz(1.0539609) q[3];
sx q[3];
rz(-0.64095961) q[3];
sx q[3];
rz(-0.62957233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

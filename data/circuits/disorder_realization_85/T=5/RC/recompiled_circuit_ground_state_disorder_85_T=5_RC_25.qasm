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
rz(0.036845358) q[1];
sx q[1];
rz(-0.52739066) q[1];
sx q[1];
rz(-1.3091458) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9835684) q[0];
sx q[0];
rz(-2.3166172) q[0];
sx q[0];
rz(1.0452861) q[0];
x q[1];
rz(2.6845161) q[2];
sx q[2];
rz(-1.0731878) q[2];
sx q[2];
rz(0.08809419) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6908474) q[1];
sx q[1];
rz(-1.8826797) q[1];
sx q[1];
rz(-1.3203095) q[1];
rz(-2.3917624) q[3];
sx q[3];
rz(-0.24838003) q[3];
sx q[3];
rz(-1.6102143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5347791) q[2];
sx q[2];
rz(-0.35627347) q[2];
sx q[2];
rz(-2.479539) q[2];
rz(0.74696294) q[3];
sx q[3];
rz(-0.91619879) q[3];
sx q[3];
rz(-1.0158739) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5544283) q[0];
sx q[0];
rz(-2.9968408) q[0];
sx q[0];
rz(-1.1821049) q[0];
rz(-0.74554044) q[1];
sx q[1];
rz(-0.59919557) q[1];
sx q[1];
rz(0.10339698) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7000843) q[0];
sx q[0];
rz(-1.9943976) q[0];
sx q[0];
rz(-1.3157428) q[0];
rz(-2.75704) q[2];
sx q[2];
rz(-2.1332624) q[2];
sx q[2];
rz(-0.99316521) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89241806) q[1];
sx q[1];
rz(-1.9016445) q[1];
sx q[1];
rz(-0.82280092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17354266) q[3];
sx q[3];
rz(-2.3056917) q[3];
sx q[3];
rz(-0.65803448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.55249247) q[2];
sx q[2];
rz(-1.2673667) q[2];
sx q[2];
rz(2.6972771) q[2];
rz(1.0537423) q[3];
sx q[3];
rz(-0.070662347) q[3];
sx q[3];
rz(1.1963371) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.062773) q[0];
sx q[0];
rz(-1.8906931) q[0];
sx q[0];
rz(-2.479082) q[0];
rz(1.484681) q[1];
sx q[1];
rz(-2.1187481) q[1];
sx q[1];
rz(2.9248765) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.076336) q[0];
sx q[0];
rz(-0.88533339) q[0];
sx q[0];
rz(2.5274967) q[0];
rz(-pi) q[1];
rz(2.4091085) q[2];
sx q[2];
rz(-1.9062796) q[2];
sx q[2];
rz(1.5200638) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3308476) q[1];
sx q[1];
rz(-0.89850366) q[1];
sx q[1];
rz(0.051203392) q[1];
rz(-pi) q[2];
rz(-2.8496938) q[3];
sx q[3];
rz(-0.51288285) q[3];
sx q[3];
rz(-1.8512902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66389877) q[2];
sx q[2];
rz(-0.52637664) q[2];
sx q[2];
rz(2.9065175) q[2];
rz(-2.7663686) q[3];
sx q[3];
rz(-0.90685654) q[3];
sx q[3];
rz(2.4231329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5821871) q[0];
sx q[0];
rz(-0.087787293) q[0];
sx q[0];
rz(-2.1413595) q[0];
rz(-1.2031215) q[1];
sx q[1];
rz(-1.5088046) q[1];
sx q[1];
rz(2.5968754) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7398374) q[0];
sx q[0];
rz(-0.71439636) q[0];
sx q[0];
rz(-0.06846662) q[0];
rz(-2.3096309) q[2];
sx q[2];
rz(-2.1429981) q[2];
sx q[2];
rz(-1.9257279) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.96918584) q[1];
sx q[1];
rz(-1.6891857) q[1];
sx q[1];
rz(1.7612996) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29981837) q[3];
sx q[3];
rz(-0.97917367) q[3];
sx q[3];
rz(-3.0772131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7228221) q[2];
sx q[2];
rz(-0.79221574) q[2];
sx q[2];
rz(-0.79696068) q[2];
rz(-0.289251) q[3];
sx q[3];
rz(-2.1713493) q[3];
sx q[3];
rz(-2.7204035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5274984) q[0];
sx q[0];
rz(-0.088936381) q[0];
sx q[0];
rz(1.8726789) q[0];
rz(-2.1753963) q[1];
sx q[1];
rz(-1.393001) q[1];
sx q[1];
rz(0.46636811) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9357929) q[0];
sx q[0];
rz(-1.1498435) q[0];
sx q[0];
rz(2.5942786) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63299243) q[2];
sx q[2];
rz(-2.5883753) q[2];
sx q[2];
rz(-2.4621682) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6154229) q[1];
sx q[1];
rz(-0.96994441) q[1];
sx q[1];
rz(-0.38819617) q[1];
x q[2];
rz(-1.0859231) q[3];
sx q[3];
rz(-2.4182662) q[3];
sx q[3];
rz(1.5679899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7352778) q[2];
sx q[2];
rz(-0.80253989) q[2];
sx q[2];
rz(0.80424133) q[2];
rz(0.4869701) q[3];
sx q[3];
rz(-2.099497) q[3];
sx q[3];
rz(0.0074726661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2020579) q[0];
sx q[0];
rz(-2.845293) q[0];
sx q[0];
rz(-2.1837088) q[0];
rz(-1.6288039) q[1];
sx q[1];
rz(-1.0572546) q[1];
sx q[1];
rz(1.5527976) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88547546) q[0];
sx q[0];
rz(-1.4909571) q[0];
sx q[0];
rz(0.022932963) q[0];
x q[1];
rz(1.2852766) q[2];
sx q[2];
rz(-1.3034045) q[2];
sx q[2];
rz(-1.3533139) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5918436) q[1];
sx q[1];
rz(-1.3938245) q[1];
sx q[1];
rz(-0.66284499) q[1];
rz(-pi) q[2];
rz(-0.62726043) q[3];
sx q[3];
rz(-1.2377487) q[3];
sx q[3];
rz(-1.2675389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9396886) q[2];
sx q[2];
rz(-1.0558015) q[2];
sx q[2];
rz(0.73516694) q[2];
rz(-0.9489263) q[3];
sx q[3];
rz(-2.2831423) q[3];
sx q[3];
rz(2.2775876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.3801124) q[0];
sx q[0];
rz(-2.1367456) q[0];
sx q[0];
rz(0.6066221) q[0];
rz(0.78423777) q[1];
sx q[1];
rz(-1.3332858) q[1];
sx q[1];
rz(1.3444791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30596581) q[0];
sx q[0];
rz(-0.84442645) q[0];
sx q[0];
rz(2.4854922) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44353087) q[2];
sx q[2];
rz(-0.54033989) q[2];
sx q[2];
rz(-2.4176554) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0095876) q[1];
sx q[1];
rz(-0.46087206) q[1];
sx q[1];
rz(0.51746093) q[1];
rz(-pi) q[2];
rz(1.7502039) q[3];
sx q[3];
rz(-2.6699237) q[3];
sx q[3];
rz(-0.52667945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6091696) q[2];
sx q[2];
rz(-2.6340941) q[2];
sx q[2];
rz(1.3896821) q[2];
rz(-2.9621647) q[3];
sx q[3];
rz(-0.98589412) q[3];
sx q[3];
rz(2.7348886) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0591902) q[0];
sx q[0];
rz(-0.32408369) q[0];
sx q[0];
rz(-2.3865336) q[0];
rz(-0.45626196) q[1];
sx q[1];
rz(-1.7165963) q[1];
sx q[1];
rz(1.0770575) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9108601) q[0];
sx q[0];
rz(-2.5344779) q[0];
sx q[0];
rz(2.7348628) q[0];
x q[1];
rz(0.12142373) q[2];
sx q[2];
rz(-2.3540263) q[2];
sx q[2];
rz(-1.9483669) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3216599) q[1];
sx q[1];
rz(-2.0256307) q[1];
sx q[1];
rz(-1.8881892) q[1];
x q[2];
rz(-2.9829426) q[3];
sx q[3];
rz(-1.0358184) q[3];
sx q[3];
rz(-1.8381929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.72851744) q[2];
sx q[2];
rz(-2.04144) q[2];
sx q[2];
rz(-2.4033578) q[2];
rz(1.5089367) q[3];
sx q[3];
rz(-1.8176707) q[3];
sx q[3];
rz(-1.1608605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4974834) q[0];
sx q[0];
rz(-1.4583541) q[0];
sx q[0];
rz(0.64895502) q[0];
rz(2.451918) q[1];
sx q[1];
rz(-0.70976218) q[1];
sx q[1];
rz(2.7241657) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48106347) q[0];
sx q[0];
rz(-1.821035) q[0];
sx q[0];
rz(-2.5323243) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.242495) q[2];
sx q[2];
rz(-0.30825492) q[2];
sx q[2];
rz(-2.2876842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46681225) q[1];
sx q[1];
rz(-1.9793452) q[1];
sx q[1];
rz(1.684434) q[1];
x q[2];
rz(0.77885742) q[3];
sx q[3];
rz(-0.53611272) q[3];
sx q[3];
rz(2.6440309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4735585) q[2];
sx q[2];
rz(-1.7201951) q[2];
sx q[2];
rz(-0.045698015) q[2];
rz(-0.33347305) q[3];
sx q[3];
rz(-2.5662751) q[3];
sx q[3];
rz(-3.0157109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4561653) q[0];
sx q[0];
rz(-2.6021155) q[0];
sx q[0];
rz(-1.217655) q[0];
rz(-2.894891) q[1];
sx q[1];
rz(-1.8496937) q[1];
sx q[1];
rz(0.47952476) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96146255) q[0];
sx q[0];
rz(-0.96213522) q[0];
sx q[0];
rz(-2.6070239) q[0];
x q[1];
rz(-0.87534753) q[2];
sx q[2];
rz(-1.1743059) q[2];
sx q[2];
rz(-3.0975395) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.785276) q[1];
sx q[1];
rz(-2.0802167) q[1];
sx q[1];
rz(-1.2356349) q[1];
rz(-0.44879313) q[3];
sx q[3];
rz(-2.6330607) q[3];
sx q[3];
rz(-2.3859346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.55629998) q[2];
sx q[2];
rz(-1.8712529) q[2];
sx q[2];
rz(1.4053819) q[2];
rz(-2.4893238) q[3];
sx q[3];
rz(-0.26572078) q[3];
sx q[3];
rz(-0.71776596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9115059) q[0];
sx q[0];
rz(-1.7510887) q[0];
sx q[0];
rz(-0.44816309) q[0];
rz(-0.74408342) q[1];
sx q[1];
rz(-0.71744812) q[1];
sx q[1];
rz(1.7070028) q[1];
rz(1.4135398) q[2];
sx q[2];
rz(-2.1485211) q[2];
sx q[2];
rz(-1.1358144) q[2];
rz(-2.7884095) q[3];
sx q[3];
rz(-1.0241057) q[3];
sx q[3];
rz(3.1288341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-2.5833997) q[0];
sx q[0];
rz(-0.71891958) q[0];
sx q[0];
rz(-0.67607003) q[0];
rz(-2.8683635) q[1];
sx q[1];
rz(-0.9571119) q[1];
sx q[1];
rz(2.2303384) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34150037) q[0];
sx q[0];
rz(-1.5507292) q[0];
sx q[0];
rz(-1.5930727) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5780294) q[2];
sx q[2];
rz(-0.94524318) q[2];
sx q[2];
rz(0.48239732) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76409066) q[1];
sx q[1];
rz(-1.6030178) q[1];
sx q[1];
rz(-0.5733787) q[1];
rz(-0.1333255) q[3];
sx q[3];
rz(-0.76078868) q[3];
sx q[3];
rz(2.4530792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91997826) q[2];
sx q[2];
rz(-2.6875434) q[2];
sx q[2];
rz(0.81217074) q[2];
rz(-3.0657366) q[3];
sx q[3];
rz(-2.4935738) q[3];
sx q[3];
rz(-2.5465452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49038637) q[0];
sx q[0];
rz(-2.6631329) q[0];
sx q[0];
rz(1.2745717) q[0];
rz(1.5050585) q[1];
sx q[1];
rz(-0.36205629) q[1];
sx q[1];
rz(-1.9777745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464435) q[0];
sx q[0];
rz(-2.4730485) q[0];
sx q[0];
rz(-0.77495452) q[0];
rz(-pi) q[1];
rz(-0.50588079) q[2];
sx q[2];
rz(-1.1846647) q[2];
sx q[2];
rz(-2.1979111) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.84097505) q[1];
sx q[1];
rz(-1.779044) q[1];
sx q[1];
rz(-2.5131701) q[1];
rz(0.31812654) q[3];
sx q[3];
rz(-1.162718) q[3];
sx q[3];
rz(-2.7216743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.85094467) q[2];
sx q[2];
rz(-3.1182351) q[2];
sx q[2];
rz(-1.5604431) q[2];
rz(0.23049878) q[3];
sx q[3];
rz(-2.0932308) q[3];
sx q[3];
rz(-2.7010664) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5480963) q[0];
sx q[0];
rz(-0.31262147) q[0];
sx q[0];
rz(-3.0189959) q[0];
rz(-0.7971881) q[1];
sx q[1];
rz(-1.0120579) q[1];
sx q[1];
rz(0.0880934) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49876346) q[0];
sx q[0];
rz(-2.3173712) q[0];
sx q[0];
rz(2.7902664) q[0];
rz(1.2253907) q[2];
sx q[2];
rz(-2.900335) q[2];
sx q[2];
rz(3.0221107) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.28411394) q[1];
sx q[1];
rz(-1.4339674) q[1];
sx q[1];
rz(1.72422) q[1];
rz(-2.6016079) q[3];
sx q[3];
rz(-1.267588) q[3];
sx q[3];
rz(1.4410415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2587571) q[2];
sx q[2];
rz(-1.5587403) q[2];
sx q[2];
rz(1.3239512) q[2];
rz(-2.6435408) q[3];
sx q[3];
rz(-2.1784454) q[3];
sx q[3];
rz(-2.3948885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8525456) q[0];
sx q[0];
rz(-2.9883224) q[0];
sx q[0];
rz(-2.9943941) q[0];
rz(0.64951253) q[1];
sx q[1];
rz(-1.5107061) q[1];
sx q[1];
rz(-1.46896) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8181179) q[0];
sx q[0];
rz(-2.9184353) q[0];
sx q[0];
rz(-1.5384307) q[0];
rz(-pi) q[1];
rz(-1.1403667) q[2];
sx q[2];
rz(-1.7375611) q[2];
sx q[2];
rz(-2.8362136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.34079724) q[1];
sx q[1];
rz(-2.1873475) q[1];
sx q[1];
rz(2.7113797) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1280414) q[3];
sx q[3];
rz(-2.1078133) q[3];
sx q[3];
rz(-0.90195105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9327675) q[2];
sx q[2];
rz(-2.6471477) q[2];
sx q[2];
rz(0.24173582) q[2];
rz(0.73511165) q[3];
sx q[3];
rz(-1.7416411) q[3];
sx q[3];
rz(-2.476695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6786137) q[0];
sx q[0];
rz(-2.0942056) q[0];
sx q[0];
rz(0.12272923) q[0];
rz(-2.2562476) q[1];
sx q[1];
rz(-2.131999) q[1];
sx q[1];
rz(0.56040323) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7357702) q[0];
sx q[0];
rz(-0.54600958) q[0];
sx q[0];
rz(-1.9178792) q[0];
rz(-pi) q[1];
rz(0.63626115) q[2];
sx q[2];
rz(-0.36379746) q[2];
sx q[2];
rz(-2.3788102) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.082277) q[1];
sx q[1];
rz(-2.6059285) q[1];
sx q[1];
rz(-1.0470864) q[1];
x q[2];
rz(-2.3123829) q[3];
sx q[3];
rz(-2.5253339) q[3];
sx q[3];
rz(-0.63502888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.6178599) q[2];
sx q[2];
rz(-2.3475519) q[2];
sx q[2];
rz(2.9368371) q[2];
rz(-2.2023885) q[3];
sx q[3];
rz(-1.771628) q[3];
sx q[3];
rz(0.088976629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8646249) q[0];
sx q[0];
rz(-0.11700103) q[0];
sx q[0];
rz(-2.4846039) q[0];
rz(0.14343801) q[1];
sx q[1];
rz(-1.5564432) q[1];
sx q[1];
rz(2.4030446) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.448682) q[0];
sx q[0];
rz(-1.686272) q[0];
sx q[0];
rz(0.67489745) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0623018) q[2];
sx q[2];
rz(-1.3710944) q[2];
sx q[2];
rz(-2.8815567) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9061028) q[1];
sx q[1];
rz(-0.74550438) q[1];
sx q[1];
rz(-1.0968614) q[1];
rz(-2.1288217) q[3];
sx q[3];
rz(-0.9415365) q[3];
sx q[3];
rz(2.2037639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2584381) q[2];
sx q[2];
rz(-2.482735) q[2];
sx q[2];
rz(0.0030041791) q[2];
rz(-2.352412) q[3];
sx q[3];
rz(-0.3843669) q[3];
sx q[3];
rz(-1.9743617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.060519144) q[0];
sx q[0];
rz(-0.34927148) q[0];
sx q[0];
rz(-0.14807598) q[0];
rz(0.5332467) q[1];
sx q[1];
rz(-2.8956469) q[1];
sx q[1];
rz(-1.6291133) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027372885) q[0];
sx q[0];
rz(-0.54410579) q[0];
sx q[0];
rz(-2.6065488) q[0];
rz(-pi) q[1];
rz(2.523411) q[2];
sx q[2];
rz(-0.91518171) q[2];
sx q[2];
rz(-2.5822292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.032669205) q[1];
sx q[1];
rz(-1.9100238) q[1];
sx q[1];
rz(3.0134588) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.809172) q[3];
sx q[3];
rz(-2.4052103) q[3];
sx q[3];
rz(0.045698085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6988354) q[2];
sx q[2];
rz(-1.7008702) q[2];
sx q[2];
rz(3.0485145) q[2];
rz(-0.062151521) q[3];
sx q[3];
rz(-2.744894) q[3];
sx q[3];
rz(-1.9183581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1208948) q[0];
sx q[0];
rz(-0.59497213) q[0];
sx q[0];
rz(2.4769532) q[0];
rz(-1.7918034) q[1];
sx q[1];
rz(-2.2638958) q[1];
sx q[1];
rz(2.452449) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1375547) q[0];
sx q[0];
rz(-1.8144286) q[0];
sx q[0];
rz(1.6904657) q[0];
x q[1];
rz(-0.87920636) q[2];
sx q[2];
rz(-1.9636646) q[2];
sx q[2];
rz(0.61516064) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0943567) q[1];
sx q[1];
rz(-2.7942889) q[1];
sx q[1];
rz(-2.7014665) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1647968) q[3];
sx q[3];
rz(-1.6305123) q[3];
sx q[3];
rz(-2.7415558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.050921116) q[2];
sx q[2];
rz(-0.65616578) q[2];
sx q[2];
rz(-1.7655168) q[2];
rz(-1.225166) q[3];
sx q[3];
rz(-0.71198946) q[3];
sx q[3];
rz(0.041697748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(3.0390778) q[0];
sx q[0];
rz(-0.044487655) q[0];
sx q[0];
rz(0.59120375) q[0];
rz(-2.8996331) q[1];
sx q[1];
rz(-1.0958902) q[1];
sx q[1];
rz(2.718149) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71899881) q[0];
sx q[0];
rz(-2.7072577) q[0];
sx q[0];
rz(2.3584764) q[0];
rz(-pi) q[1];
rz(2.9482163) q[2];
sx q[2];
rz(-2.6453553) q[2];
sx q[2];
rz(-1.9691182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.075942) q[1];
sx q[1];
rz(-1.0417176) q[1];
sx q[1];
rz(-0.34543646) q[1];
x q[2];
rz(2.5245111) q[3];
sx q[3];
rz(-1.9071731) q[3];
sx q[3];
rz(1.468889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.55687904) q[2];
sx q[2];
rz(-0.60606474) q[2];
sx q[2];
rz(-1.1663743) q[2];
rz(-1.1276468) q[3];
sx q[3];
rz(-0.96846628) q[3];
sx q[3];
rz(0.52496547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0135076) q[0];
sx q[0];
rz(-2.7078244) q[0];
sx q[0];
rz(0.67310131) q[0];
rz(0.85064864) q[1];
sx q[1];
rz(-1.7885845) q[1];
sx q[1];
rz(-0.32352111) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7009525) q[0];
sx q[0];
rz(-1.8201168) q[0];
sx q[0];
rz(1.9560157) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0925301) q[2];
sx q[2];
rz(-1.8169122) q[2];
sx q[2];
rz(0.1637295) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9266745) q[1];
sx q[1];
rz(-1.2387382) q[1];
sx q[1];
rz(-2.3658793) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5629285) q[3];
sx q[3];
rz(-0.98697013) q[3];
sx q[3];
rz(1.5848094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2181776) q[2];
sx q[2];
rz(-0.18109334) q[2];
sx q[2];
rz(3.0695445) q[2];
rz(-2.1273023) q[3];
sx q[3];
rz(-2.225596) q[3];
sx q[3];
rz(-2.5924276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2608248) q[0];
sx q[0];
rz(-1.4243955) q[0];
sx q[0];
rz(-2.0212174) q[0];
rz(-2.8063759) q[1];
sx q[1];
rz(-0.64246476) q[1];
sx q[1];
rz(-1.1186218) q[1];
rz(1.3041244) q[2];
sx q[2];
rz(-2.647807) q[2];
sx q[2];
rz(-0.37995445) q[2];
rz(-1.5453962) q[3];
sx q[3];
rz(-1.7138154) q[3];
sx q[3];
rz(-2.9830767) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

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
rz(-0.24124423) q[0];
sx q[0];
rz(-0.21114199) q[0];
sx q[0];
rz(-2.733732) q[0];
rz(1.3762228) q[1];
sx q[1];
rz(-1.9150182) q[1];
sx q[1];
rz(3.0768375) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5714343) q[0];
sx q[0];
rz(-1.6862094) q[0];
sx q[0];
rz(0.0053216318) q[0];
rz(-pi) q[1];
rz(0.98528905) q[2];
sx q[2];
rz(-1.3006214) q[2];
sx q[2];
rz(-2.5463242) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.25542361) q[1];
sx q[1];
rz(-1.5597289) q[1];
sx q[1];
rz(3.0003285) q[1];
rz(2.9969941) q[3];
sx q[3];
rz(-0.62718348) q[3];
sx q[3];
rz(0.36563624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.18225081) q[2];
sx q[2];
rz(-1.094341) q[2];
sx q[2];
rz(1.135929) q[2];
rz(-0.7302537) q[3];
sx q[3];
rz(-1.7466931) q[3];
sx q[3];
rz(1.5203389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092904329) q[0];
sx q[0];
rz(-0.95656675) q[0];
sx q[0];
rz(3.0862869) q[0];
rz(2.7703908) q[1];
sx q[1];
rz(-0.35938811) q[1];
sx q[1];
rz(0.41195437) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13875554) q[0];
sx q[0];
rz(-1.635972) q[0];
sx q[0];
rz(2.6145593) q[0];
rz(1.9048334) q[2];
sx q[2];
rz(-1.8530117) q[2];
sx q[2];
rz(-1.7441234) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1671092) q[1];
sx q[1];
rz(-2.2979263) q[1];
sx q[1];
rz(1.760046) q[1];
rz(-0.92276908) q[3];
sx q[3];
rz(-1.2639746) q[3];
sx q[3];
rz(1.1359147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5281795) q[2];
sx q[2];
rz(-1.4643022) q[2];
sx q[2];
rz(-2.2589653) q[2];
rz(1.4334076) q[3];
sx q[3];
rz(-1.9288758) q[3];
sx q[3];
rz(0.17361704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16647896) q[0];
sx q[0];
rz(-0.011703165) q[0];
sx q[0];
rz(-0.56667462) q[0];
rz(0.4176248) q[1];
sx q[1];
rz(-1.6028701) q[1];
sx q[1];
rz(1.8142987) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1786363) q[0];
sx q[0];
rz(-1.6554852) q[0];
sx q[0];
rz(0.87015193) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1955117) q[2];
sx q[2];
rz(-1.6951121) q[2];
sx q[2];
rz(0.92783463) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3440888) q[1];
sx q[1];
rz(-0.3080536) q[1];
sx q[1];
rz(-1.5665685) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0653082) q[3];
sx q[3];
rz(-2.3479241) q[3];
sx q[3];
rz(2.6738559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3920445) q[2];
sx q[2];
rz(-1.0907402) q[2];
sx q[2];
rz(-2.2331179) q[2];
rz(0.79489094) q[3];
sx q[3];
rz(-1.1029693) q[3];
sx q[3];
rz(-3.095043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4817151) q[0];
sx q[0];
rz(-2.2292723) q[0];
sx q[0];
rz(-2.4758441) q[0];
rz(-2.2944229) q[1];
sx q[1];
rz(-2.8916292) q[1];
sx q[1];
rz(2.7319103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9315791) q[0];
sx q[0];
rz(-2.0635475) q[0];
sx q[0];
rz(-1.3973622) q[0];
rz(1.1944067) q[2];
sx q[2];
rz(-0.79845286) q[2];
sx q[2];
rz(1.2898766) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.078881513) q[1];
sx q[1];
rz(-2.1797175) q[1];
sx q[1];
rz(2.8669861) q[1];
rz(2.014854) q[3];
sx q[3];
rz(-1.5967301) q[3];
sx q[3];
rz(1.3249205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4117671) q[2];
sx q[2];
rz(-0.7080141) q[2];
sx q[2];
rz(1.7626308) q[2];
rz(1.1369368) q[3];
sx q[3];
rz(-1.7863019) q[3];
sx q[3];
rz(2.304145) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9627422) q[0];
sx q[0];
rz(-1.0609635) q[0];
sx q[0];
rz(-2.7233997) q[0];
rz(-1.4232945) q[1];
sx q[1];
rz(-1.9530674) q[1];
sx q[1];
rz(-1.4994015) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2432775) q[0];
sx q[0];
rz(-2.0137798) q[0];
sx q[0];
rz(1.2994631) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4046481) q[2];
sx q[2];
rz(-0.43790753) q[2];
sx q[2];
rz(2.2543668) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.027335701) q[1];
sx q[1];
rz(-1.3523977) q[1];
sx q[1];
rz(3.0824003) q[1];
x q[2];
rz(-1.4410327) q[3];
sx q[3];
rz(-1.5059744) q[3];
sx q[3];
rz(-0.45634064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8666009) q[2];
sx q[2];
rz(-1.9860705) q[2];
sx q[2];
rz(-2.9633813) q[2];
rz(-0.37402672) q[3];
sx q[3];
rz(-2.1686797) q[3];
sx q[3];
rz(-1.9988352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6571534) q[0];
sx q[0];
rz(-1.3244119) q[0];
sx q[0];
rz(1.4763747) q[0];
rz(2.5560675) q[1];
sx q[1];
rz(-0.89591566) q[1];
sx q[1];
rz(-0.73385986) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7134345) q[0];
sx q[0];
rz(-2.1913986) q[0];
sx q[0];
rz(1.0125068) q[0];
x q[1];
rz(-0.30690212) q[2];
sx q[2];
rz(-1.8820805) q[2];
sx q[2];
rz(-2.5800623) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.1230525) q[1];
sx q[1];
rz(-0.55852671) q[1];
sx q[1];
rz(-1.9473507) q[1];
x q[2];
rz(1.1185944) q[3];
sx q[3];
rz(-1.7590932) q[3];
sx q[3];
rz(-1.888243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.057498589) q[2];
sx q[2];
rz(-1.7666662) q[2];
sx q[2];
rz(-1.499184) q[2];
rz(1.5205787) q[3];
sx q[3];
rz(-0.64783827) q[3];
sx q[3];
rz(2.9852941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3868189) q[0];
sx q[0];
rz(-2.2235625) q[0];
sx q[0];
rz(1.3936438) q[0];
rz(-0.58760324) q[1];
sx q[1];
rz(-1.5005451) q[1];
sx q[1];
rz(-0.29744068) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35703941) q[0];
sx q[0];
rz(-1.926657) q[0];
sx q[0];
rz(2.3054992) q[0];
x q[1];
rz(0.018902773) q[2];
sx q[2];
rz(-1.3002965) q[2];
sx q[2];
rz(-2.8342385) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45876452) q[1];
sx q[1];
rz(-1.5611083) q[1];
sx q[1];
rz(0.89114706) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24516205) q[3];
sx q[3];
rz(-0.94221598) q[3];
sx q[3];
rz(-0.50776428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3191159) q[2];
sx q[2];
rz(-1.7812704) q[2];
sx q[2];
rz(2.759867) q[2];
rz(-0.60025674) q[3];
sx q[3];
rz(-0.67265066) q[3];
sx q[3];
rz(0.265358) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47939077) q[0];
sx q[0];
rz(-1.6265765) q[0];
sx q[0];
rz(1.8040682) q[0];
rz(-3.1321101) q[1];
sx q[1];
rz(-0.9597221) q[1];
sx q[1];
rz(-1.7705852) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4776896) q[0];
sx q[0];
rz(-0.36720095) q[0];
sx q[0];
rz(2.7504671) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2393537) q[2];
sx q[2];
rz(-1.9752562) q[2];
sx q[2];
rz(-0.33360937) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7035455) q[1];
sx q[1];
rz(-2.3444972) q[1];
sx q[1];
rz(0.17532562) q[1];
rz(-0.41402702) q[3];
sx q[3];
rz(-1.4883057) q[3];
sx q[3];
rz(0.02007881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3235772) q[2];
sx q[2];
rz(-0.90138268) q[2];
sx q[2];
rz(0.14249194) q[2];
rz(2.2583708) q[3];
sx q[3];
rz(-1.3775237) q[3];
sx q[3];
rz(-1.4628598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.091081) q[0];
sx q[0];
rz(-0.091878042) q[0];
sx q[0];
rz(1.2485414) q[0];
rz(2.7342791) q[1];
sx q[1];
rz(-1.3469603) q[1];
sx q[1];
rz(1.1936845) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3486987) q[0];
sx q[0];
rz(-1.3137709) q[0];
sx q[0];
rz(-1.2299163) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75908179) q[2];
sx q[2];
rz(-2.1164701) q[2];
sx q[2];
rz(2.1496043) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15577573) q[1];
sx q[1];
rz(-2.0627506) q[1];
sx q[1];
rz(1.953275) q[1];
rz(3.1415588) q[3];
sx q[3];
rz(-1.8305998) q[3];
sx q[3];
rz(2.3855748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0672037) q[2];
sx q[2];
rz(-0.95778242) q[2];
sx q[2];
rz(-2.0286782) q[2];
rz(-3.0985966) q[3];
sx q[3];
rz(-1.4837416) q[3];
sx q[3];
rz(-0.23785166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064903108) q[0];
sx q[0];
rz(-1.7273128) q[0];
sx q[0];
rz(-2.9458556) q[0];
rz(1.9211357) q[1];
sx q[1];
rz(-2.7595322) q[1];
sx q[1];
rz(0.97741309) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.138151) q[0];
sx q[0];
rz(-1.312698) q[0];
sx q[0];
rz(-1.4679704) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9896792) q[2];
sx q[2];
rz(-1.8679233) q[2];
sx q[2];
rz(-0.10694416) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4949242) q[1];
sx q[1];
rz(-2.1416602) q[1];
sx q[1];
rz(-2.5350556) q[1];
rz(-pi) q[2];
rz(0.39267003) q[3];
sx q[3];
rz(-1.304198) q[3];
sx q[3];
rz(0.50470282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1310554) q[2];
sx q[2];
rz(-0.26630339) q[2];
sx q[2];
rz(-1.6923426) q[2];
rz(2.2629755) q[3];
sx q[3];
rz(-1.0613469) q[3];
sx q[3];
rz(-1.7884458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6286248) q[0];
sx q[0];
rz(-1.1911363) q[0];
sx q[0];
rz(-1.3918899) q[0];
rz(-1.3799008) q[1];
sx q[1];
rz(-1.5329783) q[1];
sx q[1];
rz(3.0402532) q[1];
rz(-0.78892563) q[2];
sx q[2];
rz(-0.25392214) q[2];
sx q[2];
rz(1.8352652) q[2];
rz(2.8158549) q[3];
sx q[3];
rz(-0.67254638) q[3];
sx q[3];
rz(-1.2253958) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];

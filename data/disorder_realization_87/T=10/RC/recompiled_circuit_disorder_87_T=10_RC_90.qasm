OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23624578) q[0];
sx q[0];
rz(-2.4155004) q[0];
sx q[0];
rz(0.2015764) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(-0.5402686) q[1];
sx q[1];
rz(-0.93710605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1256589) q[0];
sx q[0];
rz(-2.4056245) q[0];
sx q[0];
rz(0.65471162) q[0];
x q[1];
rz(-2.4891698) q[2];
sx q[2];
rz(-2.4180275) q[2];
sx q[2];
rz(-3.0384118) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9846676) q[1];
sx q[1];
rz(-2.942454) q[1];
sx q[1];
rz(1.9763293) q[1];
rz(-1.3017544) q[3];
sx q[3];
rz(-1.6645414) q[3];
sx q[3];
rz(-2.4274488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9822838) q[2];
sx q[2];
rz(-0.26580492) q[2];
sx q[2];
rz(3.0554331) q[2];
rz(-2.384095) q[3];
sx q[3];
rz(-2.3870654) q[3];
sx q[3];
rz(-1.0936201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57698292) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(-1.3264867) q[0];
rz(-1.8857229) q[1];
sx q[1];
rz(-1.5763667) q[1];
sx q[1];
rz(-0.27145162) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46389929) q[0];
sx q[0];
rz(-1.9672183) q[0];
sx q[0];
rz(1.0248313) q[0];
rz(-0.51606744) q[2];
sx q[2];
rz(-1.8415673) q[2];
sx q[2];
rz(2.5004636) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8703692) q[1];
sx q[1];
rz(-1.1961812) q[1];
sx q[1];
rz(0.80925525) q[1];
x q[2];
rz(0.80941697) q[3];
sx q[3];
rz(-1.8339694) q[3];
sx q[3];
rz(-1.1337048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0945956) q[2];
sx q[2];
rz(-1.3826028) q[2];
sx q[2];
rz(2.5644152) q[2];
rz(-2.2180637) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(1.4097759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24401027) q[0];
sx q[0];
rz(-2.0799461) q[0];
sx q[0];
rz(2.999021) q[0];
rz(-1.3525195) q[1];
sx q[1];
rz(-1.0357772) q[1];
sx q[1];
rz(-2.9325063) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3269743) q[0];
sx q[0];
rz(-2.1332392) q[0];
sx q[0];
rz(2.4787865) q[0];
x q[1];
rz(-1.1582583) q[2];
sx q[2];
rz(-0.87450714) q[2];
sx q[2];
rz(1.0227026) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0774539) q[1];
sx q[1];
rz(-1.1541379) q[1];
sx q[1];
rz(2.0616848) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7861869) q[3];
sx q[3];
rz(-1.8672018) q[3];
sx q[3];
rz(0.6682369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3645939) q[2];
sx q[2];
rz(-2.2475188) q[2];
sx q[2];
rz(-2.7374632) q[2];
rz(-1.8557619) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(2.9742441) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7053112) q[0];
sx q[0];
rz(-0.10665882) q[0];
sx q[0];
rz(2.1787815) q[0];
rz(-0.46936938) q[1];
sx q[1];
rz(-2.5517187) q[1];
sx q[1];
rz(-0.00096360047) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9199333) q[0];
sx q[0];
rz(-0.98826212) q[0];
sx q[0];
rz(2.4525053) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12400603) q[2];
sx q[2];
rz(-1.5713072) q[2];
sx q[2];
rz(2.309547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0469989) q[1];
sx q[1];
rz(-1.7013036) q[1];
sx q[1];
rz(-0.99884896) q[1];
rz(-pi) q[2];
rz(3.0022995) q[3];
sx q[3];
rz(-2.5677498) q[3];
sx q[3];
rz(1.3822671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0659539) q[2];
sx q[2];
rz(-1.5116189) q[2];
sx q[2];
rz(0.11165079) q[2];
rz(2.3305437) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(0.013899175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1902996) q[0];
sx q[0];
rz(-0.90068978) q[0];
sx q[0];
rz(0.35650373) q[0];
rz(-0.50645343) q[1];
sx q[1];
rz(-1.3566147) q[1];
sx q[1];
rz(2.8809663) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60359943) q[0];
sx q[0];
rz(-1.7126181) q[0];
sx q[0];
rz(1.5285368) q[0];
rz(-pi) q[1];
rz(-0.25116253) q[2];
sx q[2];
rz(-1.3248982) q[2];
sx q[2];
rz(-0.22566251) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.063559859) q[1];
sx q[1];
rz(-1.5531925) q[1];
sx q[1];
rz(2.5672008) q[1];
x q[2];
rz(0.14560933) q[3];
sx q[3];
rz(-1.429261) q[3];
sx q[3];
rz(-1.6573997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9115209) q[2];
sx q[2];
rz(-1.6610961) q[2];
sx q[2];
rz(0.22932209) q[2];
rz(0.54245943) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(2.1054161) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3689573) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(-1.649958) q[0];
rz(1.0391327) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(1.7274436) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6468069) q[0];
sx q[0];
rz(-1.2959058) q[0];
sx q[0];
rz(-1.5216212) q[0];
rz(-pi) q[1];
rz(-2.6944689) q[2];
sx q[2];
rz(-2.5782055) q[2];
sx q[2];
rz(2.2396357) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1097088) q[1];
sx q[1];
rz(-2.0307396) q[1];
sx q[1];
rz(2.1214532) q[1];
rz(-pi) q[2];
rz(-3.066091) q[3];
sx q[3];
rz(-1.3550948) q[3];
sx q[3];
rz(-2.2942033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.002939) q[2];
sx q[2];
rz(-0.61855519) q[2];
sx q[2];
rz(3.1138528) q[2];
rz(2.6489143) q[3];
sx q[3];
rz(-1.490482) q[3];
sx q[3];
rz(-1.3940575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.7572927) q[0];
sx q[0];
rz(-1.5526271) q[0];
sx q[0];
rz(-0.12167715) q[0];
rz(-1.9901468) q[1];
sx q[1];
rz(-0.45184389) q[1];
sx q[1];
rz(2.7391403) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50169045) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(2.7946266) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5791113) q[2];
sx q[2];
rz(-1.6168211) q[2];
sx q[2];
rz(-0.41551057) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9060045) q[1];
sx q[1];
rz(-1.8530122) q[1];
sx q[1];
rz(-1.3868335) q[1];
rz(-3.1137755) q[3];
sx q[3];
rz(-1.9189315) q[3];
sx q[3];
rz(-2.5607525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1272614) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(-2.2793615) q[2];
rz(0.47752738) q[3];
sx q[3];
rz(-1.1815485) q[3];
sx q[3];
rz(2.2235218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.35571337) q[0];
sx q[0];
rz(-2.9861351) q[0];
sx q[0];
rz(-2.4086337) q[0];
rz(-0.14006242) q[1];
sx q[1];
rz(-2.1439794) q[1];
sx q[1];
rz(1.0345116) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1326133) q[0];
sx q[0];
rz(-1.4011369) q[0];
sx q[0];
rz(-1.3404113) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6413692) q[2];
sx q[2];
rz(-1.6732054) q[2];
sx q[2];
rz(0.18825738) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.58971436) q[1];
sx q[1];
rz(-1.0052181) q[1];
sx q[1];
rz(-3.0977071) q[1];
rz(-pi) q[2];
rz(-1.6789867) q[3];
sx q[3];
rz(-1.9868317) q[3];
sx q[3];
rz(-0.070852208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90116477) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(-2.365716) q[2];
rz(0.72426978) q[3];
sx q[3];
rz(-0.80562076) q[3];
sx q[3];
rz(1.7075214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6999321) q[0];
sx q[0];
rz(-0.76403809) q[0];
sx q[0];
rz(3.124776) q[0];
rz(-0.018521221) q[1];
sx q[1];
rz(-1.3341981) q[1];
sx q[1];
rz(2.3628078) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9176365) q[0];
sx q[0];
rz(-2.7178239) q[0];
sx q[0];
rz(-2.5013322) q[0];
x q[1];
rz(1.0672827) q[2];
sx q[2];
rz(-0.52769606) q[2];
sx q[2];
rz(0.50349456) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15021819) q[1];
sx q[1];
rz(-2.7017936) q[1];
sx q[1];
rz(-0.99779731) q[1];
x q[2];
rz(1.2905144) q[3];
sx q[3];
rz(-2.9119769) q[3];
sx q[3];
rz(1.5667825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6918216) q[2];
sx q[2];
rz(-1.3042973) q[2];
sx q[2];
rz(-0.71643913) q[2];
rz(1.4572432) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(-1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.107782) q[0];
sx q[0];
rz(-0.53492117) q[0];
sx q[0];
rz(-2.9901436) q[0];
rz(2.9653213) q[1];
sx q[1];
rz(-1.9468032) q[1];
sx q[1];
rz(0.72296468) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82397126) q[0];
sx q[0];
rz(-2.9569607) q[0];
sx q[0];
rz(2.1872107) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1935913) q[2];
sx q[2];
rz(-1.4133487) q[2];
sx q[2];
rz(-0.27062705) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.42230095) q[1];
sx q[1];
rz(-2.0740168) q[1];
sx q[1];
rz(-3.0784688) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39977269) q[3];
sx q[3];
rz(-2.3120566) q[3];
sx q[3];
rz(-1.5012036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.67939776) q[2];
sx q[2];
rz(-1.2827736) q[2];
sx q[2];
rz(-0.52552137) q[2];
rz(-2.8578791) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(0.39653683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5491966) q[0];
sx q[0];
rz(-2.017673) q[0];
sx q[0];
rz(2.8155433) q[0];
rz(2.0422968) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(-2.2407871) q[2];
sx q[2];
rz(-1.9566036) q[2];
sx q[2];
rz(-1.4835139) q[2];
rz(1.3463734) q[3];
sx q[3];
rz(-1.8891469) q[3];
sx q[3];
rz(-0.3286152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
